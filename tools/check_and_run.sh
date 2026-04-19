#!/bin/bash
# Called by the launchd daemon every few minutes.
# Checks for a reply to the latest experiment proposal email.
# If AUTHORISE → runs the experiment.
# If FEEDBACK  → saves to pending_feedback.txt and sends an acknowledgment.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="$REPO_ROOT/.venv/bin/python"
TOOLS="$REPO_ROOT/tools"
MSGID_FILE="$TOOLS/.last_proposal_msgid"
FEEDBACK_FILE="$TOOLS/pending_feedback.txt"
THREAD_MSGID_FILE="$TOOLS/.proposal_thread_msgid"
LOG="$TOOLS/daemon.log"

# Locate the Claude CLI — prefer a system install, fall back to the VSCode
# extension binary (fragile: breaks on extension updates).
# Run `brew install node && npm install -g @anthropic-ai/claude-code` for a
# stable system install.
CLAUDE_BIN="$(which claude 2>/dev/null \
    || find "$HOME/.vscode/extensions" -name "claude" -type f 2>/dev/null | sort -r | head -1 \
    || true)"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S')  $*" >> "$LOG"; }

# Keep log bounded — trim to last 500 lines if it exceeds 1000
if [ -f "$LOG" ] && [ "$(wc -l < "$LOG")" -gt 1000 ]; then
    tail -500 "$LOG" > "$LOG.tmp" && mv "$LOG.tmp" "$LOG"
fi

if [ ! -f "$MSGID_FILE" ]; then
    exit 0  # nothing pending — silent exit
fi

# If feedback was already acked and the revised proposal was sent (msgid file
# updated), the old user-reply email is still in the inbox, so check_reply.py
# would detect FEEDBACK again on the next poll.  Guard against this by exiting
# early when pending_feedback.txt exists — /propose-experiment removes it once
# the revised proposal email is sent.
if [ -f "$FEEDBACK_FILE" ]; then
    exit 0
fi

RESULT=$("$PYTHON" "$TOOLS/check_reply.py" --msg-id-file "$MSGID_FILE" 2>>"$LOG") || true

case "$RESULT" in
    NO_REPLY)
        ;;  # silent — don't clutter the log on every poll

    AUTHORISE)
        log "AUTHORISE received. Running experiment."
        EXPERIMENT=$(printf "%02d" "$(cat "$TOOLS/.last_proposal_experiment" 2>/dev/null || echo "2")")
        # Capture proposal Message-ID for threading before cleanup
        PROPOSAL_MSGID=$(cat "$MSGID_FILE" 2>/dev/null || true)
        # Clean up before running so a re-trigger can't happen mid-run.
        # Also clear the thread root — the next proposal starts a fresh chain.
        rm -f "$MSGID_FILE" "$TOOLS/.last_proposal_experiment" "$THREAD_MSGID_FILE"

        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
            --body "Experiment $EXPERIMENT has started on the Mac mini. I'll email you when it's done and propose the next one." \
            --in-reply-to "$PROPOSAL_MSGID"

        cd "$REPO_ROOT/models"
        # Send a failure email if run.sh crashes so the loop doesn't die silently.
        if ! bash "experiments/$EXPERIMENT/run.sh" >> "$LOG" 2>&1; then
            log "Experiment $EXPERIMENT FAILED. Sending failure notification."
            "$PYTHON" "$TOOLS/send_email.py" \
                --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
                --body "Experiment $EXPERIMENT failed on the Mac mini. Check tools/daemon.log for the error. The loop is paused — open Claude Code and run /propose-experiment once you've diagnosed the issue." \
                --in-reply-to "$PROPOSAL_MSGID" || true
            exit 1
        fi
        log "Experiment $EXPERIMENT complete."

        # Automatically propose the next experiment using Claude CLI
        if [ -n "$CLAUDE_BIN" ]; then
            log "Invoking Claude to propose next experiment..."
            cd "$REPO_ROOT"
            "$CLAUDE_BIN" --print \
                --permission-mode bypassPermissions \
                "$(cat .claude/commands/propose-experiment.md)" \
                >> "$LOG" 2>&1 || true
            # Verify Claude actually armed the daemon — exit-code alone is not reliable.
            # A successful proposal always writes .last_proposal_msgid.
            if [ -f "$MSGID_FILE" ]; then
                log "Next experiment proposed and email sent."
            else
                log "Claude proposal step did not complete (no msgid file). Sending fallback email."
                "$PYTHON" "$TOOLS/send_email.py" \
                    --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
                    --body "Experiment $EXPERIMENT has finished but the automatic proposal failed. Open Claude Code and run /propose-experiment to review results and propose the next experiment." \
                    --in-reply-to "$PROPOSAL_MSGID"
            fi
        else
            log "Claude CLI not found. Open Claude Code and run /propose-experiment to propose the next experiment."
            "$PYTHON" "$TOOLS/send_email.py" \
                --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
                --body "Experiment $EXPERIMENT has finished. Open Claude Code and run /propose-experiment to review the results and propose the next experiment." \
                --in-reply-to "$PROPOSAL_MSGID"
        fi
        ;;

    FEEDBACK)
        log "Feedback received (content not logged — invoking Claude to revise proposal)."
        # Write flag so /propose-experiment knows to fetch pending feedback.
        touch "$FEEDBACK_FILE"
        # Guarantee cleanup even if set -e fires before we reach the explicit rm below
        # (e.g. if the ack send_email.py call fails with a transient SMTP error).
        trap 'rm -f "$FEEDBACK_FILE"' EXIT

        PROPOSAL_MSGID=$(cat "$MSGID_FILE" 2>/dev/null || true)
        EXPERIMENT=$(printf "%02d" "$(cat "$TOOLS/.last_proposal_experiment" 2>/dev/null || echo "0")")

        # On the first feedback round, record the original proposal as the
        # thread root so all revisions (and their acks) stay in one email chain.
        if [ ! -f "$THREAD_MSGID_FILE" ]; then
            echo "$PROPOSAL_MSGID" > "$THREAD_MSGID_FILE"
        fi
        THREAD_MSGID=$(cat "$THREAD_MSGID_FILE")

        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
            --body "Got it. I'll review your feedback and send a revised proposal shortly." \
            --in-reply-to "$THREAD_MSGID"

        if [ -n "$CLAUDE_BIN" ]; then
            log "Invoking Claude to revise proposal based on feedback..."
            cd "$REPO_ROOT"
            # Timeout after 30 min so a hung Claude process doesn't block the daemon
            # indefinitely (pending_feedback.txt would persist until the process exits).
            # macOS lacks GNU `timeout`, so emulate it with a background watchdog.
            "$CLAUDE_BIN" --print \
                --permission-mode bypassPermissions \
                "$(cat .claude/commands/propose-experiment.md)" \
                >> "$LOG" 2>&1 &
            CLAUDE_PID=$!
            ( sleep 1800; kill -0 "$CLAUDE_PID" 2>/dev/null && kill "$CLAUDE_PID" ) &
            WATCHDOG_PID=$!
            wait "$CLAUDE_PID" 2>/dev/null || true
            kill "$WATCHDOG_PID" 2>/dev/null || true

            # Always remove the flag here — never rely on Claude to clean it up.
            # If Claude succeeded, the flag is redundant. If it failed, leaving it
            # would silently block the daemon from polling on every future run.
            rm -f "$FEEDBACK_FILE"

            if [ -f "$MSGID_FILE" ] && [ "$(cat "$MSGID_FILE")" != "$PROPOSAL_MSGID" ]; then
                log "Revised proposal sent successfully."
            else
                # Claude failed to send a revised proposal. Update the msgid file
                # to point to a fallback email so the next poll doesn't re-detect
                # the same user-reply and send another ack.
                log "Claude did not complete the revised proposal. Sending fallback email."
                "$PYTHON" "$TOOLS/send_email.py" \
                    --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
                    --body "Automatic proposal revision failed. Open Claude Code and run /propose-experiment to revise the proposal." \
                    --in-reply-to "$THREAD_MSGID" \
                    --save-id "$MSGID_FILE"
            fi
        else
            rm -f "$FEEDBACK_FILE"
            log "Claude CLI not found. Sending fallback email."
            "$PYTHON" "$TOOLS/send_email.py" \
                --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
                --body "Claude CLI not found. Open Claude Code and run /propose-experiment to revise the proposal." \
                --in-reply-to "$THREAD_MSGID" \
                --save-id "$MSGID_FILE"
        fi
        ;;

    *)
        log "Unexpected reply check output: $RESULT"
        ;;
esac
