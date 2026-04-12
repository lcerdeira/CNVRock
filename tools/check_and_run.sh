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

RESULT=$("$PYTHON" "$TOOLS/check_reply.py" --msg-id-file "$MSGID_FILE" 2>>"$LOG") || true

case "$RESULT" in
    NO_REPLY)
        ;;  # silent — don't clutter the log on every poll

    AUTHORISE)
        log "AUTHORISE received. Running experiment."
        EXPERIMENT=$(cat "$TOOLS/.last_proposal_experiment" 2>/dev/null || echo "02")
        # Capture proposal Message-ID for threading before cleanup
        PROPOSAL_MSGID=$(cat "$MSGID_FILE" 2>/dev/null || true)
        # Clean up before running so a re-trigger can't happen mid-run
        rm -f "$MSGID_FILE" "$TOOLS/.last_proposal_experiment"

        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
            --body "Experiment $EXPERIMENT has started on the Mac mini. I'll email you when it's done and propose the next one." \
            --in-reply-to "$PROPOSAL_MSGID"

        cd "$REPO_ROOT/models"
        bash "experiments/$EXPERIMENT/run.sh" >> "$LOG" 2>&1
        log "Experiment $EXPERIMENT complete."

        # Automatically propose the next experiment using Claude CLI
        if [ -n "$CLAUDE_BIN" ]; then
            log "Invoking Claude to propose next experiment..."
            cd "$REPO_ROOT"
            "$CLAUDE_BIN" --print \
                "$(cat .claude/commands/propose-experiment.md)" \
                >> "$LOG" 2>&1 \
                && log "Next experiment proposed and email sent." \
                || log "Claude proposal step failed — open Claude Code and run /propose-experiment manually."
        else
            log "Claude CLI not found. Open Claude Code and run /propose-experiment to propose the next experiment."
            "$PYTHON" "$TOOLS/send_email.py" \
                --subject "Re: CNV Experiment $EXPERIMENT Proposal" \
                --body "Experiment $EXPERIMENT has finished. Open Claude Code and run /propose-experiment to review the results and propose the next experiment." \
                --in-reply-to "$PROPOSAL_MSGID"
        fi
        ;;

    FEEDBACK)
        log "Feedback received (content not logged — open Claude Code to review)."
        # Write a flag file so /propose-experiment knows to read pending feedback.
        # The feedback itself stays in your Gmail — it is never written to disk here.
        touch "$FEEDBACK_FILE"
        PROPOSAL_MSGID=$(cat "$MSGID_FILE" 2>/dev/null || true)
        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "Re: CNV Experiment $(cat "$TOOLS/.last_proposal_experiment" 2>/dev/null) Proposal" \
            --body "Got it. I'll review your feedback and send a revised proposal shortly. Open Claude Code to proceed." \
            --in-reply-to "$PROPOSAL_MSGID"
        # Leave msgid file in place until the new proposal replaces it
        ;;

    *)
        log "Unexpected reply check output: $RESULT"
        ;;
esac
