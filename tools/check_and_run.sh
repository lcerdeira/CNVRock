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

if [ ! -f "$MSGID_FILE" ]; then
    log "No proposal pending (no msgid file). Exiting."
    exit 0
fi

log "Checking for reply..."
RESULT=$("$PYTHON" "$TOOLS/check_reply.py" --msg-id-file "$MSGID_FILE" 2>>"$LOG") || true

case "$RESULT" in
    NO_REPLY)
        log "No reply yet."
        ;;

    AUTHORISE)
        log "AUTHORISE received. Running experiment."
        EXPERIMENT=$(cat "$TOOLS/.last_proposal_experiment" 2>/dev/null || echo "02")
        # Clean up before running so a re-trigger can't happen mid-run
        rm -f "$MSGID_FILE" "$TOOLS/.last_proposal_experiment"

        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "CNV Experiment $EXPERIMENT — Running" \
            --body "Experiment $EXPERIMENT has started on the Mac mini. I'll email you when it's done and propose the next one."

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
                --subject "CNV Experiment $EXPERIMENT — Complete" \
                --body "Experiment $EXPERIMENT has finished. Open Claude Code and run /propose-experiment to review the results and propose the next experiment."
        fi
        ;;

    FEEDBACK)
        log "Feedback received (content not logged — open Claude Code to review)."
        # Write a flag file so /propose-experiment knows to read pending feedback.
        # The feedback itself stays in your Gmail — it is never written to disk here.
        touch "$FEEDBACK_FILE"
        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "CNV — Feedback received" \
            --body "Got it. I'll review your feedback and send a revised proposal shortly. Open Claude Code to proceed."
        # Leave msgid file in place until the new proposal replaces it
        ;;

    *)
        log "Unexpected reply check output: $RESULT"
        ;;
esac
