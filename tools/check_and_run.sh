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
        # Determine which experiment the proposal was for from a companion file
        EXPERIMENT=$(cat "$TOOLS/.last_proposal_experiment" 2>/dev/null || echo "02")
        cd "$REPO_ROOT/models"
        bash "experiments/$EXPERIMENT/run.sh" >> "$LOG" 2>&1
        log "Experiment $EXPERIMENT complete."
        # Send confirmation email
        "$PYTHON" "$TOOLS/send_email.py" \
            --subject "CNV Experiment $EXPERIMENT — Started" \
            --body "Experiment $EXPERIMENT is now running on the Mac mini. I will update you when it completes."
        # Clean up so the daemon doesn't re-trigger
        rm -f "$MSGID_FILE" "$TOOLS/.last_proposal_experiment"
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
