# check-reply

Check whether a reply to the latest experiment proposal email has arrived, and act on it.

## Steps

1. **Check for a pending proposal**
   If `tools/.last_proposal_msgid` does not exist, say "No proposal is currently pending" and stop.

2. **Check for a reply**
   Run:
   ```bash
   .venv/bin/python tools/check_reply.py --msg-id-file tools/.last_proposal_msgid
   ```

3. **Handle the result**

   **If NO_REPLY:** Say "No reply yet. The daemon will check again in 5 minutes."

   **If AUTHORISE:**
   - Read `tools/.last_proposal_experiment` to get the experiment number.
   - Run the experiment:
     ```bash
     bash models/experiments/<N>/run.sh
     ```
   - When complete, send a summary email:
     ```bash
     .venv/bin/python tools/send_email.py \
         --subject "CNV Experiment <N> — Complete" \
         --body "Experiment <N> has finished. Run /propose-experiment to review results and propose the next one."
     ```
   - Remove `tools/.last_proposal_msgid` and `tools/.last_proposal_experiment`.
   - Update `models/experiments/<N>/README.md` with "Status: Running / Complete".

   **If FEEDBACK:**
   - Fetch the feedback text (only to this session, not to disk):
     ```bash
     .venv/bin/python tools/check_reply.py \
         --msg-id-file tools/.last_proposal_msgid --print-body
     ```
   - Read the feedback, implement the requested changes to the experiment config and README.
   - Then run `/propose-experiment` to send a revised proposal.
