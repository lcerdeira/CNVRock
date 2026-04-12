# propose-experiment

Analyse the latest experiment results and propose the next experiment. Set up the experiment folder, send an email proposal, and arm the reply-checking daemon.

**Important:** You may be running autonomously from the launchd daemon (via `claude --print`). In that case there is no user present to approve actions. Do not ask for confirmation before writing files, creating folders, or sending email — all necessary permissions are pre-granted in `~/.claude/settings.json`. Proceed directly through all steps.

## Steps

1. **Read context**
   - Read `models/experiments/` to find the highest-numbered existing experiment (call it N).
   - Read `data/results/<out_dir_of_experiment_N>/evaluation.txt` for that experiment.
   - Read `models/experiments/N/config.yaml` and `models/experiments/N/README.md`.
   - If `tools/pending_feedback.txt` exists (a flag file, not the content), fetch the feedback by running:
     ```bash
     .venv/bin/python tools/check_reply.py --msg-id-file tools/.last_proposal_msgid --print-body
     ```
     This returns the reply body once and nowhere else — it is not stored on disk.

2. **Analyse results**
   Interpret the evaluation using the GUIDANCE section at the top of evaluation.txt:
   - Which genes have the worst FNR?
   - What do the FN p50 CRR values tell you? (>> 1.0 = HMM issue, ≈ 1.0 = weak signal)
   - Are there patterns by population?
   - What does high missingness delta suggest?
   Identify the single most impactful change to make.

3. **Update experiment N's README with actual results**
   Update `models/experiments/N/README.md` to record what actually happened:
   - Change `Status: Proposed` → `Status: Complete <date>`
   - Add an **Actual outcome** section with the key metrics table (Gene / FNR / PPV / MCC)
   - Note where predictions matched and where they diverged, and why
   - Add a one-line pointer to the next experiment: `→ See experiment N+1`

4. **Design experiment N+1**
   - Keep changes minimal — one or two parameter knobs at most unless the analysis clearly demands more.
   - Reuse existing versioned components unless the algorithm itself must change.
   - If the architecture or a versioned component must change, create a new numbered file (e.g. `models/hmm/02_*.py`).
   - If only parameters change, just update the config.

5. **Create the experiment folder**
   - Copy `models/experiments/N/` to `models/experiments/N+1/`.
   - Update `config.yaml` with the new parameters, including comments explaining the rationale for each change.
   - Update `run.sh` — if reusing the existing checkpoint, invoke `wrap_up.py` with the checkpoint path instead of `train.py`.
   - Write `README.md` with: hypothesis, what changed, what you expect to happen.
   - If this is a **revised proposal** (i.e. `tools/pending_feedback.txt` existed), add a **Proposal history** section to the README:
     - Summarise what your original proposal was (parameters and rationale).
     - Summarise the feedback received (paraphrase — do not paste verbatim; feedback is never stored on disk).
     - Describe what changed between the original and revised proposal, and why.

6. **Write the email (≤ 100 lines)**
   Write to a temp file `tools/.proposal_body.txt` with this structure:
   ```
   Hi Chiyun,

   Experiment N results summary:
   <2-3 bullet points on key metrics>

   My diagnosis:
   <1 paragraph — what the numbers tell us>

   Proposed experiment N+1: <one-line title>
   Changes from N:
   <bullet list of config changes with rationale>

   What I expect:
   <1 paragraph>

   What we could do instead:
   <1-2 alternative approaches if this one doesn't work>

   Run time estimate: <rough estimate based on whether training is needed>

   Reply AUTHORISE to run this on the Mac mini, or reply with feedback
   and I'll revise the proposal.

   Claude
   ```

7. **Send the email and arm the daemon**
   Run:
   ```bash
   .venv/bin/python tools/send_email.py \
       --subject "CNV Experiment N+1 Proposal" \
       --body @tools/.proposal_body.txt \
       --save-id tools/.last_proposal_msgid
   echo "N+1" > tools/.last_proposal_experiment
   ```
   Then, if the daemon is not yet installed:
   ```bash
   bash tools/install_daemon.sh
   ```
   If already installed, it will already be polling.

8. **Update README.md** in `models/experiments/N+1/` with the proposal summary and date.

9. **Clean up** `tools/pending_feedback.txt` if it existed (feedback has now been acted on).
