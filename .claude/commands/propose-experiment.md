# propose-experiment

Analyse the latest experiment results and propose the next experiment. Set up the experiment folder, send an email proposal, and arm the reply-checking daemon.

**Important:** You may be running autonomously from the launchd daemon (via `claude --print`). In that case there is no user present to approve actions. Do not ask for confirmation before writing files, creating folders, or sending email — all necessary permissions are pre-granted in `~/.claude/settings.json`. Proceed directly through all steps.

## Steps

1. **Read context**
   - Read `models/experiments/` to find the highest-numbered existing experiment (call it N).
   - **Determine whether experiment N has been run:** check if `data/results/<out_dir_of_N>/evaluation.txt` exists. Record this as `N_has_run` (true/false). You can find `out_dir` in `models/experiments/N/config.yaml`.
   - If `N_has_run` is true: the last completed experiment is N, and the next proposal will be N+1.
   - If `N_has_run` is false: experiment N is proposed but not yet run. The last completed experiment is N-1 (find its `out_dir` from `models/experiments/N-1/config.yaml`).
   - Read the last completed experiment's `evaluation.txt`, `config.yaml`, and `README.md`.
   - If `tools/pending_feedback.txt` exists (a flag file, not the content), fetch the feedback by running:
     ```bash
     .venv/bin/python tools/check_reply.py --msg-id-file tools/.last_proposal_msgid --print-body
     ```
     This returns the reply body once and nowhere else — it is not stored on disk.

2. **Analyse results**
   Interpret the last completed experiment's evaluation using the GUIDANCE section at the top of evaluation.txt:
   - Which genes have the worst FNR?
   - What do the FN p50 CRR values tell you? (>> 1.0 = HMM issue, ≈ 1.0 = weak signal)
   - Are there patterns by population?
   - What does high missingness delta suggest?
   Identify the single most impactful change to make.

3. **Update the last completed experiment's README with actual results**
   Update `models/experiments/<last_completed>/README.md`:
   - Change `Status: Proposed` → `Status: Complete <date>`
   - Add an **Actual outcome** section with the key metrics table (Gene / FNR / PPV / MCC)
   - Note where predictions matched and where they diverged, and why
   - Add a one-line pointer to the next experiment (if known)

   **Do not skip this step.** If the last completed experiment's README already has an Actual outcome section, verify it is accurate and up to date; update it if needed.

4. **Design the experiment**

   **CRITICAL — EXPERIMENT NUMBERING RULE. Read before proceeding:**

   - **If this is NOT a feedback revision** (no `pending_feedback.txt` existed) **OR if `N_has_run` is true**: design a new experiment N+1. Create a new folder `models/experiments/N+1/`. This is the normal post-completion path.
   - **If this IS a feedback revision** (`pending_feedback.txt` existed) **AND `N_has_run` is false**: you are revising a proposed-but-not-yet-run experiment. **Do NOT create experiment N+1. Do NOT change the experiment number.** Modify `models/experiments/N/` in place (update `config.yaml`, `run.sh`, `README.md`). The experiment folder and number stay as N. Proceed to step 5 with N.

   Design principles (for whichever case applies):
   - Keep changes minimal — one or two parameter knobs at most unless the analysis clearly demands more.
   - Reuse existing versioned components unless the algorithm itself must change.
   - If the architecture or a versioned component must change, create a new numbered file (e.g. `models/hmm/02_*.py`).
   - If only parameters change, just update the config.

5. **Create or update the experiment folder**
   - If creating new (N+1): copy `models/experiments/N/` to `models/experiments/N+1/`.
   - If revising in place (N): modify the existing `models/experiments/N/` directly.
   - Update `config.yaml` with the new parameters, including comments explaining the rationale for each change.
   - **Update `out_dir` and names to reflect the actual approach.** If the experiment direction changed significantly (e.g. the original proposed a parameter tweak but the revision pivots to a new architecture), rename `out_dir` in `config.yaml` and the experiment title in `README.md` to describe what the experiment actually does — not the original suggestion. Stale names like `08_higher_confidence` for an experiment that trains a new VAE with dropout are confusing.
   - Update `run.sh` — if reusing the existing checkpoint, invoke `wrap_up.py` with the checkpoint path instead of `train.py`.
   - Write `README.md` with: hypothesis, what changed, what you expect to happen.
   - If this is a **revised proposal** (i.e. `tools/pending_feedback.txt` existed), add or update a **Proposal history** section to the README:
     - Summarise what the previous proposal was (parameters and rationale).
     - Summarise the feedback received (paraphrase — do not paste verbatim; feedback is never stored on disk).
     - Describe what changed between the previous and revised proposal, and why.

6. **Write the email (≤ 100 lines)**
   Write to a temp file `tools/.proposal_body.txt` with this structure:
   ```
   Hi Chiyun,

   Experiment N results summary:
   <2-3 bullet points on key metrics>

   My diagnosis:
   <1 paragraph — what the numbers tell us>

   Proposed experiment <N>: <one-line title>
   Changes from <previous>:
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

   **CRITICAL — EMAIL THREADING RULE. This is mandatory, not optional:**

   If `tools/.proposal_thread_msgid` exists, you are in a feedback revision chain. You **MUST** thread the email under the original proposal — regardless of whether the experiment number or direction changed. Use:
   ```bash
   .venv/bin/python tools/send_email.py \
       --subject "Re: CNV Experiment <N> Proposal" \
       --body @tools/.proposal_body.txt \
       --save-id tools/.last_proposal_msgid \
       --in-reply-to "$(cat tools/.proposal_thread_msgid)"
   ```
   **Do NOT send a fresh email. Do NOT omit `--in-reply-to`.** The user must receive the revised proposal in the same thread as the original proposal and previous acks. Starting a new thread is a bug.

   If `tools/.proposal_thread_msgid` does NOT exist (fresh proposal after an experiment completes):
   ```bash
   .venv/bin/python tools/send_email.py \
       --subject "CNV Experiment <N> Proposal" \
       --body @tools/.proposal_body.txt \
       --save-id tools/.last_proposal_msgid
   ```
   Then:
   ```bash
   echo "<N>" > tools/.last_proposal_experiment
   ```
   If the daemon is not yet installed:
   ```bash
   bash tools/install_daemon.sh
   ```
   If already installed, it will already be polling.

8. **Update README.md** in `models/experiments/<proposed_experiment>/` with the proposal summary and date.

9. **Clean up** `tools/pending_feedback.txt` if it existed.
   (When running via the daemon, the daemon removes this file after Claude exits.
   When running interactively, remove it here so the daemon doesn't stay blocked.)
   Do NOT remove `tools/.proposal_thread_msgid` — it must persist across multiple
   feedback rounds so all revisions stay in the same email chain. The daemon clears
   it when an AUTHORISE is received.
