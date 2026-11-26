---
alwaysApply: true
---
# Complete execution ledger

Maintain `tmp/EXECUTION_LOG.md` as the chronological, project-local audit log for this task. This file is intentionally ignored scratch output; do not edit the approved plan file.

For every substantive shell command, workflow/test run, generated artifact, or file modification:

1. **Before execution**, append a run entry containing a unique run ID, timestamp, purpose, exact command or tool action, working directory, input files/data/parameters, and expected outputs.
2. **Immediately after execution**, update the same entry with exit/result status, observed stdout/stderr or tool result, actual output/artifact paths, and relevant checksums or row/test counts when applicable. Do not record a run as successful unless its result was observed.
3. Record the **next step**, including any failure diagnosis or follow-up needed to reproduce or resume the work.
4. Preserve enough detail to reconstruct the run. For long output, store the complete output in a separate ignored file under `tmp/run_logs/` and link it from the entry; include the unabbreviated failure or success summary in the ledger. Never silently omit a failed command.
5. Itemize each action in a parallel tool batch. File reads used only for inspection may be grouped by purpose, but every write, code-generation step, download, environment operation, test, benchmark, and data-processing run must have its own entry. Updates to the ledger itself are part of the associated run and do not require separate recursive entries.
6. Never record credentials, tokens, or other secrets. Do not modify bulk `data`, `db`, or `results` merely to create logs.
7. At phase boundaries, add a concise checkpoint summarizing completed run IDs, generated artifacts, unresolved failures, and the exact next action.

The assistant must follow this ledger process in the current task even if rule activation normally begins with a later request.
