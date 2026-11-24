---
name: comsol-batch-licence-behaviour
description: Measured COMSOL 6.3 licence checkouts on the Harvard FlexLM server for batch vs class-file runs, plus batch-mode gotchas (security sandbox, exit code 0 on error)
metadata:
  type: project
---

Measured 2026-11-23 on 27005@license.seas.harvard.edu with COMSOL 6.3 build 420 and python-scripts/comsol/omc_unitcell_rib.mph (tool: python-scripts/scripts/run_comsol_batch.sh --watch-license):

- `comsolbatch -inputfile X.mph -study std2 -usebatchlic` takes COMSOLBATCH, WAVEOPTICSBATCH, MEMSBATCH, CADIMPORTBATCH, CLUSTERNODE (+ COMSOLUSER, CADIMPORTUSER). No interactive WAVEOPTICS/RF/COMSOL seat. The batch pool is selected by `-usebatchlic`, not by batch mode itself.
- MEMS and CADIMPORT are checked out even for an optics-only study, because the template lists them (`comsolbatch -checklicense file.mph` prints the list without checking anything out). MEMS has only 2 seats.
- A compiled Java class as `-inputfile` cannot use `-usebatchlic` ("Input file does not support batch licenses"); it takes interactive COMSOL, WAVEOPTICS, MEMS, CADIMPORT while it evaluates results.
- COMSOL batch exits 0 after an error banner (`/*****Error********/` in the log), so the log must be checked.
- Java class files run by comsolbatch cannot write files (security preference "File system access"); print to stdout instead, it reaches the caller.
- The user rule: never hold a COMSOL seat while idle (only 3 WAVEOPTICS seats, shared).

**Why:** the user is often locked out of the 3 WAVEOPTICS seats; the *BATCH pool is normally unused.
**How to apply:** route unattended solves through `-usebatchlic`; keep any step that needs interactive seats short and retryable; do not propose keep-alive sessions. Re-measure before relying on it after a COMSOL upgrade or licence renewal.
