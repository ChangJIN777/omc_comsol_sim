---
name: feedback-fab-check-separate
description: Fabricability is its own pre-solver gate (reason codes, YAML thresholds), separate from geometric validity; cite threshold sources
metadata:
  type: feedback
---

For new cell types, keep fabricability as a separate pre-solver gate with its own reason codes and YAML-configurable thresholds, apart from geometric build-validity rules. Trace each threshold to a repo file:line, and mark unsourced ones TODO(confirm).

**Why:** the user explicitly asked for a "check to see if the structure is fabricable" on 2026-11-23, on top of the template's validity rules.

**How to apply:** in any new optimizer or cell-type integration, record failed fab rules in the result record and save those candidates like infeasible ones.

Related: [[project-python-rib-integration]].
