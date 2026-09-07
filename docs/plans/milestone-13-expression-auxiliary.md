# Milestone 13: Expression auxiliary (after OOF)

> **Alias:** historical **7G″**. Index: [`MILESTONE_INDEX.md`](MILESTONE_INDEX.md).
>
> **Renumber note (2026-09-07):** previously stubbed as Milestone **12**; now
> runs **after** Milestone **12** final OOF. Detail sketch:
> [`milestone-7g-double-prime-expression-auxiliary.md`](milestone-7g-double-prime-expression-auxiliary.md).

**Status:** `deferred`

## Intent

After Milestone **12** produces OOF (or an explicit interim pretrained
checkpoint from **10d**), either:

1. **Continue training** the shared methylation encoder with an expression
   prediction head / auxiliary loss, or
2. **Finetune** a frozen or lightly updated encoder for gene-expression
   prediction on matched methylation–RNA samples.

Goal: test whether expression supervision improves biologically interpretable
promoter/body aggregation beyond phenotype-only training — not a copy of
RSMethy-Net.

## Required plan steps (data first)

1. **Download expression data** into `$MBS_DATA_ROOT` (study-grouped, documented
   source URLs, checksums). Prefer cohorts with known methylation overlap
   (e.g. TCGA / Hub studies with RNA-seq or expression arrays).
2. **Catalog + overlap census** — sample×gene expression table, platform IDs,
   study IDs; intersection with canonical methylation matrices; ≥1k usable
   labeled samples bar before GPU.
3. **Split hygiene** — study- or cancer-grouped folds; never random
   within-cancer CV as the primary claim.
4. **Train / finetune** — continue-train vs finetune arms; mask missing RNA;
   residualized expression control; report phenotype + expression metrics.
5. **Report** under `reports/inspection/stage0_7g_double_prime_expression/`.

## Non-goals

- Not a gate for Milestones **10**, **11**, or **12**.
- Do not block final OOF on expression download.
- Do not replace CascadeDeepSet with per-gene fixed-window encoders.

## Done when

- Expression download + manifests committed under data roots / inspection
  reports.
- Overlap census passes the sample bar (or documents why a smaller pilot is
  accepted).
- Continue-train and/or finetune arms reported with study-held-out metrics.
