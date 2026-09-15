# Reproducing the manuscript

Every table, figure and headline number in the paper is produced by the scripts listed below.
Intermediate tables are written to `reports/publication/`; figures to `figures/`.

## Environment

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e .
pip install evo2 boltz          # GPU inference engines
docker compose up               # Evo 2 + Boltz-2 servers (see docker-compose.yml)
```

Inputs (reference genome, wild-type structures, ligand) are downloaded automatically by
`scripts/init_project.py` when `data_autofetch.enabled: true`.

## Pipeline: building BDQ-SatDB

| Step | Command | What it produces |
|---|---|---|
| 1 | `python scripts/init_project.py` | Every single-nucleotide substitution over the CDS and the defined regulatory windows of the 10 target genes; per-gene mutational space (Table S1) |
| 2 | `python run_pipeline.py --stage evo` | Evo 2 differential triage (ΔEvo) → `LETHAL_SKIP` / `LOSS_OF_FUNCTION` / `BOLTZ_READY` |
| 3 | `python run_pipeline.py --stage boltz` | Boltz-2 structural re-examination of `BOLTZ_READY` candidates (pTM, ipTM, pLDDT, stability, clash) |
| 4 | `python run_pipeline.py --stage analyze` | Rule layer: regulator–effector epistasis, multi-gene synergy rules, per-variant interpretation and CSV export |

## Manuscript item → script

| Manuscript item | Script(s) | Result table |
|---|---|---|
| Table S1, gene classification and mutational space | `scripts/init_project.py`, `config/config.yaml` | per-gene counts in `reports/publication/` |
| Saturation library totals (49,410 = 36,990 coding + 12,420 regulatory) | `scripts/export_satdb_overview_stats.py` | — |
| Coding-region resistance candidates and tiers (1,643; 155 / 134 / 1,354) | `scripts/export_saturation_R_tiers.py` | `saturation_R_tiers_20260914.csv` |
| Catalogue-absent hypotheses (virtual mutagenesis tables) | `scripts/export_virtual_mutation_publication_tables.py` | `virtual_mutation_*.csv`, `regulatory_example_assessment.csv` |
| WHO catalogue benchmark (1,465 entries; 461 primary evaluation entries) | `scripts/prepare_who_evaluation_benchmark.py`, `scripts/prepare_who_resistant_benchmark.py`, `scripts/evaluate_who_catalog.py`, `scripts/detailed_who_analysis.py` | `who_bdq_per_row_benchmark*.csv` |
| WHO summary / stratified tables (Table 1) | `scripts/export_who_bdq_publication_tables.py`, `scripts/export_who_bdq_unified_benchmark_tables.py`, `scripts/export_who_bdq_stratified_tables.py`, `scripts/export_who_bdq_lof_augmented_tables.py` | `who_bdq_*summary*.csv`, `who_bdq_stratified_*.csv` |
| CRyPTIC performance table and confidence intervals (Table S3) | `scripts/export_cryptic_manuscript_performance_table.py`, `scripts/export_cryptic_manuscript_performance_ci.py` | `cryptic_manuscript_performance_*.csv` |
| Benchmarking against TB-Profiler and catalogue matching | `scripts/benchmark_tbprofiler.py`, `scripts/benchmark_catalog_match.py` | — |
| Paired / McNemar comparisons on the common isolate set | `scripts/export_cryptic_common_sample_paired_stats.py`, `scripts/benchmark_mcnemar.py` | `cryptic_common_sample_paired_stats.csv` |
| Figure 2A–E (CRyPTIC validation panels) | `scripts/plot_cryptic_validation_figures.py`, `scripts/make_fig2_combined.py` | `figures/fig2_combined.png/pdf` |
| Figure 2B/C (false-negative taxonomy, false-positive composition) | `scripts/export_cryptic_fn_taxonomy.py`, `scripts/export_cryptic_fn_mechanism_expansion_audit.py` | `cryptic_fn_taxonomy.csv`, `cryptic_fn_mic_expansion_audit.csv` |
| Figure 2D (mtrB combination enrichment) | `scripts/export_cryptic_mtrb_enrichment.py` | `cryptic_mtrb_enrichment_summary.csv`, `cryptic_mtrb_combo_summary.csv` |
| Figure 2E (supervised genotype-only ceiling) | `scripts/export_cryptic_supervised_upper_bound.py`, `scripts/run_cryptic_supervised_site_holdout.py`, `scripts/summarize_cryptic_supervised_robustness.py` | `cryptic_supervised_*.csv` |
| Figure 2F/G (MIC-level false positives; sensitivity by MIC) | `scripts/plot_fp_mic_shift.py`, `scripts/export_cryptic_five_level_calibration.py` | `cryptic_five_level_*.csv` |
| Indeterminate / unknown-isolate sensitivity analysis | `scripts/export_cryptic_indeterminate_sensitivity.py` | `cryptic_indeterminate_sensitivity_*.csv` |
| Table S2 and the internal cohort section | `scripts/export_internal_cohort_tables.py` | `variant_calls_final.csv`, `isolate_calls_final.csv` |
| Internal cohort figures (five-category re-interpretation) | `scripts/plot_local_cohort_figure.py`, `scripts/plot_local_alluvial_figure.py`, `scripts/refresh_local_five_category.py` | `figures/` |
| Rule matrix used by reviewers | `ANALYSIS_RULE_MATRIX.md` | — |

Shared helpers: `scripts/cryptic_bdq_validation_common.py`, `scripts/rule_postprocess.py`.

## Notes on interpretation rules

- Variant-level calls: `R` (resistance), `S` (susceptibility), `INDET` (direction not
  determined). Regulatory-window and exploratory-gene variants are recorded as `INDET` and do
  not determine isolate-level calls.
- Isolate-level aggregation is resistance-dominant: any `R` call → resistant; otherwise any `S`
  call → susceptible; isolates with no interpretable call remain unknown. In the internal
  cohort every isolate carried at least one susceptibility call, so no isolate remained
  indeterminate.
- Structural-integrity and sensitization genes can only be sensitizing or neutral, never
  resistance-conferring; their non-synonymous and indel variants are therefore interpreted as
  susceptible.
