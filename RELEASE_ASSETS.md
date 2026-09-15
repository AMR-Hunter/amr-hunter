# Release assets

Files too large to keep in the git tree are attached to the GitHub release
(tag `v1.0.0`). Download them and place them next to `reports/publication/` to reproduce the
manuscript numbers.

| Asset | Size | Contents |
|---|---|---|
| `amr_hunter_db_final.db.gz` | ~370 MB uncompressed | **BDQ-SatDB**: whole-site saturation-mutagenesis database for the 10 target genes — per-variant Evo 2 scores, Boltz-2 confidence metrics, functional state, evidence code and final interpretation |
| `cryptic_supervised_cv_predictions.csv` | 15.8 MB | Per-isolate, per-repeat cross-validation predictions of the supervised reference model |
| `cryptic_candidate_nested_cv_predictions.csv` | 13.0 MB | Nested cross-validation predictions for the candidate-priority analysis |
| `external_bdq_phenotype_validation_variant_calls.csv` | 8.3 MB | CRyPTIC per-isolate variant calls with AMR-Hunter interpretation |
| `final_external_bdq_phenotype_validation_variant_calls.csv` | 8.2 MB | Same table after the final rule-layer pass |
| `rule_layer_no_gate_variant_calls.csv` | 8.2 MB | Rule-layer-only ablation (no evolutionary gate) variant calls |

SQLite schema of BDQ-SatDB: `genes`, `mutations`, `drugs`, `species`, `boltz_cache`,
`synergy_events`, `synergy_event_mutations`. Key columns of `mutations`:

```
gene_id, pos, ref, alt, aa_change, region_type, region_name,
evo_model, evo_mt_score, evo_delta,
boltz_binding_affinity, boltz_stability_score, boltz_clash,
boltz_iptm, boltz_plddt, boltz_ptm, boltz_complex_energy,
status, functional_state, resistance_phenotype, final_interpretation, evidence_code
```

Reproducing the database from scratch requires the GPU inference services (Evo 2 7B and
Boltz-2) and roughly 1,700 GPU-hours for the structural stage; the archived database is the
practical route for re-analysis.

## Smaller tables

All remaining result tables (72 files, ~24 MB) are version-controlled inside the repository
under `reports/publication/`.
