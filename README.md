# AMR-Hunter

**A mechanism-driven framework for prioritizing bedaquiline resistance-associated variants in
the *Mycobacterium tuberculosis* complex.**

AMR-Hunter interprets bedaquiline (BDQ) variants without phenotype-labelled model fitting
(*zero-shot*). It translates established BDQ resistance biology into explicit, auditable
interpretation rules and combines them with sequence-level and structure-level evidence:

| Level | What happens | Engine |
|---|---|---|
| 1. Evolutionary-fitness triage | mutant vs. wild-type sequence likelihood, ΔEvo | Evo 2 (7B) |
| 2. Structural re-examination | pTM / ipTM / pLDDT, stability and steric-clash scores | Boltz-2 |
| 3. Pathway inference | regulator–effector epistasis and multi-gene synergy rules | rule layer |
| 4. Integrated decision | R / S / indeterminate per variant, aggregated per isolate | rule layer |

Reference: *AMR-Hunter: A Mechanism-Driven Framework for Prioritizing Bedaquiline
Resistance-Associated Variants in the Mycobacterium tuberculosis Complex* — see
`CITATION.cff`.

## Target genes

The framework operates on ten BDQ resistance-associated genes, grouped by the direction in
which loss of function (LoF) shifts BDQ susceptibility:

| Mechanistic class | Genes | Interpretation of LoF |
|---|---|---|
| Target modification (POSITIVE) | `atpE` | resistance only for binding-interface substitutions; LoF → susceptible |
| LoF-mediated resistance (NEGATIVE) | `Rv0678`, `pepQ` | resistant |
| Functional-integrity dependence (STRUCTURAL) | `mmpS5`, `mmpL5` | susceptible (collapsed efflux) |
| LoF-mediated sensitization (SENSITIZATION) | `mtrA`, `mtrB`, `lpqB` | susceptible (hypersensitive) |
| Exploratory (EXPLORATORY) | `Rv1979c`, `glpK` | recorded, no definitive call |

The complete rule set, including the Rv0678 de-repression / efflux-collapse re-sensitization
rule, is defined in `config/config.yaml`.

## Repository layout

```
main.py                 command-line entry point (single-gene / batch runs)
run_pipeline.py         staged pipeline runner (init -> evo -> boltz -> analyze)
auto_run.py             unattended batch driver
dashboard_server.py     local web console
core/                   rule engine, consequence logic, DB, generators, model clients
services/               Evo 2 and Boltz-2 inference servers (GPU)
config/config.yaml      target genes, regulatory windows, thresholds, scenarios
scripts/                data preparation, benchmark, export and figure scripts
tests/                  unit tests for consequence logic and rule engine
web/                    dashboard front-end assets
```

## Installation

Requires Python ≥ 3.11.

```bash
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e .
```

The heavy inference engines are **not** installed by `pyproject.toml`; they run as separate
services on a GPU host:

```bash
pip install evo2          # Arc Institute Evo 2
pip install boltz         # Boltz-2
```

`services/evo2_server.py` and `services/boltz2_server.py` expose them over HTTP; ports and
model names are configured in `config/config.yaml`.

## Required inputs

Inputs are fetched automatically when `data_autofetch.enabled: true` in `config/config.yaml`:

| Asset | Source |
|---|---|
| H37Rv reference genome (NC_000962.3) | NCBI Entrez |
| wild-type structures for each target gene | RCSB PDB (7D00, 4NB5, 2GWR) |
| bedaquiline ligand SDF | PubChem |

Paths are resolved relative to `config/config.yaml`.

## Running the pipeline

```bash
python scripts/init_project.py            # 1) build the variant library from the annotation
python run_pipeline.py --stage evo        # 2) sequence-level triage
python run_pipeline.py --stage boltz      # 3) structural re-examination of BOLTZ_READY candidates
python run_pipeline.py --stage analyze    # 4) rule layer: epistasis + synergy + exports
```

Both inference servers can be started with `docker compose up`.

## Reproducing the manuscript

`REPRODUCE.md` maps every table, figure and headline number in the paper to the script that
generates it. Intermediate result tables are in `reports/publication/`; the
saturation-mutagenesis database (BDQ-SatDB) and the largest prediction tables are attached
to the GitHub release — see `RELEASE_ASSETS.md`.

## Data availability

- BDQ-SatDB: whole-site saturation-mutagenesis database for the ten target genes.
- CRyPTIC per-isolate and per-variant call tables and the supervised-model prediction tables.
- The Mtb H37Rv reference genome (NC_000962.3) and the WHO catalogue of mutations are public
  resources and are not redistributed here. Third-party published material is cited, not
  redistributed.

## License

Apache License 2.0 — see `LICENSE` and `NOTICE`.
