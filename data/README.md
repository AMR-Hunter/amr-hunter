# Data inputs

Nothing in this directory is version-controlled: the required inputs are public and
downloaded automatically.

Run:

```bash
python scripts/init_project.py
```

with `data_autofetch.enabled: true` in `config/config.yaml`. The script checks for each file
listed below and downloads anything missing.

| Path | Asset | Source |
|---|---|---|
| `reference/mtb_h37rv.gbk` | *M. tuberculosis* H37Rv reference genome, NC_000962.3 | NCBI Entrez |
| `structures/atpE_wt.pdb`, `atpB`-free panel | wild-type target structures | RCSB PDB 7D00 (atpE), 4NB5 (Rv0678, mmpS5, mmpL5, Rv1979c, pepQ, glpK, lpqB), 2GWR (mtrA, mtrB) |
| `ligands/bedaquiline.sdf` | bedaquiline ligand | PubChem |

The reference genome supplies CDS coordinates (gene / locus_tag, start, end, strand) used to
enumerate the saturation library and to define the regulatory windows listed in
`config/config.yaml`.

## External validation data (not redistributed)

The CRyPTIC reuse cohort used for external validation is a public resource with its own
accession and data-sharing terms; obtain it from the CRyPTIC consortium release and place the
per-isolate VCFs under `data/external/cryptic/`. The derived per-isolate and per-variant call
tables produced from it are archived in `reports/publication/` and in the release assets.

The WHO catalogue of mutations in *M. tuberculosis* complex is distributed by the WHO and is
not redistributed here.
