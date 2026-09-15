import sqlite3
from scripts.refresh_gene_role_migration import refresh_database

MINIMAL_SCHEMA = "\nCREATE TABLE genes (\n    id INTEGER PRIMARY KEY AUTOINCREMENT,\n    species_id INTEGER NOT NULL DEFAULT 1,\n    name TEXT NOT NULL,\n    logic_type TEXT DEFAULT 'POSITIVE',\n    scenario TEXT,\n    pathway_group TEXT,\n    sequence TEXT,\n    UNIQUE(species_id, name)\n);\nCREATE TABLE mutations (\n    id INTEGER PRIMARY KEY AUTOINCREMENT,\n    gene_id INTEGER NOT NULL,\n    pos INTEGER,\n    ref TEXT,\n    alt TEXT,\n    aa_change TEXT,\n    region_type TEXT NOT NULL DEFAULT 'CDS',\n    status TEXT NOT NULL DEFAULT 'PENDING',\n    evo_delta REAL,\n    resistance_phenotype TEXT,\n    functional_state TEXT,\n    evidence_code TEXT,\n    final_interpretation TEXT,\n    synergy_tags TEXT,\n    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP\n);\n"


def _build_historical_db(tmp_path) -> sqlite3.Connection:
    conn = sqlite3.connect(tmp_path / "test.db")
    conn.executescript(MINIMAL_SCHEMA)
    conn.execute("INSERT INTO genes (name, logic_type) VALUES ('mtrA', 'NEGATIVE')")
    conn.execute("INSERT INTO genes (name, logic_type) VALUES ('mtrB', 'NEGATIVE')")
    conn.execute("INSERT INTO genes (name, logic_type) VALUES ('Rv1979c', 'NEGATIVE')")
    conn.execute("INSERT INTO genes (name, logic_type) VALUES ('glpK', 'NEGATIVE')")
    conn.execute(
        "INSERT INTO mutations (gene_id, aa_change, status, resistance_phenotype, functional_state, evidence_code) VALUES (1, '229*', 'LOSS_OF_FUNCTION', 'RESISTANT', 'LOSS_OF_FUNCTION', 'NEGATIVE_LOF_SIGNATURE')"
    )
    conn.execute(
        "INSERT INTO mutations (gene_id, aa_change, status, resistance_phenotype, functional_state, evidence_code) VALUES (2, '310*', 'LOSS_OF_FUNCTION', 'RESISTANT', 'LOSS_OF_FUNCTION', 'NEGATIVE_LOF_SIGNATURE')"
    )
    conn.execute(
        "INSERT INTO mutations (gene_id, aa_change, status, resistance_phenotype, functional_state, evidence_code) VALUES (3, '100*', 'LOSS_OF_FUNCTION', 'RESISTANT', 'LOSS_OF_FUNCTION', 'NEGATIVE_LOF_SIGNATURE')"
    )
    conn.execute(
        "INSERT INTO mutations (gene_id, aa_change, status, resistance_phenotype, functional_state, evidence_code) VALUES (4, '200*', 'LOSS_OF_FUNCTION', 'RESISTANT', 'LOSS_OF_FUNCTION', 'NEGATIVE_LOF_SIGNATURE')"
    )
    conn.commit()
    return conn


def test_refresh_database_migrates_roles_and_decisions(tmp_path):
    conn = _build_historical_db(tmp_path)
    conn.close()
    genes_updated, mutations_updated = refresh_database(tmp_path / "test.db")
    assert genes_updated == 2
    assert mutations_updated == 4
    conn = sqlite3.connect(tmp_path / "test.db")
    conn.row_factory = sqlite3.Row
    logic_by_gene = dict(conn.execute("SELECT name, logic_type FROM genes").fetchall())
    assert logic_by_gene["mtrA"] == "STRUCTURAL"
    assert logic_by_gene["mtrB"] == "STRUCTURAL"
    rows = {
        row["name"]: (
            row["resistance_phenotype"],
            row["functional_state"],
            row["evidence_code"],
        )
        for row in conn.execute(
            "\n            SELECT g.name, m.resistance_phenotype, m.functional_state, m.evidence_code\n            FROM mutations m JOIN genes g ON m.gene_id = g.id\n            "
        ).fetchall()
    }
    assert rows["mtrA"] == ("SENSITIVE", "LOSS_OF_FUNCTION", "STRUCTURAL_LOF_SIGNATURE")
    assert rows["mtrB"] == ("SENSITIVE", "LOSS_OF_FUNCTION", "STRUCTURAL_LOF_SIGNATURE")
    assert rows["Rv1979c"] == ("UNKNOWN", "LOSS_OF_FUNCTION", "EXPLORATORY_ONLY")
    assert rows["glpK"] == ("UNKNOWN", "LOSS_OF_FUNCTION", "EXPLORATORY_ONLY")
    conn.close()
