from __future__ import annotations
import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.cryptic_bdq_validation_common import (
    DEFAULT_CRYPTIC_METADATA,
    DEFAULT_METADATA_INSPECTION_OUTPUT,
    REQUIRED_METADATA_FIELDS,
    read_csv_rows,
    require_fields,
    summarize_metadata_rows,
    write_metadata_summary,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect CRyPTIC reuse metadata for BDQ phenotype availability."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_CRYPTIC_METADATA,
        help="CRyPTIC reuse metadata CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_METADATA_INSPECTION_OUTPUT,
        help="Output CSV for metadata availability summary.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        rows, fieldnames = read_csv_rows(args.input)
        require_fields(fieldnames, REQUIRED_METADATA_FIELDS)
        summary_rows = summarize_metadata_rows(rows)
        write_metadata_summary(args.output, summary_rows)
    except FileNotFoundError:
        print(f"CRyPTIC metadata file not found: {args.input}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    feasibility = next(
        (row for row in summary_rows if row["metric"] == "validation_feasibility"), None
    )
    print(f"metadata_summary={args.output}")
    if feasibility is not None:
        print(f"validation_feasibility={feasibility['value']}")
        if feasibility.get("notes"):
            print(feasibility["notes"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
