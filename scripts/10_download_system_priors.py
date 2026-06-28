from __future__ import annotations

import argparse
from pathlib import Path

import requests


ENRICHR_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary"

LIBRARIES = {
    "reactome": "Reactome_Pathways_2024",
    "hallmark": "MSigDB_Hallmark_2020",
    "tf_target": "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
    "ppi_hub": "PPI_Hub_Proteins",
}


def download_library(library_name: str, output: Path) -> None:
    if output.exists() and output.stat().st_size > 0:
        return
    response = requests.get(
        ENRICHR_URL,
        params={"mode": "text", "libraryName": library_name},
        timeout=90,
    )
    response.raise_for_status()
    output.write_text(response.text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="data/priors")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for short_name, library_name in LIBRARIES.items():
        output = out_dir / f"{short_name}.gmt"
        download_library(library_name, output)
        print(f"{short_name}: {library_name} -> {output}")


if __name__ == "__main__":
    main()
