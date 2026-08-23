from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FRESH_DIR = ROOT / "fresh_data"
FRESH_DIR.mkdir(exist_ok=True)

GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/"
    "Gencode_human/release_50/gencode.v50.annotation.gtf.gz"
)

GTF_FILE = FRESH_DIR / "gencode.v50.annotation.gtf.gz"
GENES_FILE = FRESH_DIR / "genes.csv"
INTRONS_FILE = FRESH_DIR / "introns.csv"
OUTPUT_FILE = FRESH_DIR / "nested_genes_fresh.csv"
METADATA_FILE = FRESH_DIR / "metadata.json"


def download_annotation() -> None:
    if GTF_FILE.exists():
        print(f"Using existing annotation: {GTF_FILE}")
        return

    print("Downloading current GENCODE Release 50 annotation...")
    request = urllib.request.Request(
        GTF_URL,
        headers={"User-Agent": "human-nested-gene-finder/1.0"},
    )

    with urllib.request.urlopen(request, timeout=180) as response:
        with GTF_FILE.open("wb") as output:
            shutil.copyfileobj(response, output)

    print("Download complete.")


def parse_attributes(text: str) -> dict[str, str]:
    result = {}

    for item in text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue

        parts = item.split(" ", 1)
        if len(parts) != 2:
            continue

        key, value = parts
        result[key] = value.strip('"')

    return result


def parse_gtf() -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = []
    exons_by_transcript = defaultdict(list)

    with gzip.open(GTF_FILE, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chromosome = fields[0]
            feature = fields[2]

            if chromosome not in {
                *[f"chr{i}" for i in range(1, 23)],
                "chrX",
                "chrY",
            }:
                continue

            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            gene_id = attrs.get("gene_id")
            gene_name = attrs.get("gene_name", "")
            gene_type = attrs.get(
                "gene_type",
                attrs.get("gene_biotype", "")
            )
            transcript_id = attrs.get("transcript_id")

            if feature == "gene" and gene_id:
                genes.append({
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                })

            elif feature == "exon" and transcript_id:
                exons_by_transcript[transcript_id].append(
                    {
                        "chromosome": chromosome,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "gene_id": gene_id,
                        "transcript_id": transcript_id,
                    }
                )

    genes_df = pd.DataFrame(genes)

    # Reconstruct transcript introns.
    introns = []

    for transcript_id, exons in exons_by_transcript.items():
        if len(exons) < 2:
            continue

        exons.sort(key=lambda x: x["start"])

        for left, right in zip(exons, exons[1:]):
            intron_start = left["end"] + 1
            intron_end = right["start"] - 1

            if intron_start > intron_end:
                continue

            introns.append(
                {
                    "chromosome": left["chromosome"],
                    "start": intron_start,
                    "end": intron_end,
                    "strand": left["strand"],
                    "gene_id": left["gene_id"],
                    "transcript_id": transcript_id,
                }
            )

    introns_df = pd.DataFrame(introns)

    genes_df.to_csv(GENES_FILE, index=False)
    introns_df.to_csv(INTRONS_FILE, index=False)

    print(f"Genes: {len(genes_df):,}")
    print(f"Transcripts with introns: {len(exons_by_transcript):,}")
    print(f"Introns: {len(introns_df):,}")

    return genes_df, introns_df


def find_nested_genes(
    genes_df: pd.DataFrame,
    introns_df: pd.DataFrame,
) -> pd.DataFrame:

    nested = []

    for chromosome, gene_group in genes_df.groupby("chromosome"):
        intron_group = introns_df[
            introns_df["chromosome"] == chromosome
        ]

        if intron_group.empty:
            continue

        for gene in gene_group.itertuples(index=False):
            matches = intron_group[
                (intron_group["start"] <= gene.start)
                & (intron_group["end"] >= gene.end)
                & (intron_group["gene_id"] != gene.gene_id)
            ]

            for intron in matches.itertuples(index=False):
                nested.append(
                    {
                        "host_gene_id": intron.gene_id,
                        "nested_gene_id": gene.gene_id,
                        "nested_gene_name": gene.gene_name,
                        "nested_gene_type": gene.gene_type,
                        "chromosome": chromosome,
                        "nested_start": gene.start,
                        "nested_end": gene.end,
                        "intron_start": intron.start,
                        "intron_end": intron.end,
                        "intron_transcript_id": intron.transcript_id,
                    }
                )

    result = pd.DataFrame(nested)

    if result.empty:
        result.to_csv(OUTPUT_FILE, index=False)
        return result

    result = result.drop_duplicates(
        subset=["host_gene_id", "nested_gene_id"]
    )

    # Final strict validation.
    valid = (
        (result["intron_start"] <= result["nested_start"])
        & (result["intron_end"] >= result["nested_end"])
        & (result["host_gene_id"] != result["nested_gene_id"])
    )

    result = result.loc[valid].copy()

    result.to_csv(OUTPUT_FILE, index=False)

    print(f"Fresh nested relationships: {len(result):,}")
    print(f"Unique host genes: {result['host_gene_id'].nunique():,}")
    print(f"Unique nested genes: {result['nested_gene_id'].nunique():,}")

    return result


def write_metadata(result: pd.DataFrame) -> None:
    metadata = {
        "assembly": "GRCh38.p14",
        "annotation": "GENCODE Release 50",
        "ensembl_release": 116,
        "source_url": GTF_URL,
        "relationships": int(len(result)),
        "unique_host_genes": (
            int(result["host_gene_id"].nunique())
            if not result.empty else 0
        ),
        "unique_nested_genes": (
            int(result["nested_gene_id"].nunique())
            if not result.empty else 0
        ),
        "coordinate_rule": (
            "intron_start <= nested_start and "
            "intron_end >= nested_end"
        ),
    }

    METADATA_FILE.write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    print(f"Metadata: {METADATA_FILE}")


def main() -> None:
    print("=" * 60)
    print("HUMAN NESTED GENE FINDER — FRESH LIST")
    print("=" * 60)

    download_annotation()

    genes_df, introns_df = parse_gtf()

    result = find_nested_genes(genes_df, introns_df)

    write_metadata(result)

    print("=" * 60)
    print("FRESH LIST COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
