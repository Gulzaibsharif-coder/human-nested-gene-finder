from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import urllib.request
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
GENES_BED = FRESH_DIR / "genes.bed"
EXONS_BED = FRESH_DIR / "exons.bed"
INTRONS_BED = FRESH_DIR / "introns.bed"
NESTED_BED = FRESH_DIR / "nested_matches.bed"
OUTPUT_FILE = FRESH_DIR / "nested_genes_fresh.csv"
METADATA_FILE = FRESH_DIR / "metadata.json"


CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def download_annotation() -> None:
    if GTF_FILE.exists():
        print(f"Using existing annotation: {GTF_FILE}")
        return

    print("Downloading GENCODE Release 50...")

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

        key, value = item.split(" ", 1)
        result[key] = value.strip('"')

    return result


def build_intervals() -> None:
    print("Parsing GENCODE and building genomic intervals...")

    genes = {}
    transcripts = {}

    with gzip.open(GTF_FILE, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chromosome = fields[0]

            if chromosome not in CHROMS:
                continue

            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")

            if feature == "gene" and gene_id:
                genes[gene_id] = {
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_name": attrs.get("gene_name", ""),
                    "gene_type": attrs.get(
                        "gene_type",
                        attrs.get("gene_biotype", ""),
                    ),
                }

            elif feature == "exon" and transcript_id and gene_id:
                transcripts.setdefault(transcript_id, []).append(
                    (
                        chromosome,
                        start,
                        end,
                        strand,
                        gene_id,
                    )
                )

    print(f"Genes parsed: {len(genes):,}")
    print(f"Transcripts with exons: {len(transcripts):,}")

    with GENES_BED.open("w", encoding="utf-8") as out:
        for gene_id, gene in genes.items():
            # BED is 0-based, half-open.
            bed_start = gene["start"] - 1
            bed_end = gene["end"]

            out.write(
                "\t".join(
                    [
                        gene["chromosome"],
                        str(bed_start),
                        str(bed_end),
                        gene_id,
                        gene["gene_name"],
                        gene["gene_type"],
                        gene["strand"],
                    ]
                )
                + "\n"
            )

    with EXONS_BED.open("w", encoding="utf-8") as out:
        for transcript_id, exons in transcripts.items():
            for chromosome, start, end, strand, gene_id in exons:
                bed_start = start - 1
                bed_end = end

                out.write(
                    "\t".join(
                        [
                            chromosome,
                            str(bed_start),
                            str(bed_end),
                            transcript_id,
                            gene_id,
                            strand,
                        ]
                    )
                    + "\n"
                )

    # Construct introns from consecutive exons within each transcript.
    with INTRONS_BED.open("w", encoding="utf-8") as out:
        intron_count = 0

        for transcript_id, exons in transcripts.items():
            if len(exons) < 2:
                continue

            exons.sort(key=lambda item: item[1])

            for left, right in zip(exons, exons[1:]):
                chromosome = left[0]
                left_end = left[2]
                right_start = right[1]

                intron_start = left_end + 1
                intron_end = right_start - 1

                if intron_start > intron_end:
                    continue

                # BED conversion.
                bed_start = intron_start - 1
                bed_end = intron_end

                out.write(
                    "\t".join(
                        [
                            chromosome,
                            str(bed_start),
                            str(bed_end),
                            f"{left[4]}",
                            transcript_id,
                            left[3],
                        ]
                    )
                    + "\n"
                )

                intron_count += 1

    print(f"Introns reconstructed: {intron_count:,}")


def run_bedtools() -> None:
    print("Running BEDTools containment analysis...")

    command = [
        "bedtools",
        "intersect",
        "-a",
        str(GENES_BED),
        "-b",
        str(INTRONS_BED),
        "-wa",
        "-wb",
        "-f",
        "1.0",
    ]

    with NESTED_BED.open("w", encoding="utf-8") as output:
        subprocess.run(
            command,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )

    print(f"BEDTools output: {NESTED_BED}")


def build_result() -> pd.DataFrame:
    print("Building final nested-gene table...")

    gene_columns = [
        "chromosome",
        "gene_start_0",
        "gene_end",
        "nested_gene_id",
        "nested_gene_name",
        "nested_gene_type",
        "nested_strand",
    ]

    intron_columns = [
        "chromosome_2",
        "intron_start_0",
        "intron_end",
        "host_gene_id",
        "transcript_id",
        "host_strand",
    ]

    df = pd.read_csv(
        NESTED_BED,
        sep="\t",
        header=None,
        names=gene_columns + intron_columns,
    )

    # Remove impossible self-nesting relationships.
    df = df[df["nested_gene_id"] != df["host_gene_id"]].copy()

    # Convert BED coordinates back to 1-based inclusive coordinates.
    df["nested_start"] = df["gene_start_0"] + 1
    df["nested_end"] = df["gene_end"]

    df["intron_start"] = df["intron_start_0"] + 1
    df["intron_end"] = df["intron_end"]

    result = df[
        [
            "host_gene_id",
            "nested_gene_id",
            "nested_gene_name",
            "nested_gene_type",
            "chromosome",
            "nested_start",
            "nested_end",
            "intron_start",
            "intron_end",
            "transcript_id",
        ]
    ].copy()

    result = result.drop_duplicates(
        subset=["host_gene_id", "nested_gene_id"]
    )

    # Final explicit validation.
    valid = (
        (result["intron_start"] <= result["nested_start"])
        & (result["intron_end"] >= result["nested_end"])
        & (result["host_gene_id"] != result["nested_gene_id"])
    )

    result = result.loc[valid].copy()

    result.to_csv(OUTPUT_FILE, index=False)

    print(f"Final relationships: {len(result):,}")
    print(f"Unique hosts: {result.host_gene_id.nunique():,}")
    print(f"Unique nested genes: {result.nested_gene_id.nunique():,}")

    return result


def write_metadata(result: pd.DataFrame) -> None:
    metadata = {
        "assembly": "GRCh38.p14",
        "annotation": "GENCODE Release 50",
        "ensembl_release": 116,
        "source_url": GTF_URL,
        "relationships": int(len(result)),
        "unique_host_genes": int(result.host_gene_id.nunique()),
        "unique_nested_genes": int(result.nested_gene_id.nunique()),
        "chromosomes": int(result.chromosome.nunique()),
        "tool": "BEDTools intersect",
        "containment_fraction": "1.0",
        "generated_by": "scripts/generate_fresh_list_bedtools.py",
    }

    METADATA_FILE.write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def main() -> None:
    print("=" * 64)
    print("HUMAN NESTED GENE FINDER — FRESH BEDTOOLS PIPELINE")
    print("=" * 64)

    download_annotation()
    build_intervals()
    run_bedtools()
    result = build_result()
    write_metadata(result)

    print("=" * 64)
    print("FRESH ANALYSIS COMPLETE")
    print("=" * 64)


if __name__ == "__main__":
    main()
