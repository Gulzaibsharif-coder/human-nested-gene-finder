import os
import pandas as pd

from parse_gtf import read_gtf


def build_introns(exons):

    introns = []

    for transcript_id, group in exons.groupby("transcript_id"):

        group = group.sort_values("Start")

        exon_rows = list(group.itertuples(index=False))

        if len(exon_rows) < 2:
            continue

        for i in range(len(exon_rows) - 1):

            current_exon = exon_rows[i]
            next_exon = exon_rows[i + 1]

            intron_start = current_exon.End + 1
            intron_end = next_exon.Start - 1

            if intron_start > intron_end:
                continue

            introns.append({
                "Chromosome": current_exon.Chromosome,
                "Start": intron_start,
                "End": intron_end,
                "Strand": current_exon.Strand,
                "gene_id": current_exon.gene_id,
                "transcript_id": transcript_id
            })

    return pd.DataFrame(introns)


if __name__ == "__main__":

    genes, transcripts, exons = read_gtf()

    introns = build_introns(exons)

    os.makedirs("results", exist_ok=True)

    introns.to_csv(
        "results/introns.csv",
        index=False
    )

    print("Genes:", len(genes))
    print("Transcripts:", len(transcripts))
    print("Exons:", len(exons))
    print("Introns:", len(introns))
