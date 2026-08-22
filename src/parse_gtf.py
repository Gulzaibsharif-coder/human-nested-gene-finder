import gzip
import re
import pandas as pd

GTF_FILE = "data/gencode.v50.annotation.gtf.gz"


def parse_attributes(attribute_string):
    attributes = {}

    for match in re.finditer(
        r'(\S+)\s+"([^"]*)"',
        attribute_string
    ):
        key = match.group(1)
        value = match.group(2)
        attributes[key] = value

    return attributes


def read_gtf():

    genes = []
    transcripts = []
    exons = []

    with gzip.open(GTF_FILE, "rt") as file:

        for line in file:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chromosome = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            attributes = parse_attributes(fields[8])

            gene_id = attributes.get("gene_id")
            gene_name = attributes.get("gene_name")
            gene_type = attributes.get("gene_type")
            transcript_id = attributes.get("transcript_id")

            if feature == "gene":

                genes.append({
                    "Chromosome": chromosome,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type
                })

            elif feature == "transcript":

                transcripts.append({
                    "Chromosome": chromosome,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "gene_id": gene_id,
                    "transcript_id": transcript_id
                })

            elif feature == "exon":

                exons.append({
                    "Chromosome": chromosome,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "gene_id": gene_id,
                    "transcript_id": transcript_id
                })

    return (
        pd.DataFrame(genes),
        pd.DataFrame(transcripts),
        pd.DataFrame(exons)
    )


if __name__ == "__main__":

    genes, transcripts, exons = read_gtf()

    print("Genes:", len(genes))
    print("Transcripts:", len(transcripts))
    print("Exons:", len(exons))
