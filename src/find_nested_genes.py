import pandas as pd

GENES_FILE = "results/genes.csv"
INTRONS_FILE = "results/introns.csv"
OUTPUT_FILE = "results/nested_genes.csv"

genes = pd.read_csv(GENES_FILE)
introns = pd.read_csv(INTRONS_FILE)

print("Genes:", len(genes))
print("Introns:", len(introns))

nested = []

for chromosome, gene_group in genes.groupby("Chromosome"):
    intron_group = introns[introns["Chromosome"] == chromosome]

    for _, gene in gene_group.iterrows():
        matches = intron_group[
            (intron_group["Start"] <= gene["Start"]) &
            (intron_group["End"] >= gene["End"]) &
            (intron_group["gene_id"] != gene["gene_id"])
        ]

        for _, intron in matches.iterrows():
            nested.append({
                "host_gene_id": intron["gene_id"],
                "nested_gene_id": gene["gene_id"],
                "nested_gene_name": gene["gene_name"],
                "chromosome": chromosome,
                "nested_start": gene["Start"],
                "nested_end": gene["End"],
                "intron_start": intron["Start"],
                "intron_end": intron["End"]
            })

result = pd.DataFrame(nested)

result.to_csv(OUTPUT_FILE, index=False)

print("Nested genes:", len(result))
print("Saved:", OUTPUT_FILE)
