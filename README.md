# Human Nested Gene Finder

A reproducible bioinformatics pipeline for identifying human genes whose genomic coordinates are fully contained within intronic regions of other annotated genes.

## Reference

- Genome assembly: GRCh38.p14
- GENCODE annotation: Release 49
- Annotation region: human reference chromosomes (CHR)
- Analysis definition: a nested gene is a gene whose complete genomic interval is contained within an intron reconstructed from an annotated transcript.

## Pipeline

1. Parse the GENCODE GTF annotation.
2. Extract genes, transcripts, and exons.
3. Reconstruct transcript-level introns from consecutive exon boundaries.
4. Compare gene intervals with intron intervals on the same chromosome.
5. Retain genes completely contained within another transcript's intron.
6. Validate coordinate containment.
7. Remove duplicate host–nested gene pairs.
8. Characterize nested-gene architecture.

## Results

The analysis produced:

- 78,691 genes
- 507,365 transcripts
- 3,673,949 exons
- 3,166,584 reconstructed introns
- 31,300 validated host–nested gene relationships
- 14,324 unique host genes
- 26,894 unique nested genes
- 31,300/31,300 relationships passed coordinate validation
- 0 invalid coordinate relationships
- 0 exact duplicate rows
- 0 duplicate host–nested gene pairs

The median host contained one nested gene, while the maximum observed host contained 89 distinct nested genes.

Among nested relationships, lncRNA was the most frequent nested-gene category, followed by processed pseudogenes.

## Reproducibility

Install dependencies:

```bash
pip install -r requirements.txt