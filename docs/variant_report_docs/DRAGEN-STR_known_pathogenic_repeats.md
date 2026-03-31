# DRAGEN-STR repeat genotypes for known pathogenic loci

Madeline Couse

**Version 2026-03**

## Changelog

Adapted from the [crg2-pacbio pipeline](https://github.com/ccmbioinfo/crg2-pacbio) for GRCh38 DRAGEN v4.4 genomes.
Sample support fields are DRAGEN FORMAT fields (`motif_count`, `REPCI`, `ADSP`, `ADFL`, `ADIR`) rather than PacBio TRGT motif sequence, motif structure, and methylation fields.

## Summary

DRAGEN-STR (based on ExpansionHunter) genotypes disease-causing repeats located in AFF2, AR, ARX_1, ARX_2, ATN1, ATXN1, ATXN10, ATXN2, ATXN3, ATXN7, ATXN8OS, BEAN1, C9ORF72, CACNA1A, CBL, CNBP, COMP, CSTB, DAB1, DIP2B, DMD, DMPK, EIF4A3, FMR1, FOXL2, FXN, GIPC1, GLS, HOXA13_1, HOXA13_2, HOXA13_3, HOXD13, HTT, JPH3, LRP12, MARCHF6, NIPA1, NOP56, NOTCH2NLC, NUTM2B-AS1, PABPN1, PHOX2B, PPP2R2B, PRDM12, PRNP, RAPGEF2, RFC1, RUNX2, SAMD12, SOX3, STARD7, TBP, TBX1, TCF4, TNRC6A, VWA1, XYLT1, YEATS2, ZIC2 and ZIC3 genes. We create a CSV report detailing repeat sizes with disease thresholds from the [STRchive BED file](https://strchive.org/_astro/STRchive-disease-loci.hg38.CE2vK2zA.bed). Descriptions for the report columns are listed in the table below.

Suggestions for filtering and interpretation

  - Filter out loci that are not expanded in the proband: &lt;proband_ID&gt;_DISEASE_PREDICTION != 'FALSE'

NB:

The pathogenicity status of some repeats might depend on the presence of sequence interruptions or motif changes that DRAGEN-STR does not call.

At the *VWA1* locus, any deviation from two repeat copies is thought to be pathogenic, i.e. contractions or expansions are pathogenic (STRchive).

For detailed descriptions of repeat loci, including descriptions, pathogenic expansion repeat ranges, prevalence, age of onset, references and more, please refer to the [STRchive loci page](https://strchive.org/loci/).

## Column descriptions

| **Column** | **Comment** | **Source** | **Example** |
|---|---|---|---|
| CHROM | Chromosome | DRAGEN VCF | chrX |
| POS | Position| DRAGEN VCF | 147912050 |
| REF_REPEAT_COUNT | Number of repeat units spanned by the repeat in the reference | DRAGEN VCF | 20 |
| REF_LEN_BP | Reference length in bp | DRAGEN VCF | 60 | 
| MOTIF | Repeat motif | DRAGEN VCF | CGG | 
| GENE | Gene associated with repeat | DRAGEN VCF | FMR1 | 
| DISORDER | Disorder associated with repeat | [STRchive](https://strchive.org/) | FRAXA: fragile X syndrome |
| DISEASE_THRESHOLD | Number of repeat units at/above which an expansion is considered disease-causing | [STRchive](https://strchive.org/) | 201 |
| &lt;SAMPLE_ID&gt;_DISEASE_PREDICTION | If sample motif count is greater than or equal to disease threshold, TRUE; otherwise FALSE | DRAGEN VCF | FALSE |
| &lt;SAMPLE_ID&gt;_GT| Genotype | DRAGEN VCF | 1/1 |
| &lt;SAMPLE_ID&gt;_motif_count | Motif count | DRAGEN VCF | 30/30 |
| &lt;SAMPLE_ID&gt;_REPCI | Repeat confidence interval | DRAGEN VCF | 22-30/30-39 |
| &lt;SAMPLE_ID&gt;_ADSP | Number of spanning reads consistent with allele | DRAGEN VCF | 4/4 |
| &lt;SAMPLE_ID&gt;_ADFL | Number of flanking reads consistent with allele | DRAGEN VCF | 21/21 |
| &lt;SAMPLE_ID&gt;_ADIR | Number of in-repeat reads consistent with allele | DRAGEN VCF | 0/0 |
| &lt;SAMPLE_ID&gt;_coverage | Locus coverage | DRAGEN VCF | 18.37 |
