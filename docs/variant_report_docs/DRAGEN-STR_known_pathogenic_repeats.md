# DRAGEN-STR repeat genotypes for known pathogenic loci

Madeline Couse

**Version 2026-09**

## Changelog

Adapted from the [crg2-pacbio pipeline](https://github.com/ccmbioinfo/crg2-pacbio) for GRCh38 DRAGEN v4.4 genomes.
ExpansionHunter is also run on each sample CRAM so that disease loci absent from the DRAGEN STR VCF can still be reported. Sample support fields are the shared DRAGEN/ExpansionHunter FORMAT fields (`motif_count`, `REPCI`, `ADSP`, `ADFL`, `ADIR`).

## Summary

The report uses the STRchive ExpansionHunter catalog (https://github.com/dashnowlab/STRchive/releases/download/v2.26.1/STRchive-disease-loci-v2.26.1.hg19.expansionhunter.json) and applies the following order for each locus:

1. Use the disease-labelled DRAGEN call when present.
2. Otherwise, use an unlabelled DRAGEN call when its interval and repeat motif exactly match the reportable STRchive component.
3. Otherwise, use the ExpansionHunter call.

For a multi-component ExpansionHunter locus, the reported component is selected using the hg38 coordinates and reference motif in the disease-threshold table. As in the PacBio report, its automatic disease prediction is left missing because one component count cannot be interpreted with a simple locus-wide threshold. ExpansionHunter sex is taken from the peddy sex-check output when available.

The resulting CSV reports repeat sizes and STRchive disease thresholds. Descriptions for the report columns are listed in the table below.

Suggestions for filtering and interpretation

  - Filter out loci that are not expanded in the proband: &lt;proband_ID&gt;_DISEASE_PREDICTION != 'FALSE'

NB:

The pathogenicity status of some repeats might depend on the presence of sequence interruptions or motif changes that DRAGEN-STR does not call.

At the *VWA1* locus, any deviation from two repeat copies is thought to be pathogenic, i.e. contractions or expansions are pathogenic (STRchive).

For detailed descriptions of repeat loci, including descriptions, pathogenic expansion repeat ranges, prevalence, age of onset, references and more, please refer to the [STRchive loci page](https://strchive.org/loci/).

## Column descriptions

| **Column** | **Comment** | **Source** | **Example** |
|---|---|---|---|
| CHROM | Chromosome | DRAGEN VCF with ExpansionHunter Fallback | chrX |
| POS | Position| DRAGEN VCF with ExpansionHunter Fallback | 147912050 |
| REF_REPEAT_COUNT | Number of repeat units spanned by the repeat in the reference | DRAGEN VCF with ExpansionHunter Fallback | 20 |
| REF_LEN_BP | Reference length in bp | DRAGEN VCF with ExpansionHunter Fallback | 60 |
| MOTIF | Repeat motif | DRAGEN VCF with ExpansionHunter Fallback | CGG |
| GENE | Gene associated with repeat | STRchive/DRAGEN mapping | FMR1 |
| DISORDER | Disorder associated with repeat | [STRchive](https://strchive.org/) | FRAXA: fragile X syndrome |
| DISEASE_THRESHOLD | Number of repeat units at/above which an expansion is considered disease-causing | [STRchive](https://strchive.org/) | 201 |
| &lt;SAMPLE_ID&gt;_DISEASE_PREDICTION | If sample motif count is greater than or equal to disease threshold, TRUE; otherwise FALSE | DRAGEN VCF with ExpansionHunter Fallback | FALSE |
| &lt;SAMPLE_ID&gt;_GT| Genotype | DRAGEN VCF with ExpansionHunter Fallback | 1/1 |
| &lt;SAMPLE_ID&gt;_motif_count | Motif count | DRAGEN VCF with ExpansionHunter Fallback | 30/30 |
| &lt;SAMPLE_ID&gt;_REPCI | Repeat confidence interval | DRAGEN VCF with ExpansionHunter Fallback | 22-30/30-39 |
| &lt;SAMPLE_ID&gt;_ADSP | Number of spanning reads consistent with allele | DRAGEN VCF with ExpansionHunter Fallback | 4/4 |
| &lt;SAMPLE_ID&gt;_ADFL | Number of flanking reads consistent with allele | DRAGEN VCF with ExpansionHunter Fallback | 21/21 |
| &lt;SAMPLE_ID&gt;_ADIR | Number of in-repeat reads consistent with allele | DRAGEN VCF with ExpansionHunter Fallback | 0/0 |
| &lt;SAMPLE_ID&gt;_coverage | Locus coverage | DRAGEN VCF with ExpansionHunter Fallback | 18.37 |
