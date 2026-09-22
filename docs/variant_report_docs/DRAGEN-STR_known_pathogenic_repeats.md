# DRAGEN-STR repeat genotypes for known pathogenic loci

Madeline Couse

**Version 2026-09**

## Changelog

Adapted from the [crg2-pacbio pipeline](https://github.com/ccmbioinfo/crg2-pacbio) for GRCh38 DRAGEN v4.4 genomes.
ExpansionHunter is run on each sample CRAM and is now the sole source of STR calls. Disease prediction uses current STRchive ranges and reports `BENIGN`, `INTERMEDIATE`, `PATHOGENIC`, `UNKNOWN`, or `MISSING` instead of `TRUE` or `FALSE`.

## Summary

The report uses the [STRchive v2.26.1 hg38 ExpansionHunter catalog](https://github.com/dashnowlab/STRchive/blob/v2.26.1/data/catalogs/STRchive-disease-loci.hg38.expansionhunter.json) and current STRchive regions, motifs, locus names, and thresholds. Four additional loci from the previous expanded catalog are retained, giving 86 reported loci.

For a multi-component locus, the report selects the pathogenic component using its hg38 coordinates and motif. Cyclic motif rotations and reverse complements are treated as equivalent when matching components. ExpansionHunter sex is taken from the peddy sex-check output.

The resulting CSV reports repeat sizes and classifies the called `motif_count` using STRchive benign, intermediate, and pathogenic ranges. Exact, contraction, and multiple-range rules are supported where they can be represented by repeat count. Descriptions for the report columns are listed below.

Suggestions for filtering and interpretation

  - Filter for pathogenic loci in the proband: &lt;proband_ID&gt;_DISEASE_PREDICTION == 'PATHOGENIC'

NB:

The pathogenicity status of some repeats might depend on sequence interruptions or motif changes that ExpansionHunter does not fully resolve.

Some VNTRs, including *MUC1* and *CEL*, cannot be interpreted reliably from repeat count alone and are reported as `UNKNOWN`.

At the *VWA1* locus, any deviation from two repeat copies is thought to be pathogenic, i.e. contractions or expansions are pathogenic (STRchive).

Updated STRchive regions and count definitions can change classification. In validation, 4 STRs were called `PATHOGENIC` across the tested family after previously being reported as `FALSE`. Pathogenic calls should therefore be taken with caution.

For detailed descriptions of repeat loci, including descriptions, pathogenic expansion repeat ranges, prevalence, age of onset, references and more, please refer to the [STRchive loci page](https://strchive.org/loci/).

## Column descriptions

| **Column** | **Comment** | **Source** | **Example** |
|---|---|---|---|
| CHROM | Chromosome | ExpansionHunter VCF | chrX |
| POS | Position| ExpansionHunter VCF | 147912049 |
| REF_REPEAT_COUNT | Number of repeat units spanned by the repeat in the reference | ExpansionHunter VCF | 20 |
| REF_LEN_BP | Reference length in bp | ExpansionHunter VCF | 60 |
| MOTIF | Reportable repeat motif | STRchive/ExpansionHunter VCF | CGG |
| GENE | Gene associated with repeat | STRchive | FMR1 |
| DISORDER | Disorder associated with repeat | [STRchive](https://strchive.org/) | FRAXA: fragile X syndrome |
| DISEASE_THRESHOLD | Primary STRchive disease threshold; special loci may use exact, contraction, or multiple-range rules | [STRchive](https://strchive.org/) | 201 |
| &lt;SAMPLE_ID&gt;_DISEASE_PREDICTION | Classification of the called motif count against the configured STRchive ranges | ExpansionHunter/STRchive | PATHOGENIC |
| &lt;SAMPLE_ID&gt;_GT| Genotype | ExpansionHunter VCF | 1/1 |
| &lt;SAMPLE_ID&gt;_motif_count | Motif count | ExpansionHunter VCF | 30/30 |
| &lt;SAMPLE_ID&gt;_REPCI | Repeat confidence interval | ExpansionHunter VCF | 22-30/30-39 |
| &lt;SAMPLE_ID&gt;_FILTER | ExpansionHunter record filter | ExpansionHunter VCF | PASS |
| &lt;SAMPLE_ID&gt;_SO | Read-support type for each allele | ExpansionHunter VCF | SPANNING/FLANKING |
| &lt;SAMPLE_ID&gt;_ADSP | Number of spanning reads consistent with allele | ExpansionHunter VCF | 4/4 |
| &lt;SAMPLE_ID&gt;_ADFL | Number of flanking reads consistent with allele | ExpansionHunter VCF | 21/21 |
| &lt;SAMPLE_ID&gt;_ADIR | Number of in-repeat reads consistent with allele | ExpansionHunter VCF | 0/0 |
| &lt;SAMPLE_ID&gt;_coverage | Locus coverage | ExpansionHunter VCF | 18.37 |
| STRCHIVE_URL | Link to the STRchive locus page, labelled with the target variant ID | [STRchive](https://strchive.org/) | `=HYPERLINK("https://strchive.org/loci/fxs_fmr1/","FXS_FMR1")` |
