# DRAGEN-STR repeat genotypes for known pathogenic loci

Madeline Couse

**Version 2026-09**

## Changelog

Adapted from the [crg2-pacbio pipeline](https://github.com/ccmbioinfo/crg2-pacbio) for GRCh38 DRAGEN pipeline.

The STR workflow has been updated to run [ExpansionHunter v5.0.0](https://github.com/Illumina/ExpansionHunter/releases/tag/v5.0.0), instead of using DRAGEN repeat VCFs, on each sample CRAM. The report now uses repeat definitions and ranges and reports `BENIGN`, `INTERMEDIATE`, `PATHOGENIC`, `UNKNOWN`, or `MISSING` instead of a threshold-only `TRUE` or `FALSE` result.

The report also includes the ExpansionHunter filter, allele support type, classification ranges, locus-specific interpretation notes, and links to the corresponding STRchive loci.

## Summary

ExpansionHunter genotypes disease-associated repeats at 86 configured loci. The CSV report details the associated gene and disorder, repeat size, disease prediction, confidence interval, quality information, and read support for every sample in the family. Descriptions of the report columns are listed in the table below.

The hg38 catalog and locus information use [STRchive v2.26.1](https://github.com/dashnowlab/STRchive/releases/tag/v2.26.1) as the base. The [STRchive ExpansionHunter catalog](https://github.com/dashnowlab/STRchive/releases/download/v2.26.1/STRchive-disease-loci-v2.26.1.hg38.expansionhunter.json) is supplemented at selected loci with definitions from a pinned version of the Broad Institute [`str-analysis` hg38 catalog](https://github.com/broadinstitute/str-analysis/blob/6490e03a81187795ed37f452ad5e64ca36e4b53a/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json). Four legacy supplemental loci absent from current STRchive are also retained.

Broad definitions are used where they better represent the hg38 locus. The main reasons are that the STRchive interval includes extra bases that can shift the repeat count, the motif phase does not match the start of the hg38 interval, or a compound or interrupted repeat is represented more completely by the Broad definition. STRchive remains the primary source for locus and disease information. A Broad definition changes how the repeat is represented and counted but does not by itself change its interpretation.

At loci with more than one repeat motif, the report shows the component relevant to the associated disorder. The `NOTE` column identifies loci where the displayed count has an unusual convention, the locus has a complex repeat structure, or a source difference or interpretation caveat requires additional context.

The `<SAMPLE_ID>_DISEASE_PREDICTION` columns report `PATHOGENIC` when a called allele is in a pathogenic range, `INTERMEDIATE` for an intermediate, reduced-penetrance, or premutation range, and `BENIGN` when the called alleles are in the benign range. `UNKNOWN` indicates that the call cannot be interpreted confidently from repeat count alone or does not fall within a configured range. `MISSING` indicates that no usable repeat-count call is available.

Suggestions for filtering and interpretation

  - Filter for pathogenic loci in the proband: `<proband_ID>_DISEASE_PREDICTION == 'PATHOGENIC'`.
  - To retain all results that may require review: `<proband_ID>_DISEASE_PREDICTION != 'BENIGN'`.
  - Consider `REPCI`, `FILTER`, `SO`, and the spanning, flanking, and in-repeat read counts together with the disease prediction.

NB:

The pathogenicity status of some repeats might depend on sequence interruptions or motif changes that ExpansionHunter does not fully resolve from repeat count alone.

Some VNTRs, including *MUC1* and *CEL*, cannot be interpreted reliably from repeat count alone and are reported as `UNKNOWN`.

At the *VWA1* locus, any deviation from two repeat copies is considered pathogenic; this includes both contractions and expansions.

The disease prediction is a screening aid and should not be treated as a standalone interpretation.

For detailed descriptions of repeat loci, including pathogenic repeat ranges, prevalence, age of onset, supporting evidence, and references, refer to the [STRchive loci page](https://strchive.org/loci/).

## Column descriptions

| **Column** | **Comment** | **Source** | **Example** |
|---|---|---|---|
| CHROM | Chromosome containing the repeat | STRchive/Broad catalog | chrX |
| POS | Start position of the reported repeat component | STRchive/Broad catalog | 147912050 |
| REF_REPEAT_COUNT | Number of repeat units in the hg38 reference allele | ExpansionHunter VCF (`REF`) | 20 |
| REF_LEN_BP | Length of the hg38 reference allele in base pairs | ExpansionHunter VCF (`RL`) | 60 |
| MOTIF | Repeat motif counted in the report | STRchive/Broad catalog | CGG |
| GENE | Gene associated with the repeat | STRchive | FMR1 |
| DISORDER | Disorder associated with the repeat | STRchive | FRAXA: fragile X syndrome |
| DISEASE_THRESHOLD | Primary disease-associated threshold shown for quick reference; consult the range columns and `NOTE` for interpretation | STRchive/Broad thresholds | 201 |
| NOTE | Locus-specific context that may affect interpretation; `.` when no additional note is required | STRchive/Broad comparison | Allele structure, read support, and the literature should be looked at carefully. |
| &lt;SAMPLE_ID&gt;_DISEASE_PREDICTION | Screening category for the sample call | ExpansionHunter VCF and STRchive/Broad thresholds | PATHOGENIC |
| &lt;SAMPLE_ID&gt;_GT | Genotype | ExpansionHunter VCF | 1/1 |
| &lt;SAMPLE_ID&gt;_motif_count | Estimated repeat count for each allele | ExpansionHunter VCF (`REPCN`) | 30/30 |
| &lt;SAMPLE_ID&gt;_REPCI | Confidence interval around each allele's repeat count | ExpansionHunter VCF | 22-30/30-39 |
| &lt;SAMPLE_ID&gt;_FILTER | `PASS` or the ExpansionHunter quality flag for the call | ExpansionHunter VCF | PASS |
| &lt;SAMPLE_ID&gt;_SO | Type of read evidence supporting each allele | ExpansionHunter VCF | SPANNING/FLANKING |
| &lt;SAMPLE_ID&gt;_ADSP | Number of spanning reads consistent with each allele | ExpansionHunter VCF | 4/4 |
| &lt;SAMPLE_ID&gt;_ADFL | Number of flanking reads consistent with each allele | ExpansionHunter VCF | 21/21 |
| &lt;SAMPLE_ID&gt;_ADIR | Number of in-repeat reads consistent with each allele | ExpansionHunter VCF | 0/0 |
| &lt;SAMPLE_ID&gt;_coverage | Locus coverage | ExpansionHunter VCF (`LC`) | 18.37 |
| STRCHIVE_LOCUS_ID | STRchive locus identifier, or the retained identifier for a legacy supplemental locus | STRchive | FXS_FMR1 |
| BENIGN_RANGES | Repeat-count range or ranges associated with a benign result | STRchive/Broad thresholds | 5-44 |
| INTERMEDIATE_RANGES | Repeat-count range or ranges associated with an intermediate result; `.` if none | STRchive/Broad thresholds | 45-200 |
| PATHOGENIC_RANGES | Repeat-count range or ranges associated with a pathogenic result; `*` indicates no upper limit | STRchive/Broad thresholds | 201-* |
| TARGET_REGION | hg38 coordinates of the repeat component shown in the report | STRchive/Broad catalog | chrX:147912050-147912110 |
| TARGET_VARIANT_ID | Identifier of the repeat component shown in the report | STRchive/Broad catalog | FXS_FMR1 |
| STRCHIVE_URL | Link to the corresponding STRchive locus page; `.` if unavailable | STRchive | https://strchive.org/loci/fxs_fmr1/ |
