# Mitochondrial Report

Michelle Spivak, Madeline Couse, Anjali Jain

**Version 2026-04**

## Changelog

- Adapted from crg2 pipeline for GRCh38 DRAGEN v4.4 genomes. 

- Annotation databases updated to hg38-specific resources: gnomAD v3.1 chrM sites and hg38 ClinVar chrM VCF

- MITOMAP annotations originating from multiple source tables (confirmed mutation, disease, RNA, coding, control) are collapsed into combined report fields based on source priority to avoid redundant columns. Empty values are omitted from collapsed MITOMAP fields.

## Summary

- Mitochondrial variant calls from DRAGEN whole genome sequencing are normalised and decomposed with `bcftools norm` and annotated using `mity report`(https://github.com/KCCG/mity) v2.0.1. MITOMAP annotation tables were bundled with mity v2.0.1.

- Note: variants are reported in left-normalised, decomposed format: multi-base indels and complex variants are split by `bcftools norm`. Multi-nucleotide variants may be lost during decomposition and excluded from the final report. 

## Recommended filters for analysis

  - 1. gnomAD_AF_hom <= 0.01
  - 2. MITOMAP CONFIRMED DISEASE != ‘.’
  - 3. Consider variant_heteroplasmy in unaffected individuals if available

## Column descriptions

| **Column** | **Comment** | **Source** | **Example** |
| --- | --- | --- | --- |
| `CHR` | Chromosome identifier | mity | `chrM` |
| `POS` | Variant position in the mitochondrial genome | mity | `591` |
| `REF` | Reference allele | mity | `C` |
| `ALT` | Alternate allele | mity | `T` |
| `HGVS` | Variant annotation using HGVS nomenclature, assuming SNP, insertion or deletion | mity | `m.591C>T` |
| `GENE/LOCUS` | The gene the variant falls in, or if not in a gene, the MITOMAP locus entry if it exists | mity | `MT-TF` |
| `GENE/LOCUS DESCRIPTION` | Short description of the gene or locus | mity | `tRNA phenylalanine` |
| `COHORT COUNT` | Number of times this variant (defined by position and alternate allele) is observed across all input VCFs | mity | `1` |
| `SAMPLE` | Name of the sample in which the variant is observed | mity | `FamID_SampleID` |
| `<SAMPLE>.VARIANT HETEROPLASMY` | Variant allele fraction: `alt_depth / (ref_depth + alt_depth)`. Values `>=0.95` are typically considered homoplasmic. | Mity (`mt_report`) | `0.200` |
| `<SAMPLE>.ALT DEPTH` | Alternate allele read depth for this sample | mity | `2000` |
| `<SAMPLE>.TOTAL SAMPLE DEPTH` | Sum of `ref_depth` and `alt_depth` for this sample | mity | `25000` |
| `MITOMAP AA CHANGE` | Amino acid or functional change annotation from MITOMAP. For RNA variants this may show the affected tRNA; for coding variants it may show amino acid change. | Mity (`mt_report`) | `Confirmed: tRNA Phe; Coding: tRNA` |
| `clinvar_significance` | ClinVar pathogenicity classification | ClinVar via vcfanno | `Likely_pathogenic` |
| `clinvar_status` | ClinVar review status reflecting evidence/review tier. | ClinVar via vcfanno | `reviewed_by_expert_panel` |
| `gnomAD_AC_hom` | Allele count for variants with heteroplasmy level `>=0.95` in gnomAD chrM | gnomAD via vcfanno | `3` |
| `gnomAD_AC_het` | Allele count for variants with heteroplasmy `>=0.10` and `<0.95` in gnomAD chrM | gnomAD via vcfanno | `6` |
| `gnomAD_AF_hom` | Allele frequency for variants with heteroplasmy `>=0.95` | gnomAD via vcfanno | `0.7157` |
| `gnomAD_AF_het` | Allele frequency for variants with heteroplasmy `>=0.10` and `<0.95` | gnomAD via vcfanno | `0.00010683` |
| `gnomAD_max_hl` | Maximum heteroplasmy level observed across all gnomAD samples for this variant | gnomAD via vcfanno | `0.2409999966621400` |
| `MITOMAP LOCUS` | MITOMAP locus annotation. | Mity | `MT-TF` |
| `MITOMAP CONFIRMED DISEASE` | Disease name from the MITOMAP Confirmed Mutations table only. This is the primary curated disease-association field for clinically confirmed variants. Empty if the variant is not present in the Confirmed Mutations table. | Mity | `Gitelman-like syndrome` |
| `MITOMAP ADDITIONAL REPORTED DISEASE` | Additional reported disease association from the MITOMAP RNA Mutations table or Disease table, shown when different from the Confirmed disease association(s). | Mity (`mt_report`) | `Gitelman-like-syndrome` |
| `MITOMAP SIGNIFICANCE SUMMARY` | A derived field created by collapsing MITOMAP significance/status values from the Confirmed Mutations, RNA, and Disease status fields into a single summary category. | Mity (`mt_report`) | `Cfrm[LP]` |
| `MITOMAP HOMOPLASMY` | Whether the variant has been observed homoplasmic in published MITOMAP disease cases, shown here as a collapsed/source-labeled field. | Mity | `RNA: +; Disease: +` |
| `MITOMAP HETEROPLASMY` | Whether the variant has been observed heteroplasmic in published MITOMAP disease cases, shown here as a collapsed/source-labeled field. | Mity | `RNA: -; Disease: -` |
| `MITOMAP DISEASE PUBMED IDS` | PubMed IDs associated with the variant from the MITOMAP Disease table. | Mity (`mt_report`) | `.` |
| `MITOMAP MUTATIONS RNA MITOTIP` | Pathogenicity call from the MITOMAP RNA mutations table. This reflects MITOMAP’s disease-report-based RNA interpretation and is often more reliable than the computational MitoTIP interpretation for already-known MITOMAP variants. | Mity | `Pathogenic` |
| `MITOTIP SCORE` | Computational pathogenicity score from the MitoTIP model for tRNA variants; higher scores indicate greater predicted pathogenicity. | Mity | `3.16860` |
| `MITOTIP SCORE INTERPRETATION` | Qualitative interpretation of `MITOTIP SCORE` | Mity | `likely-benign` |
| `GB_COUNT` | Number of MITOMAP full-length GenBank sequences (`>15.4 kb`) carrying this variant. | Mity | `0` |
| `GB_PERCENTAGE` | Frequency of this variant in MITOMAP's GenBank full-length sequence set (`COUNT / total sequences`). | Mity | `0.0` |
| `MITOMAP # REFERENCES` | Publication count for the variant, labeled by source where relevant, such as RNA or Coding. | Mity (`mt_report`) | `RNA: 2` |
| `MITOMAP STATUS` | An indication of the strength of evidence supporting the MITOMAP annotation. The strongest evidence is "Confirmed"| Mity | `Reported_but_unconfirmed` |
| `PHYLOTREE MUT` | Phylotree Mutation | Mity | `A73G!` |
| `PHYLOTREE HAPLOTYPE` | The haplotypes that the variant is known to contribute to. See `http://www.phylotree.org/index.htm` | Mity | `H13a2b3` |

## Column changes from previous report version

| **Previous column** | **Comment** |  **Source** |**Description of Change** |
| --- | --- | --- | --- |
| `disease_amino_acid_change_mitomap` | Amino acid change in the variant if applicable | Mity | Replaced by `MITOMAP AA CHANGE` with source prioritization |
| `allele_frequency_mitomap` | MITOMAP derives this from 32,059 GenBank full length sequences, with size greater than 15.4kbp (range 0-1) | Mity | Replaced by `GB_COUNT` |
| `GenBank_frequency_mitomap` | The frequency of GenBank records containing this variant. This is a proxy for the population frequency of the variant. As expected, none of the confirmed pathogenic variants in MITOMAP have GenBank frequency >0.0% | Mity | Replaced by `GB_PERCENTAGE` |
| `codon_position_mitomap` | Position of the codon where the change occurs | Mity | Removed, considered uninformative |
| `codon_number_mitomap` | Number of the codon that is affected | Mity | Removed, considered uninformative |
| `RNA_mitomap` | Affected RNA (tRNA, 12S rRNA, 16S rRNA, noncoding MT-TF precursor etc) | Mity | Removed, redundant with `GENE/LOCUS DESCRIPTION` |
| `commercial_panels` | This indicates whether a variant is commonly tested in a number of commercial/academic testing laboratories. | Mity | No longer available |
| `MitoTip_percentile` | Where this variant’s MitoTip score sits in relation to all of the possible variants in tRNA. `100%` is the most pathogenic, `0%` is a score of zero. | Mity | Removed, redundant |
| `anticodon` | `TRUE` if the variant falls in an anticodon of a tRNA | Mity | Removed, uninformative |
| `MGRB_frequency` | Represents the frequency of this variant in the Medical Genomics Reference Bank (MGRB) cohort | Mity | Removed MGRB database information due to limited usability |
| `MGRB_FILTER` | Represents the `FILTER` status of the variant in the MGRB cohort using GATK HaplotypeCaller | Mity | Removed MGRB database information due to limited usability |
| `MGRB_AC` | The allele count in the MGRB cohort; i.e. the number of individuals with this variant. This is not accurate for the MT, as it assumes each individual is diploid. Interpret this value as `2x` the number of individuals in the MGRB control cohort. | Mity | Removed MGRB database information due to limited usability |
| `MGRB_AN` | The total number of alleles called in the MGRB cohort. This is not accurate for the MT, as it assumes each individual is diploid. Interpret this value as `2x` the number of individuals in the MGRB control cohort. | Mity | Removed MGRB database information due to limited usability |