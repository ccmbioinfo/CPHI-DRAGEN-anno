# Small-variant report update (CRE/GEMINI to slivar)

**Version 2026-08**

## Overview

Small-variant reporting in the DRAGEN short-read and PacBio long-read pipelines is being updated from CRE/[GEMINI](https://gemini.readthedocs.io/en/latest/)/R to [slivar](https://github.com/brentp/slivar) and Python. The reasons for this change are:

- GEMINI and related tools are no longer actively maintained
- The CRE/GEMINI/R workflow does not provide special handling for disease-associated non-coding genes
- slivar uses considerably fewer computational resources than GEMINI
- more of the codebase can be kept in Python and maintained consistently
- variant filtering is easier to inspect and control directly
- transcript and consequence selection can be defined explicitly in the pipeline code

In the current workflow, CRE/GEMINI loads each VEP/vcfanno-annotated small-variant VCF into a database, where [geneimpacts](https://github.com/brentp/geneimpacts) selects one primary transcript and consequence per variant. Most gene- or consequence-dependent filtering and transcript fields use that annotation, and an R script formats the final report.

The updated workflow filters the annotated VCF directly with slivar without an intermediate database, allowing all VEP annotations to be considered when gene or consequence context matters. A Python script selects a primary annotation later for some report fields, while the new `_all` columns summarize the additional gene, transcript, and consequence context. The additional context is particularly important for variants that overlap both non-coding and protein-coding genes, for example RNU4-2 and SIRT4, since protein-coding transcripts, unless specified in the curated disease-associated non-coding gene list (see Appendix), are always prioritized over non-coding transcripts.

## Variant filtering

The slivar filtering criteria remain very similar to those used by CRE/GEMINI and described in the [DRAGEN small-variant report documentation](https://github.com/ccmbioinfo/CPHI-DRAGEN-anno/blob/main/docs/variant_report_docs/DRAGEN_small_variant_report.md).

To assess filtering parity, final CRE/GEMINI and slivar reports from one DRAGEN and three PacBio datasets were compared using unique Position + Ref + Alt keys. No variants were found only in the CRE/GEMINI reports; all observed differences were variants additionally retained by slivar.

The table reports slivar-only variants as a percentage of all unique variants in each final slivar report (slivar-only / total slivar). N/A means that no comparable report was available.

| **Report** | **DRAGEN Sample A** | **PacBio Sample B** | **PacBio Sample C** | **PacBio Sample D** |
| --- | --- | --- | --- | --- |
| Coding | 1.6% | 3.3% | 2.5% | 2.1% |
| High-impact | 10.0% | 12.3% | 10.7% | N/A |
| Denovo | 0% | N/A | N/A | N/A |
| Panel | \<0.1% | 3.0% | 2.7% | \<0.1% |
| Panel-flank | \<0.1% | 3.1% | 2.6% | \<0.1% |

### Why slivar retains additional variants

The additional variants are retained in the slivar workflow for four broad reasons:

- **Final CRE depth filter:** Variants can enter a report through different filtering paths. The main path requires an alternate-read depth (AD) of at least three, while the ClinVar rescue path accepts one or two alternate reads. CRE later (erroneously) reapplies AD \>= 3 to the combined report, removing some rescued variants. Slivar keeps the threshold assigned to each path and does not apply a final blanket depth filter.
- **Partial PacBio allele depth:** PacBio VCFs can contain partial values such as AD=.,4. These values occur most often when genotype quality is missing; however this is not always the case. During later CRE/R report processing, the partial value can be interpreted as missing and converted to zero, causing the final depth filter to remove the variant. Slivar reads the available alternate count directly. For example, PacBio Sample B contains a variant with GT:AD:DP=./1:.,4:11; CRE omitted this variant, while slivar retained it with an alternate depth of four.
- **CRE selected a LOW consequence:** geneimpacts prefers annotations it classifies as coding before comparing consequence severity. It can therefore select a LOW consequence on an exact protein_coding transcript even when another annotation has an impactful consequence on a pseudogene-related, protein_coding_LoF, or nonsense_mediated_decay transcript. CRE coding-report filtering excludes the selected LOW annotation; slivar considers all annotations and retains the variant when any annotation has an impactful consequence.
- **CRE selected a blank gene:** The CRE high-impact report removes a variant when its selected annotation has no gene symbol. Slivar can retain it when another annotation for the same variant has a usable symbol.

The following table uses PacBio Sample B as a representative breakdown. Categories are mutually exclusive within each report and identify the first point at which CRE lost the variant.

| **PacBio Sample B category** | **Coding** | **High-impact** | **Panel** | **Panel-flank** |
| --- | --- | --- | --- | --- |
| Total additional variants | 194 (3.3%) | 893 (12.3%) | 823 (3.0%) | 2,720 (3.1%) |
| Partial AD | 68 | 376 | 821 | 2,712 |
| CRE final AD \>= 3 filter | 25 | 3 | 2 | 8 |
| CRE-selected LOW consequence | 101 | 0 | 0 | 0 |
| CRE-selected blank gene | 0 | 514 | 0 | 0 |

## Transcript and consequence selection

Transcript selection is the main area of difference between the CRE/GEMINI and slivar reports.

### Previous geneimpacts selection

During database creation, CRE/GEMINI uses geneimpacts to select one primary annotation for each variant. Most transcript-dependent filtering and single-value report fields are based on that annotation.

In broad terms, geneimpacts compares annotations in this order:

1. **Pseudogene status:** non-pseudogene annotations are preferred over pseudogene annotations.

2. **Coding status:** prefer protein_coding annotations with an exonic, non-UTR consequence. However, if an annotation has a splice-related consequence, skip this step and continue to step 3.

3. **Broad consequence severity:** consequences are assigned to the three groups HIGH, MED, and LOW using consequence mapping bundled with geneimpacts.

4. **Biotype:** when severity ties, an exact protein_coding biotype is preferred, followed by processed_transcript.

5. **Final tie-breakers:** SIFT and PolyPhen scores, followed by the internal consequence order, are used to resolve further ties.

Remaining ties can depend on input order, contributing to cases where CRE selects a blank gene symbol even though another similarly ranked annotation had a gene.

### Current slivar/Python selection

Where filtering depends on a gene or consequence, slivar considers all VEP annotations. A Python script later selects one primary annotation for most single-value report columns (for example, `Gene`), while the `_all` columns (for example, `Gene_all`) summarize the additional context.

The following priorities are applied in order. A later priority is considered only when the earlier priorities tie:

1. **Curated non-coding gene:** a list of disease associated RNA and mitochondrial genes are preferred over overlapping protein coding genes. Adapted from Illumina Emedgene documentation (see Appendix).

2. **Pseudogene status:** non-pseudogene annotations are preferred. This means a lower-impact non-pseudogene annotation can be displayed ahead of an impactful pseudogene annotation.

3. **Impactful status:** among annotations with the same pseudogene status, consequences above IMPACTFUL_CUTOFF are preferred. This approximately matches the former HIGH/MED versus LOW boundary; the full order is provided in the appendix.

4. **Biotype:** an exact protein_coding biotype is preferred, followed by related coding biotypes such as protein_coding_LoF, nonsense_mediated_decay, and non_stop_decay. Other biotypes follow.

5. **Consequence rank:** the highest-ranked consequence in the Appendix is preferred. If one annotation contains multiple consequences, its highest-ranked consequence is used.

6. **Gene-symbol quality:** conventional named symbols are preferred over LOC/LINC, small-RNA-style (MIR, RNU, SNORD, SNORA, SCARNA, Y_RNA), Ensembl-style, and clone-style (AC, AL, AP, RP) symbols. These are still preferred over a missing symbol.

7. **Transcript-quality flags:** MANE Select, MANE Plus Clinical, and then VEP canonical status are used only after the preceding priorities tie.

8. **Final tie-breakers:** processed or other transcript classes, genic status, transcript ID, and original VEP order are used for final ties.

## Report columns affected by transcript selection

The primary transcript selection can directly change Gene, Ensembl_gene_id, Ensembl_transcript_id, Variation, Info, Refseq_change, coding consequence fields, and transcript-level gnomAD constraint values. It can indirectly change gene-dependent Burden, HPO, and compound-het fields. ACMG secondary findings now check Gene_all, reducing dependence on the primary gene.

## New and replaced columns

Most established columns remain. The additions below provide context that was previously reduced to one GEMINI-selected gene/transcript annotation.

| **Column** | **Comment** |
| --- | --- |
| CSQ_biotype, CSQ_impact | New columns showing the biotype and VEP impact category of the selected primary annotation. |
| Gene_all, Ensembl_gene_id_all, Ensembl_transcript_id_all, Variation_all, Info_all, Refseq_change_all | New companions to the corresponding primary columns; preserve all gene, transcript, and consequence. |
| CSQ_biotype_all, CSQ_impact_all, Sift_score_all, Polyphen_score_all, Constraint_all | New transcript-context summaries; individual constraint columns still describe the primary transcript. |
| Gene_description_all, omim_phenotype_all, omim_inheritance_all, Orphanet_all | Replace the respective Gene_description, omim_phenotype, omim_inheritance, and Orphanet columns and collect annotations across Gene_all. |

## Effect on compound-het annotation

Compound-het annotations are associated with report variants using the selected Ensembl_gene_id. Because CRE/GEMINI and slivar may select different transcript annotations for the same variant, they may also assign it to different Ensembl gene IDs. This can change CH_status and CH_variant_types even though the underlying variant call has not changed.

However, all three PacBio positive controls produced the expected compound-het results in the slivar reports. PacBio Sample C remains UNKNOWN, as expected, because the duplication is unphased:

<table>
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr>
<th><strong>Sample</strong></th>
<th><strong>Gene</strong></th>
<th><strong>Variant</strong></th>
<th><strong>Expected CH status</strong></th>
<th><strong>slivar CH status</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td>PacBio Sample A</td>
<td>GNPNAT1</td>
<td>Small-indel/<br />
SNV</td>
<td>TRUE</td>
<td>TRUE</td>
</tr>
<tr>
<td>PacBio Sample B</td>
<td>CEP152</td>
<td>SNV/<br />
SNV</td>
<td>TRUE</td>
<td>TRUE</td>
</tr>
<tr>
<td>PacBio Sample C</td>
<td>TGM1</td>
<td style="text-align: left;">SV/SNV</td>
<td>UNKNOWN</td>
<td>UNKNOWN</td>
</tr>
</tbody>
</table>

## Validation

We ran the new slivar-based workflow on two families with DRAGEN short-read sequencing and two families with PacBio long-read sequencing to validate the filtering and annotation protocols. All diagnostic/potentially diagnostic variants were present in the slivar reports and annotated appropriately.

## Appendix: slivar consequence order

CRE/GEMINI primarily compared consequences using the broad HIGH, MED, and LOW groups defined by geneimpacts. Slivar requires a fully ordered list. The current list was designed to match the geneimpacts boundary as closely as possible: consequences historically treated as HIGH or MED are above IMPACTFUL_CUTOFF, while consequences historically treated as LOW are below it. The explicit order within those groups is based on the default [slivar-ordering](https://github.com/brentp/slivar/blob/master/src/slivarpkg/default-order.txt).

chromosome_number_variation\
transcript_ablation\
splice_acceptor\
splice_donor\
stop_gained\
exon_loss\
frameshift\
stop_lost\
start_lost\
transcript_amplification\
gene_fusion\
protein_protein_contact\
structural_interaction_variant\
bidirectional_gene_fusion\
start_retained\
rare_amino_acid\
feature_fusion\
feature_ablation\
disruptive_inframe_insertion\
disruptive_inframe_deletion\
inframe_insertion\
inframe_deletion\
missense\
conservative_inframe_insertion\
conservative_inframe_deletion\
protein_altering\
inframe_altering\
splice_region\
5_prime_utr_truncation\
5_prime_UTR_truncation\
initiator_codon\
duplication\
inversion\
exon_region\
regulatory_region_ablation\
splice_donor_5th_base_variant\
splice_donor_region_variant\
splice_polypyrimidine_tract_variant\
IMPACTFUL_CUTOFF\
incomplete_terminal_codon\
stop_retained\
3_prime_utr_truncation\
3_prime_UTR_truncation\
synonymous\
gene\
exon\
coding_sequence\
mature_miRNA\
mature_mirna\
5_prime_UTR_premature_start_codon_variant\
5_prime_utr_premature_start_codon\
5_prime_UTR_premature_start_codon_gain_variant\
5_prime_utr_premature_start_codon_gain_variant\
5_prime_UTR\
5_prime_utr\
3_prime_UTR\
3_prime_utr\
coding_transcript\
miRNA\
non_coding_transcript_exon\
non_coding_exon\
non_canonical_start_codon\
nc_transcript\
conserved_intron\
GENIC_CUTOFF\
intron\
NMD_transcript\
nmd_transcript\
non_coding_transcript\
non_coding\
transcript\
upstream_gene\
sequence_feature\
sequence\
downstream_gene\
intragenic\
TFBS_ablation\
TFBS_amplification\
TF_binding_site\
tf_binding_site\
regulatory_region_amplification\
feature_elongation\
regulatory_region\
conserved_intergenic\
feature_truncation\
intergenic\
intergenic_region

## Appendix: curated disease associated RNA genes

Curated RNA and non-coding disease genes, adapted from [Illumina emedgene's](https://help.connected.illumina.com/emedgene/emedgene-analyze-manual/tertiary-analysis-pipeline/how_does_emedgene_analyze_prioritize_transcripts) prioritized RNA-gene list, with the addition of RNU6-1, RNU6-2, RNU6-8, RNU6-9, RNU6ATAC. These small non-coding disease genes are usually annotated by VEP under an overlapping protein-coding gene, so they lose the default primary transcript selection. Listing a gene here promotes it to the primary report gene when a variant overlaps it.

ATXN8OS, CHASERR, H19, KCNQ1OT1, MIAT, MIR17HG, MIR137, MIR140, MIR184, MIR19B1, MIR204, MIR2861, MIR4718, MIR605, MIR96, MIR99A, RMRP, RNU2-2, RNU4-2, RNU4ATAC, RNU5A-1, RNU5B-1, RNU6-1, RNU6-2, RNU6-8, RNU6-9, RNU6ATAC, RNU7-1, RNU12, SNORA31, SNORD116-1, SNORD118, TERC, TRU-TCA1-1, MT-RNR1, MT-RNR2, MT-TA, MT-TC, MT-TD, MT-TE, MT-TF, MT-TG, MT-TH, MT-TI, MT-TK, MT-TL1, MT-TL2, MT-TM, MT-TN, MT-TP, MT-TQ, MT-TR, MT-TS1, MT-TS2, MT-TT, MT-TV, MT-TW, MT-TY
