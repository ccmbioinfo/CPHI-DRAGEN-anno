# Annotating compound heterozygous variants across variant types in short-read WGS

Madeline Couse

**Version 2026-04**

## Changelog

Adapted from the long-read [crg2-pacbio pipeline](https://github.com/ccmbioinfo/crg2-pacbio) for GRCh38 DRAGEN v4.4 short-read genomes. Note that variant phasing here will mostly rely on the availability of parental genotypes, as read-only variant phasing is only possible across very short distances with short-read WGS.

## Variant filtering and preparation for phasing

### Sequence variants

For SNVs/indels, variants are first filtered by the gnomAD filtering allele frequency (less than 1%). We are interested in rare nonsynonymous variants as well as rare synonymous/noncoding (deep intronic, UTR) variants that may be damaging. So, when determining gene compound heterozygous (CH) status in the ‘wgs.coding’ report, we also consider rare genic synonymous/noncoding variants with the following filter:

  - SpliceAI \>= 0.2 OR CADD \>= 10 OR abs(promoterAI) \>= 0.1 OR indel without CADD score (indels do not receive SpliceAI scores, and only common indels are scored by CADD)

These synonymous/noncoding variants do not appear in the ‘wgs.coding’ report (unless they are in ClinVar) but are still considered in the CH status determination. Note that though variants that VEP has annotated as having the impact ‘intergenic\_variant’ are not considered in the CH status determination, those annotated as ‘upstream\_gene\_variant’ and ‘downstream\_gene\_variant’ are. VEP annotates variants within 5kb upstream or downstream as ‘upstream\_gene\_variant’ and ‘downstream\_gene\_variant’ respectively.

Also note that if an SNV/indel overlaps multiple genes, the annotated gene is the gene with the higher impact. For example, if a variant is intronic in one gene and coding in the other, the variant will be annotated against the latter. If it is a tie between impact severities, the gene is effectively chosen randomly. These rules are imposed by the [gemini database](https://gemini.readthedocs.io/en/latest/index.html) we use to store and query variants.

### SVs

The SV report contains all SVs, so we work directly with the report. As with sequence variants, we do not consider intergenic variants in CH status determination, with the exception of those within 5kb upstream or downstream of a gene. We also exclude variants with \>1% filtering allele frequency in gnomAD.

### CNVs

Like SVs, we work directly with the CNV report. As with sequence variants, we do not consider intergenic variants in CH status determination, with the exception of those within 5kb upstream or downstream of a gene.

## Determining CH status

### CH status using proband-only short-read phasing

With short-read sequencing, we can only phase variants in the same read or read-pair. 

  - Consider only heterozygous variants in proband
  - If there is only one variant in a gene, the CH status of the gene is ‘False’.
  - If there are two or more variants in a gene, and neither is phased, the CH status of the gene is ‘Unknown’.
  - If there are two or more phased variants in one gene belonging to the same phase block, and the variants are on different haplotypes (e.g. 1|0 and 0|1 genotypes), the CH status of the gene is ‘True’. If the variants are on the same haplotype (e.g. 1|0 and 1|0 genotypes), the CH status is ‘False’.
  - If there are two or more variants in a gene belonging to more than one phase block, the CH status may be ‘True’ if there are at least two variants in trans in at least one phase block, or ‘Unknown’ otherwise (as we cannot know the relative phasing between different phase blocks).

### CH status using parental genotypes

  - Consider only heterozygous variants in proband where:
    
      - Neither parent is homozygous for the alternate allele
    
      - Both parents are not heterozygous
  - If mother is heterozygous and father is reference, variant belongs to maternal haplotype
  - If father is heterozygous and mother is reference, variant belongs to paternal haplotype
  - If both parents are homozygous reference or one or both parent genotypes are missing, haplotype is unknown

If there is at least one variant on the maternal haplotype and one variant on the paternal haplotype in a gene, the CH status of that gene is ‘True’.

If there is at least one variant on a parental haplotype and at least one on an unknown haplotype, the CH status of the gene is ‘Unknown’.

If there are multiple variants with unknown haplotype, the CH status of the gene is ‘Unknown’.

### CH status determination: implementation

CH algorithm inputs:

  - A table of all medium and high impact variants in the family\*
  - A table of all low impact variants in the family\*
  - Family ‘wgs.coding’ report
  - Family SV report
  - Family CNV report
  - A table with Ensembl, HGNC, and NCBI gene IDs
  - Family pedigree file

\* See documentation for variant impact classification [here](https://gemini.readthedocs.io/en/latest/content/database_schema.html#details-of-the-impact-and-impact-severity-columns).

CH algorithm outputs:

  - ‘wgs.coding’ report CSV with two additional columns, ‘CH\_variant\_types’ and ‘CH\_status’, that indicate if a gene contains at least one pair of heterozygous variants in *trans* in the proband, and what type of variants are in the gene.
  - SV report CSV with two additional columns, ‘CH\_variant\_types’ and ‘CH\_status’, as above.
  - CNV report CSV with two additional columns, ‘CH\_variant\_types’ and ‘CH\_status’, as above.
  - A ‘\<family\>\_compound\_het\_variants’ CSV that comprises all sequence variants, SVs, and CNVs that are heterozygous in the proband, with gene CH status annotated (see Table 1). Sequence variants are annotated with gnomAD AF, CADD score, SpliceAI score, and ClinVar. Refer to this table to view CH status for all genes. If there is only one sequence variant in a gene that is annotated as CH in the ‘wgs.coding’ report, refer to this table to see additional variant(s) in the gene (e.g a deep intronic variant that would not be included in the ‘wgs.coding’ report). Note – SVs and CNVs that overlap multiple genes will be spread over multiple rows, with one gene per row.

#### CH algorithm

Sequence variant processing

1.  Extract sequence variants from gemini db.

2.  Filter low impact sequence variants\*:
    
    1.  Variant is not intergenic
    
    2.  SpliceAI \>= 0.2 OR CADD \>= 10 OR abs(promoterAI) \>= 0.1 OR indel without CADD score (indels do not receive SpliceAI scores, and only common indels are scored by CADD)

3.  Index sequence variants by chromosome, position, ref allele, alt allele

4.  Create a sequence variant table in long format with Ensembl gene ID and per-sample genotype information (genotype, zygosity, phase set ID).

\*low impact variants present in ClinVar are exempt from these filters.

SV processing

5.  Index SVs by chromosome, start position, end position, SVTYPE, and pbsv SV ID.

6.  Filter for high impact SVs: filter out intergenic SVs and SVs with a gnomAD allele frequency less than 1%.

7.  Create an SV table in long format with Ensembl gene ID and per-sample genotype information (genotype, zygosity, phase set ID). Split variants by gene so that an SV that overlaps multiple genes is represented by multiple rows with one Ensembl gene ID per row.

CNV processing

8.  Repeat steps 5-7 for CNVs.

9.  Concatenate sequence, SV, and CNV tables.

CH status determination

10. Identify genes with CH status in proband using only proband long-read phased genotypes and phase block IDs from the variant table created in step 9.

11. For genes with unknown CH status after step 9, attempt to determine CH status using parental variant genotypes if available.

12. Export proband heterozygous variants with CH gene annotations to a CSV in the format described in Table 1 below.

Annotate variant reports with CH status

13. Add per-gene CH status and variant type(s) to small variant reports, SV, and CNV reports as ‘CH\_variant\_types’ and ‘CH\_status’ columns. For SVs/CNVs overlapping multiple genes, the per-gene status and variant type(s) will be delimited by the ‘|’ symbol. If an SV is annotated as ‘intergenic\_variant’ by SnpEff, or is present in gnomAD at greater than 1% frequency, the CH\_status is ‘.’. Note that rare variants annotated against genes with an ‘upstream\_variant’ or ‘downstream\_variant’ impact ARE included in the CH status determination.

| Column | Description | Example |
| --- | --- | --- |
| Variant_id | For sequence variants: chr-pos-ref-alt<br>For SVs: chr-start-end-SVTYPE-SVID<br>For CNVs: chr-start-end-SVTYPE | 15-48781193-T-C |
| Variant_type | One of: sequence_variant, SV, CNV | sequence_variant |
| Gene_name | HGNC gene symbol | CEP152 |
| Ensembl_gene_id | Ensembl gene ID | ENSG00000103995 |
| CH_status | Compound heterozygous status:<br>TRUE if there are at least two variants in the gene in *trans*<br>FALSE if there is only one variant in the gene, or all variants are in *cis*<br>UNKNOWN if there are at least two variants in a gene but phase cannot be determined (e.g. parental data is not available, and the variants reside in different phase blocks in the gene) | TRUE |
| `Zygosity_<sample>` | Variant zygosity for this particular sample. | het |
| Variation | Most severe impact of sequence variant. Only available for sequence variants in this table. | stop_gained |
| Clinvar | Clinvar pathogenic status. Only available for sequence variants in this table. | Pathogenic |
| Gnomad_af_grpmax | Maximum allele frequency across genetic ancestry groups in joint (WES and WGS) dataset. Only available for sequence variants in this table. | 0.000946139974985 |
| Cadd_score | CADD score. Only available for sequence variants in this table. | 27.9 |
| SpliceAI_score | SpliceAI score. Only available for sequence variants in this table. | C\|CEP152\|0.00\|0.00\|0.00\|0.00\|-1\|-7\|16\|-7 |
| promoterAI_score | promoterAI score. Only available for sequence variants in this table. | . |
| Nucleotide_change_ensembl | Ensembl HGVS. Only available for sequence variants in this table. | ENST00000325747.9:c.1755T\>G |
| `GT_abstracted_<sample>` | Genotype for phased variants is represented as ‘ref\|alt’ or ‘alt\|ref’ to facilitate a quick visual check for compound heterozygosity. For unphased variants, genotypes are left as is (e.g. T/A). | ref\|alt |
| `PS_<sample>` | Phase set ID. Note that this will be missing for the majority of variants called from short-read sequencing data. Variants that belong to the same phase block (ie have the same PS ID) have phased genotypes that can be interpreted relative to one another. So, variants with genotypes of ‘ref\|alt’ and ‘alt\|ref’ are in *trans*, while variants with genotypes of ‘ref\|alt’ and ‘ref\|alt’ are in *cis*. If variants reside in different phase blocks, phased genotypes *cannot* be interpreted relative to one another. | 48686679 |

Table 1. Family compound heterozygosity table.
