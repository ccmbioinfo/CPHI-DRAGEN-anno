# Annotating and filtering SV and CNV calls from DRAGEN whole genome sequencing

Rohan Khan, Madeline Couse

**Version 2026-03**

## Changelog

Changes from the previous `DRAGEN_SV_CNV_report.md` version

Differences from `crg2-pacbio`

  - This workflow annotates DRAGEN SV and CNV VCFs rather than PacBio pbsv and HiFiCNV calls
  - Population annotations use gnomAD SV, DGV, and a DRAGEN-derived 1000 Genomes SV database rather than PacBio-specific C4R, TG, and CoLoRSdb resources
  - Sample support fields are DRAGEN FORMAT fields (`GT`, `PR`, `SR`, `VF`, `GQ`, `FS`, `DQ`, `DN`, `CN`) rather than PacBio long-read genotype/read-depth and phasing fields
  - `CPHI-DRAGEN-anno` converts inversion-like BND records to `INV` before SV annotation

## Summary

SV and CNV reports are generated from joint-genotyped or singleton DRAGEN VCFs.

Annotate unfiltered SV or CNV calls from DRAGEN sequencing output using [SnpEff](https://pcingola.github.io/SnpEff/) and [AnnotSV](https://github.com/lgmgeo/AnnotSV). SnpEff predicts gene impacts (e.g. if a variant results in a frameshift, or falls in an intron). AnnotSV provides transcript overlap, repeat, regulatory, and gene constraint annotations that assist in clinical interpretation. Gene impact annotations used in the final report are taken from SnpEff rather than AnnotSV gene annotations because SnpEff incorporates 5kb upstream and downstream of genes and would therefore capture, say, an SV impacting a gene promoter.

For structural variants, `CPHI-DRAGEN-anno` converts DRAGEN inversion-like BND records to `INV` before annotation so that inversions are represented more consistently downstream.

[SnpEff](https://pcingola.github.io/SnpEff/) 5.4.0a for gene impact annotation

  - Annotates impact on every transcript overlapped by a variant which AnnotSV annotates

[AnnotSV](https://github.com/lgmgeo/AnnotSV) 3.1.1 for annotation of:

  - Repeats
  - Gene constraint metrics
  - Regulatory and regional annotations used in the final report

Custom python script for parsing and further annotation

  - Merge AnnotSV `split` and `full` annotations
  - Pull gene impact annotations from SnpEff annotated VCF, and merge these with AnnotSV annotations
  - Add exon overlap counts, Ensembl CDS overlap, SV/CNV start and end Ensembl gene overlaps
  - Add optional proband-specific HPO terms if an HPO file is provided in the config
  - Add OMIM
  - Add ClinGen annotations
  - Add SV/CNV frequencies from gnomAD SV, DGV Gold Standard, and a 1000 Genomes SV database
  - Add overlap with Adotto tandem repeat regions used for repeat analysis

## Population databases used in report generation

  - gnomAD SV v4.1
  - DGV Gold Standard CNV/SV set for hg38
  - 1000 Genomes SV database

Matching logic used in report generation:

  - Insertions and BNDs are matched to population SVs within 50 bp
  - Insertions must also have at least 80% size similarity
  - Deletions, duplications, and inversions are matched using at least 80% reciprocal overlap
  - Sample and database variants must have the same `SVTYPE`
  - CNVs follow the DEL/DUP reciprocal-overlap matching path

## 1000 Genomes SV database

The 1000 Genomes SV database used in report generation was created from [DRAGEN SV callsets v4.2.7](https://registry.opendata.aws/ilmn-dragen-1kgp/) using the following workflow:

  - Download per-sample DRAGEN SV VCFs from the public 1000 Genomes DRAGEN S3 bucket
  - Convert inversion-like BND records to `INV` using a [`bnd_to_inv.py`](https://github.com/srbehera/DRAGEN_Analysis/blob/main/convertInversion.py) script so inversion representation is more consistent before aggregation
  - Merge the per-sample inversion-fixed VCFs with `bcftools merge`
  - Retain PASS calls and GT-only sample data, then adjust confidence interval fields
  - Split the merged callset into INS and DEL/DUP/INV groups; BNDs are separated earlier and are not carried into the final collapsed frequency database
  - Collapse similar variants with `truvari collapse`
  - INS collapse uses a tighter positional threshold (`-s 50`)
  - DEL/DUP/INV collapse uses a broader positional threshold (`-r 500`)
  - Add allele count and allele frequency tags with `bcftools +fill-tags`
  - Convert the final collapsed VCF to a tabular database with chromosome, coordinates, SV type, SV length, ID, allele frequency, allele count, and homozygous alternate count

## Recommended filters for analysis

  - `gnomad_maxAF < 0.01`
  - `DGV_maxAF < 0.01`
  - `1000G_maxAF < 0.01`

For analysis, we recommend considering any exon-overlapping SVs or CNVs seen in the proband:

  - `sampleID_zyg` uncheck `-` and `./.`

And also any genic SV/CNV either in OMIM and/or in genes matching HPO terms:

  - `omim_phenotype` uncheck `.`
  - If an HPO file was added in the config, `HPO` uncheck `nan`

Other analysis considerations:

  - Insertions/deletions within short tandem repeat (STR) regions can be difficult to place precisely. Consider the `Repeat_type_left`, `Repeat_type_right`, and `Adotto_tandem_repeat` columns when interpreting variants in repeat-rich regions.
  - For BND-derived SV records, the `INFO` column in the final SV report is extended to include breakpoint directionality information after DRAGEN BND-to-INV conversion logic is applied.
  - For SVs, use the sample-level support fields when assessing confidence in a candidate call. In particular, review `sampleID_PR_alt`, `sampleID_SR_alt`, `sampleID_VF_alt`, `sampleID_GQ`, and `sampleID_FS`.
  - For CNVs, review `sampleID_CN`, `sampleID_DQ`, and `sampleID_DN` instead of the SV read-support fields.
  - The `sampleID_GT` column reports the called genotype, where `0` is reference and `1` is alternate. The corresponding `sampleID_zyg` column simplifies this to `het`, `hom`, `-`, `./.`, or for some sex-chromosome CNVs `hemi`.
  - For very large SVs/CNVs (`|SVLEN| > 500000`), the `GENE_NAME` and `ENSEMBL_GENE` fields are condensed in the final report to keep the output manageable and Excel-compatible. In these cases, only the first and last gene entries are retained rather than the full gene list.

## Report columns

The SV and CNV reports share most columns. Columns that are SV-only or CNV-only are noted below.

| **Column\_id** | **Comment** | **Source** | **Example** |
| --- | --- | --- | --- |
| `CHROM` | Chromosome number | DRAGEN VCF | 7 |
| `POS` | Structural variant 5’ position | DRAGEN VCF | 248484 |
| `END` | Structural variant 3’ position | DRAGEN VCF | 248485 |
| `SVLEN` | The length in nucleotides of an SV or CNV. A negative `SVLEN` indicates a deletion. | DRAGEN VCF | 109 |
| `SVTYPE` | The SVTYPE: duplication (DUP), deletion (DEL), insertion (INS), inversion (INV), breakend (BND) | DRAGEN VCF | INS |
| `ID` | DRAGEN variant ID. Present in the SV report; excluded from the final CNV report. | DRAGEN VCF | DRAGEN:INS:207297:0:0:0:0:0 |
| `INFO` | VCF INFO field. For SVs, BND directionality information is retained/added where applicable. | Generated value from DRAGEN VCF INFO | END=248484;SVTYPE=INS; SVLEN=109;CIGAR=1M109I; CIPOS=0,50;HOMLEN=50; HOMSEQ=CCGCATTCATCTCTC CAGGTAGCCTGGCACGGGGAGC CCACATTCATCTC |
| `FILTER` | Filter field from VCF | DRAGEN VCF | PASS |
| `GENE_NAME` | Gene symbol(s) impacted by the SV/CNV | SnpEff | FAM20C |
| `ENSEMBL_GENE` | Ensembl gene ID(s) impacted by the SV/CNV | SnpEff | ENSG00000177706 |
| `GENE_NAME_CDS` | Gene symbol(s) whose CDS overlaps the SV/CNV | Generated value | FAM20C |
| `GENE_NAME_START` | Gene symbol(s) overlapped by the start coordinate | Generated value | FAM20C |
| `GENE_NAME_END` | Gene symbol(s) overlapped by the end coordinate | Generated value | FAM20C |
| `VARIANT` | Gene impact term(s) | SnpEff | intron\_variant; downstream\_gene\_variant |
| `IMPACT` | SnpEff impact level(s): high, moderate, low, modifier | SnpEff | MODIFIER |
| `UCSC_link` | Excel hyperlink formula to UCSC genome browser | Generated value | [UCSC\_link](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=7:248484-248485) |
| `omim_phenotype` | OMIM phenotype(s) associated with impacted gene(s) | OMIM | Raine syndrome |
| `omim_inheritance` | OMIM mode(s) of inheritance | OMIM | AR |
| `clingen_disease` | ClinGen disease name | ClinGen | complex neurodevelopmental disorder |
| `clingen_ classification` | ClinGen clinical validity classification | ClinGen | Moderate |
| `clingen_HI` | ClinGen haploinsufficiency score | ClinGen | 1.0 |
| `clingen_TS` | ClinGen triplosensitivity score | ClinGen | 0.0 |
| `clingen_region _curation` | ClinGen genomic disease region curation | ClinGen | 10p13 population region (DGV\_Gold\_Standard\_ June\_2021\_gssvL10668) |
| `clingen_region` | ClinGen genomic disease region coordinates | ClinGen | 10:13217247-13218742 |
| `sample_SV_clingen_ region_perc_overlap` | Percent of ClinGen region overlapped by sample SV/CNV | Generated value | 98.9 |
| `HPO` | HPO terms from the supplied HPO file associated with impacted genes; only present if an HPO file is added in the config | Phenotips-style input HPO table | nan |
| `sampleID_zyg` | Sample zygosity summary. For CNVs, this is derived with CNV-aware logic and may include `hemi` on sex chromosomes. | Generated value from DRAGEN VCF | het |
| `sampleID_GT` | Sample genotype from FORMAT field | DRAGEN VCF FORMAT | 0/1 |
| `sampleID_PR_alt` | Alternate paired-read support count (`PR` alt value). Present in SV report; excluded from final CNV report. | DRAGEN VCF FORMAT | 4 |
| `sampleID_SR_alt` | Alternate split-read support count (`SR` alt value). Present in SV report; excluded from final CNV report. | DRAGEN VCF FORMAT | 10 |
| `sampleID_VF_alt` | Number of DNA fragments support of the breakend mates for the alternate allele. Present in SV report; excluded from final CNV report. | DRAGEN SV VCF FORMAT | 11 |
| `sampleID_GQ` | Sample genotype quality. Present in SV report; excluded from final CNV report. | DRAGEN SV VCF FORMAT | 10 |
| `sampleID_FS` | Phred-scaled p-value using Fisher’s exact test to detect strand bias. Present in SV report; excluded from final CNV report. | DRAGEN SV VCF FORMAT | 7.395 |
| `sampleID_DQ` | De novo quality. | DRAGEN CNV VCF FORMAT | 35 |
| `sampleID_DN` | Possible values are ‘Inherited’, ‘DeNovo’ or ‘LowDQ’. Threshold for a passing de novo call is based on DQ value. | DRAGEN CNV VCF FORMAT | 0.12 |
| `sampleID_CN` | Estimated copy number. | DRAGEN CNV VCF FORMAT | 3 |
| `Tx` | Transcript(s) reported by AnnotSV split annotation | AnnotSV | NM\_020223 |
| `Frameshift` | Indicates if the CDS length is not divisible by three (yes or no) | AnnotSV | no |
| `EXONS_SPANNED` | The number of exons that a structural variant overlaps with | Generated value | 4 |
| `Nearest_SS_type` | Nearest splice site type: 5’ (donor) or 3’ (acceptor). Present in SV report; excluded from final CNV report. | AnnotSV | 5’ |
| `Dist_nearest_SS` | Distance to nearest splice site. Present in SV report; excluded from final CNV report. | AnnotSV | 102 |
| `gnomad_NAME` | Matching gnomAD SV identifier(s) | gnomAD SV | gnomAD-SV\_v3\_ DUP\_chr6\_cbd0f366 |
| `gnomad_GRPMAX_AF` | Group-max allele frequency of overlapping gnomAD SV(s) | gnomAD SV | 0.00755 |
| `gnomad_AC` | Allele count in gnomAD SV | gnomAD SV | 54 |
| `gnomad_HOM` | Homozygous alternate count in gnomAD SV | gnomAD SV | 0 |
| `gnomad_SV` | Coordinates/type of matching gnomAD SV(s), where matching is determined using the report’s SVTYPE-specific comparison logic | gnomAD SV | DUP:6:73747426-73766256 |
| `gnomad_maxAF` | Maximum allele frequency among matching gnomAD SV(s), where matching is determined using the report’s SVTYPE-specific comparison logic | gnomAD SV | 0.00755 |
| `DGV_AF` | DGV allele frequency | DGV | 0.0955; 0.0955; 0.0955 |
| `DGV_NUM_ SAMPLES_TESTED` | Number of samples tested in the DGV study | DGV | 3425; 3425; 3425 |
| `DGV_SV` | Coordinates/type of matching DGV event(s), where matching is determined using the report’s SVTYPE-specific comparison logic | DGV | DEL:10:101634196-101635451; DEL:10:101634196-101635451; DEL:10:101634196-101635451 |
| `DGV_maxAF` | Maximum allele frequency among matching DGV event(s), where matching is determined using the report’s SVTYPE-specific comparison logic | DGV | 0.002 |
| `1000G_AF` | Allele frequency in the 1000 Genomes SV database | 1000 Genomes SV database | 0.15 |
| `1000G_AC` | Allele count in the 1000 Genomes SV database | 1000 Genomes SV database | 958 |
| `1000G_nhomalt` | Number of homozygous alternate individuals in the 1000 Genomes SV database | 1000 Genomes SV database | 254 |
| `1000G_SV` | Coordinates/type of matching 1000 Genomes SV(s), where matching is determined using the report’s SVTYPE-specific comparison logic | 1000 Genomes SV database | DEL:10:103293394-103294434 |
| `1000G_maxAF` | Maximum allele frequency among matching 1000 Genomes SV(s), where matching is determined using the report’s SVTYPE-specific comparison logic | 1000 Genomes SV database | 0.15 |
| `ExAC_delZ` | Positive delZ\_ExAC (Z score) from ExAC indicate gene intolerance to deletion | AnnotSV | 0.018797621 |
| `ExAC_dupZ` | Positive dupZ\_ExAC (Z score) from ExAC indicate gene intolerance to duplication | AnnotSV | \-1.900013288 |
| `ExAC_cnvZ` | Positive cnvZ\_ExAC (Z score) from ExAC indicate gene intolerance to CNV | AnnotSV | \-1.383777638 |
| `ExAC_synZ` | Positive synZ\_ExAC (Z score) from ExAC indicate gene intolerance to synonymous variation | AnnotSV | 0.310841961 |
| `ExAC_misZ` | Positive `misZ_ExAC` indicates gene intolerance to missense variation | AnnotSV | 1.612685911 |
| `ExAC_pLI` | Score computed by ExAC indicating the probability that a gene is intolerant to a loss of function variation (Nonsense, splice acceptor/donor variants due to SNV/indel). ExAC considers pLI\>=0.9 as an extremely LoF intolerant gene | AnnotSV | 0.99411 |
| `CytoBand` | Cytogenetic band annotation | AnnotSV | q24.32 |
| `RE_gene` | Name of the genes regulated by a regulatory element overlapped with the SV to annotate. When available, the regulated gene name is detailed with associated haploinsufficiency (HI), triplosensitivity (TS), Exomiser (EX) scores, OMIM and candidate genes. | AnnotSV | ACAT1 (morbid/RE=EA\_enhancer) |
| `TAD_coordinate` | Coordinates of the TAD whose boundaries overlapped with the annotated SV (boundaries included in the coordinates) | AnnotSV | 10:108760243-109880242 |
| `ENCODE_experiment` | ENCODE experiments supporting the TAD annotation | AnnotSV | ENCFF241RXN:ENCFF869GJT |
| `Repeat_type_left` | Repeat annotation around the left breakpoint (+/- 100 bp) | AnnotSV | L1PA10 |
| `Repeat_type_right` | Repeat annotation around the right breakpoint (+/- 100 bp) | AnnotSV | L1PA10 |
| `SegDup_left` | Segmental Duplication regions coordinates around the left SV breakpoint (+/- 100bp) | AnnotSV | 10:1239022-1240368 |
| `SegDup_right` | Segmental Duplication regions coordinates around the right SV breakpoint (+/- 100bp) | AnnotSV | 10:1239022-1240368 |
| `Adotto_tandem_repeat` | Overlap with Adotto repeat catalog | Adotto repeat catalog | ID=chr10\_1234947\_1240398\_ CGGTGTTTACTCCCCTCTGCCTCC; MOTIFS=CGGTGTTTA CTCCCCTCTGCCTCC; STRUC=(CGGTGTTTACTCCCCTCTGCCTCC)n |
| `ENCODE_blacklist_ characteristics_left` | ENCODE blacklist annotation near the left breakpoint | AnnotSV | High Signal Region |
| `ENCODE_blacklist_ characteristics_right` | ENCODE blacklist annotation near the right breakpoint | AnnotSV | High Signal Region |
