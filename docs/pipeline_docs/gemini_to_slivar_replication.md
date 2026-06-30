# Replicating CRE/GEMINI small variant reporting with slivar

## Summary

This document summarizes the slivar/Python replacement for the CRE/GEMINI small
variant reporting workflow in `CPHI-DRAGEN-anno`. The replacement currently targets
these report modes:

- `coding`
- `wgs-high-impact`
- `wgs`-style `denovo`, `panel`, and `panel-flank`

The comparisons and counts in this document are based on testing one sample set.
A demo of the final reports can be found here:

CRE:
- `/hpf/largeprojects/ccmbio/students/rkhan/cphi_dragen_runs/260623_slivar_demo/reports`

slivar:
- `/hpf/largeprojects/ccmbio/students/rkhan/cphi_dragen_runs/260623_slivar_demo/reports_slivar`

The CRE workflow loads an annotated VCF into GEMINI with `vcf2db`. During
loading, `vcf2db.py` uses `geneimpacts` to choose one top transcript/effect for
the GEMINI `variants` table. Most CRE filtering then queries fields from that
selected `variants` table row with GEMINI SQL filters. The R report script
`workflow/scripts/cre/cre.vcf2db.R` applies later report-level logic, including
the final AD >= 3 filter, `high-impact` score filtering, `high-impact` cleanup
requiring a selected Gene annotation, depth/genotype cleanup, external
annotation joins, and final CSV formatting.

The slivar replacement avoids the GEMINI database. `slivar expr` mostly
replicates the GEMINI query branches directly on the VCF.
`workflow/scripts/slivar/postfilter.py` then applies additional filters that
are not covered cleanly by `slivar expr`, including sample depth/AD checks, the
final AD >= 3 gate, `high-impact` score filtering, and the `high-impact`
primary-gene-present gate. `workflow/scripts/slivar/build_report.py` formats
the final CSV and joins the same local CRE annotation tables. Compound-het input
TSVs are produced by `workflow/scripts/slivar/build_ch_tsv.py` and then passed
to the same configured downstream compound-het annotation scripts used by the
CRE path.

Based on matching variant keys with `Position + Ref + Alt`, there is very high
parity in the final reported variants for each mode. `denovo`, `panel`, and
`panel-flank` have no one-sided final variants in the current demo comparison.
`coding` has 52 slivar-only keys. `high-impact` has 73 CRE/GEMINI-only keys
and 1598 slivar-only keys.

The larger and more impactful differences come from transcript and consequence
selection. In the CRE/GEMINI workflow, `vcf2db` and `geneimpacts` choose a
single primary transcript/effect during database loading, and those selected
values in the GEMINI `variants` table are then used for most downstream
filtering, reporting, and annotation joins. In the slivar workflow, filtering
can evaluate the full VEP `CSQ` annotation set on the VCF; only later does a
Python script choose a primary CSQ for the final report. Additionally, the
slivar report keeps `_all` columns such as `Gene_all`, `Ensembl_gene_id_all`,
and `Variation_all`. These preserve the broader VEP context, but downstream
joins for compound hets and ACMG secondary findings still depend mainly on the
selected primary gene or Ensembl gene ID.

The consequence ranking also affects `coding`-mode inclusion, compound-het HIGH-MED
versus LOW placement, and which consequence is displayed in the report. The
local slivar consequence order in `workflow/scripts/slivar/default-order.txt`
is close to, but not identical to, the `geneimpacts` severity mapping in
`/Users/rohan.khan/Desktop/slivar_attempt/effect.py`.

## CRE/GEMINI Workflow

The CRE path starts from a PASS-filtered, annotated VCF and creates a GEMINI
database with `vcf2db`. `vcf2db.py` parses the VEP `CSQ` annotations with `geneimpacts` and writes two
main annotation tables:

- `variants`: one row per VCF variant, with one primary CSQ/transcript selected
  for that variant
- `variant_impacts`: one row per parsed CSQ/transcript/effect

The `variants` table is the primary table for CRE small-variant reporting.
`geneimpacts` chooses the primary CSQ/transcript for each variant based on its
own consequence severity ranking, coding/biotype handling, pseudogene
deprioritization, and related transcript tie-breaks. It also maps the selected
consequence/variation term to the custom GEMINI `impact` and `impact_severity`
values. Those severity values help determine which transcript is selected and
are later used directly in report filtering. The transcript-selection details
are described more in the transcript selection section below.

The selected primary effect/transcript is represented in `variants` table
fields used later by CRE workflow filtering and reporting, including:

- `gene`
- `ensembl_gene_id`
- `transcript`
- `impact`
- `impact_severity`
- `is_coding`

GEMINI report queries in `workflow/scripts/cre/cre.gemini2txt.vcf2db.sh`
evaluate those selected `variants` table values. For example, `coding`
filtering uses the selected `impact_severity <> 'LOW'`, while non-`coding`
modes set `severity_filter=ALL` and do not restrict the main branch by
HIGH/MED. Queries against `impact_severity`, `impact`, `gene`, and `transcript`
are therefore already downstream of `geneimpacts` transcript/effect selection.

`workflow/scripts/cre/cre.vcf2db.R` then builds the final CSV. At this stage,
`variant_impacts` is used mainly to populate report context such as the `Info`
field, while the main report rows still come from the GEMINI `variants` export.
Important report-level behavior in the R script includes:

- requiring at least one sample with final AD >= 3 before writing the report
- applying the `high-impact` score gate for `wgs.high.impact`
- removing `high-impact` rows whose selected `Variation` is `intergenic_variant`
- removing `high-impact` rows with missing selected `Gene`
- joining CRE data tables

The GEMINI database also feeds compound-het input preparation. The HIGH-MED and
LOW sequence-variant TSVs are split from the `variants` table using
`impact_severity`; the selected primary transcript/effect determines the
`Ensembl_gene_id` attached to each sequence variant. That Ensembl gene ID is
then used downstream when CH status and variant types are merged back into the
small-variant reports.

## slivar Replacement Workflow

The slivar path keeps the VCF as the source of truth and reproduces the
GEMINI/R behavior without creating a GEMINI database. The replacement has three
main stages.

1. `workflow/wrappers/slivar/wrapper.py` runs `slivar expr` to create broad
   report branches. This is the closest replacement for the GEMINI query step:
   it applies PASS-only selection, currently rejects star ALT records, and
   separates rare, rare-ClinVar, common pathogenic ClinVar, and
   `coding`-impactful candidates. The wrapper passes
   `workflow/scripts/slivar/default-order.txt` to slivar; terms above the
   `IMPACTFUL_CUTOFF` sentinel are used to set `INFO.impactful`, which is the
   slivar-side approximation of the geneimpacts HIGH/MED cutoff for `coding`
   filtering.

2. `workflow/scripts/slivar/postfilter.py` rechecks branch logic in Python,
   applies sample depth/AD thresholds that are awkward to express in
   `slivar expr`, deduplicates branch outputs, and applies the final AD >= 3
   report gate. For `wgs-high-impact`, it also applies the numeric `high-impact`
   score gate and removes records that are all intergenic or whose
   postfilter-selected primary CSQ lacks a gene symbol.

3. `workflow/scripts/slivar/build_report.py` groups records by variant key,
   parses all VEP `CSQ` records, chooses a primary CSQ for display and joins,
   formats CRE-like report columns, computes SpliceAI and `Noncoding_path_pred`
   with the same thresholds as the R script, joins the same CRE annotation
   resources, and writes `_all` columns that preserve the full CSQ context.

For compound hets, `workflow/scripts/slivar/build_ch_tsv.py` creates the
HIGH-MED and LOW sequence-variant TSVs expected by the existing compound-het
annotation code. It uses the slivar consequence order and `IMPACTFUL_CUTOFF` to
split variants between HIGH-MED and LOW, and uses slivar `coding`-mode primary
CSQ selection to choose the `Ensembl_gene_id` written to the TSV.

## Filtering Logic

The GEMINI databases are queried through different branches depending on the
report mode. `slivar expr` is used to replicate the GEMINI filtering where
possible, and anything that cannot be handled cleanly in `slivar expr` is done
in Python. Some filtering also happens in R in the CRE/GEMINI workflow; in the
slivar workflow, that logic is consolidated with the other slivar filtering
instead of being left in the report-generation step.

### `coding`

CRE/GEMINI logic:

- Main branch: `gnomad_fafmax_faf95_max <= 0.01` and HIGH/MED effect only.
  The main branch is also genotype-filtered to variants where any sample has
  `gt_alt_depths >= 3`, or all sample alt depths are `-1`.
- Rare ClinVar branch:
  `gnomad_fafmax_faf95_max <= 0.01`, ClinVar text present, and AD >= 1 or all
  sample alt depths are `-1`.
- Common pathogenic ClinVar branch:
  `gnomad_fafmax_faf95_max > 0.01`, ClinVar status not
  `no_assertion_criteria_provided`, and `lower(Clinvar) like '%pathogenic%'`.
- Although the rare ClinVar branch initially allows AD >= 1, the final report
  still requires at least one sample with AD >= 3 after `cre.vcf2db.R`
  merges/fixes depth values.

slivar logic:

- Rare impactful branch: slivar first selects PASS, non-star variants with
  `INFO.impactful` and FAF <= 0.01 or missing FAF. `INFO.impactful` is based on
  `workflow/scripts/slivar/default-order.txt`. In `postfilter.py`, the same
  branch is checked again for FAF <= 0.01 or missing FAF, `impact_severity` or
  the slivar-produced `INFO.impactful` flag, and AD >= 3 or all AD missing.
- Rare ClinVar branches: split into two branches based on missing FAF or
  FAF <= 0.01. Both require ClinVar text, and `postfilter.py` applies AD >= 1
  or all sample alt depths are `-1`. The two slivar branches are needed because
  the GEMINI query behavior treats missing values differently from the direct
  slivar expression. In `coding` mode, the allow-missing-FAF `slivar expr` branch
  also requires `INFO.genic`. This flag is set by slivar from the
  `GENIC_CUTOFF` sentinel in `default-order.txt`; it is less stringent than
  `INFO.impactful` and is used here to keep the missing-FAF ClinVar rescue in a
  gene/transcript annotation context rather than pulling arbitrary intergenic
  ClinVar records into the `coding` branch.
- Common pathogenic ClinVar branch: FAF > 0.01, ClinVar status present and not
  `no_assertion_criteria_provided`, and ClinVar text contains `pathogenic`.
- Final output requires AD >= 3.

### `wgs-high-impact`

CRE/GEMINI logic:

- Main branch: all severities are allowed, `gnomad_fafmax_faf95_max <= 0.001`,
  and any sample must have `gt_alt_depths >= 3`, or all sample alt depths must
  be `-1`.
- Rare ClinVar branch:
  `gnomad_fafmax_faf95_max <= 0.01`, ClinVar text present, and AD >= 1 or all
  sample alt depths are `-1`.
- Common pathogenic ClinVar branch:
  `gnomad_fafmax_faf95_max > 0.01`, ClinVar status not
  `no_assertion_criteria_provided`, and `lower(Clinvar) like '%pathogenic%'`.
- After GEMINI export, `cre.vcf2db.R` keeps rows with SpliceAI score >= 0.2,
  missing CADD, CADD >= 10, or promoterAI absolute score >= 0.1. It then
  removes rows whose selected `Variation` is `intergenic_variant` and rows with
  missing selected `Gene`.
- Although the rare ClinVar branch initially allows AD >= 1, the final report
  still requires at least one sample with AD >= 3 after `cre.vcf2db.R`
  merges/fixes depth values.

slivar logic:

- Rare main branch: slivar first selects PASS, non-star variants with
  FAF <= 0.001 or missing FAF. In `postfilter.py`, the same branch is checked
  again for FAF <= 0.001 or missing FAF and AD >= 3 or all AD missing.
- Rare ClinVar branches: split into two branches based on missing FAF or
  FAF <= 0.01. Both require ClinVar text, and `postfilter.py` applies AD >= 1
  or all sample alt depths are `-1`.
- Common pathogenic ClinVar branch: FAF > 0.01, ClinVar status present and not
  `no_assertion_criteria_provided`, and ClinVar text contains `pathogenic`.
- After any branch passes, `postfilter.py` applies the intended `high-impact`
  score gate, removes records where all CSQ annotations are
  `intergenic_variant`, and removes records whose postfilter-selected primary
  CSQ has no gene symbol.
- Final output requires AD >= 3.

For `high-impact` intergenic cleanup, CRE and slivar
checks are not exactly identical. CRE removes a row when the
geneimpacts-selected GEMINI `Variation` is `intergenic_variant`; slivar removes
a record only when every CSQ annotation is `intergenic_variant`. This has not
been the main source of current one-sided `high-impact` keys, but it is a logic distinction.

### `wgs`-Style `denovo`, `panel`, and `panel-flank`

These reports are subset before report filtering:

- `denovo` starts from `filtered/{family}.denovo.vcf.gz`, produced with
  `FORMAT/DN == 'DeNovo'`.
- `panel` starts from the base filtered VCF intersected with the HPO-derived
  panel BED.
- `panel-flank` starts from the base filtered VCF intersected with the same
  panel BED after flank expansion.

After that upstream subsetting, the workflow runs these reports with `wgs`-style
logic rather than the older CRE pedigree-based `type=denovo` GEMINI query.

CRE/GEMINI logic:

- Main branch: all severities are allowed, `gnomad_fafmax_faf95_max <= 0.01`,
  and any sample must have `gt_alt_depths >= 3`, or all sample alt depths must
  be `-1`.
- Rare ClinVar branch:
  `gnomad_fafmax_faf95_max <= 0.01`, ClinVar text present, and AD >= 1 or all
  sample alt depths are `-1`.
- Common pathogenic ClinVar branch:
  `gnomad_fafmax_faf95_max > 0.01`, ClinVar status not
  `no_assertion_criteria_provided`, and `lower(Clinvar) like '%pathogenic%'`.
- Final report still requires AD >= 3 after depth values are merged/fixed.

slivar logic:

- Rare main branch: slivar first selects PASS, non-star variants with
  FAF <= 0.01 or missing FAF. In `postfilter.py`, the same branch is checked
  again for FAF <= 0.01 or missing FAF and AD >= 3 or all AD missing.
- Rare ClinVar branches: split into two branches based on missing FAF or
  FAF <= 0.01. Both require ClinVar text, and `postfilter.py` applies AD >= 1
  or all sample alt depths are `-1`.
- Common pathogenic ClinVar branch: FAF > 0.01, ClinVar status present and not
  `no_assertion_criteria_provided`, and ClinVar text contains `pathogenic`.
- Final output requires AD >= 3.

Follow-up from this section:

- Fix slivar star ALT handling to include CRE's
  retained `ALT == "*"` records.
- Decide whether slivar `high-impact` intergenic cleanup should match CRE's
  selected-`Variation` behavior instead of using the all-CSQ-intergenic rule.
- Review which `postfilter.py` validation checks are still needed now that
  `slivar expr` already applies broad branch filtering. Some duplicate checks,
  such as PASS status and possibly the `coding` `INFO.genic` branch guard, may no
  longer be necessary.
- Revisit the rare ClinVar FAF handling. The current implementation uses
  separate present-FAF and allow-missing-FAF branches; check whether this split
  is still needed for all modes, and why `coding` can appear to work with a single
  branch in some comparisons.
- Looking at some `slivar` documentation, there might be a way to filter by depth
  and potentially other criteria so it does not need to be done in Python. Need
  to look into further.

## Final Variant-Key Differences

The CRE and slivar final reports were compared by variant keys based on
`Position + Ref + Alt`.

The `wgs`-mode reports, `denovo`, `panel`, and `panel-flank`, match exactly by
final variant key and have no CRE-only or slivar-only keys.

`coding` and `high-impact` have the following one-sided keys:

| Report | Shared keys | CRE/GEMINI-only | slivar-only |
| --- | ---: | ---: | ---: |
| `coding` | 3593 | 0 | 52 |
| `high-impact` | 14410 | 73 | 1598 |

For these one-sided variants, checks were done to see how `geneimpacts` would
handle each variant and which primary transcript/effect it would choose. For
all 52 slivar-only `coding` variants, the `geneimpacts`-selected transcript has
a LOW effect and no ClinVar annotation. CRE therefore does not return them from
the `coding` queries. The slivar-only `coding` set breaks down as:

- 30 variants where a pseudogene annotation has a slivar-impactful consequence,
  often splice-region-related, but `geneimpacts` deprioritizes the pseudogene
  and selects a LOW transcript.
- 10 variants with the same pseudogene-trigger pattern, but the final slivar
  report displays a protein-coding transcript because filtering and display
  selection happen in different steps.
- 7 variants where the slivar-impactful transcript has `protein_coding_LoF`;
  `geneimpacts` only treats exactly `protein_coding` as ordinary
  protein-coding for its biotype logic.
- 3 variants where the impactful consequence is on a
  `nonsense_mediated_decay` transcript, which `geneimpacts` does not treat as
  ordinary `protein_coding`.
- 2 ZNF302 indels where slivar treats `incomplete_terminal_codon_variant` as
  impactful, while `geneimpacts` classifies the selected effect as LOW.

For `high-impact`, 72 of the 73 CRE/GEMINI-only variants appear to fail the
intended `high-impact` score gate. Manual checks through the R script in RStudio
suggest this is due to a CRE R handling bug where some score fields are not
consistently normalized as numeric before filtering. For example,
`chr10:67837580:A:G` is present in the CRE/GEMINI `high-impact` report with
CADD 7.938, SpliceAI 0.0, and promoterAI 0.0005, which should not pass the
intended SpliceAI >= 0.2, missing CADD, CADD >= 10, or promoterAI >= 0.1 gate.
The only other CRE/GEMINI-only key is `chr7:104833597:A:*`, a star ALT record
excluded by slivar branch selection.

For the 1598 slivar-only `high-impact` variants, 1597 are explained by
CRE/geneimpacts selecting a top transcript/effect with a blank `Gene`. CRE then
removes those rows in `cre.vcf2db.R`. Slivar's transcript selection can choose a
different primary CSQ with a gene symbol, allowing the variants to remain. The
remaining slivar-only `high-impact` key is `chr14:18631541:C:CA`; from checks
done so far, it should have passed `geneimpacts` and GEMINI querying, but its
CRE absence is still unexplained.

Follow-up from this section:

- Investigate why `chr14:18631541:C:CA` is absent from the CRE/GEMINI
  `high-impact` report despite appearing to satisfy the expected geneimpacts and
  GEMINI query checks.
- Look for a better way to incorporate `Gene` presence into transcript
  selection without explicitly biasing toward any transcript with a gene
  annotation. This may also help with consolidating transcript ordering logic across
  `coding`, `high-impact`, and `wgs`.

## Transcript and Consequence Selection

Transcript selection is the main parity risk because CRE/GEMINI chooses a
primary effect during database loading, while slivar chooses it later from the
CSQ records that survived VCF filtering.

### CRE/GEMINI Selection

`vcf2db.py` uses `geneimpacts.Effect.top_severity()` to select the top effect.
For VEP `CSQ`, each transcript/effect is parsed into a `geneimpacts.VEP`
object. The selected effect is then written to the GEMINI `variants` table, and
CRE queries/report fields generally use those selected table values.

`geneimpacts` does this as a pairwise sort in `effect.py`; the higher-ranked
effect sorts later and is then selected. The comparison logic is:

- If one transcript annotation is from a pseudogene biotype and the other is not, 
  `geneimpacts` ranks the non-pseudogene annotation higher, before checking 
  consequence severity.
- If one effect is `coding`, prefer the coding effect unless either side is
  splicing. In splice comparisons, the code skips this coding-status shortcut
  and continues to severity.
- Compare custom severity next: HIGH > MED > LOW > UNKNOWN. Severity is the
  maximum mapped severity across the annotation's consequence terms, using the
  `IMPACT_SEVERITY` map. Separately, the displayed `top_consequence` is the
  highest-ranked consequence term within that annotation.
- If severity ties, prefer `protein_coding` biotype over non-`protein_coding`.
- If neither side was resolved by `protein_coding`, prefer
  `processed_transcript` over other biotypes.
- If still tied, compare parsed SIFT scores. The code treats an effect with a
  lower SIFT score as lower-ranked than an effect with a higher SIFT score;
  missing values use a very large fallback.
- Then compare parsed PolyPhen scores. The code treats an effect with a higher
  PolyPhen score as lower-ranked than an effect with a lower PolyPhen score;
  missing values use a very small fallback.
- Finally, break remaining ties by the maximum `IMPACT_SEVERITY_ORDER` rank
  among each effect's consequence terms.

`top_severity()` sorts effects with the highest-ranked effect last. If multiple
effects remain tied at the top, it can return a list; `vcf2db.py` then keeps the
first entry from that returned list for the `variants` table. There is no MANE
or gene-symbol-presence tie-break in the CRE/GEMINI selection path.

### slivar Selection

Final slivar report display and annotation joins use
`choose_primary_csq()` from `workflow/scripts/slivar/shared.py`.
`postfilter.py` has a separate local primary-CSQ selector for `high-impact`
missing-gene filtering. Its priorities are similar but not identical to the
shared report selector, so `high-impact` missing-gene filtering and final
display are not guaranteed to choose the same CSQ in every tie case.

For `coding`, slivar prioritizes:

- MANE Select / MANE Plus Clinical: prefer clinically curated MANE transcripts
  when present based on VEP annotation.
- Consequence order: choose the more severe consequence according to
  `default-order.txt`.
- VEP `IMPACT`: use VEP's HIGH/MODERATE/LOW/MODIFIER label as an additional
  severity tie-break.
- Canonical status: prefer VEP canonical transcripts when earlier fields tie.
- Biotype: prefer protein-coding biotype.
- Feature ID and original CSQ index: final tie-breaks are based on the feature
  ID and original CSQ index.

For `wgs-high-impact` and `wgs`, slivar prioritizes:

- MANE Select / MANE Plus Clinical: prefer clinically curated MANE transcripts
  when present based on VEP annotation.
- Biotype: prefer protein-coding transcripts earlier in the sort, and push
  pseudogene biotypes later than other non-protein-coding annotations.
- Consequence order: choose the more severe consequence within the preferred
  biotype group according to `default-order.txt`.
- Gene-symbol presence: prefer annotations with a usable gene symbol.
- Canonical status: prefer VEP canonical transcripts when earlier fields tie.
- Feature ID and original CSQ index: final tie-breaks are based on the feature
  ID and original CSQ index.

Because MANE is first in both modes, slivar can select a MANE transcript even
when another transcript has a more severe consequence. This can make the final
report display and downstream joins use a different transcript/gene than the
annotation that caused a variant to pass filtering, and it is one reason slivar
can disagree with the `geneimpacts`-selected GEMINI row even when both
workflows keep the same genomic variant.

Follow-up from this section:

- Make `postfilter.py` use the shared `choose_primary_csq()` logic, or otherwise
  centralize primary-CSQ selection, so `high-impact` filtering and final report
  display do not silently diverge.
- Revisit whether MANE should remain the first transcript-selection key or move
  behind consequence severity for closer `geneimpacts` parity.

## slivar Default Order Versus geneimpacts

`workflow/scripts/slivar/default-order.txt` is used in two places:

- by slivar to set `INFO.impactful` for consequences above `IMPACTFUL_CUTOFF`
- by Python report/CH code to rank CSQ consequences

This approximates, but does not exactly reproduce, the `geneimpacts` HIGH/MED
versus LOW classification.

Terms that are HIGH or MED in `geneimpacts` but below the current slivar
`IMPACTFUL_CUTOFF` or absent from the slivar order file include:

| Consequence term | geneimpacts severity | Current slivar order status |
| --- | --- | --- |
| `initiator_codon` | MED | Below `IMPACTFUL_CUTOFF` |
| `splice_donor_5th_base` | MED | Below `IMPACTFUL_CUTOFF` |
| `splice_polypyrimidine_tract` | MED | Below `IMPACTFUL_CUTOFF` |
| `splice_donor_region` | MED | Below `IMPACTFUL_CUTOFF` |
| `exon_region` | MED | Below `IMPACTFUL_CUTOFF` |
| `feature_fusion` | HIGH | Below `IMPACTFUL_CUTOFF` |
| `regulatory_region_ablation` | MED | Below `IMPACTFUL_CUTOFF` |
| `5_prime_UTR_truncation` | MED | Absent as written; slivar has lowercase `5_prime_utr_truncation` |
| `chromosome_number_variation` | HIGH | Absent |
| `duplication` | MED | Absent |
| `feature_ablation` | HIGH | Absent |
| `inversion` | MED | Absent |

There are also terms where slivar can be more permissive than `geneimpacts`.
For example, `incomplete_terminal_codon` is above the slivar cutoff, while
`geneimpacts` classifies `incomplete_terminal_codon_variant` as LOW. These
mapping differences affect `coding` inclusion and compound-het HIGH-MED versus
LOW placement.

Follow-up from this section:

- Decide whether `default-order.txt` should be changed to match
  `geneimpacts` HIGH/MED/LOW more closely or kept intentionally different.
- Recheck terms affected by case or suffix normalization, especially
  `5_prime_UTR_truncation` versus `5_prime_utr_truncation`.

## Compound-Het Replication

The logic behind compound het reporting does not change between slivar and CRE 
workflows. However, compound-het parity is affected by two main differences:

- whether sequence variants are placed in the HIGH-MED or LOW TSV
- which primary `Ensembl_gene_id` is selected before CH annotations are merged
  back into reports

The CRE/GEMINI workflow creates HIGH-MED and LOW sequence-variant TSVs from the
GEMINI `variants` table. The split is based on `impact_severity`, using rare
`gnomad_af_grpmax <= 0.01` variants with AD >= 3, rare ClinVar variants with
AD >= 1, and common ClinVar rescue variants. The `Ensembl_gene_id` comes from
the `vcf2db`/`geneimpacts` selected `variants` table values.

The slivar workflow uses `build_ch_tsv.py` to create the same TSV shape from
the annotated VCF. It mirrors the rare, rare-ClinVar, and common-ClinVar rescue
branches with `gnomad_af_grpmax`, splits HIGH-MED versus LOW using
`default-order.txt` and `IMPACTFUL_CUTOFF`, and uses slivar `coding`-mode
primary-CSQ selection for `Ensembl_gene_id`.

Both paths then pass the HIGH-MED and LOW TSVs, SV report, CNV report, pedigree,
gene-ID mapping tables, and sample order to the same downstream CH annotation
script `annotate_compound_hets.py`. 

When HIGH-MED and LOW are combined, both workflows account for the same genomic
variant keys. The differences are variants moving between groups:

| CH TSV group | GEMINI unique keys | slivar unique keys | Shared keys | GEMINI-only in group | slivar-only in group |
| --- | ---: | ---: | ---: | ---: | ---: |
| HIGH-MED | 3252 | 2842 | 2805 | 447 | 37 |
| LOW | 562238 | 562648 | 562201 | 37 | 447 |
| Combined HIGH-MED + LOW | 564594 | 564594 | 564594 | 0 | 0 |

Variants moved into LOW may no longer survive the downstream LOW-impact rescue
filter. That filter excludes intergenic sequence variants but keeps genic 
upstream/downstream annotations, and it keeps LOW variants only when they have 
ClinVar evidence or satisfy SpliceAI >= 0.2, CADD >= 10, promoterAI >= 0.1, or 
the indel-with-missing-score rule. 

Final `CH_status` and `CH_variant_types` are merged back into the `coding`,
`high-impact`, `denovo`, `panel`, and `panel-flank` reports by
`Ensembl_gene_id`. If CRE/GEMINI and slivar attach the same genomic variant to
different Ensembl gene IDs, CH annotations can differ even when the variant key
matches.

The table below compares shared final-report rows by `Position + Ref + Alt`.
The final column counts rows where at least one CH field differs and the
primary `Ensembl_gene_id` also differs.

| Report | Shared keys | Different `Ensembl_gene_id` | `CH_status` mismatches | `CH_variant_types` mismatches | Either CH field mismatch | Either CH mismatch with different `Ensembl_gene_id` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `coding` | 3593 | 494 | 534 | 446 | 578 | 308 |
| `high-impact` | 14410 | 709 | 1619 | 1163 | 1811 | 473 |
| `denovo` | 370 | 21 | 37 | 32 | 38 | 12 |
| `panel` | 17648 | 690 | 1104 | 916 | 1255 | 450 |
| `panel-flank` | 56502 | 3738 | 3551 | 3069 | 3968 | 1796 |

Follow-up from this section:

- Confirming cleaned up CH TSV handling for partial/full-missing genotypes, depth/AD
  formatting, phase-set formatting, genotype-quality missing values, and shared
  sample-name normalization between report generation and CH TSV generation.

## ACMG Secondary Findings

The CRE and slivar workflows use the same ACMG SF scripts:

- `workflow/scripts/add_acmg_sf_columns.py`
- `workflow/scripts/create_acmg_sf_report.py`

`add_acmg_sf_columns.py` searches for exact gene-name matches against the ACMG
SF gene list using the `Gene` column for small-variant reports. It does not 
currently inspect slivar's `Gene_all` column. Therefore, transcript-selection 
differences can change whether a variant is marked as an ACMG SF hit.

In the demo comparison, the CRE/GEMINI SF report contains 43 variants. By
genomic variant key, 42 are accounted for in the slivar SF report. The only
CRE-only genomic variant is `chr11:2628075 C>T`: CRE reports primary gene
`KCNQ1`, while slivar reports primary gene `KCNQ1OT1`, so the slivar row does
not match the ACMG SF gene list.

Follow-up from this section:

- Decide whether ACMG SF matching should also inspect slivar `Gene_all` so a
  secondary gene annotation can rescue cases where the primary gene differs.

## Known Remaining Bugs and Gaps

Most core report filtering has been replicated. The detailed follow-up items
are listed at the end of the relevant sections above. At a high level, the
remaining work clusters into:

- transcript/consequence-selection parity, including centralizing the primary
  CSQ selector and deciding whether MANE-first behavior should remain
- consequence-order parity between `default-order.txt` and `geneimpacts`
- filtering edge cases, especially star ALT handling and `high-impact`
  intergenic cleanup
- report and CH TSV formatting edge cases around missing values, genotypes,
  depth/AD, phase sets, and sample-name normalization
- ACMG SF gene matching policy for primary versus `_all` gene annotations
- code cleanup to reduce repeated validation between `slivar expr`,
  `postfilter.py`, and report generation
