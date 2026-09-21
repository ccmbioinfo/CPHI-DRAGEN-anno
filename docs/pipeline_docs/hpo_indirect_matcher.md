# Indirect HPO term matching

## Overview

The following accepts patient HPO terms, `hp.json`, and `genes_to_phenotype.txt` and returns related, non-exact HPO terms and their directly annotated genes. This recovers genes missed by an exact HPO-ID join. For example, direct gene-to-HPO matching for a patient with *bilateral sensorineural hearing impairment* excludes *CDC14A* because that gene is annotated only to *sensorineural hearing impairment*. Indirect matching recovers *CDC14A*, resulting in the inclusion of this HPO term in the annotated report for variants in *CDC14A*. Exact HPO matches remain separate and receive score `1.0`.

`hpo-toolkit==0.8.1` loads the ontology. The `hp.json` and annotation file (`genes_to_phenotype.txt`) should come from the same HPO release. The thresholds were calibrated with HPO release `2026-06-23` and retested with release `2026-09-01`, which produced the same control performance. A future release with substantial ontology or annotation changes may require recalibration and review.

Important caveats:

- **The score is not a probability; it is a heuristic for ranking indirect HPO matches.** It combines ontology distance and shared-ancestor specificity with label and annotation similarity when needed. Higher scores indicate that there may be stronger evidence of relatedness based on the matching rules. The score is intended only as an aid for prioritizing matches.
- **Indirect matching is intentionally somewhat lenient, but not indiscriminate.** The aim is to retain meaningful related findings while rejecting terms whose genes provide essentially no phenotype-based reason for inspection. Contrasting terms on the same clinical axis may therefore be retained for a reviewer to assess.

The approach was inspired by:
- [Phrank](https://github.com/meng-ma-biomedical-AI/F29_Phrank)
- [LIRICAL's](https://github.com/TheJacksonLaboratory/LIRICAL/tree/master/lirical-core/src/main/java/org/monarchinitiative/lirical/core/likelihoodratio) phenotype-to-disease likelihood ratio
- Lin, D., 1998, July. An information-theoretic definition of similarity. In Icml (Vol. 98, No. 1998, pp. 296-304).

## 1. Building the HPO graph

HPO Toolkit loads `hp.json` as a DAG. Each current (non-obsolete) term supplies its ID, name, synonyms, alternative IDs, and broader `is_a` parents, and a term may have multiple parents. The direct parent links are retained and child links are derived so graph traversal can happen in either direction for a term.

Ontology loading redirects alternative and replaced obsolete IDs to their current terms.

Direct `gene_symbol` and `disease_id` annotations from `genes_to_phenotype.txt` are stored for each `hpo_id` and propagated upward to measure term breadth. The final matcher output contains only genes directly annotated to the respective candidate term; propagated gene annotations are used for scoring but are not included as matches.

The figures here reuse one ontology neighbourhood containing terms up to four parent-child `is_a` links (edges) away from the patient term. Nodes begin white, accepted terms turn green and show their score, and rejected terms turn red (Figure 3).

## 2. Selecting nearby candidate terms

First, potential candidates are discovered by searching through a shared ancestor and considering candidates within **three total `is_a` edges**. This permits nearby ancestors, descendants, siblings, and some children of siblings. The exact patient term is excluded.

For each candidate/patient term pair, the shortest path is selected, with ties favouring a direct ancestor/descendant path or the path with the higher shared-ancestor specificity (Section 4).

## 3. Direct ancestors and descendants

A candidate on the same direct lineage as the patient term (nodes highlighted green in Figure 3) is always considered related when it is within three edges. A fixed score is used to record its direction and distance:

| Candidate relative to patient | 1 edge | 2 edges | 3 edges |
|---|---:|---:|---:|
| broader annotation ancestor | 0.85 | 0.70 | 0.55 |
| narrower annotation descendant | 0.70 | 0.55 | 0.40 |

Since all terms that are direct ancestors or descendants of a patient term are included as matches, these fixed scores are arbitrary, being reduced by `0.15` as they get further from the patient term. Candidate descendants receive a score `0.15` lower than an ancestor at the same distance. A specific patient finding strongly supports a broader gene annotation, but a broad patient finding does not confirm which narrower subtype the patient has. The descendant can still be useful, so it is retained at a lower score.

## 4. Information content and ancestor specificity

#### Information content
The remaining cross-branch (sibling/child-of-sibling) decisions use information content (IC). IC is used to determine an HPO term's annotation breadth. The calculation runs separately for genes and diseases:

$$
IC_X(t) = \frac{\log\left(N_X / n_X(t)\right)}{\log(N_X)}
$$

- $t$ is the HPO term being measured.
- $X$ is the annotation type: either `gene` or `disease`.
- $N_X$ is the total number of distinct annotations of type $X$ in the reference `genes_to_phenotype.txt`. For gene IC, it is the total number of annotated genes; for disease IC, it is the total number of annotated disease IDs.
- $n_X(t)$ is the number of distinct annotations of type $X$ associated with term $t$ after propagation. It includes annotations made directly to $t$ and annotations made to any descendant of $t$.

If a term is associated with a large proportion of all genes or diseases ($n_X(t)$ is close to $N_X$), then IC approaches `0`. If it is associated with only a small annotation set its IC approaches `1`.

#### Annotation informativeness
A term's **annotation informativeness** is the mean of its gene and disease IC values:

$$
I(t) = \frac{IC_{gene}(t) + IC_{disease}(t)}{2}
$$

#### Descendant specificity
IC measures annotation breadth, but it does not directly measure how much of the HPO graph sits below a term. For this, descendant specificity is used:

$$
D(t) = 1 - \frac{\log\left(n_{desc}(t)+1\right)}{\log(N_{HPO})}
$$

- $t$ is the HPO term being measured.
- $n_{desc}(t)$ is the number of unique live HPO terms reachable by repeatedly following child links below $t$. This includes children, grandchildren, and every deeper descendant, but not $t$ itself. A descendant reached through multiple parents is counted once.
- $N_{HPO}$ is the total number of live terms in the loaded HPO graph.
- The `+1` makes a leaf term with no descendants valid.

A leaf term therefore has $D(t)=1$. As the number of descendants grows, $D(t)$ decreases.

#### Shared-ancestor specificity
Finally, **shared-ancestor specificity** combines gene IC, disease IC, and descendant specificity:

$$
S(t) = \frac{2I(t) + D(t)}{3}
     = \frac{IC_{gene}(t) + IC_{disease}(t) + D(t)}{3}
$$

$S(t)$ is therefore the mean of three signals. For a cross-branch comparison, $t$ is the ancestor shared by the patient and candidate terms. Values near `0` describe a broad shared ancestor; values near `1` describe a narrow ancestor.

**Before any cross-branch pair is considered further, its shared ancestor must have specificity of at least `0.20`.** This prevents two findings from matching only because both sit under a very broad parent such as *Abnormal skin morphology*.

## 5. Cross-branch terms on the same clinical axis

For sibling matching, first it is determined whether or not the terms describe the same clinical feature with different qualifiers. Every combination of their names and synonyms is normalized:

| Rule | Representative examples | Use |
|---|---|---|
| structural words | `a`, `and`, `in`, `of`, `the`, `with` | Remove wording without phenotype meaning. |
| clinical modifiers | `bilateral`, `congenital`, `developmental`, `juvenile`, `mild`, `postlingual`, `prelingual`, `progressive`, `unilateral` | Ignore onset, course, severity, laterality, and distribution when identifying the clinical core. |
| assay and specimen words | `activity`, `circulating`, `concentration`, `csf`, `level`, `serum` | Prevent common laboratory wording from defining the clinical core. |
| fused axis prefixes | `hyper`, `hypo`, `macro`, `micro`, `brady`, `tachy` | Provide a second, stricter check for directions fused into one word. |
| fused axis suffixes | `cytosis`, `cytopenia`, `penia`, `philia` | Apply the same opposing-direction check as above to fused blood-count directions. |
| generic finding words | `abnormal`, `anomaly`, `malformation`, `morphology`, `defect`, `disease` | Removed later when checking whether ordinary siblings share distinctive wording. |

The first three rules reduce *developmental cataract* and *juvenile cataract* to `cataract`, and prelingual and postlingual SNHI to `sensorineural hearing impairment`. Prefix and suffix handling similarly recognizes *hyperinsulinemia* and *hypoinsulinemia* as opposing directions of the same clinical feature because both reduce to the stem `insulinemia`. Directly opposing terms can still be retained as indirect matches because different variants or mechanisms involving the same gene may produce opposite phenotypic effects.

Same-axis terms may be direct siblings or three-edge children of siblings. The distance and shared-ancestor specificity are combined in one score:

$$
score_{axis}(p,c) = \left[1-0.15(d-1)\right]S(a)
$$

Here, $p$ is the patient term, $c$ is the candidate term, $a$ is their selected shared ancestor, $d$ is the total number of edges between the two terms, and $S(a)$ is the shared-ancestor specificity (Section 4).

## 6. Other cross-branch terms

Terms not on the same clinical axis must be **direct siblings**, each one edge below the shared ancestor. More distant non-axis relationships are rejected to prevent matching across branches that diverged for meaningful reasons.

#### Clinical similarity
For each normalized name or synonym pair, clinical similarity calculates the matched fraction of tokens in both directions and keeps the smaller value. The highest value across all pairs is $C(p,c)$. This prevents a short label contained in a longer one from scoring perfectly.

Tokens may match exactly or through a near-complete stem, such as *dystrophy* and *dystrophic*. Thus `hypotonia` versus `facial hypotonia` scores `0.50`: one direction matches completely `1.0`, but the other matches only one of two words `0.5`.

#### Distinctive-word similarity
Distinctive-word similarity $W(p,c)$ repeats this comparison after removing generic finding words and must be greater than `0` for the wording-supported gate. For example, removing `abnormal` and `morphology` exposes the difference between *Abnormal mandibular symphysis morphology* and *Abnormal oral cavity morphology*.

#### Semantic similarity
Semantic similarity measures the annotation information lost by moving from the two terms to their shared ancestor, separately for genes and diseases:

$$
M_X(p,c) = \frac{2IC_X(a)}{IC_X(p)+IC_X(c)}
$$

The final semantic similarity is their mean:

$$
M(p,c) = \frac{M_{gene}(p,c)+M_{disease}(p,c)}{2}
$$

This compares annotation breadth through IC.

### Two ways for a direct sibling to pass
A non-axis direct sibling can pass through either gate:

| Gate | Requirements | Purpose and example |
|---|---|---|
| wording-supported | clinical similarity `>= 0.50`, semantic similarity `>= 0.45`, and distinctive-word similarity `> 0` | Requires wording and annotation support; retains *hypertrophic* and *dilated cardiomyopathy*. |
| annotation-profile rescue | semantic similarity `>= 0.90` and shared-ancestor specificity `>= 0.35` | Allows different wording only with exceptional annotation support; can retain *recurrent boils* and *bacterial cellulitis*. |

Candidates passing either gate receive the full sibling score:

$$
\begin{aligned}
score_{sibling}(p,c) ={}& \left[1-0.15(d-1)\right]S(a) \\
&\times \sqrt{\frac{C(p,c)+M(p,c)}{2}} \\
&\times \sqrt{0.5+0.5C(p,c)}
\end{aligned}
$$

The first square root balances clinical and semantic evidence; the second adds weight to shared wording. Square roots keep multiplication from pushing most scores close to zero.

**The final sibling score must be at least `0.15`; otherwise the pair is rejected.**

In Figures 4 and 5, S2 represents a direct sibling that does not pass either the wording-supported or annotation-profile rescue gate, so it is rejected.

## 7. Deciding whether to add the candidate's genes

Related candidates may still be too broad to add useful information to the report:

| Direct genes on the candidate term | Decision                                              |
| ---------------------------------: | ----------------------------------------------------- |
|                               none | Skip because the term cannot add report rows.         |
|                            `1–100` | Include as a bounded, low-cost expansion.             |
|                        `101–2,000` | Include only when informativeness is at least `0.22`. |
|                  more than `2,000` | Do not indirectly expand.                             |

An indirect candidate is also omitted when it is within two edges of *Phenotypic abnormality* and its propagated branch contains more than 2,000 genes. Only that broad candidate is removed; other matches from the patient term and all exact patient-term matches are retained.

The final output contains one row per directly annotated `(gene, candidate HPO term)` pair. If multiple patient terms reach the same pair, only the highest relationship score is retained.

## 8. Validation

Validation of the thresholds and expansion rules used positive, negative, and borderline relationships; repeated random samples of 2,000 HPO terms to identify broad or unrelated expansions that were not represented in the controls; and manual review of retained and rejected examples.

The curated set contains 28 required relationships, 24 potentially useful relationships that do not drive tuning, and 12 unrelated relationships. Current validation retains 27/28 required pairs and 7/24 potentially useful pairs, while rejecting all 12 unrelated pairs.

The complete terms, tiers, and scores are in [`HPO_MATCHER_CONTROLS.tsv`](HPO_MATCHER_CONTROLS.tsv). Controls were selected by exploring parts of the HPO ontology around SNHI, cataract, and rod-cone dystrophy. Additional controls were added later by randomly sampling HPO terms and reviewing the resulting matches to decide whether a match should be retained or rejected. A score of `0.00` means the pair was rejected; a positive score does not guarantee that its genes pass the reportability checks in Section 7.
