// Helper functions loaded into `slivar expr` with `--js`.
// Slivar evaluates the expressions in wrapper.py as JavaScript.

// Return true when a Slivar INFO value has usable content. The wrapper passes
// numeric frequency fields, string ClinVar/status fields, or multi-valued INFO arrays.
function present(x) {
  if (x === undefined || x === null) return false;

  if (Array.isArray(x)) {
    // INFO fields maybe arrays; accept the field if
    // any array element is present.
    for (var i = 0; i < x.length; i++) {
      if (present(x[i])) return true;
    }
    return false;
  }

  // Numeric INFO values are present when they are real numbers, not NaN.
  if (typeof x === "number") return isFinite(x);

  // For string-like values, treat the common VCF missing tokens as absent.
  var s = String(x);
  return s !== "" && s !== "." && s !== "None";
}

// Return the first alternate-allele depth from a sample's FORMAT/AD value.
function alt_depth(sample) {
  // AD is expected to look like [ref_depth, first_alt_depth] after decomposition.
  if (!("AD" in sample) || !Array.isArray(sample.AD) || sample.AD.length < 2) {
    return null;
  }

  // After decomposition the first ALT allele is the alt allele.
  var depth = sample.AD[1];
  if (depth === undefined || depth === null || !isFinite(depth) || depth < 0) {
    return null;
  }
  return depth;
}

// Per-sample AD predicate used by wrapper.py's `ad_pass` sample expression.
function passes_alt_depth(sample, threshold) {
  var depth = alt_depth(sample);
  return depth !== null && depth >= threshold;
}

// Return true when no sample in has usable first-ALT AD.
function all_alt_depths_missing() {
  // `$S` is provided by Slivar and contains all samples keyed by sample ID.
  for (var sample_id in $S) {
    if (alt_depth($S[sample_id]) !== null) return false;
  }
  return true;
}

// Case-insensitive string match used by common ClinVar rescue in report branches.
function contains_pathogenic(x) {
  if (!present(x)) return false;

  if (Array.isArray(x)) {
    // Pass if any ClinVar annotation contains "pathogenic".
    for (var i = 0; i < x.length; i++) {
      if (String(x[i]).toLowerCase().indexOf("pathogenic") !== -1) return true;
    }
    return false;
  }

  return String(x).toLowerCase().indexOf("pathogenic") !== -1;
}

// Common CH rescue only applies when the ClinVar field text exactly matches one of these labels.
function is_ch_common_clinvar_rescue(x) {
  if (!present(x)) return false;

  if (Array.isArray(x)) {
    // Multi-valued ClinVar annotations pass if any element is one of the CH rescue labels.
    for (var i = 0; i < x.length; i++) {
      if (is_ch_common_clinvar_rescue(x[i])) return true;
    }
    return false;
  }

  // This is stricter than contains_pathogenic(): substring matches such as
  // "Pathogenic/Likely_pathogenic" are not rescued for common CH candidates.
  var s = String(x);
  return s === "Pathogenic" ||
    s === "Likely_pathogenic" ||
    s === "Conflicting_interpretations_of_pathogenicity";
}
