function present(x) {
  if (x === undefined || x === null) return false;

  if (Array.isArray(x)) {
    for (var i = 0; i < x.length; i++) {
      if (present(x[i])) return true;
    }
    return false;
  }

  if (typeof x === "number") return isFinite(x);

  var s = String(x);
  return s !== "" && s !== "." && s !== "None";
}

function alt_depth(sample) {
  if (!("AD" in sample) || !Array.isArray(sample.AD) || sample.AD.length < 2) {
    return null;
  }

  var depth = sample.AD[1];
  if (depth === undefined || depth === null || !isFinite(depth) || depth < 0) {
    return null;
  }
  return depth;
}

function passes_alt_depth(sample, threshold) {
  var depth = alt_depth(sample);
  return depth !== null && depth >= threshold;
}

function all_alt_depths_missing() {
  for (var sample_id in $S) {
    if (alt_depth($S[sample_id]) !== null) return false;
  }
  return true;
}

function contains_pathogenic(x) {
  if (!present(x)) return false;

  if (Array.isArray(x)) {
    for (var i = 0; i < x.length; i++) {
      if (String(x[i]).toLowerCase().indexOf("pathogenic") !== -1) return true;
    }
    return false;
  }

  return String(x).toLowerCase().indexOf("pathogenic") !== -1;
}

function is_ch_common_clinvar_rescue(x) {
  if (!present(x)) return false;

  if (Array.isArray(x)) {
    for (var i = 0; i < x.length; i++) {
      if (is_ch_common_clinvar_rescue(x[i])) return true;
    }
    return false;
  }

  var s = String(x);
  return s === "Pathogenic" ||
    s === "Likely_pathogenic" ||
    s === "Conflicting_interpretations_of_pathogenicity";
}
