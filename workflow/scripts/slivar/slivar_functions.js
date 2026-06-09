function present(x) {
  if (x === undefined || x === null) return false;

  if (Array.isArray(x)) {
    for (var i = 0; i < x.length; i++) {
      var s = String(x[i]);
      if (s !== "" && s !== "." && s !== "None") return true;
    }
    return false;
  }

  var s = String(x);
  return s !== "" && s !== "." && s !== "None";
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
