#!/usr/bin/env bash
#SBATCH --job-name=doc-convert
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /hpf/largeprojects/ccm_dccforge/dccdipg/Common/conda_envs/pandoc

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input.md|input.docx>"
  exit 1
fi

input="$1"

if [ ! -f "$input" ]; then
  echo "Error: file not found: $input"
  exit 1
fi

dir="$(cd "$(dirname "$input")" && pwd)"
base="$(basename "$input")"
name="${base%.*}"
ext="${base##*.}"

reference_doc="$dir/pandoc-reference.docx"

normalize_md_lists() {
  local file="$1"
  local tmp
  tmp="$(mktemp)"

  awk '
  function is_bullet(line) {
    return line ~ /^[[:space:]]*[-+*][[:space:]]+/
  }

  {
    lines[NR] = $0
  }

  END {
    for (i = 1; i <= NR; i++) {
      prev_line = (i > 1 ? lines[i-1] : "")
      curr_line = lines[i]
      next_line = (i < NR ? lines[i+1] : "")

      if (curr_line == "" && is_bullet(prev_line) && is_bullet(next_line)) {
        continue
      }

      print curr_line
    }
  }' "$file" > "$tmp"

  mv "$tmp" "$file"
}

normalize_md_tables() {
  local file="$1"
  local tmp
  tmp="$(mktemp)"

  awk '
  function trim(text) {
    sub(/^[[:space:]]+/, "", text)
    sub(/[[:space:]]+$/, "", text)
    return text
  }

  function is_table_line(line) {
    return line ~ /^[[:space:]]*\|/ && line ~ /\|[[:space:]]*$/
  }

  function is_separator_cell(cell) {
    cell = trim(cell)
    return cell ~ /^:?-+:?$/
  }

  function normalize_table_line(line,    n, i, cell, out, separator_row) {
    n = split(line, parts, "|")
    if (n < 3) {
      return line
    }

    separator_row = 1
    for (i = 2; i < n; i++) {
      cell = trim(parts[i])
      if (!is_separator_cell(cell)) {
        separator_row = 0
        break
      }
    }

    out = "|"
    for (i = 2; i < n; i++) {
      cell = trim(parts[i])
      if (separator_row) {
        cell = "---"
      }
      out = out " " cell " |"
    }

    return out
  }

  {
    if (is_table_line($0)) {
      print normalize_table_line($0)
    } else {
      print $0
    }
  }' "$file" > "$tmp"

  mv "$tmp" "$file"
}

case "$ext" in
  md)
    output="$dir/$name.docx"
    if [ -f "$reference_doc" ]; then
      pandoc "$input" -o "$output" --reference-doc="$reference_doc"
    else
      pandoc "$input" -o "$output"
    fi
    echo "Wrote $output"
    ;;
  docx)
    output="$dir/$name.md"
    pandoc -f docx -t gfm --wrap=none "$input" -o "$output"
    normalize_md_lists "$output"
    normalize_md_tables "$output"
    echo "Wrote $output"
    ;;
  *)
    echo "Error: input must end in .md or .docx"
    exit 1
    ;;
esac
