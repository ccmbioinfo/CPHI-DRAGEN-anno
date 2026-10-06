#!/usr/bin/env python3

"""Build a static family STR report from annotated calls and REViewer SVGs.

FlipBook expects one image directory per table row and uses a metadata TSV for
both the overview table and locus pages. This script reshapes the pipeline's
per-sample REViewer output into that layout, then adds small report-specific
HTML and CSS fragments to the generated site.

"""

import csv
import html
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import quote


def safe_name(value):
    """Return a name safe to use for FlipBook directories and SVG files."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")


def esc(value):
    """Escape metadata for HTML and keep the pipeline's '.' missing-value style."""
    return html.escape(str(value if value not in (None, "") else "."), quote=True)


def line(label, value):
    """Format one labelled value inside a metadata or sample cell."""
    return (
        f'<b class="str-label">{html.escape(label)}:</b> '
        f'<span class="str-value">{esc(value)}</span>'
    )


def sequence(value):
    # Soft breaks keep long motifs and locus structures from widening the table.
    chunks = re.findall(r".{1,12}", str(value))
    return f'<span class="str-sequence">{"<wbr>".join(esc(chunk) for chunk in chunks)}</span>'


def make_metadata(report_rows, catalog, samples):
    """Translate the annotated report into FlipBook metadata rows."""

    # FlipBook treats Path as the image-directory name. Every other field is
    # rendered as HTML in both the overview table and the corresponding page.
    metadata = []
    for report in report_rows:
        locus = report["STRCHIVE_LOCUS_ID"]
        url = report.get("STRCHIVE_URL", ".")
        link = (
            ""
            if url == "."
            else f'<br><a href="{esc(url)}" target="_blank" rel="noopener noreferrer">STRchive ↗</a>'
        )
        row = {
            "Path": safe_name(locus),
            "Gene": esc(report.get("GENE")) + link,
            "Disorder": esc(report.get("DISORDER")),
            "Motif": sequence(report.get("MOTIF")),
            "Locus structure": sequence(catalog[locus]["LocusStructure"]),
            "Ranges": '<div class="str-ranges">' + "".join(
                f"<span>{value}</span>"
                for value in (
                    line("Benign", report.get("BENIGN_RANGES")),
                    line("Intermediate", report.get("INTERMEDIATE_RANGES")),
                    line("Pathogenic", report.get("PATHOGENIC_RANGES")),
                )
            ) + "</div>",
            "Target": '<div class="str-target">' + "".join(
                f"<span>{value}</span>"
                for value in (
                    line("Region", report.get("TARGET_REGION")),
                    line("Variant", report.get("TARGET_VARIANT_ID")),
                )
            ) + "</div>",
            "Note": f'<div class="str-note">{esc(report.get("NOTE"))}</div>',
        }
        for sample in samples:
            prediction = report.get(f"{sample}_DISEASE_PREDICTION", ".")
            prediction_class = re.sub(
                r"[^a-z0-9]+", "-", str(prediction).lower()
            ).strip("-")
            # Each sample is one FlipBook column containing its primary call
            # fields followed by the three REViewer evidence categories.
            calls = "".join(
                f"<span>{value}</span>"
                for value in (
                    line("GT", report.get(f"{sample}_GT")),
                    line("Repeat", report.get(f"{sample}_motif_count")),
                    line("CI", report.get(f"{sample}_REPCI")),
                    line("Filter", report.get(f"{sample}_FILTER")),
                    line("Coverage", report.get(f"{sample}_coverage")),
                    line("Support", report.get(f"{sample}_SO")),
                )
            )
            evidence = "".join(
                f"<span>{value}</span>"
                for value in (
                    line("Spanning", report.get(f"{sample}_ADSP")),
                    line("Flanking", report.get(f"{sample}_ADFL")),
                    line("In-repeat", report.get(f"{sample}_ADIR")),
                )
            )
            row[sample] = (
                '<div class="str-sample-card"><div class="str-card-badges">'
                + f'<span class="str-badge prediction-{prediction_class}">{esc(prediction)}</span>'
                + '</div><div class="str-call-grid">'
                + calls
                + '</div><div class="str-evidence"><b>Evidence</b><div>'
                + evidence
                + "</div></div></div>"
            )
        metadata.append(row)
    return metadata


def write_headers(image_root, family, locus_count, css):
    """Write the HTML fragments FlipBook injects into its generated pages."""

    # FlipBook copies these specially named fragments into the generated main
    # and locus pages. Embedding the CSS keeps the static report self-contained.
    title = esc(family)

    # The overview keeps only fields useful for scanning a family. DataTables
    # still provides sorting and paging; this adds the report-specific filter.
    main_header = f"""
<style>{css}</style>
<!-- Inserted above FlipBook's generated DataTable on the overview page. -->
<div class="str-report-header">
  <h1>{title} STR review</h1>
  <p>{locus_count} disease-associated loci. Select a row to view the family REViewer plots.</p>
  <label for="str-table-filter"><b>Search loci and family calls</b></label><br>
  <input id="str-table-filter" type="search" placeholder="Gene, disorder, locus, prediction, sample…">
</div>
<script>
document.addEventListener('DOMContentLoaded', function () {{
  const table = window.jQuery('#data-table').DataTable();
  table.columns().every(function () {{
    const title = this.header().textContent.trim();
    // Path contains FlipBook's direct link to the generated locus page.
    if (title === 'Path') {{
      this.visible(true);
      this.header().textContent = 'Locus page';
    }}
    // These fields remain on locus pages but would make the overview too wide.
    if (['Locus structure', 'Ranges', 'Target', 'Note'].includes(title)) this.visible(false);
  }});
  document.getElementById('str-table-filter').addEventListener('input', function () {{
    table.search(this.value).draw();
  }});
}});
</script>
""".strip()

    # FlipBook emits all locus metadata and images as one vertical sequence.
    # This fragment groups them into compact metadata grids and plot sections.
    data_header = f"""
<style>{css}</style>
<!-- Inserted above the metadata and images on every FlipBook locus page. -->
<div class="str-report-header">
  <h1>{title} STR review</h1>
  <p>Family calls and evidence are shown above the sample-ordered REViewer plots.</p>
</div>
<script>
document.addEventListener('DOMContentLoaded', function () {{
  const column = document.querySelector('.ui.grid > .row > .fourteen.wide.column');
  // FlipBook writes metadata as sibling inline blocks. Split those blocks into
  // a locus grid and a sample grid without changing the generated values.
  const items = Array.from(column.querySelectorAll(':scope > div[style*="display: inline-block"]'));
  const locusGrid = document.createElement('div');
  const sampleGrid = document.createElement('div');
  locusGrid.className = 'str-locus-grid';
  sampleGrid.className = 'str-sample-grid';
  column.insertBefore(locusGrid, items[0]);
  column.insertBefore(sampleGrid, items[0]);
  items.forEach(function (item) {{
    const keyElement = item.firstElementChild;
    const key = keyElement.textContent.trim().toLowerCase();
    // FlipBook inserts a text-node colon between each generated key and value.
    if (keyElement.nextSibling && keyElement.nextSibling.nodeType === Node.TEXT_NODE) {{
      keyElement.nextSibling.textContent = keyElement.nextSibling.textContent.replace(/^\\s*:\\s*/, '');
    }}
    item.classList.add('str-meta-item', 'str-meta-' + key.replace(/[^a-z0-9]+/g, '-'));
    (item.querySelector('.str-sample-card') ? sampleGrid : locusGrid).appendChild(item);
  }});

  // Each plot is emitted as four siblings: divider, anchor, caption and image.
  // Wrap them so plots can be centered and horizontally scrolled when zoomed.
  const sections = Array.from(document.querySelectorAll('img.img-default-view')).map(function (image) {{
    const heading = image.previousElementSibling;
    heading.classList.add('str-figure-caption');
    const anchor = heading.previousElementSibling;
    const divider = anchor.previousElementSibling;
    const wrapper = document.createElement('div');
    wrapper.className = 'str-image-section';
    divider.parentNode.insertBefore(wrapper, divider);
    wrapper.append(divider, anchor, heading, image);
    return wrapper;
  }}).sort(function (left, right) {{
    // Staged filenames begin with sample order, so lexical order is pedigree order.
    return left.querySelector('img').getAttribute('src').localeCompare(
      right.querySelector('img').getAttribute('src')
    );
  }});
  sections.forEach(function (section, index) {{
    // Rebuild FlipBook's section targets after sorting the sample plots.
    section.querySelector('a').setAttribute('name', 'section' + (index + 1));
    section.querySelector('img').setAttribute('tabindex', index + 1);
    column.appendChild(section);
  }});
}});
</script>
""".strip()
    (image_root / "flipbook_main_page_header.html").write_text(main_header + "\n")
    (image_root / "flipbook_data_page_header.html").write_text(data_header + "\n")


def build_site(job):
    """Stage REViewer plots, run FlipBook, and write the report launcher."""

    with open(job.input.report, newline="") as handle:
        report_rows = list(csv.DictReader(handle))
    with open(job.input.samples_tsv, newline="") as handle:
        samples = [row["sample"] for row in csv.DictReader(handle, delimiter="\t")]

    reviewer_dirs = [Path(path) for path in job.input.reviewer_dirs]
    with open(job.input.catalog) as handle:
        catalog = {row["LocusId"]: row for row in json.load(handle)}

    output_dir = Path(job.output.metadata).parent
    image_root = output_dir / "_flipbook"
    shutil.rmtree(image_root, ignore_errors=True)
    image_root.mkdir(parents=True)

    # FlipBook makes one page per directory. Filename prefixes preserve
    # the sample order from the pedigree TSV when FlipBook sorts the SVGs.
    for report in report_rows:
        locus = report["STRCHIVE_LOCUS_ID"]
        locus_dir = image_root / safe_name(locus)
        locus_dir.mkdir()
        image_count = 0
        for order, (sample, reviewer_dir) in enumerate(
            zip(samples, reviewer_dirs), start=1
        ):
            prefix = f"{order:02d}_{safe_name(sample)}"
            source_dir = reviewer_dir / safe_name(locus)
            source = next(source_dir.glob(f"{safe_name(sample)}*.svg"), None)
            if source:
                destination = locus_dir / f"{prefix}.svg"
                shutil.copy2(source, destination)
                image_count += 1

        if image_count == 0:
            # A non-empty directory is required for FlipBook to retain the locus.
            (locus_dir / "00_no_reviewer_svg.svg").write_text(
                '<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="300">\n'
                '<rect width="100%" height="100%" fill="#f6f8fa" stroke="#d0d7de"/>\n'
                f'<text x="40" y="55" font-family="sans-serif" font-size="25" font-weight="bold">{esc(locus)}</text>\n'
                '<text x="40" y="92" font-family="sans-serif" font-size="15">No REViewer SVG was generated for this family.</text>\n'
                "</svg>\n"
            )

    # Create TSV for FlipBook's input.
    metadata = make_metadata(report_rows, catalog, samples)
    metadata_path = Path(job.output.metadata)
    with open(metadata_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(metadata[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(metadata)

    css = (Path(job.scriptdir) / "str_reviewer.css").read_text()
    write_headers(image_root, output_dir.name, len(report_rows), css)
    # FlipBook creates a portable static site under image_root/flipbook_html.
    log_path = Path(job.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log:
        subprocess.check_call(
            [
                sys.executable,
                "-m",
                "flipbook",
                "--generate-static-website",
                "-m",
                str(metadata_path),
                str(image_root),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    # Keep only the finished static site; the copied SVG staging tree is no
    # longer needed after FlipBook has populated flipbook_html.
    shutil.move(image_root / "flipbook_html", Path(job.output.site))
    shutil.rmtree(image_root)

    # The report target stays under reports/, while the larger site remains with
    # the STR working files. This tiny page redirects browsers to that site.
    url = quote("../" + Path(job.output.site).as_posix() + "/index.html", safe="/")
    launcher = Path(job.output.launcher)
    launcher.parent.mkdir(parents=True, exist_ok=True)
    launcher.write_text(f"""<!doctype html>
<html lang="en"><head>
  <meta charset="utf-8">
  <meta http-equiv="refresh" content="0; url={url}">
  <title>{esc(output_dir.name)} STR review</title>
</head><body>
  <p>Opening the <a href="{url}">{esc(output_dir.name)} STR review</a>…</p>
</body></html>
""")


if __name__ == "__main__":
    build_site(snakemake)
