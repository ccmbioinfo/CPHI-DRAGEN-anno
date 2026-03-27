# Document Conversion

Use `doc_conversion.sh` to convert documents between Markdown and Word format in this directory.

## Usage

Submit the script through Slurm with one input file:

```bash
sbatch doc_conversion.sh DRAGEN_SV_CNV_report_March_2026.md
sbatch doc_conversion.sh DRAGEN_SV_CNV_report_March_2026.docx
```

The script should:

- load the `pandoc` module
- convert `.md` to `.docx` using `pandoc-reference.docx`
- convert `.docx` to `.md` using `pandoc -f docx -t gfm --wrap=none`
- normalize Markdown formatting after `docx -> md` conversion

## Notes

- `pandoc-reference.docx` is only required for `.md -> .docx`
- `pandoc-reference.docx` is a Word template/reference file that Pandoc uses to style generated `.docx` documents. It controls formatting such as fonts, heading styles, table borders, and spacing.
- To change how generated Word documents look, open `pandoc-reference.docx` in Word, update the styles or formatting you want to use, and save it in place. Later `.md -> .docx` conversions will pick up those changes.
- output files are written beside the input file
- the generated Markdown is normalized to reduce extra spacing in lists and compact table formatting after conversion to make it look a little nicer. 
