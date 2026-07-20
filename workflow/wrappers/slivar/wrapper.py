from textwrap import dedent

from snakemake.shell import shell

shell.executable("bash")

mode = snakemake.params.mode
vcf = snakemake.input.vcf
js = snakemake.params.js
consequence_order_file = snakemake.params.consequence_order_file
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def run_slivar_expr(output_vcf_gz, info_expr, min_alt_depth=None):
    output_vcf = str(output_vcf_gz).removesuffix(".gz")
    info_expr = dedent(info_expr).strip()

    # Slivar sets INFO.impactful for consequences ranked before IMPACTFUL_CUTOFF in this file.
    slivar_env = f"SLIVAR_IMPACTFUL_ORDER={consequence_order_file} "

    # Slivar first applies --info to the variant-level fields in each record. For records
    # that pass, --sample-expr evaluates the expression once per sample; the name before
    # the colon becomes an INFO field containing the IDs of samples where it returned true.
    # With --pass-only, a sample expression must contain at least one passing sample ID.
    sample_expr = ""
    if min_alt_depth is not None:
        # Here, INFO/ad_pass lists samples meeting the ALT-depth cutoff. If every sample
        # lacks usable AD, the second condition passes all samples so the record is retained.
        sample_expr = (
            "--sample-expr "
            f"'ad_pass:passes_alt_depth(sample, {min_alt_depth}) || all_alt_depths_missing()' "
        )

    shell(
        f"""({slivar_env}slivar expr \
        --vcf {vcf} \
        --js {js} \
        {sample_expr}\
        --info '{info_expr}' \
        --pass-only \
        -o {output_vcf} \
        && bgzip -f {output_vcf} \
        && tabix -f {output_vcf}.gz) {log}"""
    )


# Report modes share three Slivar branches: rare main, rare ClinVar, and common ClinVar rescue.
# The branches can overlap, so postfilter.py merges and deduplicates their output.
if mode in {"coding", "wgs-high-impact", "denovo", "panel", "panel-flank"}:
    rare_main_max_af = "0.01"
    if mode == "wgs-high-impact":
        rare_main_max_af = "0.001"

    rare_main_expr = dedent(
        f"""
        variant.FILTER == "PASS" &&
        variant.ALT[0] != "*" &&
        (
          !("gnomad_fafmax_faf95_max" in INFO) ||
          !present(INFO.gnomad_fafmax_faf95_max) ||
          INFO.gnomad_fafmax_faf95_max <= {rare_main_max_af}
        )
        """
    ).strip()

    if mode == "coding":
        # Coding requires a consequence ranked before IMPACTFUL_CUTOFF in the order file.
        rare_main_expr = "INFO.impactful &&\n" + rare_main_expr

    # Keep PASS, non-star variants with missing/low FAF and at least one sample with AD >= 3
    # (or family-wide missing AD); coding also requires INFO.impactful and high-impact FAF <= 0.001.
    run_slivar_expr(
        snakemake.output.rare_main,
        rare_main_expr,
        min_alt_depth=3,
    )

    # Keep PASS, non-star variants with missing/low FAF, a ClinVar value, and at least one
    # sample with AD >= 1 (or family-wide missing AD).
    run_slivar_expr(
        snakemake.output.rare_clinvar,
        """
        variant.FILTER == "PASS" &&
        variant.ALT[0] != "*" &&
        (
          !("gnomad_fafmax_faf95_max" in INFO) ||
          !present(INFO.gnomad_fafmax_faf95_max) ||
          INFO.gnomad_fafmax_faf95_max <= 0.01
        ) &&
        (
          (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
          (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
          (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
        )
        """,
        min_alt_depth=1,
    )

    # Rescue PASS, non-star variants with FAF > 0.01, an asserted pathogenic ClinVar
    # interpretation, and no sample AD requirement.
    run_slivar_expr(
        snakemake.output.common_pathogenic_clinvar,
        """
        variant.FILTER == "PASS" &&
        variant.ALT[0] != "*" &&
        ("gnomad_fafmax_faf95_max" in INFO) &&
        present(INFO.gnomad_fafmax_faf95_max) &&
        INFO.gnomad_fafmax_faf95_max > 0.01 &&
        ("clinvar_status" in INFO) &&
        present(INFO.clinvar_status) &&
        INFO.clinvar_status != "no_assertion_criteria_provided" &&
        (
          (("clinvar_pathogenic" in INFO) && contains_pathogenic(INFO.clinvar_pathogenic)) ||
          (("clinvar_sig" in INFO) && contains_pathogenic(INFO.clinvar_sig)) ||
          (("clinvar_sig_conf" in INFO) && contains_pathogenic(INFO.clinvar_sig_conf))
        )
        """,
    )

# Keep compound-het candidates with missing/low group-max AF or an exact common ClinVar
# rescue label; build_ch_tsv.py later applies AD thresholds and assigns HIGH-MED/LOW.
elif mode == "compound-hets":
    run_slivar_expr(
        snakemake.output.candidates,
        """
        variant.FILTER == "PASS" &&
        variant.ALT[0] != "*" &&
        (
          (
            !("gnomad_af_grpmax" in INFO) ||
            !present(INFO.gnomad_af_grpmax) ||
            INFO.gnomad_af_grpmax <= 0.01
          ) ||
          (
            ("gnomad_af_grpmax" in INFO) &&
            present(INFO.gnomad_af_grpmax) &&
            INFO.gnomad_af_grpmax > 0.01 &&
            ("clinvar_status" in INFO) &&
            present(INFO.clinvar_status) &&
            INFO.clinvar_status != "no_assertion_criteria_provided" &&
            (
              (("clinvar_pathogenic" in INFO) && is_ch_common_clinvar_rescue(INFO.clinvar_pathogenic)) ||
              (("clinvar_sig" in INFO) && is_ch_common_clinvar_rescue(INFO.clinvar_sig)) ||
              (("clinvar_sig_conf" in INFO) && is_ch_common_clinvar_rescue(INFO.clinvar_sig_conf))
            )
          )
        )
        """,
    )

else:
    raise ValueError(f"Unsupported slivar mode: {mode}")
