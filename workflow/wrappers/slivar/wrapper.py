from snakemake.shell import shell

shell.executable("bash")

mode = snakemake.params.mode
vcf = snakemake.input.vcf
js = snakemake.params.js
custom_order = snakemake.params.order
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def run(output_vcf, expr, min_alt_depth=None, extra_env=None):
    env = f"SLIVAR_IMPACTFUL_ORDER={custom_order} "
    if extra_env:
        env = f"{env}{extra_env}"
    sample_expr = ""
    if min_alt_depth is not None:
        sample_expr = (
            "--sample-expr "
            f"'ad_pass:passes_alt_depth(sample, {min_alt_depth}) || all_alt_depths_missing()' "
        )
    shell(
        f"""({env}slivar expr \
        --vcf {vcf} \
        --js {js} \
        {sample_expr}\
        --info '{expr}' \
        --pass-only \
        -o {output_vcf} \
        && bgzip -f {output_vcf} \
        && tabix -f {output_vcf}.gz) {log}"""
    )


if mode == "coding":
    run(
        snakemake.output.rare_impactful[:-3],
        """INFO.impactful &&
variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
(
  !("gnomad_fafmax_faf95_max" in INFO) ||
  !present(INFO.gnomad_fafmax_faf95_max) ||
  INFO.gnomad_fafmax_faf95_max <= 0.01
)""",
        min_alt_depth=3,
    )
    run(
        snakemake.output.rare_clinvar[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
        min_alt_depth=1,
    )
    run(
        snakemake.output.common_pathogenic_clinvar[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
    )
elif mode == "wgs-high-impact":
    run(
        snakemake.output.rare_main[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
(
  !("gnomad_fafmax_faf95_max" in INFO) ||
  !present(INFO.gnomad_fafmax_faf95_max) ||
  INFO.gnomad_fafmax_faf95_max <= 0.001
)""",
        min_alt_depth=3,
    )
    run(
        snakemake.output.rare_clinvar[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
        min_alt_depth=1,
    )
    run(
        snakemake.output.common_pathogenic_clinvar[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
    )
elif mode == "wgs":
    run(
        snakemake.output.rare_main[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
(
  !("gnomad_fafmax_faf95_max" in INFO) ||
  !present(INFO.gnomad_fafmax_faf95_max) ||
  INFO.gnomad_fafmax_faf95_max <= 0.01
)""",
        min_alt_depth=3,
    )
    run(
        snakemake.output.rare_clinvar[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
        min_alt_depth=1,
    )
    run(
        snakemake.output.common_pathogenic_clinvar[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
    )
elif mode == "compound-hets":
    run(
        snakemake.output.candidates[:-3],
        """variant.FILTER == "PASS" &&
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
)""",
    )
else:
    raise ValueError(f"Unsupported slivar mode: {mode}")
