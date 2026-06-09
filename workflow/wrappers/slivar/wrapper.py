from snakemake.shell import shell

shell.executable("bash")

mode = snakemake.params.mode
vcf = snakemake.input.vcf
js = snakemake.params.js
custom_order = snakemake.params.order
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def run(output_vcf, expr, extra_env=""):
    shell(
        f"""({extra_env}slivar expr \
        --vcf {vcf} \
        --js {js} \
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
  (("gnomad_fafmax_faf95_max" in INFO) && INFO.gnomad_fafmax_faf95_max < 0.01) ||
  !("gnomad_fafmax_faf95_max" in INFO)
)""",
    )
    run(
        snakemake.output.rare_clinvar[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
("gnomad_fafmax_faf95_max" in INFO) &&
INFO.gnomad_fafmax_faf95_max < 0.01 &&
(
  (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
)""",
    )
    run(
        snakemake.output.common_pathogenic_clinvar[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
("gnomad_fafmax_faf95_max" in INFO) &&
INFO.gnomad_fafmax_faf95_max >= 0.01 &&
("clinvar_status" in INFO) &&
present(INFO.clinvar_status) &&
INFO.clinvar_status != "no_assertion_criteria_provided" &&
(
  (("clinvar_pathogenic" in INFO) && contains_pathogenic(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && contains_pathogenic(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && contains_pathogenic(INFO.clinvar_sig_conf))
)""",
    )
    run(
        snakemake.output.rare_clinvar_allow_missing_faf[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
INFO.genic &&
(
  (
    ("gnomad_fafmax_faf95_max" in INFO) &&
    INFO.gnomad_fafmax_faf95_max < 0.01
  ) ||
  !("gnomad_fafmax_faf95_max" in INFO)
) &&
(
  (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
)""",
        extra_env=f"SLIVAR_IMPACTFUL_ORDER={custom_order} ",
    )
elif mode == "wgs-high-impact":
    run(
        snakemake.output.rare_main[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
(
  (("gnomad_fafmax_faf95_max" in INFO) && INFO.gnomad_fafmax_faf95_max <= 0.001) ||
  !("gnomad_fafmax_faf95_max" in INFO)
)""",
    )
    run(
        snakemake.output.rare_clinvar[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
("gnomad_fafmax_faf95_max" in INFO) &&
INFO.gnomad_fafmax_faf95_max <= 0.01 &&
(
  (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
)""",
    )
    run(
        snakemake.output.rare_clinvar_allow_missing_faf[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
(
  (
    ("gnomad_fafmax_faf95_max" in INFO) &&
    INFO.gnomad_fafmax_faf95_max <= 0.01
  ) ||
  !("gnomad_fafmax_faf95_max" in INFO)
) &&
(
  (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
)""",
    )
    run(
        snakemake.output.common_pathogenic_clinvar[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
("gnomad_fafmax_faf95_max" in INFO) &&
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
  (("gnomad_fafmax_faf95_max" in INFO) && INFO.gnomad_fafmax_faf95_max <= 0.01) ||
  !("gnomad_fafmax_faf95_max" in INFO)
)""",
    )
    run(
        snakemake.output.rare_clinvar[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
("gnomad_fafmax_faf95_max" in INFO) &&
INFO.gnomad_fafmax_faf95_max <= 0.01 &&
(
  (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
)""",
    )
    run(
        snakemake.output.rare_clinvar_allow_missing_faf[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
(
  (
    ("gnomad_fafmax_faf95_max" in INFO) &&
    INFO.gnomad_fafmax_faf95_max <= 0.01
  ) ||
  !("gnomad_fafmax_faf95_max" in INFO)
) &&
(
  (("clinvar_pathogenic" in INFO) && present(INFO.clinvar_pathogenic)) ||
  (("clinvar_sig" in INFO) && present(INFO.clinvar_sig)) ||
  (("clinvar_sig_conf" in INFO) && present(INFO.clinvar_sig_conf))
)""",
    )
    run(
        snakemake.output.common_pathogenic_clinvar[:-3],
        """variant.FILTER == "PASS" &&
variant.ALT[0] != "*" &&
("gnomad_fafmax_faf95_max" in INFO) &&
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
else:
    raise ValueError(f"Unsupported slivar mode: {mode}")
