import os

def get_filt_vcf(wildcards):
    if wildcards.p == "coding":
        return "filtered/{family}.vcf.gz"
    elif wildcards.p == "wgs-high-impact":
        return "filtered/{family}.vcf.gz"
    elif wildcards.p == "denovo":
        return "filtered/{family}.denovo.vcf.gz"
    else:
        return "filtered/{p}/{family}.{p}.vcf.gz".format(p=wildcards.p,family=family)


slivar_script_dir = os.path.join(workflow.basedir, "scripts", "slivar")
cre_data_dir = os.path.join(workflow.basedir, "scripts", "cre", "data")


def get_slivar_report_table_inputs():
    return {
        "gene_descriptions": os.path.join(cre_data_dir, "ensembl_w_description.txt"),
        "omim": os.path.join(cre_data_dir, "OMIM_hgnc_join_omim_phenos_2025-07-10.tsv"),
        "orphanet": os.path.join(cre_data_dir, "orphanet.txt"),
        "constraint": os.path.join(cre_data_dir, "gnomad_scores_transcript_level_v4.1.1.csv"),
        "imprinting": os.path.join(cre_data_dir, "imprinting.txt"),
        "pseudoautosomal": os.path.join(cre_data_dir, "pseudoautosomal.txt"),
    }


rule input_prep:
    input:
        vcf=get_small_variant_vcf
    params:
        outdir="filtered"
    output:
        "filtered/{family}.vcf.gz"
    wildcard_constraints:
        family = "(?!.*panel|.*coding|.*denovo).*"
    log:
        "logs/input_prep/{family}.log"
    conda:
        "../envs/common.yaml"
    shell:
        '''
        ( bcftools view -e 'CHROM~"alt" || CHROM~"random" || CHROM~"Un"' -O z -o {output} {input.vcf} ) > {log} 2>&1
        '''

rule vt:
    input: get_filt_vcf # (vcf, bcf, or vcf.gz)
    output:
        temp("filtered/{p}/{family}.{p}.uniq.normalized.decomposed.vcf"),
    params:
        ref=config["ref"]["genome"],
    log:
        "logs/vt/{family}.vt.{p}.uniq.normalized.decomposed.log"
    wrapper:
        get_wrapper_path("vt")

rule pass:
    input:
        "{prefix}.{ext}"
    output:
        temp("{prefix}.pass.{ext,(vcf|vcf\.gz)}")
    threads: 6
    resources:
        mem=lambda wildcards, threads: threads * 2
    params: 
        filter = "-f PASS"
    wrapper:
        get_wrapper_path("bcftools", "view")

rule vep:
    input:
        "filtered/{p}/{family}.{p}.uniq.normalized.decomposed.pass.vcf",
    output:
        temp("annotated/{p}/vep/{family}.{p}.vep.vcf"),
    log:
        "logs/vep/{family}.vep.{p}.log"
    threads: 10
    resources:
        mem_mb = 30000
    params:
        dir=config["annotation"]["vep"]["dir"],
        dir_cache=config["annotation"]["vep"]["dir_cache"],
        ref=config["ref"]["genome"],
        phyloP100way=config["annotation"]["vep"]["phyloP100way"],
    wrapper:
        get_wrapper_path("vep")

rule vcfanno:
    input:
        "annotated/{p}/vep/{family}.{p}.vep.vcf",
    output:
        "annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.no.PS.vcf.gz",
    log:
        "logs/vcfanno/{family}.vcfanno.{p}.log"
    threads: 10
    resources:
        mem_mb = 20000
    params:
        lua_script=config["annotation"]["vcfanno"]["lua_script"],
        conf=config["annotation"]["vcfanno"]["conf"],
        base_path=config["annotation"]["vcfanno"]["base_path"],
    wrapper:
        get_wrapper_path("vcfanno")

rule add_ps_field:
    input:
       vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.no.PS.vcf.gz",
    output:
        temp("annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf"),
    log:
        "logs/bcftools/{family}.add.PS.field.{p}.log"
    conda: 
        "../envs/common.yaml"
    shell:
        '''
            #Get the PS values for every sample from the FORMAT VCF field, remove the trailing comma and store the results in PS_annot.txt
            bcftools query -f '%CHROM\t%POS\t[%PS,]\n' {input.vcf} | sed 's/,*$//g' > annotated/{wildcards.p}/vcfanno/PS_annot.txt

            #BGZIP and TABIX the PS_annot.txt file
            bgzip annotated/{wildcards.p}/vcfanno/PS_annot.txt
            tabix -s1 -b2 -e2 annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz

            #create header file containing PS INFO field info
            echo -e "##INFO=<ID=PS,Number=.,Type=String,Description="Phase set">" > annotated/{wildcards.p}/vcfanno/hdr.txt

            #Annotate the VCF with the PS_annot.txt.gz file to add PS tag info to the VCF
            bcftools annotate -a annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz -h annotated/{wildcards.p}/vcfanno/hdr.txt -c CHROM,POS,INFO/PS {input.vcf} > {output}

            #Remove intermediate files  
            rm annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz annotated/{wildcards.p}/vcfanno/PS_annot.txt.gz.tbi annotated/{wildcards.p}/vcfanno/hdr.txt
        '''

rule vcf2db:
    input:
        "annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz",
    output:
         db="annotated/{p}/{family}-gemini.db",
    log:
        "logs/vcf2db/{family}.vcf2db.{p}.log"
    threads: 1
    resources:
        mem_mb = 20000
    wrapper:
        get_wrapper_path("vcf2db")

rule bgzip:
   input:
       "{prefix}.vcf"
   output:
       "{prefix}.vcf.gz"
   conda:
       "../envs/common.yaml"
   shell:
        '''
        bgzip -c {input} > {output}
        '''
        
rule tabix:
    input: 
        "{prefix}.vcf.gz"
    output: 
        "{prefix}.vcf.gz.tbi"
    log: 
        "logs/{prefix}.log"
    conda:
           "../envs/common.yaml"
    wrapper:
        get_wrapper_path("tabix")

rule allsnvreport:
    input:
        db="annotated/{p}/{family}-gemini.db",
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz"
    output:
        directory("small_variants/{p}/{family}")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/{p}/{family}.cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         database_path=config["annotation"]["cre"]["database_path"],
         ref=config["ref"]["genome"]
    shell:
         '''
         (mkdir -p {output}
         cd {output}
         ln -s ../../../{input.db} {family}-ensemble.db
         #bgzip ../../../{input.vcf} -c > {family}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s ../../../{input.vcf} {family}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {family}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {family}-gatk-haplotype-annotated-decomposed.vcf.gz {family}-ensemble-annotated-decomposed.vcf.gz
         ln -s {family}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {family}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         if [ {wildcards.p} == "coding" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} {params.cre}/cre.sh {family} 
         elif [ {wildcards.p} == "wgs-high-impact" ]; then  
         cre={params.cre} reference={params.ref} database={params.database_path} type=wgs.high.impact {params.cre}/cre.sh {family}
         else
         cre={params.cre} reference={params.ref} database={params.database_path} type=wgs {params.cre}/cre.sh {family}
         unset type
         fi;
         ) > {log} 2>&1
         '''


rule slivar_select_coding:
    input:
        vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz"
    output:
        rare_impactful="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_impactful.vcf.gz",
        rare_impactful_tbi="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_impactful.vcf.gz.tbi",
        rare_clinvar="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_clinvar.vcf.gz",
        rare_clinvar_tbi="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_clinvar.vcf.gz.tbi",
        common_pathogenic_clinvar="small_variants_slivar/coding/{family}/branches/{family}.coding.common_pathogenic_clinvar.vcf.gz",
        common_pathogenic_clinvar_tbi="small_variants_slivar/coding/{family}/branches/{family}.coding.common_pathogenic_clinvar.vcf.gz.tbi",
        rare_clinvar_allow_missing_faf="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_clinvar_allow_missing_faf.vcf.gz",
        rare_clinvar_allow_missing_faf_tbi="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_clinvar_allow_missing_faf.vcf.gz.tbi",
    log:
        "logs/slivar/{family}.coding.select.log"
    conda:
        "../wrappers/slivar/environment.yaml"
    params:
        js=os.path.join(slivar_script_dir, "slivar_functions.js"),
        order=os.path.join(slivar_script_dir, "default-order.txt"),
        mode="coding"
    wrapper:
        get_wrapper_path("slivar")


rule slivar_postfilter_coding:
    input:
        rare_impactful="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_impactful.vcf.gz",
        rare_clinvar="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_clinvar.vcf.gz",
        common_pathogenic_clinvar="small_variants_slivar/coding/{family}/branches/{family}.coding.common_pathogenic_clinvar.vcf.gz",
        rare_clinvar_allow_missing_faf="small_variants_slivar/coding/{family}/branches/{family}.coding.rare_clinvar_allow_missing_faf.vcf.gz",
    output:
        vcf="small_variants_slivar/coding/{family}/{family}.coding.post_gemini_filter.vcf",
        key_file="small_variants_slivar/coding/{family}/{family}.coding.post_gemini_filter.keys.txt",
        audit="small_variants_slivar/coding/{family}/{family}.coding.post_gemini_filter.audit.tsv",
        reasons="small_variants_slivar/coding/{family}/{family}.coding.post_gemini_filter.reasons.tsv",
        summary="small_variants_slivar/coding/{family}/{family}.coding.post_gemini_filter.summary.txt",
    log:
        "logs/slivar/{family}.coding.postfilter.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "postfilter.py"),
        out_prefix="small_variants_slivar/coding/{family}/{family}.coding",
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (python3 {params.script} \
        --mode coding \
        --rare-impactful-vcf {input.rare_impactful} \
        --rare-clinvar-vcf {input.rare_clinvar} \
        --common-pathogenic-clinvar-vcf {input.common_pathogenic_clinvar} \
        --rare-clinvar-allow-missing-faf-vcf {input.rare_clinvar_allow_missing_faf} \
        --impact-order-file {params.order} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


rule slivar_report_coding:
    input:
        vcf="small_variants_slivar/coding/{family}/{family}.coding.post_gemini_filter.vcf",
        **get_slivar_report_table_inputs()
    output:
        "reports_slivar/{family}.wgs.coding.slivar.hg38.csv"
    log:
        "logs/slivar/{family}.coding.report.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "build_report.py"),
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (mkdir -p $(dirname {output})
        python3 {params.script} \
        --mode coding \
        --vcf {input.vcf} \
        --out-csv {output} \
        --impact-order-file {params.order} \
        --gene-descriptions {input.gene_descriptions} \
        --omim {input.omim} \
        --orphanet {input.orphanet} \
        --constraint {input.constraint} \
        --imprinting {input.imprinting} \
        --pseudoautosomal {input.pseudoautosomal}) > {log} 2>&1
        """


rule slivar_select_wgs_high_impact:
    input:
        vcf="annotated/wgs-high-impact/vcfanno/{family}.wgs-high-impact.vep.vcfanno.vcf.gz"
    output:
        rare_main="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_main.vcf.gz",
        rare_main_tbi="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_main.vcf.gz.tbi",
        rare_clinvar="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_clinvar.vcf.gz",
        rare_clinvar_tbi="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_clinvar.vcf.gz.tbi",
        rare_clinvar_allow_missing_faf="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_clinvar_allow_missing_faf.vcf.gz",
        rare_clinvar_allow_missing_faf_tbi="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_clinvar_allow_missing_faf.vcf.gz.tbi",
        common_pathogenic_clinvar="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.common_pathogenic_clinvar.vcf.gz",
        common_pathogenic_clinvar_tbi="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.common_pathogenic_clinvar.vcf.gz.tbi",
    log:
        "logs/slivar/{family}.wgs-high-impact.select.log"
    conda:
        "../wrappers/slivar/environment.yaml"
    params:
        js=os.path.join(slivar_script_dir, "slivar_functions.js"),
        order=os.path.join(slivar_script_dir, "default-order.txt"),
        mode="wgs-high-impact"
    wrapper:
        get_wrapper_path("slivar")


rule slivar_postfilter_wgs_high_impact:
    input:
        rare_main="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_main.vcf.gz",
        rare_clinvar="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_clinvar.vcf.gz",
        common_pathogenic_clinvar="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.common_pathogenic_clinvar.vcf.gz",
        rare_clinvar_allow_missing_faf="small_variants_slivar/wgs-high-impact/{family}/branches/{family}.wgs-high-impact.rare_clinvar_allow_missing_faf.vcf.gz",
    output:
        vcf="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.post_high_impact_filter.vcf",
        key_file="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.post_high_impact_filter.keys.txt",
        audit="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.post_high_impact_filter.audit.tsv",
        reasons="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.post_high_impact_filter.reasons.tsv",
        summary="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.post_high_impact_filter.summary.txt",
    log:
        "logs/slivar/{family}.wgs-high-impact.postfilter.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "postfilter.py"),
        out_prefix="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact",
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (python3 {params.script} \
        --mode wgs-high-impact \
        --rare-main-vcf {input.rare_main} \
        --rare-clinvar-vcf {input.rare_clinvar} \
        --common-pathogenic-clinvar-vcf {input.common_pathogenic_clinvar} \
        --rare-clinvar-allow-missing-faf-vcf {input.rare_clinvar_allow_missing_faf} \
        --impact-order-file {params.order} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


rule slivar_report_wgs_high_impact:
    input:
        vcf="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.post_high_impact_filter.vcf",
        **get_slivar_report_table_inputs()
    output:
        "reports_slivar/{family}.wgs.high.impact.slivar.hg38.csv"
    log:
        "logs/slivar/{family}.wgs-high-impact.report.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "build_report.py"),
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (mkdir -p $(dirname {output})
        python3 {params.script} \
        --mode wgs-high-impact \
        --vcf {input.vcf} \
        --out-csv {output} \
        --impact-order-file {params.order} \
        --gene-descriptions {input.gene_descriptions} \
        --omim {input.omim} \
        --orphanet {input.orphanet} \
        --constraint {input.constraint} \
        --imprinting {input.imprinting} \
        --pseudoautosomal {input.pseudoautosomal}) > {log} 2>&1
        """


rule slivar_select_wgs:
    input:
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf.gz"
    output:
        rare_main="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz",
        rare_main_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz.tbi",
        rare_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz",
        rare_clinvar_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz.tbi",
        rare_clinvar_allow_missing_faf="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar_allow_missing_faf.vcf.gz",
        rare_clinvar_allow_missing_faf_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar_allow_missing_faf.vcf.gz.tbi",
        common_pathogenic_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz",
        common_pathogenic_clinvar_tbi="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz.tbi",
    wildcard_constraints:
        p="panel|panel-flank|denovo"
    log:
        "logs/slivar/{family}.{p}.select.log"
    conda:
        "../wrappers/slivar/environment.yaml"
    params:
        js=os.path.join(slivar_script_dir, "slivar_functions.js"),
        order=os.path.join(slivar_script_dir, "default-order.txt"),
        mode="wgs"
    wrapper:
        get_wrapper_path("slivar")


rule slivar_postfilter_wgs:
    input:
        rare_main="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_main.vcf.gz",
        rare_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar.vcf.gz",
        common_pathogenic_clinvar="small_variants_slivar/{p}/{family}/branches/{family}.{p}.common_pathogenic_clinvar.vcf.gz",
        rare_clinvar_allow_missing_faf="small_variants_slivar/{p}/{family}/branches/{family}.{p}.rare_clinvar_allow_missing_faf.vcf.gz",
    output:
        vcf="small_variants_slivar/{p}/{family}/{family}.{p}.post_wgs_filter.vcf",
        key_file="small_variants_slivar/{p}/{family}/{family}.{p}.post_wgs_filter.keys.txt",
        audit="small_variants_slivar/{p}/{family}/{family}.{p}.post_wgs_filter.audit.tsv",
        reasons="small_variants_slivar/{p}/{family}/{family}.{p}.post_wgs_filter.reasons.tsv",
        summary="small_variants_slivar/{p}/{family}/{family}.{p}.post_wgs_filter.summary.txt",
    wildcard_constraints:
        p="panel|panel-flank|denovo"
    log:
        "logs/slivar/{family}.{p}.postfilter.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "postfilter.py"),
        out_prefix="small_variants_slivar/{p}/{family}/{family}.{p}",
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (python3 {params.script} \
        --mode wgs \
        --rare-main-vcf {input.rare_main} \
        --rare-clinvar-vcf {input.rare_clinvar} \
        --common-pathogenic-clinvar-vcf {input.common_pathogenic_clinvar} \
        --rare-clinvar-allow-missing-faf-vcf {input.rare_clinvar_allow_missing_faf} \
        --impact-order-file {params.order} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


rule slivar_report_panel_wgs:
    input:
        vcf="small_variants_slivar/{p}/{family}/{family}.{p}.post_wgs_filter.vcf",
        **get_slivar_report_table_inputs()
    output:
        "reports_slivar/{family}.{p}.slivar.hg38.csv"
    wildcard_constraints:
        p="panel|panel-flank"
    log:
        "logs/slivar/{family}.{p}.report.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "build_report.py"),
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (mkdir -p $(dirname {output})
        python3 {params.script} \
        --mode wgs \
        --vcf {input.vcf} \
        --out-csv {output} \
        --impact-order-file {params.order} \
        --gene-descriptions {input.gene_descriptions} \
        --omim {input.omim} \
        --orphanet {input.orphanet} \
        --constraint {input.constraint} \
        --imprinting {input.imprinting} \
        --pseudoautosomal {input.pseudoautosomal}) > {log} 2>&1
        """


rule slivar_report_denovo_wgs:
    input:
        vcf="small_variants_slivar/denovo/{family}/{family}.denovo.post_wgs_filter.vcf",
        **get_slivar_report_table_inputs()
    output:
        "reports_slivar/{family}.wgs.denovo.slivar.hg38.csv"
    log:
        "logs/slivar/{family}.denovo.report.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "build_report.py"),
        order=os.path.join(slivar_script_dir, "default-order.txt")
    shell:
        """
        (mkdir -p $(dirname {output})
        python3 {params.script} \
        --mode wgs \
        --vcf {input.vcf} \
        --out-csv {output} \
        --impact-order-file {params.order} \
        --gene-descriptions {input.gene_descriptions} \
        --omim {input.omim} \
        --orphanet {input.orphanet} \
        --constraint {input.constraint} \
        --imprinting {input.imprinting} \
        --pseudoautosomal {input.pseudoautosomal}) > {log} 2>&1
        """


rule compare_coding_report_keys:
    input:
        gemini=f"reports/{{family}}.wgs.coding.CH{sf_suffix}.hg38.csv",
        slivar="reports_slivar/{family}.wgs.coding.slivar.hg38.csv",
    output:
        summary="reports_slivar_compare/{family}.wgs.coding.summary.tsv",
        shared="reports_slivar_compare/{family}.wgs.coding.shared.csv",
        gemini_only="reports_slivar_compare/{family}.wgs.coding.gemini_only.csv",
        slivar_only="reports_slivar_compare/{family}.wgs.coding.slivar_only.csv",
    log:
        "logs/slivar/{family}.coding.compare.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "compare_report_variant_keys.py"),
        out_prefix="reports_slivar_compare/{family}.wgs.coding"
    shell:
        """
        (python3 {params.script} \
        --gemini-report {input.gemini} \
        --slivar-report {input.slivar} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


rule compare_wgs_high_impact_report_keys:
    input:
        gemini=f"reports/{{family}}.wgs.high.impact.CH{sf_suffix}.hg38.csv",
        slivar="reports_slivar/{family}.wgs.high.impact.slivar.hg38.csv",
    output:
        summary="reports_slivar_compare/{family}.wgs.high.impact.summary.tsv",
        shared="reports_slivar_compare/{family}.wgs.high.impact.shared.csv",
        gemini_only="reports_slivar_compare/{family}.wgs.high.impact.gemini_only.csv",
        slivar_only="reports_slivar_compare/{family}.wgs.high.impact.slivar_only.csv",
    log:
        "logs/slivar/{family}.wgs-high-impact.compare.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "compare_report_variant_keys.py"),
        out_prefix="reports_slivar_compare/{family}.wgs.high.impact"
    shell:
        """
        (python3 {params.script} \
        --gemini-report {input.gemini} \
        --slivar-report {input.slivar} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


rule compare_panel_wgs_report_keys:
    input:
        gemini="reports/{family}.{p}.CH.hg38.csv",
        slivar="reports_slivar/{family}.{p}.slivar.hg38.csv",
    output:
        summary="reports_slivar_compare/{family}.{p}.summary.tsv",
        shared="reports_slivar_compare/{family}.{p}.shared.csv",
        gemini_only="reports_slivar_compare/{family}.{p}.gemini_only.csv",
        slivar_only="reports_slivar_compare/{family}.{p}.slivar_only.csv",
    wildcard_constraints:
        p="panel|panel-flank"
    log:
        "logs/slivar/{family}.{p}.compare.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "compare_report_variant_keys.py"),
        out_prefix="reports_slivar_compare/{family}.{p}"
    shell:
        """
        (python3 {params.script} \
        --gemini-report {input.gemini} \
        --slivar-report {input.slivar} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """


rule compare_denovo_wgs_report_keys:
    input:
        gemini=f"reports/{{family}}.wgs.denovo.CH{sf_suffix}.hg38.csv",
        slivar="reports_slivar/{family}.wgs.denovo.slivar.hg38.csv",
    output:
        summary="reports_slivar_compare/{family}.wgs.denovo.summary.tsv",
        shared="reports_slivar_compare/{family}.wgs.denovo.shared.csv",
        gemini_only="reports_slivar_compare/{family}.wgs.denovo.gemini_only.csv",
        slivar_only="reports_slivar_compare/{family}.wgs.denovo.slivar_only.csv",
    log:
        "logs/slivar/{family}.denovo.compare.log"
    conda:
        "../envs/slivar.yaml"
    params:
        script=os.path.join(slivar_script_dir, "compare_report_variant_keys.py"),
        out_prefix="reports_slivar_compare/{family}.wgs.denovo"
    shell:
        """
        (python3 {params.script} \
        --gemini-report {input.gemini} \
        --slivar-report {input.slivar} \
        --out-prefix {params.out_prefix}) > {log} 2>&1
        """

if config["run"]["hpo"]:

    def get_panel(wildcards):
        if wildcards.p == "panel":
            return "genes/{family}.bed"
        else:
            return "genes/{family}_{p}.bed"

    rule hpo_to_panel:
        input: 
            hpo=config["run"]["hpo"],
            ensembl=config["genes"]["ensembl"],
            refseq=config["genes"]["refseq"],
            hgnc=config["genes"]["hgnc"]
        output: 
            genes="genes/{family}.bed"
        wildcard_constraints:
            family = "(?!.*panel|.*coding|.*denovo).*"
        conda: "../envs/hpo_to_panel.yaml"
        log: "logs/hpo_to_panel/{family}.log"
        script:
            "../scripts/hpo_to_panel.py"

    rule add_flank:
        input: "genes/{family}.bed"
        output: "genes/{family}_{p}.bed"
        params: config["run"]["flank"]
        conda: "../envs/hpo_to_panel.yaml"
        shell:
            '''
            cat {input} | awk -F "\t" '{{print $1"\t"$2-{params}"\t"$3+{params}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output}
            '''

    rule intersect:
        input: 
            left="filtered/{family}.vcf.gz",
            right=get_panel
        output:
            vcf="filtered/{p}/{family}.{p}.vcf.gz"
        params:
            extra="-header"
        log: "logs/report/bedtools-{family}-{p}.log"
        wrapper:
            get_wrapper_path("bedtools", "intersect")

rule filter_denovo:
    input:
        "filtered/{family}.vcf.gz",
    output:
        "filtered/{family}.denovo.vcf.gz",
    log:
        "logs/filter_denovo/{family}.denovo.log",
    conda:
        "../envs/common.yaml"
    shell:
        '''
        (bcftools filter -i "FORMAT/DN == 'DeNovo'" -O z -o {output} {input}) > {log} 2>&1
        '''
