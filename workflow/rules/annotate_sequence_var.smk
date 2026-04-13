def get_filt_vcf(wildcards):
    if wildcards.p == "coding":
        return "filtered/{family}.vcf.gz"
    elif wildcards.p == "wgs-high-impact":
        return "filtered/{family}.vcf.gz"
    elif wildcards.p == "denovo":
        return "filtered/{family}.denovo.vcf.gz"
    else:
        return "filtered/{p}/{family}.{p}.vcf.gz".format(p=wildcards.p,family=family)

rule input_prep:
    input:
        vcf=get_sequence_var_vcf
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
#    wildcard_constraints:
#        prefix = "(?!.*panel).*"
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
        directory("sequence_variants/{p}/{family}")
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