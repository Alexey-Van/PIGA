rule all_call_lr_snv:
    input:
        f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.beagle.vcf.gz"


rule lr_ont_minimap2_map:
    input:
        ont_fq = config['lr_fastqs'],
        ref = config['reference']['CHM13']
    output:
        ont_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.ont.srt.bam"
    threads: 8
    resources:
        mem_mb = 60*1024
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.ref} {input.ont_fq} | \
            samtools sort -@ {threads} -o {output.ont_bam} -
        samtools index -@ {threads} {output.ont_bam}
        """

#How to run
# back to singularity temporarily.
# I add -B /storage:/storage for the test because our reference file is a link.
#             -B /storage:/storage \
#             -B $(pwd):/project \
#             -B $(pwd)/{params.TMP_DIR}:/tmp \
#             -B {params.REF_DIR}:/project/{params.REF_DIR} \
rule lr_ont_dv:
    input:
        ont_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.ont.srt.bam",
        ref = os.path.abspath(config['reference']['CHM13'])
    output:
        dv_vcf = temp("c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.vcf.gz"),
        dv_gvcf = "c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.g.vcf.gz"
    params:
        PMD_DIR = "c2_call_lr_snv/sample_gvcf/{sample}/pepper_margin_deepvariant"
    singularity:
        "file://scripts/call_lr_snv/pepper_deepvariant.sif"
    threads: 8
    resources:
        mem_mb = 60*1024
    shell:
        """
        rm -rf {params.PMD_DIR}
        mkdir -p {params.PMD_DIR}
        run_pepper_margin_deepvariant call_variant \
            --bam {input.ont_bam} --fasta {input.ref} \
            --output_dir {params.PMD_DIR} --output_prefix {wildcards.sample} \
            --threads {threads} --sample_name {wildcards.sample} \
            --ont_r10_q20 --gvcf
        final_vcf=$(find {params.PMD_DIR} -name '{wildcards.sample}.vcf.gz' -print -quit)
        final_gvcf=$(find {params.PMD_DIR} -name '{wildcards.sample}.g.vcf.gz' -print -quit)
        test -n "$final_vcf" && test -n "$final_gvcf"
        cp "$final_vcf" {output.dv_vcf}
        cp "$final_gvcf" {output.dv_gvcf}
        tabix -f {output.dv_gvcf}
        """

rule prepare_chr_bed:
    input:
        ref_fai = config['reference']['CHM13'] + '.fai'
    output:
        chr_bed = "c2_call_lr_snv/chr_vcf/chr_bed/chm13.{chrom}.bed"
    resources:
        mem_mb = 10*1024
    shell:
        """
        awk -v OFS='\\t' -v chrom={wildcards.chrom} '{{if($1==chrom)print $1,"0",$2}}' {input.ref_fai} > {output.chr_bed}
        """

# gvcf_map = f"c2_call_lr_snv/sample_gvcf/{config['prefix']}.deepvariant.gvcf.map",
        #ls {input.dv_gvcfs} > {output.gvcf_map}
rule lr_glnexus:
    input:
        ref = config['reference']['CHM13'],
        chr_bed = "c2_call_lr_snv/chr_vcf/chr_bed/chm13.{chrom}.bed",
        dv_gvcfs = expand("c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.g.vcf.gz", sample=samples_list)
    output:
        glnexus_bcf = f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.{{chrom}}.bcf",
        glnexus_vcf = f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.{{chrom}}.vcf.gz",
        glnexus_vcf_tbi = f"c2_call_lr_snv/chr_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.{{chrom}}.vcf.gz.tbi"
    params:
        chr_db = "c2_call_lr_snv/chr_vcf/{chrom}_db"
    threads: 16
    resources:
        max_mem_gb = 200,
        mem_mb = 200*1024
    shell:
        """

        if [ -d {params.chr_db} ];
          then rm -rf {params.chr_db}
        fi

        glnexus_cli \
            --config DeepVariant \
            --dir {params.chr_db} \
            --bed {input.chr_bed} \
            --threads {threads} \
            --mem-gbytes {resources.max_mem_gb} \
            {input.dv_gvcfs} > {output.glnexus_bcf}

        bcftools view --threads {threads} {output.glnexus_bcf} | \
            bcftools norm --threads {threads} -m -any -f {input.ref} |\
            bcftools view --threads {threads} -i "AQ>40" |\
            bcftools norm --threads {threads} -m +any |\
            bcftools view --threads {threads} -m2 -M2 -v snps |\
            bcftools norm --threads {threads} -m -any |\
            bcftools sort -o {output.glnexus_vcf}

        tabix -f {output.glnexus_vcf}
        rm -rf {params.chr_db}
        """

rule glnexus_vcf_merge:
    input:
        glnexus_vcfs = expand("c2_call_lr_snv/chr_vcf/{prefix}.deepvariant.filter.biallelic_snp.{chrom}.vcf.gz", chrom=[f'chr{i}' for i in range(1, 23)] + ['chrX'], prefix = config['prefix'])
    output:
        merged_glnexus_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.vcf.gz"
    threads: 8
    resources:
        mem_mb = 60*1024
    shell:
        """
        bcftools concat --threads {threads} -Oz -o {output.merged_glnexus_vcf} {input.glnexus_vcfs}

        tabix -f {output.merged_glnexus_vcf}
        """

rule lr_whatshap_genotype:
    input:
        merged_glnexus_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.filter.biallelic_snp.vcf.gz",
        ont_bam = "c2_call_lr_snv/sample_bam/{sample}/{sample}.ont.srt.bam",
        ref = config['reference']['CHM13']
    output:
        whatshap_tmp_snp_vcf = temp("c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.tmp.biallelic_snp.vcf.gz"),
        whatshap_snp_vcf = "c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.vcf.gz",
        whatshap_reform_snp_vcf = "c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.reform.vcf.gz"
    params:
        sex = get_sex
    threads: 8
    resources:
        mem_mb = 60*1024
    shell:
        """
        whatshap genotype -o {output.whatshap_tmp_snp_vcf} -r {input.ref} --ignore-read-groups --sample {wildcards.sample} {input.merged_glnexus_vcf} {input.ont_bam}
        tabix -f {output.whatshap_tmp_snp_vcf}

        bcftools view -s {wildcards.sample} {output.whatshap_tmp_snp_vcf} -Oz -o {output.whatshap_snp_vcf}
        tabix -f {output.whatshap_snp_vcf}

        bash scripts/call_lr_snv/whatshap_vcf_reform.sh {output.whatshap_snp_vcf} {params.sex} | bgzip -c > {output.whatshap_reform_snp_vcf}
        tabix -f {output.whatshap_reform_snp_vcf}
        """

rule lr_whatshap_vcf_merge:
    input:
        whatshap_reform_snp_vcfs = expand("c2_call_lr_snv/sample_gvcf/{sample}/{sample}.deepvariant.whatshap.biallelic_snp.reform.vcf.gz", sample = samples_list)
    output:
        merge_whatshap_snp_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.vcf.gz",
        merge_whatshap_filter_snp_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.vcf.gz",
        merge_whatshap_filter_snp_vcf_tbi = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.vcf.gz.tbi"
    threads: 4
    resources:
        mem_mb = 30*1024
    shell:
        """
        bcftools merge --threads {threads} -m none {input.whatshap_reform_snp_vcfs} -Oz -o {output.merge_whatshap_snp_vcf}
        bcftools plugin fill-tags --threads {threads} {output.merge_whatshap_snp_vcf} | \
        bcftools view --threads {threads} -e "(HWE<1e-6 && ExcHet<0.1) || F_MISSING>0.1" | \
        bcftools view --threads {threads} -i "(CHR!='chrX' && MEDIAN(GQ)>50 && AC!=0) || (CHR=='chrX' && MEDIAN(GQ)>40 && AC!=0)" -o {output.merge_whatshap_filter_snp_vcf}
        tabix {output.merge_whatshap_filter_snp_vcf}
        """

rule lr_beagle_split_vcf:
    input:
        whatshap_filter_snp_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.vcf.gz"
    output:
        chr_num_vcf = f"c2_call_lr_snv/chr_vcf/{{chrom}}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{{chrom}}.1.vcf.gz"
    params:
        prefix = f"c2_call_lr_snv/chr_vcf/{{chrom}}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{{chrom}}"
    threads: 1
    resources:
        mem_mb = 30*1024
    shell:
        """
        mkdir -p c2_call_lr_snv/lr_beagle/{wildcards.chrom}

        bcftools view {input.whatshap_filter_snp_vcf} {wildcards.chrom} | \
            java -jar scripts/call_lr_snv/splitvcf.jar {wildcards.chrom} 50000 3000 {params.prefix}
        
        """

# checkpoint for counting all of the "chr"-"num" pairs.
checkpoint generate_beagle_split_num_list:
    input:
        chr_num_vcfs = expand("c2_call_lr_snv/chr_vcf/{chrom}/{prefix}.deepvariant.whatshap.biallelic_snp.filter.{chrom}.1.vcf.gz", chrom = [f"chr{i}" for i in range(1, 23)] + ["chrX"], prefix = config['prefix'])
    output:
        beagle_split_num_list = "c2_call_lr_snv/chr_vcf/beagle_split_num.list"
    params:
        file_prefix = f"c2_call_lr_snv/chr_vcf/*/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter"
    threads: 1
    resources:
        mem_mb = 30*1024
    shell:
        """
        ls {params.file_prefix}.*.vcf.gz > {output.beagle_split_num_list}
        """


def get_chr_beagle_result_files(wildcards, prefix):
    chr_num_combination = checkpoints.generate_beagle_split_num_list.get().output[0]
    pairs = []
    with open(chr_num_combination) as f:
        for line in f:
            chrom, num = line.strip().split("/")[-1].split(".")[5:7]
            if chrom == wildcards.chrom:
                pairs.append(f"c2_call_lr_snv/chr_vcf/{chrom}/{prefix}.deepvariant.whatshap.beagle.{chrom}.{num}.vcf.gz")
    return pairs


rule lr_beagle:
    input:
        chr_num_vcf = f"c2_call_lr_snv/chr_vcf/{{chrom}}/{config['prefix']}.deepvariant.whatshap.biallelic_snp.filter.{{chrom}}.{{num}}.vcf.gz",
    output:
        chr_num_beagle_vcf = f"c2_call_lr_snv/chr_vcf/{{chrom}}/{config['prefix']}.deepvariant.whatshap.beagle.{{chrom}}.{{num}}.vcf.gz"
    params:
        prefix = f"c2_call_lr_snv/chr_vcf/{{chrom}}/{config['prefix']}.deepvariant.whatshap.beagle.{{chrom}}.{{num}}"
    threads: 8
    resources:
        mem_mb = 60*1024
    shell:
        """
        beagle \
            nthreads={threads} \
            gl={input.chr_num_vcf} \
            out={params.prefix}
        """

rule chr_num_beagle_vcf_merge:
    input:
        chr_num_beagle_vcfs = partial(get_chr_beagle_result_files, prefix = config['prefix'])
    output:
        chr_beagle_vcf = f"c2_call_lr_snv/chr_vcf/{{chrom}}/{config['prefix']}.deepvariant.whatshap.beagle.{{chrom}}.vcf.gz"
    threads: 1
    resources:
        max_mem_gb = 20,
        mem_mb = 20*1024
    run:
        vcf_num = len(input.chr_num_beagle_vcfs)

        ## if only one file on each chromosome:
        if vcf_num == 1:
          shell(f"""
            mv {input.chr_num_beagle_vcfs} {output.chr_beagle_vcf}
            tabix -f {output.chr_beagle_vcf}
                """)
        ## if several splitted files.
        else:
          shell(f"""
            java "-Xmx{resources.max_mem_gb}G" \
                -jar scripts/call_lr_snv/mergevcf.jar \
                chr{wildcards.chrom} \
                {input.chr_num_beagle_vcfs} | \
                bgzip -c > {output.chr_beagle_vcf}
            tabix -f {output.chr_beagle_vcf}
                """)


rule concat_beagle_vcf:
    input:
        chr_beagle_vcfs = expand("c2_call_lr_snv/chr_vcf/{chrom}/{prefix}.deepvariant.whatshap.beagle.{chrom}.vcf.gz", chrom = [f"chr{i}" for i in range(1, 23)] + ["chrX"], prefix = config['prefix'])
    output:
        concat_beagle_vcf = f"c2_call_lr_snv/merged_vcf/{config['prefix']}.deepvariant.whatshap.beagle.vcf.gz"
    threads: 16
    resources:
        mem_mb = 120*1024
    shell:
        """        
        bcftools concat --threads {threads} -a {input.chr_beagle_vcfs} | \
        bcftools annotate -x FORMAT/DS,FORMAT/GP | \
        bcftools plugin fill-tags --threads {threads} | \
        bcftools view --threads {threads} -i "AC!=0" -o {output.concat_beagle_vcf}
        tabix -f {output.concat_beagle_vcf}
        """
# rule xxx:
#     input:
#     output:
#     threads:
#     resources:
#         mem_mb = 
#     shell:
#         """
#         """
