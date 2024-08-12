from scripts.utils import allocated,ignore

rule somatic_tn_gridss_prep:
    input:
        ok1   = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss", "preprocess.ok"),
        ok0   = lambda wildcards: "{}/15.germline_sv_gridss/{}/_gridss/preprocess.ok".format(config['workdir'], dic_tumor_to_normal[wildcards.sample]),
    output:
        ok = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss", "prep.ok"),
    params:
        workspace = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss"),
        pre1 = lambda wildcards, input: input.ok1.replace('preprocess.ok', '%s.cram.gridss.working'%(wildcards.sample)),
        pre0 = lambda wildcards, input: input.ok0.replace('preprocess.ok', '%s.cram.gridss.working'%(dic_tumor_to_normal[wildcards.sample])),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn_gridss_prep", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn_gridss_prep", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn_gridss_prep", cluster))
    shell:
        "cd {params.workspace} \n"
        "rm -rf * \n"
        "ln -s {params.pre1} \n"
        "ln -s {params.pre0} \n"
        "touch {output.ok}"


rule somatic_tn__gridss_assemble_shards:
    input:
        cram1 = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
        ok = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss", "prep.ok"),
    output:
        ok = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss", "assemble_{shard}.ok"),
    params:
        workspace = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss"),
        bam       = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.assemble.bam"),
        shard     = "{shard}"
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__gridss_assemble_shards", "{sample}.{shard}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__gridss_assemble_shards", "{sample}.{shard}.e"),
    threads:
        int(allocated("threads", "somatic_tn__gridss_assemble_shards", cluster))
    container:
        config['container']['gridss']
    shell:
        "cd {params.workspace}"
        "  > {log.out} 2> {log.err}\n"
        "gridss "
        "  -s assemble "
        "  -r {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -t {threads}"
        "  -a {params.bam}"
        "  -w {params.workspace}"
        "  --jobnodes {config[parameter][gridss_shards]}"
        "  --jobindex {params.shard}"
        "  {input.cram0}"
        "  {input.cram1}"
        "  >> {log.out} 2>> {log.err}\n"
        "touch {output.ok}"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"


rule somatic_tn__gridss_call:
    input:
        cram1 = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
        ok =   lambda wildcards: [ join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss", "assemble_{shard}.ok".format(shard=i)) for i in range(config['parameter']['gridss_shards'])]
    output:
        bam = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.assemble.bam"),
        vcf = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.gridss.vcf"),
    params:
        workspace = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "_gridss"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__gridss_call", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__gridss_call", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__gridss_call", cluster))
    container:
        config['container']['gridss']
    shell:
        "gridss "
        "  -s assemble,call "
        "  -r {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -t {threads}"
        "  -a {output.bam}"
        "  -o {output.vcf}"
        "  -w {params.workspace}"
        "  {input.cram0}"
        "  {input.cram1}"
        "  > {log.out} 2> {log.err}\n"


rule somatic_tn__gripss:
    input:
        vcf = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.gridss.vcf"),
    output:
        vcf = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.gripss.filtered.vcf.gz"),
    params:
        dir = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}"),
        normal = lambda wildcards: dic_tumor_to_normal[wildcards.sample],
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__gripss", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__gripss", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__gripss", cluster))
    container:
        config['container']['gripss']
    shell:
        "java -Xmx32g -jar /usr/local/share/hmftools-gripss-2.4-0/gripss.jar "
        "  -sample {wildcards.sample} "
        "  -reference {params.normal} "
        "  -ref_genome_version 38 "
        "  -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "  -pon_sgl_file {config[references][hmftools]}/ref/38/sv/sgl_pon.38.bed.gz "
        "  -pon_sv_file {config[references][hmftools]}/ref/38/sv/sv_pon.38.bedpe.gz "
        "  -known_hotspot_file {config[references][hmftools]}/ref/38/sv/known_fusions.38.bedpe "
        "  -repeat_mask_file {config[references][hmftools]}/ref/38/sv/repeat_mask_data.38.fa.gz "
        "  -vcf {input.vcf} "
        "  -output_dir {params.dir}"
        "  > {log.out} 2> {log.err}\n"
