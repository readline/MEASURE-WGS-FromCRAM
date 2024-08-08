from scripts.utils import allocated,ignore

rule somatic_ss__mutect2_split:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = temp(join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.mutect2.vcf")),
    params:
        dir = join(config['workdir'], "41.somatic_ss_sv__gridss", "{sample}"),
        vcf = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.mutect2.vcf"),
        f1r2 = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.f1r2.tar.gz"),
        rom = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.read-orientation-model.tar.gz"),
        pst = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.getpileupsummaries.table"),
        st = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.segments.table"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_split", "{sample}.{chr}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_split", "{sample}.{chr}.e"),
    threads:
        int(allocated("threads", "somatic_ss__mutect2_split", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" Mutect2 "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.cram} "
        "    -L {wildcards.chr} "
        "    --panel-of-normals {config[references][pon]}/Mutect2.pon.vcf.gz "
        "    -O {params.vcf} "
        "    -germline-resource {config[references][gatksomatic]}/af-only-gnomad.hg38.vcf.gz "
        "    --f1r2-tar-gz {params.f1r2} "
        "    > {log.out} 2> {log.err}\n"
        "echo 1 >> {log.out}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" LearnReadOrientationModel "
        "    -O {params.rom} "
        "    -I {params.f1r2}"
        "    >> {log.out} 2>> {log.err}\n"
        "echo 2 >> {log.out}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" GetPileupSummaries "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.cram} "
        "    -V {config[references][gatksomatic]}/small_exac_common_3.hg38.vcf.gz "
        "    -L {wildcards.chr} "
        "    -O {params.pst} "
        "    >> {log.out} 2>> {log.err}\n"
        "echo 3 >> {log.out}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" CalculateContamination "
        "    -I {params.pst} "
        "    -tumor-segmentation {params.st}  "
        "    -O $sample.$chr.calculatecontamination.table "
        "    >> {log.out} 2>> {log.err}\n"
        "echo 4 >> {log.out}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" FilterMutectCalls "
        "    --reference {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -V {params.vcf} "
        "    --tumor-segmentation {params.st} "
        "    --ob-priors {params.rom} "
        "    -O {output.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        "echo 5 >> {log.out}\n"

rule somatic_ss__mutect2_merge:
    input:
        vcfs =lambda wildcards: \
             ["{}/41.somatic_ss_snvindel_mutect2/{}/chroms/{}.{}.mutect2.vcf".format(config['workdir'], wildcards.sample, wildcards.sample, chrid) for chrid in ['chr%d'%(i) for i in range(1,23)]+['chrX','chrY']],
    output:
        vcf  = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.vcf.gz"),
        vcfp = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}"),
        vcf  = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.vcf"),
        vcfp = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf"),
        inputvcfs=lambda wildcards, input: " ".join("-I {} ".format(in_) for in_ in input.vcfs),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__mutect2_merge", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" MergeVcfs "
        "    {params.inputvcfs} "
        "    -O {params.vcf} \n"
        "    > {log.out} 2> {log.err}\n"
        "head -n 10000 {params.vcf} |grep '^#' > {params.vcfp} 2>>{log.err} \n"
        "grep -v '^#' {params.vcf} | grep PASS >> {params.vcfp} 2>>{log.err} \n"
        "bgzip  {params.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        "bgzip  {params.vcfp} "
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcfp} "
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        "rm -rf {params.dir}/chroms"
        "    >> {log.out} 2>> {log.err}\n"


rule somatic_ss__octopus_split:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = temp(join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "chroms", "{sample}.{chr}.octopus.vcf")),
    params:
        dir = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "chroms", "{chr}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__octopus_split", "{sample}.{chr}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__octopus_split", "{sample}.{chr}.e"),
    threads:
        int(allocated("threads", "somatic_ss__octopus_split", cluster))
    container:
        config['container']['octopus']
    shell:
        "mkdir -p {params.dir}/tmp\n"
        "cd {params.dir}\n"
        "octopus "
        "    --threads {threads} "
        "    -C cancer "
        "    --working-directory {params.dir} "
        "    --temp-directory-prefix tmp "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.cram} "
        "    -o {output.vcf} "
        "    --forest-model /opt/octopus/resources/forests/germline.v0.7.4.forest "
        "    --somatic-forest-model /opt/octopus/resources/forests/somatic.v0.7.4.forest "
        "    --annotations AC AD DP "
        "    -T {wildcards.chr} "
        "    > {log.out} 2> {log.err}\n"

rule somatic_ss__octopus_merge:
    input:
        vcfs =lambda wildcards: \
             ["{}/42.somatic_ss_snvindel_octopus/{}/chroms/{}.{}.octopus.vcf".format(config['workdir'], wildcards.sample, wildcards.sample, chrid) for chrid in ['chr%d'%(i) for i in range(1,23)]+['chrX','chrY']],
    output:
        vcf  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.vcf.gz"),
        vcfp = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}"),
        vcf  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.vcf"),
        vcfp = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.pass.vcf"),
        inputvcfs=lambda wildcards, input: " ".join("-I {} ".format(in_) for in_ in input.vcfs),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__octopus_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__octopus_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__octopus_merge", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" MergeVcfs "
        "    {params.inputvcfs} "
        "    -O {params.vcf} \n"
        "    > {log.out} 2> {log.err}\n"
        "head -n 10000 {params.vcf} |grep '^#' > {params.vcfp} 2>>{log.err} \n"
        "grep -v '^#' {params.vcf} | grep PASS >> {params.vcfp} 2>>{log.err} \n"
        "bgzip  {params.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        "bgzip  {params.vcfp} "
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcfp} "
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        "rm -rf {params.dir}/chroms"
        "    >> {log.out} 2>> {log.err}\n"


rule somatic_ss__gripss:
    input:
        vcf = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.gridss.vcf"),
    output:
        vcf = join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}", "{sample}.gripss.filtered.vcf.gz"),
    params:
        dir = join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}"),
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
        "  -ref_genome_version 38 "
        "  -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "  -pon_sgl_file {config[references][hmftools]}/ref/38/sv/sgl_pon.38.bed.gz "
        "  -pon_sv_file {config[references][hmftools]}/ref/38/sv/sv_pon.38.bedpe.gz "
        "  -known_hotspot_file {config[references][hmftools]}/ref/38/sv/known_fusions.38.bedpe "
        "  -repeat_mask_file {config[references][hmftools]}/ref/38/sv/repeat_mask_data.38.fa.gz "
        "  -vcf {input.vcf} "
        "  -output_dir {params.dir}"
        "  > {log.out} 2> {log.err}\n"