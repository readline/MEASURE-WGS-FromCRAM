from scripts.utils import allocated,ignore

rule somatic_tn__mutect2_split:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
    output:
        vcf = temp(join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.mutect2.flt.vcf")),
    params:
        normal = lambda wildcards: dic_tumor_to_normal[wildcards.sample],
        vcf = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.mutect2.vcf"),
        f1r2 = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.f1r2.tar.gz"),
        rom = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.read-orientation-model.tar.gz"),
        pst = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.getpileupsummaries.table"),
        ct = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.calculatecontamination.table"),
        st = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.segments.table"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__mutect2_split", "{sample}.{chr}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__mutect2_split", "{sample}.{chr}.e"),
    threads:
        int(allocated("threads", "somatic_tn__mutect2_split", cluster))
    container:
        config['container']['gatk']
    shell:
        "normal=$(samtools view -H {input.cram0} | grep '^@RG'|  awk '/^@RG/ {{for (i=1; i<=NF; i++) if ($i ~ /^SM:/) print substr($i,4)}}'|head -1) 2> {log.err}\n"
        "echo Tumor: {wildcards.sample}.  Normal: {params.normal}, RG:$normal, 2>> {log.err}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" Mutect2 "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.cram} "
        "    -I {input.cram0} "
        "    -normal $normal"
        "    -L {wildcards.chr} "
        "    --panel-of-normals {config[references][pon]}/Mutect2.pon.vcf.gz "
        "    -O {params.vcf} "
        "    -germline-resource {config[references][gatksomatic]}/af-only-gnomad.hg38.vcf.gz "
        "    --f1r2-tar-gz {params.f1r2} "
        "    > {log.out} 2>> {log.err}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" LearnReadOrientationModel "
        "    -O {params.rom} "
        "    -I {params.f1r2}"
        "    >> {log.out} 2>> {log.err}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" GetPileupSummaries "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.cram} "
        "    -I {input.cram0} "
        "    -V {config[references][gatksomatic]}/small_exac_common_3.hg38.vcf.gz "
        "    -L {wildcards.chr} "
        "    -O {params.pst} "
        "    >> {log.out} 2>> {log.err}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" CalculateContamination "
        "    -I {params.pst} "
        "    -tumor-segmentation {params.st}  "
        "    -O {params.ct} "
        "    >> {log.out} 2>> {log.err}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" FilterMutectCalls "
        "    --reference {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -V {params.vcf} "
        "    --tumor-segmentation {params.st} "
        "    --ob-priors {params.rom} "
        "    -O {output.vcf} "
        "    >> {log.out} 2>> {log.err}\n"

rule somatic_tn__mutect2_merge:
    input:
        tbis =lambda wildcards: \
             ["{}/31.somatic_tn_snvindel_mutect2/{}/chroms/{}.{}.mutect2.flt.vcf.gz.tbi".format(config['workdir'], wildcards.sample, wildcards.sample, chrid) for chrid in ['chr%d'%(i) for i in range(1,23)]+['chrX','chrY']],
    output:
        vcf  = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "{sample}.mutect2.vcf.gz"),
        vcfp = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}"),
        normal = lambda wildcards: dic_tumor_to_normal[wildcards.sample],
        vcf0  = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "{sample}.mutect2.tmp.vcf.gz"),
        inputvcfs=lambda wildcards, input: " ".join(" {} ".format(in_.replace('.tbi','')) for in_ in input.tbis),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__mutect2_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__mutect2_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__mutect2_merge", cluster))
    container:
        config['container']['gatk']
    shell:
        "bcftools concat "
        "    -a "
        "    -O z "
        "    -o {params.vcf0} "
        "    {params.inputvcfs} "
        "    > {log.out} 2> {log.err}\n"
        "echo {wildcards.sample} > samplename.txt\n"
        "echo {params.normal} >> samplename.txt\n"
        "bcftools reheader "
        "    -s samplename.txt "
        "    -o {output.vcf} "
        "    {params.vcf0} "
        "    >> {log.out} 2>> {log.err}\n"
        "rm {params.vcf0} samplename.txt\n"
        "bcftools view "
        "    -f 'PASS' "
        "    {output.vcf}"
        "    -O z "
        "    -o {output.vcfp}"
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcfp} "
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        # "rm -rf {params.dir}/chroms"
        # "    >> {log.out} 2>> {log.err}\n"

rule somatic_tn__strelka:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
        manta = join(config['workdir'], "34.somatic_sv_manta", "{sample}", "results", "variants", "candidateSmallIndels.vcf.gz"),
    output:
        vgz = join(config['workdir'], "32.somatic_snv_strelka", "{sample}", "results", "variants", "somatic.snvs.vcf.gz"),
    params:
        dir = join(config['workdir'], "32.somatic_snv_strelka", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__strelka", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__strelka", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__strelka", cluster))
    container:
        config['container']['strelka']
    shell:
        "configureStrelkaSomaticWorkflow.py"
        "  --normalBam {input.cram0}"
        "  --tumorBam {input.cram}"
        "  --referenceFasta {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --indelCandidates {input.manta}"
        "  --runDir {params.dir}"
        "  --outputCallableRegions"
        "  > {log.out} 2> {log.err}\n"
        "cd {params.dir} \n"
        "./runWorkflow.py"
        "  -m local"
        "  -j {threads}"
        "  >> {log.out} 2>> {log.err}\n"
        "rm -rf workspace"
        "  >> {log.out} 2>> {log.err}\n"


rule somatic_tn__manta:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
    output:
        vgz = join(config['workdir'], "34.somatic_sv_manta", "{sample}", "results", "variants", "candidateSV.vcf.gz"),
        indel = join(config['workdir'], "34.somatic_sv_manta", "{sample}", "results", "variants", "candidateSmallIndels.vcf.gz"),
    params:
        dir = join(config['workdir'], "34.somatic_sv_manta", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__manta", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__manta", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__manta", cluster))
    container:
        config['container']['manta']
    shell:
        "configManta.py"
        "  --normalBam {input.cram0}"
        "  --tumorBam {input.cram}"
        "  --reference {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --runDir {params.dir}"
        "  > {log.out} 2> {log.err}\n"
        "cd {params.dir} \n"
        "./runWorkflow.py"
        "  -m local"
        "  -j {threads}"
        "  >> {log.out} 2>> {log.err}\n"
        "rm -rf workspace"
        "  >> {log.out} 2>> {log.err}\n"

rule somatic_tn__gridss_assemble_shards:
    input:
        cram1 = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
        ok1   = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss", "preprocess.ok"),
        ok0   = lambda wildcards: "{}/15.germline_sv_gridss/{}/_gridss/preprocess.ok".format(config['workdir'], dic_tumor_to_normal[wildcards.sample]),
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
