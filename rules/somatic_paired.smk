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
    max_attempts: 3
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
        sample_sid = lambda wildcards: dic_sample_to_sid[wildcards.sample],
        normal_sid = lambda wildcards: dic_sample_to_sid[dic_tumor_to_normal[wildcards.sample]],
        vcf0  = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "{sample}.mutect2.tmp.vcf.gz"),
        inputvcfs=lambda wildcards, input: " ".join(" {} ".format(in_.replace('.tbi','')) for in_ in input.tbis),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__mutect2_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__mutect2_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__mutect2_merge", cluster))
    max_attempts: 3
    container:
        config['container']['gatk']
    shell:
        "cd {params.dir} \n"
        "bcftools concat "
        "    -a "
        "    -O z "
        "    -o {params.vcf0} "
        "    {params.inputvcfs} "
        "    > {log.out} 2> {log.err}\n"
        "echo {params.sample_sid} {wildcards.sample} > samplename.txt\n"
        "echo {params.normal_sid} {params.normal} >> samplename.txt\n"
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
        "rm -rf {params.dir}/chroms"
        "    >> {log.out} 2>> {log.err}\n"

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
    max_attempts: 3
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


rule somatic_tn__muse:
    input:
        bam = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
        bam0 = lambda wildcards: "{}/02.bam/{}/{}.bam".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]), 
    output:
        vcf = join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}.MuSE.vcf"),
    params:
        dir = join(config['workdir'], "33.somatic_snv_muse", "{sample}"),
        prefix = join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__muse", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__muse", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__muse", cluster))
    max_attempts: 3
    container:
        config['container']['muse']
    shell:
        "cd {params.dir}\n"
        "ln -s {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz "
        " > {log.out} 2> {log.err}\n"
        "tabix -p vcf Homo_sapiens_assembly38.dbsnp138.vcf.gz "
        " >> {log.out} 2>> {log.err}\n"
        "/MuSE/bin/MuSE call "
        "    -f {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -O {params.prefix} "
        "    -n {threads}"
        "    {input.bam} "
        "    {input.bam0} "
        " >> {log.out} 2>> {log.err}\n"
        "/MuSE/bin/MuSE sump "
        "    -I {params.prefix}.MuSE.txt "
        "    -n {threads}"
        "    -G "
        "    -D Homo_sapiens_assembly38.dbsnp138.vcf.gz "
        "    -O {params.prefix}.MuSE.vcf"
        " >> {log.out} 2>> {log.err} \n"
        "rm Homo_sapiens_assembly38.dbsnp138.vcf.gz"
        " > {log.out} 2> {log.err}\n"

rule somatic_tn__muse_post:
    input:
        vcf = join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}.MuSE.vcf.gz.tbi"),
    output:
        vcfp  = join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}.MuSE.pass.vcf.gz"),
    params:
        dir = join(config['workdir'], "33.somatic_snv_muse", "{sample}"),
        vcf = join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}.MuSE.vcf.gz"),
        vcf0  = join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}.MuSE.tmp.vcf.gz"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__muse_post", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__muse_post", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__muse_post", cluster))
    max_attempts: 3
    container:
        config['container']['gatk']
    shell:
        "echo TUMOR {wildcards.sample} > samplename.txt\n"
        "echo NORMAL {params.normal} >> samplename.txt\n"
        "bcftools reheader "
        "    -s samplename.txt "
        "    -o {params.vcf0} "
        "    {params.vcf} "
        "    >> {log.out} 2>> {log.err} \n"
        "mv {params.vcf0} {params.vcf} \n"
        "rm samplename.txt \n"
        "bcftools view "
        "    -f 'PASS' "
        "    {params.vcf}"
        "    -O z "
        "    -o {output.vcfp}"
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcfp} "
        "    >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {params.vcf} "
        "    >> {log.out} 2>> {log.err}\n"
        
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
    max_attempts: 3
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
    max_attempts: 3
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
    max_attempts: 3
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
        int(allocated("threads", "somatic_tn__gripss", cluster))
    max_attempts: 3
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

rule somatic_tn__dellysv:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
    output:
        prebcf = join(config['workdir'], "36.somatic_tn_sv__delly", "{sample}", "{sample}.delly.pre.bcf"),
    params:
        dir = join(config['workdir'], "36.somatic_tn_sv__delly", "{sample}"),
        callbcf = join(config['workdir'], "36.somatic_tn_sv__delly", "{sample}", "{sample}.delly.bcf"),
        prebcf = join(config['workdir'], "36.somatic_tn_sv__delly", "{sample}", "{sample}.delly.pre.bcf"),
        genobcf = join(config['workdir'], "36.somatic_tn_sv__delly", "{sample}", "{sample}.delly.geno.bcf"),
        samplelist = join(config['workdir'], "36.somatic_tn_sv__delly", "{sample}", "sample.list"),
        tumorsid  = lambda wildcards: dic_sample_to_sid[wildcards.sample],
        normalsid = lambda wildcards: dic_sample_to_sid[dic_tumor_to_normal[wildcards.sample]],
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__dellysv", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__dellysv", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__dellysv", cluster))
    max_attempts: 3
    container:
        config['container']['delly']
    shell:
        "cd {params.dir} \n"
        "delly call "
        "    -g {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -x {config[references][delly]}/human.hg38.excl.tsv "
        "    -o {params.callbcf} "
        "    {input.cram} "
        "    {input.cram0}"
        "  > {log.out} 2> {log.err}\n"
        "echo '{params.tumorsid}\ttumor' > {params.samplelist}\n"
        "echo '{params.normalsid}\tcontrol' >> {params.samplelist}\n"
           "delly filter "
        "    -f somatic "
        "    -o {output.prebcf} "
        "    -s {params.samplelist} "
        "    {params.callbcf} "
        "  >> {log.out} 2>> {log.err}\n"

rule somatic_tn__dellysv_joint:
    input:
        crams = lambda wildcards: 
            ["{}/01.cram/{}/{}.cram".format(config['workdir'], sample, sample) for sample in [i[0] for i in tasks['somatic_paired']]],
        crams0 = lambda wildcards: 
            ["{}/01.cram/{}/{}.cram".format(config['workdir'], sample, sample) for sample in [i[1] for i in tasks['somatic_paired']]],
        bcfs = lambda wildcards: 
            ["{}/36.somatic_tn_sv__delly/{}/{}.delly.pre.bcf".format(config['workdir'], sample, sample) for sample in [i[0] for i in tasks['somatic_paired']]],
    output:
        premerge = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.pre.bcf"),
        genobcf = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.geno.bcf"),
    params:
        dir = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__dellysv2", "Merge.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__dellysv2", "Merge.e"),
    threads:
        int(allocated("threads", "somatic_tn__dellysv_joint", cluster))
    max_attempts: 3
    container:
        config['container']['delly']
    shell:
        "cd {params.dir} \n"
        "ls {input.crams} | while read n; do\n"
        "  echo $(samtools samples $n | cut -f1)'\ttumor'\n"
        "done > sample.list \n"
        "ls {input.crams0} | while read n; do\n"
        "  echo $(samtools samples $n | cut -f1)'\tcontrol'\n"
        "done >> sample.list \n"
        "ls {config[references][pon]}/cram/*.cram | while read n; do\n"
        "  echo $(samtools samples $n | cut -f1)'\tcontrol'\n"
        "done >> sample.list \n"
        "bcftools merge "
        "    -m id "
        "    -O b "
        "    --force-single "
        "    -o {output.premerge} "
        "    {input.bcfs}"
        "  > {log.out} 2> {log.err}\n"
        "delly call "
        "    -g {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -v {output.premerge} "
        "    -o {output.genobcf} "
        "    -x {config[references][delly]}/human.hg38.excl.tsv "
        "    {input.crams} "
        "    {input.crams0} "
        "    $(ls {config[references][pon]}/cram/*.cram)"
        "  >> {log.out} 2>> {log.err}\n"


rule somatic_tn__dellysv_post:
    input:
        genobcf = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.geno.bcf"),
    output:
        bcf = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.somatic.bcf"),
        vcf = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.somatic.vcf.gz"),
        bcfg = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.germline.bcf"),
        vcfg = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.germline.vcf.gz"),
        gvcf = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls", "Merge.delly.gvcf.gz"),
    params:
        dir = join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__dellysv_post", "Post.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__dellysv_post", "Post.e"),
    threads:
        int(allocated("threads", "somatic_tn__dellysv_post", cluster))
    max_attempts: 3
    container:
        config['container']['delly']
    shell:
        "cd {params.dir} \n"
        "delly filter "
        "    -f somatic "
        "    -o {output.bcf} "
        "    -s sample.list "
        "    {input.genobcf}"
        "  > {log.out} 2> {log.err}\n"
        "delly filter "
        "    -f germline "
        "    -o {output.bcfg} "
        "    -s sample.list "
        "    {input.genobcf}"
        "  >> {log.out} 2>> {log.err}\n"
        "bcftools view {output.bcf} | bgzip > {output.vcf}"
        "  2>> {log.err}\n"
        "tabix -p vcf {output.vcf}"
        "  >> {log.out} 2>> {log.err}\n"
        "bcftools view {output.bcfg} | bgzip > {output.vcfg}"
        "  2>> {log.err}\n"
        "tabix -p vcf {output.vcfg}"
        "  >> {log.out} 2>> {log.err}\n"
        "bcftools view {input.genobcf} | bgzip > {output.gvcf}"
        "  2>> {log.err}\n"
        "tabix -p vcf {output.gvcf}"
        "  >> {log.out} 2>> {log.err}\n"

rule somatic_tn__amber:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
    output:
        ok = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "amber", "amber.ok"),
    params:
        dir = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "amber"),
        normal = lambda wildcards: dic_tumor_to_normal[wildcards.sample],
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__amber", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__amber", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__amber", cluster))
    max_attempts: 3
    container:
        config['container']['amber']
    shell:
        "cd {params.dir} \n"
        "java -Xmx24G "
        "    -jar /usr/local/share/hmftools-amber-4.0.1-0/amber.jar "
        "    -tumor {wildcards.sample} "
        "    -tumor_bam {input.cram} "
        "    -reference {params.normal}"
        "    -reference_bam {input.cram0}"
        "    -output_dir {params.dir} "
        "    -threads {threads} "
        "    -loci {config[references][hmftools]}/ref/38/copy_number/AmberGermlineSites.38.tsv.gz "
        "    -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -ref_genome_version 38"
        "  > {log.out} 2> {log.err}\n"
        "touch {output.ok}"

rule somatic_tn__cobalt:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        cram0 = lambda wildcards: "{}/01.cram/{}/{}.cram".format(config['workdir'], dic_tumor_to_normal[wildcards.sample], dic_tumor_to_normal[wildcards.sample]),
    output:
        ok = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "cobalt", "cobalt.ok"),
    params:
        dir = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "cobalt"),
        normal = lambda wildcards: dic_tumor_to_normal[wildcards.sample],
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__cobalt", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__cobalt", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__cobalt", cluster))
    max_attempts: 3
    container:
        config['container']['cobalt']
    shell:
        "cd {params.dir} \n"
        "java -Xmx12G "
        "    -jar /usr/local/share/hmftools-cobalt-1.16-0/cobalt.jar "
        "    -tumor {wildcards.sample} "
        "    -tumor_bam {input.cram} "
        "    -reference {params.normal}"
        "    -reference_bam {input.cram0}"
        "    -output_dir {params.dir} "
        "    -threads {threads} "
        "    -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -gc_profile {config[references][hmftools]}/ref/38/copy_number/GC_profile.1000bp.38.cnp"
        "  > {log.out} 2> {log.err}\n"
        "touch {output.ok}"

rule somatic_tn__purple:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        amber = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "amber", "amber.ok"),
        cobalt = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "cobalt", "cobalt.ok"),
        mutect = join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
        gripss = join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.gripss.filtered.vcf.gz"),
    output:
        cnv = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "{sample}.purple.cnv.somatic.tsv"),
    params:
        dir = join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}"),
        normal = lambda wildcards: dic_tumor_to_normal[wildcards.sample],
        gripssraw = lambda wildcards, input: input.gripss.replace('filtered.','')
    log:
        out = join(config['pipelinedir'], "logs", "somatic_tn__purple", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_tn__purple", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_tn__purple", cluster))
    max_attempts: 3
    container:
        config['container']['purple']
    shell:
        "cd {params.dir}\n"
        "java -Xmx24G "
        "    -jar /usr/local/share/hmftools-purple-4.0.2-0/purple.jar "
        "    -tumor {wildcards.sample} "
        "    -reference {params.normal}"
        "    -amber $(dirname $(realpath {input.amber})) "
        "    -cobalt $(dirname $(realpath {input.cobalt})) "
        "    -gc_profile {config[references][hmftools]}/ref/38/copy_number/GC_profile.1000bp.38.cnp "
        "    -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -ref_genome_version 38 "
        "    -ensembl_data_dir {config[references][hmftools]}/ref/38/common/ensembl_data "
        "    -somatic_hotspots {config[references][hmftools]}/ref/38/variants/KnownHotspots.somatic.38.vcf.gz "
        "    -somatic_vcf {input.mutect} "
        "    -somatic_sv_vcf {input.gripss} "
        "    -sv_recovery_vcf {params.gripssraw} "
        "    -driver_gene_panel {config[references][hmftools]}/ref/38/common/DriverGenePanel.38.tsv "
        "    -output_dir {params.dir}"
        "  > {log.out} 2> {log.err}\n"