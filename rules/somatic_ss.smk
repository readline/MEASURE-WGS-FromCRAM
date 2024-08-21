from scripts.utils import allocated,ignore

rule somatic_ss__mutect2_split:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = temp(join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.mutect2.flt.vcf")),
    params:
        vcf = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.mutect2.vcf"),
        f1r2 = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.f1r2.tar.gz"),
        rom = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.read-orientation-model.tar.gz"),
        pst = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.getpileupsummaries.table"),
        ct = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.calculatecontamination.table"),
        st = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "chroms", "{sample}.{chr}.segments.table"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_split", "{sample}.{chr}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_split", "{sample}.{chr}.e"),
    threads:
        int(allocated("threads", "somatic_ss__mutect2_split", cluster))
    max_attempts: 3
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
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" LearnReadOrientationModel "
        "    -O {params.rom} "
        "    -I {params.f1r2}"
        "    >> {log.out} 2>> {log.err}\n"
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" GetPileupSummaries "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.cram} "
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

rule somatic_ss__mutect2_merge:
    input:
        tbis =lambda wildcards: \
             ["{}/41.somatic_ss_snvindel_mutect2/{}/chroms/{}.{}.mutect2.flt.vcf.gz.tbi".format(config['workdir'], wildcards.sample, wildcards.sample, chrid) for chrid in ['chr%d'%(i) for i in range(1,23)]+['chrX','chrY']],
    output:
        vcf  = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.vcf.gz"),
        vcfp = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}"),
        vcf0  = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.tmp.vcf.gz"),
        inputvcfs=lambda wildcards, input: " ".join(" {} ".format(in_.replace('.tbi','')) for in_ in input.tbis),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__mutect2_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__mutect2_merge", cluster))
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
        "echo {wildcards.sample} > samplename.txt\n"
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


rule somatic_ss__octopus_split:
    input:
        bam = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
    output:
        vcf = temp(join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "itvs", "{sample}.{itv}.octopus.vcf")),
    params:
        dir = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "itvs", "{itv}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__octopus_split", "{sample}.{itv}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__octopus_split", "{sample}.{itv}.e"),
    threads:
        int(allocated("threads", "somatic_ss__octopus_split", cluster))
    max_attempts: 3
    container:
        config['container']['octopus']
    shell:
        "mkdir -p {params.dir}/tmp\n"
        "cd {params.dir}\n"
        "grep ^chr {config[references][vardictitv]}/{wildcards.itv}/scattered.interval_list|"
        "    cut -f1-3 > {wildcards.itv}.bed"
        "    2> {log.err}\n"
        "octopus "
        "    --threads {threads} "
        "    -C cancer "
        "    --working-directory {params.dir} "
        "    --temp-directory-prefix tmp "
        "    -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -I {input.bam} "
        "    -o {output.vcf} "
        "    --forest-model /opt/octopus/resources/forests/germline.v0.7.4.forest "
        "    --somatic-forest-model /opt/octopus/resources/forests/somatic.v0.7.4.forest "
        "    --annotations AC AD DP "
        "    -t {wildcards.itv}.bed "
        "    > {log.out} 2> {log.err}\n"

rule somatic_ss__octopus_merge:
    input:
        tbis =lambda wildcards: \
             ["{}/42.somatic_ss_snvindel_octopus/{}/itvs/{}.{}.octopus.vcf.gz.tbi".format(config['workdir'], wildcards.sample, wildcards.sample, itv) \
             for itv in ['temp_%.4d_of_%d'%(i,config['parameter']['vardict']) for i in range(1,config['parameter']['vardict']+1)]],
    output:
        vcf  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.vcf.gz"),
        vcfp = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}"),
        vcf0 = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.tmp.vcf.gz"),
        inputvcfs=lambda wildcards, input: " ".join(" {} ".format(in_.replace('.tbi','')) for in_ in input.tbis),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__octopus_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__octopus_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__octopus_merge", cluster))
    max_attempts: 3
    container:
        config['container']['gatk']
    shell:
        "bcftools concat "
        "    -a "
        "    -O z "
        "    -o {params.vcf} "
        "    {params.inputvcfs} "
        "    > {log.out} 2> {log.err}\n"
        "echo {wildcards.sample} > samplename.txt\n"
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
        "rm -rf {params.dir}/itvs"
        "    >> {log.out} 2>> {log.err}\n"

rule somatic_ss__vardict_split:
    input:
        bam = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
    output:
        vcf = temp(join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "itvs", "{sample}.{itv}.vardict.vcf")),
    params:
        dir = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "itvs"),
        temp= "/lscratch/$SLURM_JOB_ID/{sample}.{itv}.txt"
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__vardict_split", "{sample}.{itv}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__vardict_split", "{sample}.{itv}.e"),
    threads:
        int(allocated("threads", "somatic_ss__vardict_split", cluster))
    max_attempts: 3
    container:
        config['container']['vardict']
    shell:
        "cd {params.dir}\n"
        "grep ^chr {config[references][vardictitv]}/{wildcards.itv}/scattered.interval_list|"
        "    cut -f1-3 > {wildcards.itv}.bed"
        "    2> {log.err}\n"
        "export JAVA_OPTS='\"-Xms90g\" \"-Xmx90g\"'"
        "    > {log.out} 2>> {log.err}\n"
        "vardict-java "
        "    -G {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -f 0.05 "
        "    -N {wildcards.sample} "
        "    -b {input.bam} "
        "    -th {threads} "
        "    -c 1 "
        "    -S 2 "
        "    -E 3 "
        "    {wildcards.itv}.bed |"
        "teststrandbias.R "
        "    -N {wildcards.sample} "
        "    -E "
        "    -f 0.05 > {params.temp} "
        "    2>> {log.err}\n"
        "var2vcf_valid.pl "
        "    -N {wildcards.sample} "
        "    -E "
        "    -f 0.05 "
        "    {params.temp} "
        "    > {output.vcf} "
        "    2>> {log.err}\n"

rule somatic_ss__vardict_merge:
    input:
        tbis =lambda wildcards: \
             ["{}/43.somatic_ss_snvindel_vardict/{}/itvs/{}.{}.vardict.vcf.gz.tbi".format(config['workdir'], wildcards.sample, wildcards.sample, itv) \
                for itv in ['temp_%.4d_of_%d'%(i,config['parameter']['vardict']) for i in range(1,config['parameter']['vardict']+1)]],
    output:
        vcf  = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.vcf.gz"),
        vcfp = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}"),
        inputvcfs=lambda wildcards, input: " ".join(" {} ".format(in_.replace('.tbi','')) for in_ in input.tbis),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__vardict_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__vardict_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__vardict_merge", cluster))
    max_attempts: 3
    container:
        config['container']['gatk']
    shell:
        "bcftools concat "
        "    -a "
        "    -O z "
        "    -o {output.vcf} "
        "    {params.inputvcfs} "
        "    > {log.out} 2> {log.err}\n"
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
        "rm -rf {params.dir}/itvs"
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
    max_attempts: 3
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

rule somatic_ss__gridss_annot:
    input:
        vcf = join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}", "{sample}.gripss.filtered.vcf.gz"),
        deepvariant = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.pass.vcf.gz" ),
    output:
        annot = join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}", "AnnotSV", "{sample}.gripss.pass.tsv"),
    params:
        annotdir = join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}", "AnnotSV"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__gridss_annot", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__gridss_annot", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__gridss_annot", cluster))
    max_attempts: 3
    container:
        config['container']['annotsv']
    shell:
        "AnnotSV "
        "   -SVinputFile {input.vcf} "
        "   -annotationsDir {config[references][annotsv]} "
        "   -bedtools bedtools "
        "   -bcftools bcftools "
        "   -annotationMode full "
        "   -genomeBuild GRCh38 "
        "   -includeCI 1 "
        "   -overwrite 1 "
        "   -outputFile {wildcards.sample}.gripss.pass "
        "   -outputDir {params.annotdir} "
        "   -SVinputInfo 1 "
        "   -SVminSize 50 "
        "   -overlap 70 "
        "   -snvIndelFiles {input.deepvariant} "
        "   -snvIndelPASS 1 "
        "  > {log.out} 2> {log.err}\n"

rule somatic_ss__dellysv:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        bcf = join(config['workdir'], "45.somatic_ss_sv__delly", "{sample}", "{sample}.delly.bcf"),
        vcfgz = join(config['workdir'], "45.somatic_ss_sv__delly", "{sample}", "{sample}.delly.vcf.gz"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__dellysv", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__dellysv", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__dellysv", cluster))
    max_attempts: 3
    container:
        config['container']['delly']
    shell:
        "delly call "
        "    -g {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -x {config[references][delly]}/human.hg38.excl.tsv "
        "    -o {output.bcf} "
        "    {input.cram} "
        "  > {log.out} 2> {log.err}\n"
        "bcftools view -O z -W -o {output.vcfgz} {output.bcf}"
        "  >> {log.out} 2>> {log.err}\n"

rule somatic_ss__dellysv_joint:
    input:
        crams = lambda wildcards: 
            ["{}/01.cram/{}/{}.cram".format(config['workdir'], sample, sample) for sample in tasks['somatic_ss']+[i[0] for i in tasks['somatic_paired']]],
        bcfs = lambda wildcards: 
            ["{}/45.somatic_ss_sv__delly/{}/{}.delly.bcf".format(config['workdir'], sample, sample) for sample in tasks['somatic_ss']+[i[0] for i in tasks['somatic_paired']]],
    output:
        premerge = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.pre.bcf"),
        genobcf = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.geno.bcf"),
    params:
        dir = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__dellysv_int", "Merge.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__dellysv_int", "Merge.e"),
    threads:
        int(allocated("threads", "somatic_ss__dellysv_joint", cluster))
    max_attempts: 3
    container:
        config['container']['delly']
    shell:
        "cd {params.dir} \n"
        "ls {input.crams} | while read n; do\n"
        "  echo $(samtools samples $n | cut -f1)'\ttumor'\n"
        "done > sample.list \n"
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
        "    $(ls {config[references][pon]}/cram/*.cram)"
        "  >> {log.out} 2>> {log.err}\n"

rule somatic_ss__dellysv_post:
    input:
        genobcf = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.geno.bcf"),
    output:
        bcf = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.somatic.bcf"),
        vcf = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.somatic.vcf.gz"),
        bcfg = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.germline.bcf"),
        vcfg = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.germline.vcf.gz"),
        gvcf = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.gvcf.gz"),
    params:
        dir = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls"),
        samplelist = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "sample.list"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__dellysv_post", "Post.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__dellysv_post", "Post.e"),
    threads:
        int(allocated("threads", "somatic_ss__dellysv_post", cluster))
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
        "bcftools view -O z -W -o {output.vcf} {output.bcf}"
        "  >> {log.out} 2>> {log.err}\n"
        "bcftools view -O z -W -o {output.vcfg} {output.bcfg}"
        "  >> {log.out} 2>> {log.err}\n"
        "bcftools view -O z -W -o {output.gvcf} {input.genobcf}"
        "  >> {log.out} 2>> {log.err}\n"

rule somatic_ss__dellysv_annot:
    input:
        vcf = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.somatic.vcf.gz"),
        vcfg = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls", "Merge.delly.germline.vcf.gz"),
    output:
        annot = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls",  "AnnotSV", "Merge.delly.somatic.tsv"),
        annotg = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls",  "AnnotSV", "Merge.delly.germline.tsv"),
    params:
        annotdir = join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls",  "AnnotSV"),
        annotfile1= "Merge.delly.somatic",
        annotfile2= "Merge.delly.germline",
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__dellysv_annot", "AnnotSV.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__dellysv_annot", "AnnotSV.e"),
    threads:
        int(allocated("threads", "somatic_ss__dellysv_annot", cluster))
    max_attempts: 3
    container:
        config['container']['annotsv']
    shell:
        "AnnotSV "
        "   -SVinputFile {input.vcf} "
        "   -annotationsDir {config[references][annotsv]} "
        "   -bedtools bedtools "
        "   -bcftools bcftools "
        "   -annotationMode full "
        "   -genomeBuild GRCh38 "
        "   -includeCI 1 "
        "   -overwrite 1 "
        "   -outputFile {params.annotfile1} "
        "   -outputDir {params.annotdir} "
        "   -SVinputInfo 1 "
        "   -SVminSize 50 "
        "   -overlap 70 "
        "  > {log.out} 2> {log.err}\n"
        "AnnotSV "
        "   -SVinputFile {input.vcfg} "
        "   -annotationsDir {config[references][annotsv]} "
        "   -bedtools bedtools "
        "   -bcftools bcftools "
        "   -annotationMode full "
        "   -genomeBuild GRCh38 "
        "   -includeCI 1 "
        "   -overwrite 1 "
        "   -outputFile {params.annotfile2} "
        "   -outputDir {params.annotdir} "
        "   -SVinputInfo 1 "
        "   -SVminSize 50 "
        "   -overlap 70 "
        "  >> {log.out} 2>> {log.err}\n"


rule somatic_ss__cnvkit:
    input:
        bam = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
        mosdepth = join(config['workdir'], "01.cram", "{sample}", "QC", "Mosdepth", "{sample}.mosdepth.summary.txt"),
    output:
        cns = join(config['workdir'], "46.somatic_ss_cnv__cnvkit", "{sample}", "{sample}.call.cns"),
        ok = join(config['workdir'], "46.somatic_ss_cnv__cnvkit", "{sample}", "cnvkit.ok"),
    params:
        dir = join(config['workdir'], "46.somatic_ss_cnv__cnvkit", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__cnvkit", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__cnvkit", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__cnvkit", cluster))
    max_attempts: 3
    container:
        config['container']['cnvkit']
    shell:
        "cd {params.dir} \n"
        "genderinfo=$(python {config[pipelinedir]}/scripts/getgender.py {input.mosdepth}) \n"
        "cnvkit.py "
        "   batch {input.bam} "
        "   -m wgs "
        "   -r {config[references][pon]}/CNVkit.pon.cnn "
        "   --output-dir {params.dir} "
        "   --scatter "
        "   -p {threads} "
        "   $genderinfo "
        "  > {log.out} 2> {log.err}\n"
        "touch {output.ok}"

rule somatic_ss__amber:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        ok = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "amber", "amber.ok"),
    params:
        dir = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "amber"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__amber", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__amber", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__amber", cluster))
    max_attempts: 3
    container:
        config['container']['amber']
    shell:
        "cd {params.dir} \n"
        "java -Xmx24G "
        "    -jar /usr/local/share/hmftools-amber-4.0.1-0/amber.jar "
        "    -tumor {wildcards.sample} "
        "    -tumor_bam {input.cram} "
        "    -output_dir {params.dir} "
        "    -threads {threads} "
        "    -loci {config[references][hmftools]}/ref/38/copy_number/AmberGermlineSites.38.tsv.gz "
        "    -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -ref_genome_version 38"
        "  > {log.out} 2> {log.err}\n"
        "touch {output.ok}"

rule somatic_ss__cobalt:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        ok = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "cobalt", "cobalt.ok"),
    params:
        dir = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "cobalt"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__cobalt", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__cobalt", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__cobalt", cluster))
    max_attempts: 3
    container:
        config['container']['cobalt']
    shell:
        "cd {params.dir} \n"
        "java -Xmx12G "
        "    -jar /usr/local/share/hmftools-cobalt-1.16-0/cobalt.jar "
        "    -tumor {wildcards.sample} "
        "    -tumor_bam {input.cram} "
        "    -output_dir {params.dir} "
        "    -threads {threads} "
        "    -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -tumor_only_diploid_bed {config[references][hmftools]}/ref/38/copy_number/DiploidRegions.38.bed.gz "
        "    -gc_profile {config[references][hmftools]}/ref/38/copy_number/GC_profile.1000bp.38.cnp"
        "  > {log.out} 2> {log.err}\n"
        "touch {output.ok}"

rule somatic_ss__purple:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        amber = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "amber", "amber.ok"),
        cobalt = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "cobalt", "cobalt.ok"),
        mutect = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
        gripss = join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}", "{sample}.gripss.filtered.vcf.gz"),
    output:
        cnv = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "{sample}.purple.cnv.somatic.tsv"),
    params:
        dir = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}"),
        gripssraw = lambda wildcards, input: input.gripss.replace('filtered.','')
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__purple", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__purple", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__purple", cluster))
    max_attempts: 3
    container:
        config['container']['purple']
    shell:
        "cd {params.dir}\n"
        "java -Xmx24G "
        "    -jar /usr/local/share/hmftools-purple-4.0.2-0/purple.jar "
        "    -tumor {wildcards.sample} "
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

rule somatic_ss__canvas:
    input:
        selfsc = join(config['workdir'], "05.peddy", "{sample}", "{sample}.sex_check.csv"),
        bam = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
        vcfp = join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
        vcf  = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf.gz" ),
    output:
        vcf = join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}", "CNV.vcf.gz"),
    params:
        dir = join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__canvas", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__canvas", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__canvas", cluster))
    max_attempts: 3
    container:
        config['container']['canvas']
    shell: 
        "cd {params.dir} \n"
        "rm -rf * \n"
        "mkdir ref \n"
        "cd ref \n"
        "ln -s {config[references][canvas]} canvas\n"
        "ln -s {config[references][canvas]}/Sequence Sequence\n"
        "cd {params.dir}\n"
        "python3 {config[pipelinedir]}/scripts/peddy2ploidy.py {input.selfsc} ref/ploidy.vcf"
        "  > {log.out} 2> {log.err}\n"
        "Canvas.sh"
        "    Somatic-WGS"
        "    -b {input.bam}"
        "    -n {wildcards.sample} "
        "    -o {params.dir}"
        "    -r ref/canvas/kmer.fa"
        "    -g ref/canvas/WholeGenomeFasta"
        "    -f ref/canvas/filter13.bed"
        "    --sample-b-allele-vcf={input.vcf}"
        "    --ploidy-vcf=ref/ploidy.vcf"
        "    >> {log.out} 2>> {log.err}\n"
        "rm -rf ref TempCNV*"
        "    >> {log.out} 2>> {log.err}\n"

rule somatic_ss__canvas_annot:
    input:
        vcf = join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}", "CNV.vcf.gz"),
        deepvariant = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.pass.vcf.gz" ),
    output:
        vcfp = join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}", "CNV.pass.vcf.gz"),
        annot = join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}", "AnnotSV", "Canvas.CNV.pass.tsv"),
    params:
        annotdir = join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}", "AnnotSV"),
        annotfile1= "Canvas.CNV.pass",
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__canvas_annot", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__canvas_annot", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__canvas_annot", cluster))
    max_attempts: 3
    container:
        config['container']['annotsv']
    shell:
        "bcftools view -O z -W -f PASS -o {output.vcfp} {input.vcf}"
        "AnnotSV "
        "   -SVinputFile {output.vcfp} "
        "   -annotationsDir {config[references][annotsv]} "
        "   -bedtools bedtools "
        "   -bcftools bcftools "
        "   -annotationMode full "
        "   -genomeBuild GRCh38 "
        "   -includeCI 1 "
        "   -overwrite 1 "
        "   -outputFile {params.annotfile1} "
        "   -outputDir {params.annotdir} "
        "   -SVinputInfo 1 "
        "   -SVminSize 50 "
        "   -overlap 70 "
        "   -snvIndelFiles {input.deepvariant} "
        "   -snvIndelPASS 1 "
        "  > {log.out} 2> {log.err}\n"