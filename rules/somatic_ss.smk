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
        tbis =lambda wildcards: \
             ["{}/42.somatic_ss_snvindel_octopus/{}/chroms/{}.{}.octopus.vcf.gz.tbi".format(config['workdir'], wildcards.sample, wildcards.sample, chrid) for chrid in ['chr%d'%(i) for i in range(1,23)]+['chrX','chrY']],
    output:
        vcf  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.vcf.gz"),
        vcfp = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}"),
        inputvcfs=lambda wildcards, input: " ".join(" {} ".format(in_.replace('.tbi','')) for in_ in input.tbis),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__octopus_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__octopus_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__octopus_merge", cluster))
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
        "rm -rf {params.dir}/chroms"
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
    container:
        config['container']['vardict']
    shell:
        "cd {params.dir}\n"
        "grep ^chr {config[references][gatkbundle]}/scattered_calling_intervals/{wildcards.itv}/scattered.interval_list|"
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
        vcfs =lambda wildcards: \
             ["{}/43.somatic_ss_snvindel_vardict/{}/itvs/{}.{}.vardict.vcf".format(config['workdir'], wildcards.sample, wildcards.sample, itv) \
                for itv in ['temp_%.4d_of_%d'%(i,config['parameter']['gatkitv']) for i in range(1,config['parameter']['gatkitv']+1)]],
    output:
        vcf  = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.vcf.gz"),
        vcfp = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.pass.vcf.gz"),
    params:
        dir  = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}"),
        vcf  = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.vcf"),
        vcfp = join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.pass.vcf"),
        inputvcfs=lambda wildcards, input: " ".join("-I {} ".format(in_) for in_ in input.vcfs),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__vardict_merge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__vardict_merge", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__vardict_merge", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID\" MergeVcfs "
        "    {inputvcfs} "
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

rule somatic_ss__dellysv:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcfgz = join(config['workdir'], "45.somatic_ss_sv__delly", "{sample}", "{sample}.delly.vcf.gz"),
    params:
        vcf = join(config['workdir'], "45.somatic_ss_sv__delly", "{sample}", "{sample}.delly.vcf"),
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__dellysv", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__dellysv", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__dellysv", cluster))
    container:
        config['container']['delly']
    shell:
        "delly call "
        "    -g {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -x {config[references][delly]}/human.hg38.excl.tsv "
        "    {input.cram}"
        "    > {params.vcf} |"
        "bgzip > {output.vcfgz}"
        "    2> {log.err}\n"
        "tabix {output.vcfgz}"
        "  > {log.out} 2>> {log.err}\n"

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
    container:
        config['container']['amber']
    shell:
        "cd {params.dir}"
        "java -Xmx24G "
        "    -jar /usr/local/share/hmftools-amber-4.0-1/amber.jar "
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
    container:
        config['container']['cobalt']
    shell:
        "cd {params.dir}"
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
        vcf = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "M00001_TOD_BM.purple.sv.vcf.gz"),
    params:
        dir = join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}"),
        inputvcfs=lambda wildcards, input: " ".join("-I {} ".format(in_) for in_ in input.vcfs),
        gripssraw = lambda wildcards, input: input.gripss.replace('filtered.','')
    log:
        out = join(config['pipelinedir'], "logs", "somatic_ss__purple", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "somatic_ss__purple", "{sample}.e"),
    threads:
        int(allocated("threads", "somatic_ss__purple", cluster))
    container:
        config['container']['purple']
    shell:
        "cd {params.dir}"
        "java -Xmx24G "
        "    -jar /usr/local/share/hmftools-purple-4.0.2-0/purple.jar "
        "    -tumor M00001_TOD_BM "
        "    -amber $(dirname $(realpath {input.amber})) "
        "    -cobalt $(dirname $(realpath {input.cobalt})) "
        "    -gc_profile {config[references][hmftools]}/ref/38/copy_number/GC_profile.1000bp.38.cnp "
        "    -ref_genome {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta "
        "    -ref_genome_version 38 "
        "    -ensembl_data_dir {config[references][hmftools]}/ref/38/common/ensembl_data "
        "    -somatic_vcf {input.mutect} "
        "    -structural_vcf {input.gripss} "
        "    -sv_recovery_vcf {params.gripssraw} "
        "    -run_drivers "
        "    -driver_gene_panel {config[references][hmftools]}/ref/38/common/DriverGenePanel.38.tsv "
        "    -output_dir {output.vcf}"
        "  > {log.out} 2> {log.err}\n"
