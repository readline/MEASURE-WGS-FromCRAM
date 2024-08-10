from scripts.utils import allocated,ignore

itv50 = ['%.4d'%i for i in range(1,config['parameter']['gatkitv']+1)]

def generate_germline__gatk_hcmerge_input_gvcf(wildcards):
    return expand(join(config['workdir'], "10.germline_snv_gatk", "{sample}", "itvs", "{sample}.{itv}.vcf.gz"), sample=wildcards.sample, itv=itv50)

def generate_germline__gatk_hcmerge_input_bam(wildcards):
    return expand(join(config['workdir'], "10.germline_snv_gatk", "{sample}", "itvs", "{sample}.{itv}.vcf.bam"), sample=wildcards.sample, itv=itv50)


rule germline__deepvariant1:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        dv1a = temp( join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "gvcf.tfrecord-{itv}-of-%s.gz"%(itvdv[-1]) ) ),
        dv1b = temp( join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "make_examples.tfrecord-{itv}-of-%s.gz"%(itvdv[-1]) ) ),
        dv1c = temp( join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "make_examples.tfrecord-{itv}-of-%s.gz.example_info.json"%(itvdv[-1]) ) ),
    params:
        folder = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp"),
        itv="{itv}",
        dv1a = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "gvcf.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']) ),
        dv1b = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "make_examples.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']) ),
    log: 
        out = join(config['pipelinedir'], "logs", "germline__deepvariant1", "{sample}.{itv}.o"),
        err = join(config['pipelinedir'], "logs", "germline__deepvariant1", "{sample}.{itv}.e"),
    threads:
        int(allocated("threads", "germline__deepvariant1", cluster))
    container:
        config['container']['deepvariant']
    shell:
        "/opt/deepvariant/bin/make_examples"
        "  --mode calling"
        "  --ref {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --reads {input.cram}"
        "  --gvcf {params.dv1a}"
        "  --examples {params.dv1b}"
        "  --channels insert_size"
        "  --task {params.itv}"
        " > {log.out} 2> {log.err}\n"


rule germline__deepvariant2:
    input:
        dv1a = lambda wildcards: [join( config['workdir'], "11.germline_snv_deepvariant", wildcards.sample, "temp", 
                                  "gvcf.tfrecord-{itv}-of-{all}.gz".format(itv=i, all=itvdv[-1]) ) for i in itvdv[:-1]],
        dv1b = lambda wildcards: [join( config['workdir'], "11.germline_snv_deepvariant", wildcards.sample, "temp", 
                                  "make_examples.tfrecord-{itv}-of-{all}.gz".format(itv=i, all=itvdv[-1]) ) for i in itvdv[:-1]],
        dv1c = lambda wildcards: [join( config['workdir'], "11.germline_snv_deepvariant", wildcards.sample, "temp", 
                                  "make_examples.tfrecord-{itv}-of-{all}.gz.example_info.json".format(itv=i, all=itvdv[-1]) ) for i in itvdv[:-1]],
    output:
        status = temp( join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "step2.ok" ) ),
    params:
        call = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "call_variants_output.tfrecord.gz" ),
        folder = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp"),
        dv1b = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", 
                                  "make_examples.tfrecord@%d.gz" )%(config['parameter']['deepvariant_shards'] ),
        sif    = config['container']['deepvariant-gpu'],
        bind   = "/vf,/spin1,/data,/fdb,/gpfs,/lscratch,/home",
    log: 
        out = join(config['pipelinedir'], "logs", "germline__deepvariant2", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__deepvariant2", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__deepvariant2", cluster))
    shell:
        "singularity exec -B {params.bind} --nv {params.sif} "
        "/opt/deepvariant/bin/call_variants"
        "  --outfile {params.call}"
        "  --examples {params.dv1b}"
        "  --checkpoint /opt/models/wgs"
        " > {log.out} 2> {log.err}\n"
        "touch {output.status}"


rule germline__deepvariant3:
    input:
        dv1a = lambda wildcards: [join( config['workdir'], "11.germline_snv_deepvariant", wildcards.sample, "temp", 
                                  "gvcf.tfrecord-{itv}-of-{all}.gz".format(itv=i, all=itvdv[-1]) ) for i in itvdv[:-1]],
        status = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "step2.ok" ),
    output:
        vcf  = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf.gz" ),
        gvcf = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.gvcf.gz" ),
    params:
        folder = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp"),
        call   = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "call_variants_output.tfrecord.gz" ),
        vcf    = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf" ),
        gvcf   = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.gvcf" ),
        dv1a   = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "temp", "gvcf.tfrecord@%d.gz"%(config['parameter']['deepvariant_shards']) ),
        sif    = config['container']['deepvariant-gpu'],
        bind   = "/vf,/spin1,/data,/fdb,/gpfs,/lscratch,/home",
    log: 
        out = join(config['pipelinedir'], "logs", "germline__deepvariant3", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__deepvariant3", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__deepvariant3", cluster))
    shell:
        "singularity exec -B {params.bind} --nv {params.sif} "
        "/opt/deepvariant/bin/postprocess_variants"
        "    --ref {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "    --infile {params.call}"
        "    --nonvariant_site_tfrecord_path {params.dv1a}"
        "    --outfile {params.vcf}"
        "    --gvcf_outfile {params.gvcf}"
        " > {log.out} 2> {log.err}\n"
        "singularity exec -B {params.bind} --nv {params.sif} "
        "bgzip {params.vcf}"
        " >> {log.out} 2>> {log.err}\n"
        "singularity exec -B {params.bind} --nv {params.sif} "
        "bgzip {params.gvcf}"
        " >> {log.out} 2>> {log.err}\n"
        "singularity exec -B {params.bind} --nv {params.sif} "
        "tabix -p vcf {output.vcf}"
        " >> {log.out} 2>> {log.err}\n"
        "singularity exec -B {params.bind} --nv {params.sif} "
        "tabix -p vcf  {output.gvcf}"
        " >> {log.out} 2>> {log.err}\n"


rule germline__gatk_hcitv:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        gvcf = temp(join(config['workdir'], "10.germline_snv_gatk", "{sample}", "itvs", "{sample}.{itv}.vcf.gz")),
        bam  = temp(join(config['workdir'], "10.germline_snv_gatk", "{sample}", "itvs", "{sample}.{itv}.vcf.bam")),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gatk_hcitv", "{sample}.{itv}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gatk_hcitv", "{sample}.{itv}.e"),
    threads:
        int(allocated("threads", "germline__gatk_hcitv", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" "
        "    HaplotypeCaller" 
        "    --reference {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta" 
        "    --input {input.cram} "
        "    --output {output.gvcf} "
        "    --dbsnp {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
        "    --intervals {config[references][gatkbundle]}/scattered_calling_intervals/temp_{wildcards.itv}_of_50/scattered.interval_list"
        "    --bam-output {output.bam}"
        "    --tmp-dir /lscratch/$SLURM_JOB_ID/"
        "    -ERC GVCF"
        " >> {log.out} 2>> {log.err}\n"


rule germline__gatk_hcmerge:
    input:
        gvcf = generate_germline__gatk_hcmerge_input_gvcf,
        bam  = generate_germline__gatk_hcmerge_input_bam,
    params:
        tmpbam="/lscratch/$SLURM_JOB_ID/{sample}.vcf.bamout.bam",
        inputvcfs=lambda wildcards, input: " ".join("-I {} ".format(in_) for in_ in input.gvcf),
        inputbams=lambda wildcards, input: " ".join(" {} ".format(in_) for in_ in input.bam),
    output:
        gvcf = join(config['workdir'], "10.germline_snv_gatk", "{sample}", "{sample}.gvcf.gz"),
        gbam = join(config['workdir'], "10.germline_snv_gatk", "{sample}", "{sample}.gvcf.bam"),
        gmet = join(config['workdir'], "10.germline_snv_gatk", "{sample}", "{sample}.gvcf.gz.variant_calling_summary_metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gatk_hcmerge", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gatk_hcmerge", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__gatk_hcmerge", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" "
        "    MergeVcfs "
        "    -O {output.gvcf}"
        "    {params.inputvcfs}"
        " >> {log.out} 2>> {log.err}\n"
        "samtools merge"
        "    -@ {threads} "
        "    -o {params.tmpbam} "
        "    {params.inputbams}"
        " >> {log.out} 2>> {log.err}\n"
        "samtools sort -@ {threads}"
        "    -T /lscratch/$SLURM_JOB_ID/{wildcards.sample} "
        "    -o {output.gbam} "
        "    {params.tmpbam} "
        " >> {log.out} 2>> {log.err}\n"
        "gatk --java-options \"-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms20G -Xmx20G -XX:ParallelGCThreads=2\" "
        "    CollectVariantCallingMetrics"
        "    -I {output.gvcf}"
        "    -O {output.gvcf}"
        "    --DBSNP {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
        "    -SD {config[references][gatkbundlesup]}/Homo_sapiens_assembly38.dict"
        "    --GVCF_INPUT"
        "    --THREAD_COUNT {threads}"
        "    -TI {config[references][gatkbundle]}/wgs_calling_regions.hg38.interval_list"
        " >> {log.out} 2>> {log.err}\n"


rule germline__gdbimport:
    input:
        gvcf= lambda wildcards: ["{}/10.germline_snv_gatk/{}/{}.gvcf.gz".format(config['workdir'], sample, sample) for sample in dic_sample_to_runs],
    params:
        interval = "{}/scattered_calling_intervals/temp_{}_of_50/scattered.interval_list".format(config['references']['gatkbundle'], "{itv}"),
        inputvcfs=lambda wildcards, input: " ".join("-V {} ".format(in_) for in_ in input.gvcf),
    output:
        itvvcf = temp(join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.itv_{itv}.vcf.gz")),
        itvmf  = temp(join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.itv_{itv}.mf.vcf.gz")),
        itvso  = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.itv_{itv}.siteonly.vcf.gz"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gdbimport", "itv_{itv}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gdbimport", "itv_{itv}.e"),
    threads:
        int(allocated("threads", "germline__gdbimport", cluster))
    container:
        config['container']['gatk']
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenomicsDBImport \
            {params.inputvcfs} \
            -L {params.interval} \
            --genomicsdb-workspace-path /lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv} \
            --merge-input-intervals \
            --consolidate >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenotypeGVCFs \
            -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta \
            -O {output.itvvcf} \
            -D {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
            -G StandardAnnotation -G AS_StandardAnnotation \
            --allow-old-rms-mapping-quality-annotation-data \
            --merge-input-intervals \
            -V gendb:///lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv} \
            -L {params.interval}   >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantFiltration \
            --filter-expression "ExcessHet>54.69" \
            --filter-name ExcessHet \
            -V {output.itvvcf} \
            -O {output.itvmf} >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            MakeSitesOnlyVcf \
             -I {output.itvmf} \
             -O {output.itvso} >> {log.out} 2>> {log.err}
        """
            
            
rule germline__genotyping:
    input:
        mfitvs=lambda wildcards: \
             ["{}/10.germline_snv_gatk/VQSR/Merge.itv_{}.mf.vcf.gz".format(config['workdir'], '%.4d'%itv) for itv in range(1, config['parameter']['gatkitv']+1)],
        soitvs=lambda wildcards: \
             ["{}/10.germline_snv_gatk/VQSR/Merge.itv_{}.siteonly.vcf.gz".format(config['workdir'], '%.4d'%itv) for itv in range(1, config['parameter']['gatkitv']+1)],
    output:
        sovcf = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.sites_only.vcf.gz"),
        mfvcf = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.mf.vcf.gz"),
        mir   = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.indels.recal"),
        mit   = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.indels.tranches"),
        msr   = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.snps.recal"),
        mst   = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.snps.tranches"),
        msmr  = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "Merge.snps.model.report"),
        mirv  = join(config['workdir'], "10.germline_snv_gatk", "VQSR", "tmp.indel.recalibrated.vcf"),
        vqsr  = join(config['workdir'], "10.germline_snv_gatk", "Merge.flt.vqsr.vcf.gz"),
    params:
        mfinputs=lambda wildcards, input: " ".join("--input {} ".format(in_) for in_ in input.mfitvs),
        soinputs=lambda wildcards, input: " ".join("--input {} ".format(in_) for in_ in input.soitvs),
    log:
        out = join(config['pipelinedir'], "logs", "germline__genotyping.o"),
        err = join(config['pipelinedir'], "logs", "germline__genotyping.e"),
    threads:
        int(allocated("threads", "germline__genotyping", cluster))
    container:
        config['container']['gatk']
    shell:
        """
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GatherVcfsCloud \
          --ignore-safety-checks \
          --gather-type BLOCK \
          --output {output.sovcf} \
          {params.soinputs}  >> {log.out} 2>> {log.err}
        gatk --java-options "-Xmx16g -Xms16g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
          GatherVcfsCloud \
          --ignore-safety-checks \
          --gather-type BLOCK \
          --output {output.mfvcf} \
          {params.mfinputs}  >> {log.out} 2>> {log.err}
        tabix -p vcf {output.sovcf}
        tabix -p vcf {output.mfvcf}
        gatk --java-options -Xms24g \
          VariantRecalibrator \
          -V {output.sovcf} \
          -O {output.mir} \
          --tranches-file {output.mit} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
          -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
          -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP\
          --use-allele-specific-annotations \
          -mode INDEL \
          --max-gaussians 4 \
          --resource:mills,known=false,training=true,truth=true,prior=12 {config[references][gatkbundle]}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
          --resource:axiomPoly,known=false,training=true,truth=false,prior=10 {config[references][gatkbundle]}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
          --resource:dbsnp,known=true,training=false,truth=false,prior=2 {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz >> {log.out} 2>> {log.err}
        gatk --java-options -Xms50g \
          VariantRecalibrator \
          -V {output.sovcf} \
          -O {output.msr} \
          --tranches-file {output.mst} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 \
          -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
          -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
          --use-allele-specific-annotations \
          -mode SNP \
          --sample-every-Nth-variant 10 \
          --output-model {output.msmr} \
          --max-gaussians 6 \
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {config[references][gatkbundle]}/hapmap_3.3.hg38.vcf.gz \
          -resource:omni,known=false,training=true,truth=true,prior=12 {config[references][gatkbundle]}/1000G_omni2.5.hg38.vcf.gz \
          -resource:1000G,known=false,training=true,truth=false,prior=10 {config[references][gatkbundle]}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz >> {log.out} 2>> {log.err}
        gatk --java-options -Xms80g \
          VariantRecalibrator \
          -V {output.sovcf} \
          -O {output.msr} \
          --tranches-file {output.mst} \
          --trust-all-polymorphic \
          -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 \
          -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
          -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
          --use-allele-specific-annotations \
          -mode SNP \
          --input-model {output.msmr} \
          --max-gaussians 6 \
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {config[references][gatkbundle]}/hapmap_3.3.hg38.vcf.gz \
          -resource:omni,known=false,training=true,truth=true,prior=12 {config[references][gatkbundle]}/1000G_omni2.5.hg38.vcf.gz \
          -resource:1000G,known=false,training=true,truth=false,prior=10 {config[references][gatkbundle]}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {config[references][gatkbundle]}/Homo_sapiens_assembly38.dbsnp138.vcf.gz >> {log.out} 2>> {log.err}
        gatk --java-options -Xms5g \
          ApplyVQSR \
          -O {output.mirv} \
          -V {output.mfvcf} \
          --recal-file {output.mir} \
          --use-allele-specific-annotations \
          --tranches-file {output.mit} \
          --truth-sensitivity-filter-level 95.0 \
          --create-output-variant-index true \
          -mode INDEL >> {log.out} 2>> {log.err}
        gatk --java-options -Xms5g \
          ApplyVQSR \
          -O {output.vqsr} \
          -V {output.mirv} \
          --recal-file {output.msr} \
          --use-allele-specific-annotations \
          --tranches-file {output.mst} \
          --truth-sensitivity-filter-level 99.7 \
          --create-output-variant-index true \
          -mode SNP >> {log.out} 2>> {log.err}
        """

rule germline__peddy:
    input:
        vcf  = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf.gz" ),
        rgsm = join( config['workdir'], "01.cram", "{sample}", "{sample}.cram.samples" ),
    output:
        selfsm = join(config['workdir'], "05.peddy", "{sample}", "{sample}.peddy.ped"),
        selfsc = join(config['workdir'], "05.peddy", "{sample}", "{sample}.sex_check.csv"),
    log: 
        out = join(config['pipelinedir'], "logs", "germline__peddy", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__peddy", "{sample}.e"),
    params:
        pid = lambda wildcards: dic_sample_to_patient[wildcards.sample],
        sid = "{sample}",
        folder = join(config['workdir'], "05.peddy", "{sample}"),
        prefix = "{sample}",
    threads:
        int(allocated("threads", "germline__peddy", cluster))
    container:
        config['container']['peddy']
    shell:
        "cd {params.folder}\n"
        "sm=$(cat {input.rgsm}|head -n 1)\n"
        "echo '{params.pid}\t'$sm'\t0\t0\t0\t-9' > peddy.ped\n"
        "peddy"
        "  -p {threads}"
        "  --plot"
        "  --prefix {params.sid}"
        "  {input.vcf}"
        "  peddy.ped"
        " >> {log.out} 2>> {log.err}\n"


rule germline__strelka:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        manta = join(config['workdir'], "13.germline_sv_manta", "{sample}", "results", "variants", "candidateSmallIndels.vcf.gz"),
    output:
        vgz = join(config['workdir'], "12.germline_snv_strelka", "{sample}", "results", "variants", "variants.vcf.gz"),
    params:
        dir = join(config['workdir'], "12.germline_snv_strelka", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__strelka", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__strelka", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__strelka", cluster))
    container:
        config['container']['strelka']
    shell:
        "configureStrelkaGermlineWorkflow.py"
        "  --bam {input.cram}"
        "  --referenceFasta {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --indelCandidates {input.manta}"
        "  --runDir {params.dir}"
        "  > {log.out} 2> {log.err}\n"
        "cd {params.dir} \n"
        "./runWorkflow.py"
        "  -m local"
        "  -j {threads}"
        "  >> {log.out} 2>> {log.err}\n"
        "rm -rf workspace"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__manta:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vgz = join(config['workdir'], "13.germline_sv_manta", "{sample}", "results", "variants", "candidateSV.vcf.gz"),
        indel = join(config['workdir'], "13.germline_sv_manta", "{sample}", "results", "variants", "candidateSmallIndels.vcf.gz"),
    params:
        dir = join(config['workdir'], "13.germline_sv_manta", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__manta", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__manta", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__manta", cluster))
    container:
        config['container']['manta']
    shell:
        "configManta.py"
        "  --bam {input.cram}"
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


rule germline__tiddit:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = join(config['workdir'], "14.germline_sv_tiddit", "{sample}", "{sample}.tiddit.vcf"),
    params:
        dir    = join(config['workdir'], "14.germline_sv_tiddit", "{sample}"),
        prefix = join(config['workdir'], "14.germline_sv_tiddit", "{sample}", "{sample}.tiddit"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__tiddit", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__tiddit", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__tiddit", cluster))
    container:
        config['container']['tiddit']
    shell:
        "cd {params.dir}\n"
        "rm -rf *\n"
        "tiddit"
        "  --sv"
        "  --threads {threads}"
        "  --bam {input.cram}"
        "  --ref {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -o {params.prefix}"
        "  > {log.out} 2> {log.err}\n"


rule germline__gridss_preprocess:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        ok = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss", "preprocess.ok"),
    params:
        workspace = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gridss_preprocess", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gridss_preprocess", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__gridss_preprocess", cluster))
    container:
        config['container']['gridss']
    shell:
        "cd {params.workspace} > {log.out} 2> {log.err}\n"
        "gridss "
        "  -s preprocess "
        "  -r {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -t {threads} "
        "  -w {params.workspace} "
        "  {input.cram}"
        "  -b {config[references][blacklist]}/Merge.exclude.bed"
        "  >> {log.out} 2>> {log.err}\n"
        "touch {output.ok}"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__gridss_assemble:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        ok   = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss", "preprocess.ok"),
    output:
        ok = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss", "assemble_{shard}.ok"),
    params:
        workspace = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss"),
        bam       = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.assemble.bam"),
        shard     = "{shard}"
    log:
        out = join(config['pipelinedir'], "logs", "germline__gridss_assemble", "{sample}.{shard}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gridss_assemble", "{sample}.{shard}.e"),
    threads:
        int(allocated("threads", "germline__gridss_assemble", cluster))
    container:
        config['container']['gridss']
    shell:
        "gridss "
        "  -s assemble "
        "  -r {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -t {threads}"
        "  -a {params.bam}"
        "  --jobnodes {config[parameter][gridss_shards]}"
        "  --jobindex {params.shard}"
        "  {input.cram}"
        "  -w {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"
        "touch {output.ok}"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.workspace}"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__gridss_call:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        ok =   lambda wildcards: [ join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss", "assemble_{shard}.ok".format(shard=i)) for i in range(config['parameter']['gridss_shards'])]
    output:
        bam = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.assemble.bam"),
        vcf = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.gridss.vcf"),
    params:
        workspace = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "_gridss"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gridss_call", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gridss_call", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__gridss_call", cluster))
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
        "  {input.cram}"
        "  > {log.out} 2> {log.err}\n"

rule germline__gridss_virusbreakend:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = join(config['workdir'], "16.germline_sv_virusbreakend", "{sample}", "{sample}.virusbreakend.vcf"),
    params:
        dir = join(config['workdir'], "16.germline_sv_virusbreakend", "{sample}", "tmp"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gridss_virusbreakend", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gridss_virusbreakend", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__gridss_virusbreakend", cluster))
    container:
        config['container']['gridss']
    shell:
        "mkdir -p {params.dir} \n"
        "cd {params.dir} \n"
        "ln -s {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta \n"
        "virusbreakend "
        "  -r Homo_sapiens_assembly38.fasta "
        "  -j /usr/local/share/gridss-2.13.2-3/gridss.jar "
        "  -o {output.vcf} "
        "  --db {config[references][virusbreakend]} "
        "  {input.cram}"
        "  > {log.out} 2> {log.err}\n"
        "cd .. \n"
        "rm -rf tmp"

rule germline__gripss_germline:
    input:
        vcf = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.gridss.vcf"),
    output:
        vcf = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "gripss_germline", "{sample}.gripss.filtered.vcf.gz"),
    params:
        dir = join(config['workdir'], "15.germline_sv_gridss", "{sample}", "gripss_germline"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__gripss_germline", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__gripss_germline", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__gripss_germline", cluster))
    container:
        config['container']['gripss']
    shell:
        "java -Xmx32g -jar /usr/local/share/hmftools-gripss-2.4-0/gripss.jar "
        "  -germline "
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

rule germline__canvas:
    input:
        bam    = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
        selfsc = join(config['workdir'], "05.peddy", "{sample}", "{sample}.sex_check.csv"),
        vgz    = join(config['workdir'], "12.germline_snv_strelka", "{sample}", "results", "variants", "variants.vcf.gz"),
    output:
        vcf = join(config['workdir'], "17.germline_cnv_canvas", "{sample}", "CNV.vcf.gz"),
    params:
        script = join(config['pipelinedir'], "scripts", "peddy2ploidy.py"),
        prefix = join(config['workdir'], "17.germline_cnv_canvas", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__canvas", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__canvas", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__canvas", cluster))
    container:
        config['container']['canvas']
    shell:
        "cd {params.prefix}\n"
        "rm -rf *\n"
        "mkdir ref\n"
        "cd ref\n"
        "ln -s {config[references][canvas]} canvas\n"
        "ln -s {config[references][canvas]}/Sequence Sequence\n"
        "cd {params.prefix}\n"
        "python3 {params.script} {input.selfsc} ref/ploidy.vcf"
        "  > {log.out} 2> {log.err}\n"
        "Canvas.sh"
        "  SmallPedigree-WGS"
        "  -b {input.bam}"
        "  -o {params.prefix}"
        "  -r ref/canvas/kmer.fa"
        "  -g ref/canvas/WholeGenomeFasta"
        "  -f ref/canvas/filter13.bed"
        "  --sample-b-allele-vcf={input.vgz}"
        "  --ploidy-vcf=ref/ploidy.vcf"
        "  >> {log.out} 2>> {log.err}\n"
        "rm -rf ref TempCNV"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__melt_ins:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Insert", "LINE1.final_comp.vcf"),
    params:
        prefix = join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Insert"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__melt_ins", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__melt_ins", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__melt_ins", cluster))
    container:
        config['container']['melt']
    shell:
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Single "
        "  -a "
        "  -c {threads} "
        "  -h {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -bamfile {input.cram} "
        "  -n /opt/MELT/add_bed_files/Hg38/Hg38.genes.bed "
        "  -t /opt/MELT/me_refs/transposon_file_list.txt "
        "  -w {params.prefix}"
        "  > {log.out} 2> {log.err}\n"


rule germline__melt_del1:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_LINE1", "DEL.final_comp.vcf"),
    params:
        prefix = join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_LINE1"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__melt_del1", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__melt_del1", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__melt_del1", cluster))
    container:
        config['container']['melt']
    shell:
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Genotype"
        "  -bamfile {input.cram} "
        "  -w {params.prefix}/workspace"
        "  -bed /opt/MELT/add_bed_files/Hg38/LINE1.deletion.bed "
        "  -h {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.prefix}/workspace/*.del.tsv > /lscratch/$SLURM_JOB_ID/del_LINE1.output.list"
        "  2>> {log.err}\n"
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Merge"
        "  -mergelist /lscratch/$SLURM_JOB_ID/del_LINE1.output.list"
        "  -bed /opt/MELT/add_bed_files/Hg38/LINE1.deletion.bed"
        "  -h {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -o {params.prefix}"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__melt_del2:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        vcf = join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_AluY", "DEL.final_comp.vcf"),
    params:
        prefix = join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_AluY"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__melt_del2", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__melt_del2", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__melt_del2", cluster))
    container:
        config['container']['melt']
    shell:
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Genotype"
        "  -bamfile {input.cram} "
        "  -w {params.prefix}/workspace"
        "  -bed /opt/MELT/add_bed_files/Hg38/AluY.deletion.bed "
        "  -h {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.prefix}/workspace/*.del.tsv > /lscratch/$SLURM_JOB_ID/del_AluY.output.list"
        "  2>> {log.err}\n"
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Merge"
        "  -mergelist /lscratch/$SLURM_JOB_ID/del_AluY.output.list"
        "  -bed /opt/MELT/add_bed_files/Hg38/AluY.deletion.bed"
        "  -h {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -o {params.prefix}"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__msi_msisensorpro:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        prefix = join(config['workdir'], "23.germline_msi_msisensorpro", "{sample}", "{sample}.msisensorpro"),
    log:
        out = join(config['pipelinedir'], "logs", "germline__msi_msisensorpro", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__msi_msisensorpro", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__msi_msisensorpro", cluster))
    container:
        config['container']['msisensor-pro']
    shell:
        "msisensor-pro "
        "  pro"
        "  -d {config[references][pon]}/MSIsensorpro.workspace/baseline/ref_scan_baseline"
        "  -g {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  -t {input.cram} "
        "  -o {output.prefix}"
        "  -b {threads}"
        "  >> {log.out} 2>> {log.err}\n"


rule germline__hlala:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        txt = join(config['workdir'], "25.germline_hla-la", "{sample}", "hla", "summaryStatistics.txt"),
    params: 
        sample  = "{sample}",
        outdir = join(config['workdir'], "25.germline_hla-la"),
        graph  = config['references']['hla-la'],
        sif    = config['container']['hla-la'],
        bind   = "/vf,/spin1,/data,/fdb,/gpfs,/lscratch,/home,%s/graphs:/usr/local/opt/hla-la/graphs"%(config['references']['hla-la'])
    log:
        out = join(config['pipelinedir'], "logs", "germline__hlala", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "germline__hlala", "{sample}.e"),
    threads:
        int(allocated("threads", "germline__hlala", cluster))
    shell:
        "singularity exec -B {params.bind} {params.sif} "
        "  HLA-LA.pl "
        "  --BAM {input.cram} "
        "  --graph PRG_MHC_GRCh38_withIMGT "
        "  --sampleID {params.sample} "
        "  --maxThreads {threads} "
        "  --workingDir {params.outdir}"
        "  >> {log.out} 2>> {log.err}\n"

