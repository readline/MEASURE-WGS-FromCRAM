from scripts.utils import allocated,ignore

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
        "  --ref {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "    --ref {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
        "    --infile {params.call}"
        "    --nonvariant_site_tfrecord_path {params.dv1a}"
        "    --outfile {params.vcf}"
        "    --gvcf_outfile {params.gvcf}"
        " > {log.out} 2> {log.err}\n"
        "bgzip {params.vcf}"
        " >> {log.out} 2>> {log.err}\n"
        "bgzip {params.gvcf}"
        " >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf {output.vcf}"
        " >> {log.out} 2>> {log.err}\n"
        "tabix -p vcf  {output.gvcf}"
        " >> {log.out} 2>> {log.err}\n"


rule germline__peddy:
    input:
        vcf  = join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf.gz" ),
        rgsm = join( config['workdir'], "01.cram", "{sample}", "{sample}.cram.samples" ),
    output:
        selfsm = join(config['workdir'], "05.peddy", "{sample}", "{sample}.peddy.ped"),
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
        "  --referenceFasta {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "  --reference {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "  --ref {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "  -r {config[references][gridss][path]}/Homo_sapiens_assembly38.fasta"
        "  -t {threads} "
        "  -w {params.workspace} "
        "  {input.cram}"
        "  -b {config[references][blacklist][path]}/Merge.exclude.bed"
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
        "  -r {config[references][gridss][path]}/Homo_sapiens_assembly38.fasta"
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
        "  -r {config[references][gridss][path]}/Homo_sapiens_assembly38.fasta"
        "  -t {threads}"
        "  -a {output.bam}"
        "  -o {output.vcf}"
        "  -w {params.workspace}"
        "  {input.cram}"
        "  > {log.out} 2> {log.err}\n"


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
        "  -h {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "  -h {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.prefix}/workspace/*.del.tsv > /lscratch/$SLURM_JOB_ID/del_LINE1.output.list"
        "  2>> {log.err}\n"
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Merge"
        "  -mergelist /lscratch/$SLURM_JOB_ID/del_LINE1.output.list"
        "  -bed /opt/MELT/add_bed_files/Hg38/LINE1.deletion.bed"
        "  -h {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "  -h {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
        "  >> {log.out} 2>> {log.err}\n"
        "ls {params.prefix}/workspace/*.del.tsv > /lscratch/$SLURM_JOB_ID/del_AluY.output.list"
        "  2>> {log.err}\n"
        "java -Xmx68g -jar /opt/MELT/MELT.jar "
        "  Deletion-Merge"
        "  -mergelist /lscratch/$SLURM_JOB_ID/del_AluY.output.list"
        "  -bed /opt/MELT/add_bed_files/Hg38/AluY.deletion.bed"
        "  -h {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        "  -d {config[references][PanelOfNormal][path]}/wgs_nygc112.v1/msisensorpro/Homo_sapiens_assembly38.fasta.scan.list_baseline "
        "  -g {config[references][gatkbundle][path]}/Homo_sapiens_assembly38.fasta"
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
        graph  = config['references']['hla-la']['path'],
        sif    = config['container']['hla-la'],
        bind   = "/vf,/spin1,/data,/fdb,/gpfs,/lscratch,/home,%s/graphs:/usr/local/opt/hla-la/graphs"%(config['references']['hla-la']['path'])
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

