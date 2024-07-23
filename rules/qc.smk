from scripts.utils import allocated,ignore

rule cramqc__flagstat:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        flagstat = join(config['workdir'], "01.cram", "{sample}", "QC", "Flagstat", "Flagstat.json"),
    log: 
        out = join(config['pipelinedir'], "logs", "cramqc__flagstat", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__flagstat", "{sample}.e"),
    params:
        folder = join(config['workdir'], "01.cram", "{sample}", "QC", "Flagstat",),
    threads:
        int(allocated("threads", "cramqc__flagstat", cluster))
    container:
        config['container']['samtools']
    shell:
        """
        cd {params.folder}
        samtools flagstat -@ {threads} -O json {input.cram} > {output.flagstat} 2>{log.err}
        """

rule cramqc__collect_quality_yield_metrics:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        metrix = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "QualityYield.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_quality_yield_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_quality_yield_metrics", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__collect_quality_yield_metrics", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms2000m -Xmx3000m\" "
        "  CollectQualityYieldMetrics"
        "  -I {input.cram}"
        "  -O {output.metrix}"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  > {log.out} 2> {log.err}\n"

rule cramqc__collect_wgs_metrics:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        metrix   = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "WGS.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_wgs_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_wgs_metrics", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__collect_wgs_metrics", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms2000m -Xmx3000m\" "
        " CollectWgsMetrics"
        "  -I {input.cram}"
        "  -O {output.metrix}"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --INCLUDE_BQ_HISTOGRAM true"
        "  --INTERVALS {config[references][gatkbundle]}/wgs_coverage_regions.hg38.interval_list"
        "  --VALIDATION_STRINGENCY SILENT"
        "  --USE_FAST_ALGORITHM true"
        "  --CREATE_MD5_FILE true"
        "  --READ_LENGTH 151"
        "  > {log.out} 2> {log.err}\n"

rule cramqc__collect_all_reads_multiple_metrics:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        status   = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "AllReadsMultiple.ok"),
    params:
        prefix = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "AllReadsMultiple.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_all_reads_multiple_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_all_reads_multiple_metrics", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__collect_all_reads_multiple_metrics", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "  CollectMultipleMetrics"
        "  -I {input.cram}"
        "  -O {params.prefix}"
        "  --ASSUME_SORTED true"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --VALIDATION_STRINGENCY LENIENT"
        "  --PROGRAM null "
        "  --PROGRAM CollectBaseDistributionByCycle "
        "  --PROGRAM CollectInsertSizeMetrics "
        "  --PROGRAM MeanQualityByCycle "
        "  --PROGRAM QualityScoreDistribution "
        "  --METRIC_ACCUMULATION_LEVEL null "
        "  --METRIC_ACCUMULATION_LEVEL ALL_READS "
        " > {log.out} 2> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__collect_read_groups_multiple_metrics:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        status   = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "ReadGroupsMultiple.ok"),
    params:
        prefix = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "ReadGroupsMultiple.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_read_groups_multiple_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_read_groups_multiple_metrics", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__collect_read_groups_multiple_metrics", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "  CollectMultipleMetrics"
        "  -I {input.cram}"
        "  -O {params.prefix}"
        "  --ASSUME_SORTED true"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --VALIDATION_STRINGENCY LENIENT"
        "  --PROGRAM null "
        "  --PROGRAM CollectBaseDistributionByCycle "
        "  --PROGRAM CollectInsertSizeMetrics "
        "  --PROGRAM CollectAlignmentSummaryMetrics "
        "  --PROGRAM QualityScoreDistribution "
        "  --METRIC_ACCUMULATION_LEVEL ALL_READS "
        "  --METRIC_ACCUMULATION_LEVEL READ_GROUP "
        " > {log.out} 2> {log.err}\n"
        "touch {output.status}"
        "  >> {log.out} 2>> {log.err}\n"

rule cramqc__collect_aggregation_metrics:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        status   = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "SM_LB_Aggregation.ok"),
    params:
        prefix = join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "SM_LB_Aggregation.metrics"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__collect_aggregation_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__collect_aggregation_metrics", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__collect_aggregation_metrics", cluster))
    container:
        config['container']['gatk']
    shell:
        "gatk --java-options \"-Xms5000m -Xmx6500m\" "
        "  CollectMultipleMetrics"
        "  -I {input.cram}"
        "  -O {params.prefix}"
        "  --ASSUME_SORTED true"
        "  --CREATE_MD5_FILE true"
        "  -R {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --PROGRAM null "
        "  --PROGRAM CollectAlignmentSummaryMetrics "
        "  --PROGRAM CollectInsertSizeMetrics "
        "  --PROGRAM CollectSequencingArtifactMetrics "
        "  --PROGRAM QualityScoreDistribution "
        "  --METRIC_ACCUMULATION_LEVEL null "
        "  --METRIC_ACCUMULATION_LEVEL SAMPLE "
        "  --METRIC_ACCUMULATION_LEVEL LIBRARY "
        " >> {log.out} 2>> {log.err}\n"
        "touch {output.status}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__mosdepth:
    input:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        status = join(config['workdir'], "01.cram", "{sample}", "QC", "Mosdepth", "{sample}.mosdepth.summary.txt"),
    params:
        prefix = join(config['workdir'], "01.cram", "{sample}", "QC", "Mosdepth", "{sample}"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__mosdepth", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__mosdepth", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__mosdepth", cluster))
    container:
        config['container']['mosdepth']
    shell:
        "export MOSDEPTH_Q0=NO_COVERAGE \n"
        "export MOSDEPTH_Q1=LOW_COVERAGE \n"
        "export MOSDEPTH_Q2=CALLABLE \n"
        "export MOSDEPTH_Q3=HIGH_COVERAGE \n"
        "export MOSDEPTH_Q4=HIGH_COVERAGE_300 \n"
        "export MOSDEPTH_Q5=HIGH_COVERAGE_1000 \n"
        "mosdepth"
        "    --threads {threads}"
        "    --fasta {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "    --quantize 0:1:5:150:300:1000:"
        "    {params.prefix}"
        "    {input.cram}"
        " >> {log.out} 2>> {log.err}\n"

rule cramqc__summary_metrics:
    input:
        join(config['workdir'], "01.cram", "{sample}", "QC", "Flagstat", "Flagstat.json"),
        join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "QualityYield.metrics"),
        join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "WGS.metrics"),
        join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "AllReadsMultiple.ok"),
        join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "ReadGroupsMultiple.ok"),
        join(config['workdir'], "01.cram", "{sample}", "QC", "Metrics", "SM_LB_Aggregation.ok"),
        join(config['workdir'], "01.cram", "{sample}", "QC", "Mosdepth", "{sample}.mosdepth.summary.txt"),
    output:
        status   = join(config['workdir'], "01.cram", "{sample}", "QC", "Summary.ok"),
    log:
        out = join(config['pipelinedir'], "logs", "cramqc__summary_metrics", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cramqc__summary_metrics", "{sample}.e"),
    threads:
        int(allocated("threads", "cramqc__summary_metrics", cluster))
    shell:
        "touch {output.status}"
