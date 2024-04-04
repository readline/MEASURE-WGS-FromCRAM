from scripts.utils import allocated,ignore

rule cramqc__flagstat:
    input:
        cram     = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram"),
    output:
        flagstat = join(config['workdir'], "01.cram", "{sm}", "QC_flagstat", "flagstat"),
    log: 
        out = join(config['pipelinedir'], "logs", "cram__to_bam", "{sm}.o"),
        err = join(config['pipelinedir'], "logs", "cram__to_bam", "{sm}.e"),
    params:
        folder = join(config['workdir'], "02.bam", "{sm}"),
    threads:
        int(allocated("threads", "cram__to_bam", cluster))
    container:
        config['container']['samtools']
    shell:
        """
        cd {params.folder}
        samtools view -@ {threads} --bam -o {output.bam} {input.bam} > {log.out} 2>{log.err}
        samtools index -@ {threads} {output.bam} >> {log.out} 2>>{log.err}
        """