from scripts.utils import allocated,ignore

rule cram__check_md5:
    input:
        cram = lambda wildcards: dic_run[wildcards.run]['CRAM'],
    output:
        check   = join(config['workdir'], "01.cram", "md5check_run", "{run}.check"),
    log:
        out = join(config['pipelinedir'], "logs", "cram__check_md5", "{run}.o"),
        err = join(config['pipelinedir'], "logs", "cram__check_md5", "{run}.e"),
    params:
        md5 = lambda wildcards: dic_run[wildcards.run]['MD5'],
    threads:
        int(allocated("threads", "cram__check_md5", cluster))
    max_attempts: 3
    container:
        config['container']['samtools']
    shell:
        """
        # Check if the file exists
        if [ ! -f "{input.cram}" ]; then
            echo "Error: File '{input.cram}' not found."
            exit 1
        fi

        # Calculate the MD5 hash of the file
        actual_md5=$(md5sum {input.cram} | awk '{{ print $1 }}')

        # Compare the calculated MD5 hash with the expected one
        if [ "$actual_md5" = "{params.md5}" ]; then
            touch {output.check}
            echo "MD5 matched: $actual_md5" > {log.out}
        else
            echo "MD5 hash does not match the expected value." > {log.out}
            echo "Actural md5: $(actual_md5)" >> {log.out}
            echo "Listed md5:  {params.md5}" >> {log.out}
            exit 1
        fi
        """

rule cram_merge_sample:
    input:
        check = lambda wildcards: [join(config['workdir'], "01.cram", "md5check_run", "%s.check"%(run)) for run in dic_sample_to_runs[wildcards.sample]],
        crams = lambda wildcards: [dic_run[run]['CRAM'] for run in dic_sample_to_runs[wildcards.sample]],
    output:
        cram = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
        rgsm = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram.samples"),
    log: 
        out = join(config['pipelinedir'], "logs", "merge_sample_cram", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "merge_sample_cram", "{sample}.e"),
    threads:
        int(allocated("threads", "cram_merge_sample", cluster))
    max_attempts: 3
    container:
        config['container']['samtools']
    shell:
        """
        if [ $(ls -1 {input.crams} | wc -l) -eq 1 ]; then
            ln -s {input.crams} {output.cram} \
               > {log.out} 2> {log.err}
            samtools index -@ {threads} {output.cram}
            samtools samples {output.cram} | cut -f1 >> {output.rgsm}
        else
            samtools merge -@ {threads} --write-index -o {output.cram} {input.crams}  2> {log.err}
            samtools samples {output.cram} | cut -f1 >> {output.rgsm}
        fi
        """

rule cram__to_bam:
    input:
        cram   = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        bam   = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
        bai   = join(config['workdir'], "02.bam", "{sample}", "{sample}.bam.bai"),
    log: 
        out = join(config['pipelinedir'], "logs", "cram__to_bam", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cram__to_bam", "{sample}.e"),
    params:
        folder = join(config['workdir'], "02.bam", "{sample}"),
    threads:
        int(allocated("threads", "cram__to_bam", cluster))
    max_attempts: 3
    container:
        config['container']['samtools']
    shell:
        """
        cd {params.folder}
        samtools view -@ {threads} -O BAM -T {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta -o {output.bam} {input.cram} > {log.out} 2>{log.err}
        samtools index -@ {threads} {output.bam} >> {log.out} 2>>{log.err}
        """

rule cram__to_chrm:
    input:
        cram   = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        bam   = join(config['workdir'], "03.chrM", "{sample}", "{sample}.chrM.bam"),
        bai   = join(config['workdir'], "03.chrM", "{sample}", "{sample}.chrM.bam.bai"),
    log: 
        out = join(config['pipelinedir'], "logs", "cram__to_chrm", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cram__to_chrm", "{sample}.e"),
    params:
        folder = join(config['workdir'], "03.chrM", "{sample}"),
    threads:
        int(allocated("threads", "cram__to_chrm", cluster))
    max_attempts: 3
    container:
        config['container']['samtools']
    shell:
        """
        cd {params.folder}
        samtools view -@ {threads} --bam -o {output.bam} {input.cram} chrM > {log.out} 2>{log.err}
        samtools index -@ {threads} {output.bam} >> {log.out} 2>>{log.err}
        """


rule cram__verifybamid2:
    input:
        cram   = join(config['workdir'], "01.cram", "{sample}", "{sample}.cram"),
    output:
        selfsm   = join(config['workdir'], "04.verifybamid2", "{sample}", "{sample}.selfSM"),
    log: 
        out = join(config['pipelinedir'], "logs", "cram__verifybamid2", "{sample}.o"),
        err = join(config['pipelinedir'], "logs", "cram__verifybamid2", "{sample}.e"),
    params:
        folder = join(config['workdir'], "04.verifybamid2", "{sample}"),
        prefix = "{sample}"
    threads:
        int(allocated("threads", "cram__verifybamid2", cluster))
    max_attempts: 3
    container:
        config['container']['verifybamid2']
    shell:
        "cd {params.folder}\n"
        "verifybamid2"
        "  --SVDPrefix /usr/local/share/verifybamid2-2.0.1-10/resource/1000g.phase3.100k.b38.vcf.gz.dat"
        "  --BamFile {input.cram}"
        "  --Reference {config[references][gatkbundle]}/Homo_sapiens_assembly38.fasta"
        "  --Output {params.prefix}"
        "  --NumThread {threads}"
        "  >{log.out} 2>{log.err}"