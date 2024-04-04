from scripts.utils import allocated,ignore

rule cram__check_md5:
    input:
        cram = lambda wildcards: samples[wildcards.sm]['CRAM'],
    output:
        cram   = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram"),
        crai   = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram.crai"),
        md5    = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram.pass.md5"),
        rgsm   = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram.rg_sm"),
    log: 
        out = join(config['pipelinedir'], "logs", "cram__check_md5", "{sm}.o"),
        err = join(config['pipelinedir'], "logs", "cram__check_md5", "{sm}.e"),
    params:
        folder = join(config['workdir'], "01.cram", "{sm}"),
        md5 = lambda wildcards: samples[wildcards.sm]['MD5'],
    threads:
        int(allocated("threads", "cram__check_md5", cluster))
    container:
        config['container']['samtools']
    shell:
        """
        cd {params.folder}
        # Check if the file exists
        if [ ! -f "{input.cram}" ]; then
            echo "Error: File '{input.cram}' not found."
            exit 1
        fi

        # Calculate the MD5 hash of the file
        actual_md5=$(md5sum {input.cram} | awk '{{ print $1 }}')

        # Compare the calculated MD5 hash with the expected one
        if [ "$actual_md5" = "{params.md5}" ]; then
            ln -s {input.cram} {output.cram}
            samtools index -@ {threads} {output.cram}
            echo "$actual_md5" > {output.md5}
            samtools samples {output.cram} | cut -f1 > {output.rgsm}
        else
            echo "MD5 hash does not match the expected value."
            echo "Actural md5: $(actual_md5)"
            echo "Listed md5:  {params.md5}"
            exit 1
        fi
        """

rule cram__to_bam:
    input:
        cram   = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram"),
    output:
        bam   = join(config['workdir'], "02.bam", "{sm}", "{sm}.bam"),
        bai   = join(config['workdir'], "02.bam", "{sm}", "{sm}.bam.bai"),
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
        samtools view -@ {threads} --bam -o {output.bam} {input.cram} > {log.out} 2>{log.err}
        samtools index -@ {threads} {output.bam} >> {log.out} 2>>{log.err}
        """

rule cram__to_chrm:
    input:
        cram   = join(config['workdir'], "01.cram", "{sm}", "{sm}.cram"),
    output:
        bam   = join(config['workdir'], "03.chrM", "{sm}", "{sm}.chrM.bam"),
        bai   = join(config['workdir'], "03.chrM", "{sm}", "{sm}.chrM.bam.bai"),
    log: 
        out = join(config['pipelinedir'], "logs", "cram__to_chrm", "{sm}.o"),
        err = join(config['pipelinedir'], "logs", "cram__to_chrm", "{sm}.e"),
    params:
        folder = join(config['workdir'], "03.chrM", "{sm}"),
    threads:
        int(allocated("threads", "cram__to_chrm", cluster))
    container:
        config['container']['samtools']
    shell:
        """
        cd {params.folder}
        samtools view -@ {threads} --bam -o {output.bam} {input.cram} chrM > {log.out} 2>{log.err}
        samtools index -@ {threads} {output.bam} >> {log.out} 2>>{log.err}
        """