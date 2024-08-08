rule vcf_compress:
    input:
        vcf = "{path}/{file}.vcf"
    output:
        tbi = "{path}/{file}.vcf.gz.tbi"
    params:
        vcfgz = "{path}/{file}.vcf.gz"
    threads:
        2
    container:
        config['container']['samtools']
    shell:
        "bgzip -c {input.vcf} > {params.vcfgz}\n"
        "tabix -p vcf {params.vcfgz} "