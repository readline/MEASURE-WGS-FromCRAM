rule interval_prep:
    input:
        interval = join("config[references][gatkbundle]", "scattered_calling_intervals", "{itv}", "scattered.interval_list"),
    output:
        bed   = join(config['workdir'], "09.supp_data", "scattered_calling_intervals", "{itv}.bed"),
        bedgz = join(config['workdir'], "09.supp_data", "scattered_calling_intervals", "{itv}.bed.gz"),
    params:
        dir = join(config['workdir'], "09.supp_data", "scattered_calling_intervals")
    log:
        out = join(config['pipelinedir'], "logs", "interval_prep", "{itv}.o"),
        err = join(config['pipelinedir'], "logs", "interval_prep", "{itv}.e"),
    threads:
        int(allocated("threads", "interval_prep", cluster))
    container:
        config['container']['gatk']
    shell:
        "cd {params.dir} \n"
        "grep ^chr {input.interval} | cut -f1-3 > {output.bed} \n"
        "bgzip -c {output.bed} > {output.bedgz} \n"
        "tabix -p bed {output.bedgz}"