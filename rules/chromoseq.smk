from scripts.utils import allocated,ignore
# Filter Somatic SV calls
rule manta_somatic_filter:
    """
    Data processing step to call somatic structural variants using manta. 
    Manta is optimized for optimized for analysis of germline variation 
    in small sets of individuals and somatic variation in tumor/normal 
    sample pairs. Manta's Github Repo: https://github.com/Illumina/manta
    @Input:
        Single-sample VCF file with called somatic structural variants (scatter)
    @Output:
        Filtered single-sample VCF file with called somatic structural variants
    """
    input: 
        vcf  = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.vcf.gz"),
    output:
        samplevcf  = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.sample.vcf.gz"),
        flagged    = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.flagged.vcf.gz"),
        filtered   = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.filtered.vcf"),
    params: 
        rname       = "manta_somatic_filter",
        sample      = "{name}",
        outdir      = join(workpath, "MANTA", "somatic", "{name}"),
        filtermanta = join(workpath, "workflow", "scripts", "FilterManta.pl"),
        genome      = config['references']['GENOME'],
        filter_ref  = config['references']['MANTA_FILTER_CHROMOSEQ_TRANSLOCATION'],
        memory      = allocated("mem", "manta_somatic_filter", cluster).rstrip('G'),
    threads: int(allocated("threads", "manta_somatic_filter", cluster))
    container: config['images']['genome-seek_sv']
    envmodules: 
        config['tools']['bcftools'],
        config['tools']['samtools'],
        config['tools']['bedtools'],
        config['tools']['svtools'],
        config['tools']['minimap2'],
        config['tools']['perl'],
    shell: """
    # Flag and filter SVs based on the
    # following: read support, contig 
    # re-mapping, and allele fraction. 
    # Also removes non-PASS SVs from
    # somatic SV callset.
    bcftools view \\
        -s {params.sample} \\
        -O z \\
        -o {output.samplevcf} \\
        {input.vcf}
    # Script to filter manta results 
    # from chromoseq pipeline,
    # https://github.com/genome/docker-basespace_chromoseq
    # Needs paths to its depedencies:
   bedtools_path="$(type -P bedtools)"
    minimap2_path="$(type -P minimap2)"
    svtools_path="$(type -P svtools)"
    echo "Paths of all dependencies... ${{bedtools_path}}:${{minimap2_path}}:${{svtools_path}}"
    perl {params.filtermanta} \\
        -r {params.genome} \\
        -k {params.filter_ref} \\
        -b "${{bedtools_path}}" \\
        -p "${{minimap2_path}}" \\
        -s "${{svtools_path}}" \\
        -t {params.outdir} \\
        {output.samplevcf} \\
        {output.flagged}
    # Remove non-PASS-ing SVs
    bcftools view \\
        -O v \\
        -i 'FILTER=="PASS"' \\
        {output.flagged} \\
    > {output.filtered}
    """

#Start Chromoseq Rules

rule duphold:
        input:
            svs = join(workpath, "MANTA", "somatic", "{name}", "results", "variants", "somaticSV.vcf.gz"),
            bam = join(workpath, "BAM", "{name}.recal.bam"),
        output:
            sv = join(workpath, "chromoseq", "{name}", "{name}.tumorSV.vcf"),
            fixed = join(workpath, "chromoseq", "{name}", "fixed.vcf"),
            sv_noextra = join(workpath, "chromoseq", "{name}", "{name}.noextra.vcf.gz")
        params:
            genome = config['references']['GENOME'],
            sample = "{name}",
            rname = "duphold"
        shell:  """
                mkdir -p chromoseq
                mkdir -p chromoseq/{params.sample}
                zcat {input.svs} | sed 's/DUP:TANDEM/DUP/g' > {output.fixed} && /data/OpenOmics/references/genome-seek/chromoseq/duphold -f {params.genome} -v {output.fixed} -b {input.bam} -t 8 -o {output.sv} && grep -v -E >
                """
