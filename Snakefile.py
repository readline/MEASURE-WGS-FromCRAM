import os,sys
import yaml
from os.path import join
import pandas as pd
from scripts.load import samplesheet
from scripts.utils import (allocated, ignore)

# Load config file
configfile: 'config.yaml'
workdir: config['workdir']

pipedir = config['pipelinedir']
print(config)

# Load cluster config
with open(join(config['pipelinedir'], 'cluster.yaml')) as infile:
    cluster = yaml.safe_load(infile)
# Load sample sheet
dic_patient_to_eventsamples, dic_sample_to_runs, dic_run, dic_sample_to_patient, dic_sample_to_sid, dic_tumor_to_normal, tasks = \
    samplesheet(join(config['pipelinedir'],'sample_sheet.txt'))

snakedir = os.getcwd()
print('Runtime dir:', snakedir)

# Init interval lists
itvdv = ['%.5d'%(itv) for itv in range(int(config['parameter']['deepvariant_shards'])+1)]

# include rules
include: join(config['pipelinedir'], "rules", "cram.smk")
include: join(config['pipelinedir'], "rules", "qc.smk")
include: join(config['pipelinedir'], "rules", "germline.smk")
include: join(config['pipelinedir'], "rules", "somatic_ss.smk")
include: join(config['pipelinedir'], "rules", "somatic_paired.smk")

# include wrappers
include: join(config['pipelinedir'], "rules", "wrappers", "vcf2tbi.smk")


def configure_all(wildcards):
    inputlist = [expand(
                    join(config['workdir'], "01.cram", "{sample}", "QC", "Summary.ok"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "02.bam", "{sample}", "{sample}.bam"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "03.chrM", "{sample}", "{sample}.chrM.bam"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "04.verifybamid2", "{sample}", "{sample}.selfSM"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "01.cram", "{sample}", "QC", "Mosdepth", "{sample}.mosdepth.summary.txt"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf.gz" ),
                    sample=dic_sample_to_runs,
                ),
                join(config['workdir'], "10.germline_snv_gatk", "Merge.flt.vqsr.vcf.gz"),
                expand(
                    join(config['workdir'], "05.peddy", "{sample}", "{sample}.peddy.ped"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "12.germline_snv_strelka", "{sample}", "results", "variants", "variants.vcf.gz"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "13.germline_sv_manta", "{sample}", "results", "variants", "diploidSV.vcf.gz"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "13.germline_sv_manta", "{sample}", "AnnotSV", "{sample}.manta.germline.tsv"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "14.germline_sv_tiddit", "{sample}", "{sample}.tiddit.vcf"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "14.germline_sv_tiddit", "{sample}", "AnnotSV", "{sample}.tiddit.pass.tsv"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.gridss.vcf"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "15.germline_sv_gridss", "{sample}", "gripss_germline", "AnnotSV", "{sample}.gripss.pass.tsv"),
                    sample=tasks['germline'],
                ),
                expand(
                    join(config['workdir'], "16.germline_sv_virusbreakend", "{sample}", "{sample}.virusbreakend.vcf.gz"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "17.germline_cnv_canvas", "{sample}", "CNV.vcf.gz"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "17.germline_cnv_canvas", "{sample}", "AnnotSV", "Canvas.CNV.pass.tsv"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Insert", "AnnotSV", "{sample}.LINE1.final_comp.tsv"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_LINE1", "AnnotSV", "{sample}.DEL.final_comp.tsv"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_AluY", "AnnotSV", "{sample}.DEL.final_comp.tsv"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "23.germline_msi_msisensorpro", "{sample}", "{sample}.msisensorpro"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "25.germline_hla-la", "{sample}", "hla", "summaryStatistics.txt"),
                    sample=dic_sample_to_runs,
                ),
                expand(
                    join(config['workdir'], "41.somatic_ss_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
                    sample=tasks['somatic_ss'],
                ),
                expand(
                    join(config['workdir'], "42.somatic_ss_snvindel_octopus", "{sample}", "{sample}.octopus.pass.vcf.gz"),
                    sample=tasks['somatic_ss'],
                ),
                expand(
                    join(config['workdir'], "43.somatic_ss_snvindel_vardict", "{sample}", "{sample}.vardict.pass.vcf.gz"),
                    sample=tasks['somatic_ss']+[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "44.somatic_ss_sv__gridss", "{sample}",  "AnnotSV", "{sample}.gripss.pass.tsv"),
                    sample=tasks['somatic_ss'],
                ),
                expand(
                    join(config['workdir'], "45.somatic_ss_sv__delly", "JointCalls",  "AnnotSV", "Merge.delly.somatic.tsv"),
                ),
                expand(
                    join(config['workdir'], "46.somatic_ss_cnv__cnvkit", "{sample}", "cnvkit.ok"),
                    sample=tasks['somatic_ss']+[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "47.somatic_ss_cnv__purple", "{sample}", "{sample}.purple.cnv.somatic.tsv"),
                    sample=tasks['somatic_ss']+[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "48.somatic_ss_cnv__canvas", "{sample}", "AnnotSV", "Canvas.CNV.pass.tsv"),
                    sample=tasks['somatic_ss']+[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "31.somatic_tn_snvindel_mutect2", "{sample}", "{sample}.mutect2.pass.vcf.gz"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "32.somatic_snv_strelka", "{sample}", "results", "variants", "somatic.snvs.vcf.gz"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "33.somatic_snv_muse", "{sample}", "{sample}.MuSE.pass.vcf.gz"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "34.somatic_sv_manta", "{sample}", "AnnotSV", "{sample}.manta.somatic.tsv"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "{sample}.gripss.filtered.vcf.gz"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "35.somatic_tn_sv__gridss", "{sample}", "AnnotSV", "{sample}.gripss.somatic.tsv"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                ),
                expand(
                    join(config['workdir'], "38.somatic_tn_cnv__purple", "{sample}", "{sample}.purple.cnv.somatic.tsv"),
                    sample=[i[0] for i in tasks['somatic_paired']],
                )
    ]
    if len(tasks['somatic_paired']) > 0:
        inputlist.append(join(config['workdir'], "36.somatic_tn_sv__delly", "JointCalls",  "AnnotSV", "Merge.delly.somatic.tsv"))
    return inputlist

rule all:
    input:
        unpack(configure_all)