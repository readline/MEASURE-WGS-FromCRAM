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
dic_patient_to_eventsamples, dic_sample_to_runs, dic_run, dic_sample_to_patient, tasks = samplesheet(join(config['pipelinedir'],'sample_sheet.txt'))

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



rule all:
    input:
        expand(
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
            join( config['workdir'], "11.germline_snv_deepvariant", "{sample}", "{sample}.deepvariant.vcf.gz" ),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "05.peddy", "{sample}", "{sample}.peddy.ped"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "12.germline_snv_strelka", "{sample}", "results", "variants", "variants.vcf.gz"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "13.germline_sv_manta", "{sample}", "results", "variants", "candidateSV.vcf.gz"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "14.germline_sv_tiddit", "{sample}", "{sample}.tiddit.vcf"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "15.germline_sv_gridss", "{sample}", "{sample}.gridss.vcf"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Insert", "LINE1.final_comp.vcf"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_LINE1", "DEL.final_comp.vcf"),
            sample=dic_sample_to_runs,
        ),
        expand(
            join(config['workdir'], "21.germline_mei_melt", "{sample}", "ME_Deletion_AluY", "DEL.final_comp.vcf"),
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