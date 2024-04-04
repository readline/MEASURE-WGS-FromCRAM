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
samples, patients, tasks = samplesheet(join(config['pipelinedir'],'sample_sheet.txt'))
# include rules
include: join(config['pipelinedir'], "rules", "cram.smk")
include: join(config['pipelinedir'], "rules", "qc.smk")
include: join(config['pipelinedir'], "rules", "germline.smk")
include: join(config['pipelinedir'], "rules", "somatic_ss.smk")
include: join(config['pipelinedir'], "rules", "somatic_paired.smk")

snakedir = os.getcwd()
print('Runtime dir:', snakedir)

rule all:
    input:
        expand(
            join(config['workdir'], "01.cram", "{sm}", "{sm}.cram"),
            sm=samples,
        ),
        expand(
            join(config['workdir'], "02.bam", "{sm}", "{sm}.bam"),
            sm=samples,
        ),
        expand(
            join(config['workdir'], "03.chrM", "{sm}", "{sm}.chrM.bam"),
            sm=samples,
        ),