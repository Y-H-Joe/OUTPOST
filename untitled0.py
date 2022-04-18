#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 14:42:52 2022

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####

=================================== input =====================================

=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================

####=======================================================================####
"""
import pandas as pd
import os
import re
import sys
import GEMINI as ge

# %% settings

wkdir = "/analysis1/yihang_analysis/pipeline"
os.chdir(wkdir)

bwa = 'bwa'
threads = 96
samtools = 'samtools'
kaiju = 'kaiju'
python3 = 'python3'
Rscript = 'Rscript'
kaiju_addTaxonNames = 'kaiju-addTaxonNames'
kaiju_nodes = '/home/yihang/software/kaijudb_nr_euk/nodes.dmp'
kaiju_fmi = '/home/yihang/software/kaijudb_nr_euk/nr_euk/kaiju_db_nr_euk.fmi'
kaiju_names = '/home/yihang/software/kaijudb_nr_euk/names.dmp'
taxa_level = ['taxaID','superkingdom','phylum','class','order','family','genus','species']
bamToFastq  ='bamToFastq'
samtools = 'samtools'
seqkit = 'seqkit'
humann = '/home/yihang/anaconda3/envs/pipeline2/bin/humann3'
humann_renorm_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_renorm_table'
humann_regroup_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_regroup_table'
humann_join_tables = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_join_tables'
humann_split_stratified_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_split_stratified_table'
humann_rename_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_rename_table'


config="GEMINI/config.tsv"
df_config=pd.read_csv(config,sep='\t')

sample_list, fq_list, bam_list, assembly, assembly_dir, comparison_dict = ge.check_config(df_config)
bam_basename=[os.path.basename(x) for x in bam_list]

rule all:
    input:
        "humann.done",
        f"{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm",
        "idxstats.done",
        expand("{assembly}/taxa_analysis/temp/{bam_basename}.idxstats",assembly = assembly, bam_basename = bam_basename),
        expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.named.tsv",
               group = comparison_dict.keys(), assembly = assembly),
        expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.named.tsv",
               group = comparison_dict.keys(), assembly = assembly),
        expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.named.tsv",
               group = comparison_dict.keys(), assembly = assembly),
        expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv",
               group = comparison_dict.keys(), assembly = assembly),
        expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.named.tsv",
               group = comparison_dict.keys(), assembly = assembly)
# %% bam
"""
def return_bam(wildcards):
     bam=bam_list[sample_list.index(wildcards.sample)]
     return bam

rule bam2fastq:
    input:
        return_bam
    output:
        "data/reads/{sample}.fq.gz",
        "bam2fastq.done"
    log:
        "log/bam2fastq.{sample}.log"
    benchmark:
        "benchmark/bam2fastq.{sample}.benchmark"
    threads: 2
    run:
        shell("{bamToFastq} -i {input} -fq data/reads/{wildcards.sample}.fq.tmp")
        shell("cat data/reads/{fq}.tmp | {seqkit} rmdup --threads {threads} -s -o data/reads/{wildcards.sample}.fq.gz")
        shell("rm data/reads/{wildcards.sample}.fq.tmp")
        shell("touch bam2fastq.done")
"""

# %% humann
def return_fq(wildcards):
    return fq_list[sample_list.index(wildcards.sample)]

rule humann_init:
    input:
        return_fq
    output:
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv"
    log:
        "{assembly}/log/humann_init.{sample}.log"
    benchmark:
        "{assembly}/benchmark/humann_init.{sample}.benchmark"
    threads: threads
    run:
        shell("{humann} --input {input}  --output {assembly}/metabolism_analysis/humann3/ori_results/ "
        " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
        " --threads {threads} --memory-use maximum")


rule humann_annotate:
    input:
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv"
    output:
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance_cpm.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage_cpm.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance_relab.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage_relab.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_level4ec.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_level4ec.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_ko.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_ko.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_eggnog.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_eggnog.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_pfam.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_pfam.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_rxn.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_rxn.tsv"
    log:
        "{assembly}/log/humann_annotate.{sample}.log"
    benchmark:
        "{assembly}/benchmark/humann_annotate.{sample}.benchmark"
    run:
        # normalize
        shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --units relab")
        shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --units cpm")
        shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathabundance.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathabundance_relab.tsv "
              " --units relab")
        shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathabundance.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathabundance_cpm.tsv "
              " --units cpm")
        shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathcoverage.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathcoverage_relab.tsv "
              " --units relab")
        shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathcoverage.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_pathcoverage_cpm.tsv "
              " --units cpm")
        # regroup
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_level4ec.tsv --groups uniref90_level4ec")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab_level4ec.tsv --groups uniref90_level4ec")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_ko.tsv --groups uniref90_ko")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab_ko.tsv --groups uniref90_ko")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_eggnog.tsv --groups uniref90_eggnog")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab_eggnog.tsv --groups uniref90_eggnog")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_pfam.tsv --groups uniref90_pfam")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab_pfam.tsv --groups uniref90_pfam")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_rxn.tsv --groups uniref90_rxn")
        shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output {assembly}/metabolism_analysis/humann3/ori_results/{wildcards.sample}_genefamilies_relab_rxn.tsv --groups uniref90_rxn")

rule humann_group:
    input:
# =============================================================================
#         expand("{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_level4ec.tsv",sample = sample_list, assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_ko.tsv",sample = sample_list, assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_eggnog.tsv",sample = sample_list, assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_rxn.tsv",sample = sample_list, assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_pfam.tsv",sample = sample_list, assembly = assembly)
# =============================================================================
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_level4ec.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_ko.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_eggnog.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_rxn.tsv",
        "{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_pfam.tsv"
    output:
        directory("{assembly}/metabolism_analysis/humann3/{group}")
    run:
        for group,samples in comparison_dict.items():
            os.makedirs(f"{assembly}/metabolism_analysis/humann3/{group}", exist_ok=True)
            for sample in samples:
                os.system(f"cp {assembly}/metabolism_analysis/humann3/ori_results/{sample}*tsv  {assembly}/metabolism_analysis/humann3/{group} ")

rule humann_output:
    input:
        ("{assembly}/metabolism_analysis/humann3/{group}", group = comparison_dict.keys(), assembly = assembly)
    output:
        "humann.done",
# =============================================================================
#         expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.named.tsv",
#                group = comparison_dict.keys(), assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.named.tsv",
#                group = comparison_dict.keys(), assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.named.tsv",
#                group = comparison_dict.keys(), assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv",
#                group = comparison_dict.keys(), assembly = assembly),
#         expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.named.tsv",
#                group = comparison_dict.keys(), assembly = assembly)
# =============================================================================
        "{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.named.tsv",
        "{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.named.tsv",
        "{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.named.tsv",
        "{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv",
        "{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.named.tsv"
    run:
        for group in comparison_dict.keys():
            # join tables for each group
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  "--file_name _genefamilies_relab.tsv ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  "--file_name _genefamilies_cpm.tsv ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathabundance_relab.tsv "
                  "--file_name _pathabundance_relab.tsv ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathabundance_cpm.tsv "
                  "--file_name _pathabundance_cpm.tsv ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathcoverage_relab.tsv "
                  "--file_name _pathcoverage_relab.tsv ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathcoverage_cpm.tsv "
                  "--file_name _pathcoverage_cpm.tsv ")
            # regroup tables for each group
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec.tsv --groups uniref90_level4ec")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko.tsv --groups uniref90_ko")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog.tsv --groups uniref90_eggnog")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam.tsv --groups uniref90_pfam")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn.tsv --groups uniref90_rxn")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec.tsv --groups uniref90_level4ec")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko.tsv --groups uniref90_ko")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog.tsv --groups uniref90_eggnog")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam.tsv --groups uniref90_pfam")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn.tsv --groups uniref90_rxn")
            # stratify tables for each group
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn.tsv "
                  "-o {assembly}/metabolism_analysis/humann3/output/")
            # rename tables for each group
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.tsv "
                  " --names kegg-orthology -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.tsv "
                  " --names eggnog -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.tsv "
                  " --names ec -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.tsv "
                  " --names pfam -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.tsv "
                  " --names metacyc-rxn -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko_unstratified.tsv "
                  " --names kegg-orthology -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog_unstratified.tsv "
                  " --names eggnog -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec_unstratified.tsv "
                  " --names ec -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam_unstratified.tsv "
                  " --names pfam -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam_unstratified.named.tsv")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.tsv "
                  " --names metacyc-rxn -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.named.tsv")
            shell("touch humann.done")

# %% taxonomy

rule kaiju_annotate:
    input:
        assembly_dir
    output:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref"
    threads: threads
    log:
        "{assembly}/log/kaiju_annotate.log"
    benchmark:
        "{assembly}/benchmark/kaiju_annotate.benchmark"
    shell:
        "{kaiju} -t {kaiju_nodes} -v -f {kaiju_fmi} -z {threads} -i {input}  -o {output}"

rule kaiju_addTaxonNames:
    input:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref"
    output:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm"
    log:
        "{assembly}/log/kaiju_addTaxonNames.log"
    benchmark:
        "{assembly}/benchmark/kaiju_addTaxonNames.benchmark"
    threads: threads
    shell:
        "{kaiju_addTaxonNames} -i {input} -o {output} -t  {kaiju_nodes} -n {kaiju_names} "
        " -v -r superkingdom,phylum,class,order,family,genus,species"

def return_bam(wildcards):
     bam = bam_list[sample_list.index(wildcards.sample)]
     return bam

rule idxstats:
    input:
        return_bam
    output:
        "{assembly}/taxa_analysis/temp/{bam_basename}.idxstats",
        "idxstats.done"
    log:
        "{assembly}/log/idxstats.{bam_basename}.log"
    benchmark:
        "{assembly}/benchmark/idxstats.{bam_basename}.benchmark"
    threads: threads
    run:
        for basename in bam_basename:
            os.system(f"{samtools} idxstats --threads {threads} {input} > {assembly}/taxa_analysis/temp/{basename}.idxstats ")
        os.system("touch idxstats.done")




rule paste_counts_table_together:
    input:
        expand("bam_counts/{bam_basename}.idxstats",bam_basename=bam_basename)
    output:
        expand("bam_counts/{group}.counts.pure",group=group)
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    run:
        with open("code/namefile.paste_counts_table_together.txt",'w') as f:
            for i in df_config.index:
                f.write("bam_counts/{bam_basename}.idxstats".format(bam_basename=bam_basename[i]))
                f.write('\n')
        for g in group:
            shell("{python3} code/paste_counts_table_together.py 3 code/namefile.paste_counts_table_together.txt "
             " bam_counts/{g}.counts.pure".format(python3=python3,g=g))
            shell("sed -i '$d' bam_counts/{g}.counts.pure".format(g=g))

rule format_kaiju_output_to_tab_seperated:
    input:
        expand("kaiju/{group}_kaiju.ref.nm",group=group)
    output:
        expand("kaiju/{group}_kaiju.ref.nm.tsv",group=group)
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    run:
        for g in group:
            shell("{python3} code/format_kaiju_output_to_tab_seperated.py "
            " kaiju/{g}_kaiju.ref.nm ".format(python3=python3,g=g))

rule merge_counts_pure_and_kaiju:
    input:
        counts_dp=expand("bam_counts/{group}.counts.pure",group=group),
        kaiju_dp=expand("kaiju/{group}_kaiju.ref.nm.tsv",group=group)
    output:
        expand("bam_counts/{group}.counts.taxa.tsv",group=group)
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    run:
        for g in group:
            shell("{Rscript} code/merge_counts_and_kaiju.R bam_counts/{g}.counts.pure "
            " kaiju/{g}_kaiju.ref.nm.tsv bam_counts/{g}.counts.taxa.tsv".format(Rscript=Rscript,g=g))

rule prepare_contig_table_from_counts_table: ## need to modify the header of contig_table
    input:
        expand("bam_counts/{group}.counts.pure",group=group)
    output:
        expand("MAG/{group}.contig_table.tsv",group=group)
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    run:
        for g in group:
            shell("{python3} code/prepare_contig_table_from_counts_table.py "
            " bam_counts/{g}.counts.taxa.tsv {group1} {group2} "
            " MAG/{g}.contig_table.tsv".format(python3=python3,g=g,group1=group1,group2=group2))

rule counts_table2rel_abun:
    input:
        expand("bam_counts/{group}.counts.taxa.tsv",group=group)
    output:
        ["taxa_abun/rel_abun/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)],
        ["taxa_abun/rel_abun/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)]
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    run:
        for g in group:
            prefix="taxa_abun/rel_abun/{group1_name}_vs_{group2_name}.{g}"\
            .format(group1_name=group1_name,group2_name=group2_name,g=g)

            shell("{python3} code/counts_table2rel_abun.py bam_counts/{g}.counts.taxa.tsv "
            " {prefix} {index} ".format(python3=python3,prefix=prefix,index=index,g=g))

rule rel_abun_utest:
    input:
        ["taxa_abun/rel_abun/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)]
    output:
        ["taxa_abun/utest/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.average_change.equal.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)],

        ["taxa_abun/utest/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.average_change.unequal.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)],

        ["taxa_abun/utest/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.u-test.two_sided.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)],

        ["taxa_abun/utest/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.rel_abun.equal.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)],

        ["taxa_abun/utest/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.rel_abun.unequal.csv"\
        .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i) for g in group for i in range(1,9)]
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    run:
        for g in group: # need to modify to fit multi-group cases
            for i in range(1,9):
                infile="taxa_abun/rel_abun/{group1_name}_vs_{group2_name}.{g}.rel_abun.{i}.rmU.csv"\
                .format(group1_name=group1_name,group2_name=group2_name,g=g,i=i)

                shell("{python3} code/rel_abun_utest.py {infile} {group1} {group1_name} {group2} {group2_name} "\
                .format(python3=python3,infile=infile,group1=group1,group1_name=group1_name,group2=group2,group2_name=group2_name))

rule extract_top_taxa_from_rel_abun_table:
    input:
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv",
        "taxa_abun/utest/sample12_rel_abun.1.rmU.euk.csv_relative_abun_unequal_horse_vs_donkey.csv"
    output:
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv.top30.csv"
        "taxa_abun/utest/sample12_rel_abun.1.rmU.euk.csv_relative_abun_unequal_horse_vs_donkey.csv.top30.csv"
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    scripts:
        "code/extract_top_taxa_by_rel_abun_from_rel_abun_table.py"




























rule scale_rel_abun_table:
    input:
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv.top30.csv"
    output:
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv.top30.csv.log2.fillmin.csv",
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv.top30.csv.log10.fillmin.csv",
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv.top30.csv.norm.fillmin.csv"
    scripts:
        "code/scale_rel_abun_table.py"

rule heatmap:
    input:
        "taxa_abun/rel_abun/sample12_rel_abun.1.rmU.euk.csv.top30.csv.log10.fillmin.csv",
    output:
        "taxa_abun/figs/sample12_rel_abun.1.rmU.euk.csv.top30.csv.log10.fillmin.csv.heatmap.pdf",
    scripts:
        "code/heatmap.py"

rule config2alphabeta_group:
    input:
        "config.tsv",
        "../taxa_abun/rel_abun/sample12_rel_abun.8.rmU.euk.csv"
    output:
        "alpha_beta_diversity_group.tsv",
        "../taxa_abun/alpha_beta_diversity/sample12_rel_abun.8.rmU.euk.csv.alpha_beta.csv"
    scripts:
        "config2alphabeta_group_20211221.py"

rule alpha_beta_diversity:
    input:
        "../taxa_abun/alpha_beta_diversity/sample12_rel_abun.8.rmU.euk.csv.alpha_beta.csv",
        "alpha_beta_diversity_group.tsv"
    output:
        "../taxa_abun/alpha_beta_diversity/..."
    scripts:
        "alpha_beta_diversity.R"

rule rel_abun2lefse:
    input:
        "../taxa_abun/rel_abun/sample12_rel_abun.8.rmU.euk.csv"
    output:
        "../lefse/sample12_rel_abun.8.rmU.euk.csv.lefse.tsv"
    scripts:
        "rel_abun2lefse.py"

rule run_lefse:
    input:
        "../lefse/sample12_rel_abun.8.rmU.euk.csv.lefse.tsv"
    output:
        "../lefse/..."
    scripts:
        "run_lefse.sh"



rule humann2rel_abun:
    input:
        "../metabolism/humann3/horsedonkey_genefamilies_uniref90names_cpm_{}_unstratified.named.tsv"
    output:
        "../metabolism/humann3/horsedonkey_genefamilies_uniref90names_cpm_{}_unstratified.named.tsv.rel_abun.csv"
    scripts:
        "humann2rel_abun.py"

rule alluvial_plot:
    input:
    output:
    scripts:

rule abricate:
    input:
    output:
    scripts:

rule boxplot:
    input:
    output:
    scripts:

rule microbeannotator:
    input:
    output:
    scripts:

rule checkM:
    input:
    output:
    scripts:

rule phylogenetic_tree:
    input:
    output:
    scripts:






