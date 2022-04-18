#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 09:47:09 2022

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
wkdir="/analysis1/yihang_analysis/pipeline"
bam_folder="/analysis1/yihang_analysis/pipeline/data/bam"

bwa='bwa'
threads=96
samtools='samtools'
kaiju='kaiju'
python3='python3'
Rscript='Rscript'
kaiju_addTaxonNames='kaiju-addTaxonNames'
kaiju_nodes='/home/yihang/software/kaijudb_nr_euk/nodes.dmp'
kaiju_fmi='/home/yihang/software/kaijudb_nr_euk/nr_euk/kaiju_db_nr_euk.fmi'
kaiju_names='/home/yihang/software/kaijudb_nr_euk/names.dmp'
bamToFastq='bamToFastq'
seqkit='seqkit'
humann='/home/yihang/anaconda3/envs/pipeline2/bin/humann3'
humann_renorm_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_renorm_table'
humann_regroup_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_regroup_table'
humann_join_tables = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_join_tables'
humann_split_stratified_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_split_stratified_table'
humann_rename_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_rename_table'

import pandas as pd
import os
import re
import sys
import mwgx_pipeline as mp
from shutil import copy

config="mwgx_pipeline/config.tsv"
df_config=pd.read_csv(config,sep='\t')

samples_list, assembly_list, assembly_loc_list, group_list, group_name_list, bam_list=mp.check_config(df_config)

print("config.tsv check pass.")

assembly="/analysis1/yihang_analysis/pipeline/data/assembly/sample_assembly.fa"
group=["hd"]
group1='[0,1,2]'
group1_name='hinny'
group2='[3,4,5]'
group2_name='horse'
index=str(["\'{x}\'".format(x=x) for x in list(df_config['samples'])]).replace(" ",'')

taxa_level=['taxaID','superkingdom','phylum','class','order','family','genus','species']
os.chdir(wkdir)
bam=list(df_config['bam_dir'])
bam_basename=[os.path.basename(x) for x in bam_list]

cpu_count=os.cpu_count()

rule all:
    input:
        expand("data/reads/{sample}.fq.gz",sample = samples_list)
    

# %%
def return_bam(wildcards):
    bam=bam_list[samples_list.index(wildcards.sample)]
    return bam
rule bam2fastq:
    input:
        return_bam
    output:
        "data/reads/{sample}.fq.gz"
    log:
        "log/bam2fastq.{sample}.log"
    benchmark:
        "benchmark/bam2fastq.{sample}.benchmark"
    threads: 2
    run:
        shell("{bamToFastq} -i {input} -fq data/reads/{wildcards.sample}.fq.tmp")
        shell("cat data/reads/{fq}.tmp | {seqkit} rmdup --threads {threads} -s -o data/reads/{wildcards.sample}.fq.gz")
        shell("rm data/reads/{wildcards.sample}.fq.tmp")

rule bam2fastq:
    input:
        bam_list
    output:
        expand("data/reads/{fq}.gz",fq=[bam+".fq" for bam in bam_basename])
    log:
        "log/bam2fastq.log"
    benchmark:
        "benchmark/bam2fastq.benchmark"
    threads: threads
    run:
        for c,b in enumerate(bam_basename):
            fq = b+".fq"
            b_abs=bam[c]
            #print(b_abs)
            shell("{bamToFastq} -i {b_abs} -fq data/reads/{fq}.tmp ")
            shell("cat data/reads/{fq}.tmp | {seqkit} rmdup --threads {threads} -s -o data/reads/{fq}.gz")
            shell("rm data/reads/{fq}.tmp")
mp.func_end("bam2fastq")

# %%
# humann_databases --download chocophlan full /path/to/databases --update-config yes
# humann_databases --download uniref uniref90_diamond /path/to/databases --update-config yes
# humann_databases --download utility_mapping full /path/to/databases --update-config yes
## if you just want to update database location
# humann_config --update database_folders nucleotide /path/to/location
# humann_config --update database_folders protein /path/to/location
# humann_config --update database_folders utility_mapping /path/to/location
rule humann:
    input:
        "data/reads/{sample}.fq.gz"
    output:
        "metabolism/humann3/ori_results/{sample}_pathabundance.tsv",
        "metabolism/humann3/ori_results/{sample}_pathcoverage.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies.tsv"
    log:
        "log/humann.{sample}.log"
    benchmark:
        "benchmark/humann.{sample}.benchmark"
    threads: 4
    run:
        shell("{humann} --input {input}  --output metabolism/humann3/ori_results/ "
        " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
        " --threads {threads} --memory-use maximum")



rule humann_upstream: ## need to combine humann3.sh -> huamnn3_downstream -> humann3_join
    input:
        "metabolism/humann3/ori_results/{sample}_pathabundance.tsv",
        "metabolism/humann3/ori_results/{sample}_pathcoverage.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies.tsv"
    output:
        "metabolism/humann3/results"
    log:
        "log/humann.log"
    benchmark:
        "benchmark/humann.benchmark"
    threads: threads
    run:
        cmd=("find $(pwd)/data/reads -name *.fq.gz | parallel --eta -j "
         " {samples_num} --load 90% --memfree 15G --noswap "
         " '{humann} --input {{}}  --output metabolism/humann3/results/ "
         " --search-mode uniref90 --diamond-options \""--block-size 10 --fast \" "
         " --threads {threads} --memory-use maximum'".\
        format(humann=humann,threads=threads,samples_num=len(samples_list)))
        os.system(cmd)
# %%
rule humann_downstream:
    input:
        "metabolism/humann3/ori_results/{sample}_pathabundance.tsv",
        "metabolism/humann3/ori_results/{sample}_pathcoverage.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies.tsv"
    output:
        "metabolism/humann3/ori_results/{sample}_pathabundance_cpm.tsv",
        "metabolism/humann3/ori_results/{sample}_pathcoverage_cpm.tsv",
        "metabolism/humann3/ori_results/{sample}_pathabundance_relab.tsv",
        "metabolism/humann3/ori_results/{sample}_pathcoverage_relab.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_cpm_level4ec.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_relab_level4ec.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_cpm_ko.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_relab_ko.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_cpm_eggnog.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_relab_eggnog.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_cpm_pfam.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_relab_pfam.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_cpm_rxn.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies_relab_rxn.tsv"
    log:
        "log/humann_downstream.{sample}.log"
    benchmark:
        "benchmark/humann_downstream.{sample}.benchmark"
    run:
        # normalize
        shell("{humann_renorm_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --units relab &>/dev/null")
        shell("{humann_renorm_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --units cpm &>/dev/null")
        shell("{humann_renorm_table} --input metabolism/humann3/ori_results/{wildcards.sample}_pathabundance.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_pathabundance_relab.tsv "
              " --units relab &>/dev/null")
        shell("{humann_renorm_table} --input metabolism/humann3/ori_results/{wildcards.sample}_pathabundance.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_pathabundance_cpm.tsv "
              " --units cpm &>/dev/null")
        shell("{humann_renorm_table} --input metabolism/humann3/ori_results/{wildcards.sample}_pathcoverage.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_pathcoverage_relab.tsv "
              " --units relab &>/dev/null")
        shell("{humann_renorm_table} --input metabolism/humann3/ori_results/{wildcards.sample}_pathcoverage.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_pathcoverage_cpm.tsv "
              " --units cpm &>/dev/null")
        # regroup
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_level4ec.tsv --groups uniref90_level4ec &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab_level4ec.tsv --groups uniref90_level4ec &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_ko.tsv --groups uniref90_ko &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab_ko.tsv --groups uniref90_ko &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_eggnog.tsv --groups uniref90_eggnog &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab_eggnog.tsv --groups uniref90_eggnog &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_pfam.tsv --groups uniref90_pfam &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab_pfam.tsv --groups uniref90_pfam &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_cpm_rxn.tsv --groups uniref90_rxn &>/dev/null")
        shell("{humann_regroup_table} --input metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab.tsv "
              " --output metabolism/humann3/ori_results/{wildcards.sample}_genefamilies_relab_rxn.tsv --groups uniref90_rxn &>/dev/null")

# %%
rule humann_regroup:
    input:
        expand("metabolism/humann3/ori_results/{sample}_pathabundance_cpm.tsv",sample = samples_list)
    output:
        "humann_regroup.done"
    run:
        for c,name in enumerate(group_name_list):
            samples = [samples_list[i] for i in group_list[c]]
            for sample in samples:
                cmd = "mkdir -p metabolism/humann3/ori_results/{group} ; cp metabolism/humann3/ori_results/{sample}*tsv metabolism/humann3/ori_results/{group}"\
                    .format(group = name, sample = sample)
                os.system(cmd)
        shell("touch humann_regroup.done")

# %%
rule humann_join: ## need to combine humann3.sh -> huamnn3_downstream -> humann3_join
    input:
        "metabolism/humann3/ori_results/{sample}_pathabundance.tsv",
        "metabolism/humann3/ori_results/{sample}_pathcoverage.tsv",
        "metabolism/humann3/ori_results/{sample}_genefamilies.tsv"
    output:
        "metabolism/humann3/results"
    log:
        "log/humann.log"
    benchmark:
        "benchmark/humann.benchmark"
    threads: threads
    run:
        cmd=("find $(pwd)/data/reads -name *.fq.gz | parallel --eta -j "
         " {samples_num} --load 90% --memfree 15G --noswap "
         " '{humann} --input {{}}  --output metabolism/humann3/results/ "
         " --search-mode uniref90 --diamond-options \""--block-size 10 --fast \" "
         " --threads {threads} --memory-use maximum'".\
        format(humann=humann,threads=threads,samples_num=len(samples_list)))





        
rule kaiju_annotate:
    input:
        assembly_loc_list
    output:
        expand("kaiju/{assembly}_kaiju.ref",assembly=assembly_loc_list)
    threads: threads
    log:
        "log/kaiju_annotate.log"
    benchmark:
        "benchmark/kaiju_annotate.benchmark"
    shell:
        "{kaiju} -t {kaiju_nodes} -v -f {kaiju_fmi} -z {threads} -i {input}  -o {output}"


rule kaiju_addTaxonNames:
    input:
        expand("kaiju/{group}_kaiju.ref",group=group)
    output:
        expand("kaiju/{group}_kaiju.ref.nm",group=group)
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    threads: threads
    params:
        kaiju=kaiju,
        kaiju_addTaxonNames=kaiju_addTaxonNames,
        kaiju_nodes=kaiju_nodes,
        kaiju_names=kaiju_names
    shell:
        "{kaiju_addTaxonNames} -i {input} -o {output} -t  {kaiju_nodes} -n {kaiju_names} "
        " -r superkingdom,phylum,class,order,family,genus,species"

rule idxstats:
    input:
        bam
    output:
        expand("bam_counts/{bam_basename}.idxstats",bam_basename=bam_basename)
    log:
        "log/{rules.current}.log"
    benchmark:
        "benchmark/{rules.current}.benchmark"
    threads: threads
    run:
        for i,o in zip(input,output):
            shell("samtools idxstats --threads {threads} {i} > {o} ")

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
    





