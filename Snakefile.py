import pandas as pd
import os
import re
import sys
import GEMINI as ge
import time
import itertools
from multiprocessing import Pool
from multiprocessing import set_start_method

# %% settings

wkdir = "/analysis1/yihang_analysis/pipeline"
os.chdir(wkdir)

bwa = 'bwa'
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

paired = False
top = 30

config="GEMINI/config.tsv"
df_config=pd.read_csv(config,sep='\t')

sample_list, fq_list, bam_list, assembly, assembly_dir, comparison_dict = ge.check_config(df_config)
bam_basename=[os.path.basename(x) for x in bam_list]
group_pair_list = list(itertools.combinations(comparison_dict.keys(),2))

rule all:
    input:
        f"{assembly}/log/humann.done",
        f"{assembly}/log/prepare_contig_table_from_counts_table.done",
        f"{assembly}/log/scale_rel_abun_table.done",
        f"{assembly}/log/alpha_beta_diversity.done"

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
rule humann_output:
    input:
        f"{assembly}/log/humann_group.done"
    output:
        f"{assembly}/log/humann.done",
        expand("{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.named.tsv",
               assembly = assembly, group = comparison_dict.keys())
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
            shell("touch {assembly}/log/humann.done")

rule humann_group:
    input:
        f"{assembly}/log/humann_annotate.done",
    output:
        f"{assembly}/log/humann_group.done"
    run:
        for group,samples in comparison_dict.items():
            os.makedirs(f"{assembly}/metabolism_analysis/humann3/{group}", exist_ok=True)
            for sample in samples:
                os.system(f"cp {assembly}/metabolism_analysis/humann3/ori_results/{sample}*tsv  {assembly}/metabolism_analysis/humann3/{group} ")
        shell("touch {assembly}/log/humann_group.done")

rule humann_annotate:
    input:
        f"{assembly}/log/humann_init.done"
    output:
        f"{assembly}/log/humann_annotate.done"
    log:
        f"{assembly}/log/humann_annotate.log"
    benchmark:
        f"{assembly}/benchmark/humann_annotate.benchmark"
    run:
        for sample in sample_list:
            # normalize
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --units relab")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --units cpm")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance_relab.tsv "
                  " --units relab")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance_cpm.tsv "
                  " --units cpm")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage_relab.tsv "
                  " --units relab")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage_cpm.tsv "
                  " --units cpm")
            # regroup
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_level4ec.tsv --groups uniref90_level4ec")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_level4ec.tsv --groups uniref90_level4ec")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_ko.tsv --groups uniref90_ko")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_ko.tsv --groups uniref90_ko")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_eggnog.tsv --groups uniref90_eggnog")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_eggnog.tsv --groups uniref90_eggnog")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_pfam.tsv --groups uniref90_pfam")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_pfam.tsv --groups uniref90_pfam")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_rxn.tsv --groups uniref90_rxn")
            shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_rxn.tsv --groups uniref90_rxn")
            shell("touch {assembly}/log/humann_annotate.done")

def return_fq(wildcards):
    return fq_list[sample_list.index(wildcards.sample)]

rule humann_init:
    input:
        fq_list
    output:
        "{assembly}/log/humann_init.done"
    log:
        "{assembly}/log/humann_init.log"
    benchmark:
        "{assembly}/benchmark/humann_init.benchmark"
    threads: 4
    run:
        os.makedirs(f"{assembly}/metabolism_analysis/humann3/ori_results/", exist_ok=True)
        for fq in fq_list:
            shell("{humann} --input {fq}  --output {assembly}/metabolism_analysis/humann3/ori_results/ "
            " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
            " --threads {threads} --memory-use maximum")
        shell("touch {assembly}/log/humann_init.done")

# %% alpha_beta_diversity
rule alpha_beta_diversity:
    input:
        "{assembly}/log/counts_table2rel_abun.done"
    output:
        "{assembly}/log/alpha_beta_diversity.done"
    log:
        "{assembly}/log/alpha_beta_diversity.log"
    benchmark:
        "{assembly}/benchmark/alpha_beta_diversity.benchmark"
    threads: 4
    run:
        for group1,group2 in group_pair_list:
            wkdir = f"{assembly}/alpha_beta_analysis/alpha_beta_{group1}_vs_{group2}"
            os.makedirs(wkdir, exist_ok=True)
            group1_new_name = [f"{sample}_{group1}" for sample in comparison_dict[group1]]
            group2_new_name = [f"{sample}_{group2}" for sample in comparison_dict[group2]]

            group1_name = comparison_dict[group1]
            group2_name = comparison_dict[group2]

            species_dp = f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.species.rmU.csv"
            genus_dp = f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.genus.rmU.csv"
            species_df = pd.read_csv(species_dp,index_col=0)
            genus_df = pd.read_csv(genus_dp,index_col=0)

            species_new_df = species_df.reindex(group1_name + group2_name)
            genus_new_df = genus_df.reindex(group1_name + group2_name)

            species_new_dp = wkdir + f"/{assembly}.taxa_counts.rel_abun.species.rmU.{group1}_vs_{group2}.csv"
            genus_new_dp = wkdir + f"/{assembly}.taxa_counts.rel_abun.genus.rmU.{group1}_vs_{group2}.csv"
            species_new_df.to_csv(species_new_dp,index=True)
            genus_new_df.to_csv(genus_new_dp,index=True)

            index_file = wkdir + "/index.txt"
            with open(index_file,'w') as w:
                for new_name in group1_new_name + group2_new_name:
                    w.write(new_name + "\n")

            group_file = wkdir + "/group.csv"
            with open(group_file,'w') as w:
                w.write("sample,Group,Others\n")
                for new_name in group1_new_name:
                    w.write(f"{new_name},{group1},A\n")
                for new_name in group2_new_name:
                    w.write(f"{new_name},{group2},A\n")

            shell("{Rscript} GEMINI/alpha_beta_diversity.R {genus_new_dp} {wkdir} {index_file} {group_file} {group1},{group2} > {log} 2>&1")
        shell("touch {assembly}/log/alpha_beta_diversity.done")





# %% taxonomy
rule scale_rel_abun_table:
    input:
        "{assembly}/log/extract_top_taxa.done"
    output:
        "{assembly}/log/scale_rel_abun_table.done"
    log:
        "{assembly}/log/scale_rel_abun_table.log"
    benchmark:
        "{assembly}/benchmark/scale_rel_abun_table.benchmark"
    run:
        dp_list = [f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.top{top}.csv"
                       for level in taxa_level]
        output_list = [f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.top{top}.fillmin.log10.csv"
                       for level in taxa_level]

        for group1,group2 in group_pair_list:
            dp_list_equal = [f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.csv"
                                      for level in taxa_level]
            output_list_equal = [f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.fillmin.log10.csv"
                                          for level in taxa_level]

            dp_list_unequal = [f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.csv"
                                        for level in taxa_level]
            output_list_unequal = [f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.fillmin.log10.csv"
                                            for level in taxa_level]

            dp_list += dp_list_equal
            dp_list += dp_list_unequal
            output_list += output_list_equal
            output_list += output_list_unequal
        dp_list = ','.join(dp_list)
        output_list = ','.join(output_list)
        shell("{python3} GEMINI/scale_rel_abun_table.py {dp_list} {output_list}")
        shell("touch {assembly}/log/scale_rel_abun_table.done")

rule extract_top_taxa:
    input:
        "{assembly}/log/rel_abun_utest.done"
    output:
        "{assembly}/log/extract_top_taxa.done"
    log:
        "{assembly}/log/extract_top_taxa.log"
    benchmark:
        "{assembly}/benchmark/extract_top_taxa.benchmark"
    threads: 4
    run:
        # top abundant
        dp_list = ','.join([f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.csv"
                            for level in taxa_level])
        output_name_list = ','.join([f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.top{top}.csv"
                            for level in taxa_level])
        error_log = os.getcwd() + f"/{assembly}/log/extract_top_taxa_top_abundant.error"
        shell("{python3} GEMINI/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
              " {top} {error_log}")
        result = ge.wait_unti_file_exists(output_name_list.split(','),error_log)
        assert result is True, "GEMINI: extract_top_taxa_top_abundant has errors. exit."

        # top unequal abundant
        output_list = []
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}", exist_ok=True)
            dp_list = ','.join([f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.csv"
                                for level in taxa_level])
            output_name_list = ','.join([
                    f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.csv"
                    for level in taxa_level])
            error_log = os.getcwd() + f"/{assembly}/log/extract_top_taxa_top_unequal_abundant.error"
            output_list.append(f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_species.rel_abun.unequal.top{top}.csv")
            shell("nohup {python3} GEMINI/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
              " {top} {error_log} > {log}_{group1}_vs_{group2}_unequal 2>&1 &")
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: extract_top_taxa_top_unequal_abundant has errors. exit."

       # top equal abundant
        output_list = []
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}", exist_ok=True)
            dp_list = ','.join([f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.csv"
                                for level in taxa_level])
            output_name_list = ','.join([
                    f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.csv"
                    for level in taxa_level])
            error_log = os.getcwd() + f"/{assembly}/log/extract_top_taxa_top_equal_abundant.error"
            output_list.append(f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_species.rel_abun.equal.top{top}.csv")
            shell("nohup {python3} GEMINI/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
              " {top} {error_log} > {log}_{group1}_vs_{group2}_equal 2>&1 &")
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: extract_top_taxa_top_equal_abundant has errors. exit."

        shell("touch {assembly}/log/extract_top_taxa.done")


rule rel_abun_utest:
    input:
        "{assembly}/log/counts_table2rel_abun.done"
    output:
        "{assembly}/log/rel_abun_utest.done"
    log:
        "{assembly}/log/rel_abun_utest.log"
    benchmark:
        "{assembly}/benchmark/rel_abun_utest.benchmark"
    threads: 4
    run:
        output_list = []
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}", exist_ok=True)

            dp_list = ','.join([f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.csv" for level in taxa_level])
            group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
            group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
            prefix_list = ','.join([f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}"
                                    for level in taxa_level])
            output = f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_species.ave_change.unequal.csv"
            output_list.append(output)
            error_log = os.getcwd() + f"/{assembly}/log/rel_abun_utest.error"
            try:
                os.remove(error_log) # clean former residual error log
            except:
                pass
            shell("nohup {python3} GEMINI/rel_abun_utest.py {dp_list} {group1_index} {group1} "
                  " {group2_index} {group2} {prefix_list} {paired} {error_log} > {log}_{group1}_vs_{group2} 2>&1 &")
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: rel_abun_utest has errors. exit."
        shell("touch {assembly}/log/rel_abun_utest.done")

rule counts_table2rel_abun:
    input:
        "{assembly}/assembly_analysis/{assembly}.taxa_counts.tsv"
    output:
        "{assembly}/log/counts_table2rel_abun.done"
    log:
        "{assembly}/log/counts_table2rel_abun.log"
    benchmark:
        "{assembly}/benchmark/counts_table2rel_abun.benchmark"
    run:
        prefix = f"{assembly}/taxa_analysis/{assembly}.taxa_counts"
        samples = ','.join(sample_list)
        shell("{python3} GEMINI/counts_table2rel_abun.py {input} {prefix} {samples}")
        shell("touch {assembly}/log/counts_table2rel_abun.done")

rule prepare_contig_table_from_counts_table:
    input:
        "{assembly}/assembly_analysis/{assembly}.taxa_counts.tsv"
    output:
        "{assembly}/log/prepare_contig_table_from_counts_table.done"
    log:
        "{assembly}/log/prepare_contig_table_from_counts_table.log"
    benchmark:
        "{assembly}/benchmark/prepare_contig_table_from_counts_table.benchmark"
    threads: 4
    run:
        output_list = []
        for group1,group2 in group_pair_list:
            group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
            group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
            samples = ','.join(sample_list)
            output = f"{assembly}/assembly_analysis/{assembly}.{group1}_vs_{group2}.contig_table.tsv"
            output_list.append(output)
            error_log = os.getcwd() + f"/{assembly}/log/prepare_contig_table_from_counts_table.error"
            try:
                os.remove(error_log) # clean former residual error log
            except:
                pass
            shell("nohup {python3} GEMINI/prepare_contig_table_from_counts_table.py {input} "
                  " {group1_index} {group2_index} {output} {samples} {error_log} > {log} 2>&1 &")
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: prepare_contig_table_from_counts_table has errors. exit."
        shell("touch {assembly}/log/prepare_contig_table_from_counts_table.done")

rule merge_counts_pure_and_kaiju:
    input:
        counts_dp = "{assembly}/taxa_analysis/temp/{assembly}.counts.pure",
        kaiju_dp = "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm.tsv"
    output:
        "{assembly}/assembly_analysis/{assembly}.taxa_counts.tsv"
    log:
        "{assembly}/log/merge_counts_pure_and_kaiju.log"
    benchmark:
        "{assembly}/benchmark/merge_counts_pure_and_kaiju.benchmark"
    run:
        os.makedirs(f"{assembly}/assembly_analysis/", exist_ok=True)
        shell("{Rscript} GEMINI/merge_counts_and_kaiju.R {input.counts_dp} {input.kaiju_dp} {assembly}/assembly_analysis/{assembly}.taxa_counts.tsv")


rule format_kaiju_output:
    input:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm"
    output:
       "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm.tsv"
    log:
        "{assembly}/log/format_kaiju_output.log"
    benchmark:
        "{assembly}/benchmark/format_kaiju_output.benchmark"
    shell:
        "{python3} GEMINI/format_kaiju_output_to_tab_seperated.py {input}"

rule kaiju_addTaxonNames:
    input:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref"
    output:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm"
    log:
        "{assembly}/log/kaiju_addTaxonNames.log"
    benchmark:
        "{assembly}/benchmark/kaiju_addTaxonNames.benchmark"
    threads: 4
    shell:
        "{kaiju_addTaxonNames} -i {input} -o {output} -t  {kaiju_nodes} -n {kaiju_names} "
        " -v -r superkingdom,phylum,class,order,family,genus,species"

rule kaiju_annotate:
    input:
        assembly_dir
    output:
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref"
    threads: 4
    log:
        "{assembly}/log/kaiju_annotate.log"
    benchmark:
        "{assembly}/benchmark/kaiju_annotate.benchmark"
    shell:
        "{kaiju} -t {kaiju_nodes} -v -f {kaiju_fmi} -z {threads} -i {input}  -o {output}"


rule paste_counts_table:
    input:
        "{assembly}/log/idxstats.done"
    output:
        "{assembly}/taxa_analysis/temp/{assembly}.counts.pure"
    run:
        counts_table_list = ','.join(expand("{assembly}/taxa_analysis/temp/{basename}.idxstats",
                                   assembly = assembly, basename = bam_basename))
        shell("{python3} GEMINI/paste_counts_table.py 3 {counts_table_list} {assembly}/taxa_analysis/temp/{assembly}.counts.pure")

rule idxstats:
    input:
        bam_list
    output:
        "{assembly}/log/idxstats.done"
    log:
        "{assembly}/log/idxstats.log"
    benchmark:
        "{assembly}/benchmark/idxstats.benchmark"
    threads: 4
    run:
        os.makedirs(f"{assembly}/taxa_analysis/temp/", exist_ok=True)
        for bam,basename in zip(bam_list,bam_basename):
            shell(f"{samtools} idxstats --threads {threads} {bam} > {assembly}/taxa_analysis/temp/{basename}.idxstats ")
        shell("touch {assembly}/log/idxstats.done")
















































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






