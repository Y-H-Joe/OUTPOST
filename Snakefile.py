import pandas as pd
import os
import sys
import GEMINI as ge
import itertools
import time

# %% command parameters
bwa = 'bwa'
kaiju = 'kaiju'
python3 = 'python3'
Rscript = 'Rscript'
kaiju_addTaxonNames = 'kaiju-addTaxonNames'
kaiju_nodes = '/home/yihang/software/kaijudb_nr_euk/nodes.dmp'
kaiju_fmi = '/home/yihang/software/kaijudb_nr_euk/nr_euk/kaiju_db_nr_euk.fmi'
kaiju_names = '/home/yihang/software/kaijudb_nr_euk/names.dmp'
samtools = 'samtools'
humann = '/home/yihang/anaconda3/envs/pipeline2/bin/humann3'
humann_renorm_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_renorm_table'
humann_regroup_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_regroup_table'
humann_join_tables = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_join_tables'
humann_split_stratified_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_split_stratified_table'
humann_rename_table = '/home/yihang/anaconda3/envs/pipeline2/bin/humann_rename_table'
lefse_format_input = 'lefse_format_input.py'
lefse_run = 'lefse_run.py'

# %% settings
taxa_level = ['taxaID','superkingdom','phylum','class','order','family','genus','species']
databases = ['rxn','eggnog','ko','level4ec','pfam']
paired = False
top = 30
cores = 32
config="GEMINI/GEMINI_config.tsv"

# %% GEMINI starts
df_config=pd.read_csv(config,sep='\t')
sample_list, fq_list, bam_list, assembly, assembly_dir, comparison_dict = ge.check_config(df_config)
bam_basename=[os.path.basename(x) for x in bam_list]
group_pair_list = list(itertools.combinations(comparison_dict.keys(),2))

rule all:
    input:
        f"{assembly}/log/heatmap_humann.done",
        f"{assembly}/log/heatmap_taxa.done",
        f"{assembly}/log/prepare_contig_table_from_counts_table.done",
        f"{assembly}/log/alpha_beta_diversity.done",
        f"{assembly}/log/lefse_humann.done",
        f"{assembly}/log/lefse_taxa.done",
        f"{assembly}/log/taxa_barplots.done",
        f"{assembly}/log/taxa_boxplot.done",
        f"{assembly}/log/process_contig_table.done"

# %% lefse
rule lefse_taxa:
    input:
        "{assembly}/log/rel_abun2lefse_taxa.done"
    output:
        "{assembly}/log/lefse_taxa.done"
    log:
        "{assembly}/log/lefse_taxa.log"
    benchmark:
        "{assembly}/benchmark/lefse_taxa.benchmark"
    run:
        wkdir = f"{assembly}/lefse_analysis/figs_taxa/"
        os.makedirs(wkdir, exist_ok=True)
        # unequal taxa to lefse barplot
        for group1,group2 in group_pair_list:
            dp_list = [f"{assembly}/lefse_analysis/taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.lefse.tsv"
                       for level in taxa_level]
            format_list = [dp.replace(".tsv",".format") for dp in dp_list]
            res_list = [dp.replace(".tsv",".res") for dp in dp_list]
            output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

            for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                try:
                    shell("{lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ")
                    shell("{lefse_run} {format_}  {res} > {log} 2>&1 ")
                    shell("{python3} GEMINI/barplots_from_lefse_res_file.py  {res} {output} > {log} 2>&1 ")
                except:
                    print("GEMINI: barplots_from_lefse_res_file: ",dp," has no features.skip.")

        # equal humann to lefse barplot
        for group1,group2 in group_pair_list:
            dp_list = [f"{assembly}/lefse_analysis/taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.lefse.tsv"
                       for level in taxa_level]
            format_list = [dp.replace(".tsv",".format") for dp in dp_list]
            res_list = [dp.replace(".tsv",".res") for dp in dp_list]
            output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

            for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                try:
                    shell("{lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ")
                    shell("{lefse_run} {format_}  {res} > {log} 2>&1 ")
                    shell("{python3} GEMINI/barplots_from_lefse_res_file.py  {res} {output} > {log} 2>&1 ")
                except:
                    print("GEMINI: barplots_from_lefse_res_file: ",dp," has no features.skip.")

        shell("touch {assembly}/log/lefse_taxa.done")

rule lefse_humann:
    input:
        "{assembly}/log/rel_abun2lefse_humann.done"
    output:
        "{assembly}/log/lefse_humann.done"
    log:
        "{assembly}/log/lefse_humann.log"
    benchmark:
        "{assembly}/benchmark/rel_abun2lefse_humann.benchmark"
    run:
        wkdir = f"{assembly}/lefse_analysis/figs_humann/"
        os.makedirs(wkdir, exist_ok=True)
        # unequal humann to lefse barplot
        for group1,group2 in group_pair_list:
            dp_list = [f"{assembly}/lefse_analysis/humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.lefse.tsv"
                       for database in databases]
            format_list = [dp.replace(".tsv",".format") for dp in dp_list]
            res_list = [dp.replace(".tsv",".res") for dp in dp_list]
            output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

            for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                try:
                    shell("{lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ")
                    shell("{lefse_run} {format_}  {res} > {log} 2>&1 ")
                    shell("{python3} GEMINI/barplots_from_lefse_res_file.py  {res} {output} > {log} 2>&1 ")
                except:
                    print("GEMINI: barplots_from_lefse_res_file: ",dp," has no features.skip.")

        # equal humann to lefse barplot
        for group1,group2 in group_pair_list:
            dp_list = [f"{assembly}/lefse_analysis/humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.lefse.tsv"
                       for database in databases]
            format_list = [dp.replace(".tsv",".format") for dp in dp_list]
            res_list = [dp.replace(".tsv",".res") for dp in dp_list]
            output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

            for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                shell("{lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ")
                shell("{lefse_run} {format_}  {res} > {log} 2>&1 ")
                shell("{python3} GEMINI/barplots_from_lefse_res_file.py  {res} {output} > {log} 2>&1 ")

        shell("touch {assembly}/log/lefse_humann.done")

rule rel_abun2lefse_taxa:
    input:
        "{assembly}/log/rel_abun_utest.done"
    output:
        "{assembly}/log/rel_abun2lefse_taxa.done"
    log:
        "{assembly}/log/rel_abun2lefse_taxa.log"
    benchmark:
        "{assembly}/benchmark/rel_abun2lefse_taxa.benchmark"
    run:
        # unequal taxa to lefse table
        for group1,group2 in group_pair_list:
            wkdir = f"{assembly}/lefse_analysis/taxa_{group1}_vs_{group2}/"
            os.makedirs(wkdir, exist_ok=True)
            dp_list = ','.join([f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.csv"
                                for level in taxa_level])
            output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                for dp in dp_list.split(',')])
            class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
            shell("{python3} GEMINI/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")
        # equal taxa to lefse table
        for group1,group2 in group_pair_list:
            wkdir = f"{assembly}/lefse_analysis/taxa_{group1}_vs_{group2}/"
            os.makedirs(wkdir, exist_ok=True)
            dp_list = ','.join([f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.csv"
                                for level in taxa_level])
            output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                for dp in dp_list.split(',')])
            class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
            shell("{python3} GEMINI/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")

        shell("touch {assembly}/log/rel_abun2lefse_taxa.done")

rule rel_abun2lefse_humann:
    input:
        "{assembly}/log/humann_utest.done"
    output:
        "{assembly}/log/rel_abun2lefse_humann.done"
    log:
        "{assembly}/log/rel_abun2lefse_humann.log"
    benchmark:
        "{assembly}/benchmark/rel_abun2lefse_humann.benchmark"
    run:
        # unequal humann to lefse table
        for group1,group2 in group_pair_list:
            wkdir = f"{assembly}/lefse_analysis/humann_{group1}_vs_{group2}/"
            os.makedirs(wkdir, exist_ok=True)
            dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.csv"
                                for database in databases])
            output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                for dp in dp_list.split(',')])
            class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
            shell("{python3} GEMINI/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")
        # equal humann to lefse table
        for group1,group2 in group_pair_list:
            wkdir = f"{assembly}/lefse_analysis/humann_{group1}_vs_{group2}/"
            os.makedirs(wkdir, exist_ok=True)
            dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.csv"
                                for database in databases])
            output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                for dp in dp_list.split(',')])
            class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
            shell("{python3} GEMINI/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")

        shell("touch {assembly}/log/rel_abun2lefse_humann.done")

# %% heatmap
rule heatmap_taxa:
    input:
        "{assembly}/log/scale_rel_abun_table.done"
    output:
        "{assembly}/log/heatmap_taxa.done"
    log:
        "{assembly}/log/heatmap_taxa.log"
    benchmark:
        "{assembly}/benchmark/heatmap_taxa.benchmark"
    run:
        wkdir = f"{assembly}/taxa_analysis/figs/"
        os.makedirs(wkdir, exist_ok=True)
        # top taxa heatmap
        data_type = 'taxa'
        for level in taxa_level:
            dp = f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.top{top}.fillmin.log10.csv"
            output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
            index_dp = f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.top{top}.fillmin.log10.index"
            with open(index_dp,'w') as w:
                for sample in sample_list:
                    w.write(sample + '\n')
            try:
                shell("{Rscript} GEMINI/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
            except: pass
        # top unequal taxa heatmap
        for group1,group2 in group_pair_list:
            for level in taxa_level:
                dp = f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top30.fillmin.log10.csv"
                output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                index_dp = f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top30.fillmin.log10.index"
                with open(index_dp,'w') as w:
                    for sample in comparison_dict[group1]:
                        w.write(f"{sample}_{group1}\n")
                    for sample in comparison_dict[group2]:
                        w.write(f"{sample}_{group2}\n")
                try:
                    shell("{Rscript} GEMINI/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                except: pass
        # top equal taxa heatmap
        for group1,group2 in group_pair_list:
            for level in taxa_level:
                dp = f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top30.fillmin.log10.csv"
                output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                index_dp = f"{assembly}/taxa_analysis/top_taxa_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top30.fillmin.log10.index"
                with open(index_dp,'w') as w:
                    for sample in comparison_dict[group1]:
                        w.write(f"{sample}_{group1}\n")
                    for sample in comparison_dict[group2]:
                        w.write(f"{sample}_{group2}\n")
                try:
                    shell("{Rscript} GEMINI/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                except: pass
        shell("touch {assembly}/log/heatmap_taxa.done")

rule heatmap_humann:
    input:
        "{assembly}/log/extract_top_humann.done"
    output:
        "{assembly}/log/heatmap_humann.done"
    log:
        "{assembly}/log/heatmap_humann.log"
    benchmark:
        "{assembly}/benchmark/heatmap_humann.benchmark"
    run:
        wkdir = f"{assembly}/metabolism_analysis/figs/"
        os.makedirs(wkdir, exist_ok=True)

        with open(f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_cpm_{databases[0]}_unstratified.named.tsv",'r') as r:
            line = r.readline()
            humann_sample_list = [x.replace('_Abundance-RPKs','') for x in line.strip().split('\t')[1:]]

        # top humann heatmap
        data_type = 'humann'
        for database in databases:
            dp = f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.csv"
            output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
            index_dp = f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.index"
            with open(index_dp,'w') as w:
                for sample in humann_sample_list:
                    w.write(sample + '\n')
            try:
                shell("{Rscript} GEMINI/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
            except: pass
        # top unequal humann heatmap
        for group1,group2 in group_pair_list:
            for database in databases:
                dp = f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.csv"
                output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                index_dp = f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.index"
                with open(index_dp,'w') as w:
                    for sample in comparison_dict[group1]:
                        w.write(f"{sample}_{group1}\n")
                    for sample in comparison_dict[group2]:
                        w.write(f"{sample}_{group2}\n")
                try:
                    shell("{Rscript} GEMINI/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                except: pass
        # top equal humann heatmap
        for group1,group2 in group_pair_list:
            for database in databases:
                dp = f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.csv"
                output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                index_dp = f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.index"
                with open(index_dp,'w') as w:
                    for sample in comparison_dict[group1]:
                        w.write(f"{sample}_{group1}\n")
                    for sample in comparison_dict[group2]:
                        w.write(f"{sample}_{group2}\n")
                try:
                    shell("{Rscript} GEMINI/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                except: pass
        shell("touch {assembly}/log/heatmap_humann.done")

# %% humann
rule extract_top_humann:
    input:
        "{assembly}/log/humann_utest.done"
    output:
        "{assembly}/log/extract_top_humann.done"
    log:
        "{assembly}/log/extract_top_humann.log"
    benchmark:
        "{assembly}/benchmark/extract_top_humann.benchmark"
    threads: cores
    run:
        # top abundant
        dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.csv"
                       for database in databases])
        output_name_list = ','.join([f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.csv"
                       for database in databases])
        error_log = os.getcwd() + f"/{assembly}/log/extract_top_humann.error"
        shell("{python3} GEMINI/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
              " {top} {error_log} > {log} 2>&1 ")
        result = ge.wait_unti_file_exists(output_name_list.split(','),error_log)
        assert result is True, "GEMINI: extract_top_humann has errors. exit."

        # top unequal abundant
        output_list = []
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}", exist_ok=True)
            dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.csv"
                       for database in databases])
            output_name_list = ','.join([f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.csv"
                    for database in databases])
            error_log = os.getcwd() + f"/{assembly}/log/extract_top_humann_top_unequal_abundant.error"
            output_list.append(f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.csv")
            shell("nohup {python3} GEMINI/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
              " {top} {error_log} > {log}_{group1}_vs_{group2}_unequal 2>&1 &")
            time.sleep(60)
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: extract_top_humann_top_unequal_abundant has errors. exit."

       # top equal abundant
        output_list = []
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}", exist_ok=True)
            dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.csv"
                       for database in databases])
            output_name_list = ','.join([f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.csv"
                    for database in databases])
            error_log = os.getcwd() + f"/{assembly}/log/extract_top_humann_top_equal_abundant.error"
            output_list.append(f"{assembly}/metabolism_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.csv")
            shell("nohup {python3} GEMINI/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
              " {top} {error_log} > {log}_{group1}_vs_{group2}_equal 2>&1 &")
            time.sleep(60)
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: extract_top_humann_top_equal_abundant has errors. exit."

        shell("touch {assembly}/log/extract_top_humann.done")


rule humann_utest:
    input:
        "{assembly}/log/humann2rel_abun.done"
    output:
        "{assembly}/log/humann_utest.done"
    log:
        "{assembly}/log/humann_utest.log"
    benchmark:
        "{assembly}/benchmark/humann_utest.benchmark"
    threads: cores
    run:
        output_list = []
        dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.csv"
                   for database in databases])
        with open(f"{assembly}/metabolism_analysis/humann3/output/allSamples_genefamilies_uniref90names_cpm_{databases[0]}_unstratified.named.tsv",'r') as r:
            line = r.readline()
            humann_sample_list = [x.replace('_Abundance-RPKs','') for x in line.strip().split('\t')[1:]]
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}", exist_ok=True)

            group1_index = ','.join([str(humann_sample_list.index(sample)) for sample in comparison_dict[group1]])
            group2_index = ','.join([str(humann_sample_list.index(sample)) for sample in comparison_dict[group2]])
            prefix_list = ','.join([f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}"
                       for database in databases])
            output = f"{assembly}/metabolism_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.{group1}_vs_{group2}.ave_change.unequal.csv"
            output_list.append(output)
            error_log = os.getcwd() + f"/{assembly}/log/humann_utest.error"
            try:
                os.remove(error_log) # clean former residual error log
            except:
                pass
            shell("nohup {python3} GEMINI/rel_abun_utest.py {dp_list} {group1_index} {group1} "
                  " {group2_index} {group2} {prefix_list} {paired} {error_log} > {log}_{group1}_vs_{group2} 2>&1 &")
            time.sleep(60)
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: rel_abun_utest has errors. exit."
        shell("touch {assembly}/log/humann_utest.done")

rule humann2rel_abun:
    input:
        "{assembly}/log/humann_output.done"
    output:
        "{assembly}/log/humann2rel_abun.done"
    log:
        "{assembly}/log/humann2rel_abun.log"
    benchmark:
        "{assembly}/benchmark/humann2rel_abun.benchmark"
    run:
        for group in list(comparison_dict.keys()) + ['allSamples']:
            dp_list = ','.join([f"{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}_unstratified.named.tsv"
                       for database in databases])
            output_list = ','.join([f"{assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.csv"
                       for database in databases])
            shell("{python3} GEMINI/humann2rel_abun.py {dp_list} {output_list} > {log} 2>&1 ")
        shell("touch {assembly}/log/humann2rel_abun.done")

rule humann_output:
    input:
        "{assembly}/log/humann_group.done"
    output:
        "{assembly}/log/humann_output.done"
    log:
        "{assembly}/log/humann_output.log"
    benchmark:
        "{assembly}/benchmark/humann_output.benchmark"
    run:
        os.makedirs(f"{assembly}/metabolism_analysis/humann3/output/", exist_ok=True)
        for group in list(comparison_dict.keys()) + ['allSamples']:
            # join tables for each group
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                  "--file_name _genefamilies_relab.tsv  > {log} 2>&1 ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                  "--file_name _genefamilies_cpm.tsv  > {log} 2>&1 ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathabundance_relab.tsv "
                  "--file_name _pathabundance_relab.tsv  > {log} 2>&1 ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathabundance_cpm.tsv "
                  "--file_name _pathabundance_cpm.tsv  > {log} 2>&1 ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathcoverage_relab.tsv "
                  "--file_name _pathcoverage_relab.tsv  > {log} 2>&1 ")
            shell("{humann_join_tables} -i {assembly}/metabolism_analysis/humann3/{group}/ "
                  " -o {assembly}/metabolism_analysis/humann3/output/{group}_pathcoverage_cpm.tsv "
                  "--file_name _pathcoverage_cpm.tsv  > {log} 2>&1 ")

            # regroup tables for each group
            for database in databases:
                shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                      " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")
                shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                      " --output {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")

            # stratify tables for each group
            for database in databases:
                shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_{database}.tsv "
                      "-o {assembly}/metabolism_analysis/humann3/output/ > {log} 2>&1 ")
                shell("{humann_split_stratified_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}.tsv "
                      "-o {assembly}/metabolism_analysis/humann3/output/ > {log} 2>&1 ")

            # rename tables for each group
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.tsv "
                  " --names kegg-orthology -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.tsv "
                  " --names eggnog -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.tsv "
                  " --names ec -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.tsv "
                  " --names pfam -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.tsv "
                  " --names metacyc-rxn -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko_unstratified.tsv "
                  " --names kegg-orthology -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog_unstratified.tsv "
                  " --names eggnog -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec_unstratified.tsv "
                  " --names ec -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam_unstratified.tsv "
                  " --names pfam -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam_unstratified.named.tsv > {log} 2>&1 ")
            shell("{humann_rename_table} -i {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.tsv "
                  " --names metacyc-rxn -o {assembly}/metabolism_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.named.tsv > {log} 2>&1 ")
        shell("touch {assembly}/log/humann_output.done")

rule humann_group:
    input:
        "{assembly}/log/humann_annotate.done",
    output:
        "{assembly}/log/humann_group.done"
    log:
        "{assembly}/log/humann_group.log"
    benchmark:
        "{assembly}/benchmark/humann_group.benchmark"
    run:
        for group,samples in comparison_dict.items():
            os.makedirs(f"{assembly}/metabolism_analysis/humann3/{group}", exist_ok=True)
            for sample in samples:
                os.system(f"cp {assembly}/metabolism_analysis/humann3/ori_results/{sample}*tsv  {assembly}/metabolism_analysis/humann3/{group} ")
        # allSamples
        os.makedirs(f"{assembly}/metabolism_analysis/humann3/allSamples", exist_ok=True)
        for sample in sample_list:
            os.system(f"cp {assembly}/metabolism_analysis/humann3/ori_results/{sample}*tsv  {assembly}/metabolism_analysis/humann3/allSamples ")
        shell("touch {assembly}/log/humann_group.done")

rule humann_annotate:
    input:
        "{assembly}/log/rename_humann_ori_output.done"
    output:
        "{assembly}/log/humann_annotate.done"
    log:
        "{assembly}/log/humann_annotate.log"
    benchmark:
        "{assembly}/benchmark/humann_annotate.benchmark"
    run:
        for sample in sample_list:
            # normalize
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                  " --units relab > {log} 2>&1 ")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                  " --units cpm > {log} 2>&1 ")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance_relab.tsv "
                  " --units relab > {log} 2>&1 ")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance_cpm.tsv "
                  " --units cpm > {log} 2>&1 ")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage_relab.tsv "
                  " --units relab > {log} 2>&1 ")
            shell("{humann_renorm_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv "
                  " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage_cpm.tsv "
                  " --units cpm > {log} 2>&1 ")
            # regroup
            for database in databases:
                shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                      " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_cpm_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")
                shell("{humann_regroup_table} --input {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                      " --output {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies_relab_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")

        shell("touch {assembly}/log/humann_annotate.done")

rule rename_humann_ori_output:
    input:
        "{assembly}/log/humann_init.done"
    output:
        "{assembly}/log/rename_humann_ori_output.done"
    log:
        "{assembly}/log/rename_humann_ori_output.log"
    benchmark:
        "{assembly}/benchmark/humann_init.benchmark"
    run:
        for sample,fq in zip(sample_list, fq_list):
            # get the humann renamed file name, taken from humann source code
            fq_humann = os.path.basename(fq)
            # Remove gzip extension if present
            if re.search('.gz$',fq_humann):
                fq_humann='.'.join(fq_humann.split('.')[:-1])
            # Remove input file extension if present
            if '.' in fq_humann:
                fq_humann='.'.join(fq_humann.split('.')[:-1])

            shell("cp {assembly}/metabolism_analysis/humann3/ori_results/{fq_humann}_genefamilies.tsv "
                  " {assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv")
            shell("cp {assembly}/metabolism_analysis/humann3/ori_results/{fq_humann}_pathabundance.tsv "
                  " {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv")
            shell("cp {assembly}/metabolism_analysis/humann3/ori_results/{fq_humann}_pathcoverage.tsv "
                  " {assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv")

            file1 = f"{assembly}/metabolism_analysis/humann3/ori_results/{fq_humann}_genefamilies.tsv"
            file2 = f"{assembly}/metabolism_analysis/humann3/ori_results/{sample}_genefamilies.tsv"
            with open(file1,'r') as r, open(file2,'w') as w:
                lines = r.readlines()
                lines[0] = "# Gene Family\t" + sample + "_Abundance-RPKs\n"
                w.write(lines)

            file1 = f"{assembly}/metabolism_analysis/humann3/ori_results/{fq_humann}_pathabundance.tsv"
            file2 = f"{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathabundance.tsv"
            with open(file1,'r') as r, open(file2,'w') as w:
                lines = r.readlines()
                lines[0] = "# Pathway\t" + sample + "_Abundance\n"
                w.write(lines)

            file1 = f"{assembly}/metabolism_analysis/humann3/ori_results/{fq_humann}_pathcoverage.tsv"
            file2 = f"{assembly}/metabolism_analysis/humann3/ori_results/{sample}_pathcoverage.tsv"
            with open(file1,'r') as r, open(file2,'w') as w:
                lines = r.readlines()
                lines[0] = "# Pathway\t" + sample + "_Coverage\n"
                w.write(lines)
        shell("touch {assembly}/log/rename_humann_ori_output.done")

rule humann_init:
    input:
        fq_list
    output:
        "{assembly}/log/humann_init.done"
    log:
        "{assembly}/log/humann_init.log"
    benchmark:
        "{assembly}/benchmark/humann_init.benchmark"
    threads: cores
    run:
        os.makedirs(f"{assembly}/metabolism_analysis/humann3/ori_results/", exist_ok=True)
        for fq in fq_list:
            shell("{humann} --input {fq}  --output {assembly}/metabolism_analysis/humann3/ori_results/ "
            " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
            " --threads {threads} --memory-use maximum > {log} 2>&1 ")
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
    threads: cores
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

rule taxa_barplots:
    input:
        "{assembly}/log/extract_top_taxa.done"
    output:
        "{assembly}/log/taxa_barplots.done"
    log:
        "{assembly}/log/taxa_barplots.log"
    benchmark:
        "{assembly}/benchmark/taxa_barplots.benchmark"
    run:
        os.makedirs(f"{assembly}/taxa_analysis/figs/", exist_ok=True)
        for level in taxa_level:
            dp = f"{assembly}/taxa_analysis/{assembly}.taxa_counts.rel_abun.{level}.rmU.top{top}.csv"
            new_dp = dp.replace(".csv",".addOthers.csv")
            output = f"{assembly}/taxa_analysis/figs/" + os.path.basename(dp).replace(".csv",".barplot.pdf")
            df = pd.read_csv(dp, index_col = 0)
            df['Others'] = 1 - df.sum(axis = 1)
            df.to_csv(new_dp)
            try:
                shell("{Rscript} GEMINI/barplots.R {new_dp} {output} > {log} 2>&1")
            except:
                pass
        shell("touch {assembly}/log/taxa_barplots.done")

rule extract_top_taxa:
    input:
        "{assembly}/log/rel_abun_utest.done"
    output:
        "{assembly}/log/extract_top_taxa.done"
    log:
        "{assembly}/log/extract_top_taxa.log"
    benchmark:
        "{assembly}/benchmark/extract_top_taxa.benchmark"
    threads: cores
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
            time.sleep(60)
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
            time.sleep(60)
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: extract_top_taxa_top_equal_abundant has errors. exit."

        shell("touch {assembly}/log/extract_top_taxa.done")

rule taxa_boxplot:
    input:
        "{assembly}/log/rel_abun_utest.done"
    output:
        "{assembly}/log/taxa_boxplot.done"
    log:
        "{assembly}/log/taxa_boxplot.log"
    benchmark:
        "{assembly}/benchmark/taxa_boxplot.benchmark"
    threads: cores
    run:
        output_list = []
        for group1,group2 in group_pair_list:
            os.makedirs(f"{assembly}/taxa_analysis/boxplot_{group1}_vs_{group2}", exist_ok=True)
            group_pair = ','.join([group1,group2])
            groups = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
            dp_list = [f"{assembly}/taxa_analysis/utest_{group1}_vs_{group2}/{assembly}.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.csv" for level in taxa_level]
            for dp in dp_list:
                output = f"{assembly}/taxa_analysis/boxplot_{group1}_vs_{group2}/{os.path.basename(dp).replace('.csv','')}"
            try:
                shell("{Rscript} GEMINI/boxplot.R {dp} {output} {groups} {group_pair} > {log} 2>&1")
            except:
                pass
        shell("touch {assembly}/log/taxa_boxplot.done")

rule rel_abun_utest:
    input:
        "{assembly}/log/counts_table2rel_abun.done"
    output:
        "{assembly}/log/rel_abun_utest.done"
    log:
        "{assembly}/log/rel_abun_utest.log"
    benchmark:
        "{assembly}/benchmark/rel_abun_utest.benchmark"
    threads: cores
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
            time.sleep(60)
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: rel_abun_utest has errors. exit."
        shell("touch {assembly}/log/rel_abun_utest.done")

rule counts_table2rel_abun:
    input:
        "{assembly}/log/merge_counts_pure_and_kaiju.done",
        real = "{assembly}/assembly_analysis/{assembly}.taxa_counts.tsv"
    output:
        "{assembly}/log/counts_table2rel_abun.done"
    log:
        "{assembly}/log/counts_table2rel_abun.log"
    benchmark:
        "{assembly}/benchmark/counts_table2rel_abun.benchmark"
    run:
        prefix = f"{assembly}/taxa_analysis/{assembly}.taxa_counts"
        samples = ','.join(sample_list)
        shell("{python3} GEMINI/counts_table2rel_abun.py {input.real} {prefix} {samples} > {log} 2>&1 ")
        shell("touch {assembly}/log/counts_table2rel_abun.done")

rule process_contig_table:
    input:
        "{assembly}/log/prepare_contig_table_from_counts_table.done",
    output:
        "{assembly}/log/process_contig_table.done"
    log:
        "{assembly}/log/process_contig_table.log"
    benchmark:
        "{assembly}/benchmark/process_contig_table.benchmark"
    run:
        for group1, group2 in group_pair_list:
            dp = f"{assembly}/assembly_analysis/{assembly}.{group1}_vs_{group2}.contig_table.tsv"
            samples = ','.join(sample_list)
            output = f"{assembly}/assembly_analysis/{assembly}.{group1}_vs_{group2}.contig_table.processed.tsv"
            shell("{python3} GEMINI/process_contig_table.py {dp} {samples} {output} > {log} 2>&1 ")
        shell("touch {assembly}/log/process_contig_table.done")

rule prepare_contig_table_from_counts_table:
    input:
        "{assembly}/log/merge_counts_pure_and_kaiju.done",
        real = "{assembly}/assembly_analysis/{assembly}.taxa_counts.tsv"
    output:
        "{assembly}/log/prepare_contig_table_from_counts_table.done"
    log:
        "{assembly}/log/prepare_contig_table_from_counts_table.log"
    benchmark:
        "{assembly}/benchmark/prepare_contig_table_from_counts_table.benchmark"
    threads: cores
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
            shell("nohup {python3} GEMINI/prepare_contig_table_from_counts_table.py {input.real} "
                  " {group1_index} {group2_index} {output} {samples} {error_log} > {log} 2>&1 &")
            time.sleep(60)
        result = ge.wait_unti_file_exists(output_list,error_log)
        assert result is True, "GEMINI: prepare_contig_table_from_counts_table has errors. exit."
        shell("touch {assembly}/log/prepare_contig_table_from_counts_table.done")

rule merge_counts_pure_and_kaiju:
    input:
        "{assembly}/log/format_kaiju_output.done",
        "{assembly}/log/paste_counts_table.done",
        counts_dp = "{assembly}/taxa_analysis/temp/{assembly}.counts.pure",
        kaiju_dp = "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm.tsv"
    output:
        "{assembly}/log/merge_counts_pure_and_kaiju.done",
        "{assembly}/assembly_analysis/{assembly}.taxa_counts.tsv"
    log:
        "{assembly}/log/merge_counts_pure_and_kaiju.log"
    benchmark:
        "{assembly}/benchmark/merge_counts_pure_and_kaiju.benchmark"
    run:
        os.makedirs(f"{assembly}/assembly_analysis/", exist_ok=True)
        shell("{Rscript} GEMINI/merge_counts_and_kaiju.R {input.counts_dp} {input.kaiju_dp} {assembly}/assembly_analysis/{assembly}.taxa_counts.tsv")
        shell("touch {assembly}/log/merge_counts_pure_and_kaiju.done")

rule format_kaiju_output:
    input:
        "{assembly}/log/kaiju_addTaxonNames.done",
        real = "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm"
    output:
        "{assembly}/log/format_kaiju_output.done",
        "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm.tsv"
    log:
        "{assembly}/log/format_kaiju_output.log"
    benchmark:
        "{assembly}/benchmark/format_kaiju_output.benchmark"
    run:
        shell("{python3} GEMINI/format_kaiju_output_to_tab_seperated.py {input.real}")
        shell("touch {assembly}/log/format_kaiju_output.done")

rule kaiju_addTaxonNames:
    input:
        "{assembly}/log/kaiju_annotate.done",
        real = "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref"
    output:
        "{assembly}/log/kaiju_addTaxonNames.done",
        real = "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref.nm"
    log:
        "{assembly}/log/kaiju_addTaxonNames.log"
    benchmark:
        "{assembly}/benchmark/kaiju_addTaxonNames.benchmark"
    threads: cores
    run:
        shell("{kaiju_addTaxonNames} -i {input.real} -o {output.real} -t  {kaiju_nodes} -n {kaiju_names} "
        " -v -r superkingdom,phylum,class,order,family,genus,species > {log} 2>&1 ")
        shell("touch {assembly}/log/kaiju_addTaxonNames.done")

rule kaiju_annotate:
    input:
        assembly_dir
    output:
        "{assembly}/log/kaiju_annotate.done",
        real = "{assembly}/taxa_analysis/kaiju/{assembly}_kaiju.ref"
    threads: cores
    log:
        "{assembly}/log/kaiju_annotate.log"
    benchmark:
        "{assembly}/benchmark/kaiju_annotate.benchmark"
    run:
        shell("{kaiju} -t {kaiju_nodes} -v -f {kaiju_fmi} -z {threads} -i {input}  -o {output.real} > {log} 2>&1 ")
        shell("touch {assembly}/log/kaiju_annotate.done")

rule paste_counts_table:
    input:
        "{assembly}/log/idxstats.done"
    output:
        "{assembly}/log/paste_counts_table.done",
        "{assembly}/taxa_analysis/temp/{assembly}.counts.pure"
    log:
        "{assembly}/log/paste_counts_table.log"
    benchmark:
        "{assembly}/benchmark/paste_counts_table.benchmark"
    run:
        counts_table_list = ','.join(expand("{assembly}/taxa_analysis/temp/{basename}.idxstats",
                                   assembly = assembly, basename = bam_basename))
        shell("{python3} GEMINI/paste_counts_table.py 3 {counts_table_list} {assembly}/taxa_analysis/temp/{assembly}.counts.pure")
        shell("touch {assembly}/log/paste_counts_table.done")

rule idxstats:
    input:
        bam_list
    output:
        "{assembly}/log/idxstats.done"
    log:
        "{assembly}/log/idxstats.log"
    benchmark:
        "{assembly}/benchmark/idxstats.benchmark"
    threads: cores
    run:
        os.makedirs(f"{assembly}/taxa_analysis/temp/", exist_ok=True)
        for bam,basename in zip(bam_list,bam_basename):
            shell(f"{samtools} idxstats --threads {threads} {bam} > {assembly}/taxa_analysis/temp/{basename}.idxstats ")
        shell("touch {assembly}/log/idxstats.done")



