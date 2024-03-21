################################################################################
###!!!DO NOT MODIFY THE FOLLOWING CODES UNLESS YOU KNOW WHAT YOU ARE DOING!!!###
################################################################################

# %% main body
import pandas as pd
import numpy as np
import os
import sys
from OUTPOST import index_genome, wait_until_file_exists, check_config
from OUTPOST.check_snakefile_config import *
import itertools
import time
import shutil

# %% prepare
try:
    df_config=pd.read_csv(OUTPOST_config,sep='\t')
except:
    sys.exit(f"OUTPOST: read {OUTPOST_config} error. Make sure file exists in tab separated format.")
sample_list, refgenome_list, r1_fq_list, r2_fq_list, se_fq_list, bam_dir_list, assembly_list, assembly_dir_list, batch_list, comparison_dict = check_config.check_config(df_config, rm_batch_effect)

bam_basename=[os.path.basename(x) for x in bam_dir_list]
group_pair_list = list(itertools.combinations(comparison_dict.keys(),2))
os.makedirs(output_dir, exist_ok = True)

if len(r1_fq_list + se_fq_list) > 0:
    fq_list = []
    for i,sample in enumerate(sample_list):
        fq_list.append(f"{output_dir}/data/{sample}_nonhostvirus.fq")
    fq_humann_list = []
    for sample,fq in zip(sample_list, fq_list):
        # get the humann renamed file name, taken from humann source code
        fq_humann = os.path.basename(fq)
        # Remove gzip extension if present
        if re.search('.gz$',fq_humann):
            fq_humann='.'.join(fq_humann.split('.')[:-1])
        # Remove input file extension if present
        if '.' in fq_humann:
            fq_humann='.'.join(fq_humann.split('.')[:-1])
        fq_humann_list.append(fq_humann)
    if assemble_contigs:
        # update assembly information
        # overwrite original assembly, because we cannot handle multiassemble vs one sample
        assembly_list = ["outpost_assembly"] * len(sample_list)
        assembly_dir_list = [f"{output_dir}/assembly_analysis/outpost_contigs/outpost_nonrd_contigs.fasta"] * len(sample_list)
process_batch_size = min(process_batch_size, len(sample_list))

# %% starts
shell("ulimit -s 65535")
rule all:
    input:
        f"{output_dir}/log/OUTPOST_report.done"

# %% report
OUTPOST_report_input = []
if len(r1_fq_list + se_fq_list) > 0:
    OUTPOST_report_input += [
    f"{output_dir}/log/metaphlan_krona.done",
    f"{output_dir}/log/metaphlan_graphlan.done",
    f"{output_dir}/log/metaphlan_heatmap.done",
    f"{output_dir}/log/metaphlan_diverisy.done",
    f"{output_dir}/log/lefse_humann.done",
    f"{output_dir}/log/heatmap_humann.done"]
if len(assembly_dir_list) > 0:
    if len(r1_fq_list + se_fq_list) > 0:
        OUTPOST_report_input += [f"{output_dir}/log/salmon_gene_quantify.done"]
        
    OUTPOST_report_input += [f"{output_dir}/log/eggmapper_gene_annotation.done",
    f"{output_dir}/log/rgi_gene_annotation.done",
    f"{output_dir}/log/gtdbtk_genome_annotation.done",
    f"{output_dir}/log/prodigal_gene_prediction.done",
    f"{output_dir}/log/metagenemark.done"]
    if (len(r1_fq_list + se_fq_list) > 0 and len(assembly_list) > 0) or (len(bam_dir_list) > 0):
        OUTPOST_report_input += [f"{output_dir}/log/biomarker_summary.done",
        f"{output_dir}/log/alpha_beta_diversity.done",
        f"{output_dir}/log/visualize_batch_effect.done",
        f"{output_dir}/log/heatmap_taxa.done",
        f"{output_dir}/log/taxa_barplots.done",
        f"{output_dir}/log/taxa_boxplot.done",
        f"{output_dir}/log/assembly_qtest.done"]
    OUTPOST_report_input += [f"{output_dir}/log/concat_kaiju_output.done"]

# rule OUTPOST_report:
#     input:
#         OUTPOST_report_input
#     output:
#         "{output_dir}/log/OUTPOST_report.done"
#     log:
#         "{output_dir}/log/OUTPOST_report.log"
#     benchmark:
#         "{output_dir}/benchmark/OUTPOST_report.benchmark"
#     run:
#         wkdir = f"{output_dir}/report/"
#         os.makedirs(wkdir, exist_ok=True)
#         os.makedirs(f"{output_dir}/taxonomy_analysis/counts_tables/", exist_ok=True)
        
#         # clean taxonomy_analysis tables
#         try:
#             shell("mv {output_dir}/taxonomy_analysis/*index {output_dir}/taxonomy_analysis/counts_tables/")
#             shell("mv {output_dir}/taxonomy_analysis/*csv {output_dir}/taxonomy_analysis/counts_tables/")
#         except:
#             pass
        
#         # generate reports
#         for group1,group2 in group_pair_list:
#             output_filename = f"{output_dir}/report/OUTPOST_report_{group1}_vs_{group2}.pdf"
#             explanation_ori = f"OUTPOST/content_original.csv"
#             visualization_ori = f"OUTPOST/visualization_example_path_original.csv"
#             resultspath = f"{assembly}"
#             group1VSgroup2 = f"{group1}_vs_{group2}"
#             shell("{python3} OUTPOST/OUTPOST_report.py  {output_filename} {assembly} {explanation_ori} {visualization_ori} {resultspath} {group1VSgroup2} > {log} 2>&1 ")

#         shell("touch {output_dir}/log/OUTPOST_report.done")

rule OUTPOST_report:
    input:
        OUTPOST_report_input
    output:
        "{output_dir}/log/OUTPOST_report.done"
    log:
        "{output_dir}/log/OUTPOST_report.log"
    benchmark:
        "{output_dir}/benchmark/OUTPOST_report.benchmark"
    run:
        wkdir = f"{output_dir}/report/"
        os.makedirs(wkdir, exist_ok=True)
        os.makedirs(f"{output_dir}/taxonomy_analysis/counts_tables/", exist_ok=True)
        
        # clean taxonomy_analysis tables
        try:
            shell("mv {output_dir}/taxonomy_analysis/*index {output_dir}/taxonomy_analysis/counts_tables/")
            shell("mv {output_dir}/taxonomy_analysis/*csv {output_dir}/taxonomy_analysis/counts_tables/")
        except:
            pass
        
        # generate reports
        shell("{python3} OUTPOST/OUTPOST_report_html.py {output_dir} {snakefile_config} {OUTPOST_config} {python3} OUTPOST/OUTPOST_report_template.html")
        
        shell("touch {output_dir}/log/OUTPOST_report.done")

# %% reads
if len(r1_fq_list + se_fq_list) > 0:
    #%% metaphlan
    rule metaphlan_krona:
        input:
            "{output_dir}/log/metaphlan_init.done"
        output:
            "{output_dir}/log/metaphlan_krona.done"
        log:
            "{output_dir}/log/metaphlan_krona.log"
        benchmark:
            "{output_dir}/benchmark/metaphlan_krona.benchmark"
        threads: 1
        run:
            krona_tsvs = []
            for i, sample in enumerate(sample_list):
                metaphlan_output = rf"{output_dir}/metaphlan_analysis/taxonomy/metaphlan_{sample}_taxa.txt"
                krona_tsv = rf"{output_dir}/metaphlan_analysis/krona/metaphlan_{sample}_krona.tsv"
                krona_tsvs.append(krona_tsv)
                shell("{python3}  OUTPOST/outpost_metaphlan2krona.py -p {metaphlan_output} -k {krona_tsv}" )
            # pass multiple tsv to krona
            krona_plot = rf"{output_dir}/metaphlan_analysis/figs/krona.html"
            # pass multiple tsv to kiimporttext will generate embedded html
            krona_tsvs_list_str = ' '.join(krona_tsvs)
            shell("{ktimporttext} {krona_tsvs_list_str}  -o {krona_plot} ")
            shell("touch {output_dir}/log/metaphlan_krona.done")
    
    
    rule metaphlan_graphlan:
        input:
            "{output_dir}/log/metaphlan_init.done"
        output:
            "{output_dir}/log/metaphlan_graphlan.done"
        log:
            "{output_dir}/log/metaphlan_graphlan.log"
        benchmark:
            "{output_dir}/benchmark/metaphlan_graphlan.benchmark"
        threads: 1
        run:
            output_merged = rf"{output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt"
            tree_file = rf"{output_dir}/metaphlan_analysis/taxonomy/merged_abundance.tree.txt"
            tree_anno = rf"{output_dir}/metaphlan_analysis/taxonomy/merged_abundance.tree.annot.txt"
            shell("{python3} {export2graphlan} --skip_rows 0 -i {output_merged} \
                  --tree {tree_file} --annotation {tree_anno} \
                  --most_abundant 1000 --abundance_threshold 0.01 --least_biomarkers 10 \
                  --annotations 1,2,3,4,5,6 --external_annotations 7")
                  
            graphlan_xml = rf"{output_dir}/metaphlan_analysis/taxonomy/merged_abundance.graphlan.xml"
            shell("python {graphlan_annotate} --annot {tree_anno} {tree_file} {graphlan_xml}")
            
            graphlan_pdf = rf"{output_dir}/metaphlan_analysis/figs/graphlan.pdf"
            shell("{graphlan} {graphlan_xml} {graphlan_pdf} --size 20 ")
            shell("touch {output_dir}/log/metaphlan_graphlan.done")
    
    
    rule metaphlan_heatmap:
        input:
            "{output_dir}/log/metaphlan_init.done"
        output:
            "{output_dir}/log/metaphlan_heatmap.done"
        log:
            "{output_dir}/log/metaphlan_heatmap.log"
        benchmark:
            "{output_dir}/benchmark/metaphlan_heatmap.benchmark"
        threads: 1
        run:
            output_merged = rf"{output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt"
            heatmap = rf"{output_dir}/metaphlan_analysis/figs/metaphlan_merge_taxa.heatmap.pdf"
            shell("{python3} OUTPOST/metaphlan_hclust2.py -i {output_merged} \
                  -o {heatmap} --skip_rows 0 --ftop 50 --f_dist_f correlation \
                 --s_dist_f euclidean --cell_aspect_ratio 1 -s --flabel_size 4 \
                 --sname_row 0 --max_flabel_len 200 --metadata_height 0.075  \
                  --minv 0.01  --slinkage average  ")
            shell("touch {output_dir}/log/metaphlan_heatmap.done")
    
    
    rule metaphlan_diverisy:
        input:
            "{output_dir}/log/metaphlan_init.done"
        output:
            "{output_dir}/log/metaphlan_diverisy.done"
        log:
            "{output_dir}/log/metaphlan_diverisy.log"
        benchmark:
            "{output_dir}/benchmark/metaphlan_diverisy.benchmark"
        threads: 1
        run:
            # beta diversity
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d beta -m bray-curtis -o {output_dir}/metaphlan_analysis/diversity")
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d beta -m jaccard -o {output_dir}/metaphlan_analysis/diversity")
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d beta -m clr -o {output_dir}/metaphlan_analysis/diversity")
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d beta -m aitchison -o {output_dir}/metaphlan_analysis/diversity")
            # alpha diversity
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d alpha -m richness -o {output_dir}/metaphlan_analysis/diversity")
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d alpha -m shannon -o {output_dir}/metaphlan_analysis/diversity")
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d alpha -m simpson -o {output_dir}/metaphlan_analysis/diversity")
            shell("{Rscript} OUTPOST/outpost_metaphlan_calculate_diversity.R \
                   -f {output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt \
                   -d alpha -m gini -o {output_dir}/metaphlan_analysis/diversity")
            shell("touch {output_dir}/log/metaphlan_diverisy.done")
    
    
    rule metaphlan_init:
        input:
            "{output_dir}/log/cat_reads.done"
        output:
            "{output_dir}/log/metaphlan_init.done"
        log:
            "{output_dir}/log/metaphlan_init.log"
        benchmark:
            "{output_dir}/benchmark/metaphlan_init.benchmark"
        threads: cores
        run:
            os.makedirs(rf"{output_dir}/metaphlan_analysis/bowtie2", exist_ok=True)
            os.makedirs(rf"{output_dir}/metaphlan_analysis/taxonomy", exist_ok=True)
            os.makedirs(rf"{output_dir}/metaphlan_analysis/figs", exist_ok=True)
            os.makedirs(rf"{output_dir}/metaphlan_analysis/diversity", exist_ok=True)
            os.makedirs(rf"{output_dir}/metaphlan_analysis/krona", exist_ok=True)
            metaphlan_output_list = []
            for i, sample in enumerate(sample_list):
                fq = rf"{output_dir}/data/{sample}_nonhostvirus.fq"
                output = rf"{output_dir}/metaphlan_analysis/taxonomy/metaphlan_{sample}_taxa.txt"
                metaphlan_output_list.append(output)
                bowtie_tmp = rf"{output_dir}/metaphlan_analysis/bowtie2/metaphlan_{sample}.bowtie2.bz2"
                if os.path.exists(bowtie_tmp):
                    shell("{metaphlan} {bowtie_tmp} --offline  --nproc {threads} \
                           --input_type bowtie2out -o {output}")
                else:
                    shell("{metaphlan} {fq} --offline --bowtie2out {bowtie_tmp}  --nproc {threads} \
                          --input_type fastq -o {output}")
            # merge
            metaphlan_outputs = " ".join(metaphlan_output_list)
            output_merged = rf"{output_dir}/metaphlan_analysis/taxonomy/metaphlan_merge_taxa.txt"
            shell("{merge_metaphlan_tables} {metaphlan_outputs} > {output_merged}")
            # reformat to tsv
            # TODO
            # reformat metaphlan taxa annotation to tsv to match OUTPOST taxonomy analysis
            if False:
                sample_list_str = ','.join(sample_list)
                shell("{python3} OUTPOST/taxa_annatation2tsv.py {output_merged} {sample_list_str} \
                      {output_dir}/metaphlan_analysis/taxonomy")
                  
            shell("touch {output_dir}/log/metaphlan_init.done")
    
    
    # %% humann
    rule lefse_humann:
        input:
            "{output_dir}/log/rel_abun2lefse_humann.done"
        output:
            "{output_dir}/log/lefse_humann.done"
        log:
            "{output_dir}/log/lefse_humann.log"
        benchmark:
            "{output_dir}/benchmark/rel_abun2lefse_humann.benchmark"
        threads: cores
        run:
            wkdir = f"{output_dir}/LDA_analysis/figs_humann/"
            os.makedirs(wkdir, exist_ok=True)
            # unequal humann to lefse barplot
            for group1,group2 in group_pair_list:
                dp_list = [f"{output_dir}/LDA_analysis/humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.lefse.tsv"
                           for database in databases]
                format_list = [dp.replace(".tsv",".format") for dp in dp_list]
                res_list = [dp.replace(".tsv",".res") for dp in dp_list]
                output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]
    
                for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                    try:
                        shell("nohup {lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ; "
                              "{lefse_run} {format_}  {res} > {log} 2>&1 ; "
                              "{python3} OUTPOST/barplots_from_lefse_res_file.py  {res} {output} {LDA_cutoff} > {log} 2>&1 & ")
                    except:
                        print(f"OUTPOST: lefse_humann: cannot process {dp}. skip.")
    
            # equal humann to lefse barplot
            for group1,group2 in group_pair_list[:-2]:
                dp_list = [f"{output_dir}/LDA_analysis/humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.lefse.tsv"
                           for database in databases]
                format_list = [dp.replace(".tsv",".format") for dp in dp_list]
                res_list = [dp.replace(".tsv",".res") for dp in dp_list]
                output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]
    
                for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                    try:
                        shell("nohup {lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ; "
                              "{lefse_run} {format_}  {res} > {log} 2>&1 ; "
                              "{python3} OUTPOST/barplots_from_lefse_res_file.py  {res} {output} {LDA_cutoff} > {log} 2>&1 & ")
                    except:
                        print(f"OUTPOST: lefse_humann: cannot process {dp}. skip.")
    
            # equal humann to lefse barplot
            for group1,group2 in group_pair_list[-2:]:
                dp_list = [f"{output_dir}/LDA_analysis/humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.lefse.tsv"
                           for database in databases]
                format_list = [dp.replace(".tsv",".format") for dp in dp_list]
                res_list = [dp.replace(".tsv",".res") for dp in dp_list]
                output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]
    
                for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                    try:
                        shell("{lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ")
                        shell("{lefse_run} {format_}  {res} > {log} 2>&1 ")
                        shell("{python3} OUTPOST/barplots_from_lefse_res_file.py  {res} {output} {LDA_cutoff} > {log} 2>&1 ")
                    except:
                        print(f"OUTPOST: lefse_humann: cannot process {dp}. skip.")
    
            shell("touch {output_dir}/log/lefse_humann.done")
    
    
    rule rel_abun2lefse_humann:
        input:
            "{output_dir}/log/humann_utest.done"
        output:
            "{output_dir}/log/rel_abun2lefse_humann.done"
        log:
            "{output_dir}/log/rel_abun2lefse_humann.log"
        benchmark:
            "{output_dir}/benchmark/rel_abun2lefse_humann.benchmark"
        run:
            # unequal humann to lefse table
            for group1,group2 in group_pair_list:
                wkdir = f"{output_dir}/LDA_analysis/humann_{group1}_vs_{group2}/"
                os.makedirs(wkdir, exist_ok=True)
                dp_list = ','.join([f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.csv"
                                    for database in databases])
                output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                    for dp in dp_list.split(',')])
                class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
                shell("{python3} OUTPOST/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")
            # equal humann to lefse table
            for group1,group2 in group_pair_list:
                wkdir = f"{output_dir}/LDA_analysis/humann_{group1}_vs_{group2}/"
                os.makedirs(wkdir, exist_ok=True)
                dp_list = ','.join([f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.csv"
                                    for database in databases])
                output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                    for dp in dp_list.split(',')])
                class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
                shell("{python3} OUTPOST/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")
    
            shell("touch {output_dir}/log/rel_abun2lefse_humann.done")
    
    
    rule heatmap_humann:
        input:
            "{output_dir}/log/scale_humann_rel_abun_table.done"
        output:
            "{output_dir}/log/heatmap_humann.done"
        log:
            "{output_dir}/log/heatmap_humann.log"
        benchmark:
            "{output_dir}/benchmark/heatmap_humann.benchmark"
        run:
            wkdir = f"{output_dir}/function_analysis/figs/"
            os.makedirs(wkdir, exist_ok=True)
    
            with open(f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_cpm_{databases[0]}_unstratified.named.tsv",'r') as r:
                line = r.readline()
                humann_sample_list = [x.replace('_Abundance-RPKs','') for x in line.strip().split('\t')[1:]]
    
            # top humann heatmap
            data_type = 'humann'
            for database in databases:
                dp = f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.fillmin.scaled.csv"
                output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                index_dp = f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.fillmin.scaled.index"
                with open(index_dp,'w') as w:
                    for sample in humann_sample_list:
                        w.write(sample + '\n')
                try:
                    shell("{Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                except: pass
            # top unequal humann heatmap
            for group1,group2 in group_pair_list:
                for database in databases:
                    dp = f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.fillmin.scaled.csv"
                    output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                    index_dp = f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.fillmin.scaled.index"
                    with open(index_dp,'w') as w:
                        for sample in comparison_dict[group1]:
                            w.write(f"{sample}_{group1}\n")
                        for sample in comparison_dict[group2]:
                            w.write(f"{sample}_{group2}\n")
                    try:
                        shell("{Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                    except: pass
            # top equal humann heatmap
            for group1,group2 in group_pair_list:
                for database in databases:
                    dp = f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.fillmin.scaled.csv"
                    output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                    index_dp = f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.fillmin.scaled.index"
                    with open(index_dp,'w') as w:
                        for sample in comparison_dict[group1]:
                            w.write(f"{sample}_{group1}\n")
                        for sample in comparison_dict[group2]:
                            w.write(f"{sample}_{group2}\n")
                    try:
                        shell("{Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                    except: pass
            shell("touch {output_dir}/log/heatmap_humann.done")
            
            
    rule scale_humann_rel_abun_table:
        input:
            "{output_dir}/log/extract_top_humann.done"
        output:
            "{output_dir}/log/scale_humann_rel_abun_table.done"
        log:
            "{output_dir}/log/scale_humann_rel_abun_table.log"
        benchmark:
            "{output_dir}/benchmark/scale_humann_rel_abun_table.benchmark"
        run:
            dp_list = [f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.csv"
                           for database in databases]
            output_list = [f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.fillmin.scaled.csv"
                           for database in databases]
    
            for group1,group2 in group_pair_list:
                dp_list_equal = [f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.csv"
                                 for database in databases]
                output_list_equal = [f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.fillmin.scaled.csv"
                                 for database in databases]
    
                dp_list_unequal = [f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.csv"
                                 for database in databases]
                output_list_unequal = [f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.fillmin.scaled.csv"
                                 for database in databases]
    
                dp_list += dp_list_equal
                dp_list += dp_list_unequal
                output_list += output_list_equal
                output_list += output_list_unequal
            dp_list = ','.join(dp_list)
            output_list = ','.join(output_list)
            try:
                shell("{python3} OUTPOST/scale_rel_abun_table.py {dp_list} {output_list}")
            except:
                # in case augment == too long
                print("OUTPOST: rule scale_rel_abun_table Argument list too long. use chunks.")
                dp_list_chunks = [dp_list.split(',')[i:i + 100] for i in range(0, len(dp_list.split(',')), 100)]
                output_list_chunks = [output_list.split(',')[i:i + 100] for i in range(0, len(output_list.split(',')), 100)]
                for dp_list_chunk,output_list_chunk in zip(dp_list_chunks, output_list_chunks):
                    dp_list_chunk = ','.join(dp_list_chunk)
                    output_list_chunk = ','.join(output_list_chunk)
                    shell("{python3} OUTPOST/scale_rel_abun_table.py {dp_list_chunk} {output_list_chunk}")
            shell("touch {output_dir}/log/scale_humann_rel_abun_table.done")
    
    
    rule extract_top_humann:
        input:
            "{output_dir}/log/humann_utest.done"
        output:
            "{output_dir}/log/extract_top_humann.done"
        log:
            "{output_dir}/log/extract_top_humann.log"
        benchmark:
            "{output_dir}/benchmark/extract_top_humann.benchmark"
        threads: cores
        run:
            # top abundant
            dp_list = ','.join([f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.csv"
                           for database in databases])
            output_name_list = ','.join([f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.top{top}.csv"
                           for database in databases])
            error_log = os.getcwd() + f"/{output_dir}/log/extract_top_humann.error"
            shell("{python3} OUTPOST/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
                  " {top} {error_log} > {log} 2>&1 ")
            result = wait_until_file_exists(output_name_list.split(','),error_log)
            assert result == True, "OUTPOST: extract_top_humann has errors. exit."
    
            # top unequal abundant
            output_list = []
            for group1,group2 in group_pair_list:
                os.makedirs(f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}", exist_ok=True)
                dp_list = ','.join([f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.csv"
                           for database in databases])
                output_name_list = ','.join([f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.csv"
                        for database in databases])
                error_log = os.getcwd() + f"/{output_dir}/log/extract_top_humann_top_unequal_abundant.error"
                output_list.append(f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.unequal.top{top}.csv")
                shell("nohup {python3} OUTPOST/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
                  " {top} {error_log} > {log}_{group1}_vs_{group2}_unequal 2>&1 &")
                time.sleep(3)
            result = wait_until_file_exists(output_list,error_log)
            assert result == True, "OUTPOST: extract_top_humann_top_unequal_abundant has errors. exit."
    
           # top equal abundant
            output_list = []
            for group1,group2 in group_pair_list:
                os.makedirs(f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}", exist_ok=True)
                dp_list = ','.join([f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.csv"
                           for database in databases])
                output_name_list = ','.join([f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.csv"
                        for database in databases])
                error_log = os.getcwd() + f"/{output_dir}/log/extract_top_humann_top_equal_abundant.error"
                output_list.append(f"{output_dir}/function_analysis/humann3/top_humann_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.{group1}_vs_{group2}.rel_abun.equal.top{top}.csv")
                shell("nohup {python3} OUTPOST/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
                  " {top} {error_log} > {log}_{group1}_vs_{group2}_equal 2>&1 &")
                time.sleep(3)
            result = wait_until_file_exists(output_list,error_log)
            assert result == True, "OUTPOST: extract_top_humann_top_equal_abundant has errors. exit."
    
            shell("touch {output_dir}/log/extract_top_humann.done")
    
    
    rule humann_utest:
        input:
            "{output_dir}/log/humann2rel_abun.done"
        output:
            "{output_dir}/log/humann_utest.done"
        log:
            "{output_dir}/log/humann_utest.log"
        benchmark:
            "{output_dir}/benchmark/humann_utest.benchmark"
        threads: cores
        run:
            output_list = []
            dp_list = ','.join([f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.csv"
                       for database in databases])
            with open(f"{output_dir}/function_analysis/humann3/output/allSamples_genefamilies_uniref90names_cpm_{databases[0]}_unstratified.named.tsv",'r') as r:
                line = r.readline()
                humann_sample_list = [x.replace('_Abundance-RPKs','') for x in line.strip().split('\t')[1:]]
            for group1,group2 in group_pair_list:
                os.makedirs(f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}", exist_ok=True)
    
                group1_index = ','.join([str(humann_sample_list.index(sample)) for sample in comparison_dict[group1]])
                group2_index = ','.join([str(humann_sample_list.index(sample)) for sample in comparison_dict[group2]])
                prefix_list = ','.join([f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.{group1}_vs_{group2}"
                           for database in databases])
                output = f"{output_dir}/function_analysis/humann3/utest_{group1}_vs_{group2}/allSamples_genefamilies_uniref90names_relab_{databases[-1]}_unstratified.named.rel_abun_format.{group1}_vs_{group2}.ave_change.unequal.csv"
                output_list.append(output)
                error_log = os.getcwd() + f"/{output_dir}/log/humann_utest.error"
                try:
                    os.remove(error_log) # clean former residual error log
                except:
                    pass
                shell("nohup {python3} OUTPOST/rel_abun_utest.py {dp_list} {group1_index} {group1} "
                      " {group2_index} {group2} {prefix_list} {paired} {two_sided} {error_log} > {log}_{group1}_vs_{group2} 2>&1 &")
                time.sleep(3)
            result = wait_until_file_exists(output_list,error_log)
            assert result == True, "OUTPOST: rel_abun_utest has errors. exit."
            shell("touch {output_dir}/log/humann_utest.done")
    
    rule humann2rel_abun:
        input:
            "{output_dir}/log/humann_output.done"
        output:
            "{output_dir}/log/humann2rel_abun.done"
        log:
            "{output_dir}/log/humann2rel_abun.log"
        benchmark:
            "{output_dir}/benchmark/humann2rel_abun.benchmark"
        run:
            for group in list(comparison_dict.keys()) + ['allSamples']:
                dp_list = ','.join([f"{output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}_unstratified.named.tsv"
                           for database in databases])
                output_list = ','.join([f"{output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}_unstratified.named.rel_abun_format.csv"
                           for database in databases])
                shell("{python3} OUTPOST/humann2rel_abun.py {dp_list} {output_list} > {log} 2>&1 ")
            shell("touch {output_dir}/log/humann2rel_abun.done")
    
    rule humann_output:
        input:
            "{output_dir}/log/humann_group.done"
        output:
            "{output_dir}/log/humann_output.done"
        log:
            "{output_dir}/log/humann_output.log"
        benchmark:
            "{output_dir}/benchmark/humann_output.benchmark"
        run:
            os.makedirs(f"{output_dir}/function_analysis/humann3/output/", exist_ok=True)
            for group in list(comparison_dict.keys()) + ['allSamples']:
                # join tables for each group
                shell("{humann_join_tables} -i {output_dir}/function_analysis/humann3/{group}/ "
                      " -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                      "--file_name _genefamilies_relab.tsv  > {log} 2>&1 ")
                shell("{humann_join_tables} -i {output_dir}/function_analysis/humann3/{group}/ "
                      " -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                      "--file_name _genefamilies_cpm.tsv  > {log} 2>&1 ")
                shell("{humann_join_tables} -i {output_dir}/function_analysis/humann3/{group}/ "
                      " -o {output_dir}/function_analysis/humann3/output/{group}_pathabundance_relab.tsv "
                      "--file_name _pathabundance_relab.tsv  > {log} 2>&1 ")
                shell("{humann_join_tables} -i {output_dir}/function_analysis/humann3/{group}/ "
                      " -o {output_dir}/function_analysis/humann3/output/{group}_pathabundance_cpm.tsv "
                      "--file_name _pathabundance_cpm.tsv  > {log} 2>&1 ")
                shell("{humann_join_tables} -i {output_dir}/function_analysis/humann3/{group}/ "
                      " -o {output_dir}/function_analysis/humann3/output/{group}_pathcoverage_relab.tsv "
                      "--file_name _pathcoverage_relab.tsv  > {log} 2>&1 ")
                shell("{humann_join_tables} -i {output_dir}/function_analysis/humann3/{group}/ "
                      " -o {output_dir}/function_analysis/humann3/output/{group}_pathcoverage_cpm.tsv "
                      "--file_name _pathcoverage_cpm.tsv  > {log} 2>&1 ")
    
                # regroup tables for each group
                for database in databases:
                    shell("{humann_regroup_table} --input {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm.tsv "
                          " --output {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")
                    shell("{humann_regroup_table} --input {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab.tsv "
                          " --output {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")
    
                # stratify tables for each group
                for database in databases:
                    shell("{humann_split_stratified_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_{database}.tsv "
                          "-o {output_dir}/function_analysis/humann3/output/ > {log} 2>&1 ")
                    shell("{humann_split_stratified_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_{database}.tsv "
                          "-o {output_dir}/function_analysis/humann3/output/ > {log} 2>&1 ")
    
                # rename tables for each group
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.tsv "
                      " --names kegg-orthology -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_ko_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.tsv "
                      " --names eggnog -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.tsv "
                      " --names ec -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_level4ec_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.tsv "
                      " --names pfam -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_pfam_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.tsv "
                      " --names metacyc-rxn -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_cpm_rxn_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko_unstratified.tsv "
                      " --names kegg-orthology -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_ko_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog_unstratified.tsv "
                      " --names eggnog -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_eggnog_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec_unstratified.tsv "
                      " --names ec -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_level4ec_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam_unstratified.tsv "
                      " --names pfam -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_pfam_unstratified.named.tsv > {log} 2>&1 ")
                shell("{humann_rename_table} -i {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.tsv "
                      " --names metacyc-rxn -o {output_dir}/function_analysis/humann3/output/{group}_genefamilies_uniref90names_relab_rxn_unstratified.named.tsv > {log} 2>&1 ")
            shell("touch {output_dir}/log/humann_output.done")
    
    rule humann_group:
        input:
            "{output_dir}/log/humann_annotate.done",
        output:
            "{output_dir}/log/humann_group.done"
        log:
            "{output_dir}/log/humann_group.log"
        benchmark:
            "{output_dir}/benchmark/humann_group.benchmark"
        run:
            for group,samples in comparison_dict.items():
                os.makedirs(f"{output_dir}/function_analysis/humann3/{group}", exist_ok=True)
                for sample in samples:
                    os.system(f"cp {output_dir}/function_analysis/humann3/ori_results/{sample}*tsv  {output_dir}/function_analysis/humann3/{group} ")
            # allSamples
            os.makedirs(f"{output_dir}/function_analysis/humann3/allSamples", exist_ok=True)
            for sample in sample_list:
                os.system(f"cp {output_dir}/function_analysis/humann3/ori_results/{sample}*tsv  {output_dir}/function_analysis/humann3/allSamples ")
            shell("touch {output_dir}/log/humann_group.done")
    
    rule humann_annotate:
        input:
            "{output_dir}/log/rename_humann_ori_output.done"
        output:
            "{output_dir}/log/humann_annotate.done"
        log:
            "{output_dir}/log/humann_annotate.log"
        benchmark:
            "{output_dir}/benchmark/humann_annotate.benchmark"
        run:
            for sample in sample_list:
                # normalize
                shell("{humann_renorm_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies.tsv "
                      " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                      " --units relab > {log} 2>&1 ")
                shell("{humann_renorm_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies.tsv "
                      " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                      " --units cpm > {log} 2>&1 ")
                shell("{humann_renorm_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_pathabundance.tsv "
                      " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_pathabundance_relab.tsv "
                      " --units relab > {log} 2>&1 ")
                shell("{humann_renorm_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_pathabundance.tsv "
                      " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_pathabundance_cpm.tsv "
                      " --units cpm > {log} 2>&1 ")
                shell("{humann_renorm_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_pathcoverage.tsv "
                      " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_pathcoverage_relab.tsv "
                      " --units relab > {log} 2>&1 ")
                shell("{humann_renorm_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_pathcoverage.tsv "
                      " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_pathcoverage_cpm.tsv "
                      " --units cpm > {log} 2>&1 ")
                # regroup
                for database in databases:
                    shell("{humann_regroup_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies_cpm.tsv "
                          " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies_cpm_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")
                    shell("{humann_regroup_table} --input {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies_relab.tsv "
                          " --output {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies_relab_{database}.tsv --groups uniref90_{database} > {log} 2>&1 ")
    
            shell("touch {output_dir}/log/humann_annotate.done")
    
    rule rename_humann_ori_output:
        input:
            "{output_dir}/log/humann_init.done"
        output:
            "{output_dir}/log/rename_humann_ori_output.done"
        log:
            "{output_dir}/log/rename_humann_ori_output.log"
        benchmark:
            "{output_dir}/benchmark/humann_init.benchmark"
        run:
            for sample, fq_humann in zip(sample_list, fq_humann_list):
                shell("cp {output_dir}/function_analysis/humann3/ori_results/{fq_humann}_genefamilies.tsv "
                      " {output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies.tsv")
                shell("cp {output_dir}/function_analysis/humann3/ori_results/{fq_humann}_pathabundance.tsv "
                      " {output_dir}/function_analysis/humann3/ori_results/{sample}_pathabundance.tsv")
                shell("cp {output_dir}/function_analysis/humann3/ori_results/{fq_humann}_pathcoverage.tsv "
                      " {output_dir}/function_analysis/humann3/ori_results/{sample}_pathcoverage.tsv")
    
                file1 = f"{output_dir}/function_analysis/humann3/ori_results/{fq_humann}_genefamilies.tsv"
                file2 = f"{output_dir}/function_analysis/humann3/ori_results/{sample}_genefamilies.tsv"
                with open(file1,'r') as r, open(file2,'w') as w:
                    lines = r.readlines()
                    lines[0] = "# Gene Family\t" + sample + "_Abundance-RPKs\n"
                    for line in lines:
                        w.write(line)
    
                file1 = f"{output_dir}/function_analysis/humann3/ori_results/{fq_humann}_pathabundance.tsv"
                file2 = f"{output_dir}/function_analysis/humann3/ori_results/{sample}_pathabundance.tsv"
                with open(file1,'r') as r, open(file2,'w') as w:
                    lines = r.readlines()
                    lines[0] = "# Pathway\t" + sample + "_Abundance\n"
                    for line in lines:
                        w.write(line)
    
                file1 = f"{output_dir}/function_analysis/humann3/ori_results/{fq_humann}_pathcoverage.tsv"
                file2 = f"{output_dir}/function_analysis/humann3/ori_results/{sample}_pathcoverage.tsv"
                with open(file1,'r') as r, open(file2,'w') as w:
                    lines = r.readlines()
                    lines[0] = "# Pathway\t" + sample + "_Coverage\n"
                    for line in lines:
                        w.write(line)
            shell("touch {output_dir}/log/rename_humann_ori_output.done")
    
    if not skip_humann_init:
        rule humann_init:
            input:
                [f"{output_dir}/log/reads_QC_PE.done"] * bool(len(r1_fq_list) > 0) +\
                [f"{output_dir}/log/reads_QC_SE.done"] * bool(len(se_fq_list) > 0)
            output:
                "{output_dir}/log/humann_init.done"
            log:
                "{output_dir}/log/humann_init.log"
            benchmark:
                "{output_dir}/benchmark/humann_init.benchmark"
            threads: cores
            run:
                def does_humann_need_to_run(fq, wkdir):
                    # get the humann renamed file name, taken from humann source code
                    fq_humann = os.path.basename(fq)
                    # Remove gzip extension if present
                    if re.search('.gz$',fq_humann):
                        fq_humann='.'.join(fq_humann.split('.')[:-1])
                    # Remove input file extension if present
                    if '.' in fq_humann:
                        fq_humann='.'.join(fq_humann.split('.')[:-1])
                        
                    genefamilies_tsv = os.path.join(wkdir,f"{fq_humann}_genefamilies.tsv")
                    pathabundance_tsv = os.path.join(wkdir,f"{fq_humann}_pathabundance.tsv")
                    pathcoverage_tsv = os.path.join(wkdir,f"{fq_humann}_pathcoverage.tsv")
                    need_to_run = False
                    if not os.path.exists(genefamilies_tsv) and \
                        not os.path.exists(pathabundance_tsv) and \
                        not os.path.exists(pathcoverage_tsv):
                        need_to_run = True
                    return need_to_run

                os.makedirs(f"{output_dir}/function_analysis/humann3/ori_results/", exist_ok=True)
                fq_list_batch = [fq_list[i:i+int(process_batch_size)] for i in range(0, len(fq_list), int(process_batch_size))]
                wkdir = rf"{output_dir}/function_analysis/humann3/ori_results/"
                for fq_batch in fq_list_batch:
                    if len(fq_batch) == 1:
                        fq = fq_batch[0]
                        need_to_run = does_humann_need_to_run(fq, wkdir)
                        if need_to_run:
                            shell("{humann} --resume --input {fq}  --output {wkdir} "
                            " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
                            " --threads {threads} --memory-use {memory_use} > {log} 2>&1 ")
                    elif len(fq_batch) > 1:
                        for fq in fq_batch[:-1]:
                            need_to_run = does_humann_need_to_run(fq, wkdir)
                            if need_to_run:
                                shell("nohup {humann} --resume --input {fq}  --output {wkdir} "
                                " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
                                " --threads {threads} --memory-use {memory_use} > {log} 2>&1 &")
                        fq = fq_batch[-1]
                        shell("{humann} --resume --input {fq}  --output {wkdir} "
                        " --search-mode uniref90 --diamond-options '--block-size 10 --fast' "
                        " --threads {threads} --memory-use {memory_use} > {log} 2>&1")
                # clean temp folder
                if clean_unnecessary:
                    try:
                        shell("rm -r {wkdir}*_humann_temp")
                    except Exception as e:
                        pass
                shell("touch {output_dir}/log/humann_init.done")
    else:
        rule skip_humann_init:
            input:
                "{output_dir}/log/cat_reads.done"
            output:
                "{output_dir}/log/humann_init.done"
            log:
                "{output_dir}/log/humann_init.log"
            benchmark:
                "{output_dir}/benchmark/humann_init.benchmark"
            run:
                humann_genefamilies = [f"{output_dir}/function_analysis/humann3/ori_results/{fq_humann}_genefamilies.tsv" for fq_humann in fq_humann_list]
                humann_pathabundance = [f"{output_dir}/function_analysis/humann3/ori_results/{fq_humann}_pathabundance.tsv" for fq_humann in fq_humann_list]
                humann_pathcoverage = [f"{output_dir}/function_analysis/humann3/ori_results/{fq_humann}_pathcoverage.tsv" for fq_humann in fq_humann_list]
                missing_humann = [file for file in humann_genefamilies + humann_pathabundance + humann_pathcoverage if not os.path.exists(file)]
                assert len(missing_humann) == 0,\
                f"OUTPOST: detected missing genefamilies/pathabundance/pathcoverage humann files under {output_dir}/function_analysis/humann3/ori_results/. cannot skip humann_init. exit. expecting {missing_humann}"
                shell("mkdir -p {output_dir}/log/")
                shell("touch {output_dir}/log/humann_init.done")
    
    
    
#%% assembly
if len(assembly_dir_list) > 0:
    if len(r1_fq_list + se_fq_list) > 0:
        #%% salmon
        rule salmon_gene_quantify:
            input:
                ["{output_dir}/log/prodigal_gene_prediction.done"] +\
                ["{output_dir}/log/metagenemark.done"] +\
                [f"{output_dir}/log/reads_QC_PE.done"] * bool(len(r1_fq_list) > 0) +\
                [f"{output_dir}/log/reads_QC_SE.done"] * bool(len(se_fq_list) > 0)
            output:
                "{output_dir}/log/salmon_gene_quantify.done"
            log:
                "{output_dir}/log/salmon_gene_quantify.log"
            benchmark:
                "{output_dir}/benchmark/salmon_gene_quantify.benchmark"
            threads: cores
            run:
                salmon_prodigal_index_dir = f"{output_dir}/assembly_analysis/quantify/index/prodigal"
                salmon_metagenemark_index_dir = f"{output_dir}/assembly_analysis/quantify/index/metagenemark"
                os.makedirs(salmon_prodigal_index_dir, exist_ok=True)
                os.makedirs(salmon_metagenemark_index_dir, exist_ok=True)
                
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        prodigal_dir = f"{output_dir}/assembly_analysis/annotation/prodigal"
                        shell("{salmon} index -t  {prodigal_dir}/{assembly}_prodigal_gene.fa \
                              --thread {threads} --index {salmon_prodigal_index_dir} ")
                        metagenemark_dir = f"{output_dir}/assembly_analysis/annotation/metagenemark"
                        shell("{salmon} index -t  {metagenemark_dir}/{assembly}_metagenemark_gene.nonrd.fa \
                              --thread {threads} --index {salmon_metagenemark_index_dir} ")

                        salmon_prodigal_quantify_dir = f"{output_dir}/assembly_analysis/quantify/quantify/prodigal_{assembly}"
                        salmon_metagenemark_quantify_dir = f"{output_dir}/assembly_analysis/quantify/quantify/metagenemark_{assembly}"
                        os.makedirs(salmon_prodigal_quantify_dir, exist_ok=True)
                        os.makedirs(salmon_metagenemark_quantify_dir, exist_ok=True)

                        for i,sample in enumerate(sample_list):
                            shell("{salmon} quant -i {salmon_prodigal_index_dir} -l A --thread {threads} --meta \
                                  -o {salmon_prodigal_quantify_dir} -r {output_dir}/data/{sample}_nonhostvirus.fq")
                            shell("{salmon} quant -i {salmon_metagenemark_index_dir} -l A --thread {threads} --meta \
                                  -o {salmon_metagenemark_quantify_dir} -r {output_dir}/data/{sample}_nonhostvirus.fq")
                        
                shell("touch {output_dir}/log/salmon_gene_quantify.done")
    

    #%% eggmapper
    rule eggmapper_gene_annotation:
        input:
            "{output_dir}/log/prodigal_gene_prediction.done",
            "{output_dir}/log/metagenemark.done"
        output:
            "{output_dir}/log/eggmapper_gene_annotation.done"
        log:
            "{output_dir}/log/eggmapper_gene_annotation.log"
        benchmark:
            "{output_dir}/benchmark/eggmapper_gene_annotation.benchmark"
        threads: cores
        run:
            eggmapper_annotate_prodigal_dir = f"{output_dir}/assembly_analysis/annotation/eggmapper/prodigal"
            eggmapper_annotate_metagenemark_dir = f"{output_dir}/assembly_analysis/annotation/eggmapper/metagenemark"
            os.makedirs(eggmapper_annotate_prodigal_dir, exist_ok=True)
            os.makedirs(eggmapper_annotate_metagenemark_dir, exist_ok=True)
            for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                if assembly != '':
                    prodigal_dir = f"{output_dir}/assembly_analysis/annotation/prodigal"
                    resume = " --resume" * bool(os.path.exists(os.path.join(eggmapper_annotate_prodigal_dir,f"{assembly}.emapper.hits")))
                    shell("{emapper} --data_dir {emapper_db} \
                          -i {prodigal_dir}/{assembly}_prodigal_prot.fa --cpu {threads} -m diamond \
                          --output_dir {eggmapper_annotate_prodigal_dir} --output {assembly} {resume}")
                    
                    metagenemark_dir = f"{output_dir}/assembly_analysis/annotation/metagenemark"
                    resume = " --resume" * bool(os.path.exists(os.path.join(eggmapper_annotate_metagenemark_dir,f"{assembly}.emapper.hits")))
                    shell("{emapper} --data_dir {emapper_db} \
                          -i {metagenemark_dir}/{assembly}_metagenemark_prot.fa --cpu {threads} -m diamond \
                          --output_dir {eggmapper_annotate_metagenemark_dir}  --output {assembly} {resume}")
                          
            shell("touch {output_dir}/log/eggmapper_gene_annotation.done")

    
    #%% rgi
    rule rgi_gene_annotation:
        input:
            "{output_dir}/log/prodigal_gene_prediction.done",
            "{output_dir}/log/metagenemark.done"
        output:
            "{output_dir}/log/rgi_gene_annotation.done"
        log:
            "{output_dir}/log/rgi_gene_annotation.log"
        benchmark:
            "{output_dir}/benchmark/rgi_gene_annotation.benchmark"
        threads: cores
        run:
            rgi_annotate_prodigal_prefix = f"{output_dir}/assembly_analysis/annotation/rgi/prodigal"
            rgi_annotate_metagenemark_prefix = f"{output_dir}/assembly_analysis/annotation/rgi/metagenemark"
            for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                if assembly != '':
                    prodigal_dir = f"{output_dir}/assembly_analysis/annotation/prodigal"
                    shell("{rgi} main -i {prodigal_dir}/{assembly}_prodigal_gene.fa -n {threads} \
                          --include_loose  --local  --clean -o {rgi_annotate_prodigal_prefix} ")
                          
                    metagenemark_dir = f"{output_dir}/assembly_analysis/annotation/metagenemark"
                    shell("{rgi} main -i {metagenemark_dir}/{assembly}_metagenemark_gene.nonrd.fa -n {threads} \
                          --include_loose  --local  --clean -o {rgi_annotate_metagenemark_prefix} ")
            shell("touch {output_dir}/log/rgi_gene_annotation.done")
    
    
    if (len(r1_fq_list + se_fq_list) > 0) and assemble_contigs:
        metagenemark_input = f"{output_dir}/log/quast_contigs.done"
        prodigal_gene_prediction_input = "{output_dir}/log/quast_contigs.done"
        genome_split_input = "{output_dir}/log/quast_contigs.done"
    else:
        metagenemark_input = assembly_dir_list
        prodigal_gene_prediction_input = assembly_dir_list
        genome_split_input = assembly_dir_list
    #%% gtdbtk
    rule genome_split:
        input:
            genome_split_input
        output:
            "{output_dir}/log/genome_split.done"
        log:
            "{output_dir}/log/genome_split.log"
        benchmark:
            "{output_dir}/benchmark/genome_split.benchmark"
        threads: 1
        run:
            for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                if assembly != '':
                    assembly_dirname = os.path.dirname(assembly_dir)
                    genome_split_dir = os.path.join(assembly_dirname, f"{assembly}_split")
                    shell("{seqkit} split -i --by-id {assembly_dir} --out-dir {genome_split_dir}")
            shell("touch {output_dir}/log/genome_split.done")
            
            
    rule gtdbtk_genome_annotation:
        input:
            "{output_dir}/log/genome_split.done"
        output:
            "{output_dir}/log/gtdbtk_genome_annotation.done"
        log:
            "{output_dir}/log/gtdbtk_genome_annotation.log"
        benchmark:
            "{output_dir}/benchmark/gtdbtk_genome_annotation.benchmark"
        threads: cores
        run:
            gtdbtk_dir = f"{output_dir}/assembly_analysis/annotation/gtdbtk"
            os.makedirs(gtdbtk_dir, exist_ok=True)
            extenstion_list, assembly_dirname_list = [], []
            for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                if assembly != '':
                    extenstion = assembly_dir.split('.')[-1]
                    # gtdbtk 要输入一个文件夹，然后处理文件夹里的所有genome
                    assembly_dirname = os.path.dirname(assembly_dir)
                    extenstion_list.append(extenstion)
                    assembly_dirname_list.append(assembly_dirname)
            for assembly, extenstion,assembly_dirname in set(zip(assembly_list, extenstion_list,assembly_dirname_list)):
                if assembly != '':
                    genome_split_dir = os.path.join(assembly_dirname, f"{assembly}_split")
                    shell("gtdbtk classify_wf --genome_dir {genome_split_dir} --out_dir {gtdbtk_dir} \
                      --extension {extenstion} --skip_ani_screen --prefix {assembly}  --cpus {threads}")
            shell("touch {output_dir}/log/gtdbtk_genome_annotation.done")
    
    
    #%% prodigal
    rule prodigal_gene_prediction:
        input:
            prodigal_gene_prediction_input
        output:
            "{output_dir}/log/prodigal_gene_prediction.done"
        log:
            "{output_dir}/log/prodigal_gene_prediction.log"
        benchmark:
            "{output_dir}/benchmark/prodigal_gene_prediction.benchmark"
        threads: 1
        run:
            prodigal_dir = f"{output_dir}/assembly_analysis/annotation/prodigal"
            os.makedirs(prodigal_dir, exist_ok=True)
            for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                if assembly != '':
                    # predict genes
                    shell("{prodigal} -i {assembly_dir} -d {prodigal_dir}/{assembly}_prodigal_gene.fa -f gff -o {prodigal_dir}/{assembly}_prodigal_gene.gff -p meta -s {prodigal_dir}/{assembly}_prodigal_gene_scores.txt -q")
                    shell("grep 'partial=00' {prodigal_dir}/{assembly}_prodigal_gene.fa | cut -f1 -d ' '| sed 's/>//' > {prodigal_dir}/{assembly}_prodigal_complete_genomes.id")
                    shell("{seqkit} grep -f {prodigal_dir}/{assembly}_prodigal_complete_genomes.id {prodigal_dir}/{assembly}_prodigal_gene.fa > {prodigal_dir}/{assembly}_prodigal_complete_genomes.fa")
                    shell("{seqkit} translate --trim {prodigal_dir}/{assembly}_prodigal_gene.fa > {prodigal_dir}/{assembly}_prodigal_prot.fa ")
            shell("touch {output_dir}/log/prodigal_gene_prediction.done")
            
                    
    # %% metagenemark
    rule metagenemark:
        input:
            metagenemark_input
        output:
            "{output_dir}/log/metagenemark.done"
        log:
            "{output_dir}/log/metagenemark.log"
        benchmark:
            "{output_dir}/benchmark/metagenemark.benchmark"
        threads: cores
        run:
            metagenemark_output_dir = f"{output_dir}/assembly_analysis/annotation/metagenemark"
            os.makedirs(metagenemark_output_dir, exist_ok=True)
            for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                if assembly != '':
                    # predict genes
                    shell("{mgm} -r -p 1 -m {mod_file} -k -K {metagenemark_output_dir}/{assembly}_metagenemark.rbs -f G -o {metagenemark_output_dir}/{assembly}_metagenemark.gtf -a -A {metagenemark_output_dir}/{assembly}_metagenemark_prot.fa -d -D {metagenemark_output_dir}/{assembly}_metagenemark_gene.fa {assembly_dir}")
                    # remove redundancy
                    shell("{cdhit} -i {metagenemark_output_dir}/{assembly}_metagenemark_gene.fa -o {metagenemark_output_dir}/{assembly}_metagenemark_gene.nonrd.fa -c {cdhit_cutoff} -n 10 -M 0 -T {threads}")
            shell("touch {output_dir}/log/metagenemark.done")
                
                
    if (len(r1_fq_list + se_fq_list) > 0 and len(assembly_list) > 0) or (len(bam_dir_list) > 0):
        # %% biomarker
        rule biomarker_summary:
            input:
                "{output_dir}/log/ancom_identification.done",
                "{output_dir}/log/rel_abun_utest.done",
                "{output_dir}/log/virulence_factors_analysis_assembly_based.done",
                "{output_dir}/log/plasmids_analysis_assembly_based.done",
                "{output_dir}/log/antibiotic_genes_analysis_assembly_based.done",
                "{output_dir}/log/lefse_taxa.done"
            output:
                "{output_dir}/log/biomarker_summary.done"
            log:
                "{output_dir}/log/biomarker_summary.log"
            benchmark:
                "{output_dir}/benchmark/biomarker_summary.benchmark"
            threads: 1
            run:
                wkdir = f"{output_dir}/biomarkers_analysis/"
                os.makedirs(wkdir, exist_ok=True)
                
                # init table
                init_tb = pd.read_csv(f'{output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv', sep = '\t', header=None)

                for level in taxa_level:
                    for group1,group2 in group_pair_list:
                        biomarker_score = [f'taxa_abun_top_{biomarker_num}',
                                           f'taxa_abun_signicant',
                                           f'LDA_larger_than_{LDA_cutoff}',
                                           f'virulence_top_{biomarker_num}',
                                           f'plasmid_top_{biomarker_num}',
                                           f'antibiotic_top_{biomarker_num}',
                                           f'ANCOM_identified']
                        taxa_list = list(set(init_tb[taxa_level.index(level)+1]))
                        biomarker_tb = pd.DataFrame(index = taxa_list, columns = biomarker_score)
                        biomarker_tb = biomarker_tb.fillna(0)
                        
                        # taxa_abun_top
                        try:
                            dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.csv"
                            tmp = pd.read_csv(dp, index_col=0)
                            row_sums = tmp.sum(axis=0)
                            sorted_row_sums =  sorted_row_sums = row_sums.sort_values(ascending=False)
                            biomarker_list = sorted_row_sums.index[:biomarker_num]
                            for index in biomarker_list: 
                                try:
                                    biomarker_tb.loc[index, f'taxa_abun_top_{biomarker_num}'] = 1
                                except:
                                    pass
                        except:
                            pass
                        # taxa_abun_signicant
                        dp = f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.u-test." +"paired."*paired + "two_sided."*two_sided + "csv"
                        tmp = pd.read_csv(dp, index_col=0)
                        for index in biomarker_tb.index: 
                            try:
                                if tmp.loc['Wilcoxon',index] == 'unequal':
                                    biomarker_tb.loc[index, 'taxa_abun_signicant'] = 1
                            except:
                                pass
                            
                        # LDA_larger_than
                        dp1 = f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.lefse.res"
                        dp2 = f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.lefse.res"
                        dp1_exists = os.path.exists(dp1)
                        dp2_exists = os.path.exists(dp2)
                        dp1 = dp1 * dp1_exists
                        dp2 = dp2 * dp2_exists
                        dp = f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.lefse.res"
                        try:
                            shell("cat {dp1} {dp2} > {dp}")
                            tmp = pd.read_csv(dp, sep = '\t', index_col=0, header = None)
                            tmp1 = tmp.loc[tmp[2].notnull()]
                            tmp2 = tmp1.loc[tmp1[3] > LDA_cutoff]
                            if tmp2.shape[0] > 0:
                                for index in tmp2.index:
                                    hit = 0
                                    # index = index.strip().split('_')[-1]
                                    for real_index in biomarker_tb.index:
                                        if isinstance(real_index, str):
                                            try:
                                                if real_index.replace('.','').replace('-','').replace(' ','').replace(':','').replace('[','').replace(']','').lower() == index.replace('_','').lower():
                                                    biomarker_tb.loc[real_index, f'LDA_larger_than_{LDA_cutoff}'] = 1
                                                    hit = 1
                                                    break
                                            except:
                                                pass
                                    if hit == 0:
                                        biomarker_tb.loc[index, f'LDA_larger_than_{LDA_cutoff}'] = 1
                        except:
                            pass
                        
                        # virulence_top
                        dp = f"{output_dir}/virulence_factors_analysis/assembly_based/genes_taxa_counts_{group1}_vs_{group2}.tsv"
                        try:
                            tmp = pd.read_csv(dp, sep = '\t')
                            level_index = taxa_level.index(level)
                            tmp1 = tmp[[level,'total']]
                            tmp2 = tmp1.groupby(level, as_index=False).sum()
                            tmp3 = tmp2.sort_values(by='total', ascending=False)
                            tmp3.index = range(tmp3.shape[0])
                            for index in tmp3.loc[:biomarker_num-1, level]:
                                try:
                                    biomarker_tb.loc[index, f'virulence_top_{biomarker_num}'] = 1
                                except:
                                    pass
                        except Exception as e:
                            print("OUTPOST warning: cannot extract virulence information for biomarker analysis.")
                            print(e)
                            

                        # plasmid_top
                        dp = f"{output_dir}/plasmids_analysis/assembly_based/genes_taxa_counts_{group1}_vs_{group2}.tsv"
                        try:
                            tmp = pd.read_csv(dp, sep = '\t')
                            level_index = taxa_level.index(level)
                            tmp1 = tmp[[level,'total']]
                            tmp2 = tmp1.groupby(level, as_index=False).sum()
                            tmp3 = tmp2.sort_values(by='total', ascending=False)
                            tmp3.index = range(tmp3.shape[0])
                            for index in tmp3.loc[:biomarker_num-1, level]:
                                try:
                                    biomarker_tb.loc[index, f'plasmid_top_{biomarker_num}'] = 1
                                except:
                                    pass
                        except Exception as e:
                            print("OUTPOST warning: cannot extract plasmid information for biomarker analysis.")
                            print(e)

                        # antibiotic_top
                        dp = f"{output_dir}/antibiotic_genes_analysis/assembly_based/genes_taxa_counts_{group1}_vs_{group2}.tsv"
                        try:
                            tmp = pd.read_csv(dp, sep = '\t')
                            level_index = taxa_level.index(level)
                            tmp1 = tmp[[level,'total']]
                            tmp2 = tmp1.groupby(level, as_index=False).sum()
                            tmp3 = tmp2.sort_values(by='total', ascending=False)
                            tmp3.index = range(tmp3.shape[0])
                            for index in tmp3.loc[:biomarker_num-1, level]:
                                try:
                                    biomarker_tb.loc[index, f'antibiotic_top_{biomarker_num}'] = 1
                                except:
                                    pass
                        except Exception as e:
                            print("OUTPOST warning: cannot extract antbiotic genes information for biomarker analysis.")
                            print(e)

                        # ANCOM_identified
                        dp = f"{output_dir}/biomarkers_analysis/ANCOM_identification/ancom_biomarker.{group1}_vs_{group2}.at_{level}.tsv"
                        tmp = pd.read_csv(dp, sep = '\t',index_col=0)
                        for index in tmp.index:
                            if tmp.loc[index,'REJECT']:
                                try:
                                    if level == 'taxaID':
                                        index = index.strip("X")
                                    biomarker_tb.loc[index,'ANCOM_identified'] = 1
                                except:
                                    pass
                        
                        # OUTPOST_biomarker_score
                        biomarker_tb['OUTPOST_biomarker_score'] = biomarker_tb.sum(axis = 1)
                        biomarker_tb = biomarker_tb.sort_values(by = 'OUTPOST_biomarker_score', ascending = False)
                        biomarker_tb.to_csv(f"{output_dir}/biomarkers_analysis/OUTPOST_biomarker_scores.{group1}_vs_{group2}.{level}.tsv",
                                            sep = '\t',header=True, index=True)
                        
                shell("touch {output_dir}/log/biomarker_summary.done")
                
                
        # %% plasmid
        rule plasmids_analysis_assembly_based:
            input:
                "{output_dir}/log/merge_counts_pure_and_kaiju.done"
            output:
                "{output_dir}/log/plasmids_analysis_assembly_based.done"
            log:
                "{output_dir}/log/plasmids_analysis_assembly_based.log"
            benchmark:
                "{output_dir}/benchmark/plasmids_analysis_assembly_based.benchmark"
            run:
                wkdir = f"{output_dir}/plasmids_analysis/assembly_based"
                os.makedirs(wkdir, exist_ok=True)
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{abricate} {assembly_dir} --db plasmidfinder --quiet > {wkdir}/{assembly}.plasmidfinder.tsv")
                        
                        taxa_level_temp = ','.join(taxa_level)
                        for group1,group2 in group_pair_list:
                            group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
                            group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
                            shell("{python3} OUTPOST/visualize_abricate.py {wkdir}/{assembly}.plasmidfinder.tsv {output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv"
                                  " {wkdir} {group1_index} {group1} {group2_index} {group2} {taxa_level_temp} True")
                        
                shell("touch {output_dir}/log/plasmids_analysis_assembly_based.done")
        
        
        # %% virulence
        rule virulence_factors_analysis_assembly_based:
            input:
                "{output_dir}/log/merge_counts_pure_and_kaiju.done"
            output:
                "{output_dir}/log/virulence_factors_analysis_assembly_based.done"
            log:
                "{output_dir}/log/virulence_factors_analysis_assembly_based.log"
            benchmark:
                "{output_dir}/benchmark/virulence_factors_analysis_assembly_based.benchmark"
            run:
                wkdir = f"{output_dir}/virulence_factors_analysis/assembly_based"
                os.makedirs(wkdir, exist_ok=True)
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{abricate} {assembly_dir} --db ecoli_vf --quiet > {wkdir}/{assembly}.ecoli_vf.tsv")
                        shell("{abricate} {assembly_dir} --db vfdb --quiet > {wkdir}/{assembly}.vfdb.tsv")
                        
                        shell("cp {wkdir}/{assembly}.vfdb.tsv {wkdir}/{assembly}.vfdb.temp.tsv")
                        shell("sed -i '1d' {wkdir}/{assembly}.vfdb.temp.tsv")
                        shell("cat {wkdir}/{assembly}.ecoli_vf.tsv {wkdir}/{assembly}.vfdb.temp.tsv > {wkdir}/{assembly}.virulence.tsv")
                        shell("rm {wkdir}/{assembly}.vfdb.temp.tsv")
                
                        taxa_level_temp = ','.join(taxa_level)
                        for group1,group2 in group_pair_list:
                            group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
                            group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
                            shell("{python3} OUTPOST/visualize_abricate.py {wkdir}/{assembly}.virulence.tsv {output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv "
                                  " {wkdir} {group1_index} {group1} {group2_index} {group2} {taxa_level_temp}")
                shell("touch {output_dir}/log/virulence_factors_analysis_assembly_based.done")

        
        
        # %% antibiotic
        rule antibiotic_genes_analysis_assembly_based:
            input:
                "{output_dir}/log/merge_counts_pure_and_kaiju.done"
            output:
                "{output_dir}/log/antibiotic_genes_analysis_assembly_based.done"
            log:
                "{output_dir}/log/antibiotic_genes_analysis_assembly_based.log"
            benchmark:
                "{output_dir}/benchmark/antibiotic_genes_analysis_assembly_based.benchmark"
            run:
                wkdir = f"{output_dir}/antibiotic_genes_analysis/assembly_based"
                os.makedirs(wkdir, exist_ok=True)
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{abricate} {assembly_dir} --db resfinder --quiet > {wkdir}/{assembly}.resfinder.tsv")
                        shell("{abricate} {assembly_dir} --db argannot --quiet > {wkdir}/{assembly}.argannot.tsv")
                        shell("{abricate} {assembly_dir} --db card --quiet > {wkdir}/{assembly}.card.tsv")
                        shell("{abricate} {assembly_dir} --db ecoh --quiet > {wkdir}/{assembly}.ecoh.tsv")
                        shell("{abricate} {assembly_dir} --db megares --quiet > {wkdir}/{assembly}.megares.tsv")
                        shell("{abricate} {assembly_dir} --db ncbi --quiet > {wkdir}/{assembly}.ncbiAMR.tsv")
        
                        shell("cp {wkdir}/{assembly}.argannot.tsv {wkdir}/{assembly}.argannot.temp.tsv")
                        shell("cp {wkdir}/{assembly}.card.tsv {wkdir}/{assembly}.card.temp.tsv")
                        shell("cp {wkdir}/{assembly}.ecoh.tsv {wkdir}/{assembly}.ecoh.temp.tsv")
                        shell("cp {wkdir}/{assembly}.megares.tsv {wkdir}/{assembly}.megares.temp.tsv")
                        shell("cp {wkdir}/{assembly}.ncbiAMR.tsv {wkdir}/{assembly}.ncbiAMR.temp.tsv")
                        shell("sed -i '1d' {wkdir}/{assembly}.argannot.temp.tsv")
                        shell("sed -i '1d' {wkdir}/{assembly}.card.temp.tsv")
                        shell("sed -i '1d' {wkdir}/{assembly}.ecoh.temp.tsv")
                        shell("sed -i '1d' {wkdir}/{assembly}.megares.temp.tsv")
                        shell("sed -i '1d' {wkdir}/{assembly}.ncbiAMR.temp.tsv")
                        shell("cat {wkdir}/{assembly}.resfinder.tsv {wkdir}/{assembly}.argannot.temp.tsv "
                              " {wkdir}/{assembly}.card.temp.tsv {wkdir}/{assembly}.ecoh.temp.tsv "
                              " {wkdir}/{assembly}.megares.temp.tsv {wkdir}/{assembly}.ncbiAMR.temp.tsv > {wkdir}/{assembly}.antibiotic.tsv")
                        shell("rm {wkdir}/{assembly}.argannot.temp.tsv")
                        shell("rm {wkdir}/{assembly}.card.temp.tsv")
                        shell("rm {wkdir}/{assembly}.ecoh.temp.tsv")
                        shell("rm {wkdir}/{assembly}.megares.temp.tsv")
                        shell("rm {wkdir}/{assembly}.ncbiAMR.temp.tsv")
                        
                        taxa_level_temp = ','.join(taxa_level)
                        for group1,group2 in group_pair_list:
                            group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
                            group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
                            shell("{python3} OUTPOST/visualize_abricate.py {wkdir}/{assembly}.antibiotic.tsv {output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv "
                                  " {wkdir} {group1_index} {group1} {group2_index} {group2} {taxa_level_temp}")
                shell("touch {output_dir}/log/antibiotic_genes_analysis_assembly_based.done")
        
        # %% ancom
        rule ancom_identification:
            input:
                "{output_dir}/log/counts_table2rel_abun.done"
            output:
                "{output_dir}/log/ancom_identification.done"
            log:
                "{output_dir}/log/ancom_identification.log"
            benchmark:
                "{output_dir}/benchmark/ancom_identification.benchmark"
            run:
                wkdir = f"{output_dir}/biomarkers_analysis/ANCOM_identification/"
                os.makedirs(wkdir, exist_ok=True)
                for level in taxa_level:
                    for group1,group2 in group_pair_list:
                        dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.csv"
                        group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
                        group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
                        shell("{Rscript} OUTPOST/fastANCOM.R {dp} {wkdir} {group1} {group2} {group1_index} {group2_index} {level}")

                shell("touch {output_dir}/log/ancom_identification.done")

                
                
        # %% alpha_beta_diversity
        rule alpha_beta_diversity:
            input:
                "{output_dir}/log/counts_table2rel_abun.done"
            output:
                "{output_dir}/log/alpha_beta_diversity.done"
            log:
                "{output_dir}/log/alpha_beta_diversity.log"
            benchmark:
                "{output_dir}/benchmark/alpha_beta_diversity.benchmark"
            threads: 1
            run:
                for group1,group2 in group_pair_list:
                    wkdir = f"{output_dir}/diversity_analysis/alpha_beta_{group1}_vs_{group2}"
                    os.makedirs(wkdir, exist_ok=True)
                    group1_new_name = [f"{sample}_{group1}" for sample in comparison_dict[group1]]
                    group2_new_name = [f"{sample}_{group2}" for sample in comparison_dict[group2]]

                    group1_name = comparison_dict[group1]
                    group2_name = comparison_dict[group2]
                    if ('species' in taxa_level) and ('genus' in taxa_level):
                        species = 'species'
                        genus = 'genus'
                    else:
                        species = taxa_level[-1]
                        genus = taxa_level[-2]
                        print(f"OUTPOST: rule alpha_beta_diversity species and genus not in {taxa_level}. use the last two {taxa_level[-1]} {taxa_level[-2]} instead.")
                    species_dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{species}.rmU.csv"
                    genus_dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{genus}.rmU.csv"
                    species_df = pd.read_csv(species_dp,index_col=0)
                    genus_df = pd.read_csv(genus_dp,index_col=0)

                    species_new_df = species_df.reindex(group1_name + group2_name)
                    genus_new_df = genus_df.reindex(group1_name + group2_name)

                    species_new_dp = wkdir + f"/all_samples.taxa_counts.rel_abun.{species}.rmU.{group1}_vs_{group2}.csv"
                    genus_new_dp = wkdir + f"/all_samples.taxa_counts.rel_abun.{genus}.rmU.{group1}_vs_{group2}.csv"

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

                    shell("{Rscript} OUTPOST/alpha_beta_diversity.R {species_new_dp} {wkdir} {index_file} {group_file} {group1},{group2} {species} > {log} 2>&1")
                    shell("{Rscript} OUTPOST/alpha_beta_diversity.R {genus_new_dp} {wkdir} {index_file} {group_file} {group1},{group2} {genus} > {log} 2>&1")
                shell("touch {output_dir}/log/alpha_beta_diversity.done")
        
        
        # %% batch_effect
        rule visualize_batch_effect:
            # the remove of batch effect was done in paste_counts_table rule
            # here we simply visualize it
            input:
                "{output_dir}/log/scale_taxa_rel_abun_table.done"
            output:
                "{output_dir}/log/visualize_batch_effect.done"
            log:
                "{output_dir}/log/visualize_batch_effect.log"
            benchmark:
                "{output_dir}/benchmark/visualize_batch_effect.benchmark"
            run:
                # visualize
                os.makedirs(f"{output_dir}/batch_effect", exist_ok=True)
                for level in taxa_level:
                    dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.csv"
                    if rm_batch_effect:
                        plot = f"{output_dir}/batch_effect/all_samples.taxa_counts.rel_abun.{level}.rmU.batch_effect_removed_PCA.pdf"
                    else:
                        plot = f"{output_dir}/batch_effect/all_samples.taxa_counts.rel_abun.{level}.rmU.PCA.pdf"
                    shell("{Rscript} OUTPOST/visualize_batch_effect.R {dp} {OUTPOST_config} {plot}")
                shell("touch {output_dir}/log/visualize_batch_effect.done")
                
                
                
        # %% taxonomy
        rule heatmap_taxa:
            input:
                "{output_dir}/log/scale_taxa_rel_abun_table.done"
            output:
                "{output_dir}/log/heatmap_taxa.done"
            log:
                "{output_dir}/log/heatmap_taxa.log"
            benchmark:
                "{output_dir}/benchmark/heatmap_taxa.benchmark"
            run:
                wkdir = f"{output_dir}/taxonomy_analysis/figs/"
                os.makedirs(wkdir, exist_ok=True)
                # top taxa heatmap
                data_type = 'taxa'
                for level in taxa_level:
                    dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.top{top}.fillmin.scaled.csv"
                    output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                    index_dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.top{top}.fillmin.scaled.index"
                    with open(index_dp,'w') as w:
                        for sample in sample_list:
                            w.write(sample + '\n')
                    try:
                        shell("{Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                    except: pass
                # top unequal taxa heatmap
                for group1,group2 in group_pair_list:
                    for level in taxa_level:
                        dp = f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.fillmin.scaled.csv"
                        output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                        index_dp = f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.fillmin.scaled.index"
                        with open(index_dp,'w') as w:
                            for sample in comparison_dict[group1]:
                                w.write(f"{sample}_{group1}\n")
                            for sample in comparison_dict[group2]:
                                w.write(f"{sample}_{group2}\n")
                        try:
                            shell("{Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                        except:
                            print(f"OUPOST error: {Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")

                # top equal taxa heatmap
                for group1,group2 in group_pair_list:
                    for level in taxa_level:
                        dp = f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.fillmin.scaled.csv"
                        output = os.path.basename(dp).replace(".csv",".heatmap.pdf")
                        index_dp = f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{30}.fillmin.scaled.index"
                        with open(index_dp,'w') as w:
                            for sample in comparison_dict[group1]:
                                w.write(f"{sample}_{group1}\n")
                            for sample in comparison_dict[group2]:
                                w.write(f"{sample}_{group2}\n")
                        try:
                            shell("{Rscript} OUTPOST/heatmap.R {dp} {wkdir} {index_dp} {data_type} {output} > {log}_taxa 2>&1")
                        except: pass
                shell("touch {output_dir}/log/heatmap_taxa.done")


        rule scale_taxa_rel_abun_table:
            input:
                "{output_dir}/log/extract_top_taxa.done"
            output:
                "{output_dir}/log/scale_taxa_rel_abun_table.done"
            log:
                "{output_dir}/log/scale_taxa_rel_abun_table.log"
            benchmark:
                "{output_dir}/benchmark/scale_taxa_rel_abun_table.benchmark"
            run:
                dp_list = [f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.top{top}.csv"
                               for level in taxa_level]
                output_list = [f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.top{top}.fillmin.scaled.csv"
                               for level in taxa_level]

                for group1,group2 in group_pair_list:
                    dp_list_equal = [f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.csv"
                                              for level in taxa_level]
                    output_list_equal = [f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.fillmin.scaled.csv"
                                                  for level in taxa_level]

                    dp_list_unequal = [f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.csv"
                                                for level in taxa_level]
                    output_list_unequal = [f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.fillmin.scaled.csv"
                                                    for level in taxa_level]

                    dp_list += dp_list_equal
                    dp_list += dp_list_unequal
                    output_list += output_list_equal
                    output_list += output_list_unequal
                dp_list = ','.join(dp_list)
                output_list = ','.join(output_list)
                try:
                    shell("{python3} OUTPOST/scale_rel_abun_table.py {dp_list} {output_list}")
                except:
                    # in case augment is too long
                    print("OUTPOST: rule scale_rel_abun_table Argument list too long. use chunks.")
                    dp_list_chunks = [dp_list.split(',')[i:i + 100] for i in range(0, len(dp_list.split(',')), 100)]
                    output_list_chunks = [output_list.split(',')[i:i + 100] for i in range(0, len(output_list.split(',')), 100)]
                    for dp_list_chunk,output_list_chunk in zip(dp_list_chunks, output_list_chunks):
                        dp_list_chunk = ','.join(dp_list_chunk)
                        output_list_chunk = ','.join(output_list_chunk)
                        shell("{python3} OUTPOST/scale_rel_abun_table.py {dp_list_chunk} {output_list_chunk}")
                shell("touch {output_dir}/log/scale_taxa_rel_abun_table.done")
        

        rule taxa_barplots:
            input:
                "{output_dir}/log/extract_top_taxa.done"
            output:
                "{output_dir}/log/taxa_barplots.done"
            log:
                "{output_dir}/log/taxa_barplots.log"
            benchmark:
                "{output_dir}/benchmark/taxa_barplots.benchmark"
            run:
                os.makedirs(f"{output_dir}/taxonomy_analysis/figs/", exist_ok=True)
                for level in taxa_level:
                    dp = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.top{top}.csv"
                    new_dp = dp.replace(".csv",".addOthers.csv")
                    output = f"{output_dir}/taxonomy_analysis/figs/" + os.path.basename(dp).replace(".csv",".barplot.pdf")
                    df = pd.read_csv(dp, index_col = 0)
                    df['Others'] = 1 - df.sum(axis = 1)
                    df.to_csv(new_dp)
                    try:
                        shell("{Rscript} OUTPOST/barplots.R {new_dp} {output} > {log} 2>&1")
                    except:
                        pass
                shell("touch {output_dir}/log/taxa_barplots.done")

        
        rule extract_top_taxa:
            input:
                "{output_dir}/log/rel_abun_utest.done"
            output:
                "{output_dir}/log/extract_top_taxa.done"
            log:
                "{output_dir}/log/extract_top_taxa.log"
            benchmark:
                "{output_dir}/benchmark/extract_top_taxa.benchmark"
            threads: cores
            run:
                # top abundant
                dp_list = ','.join([f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.csv"
                                    for level in taxa_level])
                output_name_list = ','.join([f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.top{top}.csv"
                                    for level in taxa_level])
                error_log = os.getcwd() + f"/{output_dir}/log/extract_top_taxa_top_abundant.error"
                shell("{python3} OUTPOST/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
                      " {top} {error_log}")
                result = wait_until_file_exists(output_name_list.split(','),error_log)
                assert result == True, "OUTPOST: extract_top_taxa_top_abundant has errors. exit."

                # top unequal abundant
                output_list = []
                for group1,group2 in group_pair_list:
                    os.makedirs(f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}", exist_ok=True)
                    dp_list = ','.join([f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.csv"
                                        for level in taxa_level])
                    output_name_list = ','.join([
                            f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.top{top}.csv"
                            for level in taxa_level])
                    error_log = os.getcwd() + f"/{output_dir}/log/extract_top_taxa_top_unequal_abundant.error"
                    output_list.append(f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_species.rel_abun.unequal.top{top}.csv")
                    shell("nohup {python3} OUTPOST/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
                      " {top} {error_log} > {log}_{group1}_vs_{group2}_unequal 2>&1 &")
                    time.sleep(3)
                result = wait_until_file_exists(output_list,error_log)
                assert result == True, "OUTPOST: extract_top_taxa_top_unequal_abundant has errors. exit."

               # top equal abundant
                output_list = []
                for group1,group2 in group_pair_list:
                    os.makedirs(f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}", exist_ok=True)
                    dp_list = ','.join([f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.csv"
                                        for level in taxa_level])
                    output_name_list = ','.join([
                            f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.top{top}.csv"
                            for level in taxa_level])
                    error_log = os.getcwd() + f"/{output_dir}/log/extract_top_taxa_top_equal_abundant.error"
                    output_list.append(f"{output_dir}/taxonomy_analysis/top_taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_species.rel_abun.equal.top{top}.csv")
                    shell("nohup {python3} OUTPOST/extract_top_taxa_by_rel_abun_from_rel_abun_table.py {dp_list} {output_name_list} "
                      " {top} {error_log} > {log}_{group1}_vs_{group2}_equal 2>&1 &")
                    time.sleep(3)
                result = wait_until_file_exists(output_list,error_log)
                assert result == True, "OUTPOST: extract_top_taxa_top_equal_abundant has errors. exit."

                shell("touch {output_dir}/log/extract_top_taxa.done")

        
        rule taxa_boxplot:
            input:
                "{output_dir}/log/rel_abun_utest.done"
            output:
                "{output_dir}/log/taxa_boxplot.done"
            log:
                "{output_dir}/log/taxa_boxplot.log"
            benchmark:
                "{output_dir}/benchmark/taxa_boxplot.benchmark"
            threads: 1
            run:
                output_list = []
                for group1,group2 in group_pair_list:
                    os.makedirs(f"{output_dir}/taxonomy_analysis/boxplot_{group1}_vs_{group2}", exist_ok=True)
                    group_pair = ','.join([group1,group2])
                    groups = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
                    dp_list = [f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.csv" for level in taxa_level]
                    for dp in dp_list:
                        output = f"{output_dir}/taxonomy_analysis/boxplot_{group1}_vs_{group2}/{os.path.basename(dp).replace('.csv','')}"
                        try:
                            shell("{Rscript} OUTPOST/boxplot.R {dp} {output} {groups} {group_pair} > {log} 2>&1")
                        except:
                            pass
                    # categorize boxplots
                    for level in taxa_level:
                        os.makedirs(f"{output_dir}/taxonomy_analysis/boxplot_{group1}_vs_{group2}/{level}/", exist_ok=True)
                        try:
                            shell("mv {output_dir}/taxonomy_analysis/boxplot_{group1}_vs_{group2}/*{level}.*pdf {output_dir}/taxonomy_analysis/boxplot_{group1}_vs_{group2}/{level}/")
                        except:
                            pass
                shell("touch {output_dir}/log/taxa_boxplot.done")
                
                
        rule lefse_taxa:
            input:
                "{output_dir}/log/rel_abun2lefse_taxa.done"
            output:
                "{output_dir}/log/lefse_taxa.done"
            log:
                "{output_dir}/log/lefse_taxa.log"
            benchmark:
                "{output_dir}/benchmark/lefse_taxa.benchmark"
            threads: cores
            run:
                wkdir = f"{output_dir}/LDA_analysis/figs_taxa/"
                os.makedirs(wkdir, exist_ok=True)
                # unequal taxa to lefse barplot
                for group1,group2 in group_pair_list:
                    dp_list = [f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.lefse.tsv"
                               for level in taxa_level]
                    format_list = [dp.replace(".tsv",".format") for dp in dp_list]
                    res_list = [dp.replace(".tsv",".res") for dp in dp_list]
                    output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

                    for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                        try:
                            shell("nohup {lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ; "
                                  "{lefse_run} {format_}  {res} > {log} 2>&1 ; "
                                  "{python3} OUTPOST/barplots_from_lefse_res_file.py  {res} {output} {LDA_cutoff} > {log} 2>&1 & ")
                        except:
                            print(f"OUTPOST: lefse_taxa: cannot process {dp}. skip.")

                # equal taxa to lefse barplot
                for group1,group2 in group_pair_list[:-2]:
                    dp_list = [f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.lefse.tsv"
                               for level in taxa_level]
                    format_list = [dp.replace(".tsv",".format") for dp in dp_list]
                    res_list = [dp.replace(".tsv",".res") for dp in dp_list]
                    output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

                    for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                        try:
                            shell("nohup {lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ; "
                                  "{lefse_run} {format_}  {res} > {log} 2>&1 ; "
                                  "{python3} OUTPOST/barplots_from_lefse_res_file.py  {res} {output} {LDA_cutoff} > {log} 2>&1 & ")
                        except:
                            print(f"OUTPOST: lefse_taxa: cannot process {dp}. skip.")

                # create time for calculation
                for group1,group2 in group_pair_list[-2:]:
                    dp_list = [f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.lefse.tsv"
                               for level in taxa_level]
                    format_list = [dp.replace(".tsv",".format") for dp in dp_list]
                    res_list = [dp.replace(".tsv",".res") for dp in dp_list]
                    output_list = [wkdir + os.path.basename(dp).replace(".tsv",".pdf") for dp in dp_list]

                    for dp,format_,res,output in zip(dp_list,format_list,res_list,output_list):
                        try:
                            shell("{lefse_format_input} {dp} {format_} -c 1 -s 2 -u 3 -o 1000000 > {log} 2>&1 ")
                            shell("{lefse_run} {format_}  {res} > {log} 2>&1 ")
                            shell("{python3} OUTPOST/barplots_from_lefse_res_file.py  {res} {output} {LDA_cutoff} > {log} 2>&1 ")
                        except:
                            print(f"OUTPOST: lefse_taxa: cannot process {dp}. skip.")

                shell("touch {output_dir}/log/lefse_taxa.done")
        
                
        rule rel_abun2lefse_taxa:
            input:
                "{output_dir}/log/rel_abun_utest.done"
            output:
                "{output_dir}/log/rel_abun2lefse_taxa.done"
            log:
                "{output_dir}/log/rel_abun2lefse_taxa.log"
            benchmark:
                "{output_dir}/benchmark/rel_abun2lefse_taxa.benchmark"
            threads: 1
            run:
                # unequal taxa to lefse table
                for group1,group2 in group_pair_list:
                    wkdir = f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/"
                    os.makedirs(wkdir, exist_ok=True)
                    dp_list = ','.join([f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.unequal.csv"
                                        for level in taxa_level])
                    output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                        for dp in dp_list.split(',')])
                    class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
                    shell("{python3} OUTPOST/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")
                # equal taxa to lefse table
                for group1,group2 in group_pair_list:
                    wkdir = f"{output_dir}/LDA_analysis/taxa_{group1}_vs_{group2}/"
                    os.makedirs(wkdir, exist_ok=True)
                    dp_list = ','.join([f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}.rel_abun.equal.csv"
                                        for level in taxa_level])
                    output_list = ','.join([wkdir + os.path.basename(dp).replace(".csv",".lefse.tsv")
                                        for dp in dp_list.split(',')])
                    class_ = ','.join([group1] * len(comparison_dict[group1]) + [group2] * len(comparison_dict[group2]))
                    shell("{python3} OUTPOST/rel_abun2lefse.py {dp_list} {output_list} {class_} > {log} 2>&1 ")

                shell("touch {output_dir}/log/rel_abun2lefse_taxa.done")


        rule rel_abun_utest:
            input:
                "{output_dir}/log/counts_table2rel_abun.done"
            output:
                "{output_dir}/log/rel_abun_utest.done"
            log:
                "{output_dir}/log/rel_abun_utest.log"
            benchmark:
                "{output_dir}/benchmark/rel_abun_utest.benchmark"
            threads: cores
            run:
                output_list = []
                for group1,group2 in group_pair_list:
                    os.makedirs(f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}", exist_ok=True)
                    dp_list = ','.join([f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts.rel_abun.{level}.rmU.csv" for level in taxa_level])
                    group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
                    group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
                    prefix_list = ','.join([f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_{level}"
                                            for level in taxa_level])
                    output = f"{output_dir}/taxonomy_analysis/utest_{group1}_vs_{group2}/all_samples.rel_abun.{group1}_vs_{group2}.at_species.ave_change.unequal.csv"
                    output_list.append(output)
                    error_log = os.getcwd() + f"/{output_dir}/log/rel_abun_utest.error"
                    try:
                        os.remove(error_log) # clean former residual error log
                    except:
                        pass
                    shell("nohup {python3} OUTPOST/rel_abun_utest.py {dp_list} {group1_index} {group1} "
                          " {group2_index} {group2} {prefix_list} {paired} {two_sided} {error_log} > {log}_{group1}_vs_{group2} 2>&1 &")
                    time.sleep(3)
                result = wait_until_file_exists(output_list,error_log)
                assert result == True, "OUTPOST: rel_abun_utest has errors. exit."
                shell("touch {output_dir}/log/rel_abun_utest.done")
                
                
        rule counts_table2rel_abun:
            input:
                "{output_dir}/log/merge_counts_pure_and_kaiju.done"
            output:
                "{output_dir}/log/counts_table2rel_abun.done"
            log:
                "{output_dir}/log/counts_table2rel_abun.log"
            benchmark:
                "{output_dir}/benchmark/counts_table2rel_abun.benchmark"
            run:
                input_dp = rf"{output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv"
                prefix = f"{output_dir}/taxonomy_analysis/all_samples.taxa_counts"
                samples = ','.join(sample_list)
                taxa_levels = ','.join(taxa_level)
                shell("{python3} OUTPOST/counts_table2rel_abun.py {input_dp} {prefix} {samples} {taxa_levels} > {log} 2>&1 ")
                shell("touch {output_dir}/log/counts_table2rel_abun.done")

        
        if not skip_assembly_qtest:
            rule assembly_qtest:
                input:
                    "{output_dir}/log/prepare_contig_table_from_counts_table.done",
                output:
                    "{output_dir}/log/assembly_qtest.done"
                log:
                    "{output_dir}/log/assembly_qtest.log"
                benchmark:
                    "{output_dir}/benchmark/assembly_qtest.benchmark"
                run:
                    for group1, group2 in group_pair_list:
                        dp = f"{output_dir}/assembly_analysis/taxa_counts/all_samples.{group1}_vs_{group2}.contig_table.tsv"
                        samples = ','.join(sample_list)
                        output = f"{output_dir}/assembly_analysis/taxa_counts/all_samples.{group1}_vs_{group2}.contig_table.processed.tsv"
                        shell("{python3} OUTPOST/assembly_qtest.py {dp} {samples} {output} > {log} 2>&1 ")
                    shell("touch {output_dir}/log/assembly_qtest.done")

            rule prepare_contig_table_from_counts_table:
                input:
                    "{output_dir}/log/summarize_assembly_table.done"
                output:
                    "{output_dir}/log/prepare_contig_table_from_counts_table.done"
                log:
                    "{output_dir}/log/prepare_contig_table_from_counts_table.log"
                benchmark:
                    "{output_dir}/benchmark/prepare_contig_table_from_counts_table.benchmark"
                threads: cores
                run:
                    input_dp = rf"{output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv"
                    output_list = []
                    for group1,group2 in group_pair_list:
                        group1_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group1]])
                        group2_index = ','.join([str(sample_list.index(sample)) for sample in comparison_dict[group2]])
                        samples = ','.join(sample_list)
                        output = f"{output_dir}/assembly_analysis/taxa_counts/all_samples.{group1}_vs_{group2}.contig_table.tsv"
                        output_list.append(output)
                        error_log = f"{output_dir}/log/prepare_contig_table_from_counts_table.error"
                        try:
                            os.remove(error_log) # clean former residual error log
                        except:
                            pass
                        shell("nohup {python3} OUTPOST/prepare_contig_table_from_counts_table.py {input_dp} "
                              " {group1_index} {group2_index} {output} {samples} {error_log} {qvalue} > {log} 2>&1 &")
                        time.sleep(3)
                    result = wait_until_file_exists(output_list,error_log)
                    assert result == True, "OUTPOST: prepare_contig_table_from_counts_table has errors. exit."
                    shell("touch {output_dir}/log/prepare_contig_table_from_counts_table.done")
        else:
            rule skip_assembly_qtest:
                input:
                    "{output_dir}/log/summarize_assembly_table.done"
                output:
                    "{output_dir}/log/assembly_qtest.done"
                log:
                    "{output_dir}/log/assembly_qtest.log"
                benchmark:
                    "{output_dir}/benchmark/assembly_qtest.benchmark"
                run:
                    shell("touch {output_dir}/log/assembly_qtest.done")
                    
        rule summarize_assembly_table:
            input:
                "{output_dir}/log/merge_counts_pure_and_kaiju.done"
            output:
                "{output_dir}/log/summarize_assembly_table.done"
            log:
                "{output_dir}/log/summarize_assembly_table.log"
            benchmark:
                "{output_dir}/benchmark/summarize_assembly_table.benchmark"
            run:
                os.makedirs(f"{output_dir}/assembly_analysis/", exist_ok=True)
                dp = f"{output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv"
                df = pd.read_csv(dp, sep = '\t', header = None, index_col=0)
                df.columns = taxa_level + sample_list
                summarized_df_columns = taxa_level + ['CPM_sum']
                summarized_df = pd.DataFrame(columns = summarized_df_columns)
                for sample in sample_list:
                    for level in taxa_level:
                        summarized_df.loc[sample,level] = len([x for x in list(set(df[level].dropna())) if x != 'NA'])
                    summarized_df.loc[sample,'CPM_sum'] = df[sample].sum()
                summarized_df.to_csv(f"{output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.summarized.tsv",
                                     sep = '\t',header=True, index=True)
                shell("touch {output_dir}/log/summarize_assembly_table.done")

        rule merge_counts_pure_and_kaiju:
            input:
                "{output_dir}/log/concat_kaiju_output.done",
                "{output_dir}/log/paste_counts_table.done"
            output:
                "{output_dir}/log/merge_counts_pure_and_kaiju.done"
            log:
                "{output_dir}/log/merge_counts_pure_and_kaiju.log"
            benchmark:
                "{output_dir}/benchmark/merge_counts_pure_and_kaiju.benchmark"
            run:
                os.makedirs(f"{output_dir}/assembly_analysis/taxa_counts", exist_ok=True)
                taxa_len = len(taxa_level)
                counts_dp = fr"{output_dir}/taxonomy_analysis/temp/all_samples.counts.pure"
                kaiju_dp = rf"{output_dir}/taxonomy_analysis/kaiju/all_assembly_kaiju.ref.nm.tsv"
                shell("{Rscript} OUTPOST/merge_counts_and_kaiju.R {counts_dp} {kaiju_dp} {output_dir}/assembly_analysis/taxa_counts/all_samples.taxa_counts.tsv {taxa_len}")
                shell("touch {output_dir}/log/merge_counts_pure_and_kaiju.done")

    rule concat_kaiju_output:
            input:
                "{output_dir}/log/format_kaiju_output.done"
            output:
                "{output_dir}/log/concat_kaiju_output.done"
            log:
                "{output_dir}/log/concat_kaiju_output.log"
            benchmark:
                "{output_dir}/benchmark/concat_kaiju_output.benchmark"
            run:
                kaiju_output_list = " ".join([f"{output_dir}/taxonomy_analysis/kaiju/{x}_kaiju.ref.nm.tsv" for x in set(assembly_list)])
                shell("cat {kaiju_output_list} > {output_dir}/taxonomy_analysis/kaiju/all_assembly_kaiju.ref.nm.tsv")
                shell("touch {output_dir}/log/concat_kaiju_output.done")
                
    if not skip_kaiju:
        rule format_kaiju_output:
            input:
                "{output_dir}/log/kaiju_addTaxonNames.done"
            output:
                "{output_dir}/log/format_kaiju_output.done"
            log:
                "{output_dir}/log/format_kaiju_output.log"
            benchmark:
                "{output_dir}/benchmark/format_kaiju_output.benchmark"
            run:
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{python3} OUTPOST/format_kaiju_output_to_tab_seperated.py {output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref.nm")
                        shell("touch {output_dir}/log/format_kaiju_output.done")
    
        rule kaiju_addTaxonNames:
            input:
                "{output_dir}/log/kaiju_annotate.done"
            output:
                "{output_dir}/log/kaiju_addTaxonNames.done"
            log:
                "{output_dir}/log/kaiju_addTaxonNames.log"
            benchmark:
                "{output_dir}/benchmark/kaiju_addTaxonNames.benchmark"
            threads: 1
            run:
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{kaiju_addTaxonNames} -i {output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref -o {output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref.nm -t  {kaiju_nodes} -n {kaiju_names} "
                        " -v -r superkingdom,phylum,class,order,family,genus,species > {log} 2>&1 ")
                        shell("touch {output_dir}/log/kaiju_addTaxonNames.done")
        # wait for new assembly finish if there is new assembly
        if (len(r1_fq_list + se_fq_list) > 0) and assemble_contigs:
            kaiju_annotate_input = f"{output_dir}/log/outpost_contigs.done"
        else:
            kaiju_annotate_input = [assembly_list,assembly_dir_list]
        rule kaiju_annotate:
            input:
                kaiju_annotate_input
            output:
                "{output_dir}/log/kaiju_annotate.done"
            threads: cores
            log:
                "{output_dir}/log/kaiju_annotate.log"
            benchmark:
                "{output_dir}/benchmark/kaiju_annotate.benchmark"
            run:
                os.makedirs(f"{output_dir}/taxonomy_analysis/kaiju/", exist_ok = True)
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{kaiju} -t {kaiju_nodes} -v -f {kaiju_fmi} -z {threads} -i {assembly_dir}  -o {output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref > {log} 2>&1 ")
                        shell("touch {output_dir}/log/kaiju_annotate.done")
    else:
        if (len(r1_fq_list + se_fq_list) > 0) and assemble_contigs:
            skip_kaiju_input = f"{output_dir}/log/outpost_contigs.done"
        else:
            skip_kaiju_input = [assembly_list,assembly_dir_list]
        rule skip_kaiju:
            input:
                skip_kaiju_input
            output:
                "{output_dir}/log/format_kaiju_output.done"
            log:
                "{output_dir}/log/skip_kaiju.log"
            benchmark:
                "{output_dir}/benchmark/skip_kaiju.benchmark"
            run:
                for assembly,assembly_dir in set(zip(assembly_list,assembly_dir_list)):
                    if assembly != '':
                        shell("{python3} OUTPOST/format_kaiju_output_to_tab_seperated.py {output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref.nm")
                        files = [f"{output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref",\
                                 f"{output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref.nm",\
                                 f"{output_dir}/taxonomy_analysis/kaiju/{assembly}_kaiju.ref.nm.tsv"]
                        assert all([os.path.exists(x) for x in files]),\
                        f"OUTPOST: detected missing kaiju files under {output_dir}/taxonomy_analysis/kaiju/. cannot skip kaiju. exit."
                        os.makedirs(f"{output_dir}/log", exist_ok=True)
                        shell("touch {output_dir}/log/format_kaiju_output.done")


#%% bam
if (len(r1_fq_list + se_fq_list) > 0 and len(assembly_list) > 0) or (len(bam_dir_list) > 0):
    rule paste_counts_table:
        input:
            "{output_dir}/log/idxstats.done"
        output:
            "{output_dir}/log/paste_counts_table.done"
        log:
            "{output_dir}/log/paste_counts_table.log"
        benchmark:
            "{output_dir}/benchmark/paste_counts_table.benchmark"
        run:
            counts_table_list = ','.join([f"{output_dir}/taxonomy_analysis/temp/{basename}.idxstats" for basename in bam_basename])
            samples_index = ''
            for group,samples in comparison_dict.items():
                # the sample_list order is the same as bam_list
                sample_index = '_'.join([str(sample_list.index(sample)) for sample in samples])
                samples_index = samples_index + ',' + sample_index
            # samples_index is 0_1_2,3_4_5
            batch_list_str = ','.join([str(_) for _ in batch_list])
            shell("{python3} OUTPOST/paste_counts_table.py 3 {counts_table_list} {output_dir}/taxonomy_analysis/temp/all_samples.counts.pure {batch_list_str} {rm_batch_effect} ")
            shell("touch {output_dir}/log/paste_counts_table.done")
            
    rule idxstats:
        input:
            "{output_dir}/log/align2assembly.done"
        output:
            "{output_dir}/log/idxstats.done"
        log:
            "{output_dir}/log/idxstats.log"
        benchmark:
            "{output_dir}/benchmark/idxstats.benchmark"
        threads: cores
        run:
            os.makedirs(f"{output_dir}/taxonomy_analysis/temp/", exist_ok=True)
            for bam,basename in zip(bam_dir_list,bam_basename):
                if basename != '':
                    shell(f"{samtools} idxstats --threads {threads} {bam} > {output_dir}/taxonomy_analysis/temp/{basename}.idxstats ")
                else:
                    shell(f"touch {output_dir}/taxonomy_analysis/temp/{basename}.idxstats ")
            shell("touch {output_dir}/log/idxstats.done")


if len(r1_fq_list + se_fq_list) > 0 and len(assembly_list) > 0:
    if assemble_contigs:
        index_assembly_input = [f"{output_dir}/log/outpost_contigs.done"] +\
                               [f"{output_dir}/log/reads_QC_PE.done"] * bool(len(r1_fq_list) > 0) +\
                               [f"{output_dir}/log/reads_QC_SE.done"] * bool(len(se_fq_list) > 0)
    else:
        index_assembly_input = [f"{output_dir}/log/reads_QC_PE.done"] * bool(len(r1_fq_list) > 0) +\
                               [f"{output_dir}/log/reads_QC_SE.done"] * bool(len(se_fq_list) > 0)
   
    rule index_assembly:
        input:
            index_assembly_input
        output:
            "{output_dir}/log/index_assembly.done"
        log:
            "{output_dir}/log/index_assembly.log"
        benchmark:
            "{output_dir}/benchmark/index_assembly.benchmark"
        threads: 1
        run:
            for i, sample in enumerate(sample_list):
                assembly = assembly_list[i]
                assembly_dir = assembly_dir_list[i]
                if assembly != '':
                    index_genome(assembly_dir, bwa, samtools)
            shell("touch {output_dir}/log/index_assembly.done")
                    
    rule align2assembly:
        input:
            "{output_dir}/log/index_assembly.done"
        output:
            "{output_dir}/log/align2assembly.done"
        log:
            "{output_dir}/log/align2assembly.log"
        benchmark:
            "{output_dir}/benchmark/align2assembly.benchmark"
        threads: cores
        run:
            for i, sample in enumerate(sample_list):
                nonvirus_fq = f"{output_dir}/data/{sample}_nonhostvirus.fq"
                assembly = assembly_list[i]
                assembly_dir = assembly_dir_list[i]
                align2assembly = False
                if assembly != '':
                    output_sam = f"{output_dir}/data/{sample}_align2{assembly}_unsorted.sam"
                    output_bam = f"{output_dir}/data/{sample}_align2{assembly}_unsorted.bam"
                    output_bam1 = f"{output_dir}/data/{sample}_align2{assembly}_sorted.bam"
                    # update the bam_dir_list
                    try:
                        if bam_dir_list[i] == '':
                            bam_dir_list[i] = output_bam1
                            align2assembly = True
                    except:
                        bam_dir_list.append(output_bam1)
                        align2assembly = True
                    if align2assembly == True:
                        shell("{bwa} mem -t {threads} {assembly_dir} {nonvirus_fq} > {output_sam}")
                        shell("{samtools} view -@ {threads} -bS -F 4 {output_sam} > {output_bam}")
                        shell("{samtools} sort -@ {threads} {output_bam} -o {output_bam1}")
                        shell("{samtools} index -@ {threads} {output_bam1}")
                        
                        if clean_unnecessary:
                            try:
                                os.remove(output_sam)
                                os.remove(output_bam)
                            except Exception as e:
                                pass
            # update bam_basename
            bam_basename=[os.path.basename(x) for x in bam_dir_list]
            # update bam_dir_list in config.tsv
            df_config['bam_dir'] = bam_dir_list
            df_config.to_csv(OUTPOST_config,sep='\t', index = None)
            
            shell("touch {output_dir}/log/align2assembly.done")
            
elif len(bam_dir_list) > 0:
    if (len(r1_fq_list + se_fq_list) > 0) and assemble_contigs:
        align2assembly_empty_input = f"{output_dir}/log/outpost_contigs.done"
    else:
        align2assembly_empty_input = bam_dir_list
    rule align2assembly_empty:
        input:
            align2assembly_empty_input
        output:
            "{output_dir}/log/align2assembly.done"
        log:
            "{output_dir}/log/align2assembly.log"
        benchmark:
            "{output_dir}/benchmark/align2assembly.benchmark"
        threads: 1
        run:
            shell("touch {output_dir}/log/align2assembly.done")


#%% assemble
if len(r1_fq_list + se_fq_list) > 0:
    rule quast_contigs:
        input:
            "{output_dir}/log/outpost_contigs.done"
        output:
            "{output_dir}/log/quast_contigs.done"
        log:
            "{output_dir}/log/quast_contigs.log"
        benchmark:
            "{output_dir}/benchmark/quast_contigs.benchmark"
        threads: cores
        run:
            quast_output = f"{output_dir}/assembly_analysis/quast"
            os.makedirs(quast_output, exist_ok=True)
            contig_dp = f"{output_dir}/assembly_analysis/outpost_contigs/outpost_nonrd_contigs.fasta"
            shell("{quast} {contig_dp} -t {threads} --mgm --rna-finding --circos -o {quast_output}")
            shell("touch {output_dir}/log/quast_contigs.done")
    
    rule outpost_contigs:
        input:
            "{output_dir}/log/assemble_contigs.done"
        output:
            "{output_dir}/log/outpost_contigs.done"
        log:
            "{output_dir}/log/outpost_contigs.log"
        benchmark:
            "{output_dir}/benchmark/outpost_contigs.benchmark"
        threads: cores
        run:
            os.makedirs(f"{output_dir}/assembly_analysis/outpost_contigs", exist_ok=True)
            if assembly_method == 'megahit':
                contig_dp = f"{output_dir}/assembly_analysis/contigs/final.contigs.fa"
            elif assembly_method == 'metaspades':
                contig_dp = f"{output_dir}/assembly_analysis/contigs/contigs.fasta"
            # remove redundancy
            nonrd_outpost_contigs = f"{output_dir}/assembly_analysis/outpost_contigs/outpost_nonrd_contigs.fasta"
            shell("{cdhit} -i {contig_dp} -o {nonrd_outpost_contigs} -c {cdhit_cutoff} -n 10 -M 0 -T {threads}")
            shell("cat {nonrd_outpost_contigs} | {seqkit} seq -w 60 > {nonrd_outpost_contigs}.tmp")
            shell("mv {nonrd_outpost_contigs}.tmp {nonrd_outpost_contigs}")
            shell("touch {output_dir}/log/outpost_contigs.done")
    
    rule assemble_contigs:
        input:
            "{output_dir}/log/cat_reads.done"
        output:
            "{output_dir}/log/assemble_contigs.done"
        log:
            "{output_dir}/log/assemble_contigs.log"
        benchmark:
            "{output_dir}/benchmark/assemble_contigs.benchmark"
        threads: cores
        run:
            contigs_output = f"{output_dir}/assembly_analysis/contigs"
            # os.makedirs(f"{output_dir}/assembly_analysis/figs", exist_ok=True)
            
            if assembly_method == 'megahit':
                command_megahit = "-1 " + ",".join([_ for _ in r1_fq_list if _ != '']) +\
                                  " -2 " + ",".join([_ for _ in r2_fq_list if _ != '']) +\
                                  " -r " + ",".join([_ for _ in se_fq_list if _ != ''])
                shell("{megahit} --continue -t {threads} {command_megahit} -o {contigs_output}")
            elif assembly_method == 'metaspades':
                command_spades = ""
                for r1_fq, r2_fq in zip(r1_fq_list, r2_fq_list):
                    if r1_fq != '':
                        command_spades += rf" --pe1-1 {r1_fq} --pe1-2 {r1_fq} "
                for se_fq in se_fq_list:
                    if se_fq != '':
                        command_spades += rf" --s1 {se_fq} "
                shell("{metaspades}  {command_spades} -t {threads} -o {contigs_output}")
            shell("touch {output_dir}/log/assemble_contigs.done")

    #%% reads QC
    cat_reads_input = [f"{output_dir}/log/reads_QC_PE.done"] * bool(len(r1_fq_list) > 0) +\
                      [f"{output_dir}/log/reads_QC_SE.done"] * bool(len(se_fq_list) > 0)
    if len(r1_fq_list) > 0:
        rule cat_reads:
            input:
                cat_reads_input
            output:
                "{output_dir}/log/cat_reads.done"
            log:
                "{output_dir}/log/cat_reads.log"
            benchmark:
                "{output_dir}/benchmark/cat_reads.benchmark"
            threads: 1
            run:
                for i, sample in enumerate(sample_list):
                    has_pe_reads, has_se_reads = False, False
                    r1 = r1_fq_list[i]
                    se = se_fq_list[i]
                    if r1 != '':
                        has_pe_reads = True
                    if se != '':
                        has_se_reads = True
                    
                    nonvirus_r1 = f"{output_dir}/data/{sample}_nonhostvirus_r1.fq" * has_pe_reads
                    nonvirus_r2 = f"{output_dir}/data/{sample}_nonhostvirus_r2.fq" * has_pe_reads
                    nonvirus_se = f"{output_dir}/data/{sample}_nonhostvirus_se.fq" * has_se_reads
                    nonvirus = f"{output_dir}/data/{sample}_nonhostvirus.fq"
                    shell("cat {nonvirus_r1} {nonvirus_r2} {nonvirus_se} > {nonvirus}")
                shell("touch {output_dir}/log/cat_reads.done")
                
                
        rule reads_QC_PE:
            input:
                "{output_dir}/log/fastp_PE.done"
            output:
                "{output_dir}/log/reads_QC_PE.done"
            log:
                "{output_dir}/log/reads_QC_PE.log"
            benchmark:
                "{output_dir}/benchmark/reads_QC_PE.benchmark"
            threads: cores
            run:
                # remove host contamination
                index_genome(virus_genome, bwa, samtools)
                for i, sample in enumerate(sample_list):
                    r1 = r1_fq_list[i]
                    if r1 == '':
                        continue
                    r1 = f"{output_dir}/data/{sample}.fp.r1.fastq"
                    r2 = f"{output_dir}/data/{sample}.fp.r2.fastq"
                    
                    host_genome = refgenome_list[i]
                    trimm_log = f"{output_dir}/qc/TrimSummmry_{sample}.txt"
                    nonhost_bam = f"{output_dir}/data/{sample}_onlyPE_nonhost_unsorted.bam"
                    nonvirus_bam = f"{output_dir}/data/{sample}_onlyPE_nonhostvirus_unsorted.bam"
                    nonvirus_bam1 = f"{output_dir}/data/{sample}_onlyPE_nonhostvirus_sorted.bam"
                    
                    nonhost_r1 = f"{output_dir}/data/{sample}_nonhost_r1.fq"
                    nonhost_r2 = f"{output_dir}/data/{sample}_nonhost_r2.fq"
                    nonvirus_r1 = f"{output_dir}/data/{sample}_nonhostvirus_r1.fq"
                    nonvirus_r2 = f"{output_dir}/data/{sample}_nonhostvirus_r2.fq"
                    
                    if len(refgenome_list) > 0:
                        # Index the host genome
                        index_genome(host_genome, bwa, samtools)
                        # Trimming reads with Trimmomatic
                        shell("{trimmomatic} PE -phred33 {r1} {r2} \
                         {r1}_trim.fq {r1}_trim_UnP.fq \
                         {r2}_trim.fq {r2}_trim_UnP.fq \
                         ILLUMINACLIP:{adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:32 \
                        > {trimm_log} ")
                        # Align reads to host genome and extract unmapped reads
                        shell("{bwa} mem -t {threads} {host_genome} {r1}_trim.fq {r2}_trim.fq \
                               | {samtools} view -@ {threads} -bS -f 12 -F 256  > {nonhost_bam}")
                        shell("{bamToFastq} -i {nonhost_bam} -fq {nonhost_r1} -fq2 {nonhost_r2}")
                        if rm_virus_contamination:
                            # Align to virus genome
                            # virus_aligned_sam = f"{output_dir}/data/{sample}_virus_aligned.sam"
                            shell("{bwa} mem -t {threads} {virus_genome} {nonhost_r1} {nonhost_r2} \
                                   | {samtools} view -@ {threads} -bS -f 12 -F 256 > {nonvirus_bam}")
                            shell("{samtools} sort -@ {threads} {nonvirus_bam} -o {nonvirus_bam1}")
                            shell("{samtools} index -@ {threads} {nonvirus_bam1}")
                            shell("{bamToFastq} -i {nonvirus_bam1} -fq {nonvirus_r1} -fq2 {nonvirus_r2}")
                        else:
                            shell("cp {nonhost_r1} {nonvirus_r1}")
                            shell("cp {nonhost_r2} {nonvirus_r2}")
                        # Cleanup intermediate files
                        if clean_unnecessary:
                            try:
                                os.remove(nonhost_bam)
                                os.remove(nonhost_r1)
                                os.remove(nonhost_r2)
                                if rm_virus_contamination:
                                    os.remove(nonvirus_bam)
                                os.remove(f"{output_dir}/data/{r1}_trim.fq")
                                os.remove(f"{output_dir}/data/{r2}_trim.fq")
                                os.remove(f"{output_dir}/data/{r1}_trim_Unp.fq")
                                os.remove(f"{output_dir}/data/{r2}_trim_Unp.fq")
                            except Exception as e:
                                pass
                    else:
                        if rm_virus_contamination:
                            # Align to virus genome
                            # virus_aligned_sam = f"{output_dir}/data/{sample}_virus_aligned.sam"
                            shell("{bwa} mem -t {threads} {virus_genome} {r1} {r2} \
                                   | {samtools} view -@ {threads} -bS -f 12 -F 256 > {nonvirus_bam}")
                            shell("{samtools} sort -@ {threads} {nonvirus_bam} -o {nonvirus_bam1}")
                            shell("{samtools} index -@ {threads} {nonvirus_bam1}")
                            shell("{bamToFastq} -i {nonvirus_bam1} -fq {nonvirus_r1} -fq2 {nonvirus_r2}")
                        else:
                            shell("cp {r1} {nonvirus_r1}")
                            shell("cp {r2} {nonvirus_r2}")
                shell("touch {output_dir}/log/reads_QC_PE.done")
                    
                
        rule fastp_PE:
            input:
                "{output_dir}/log/downsample_PE.done"
            output:
                "{output_dir}/log/fastp_PE.done"
            log:
                "{output_dir}/log/fastp_PE.log"
            benchmark:
                "{output_dir}/benchmark/fastp_PE.benchmark"
            threads: 1
            run:
                for i,sample in enumerate(sample_list):
                    r1 = r1_fq_list[i]
                    r2 = r2_fq_list[i]
                    if r1 == '':
                        continue
                    r1_output = f"{output_dir}/data/{sample}.fp.r1.fastq"
                    r2_output = f"{output_dir}/data/{sample}.fp.r2.fastq"
                    json = f"{output_dir}/qc/{sample}.fastp.json"
                    html = f"{output_dir}/qc/{sample}.fastp.html"
                    if downsample_reads:
                        r1_downsample = f"{output_dir}/data/{sample}.ds.r1.fastq"
                        r2_downsample = f"{output_dir}/data/{sample}.ds.r2.fastq"
                        shell("{fastp} -i {r1_downsample} -I {r2_downsample} -j {json} \
                               -h {html} -o {r1_output}  -O {r2_output} > {log} 2>&1")
                    else:
                        shell("{fastp} -i {r1} -I {r2} -j {json} -h {html} -o {r1_output} \
                               -O {r2_output} > {log} 2>&1")
                shell("touch {output_dir}/log/fastp_PE.done")
                
                
        rule downsample_PE:
            input:
                [_ for _ in r1_fq_list if _ != ''],
                [_ for _ in r2_fq_list if _ != '']
            output:
                "{output_dir}/log/downsample_PE.done"
            log:
                "{output_dir}/log/downsample_PE.log"
            benchmark:
                "{output_dir}/benchmark/downsample_PE.benchmark"
            threads: 1
            run:
                os.makedirs(f"{output_dir}/data", exist_ok=True)
                os.makedirs(f"{output_dir}/qc", exist_ok=True)
                os.makedirs(f"{output_dir}/log", exist_ok=True)
                
                for i,sample in enumerate(sample_list):
                    r1 = r1_fq_list[i]
                    r2 = r2_fq_list[i]
                    if r1 == '':
                        continue
                    # downsample reads if too large
                    if downsample_reads:
                        r1_downsample = f"{output_dir}/data/{sample}.ds.r1.fastq"
                        r2_downsample = f"{output_dir}/data/{sample}.ds.r2.fastq"
                        shell("{seqkit} head -n {downsample_reads} {r1} > {r1_downsample}")
                        shell("{seqkit} head -n {downsample_reads} {r2} > {r2_downsample}")
                shell("touch {output_dir}/log/downsample_PE.done")
    
                    
    if len(se_fq_list) > 0:
        rule reads_QC_SE:
            input:
                "{output_dir}/log/fastp_SE.done"
            output:
                "{output_dir}/log/reads_QC_SE.done",
            log:
                "{output_dir}/log/reads_QC_SE.log"
            benchmark:
                "{output_dir}/benchmark/reads_QC_SE.benchmark"
            threads: cores
            run:
                # remove host contamination
                index_genome(virus_genome, bwa, samtools)
                for i, sample in enumerate(sample_list):
                    se = se_fq_list[i]
                    if se == '':
                        continue
                    
                    se = f"{output_dir}/data/{sample}.fp.se.fastq"
                    host_genome = refgenome_list[i]
                    trimm_log = f"{output_dir}/qc/TrimSummmry_{sample}.txt"
                    nonhost_bam = f"{output_dir}/data/{sample}_onlySE_nonhost_unsorted.bam"
                    nonhost_sam = f"{output_dir}/data/{sample}_onlySE_nonhost_unsorted.sam"
                    nonvirus_bam = f"{output_dir}/data/{sample}_onlySE_nonhostvirus_unsorted.bam"
                    nonvirus_bam1 = f"{output_dir}/data/{sample}_onlySE_nonhostvirus_sorted.bam"
                    nonhost_se = f"{output_dir}/data/{sample}_nonhost_se.fq"
                    nonvirus_se = f"{output_dir}/data/{sample}_nonhostvirus_se.fq"
                    
                    index_genome(host_genome, bwa, samtools)
                    if len(refgenome_list) > 0:
                        # Index the host genome
                        # Trimming reads with Trimmomatic
                        shell("{trimmomatic} SE -phred33 {se} {se}_trim.fq \
                         ILLUMINACLIP:{adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:32 \
                        > {trimm_log} ")
                        # Align reads to host genome and extract unmapped reads
                        shell("{bwa} mem -t {threads} {host_genome} {se}_trim.fq \
                               | {samtools} view -@ {threads} -bS -f 4 > {nonhost_bam}")
                        shell("{bamToFastq} -i {nonhost_bam} -fq {nonhost_se} ")
                        # Align to virus genome
                        # virus_aligned_sam = f"{output_dir}/data/{sample}_virus_aligned.sam"
                        shell("{bwa} mem -t {threads} {virus_genome} {nonhost_se} \
                               | {samtools} view -@ {threads} -bS -f 4 > {nonvirus_bam}")
                        shell("{samtools} sort -@ {threads} {nonvirus_bam} -o {nonvirus_bam1}")
                        shell("{samtools} index -@ {threads} {nonvirus_bam1}")
                        shell("{bamToFastq} -i {nonvirus_bam1} -fq {nonvirus_se}")
                        
                        # Cleanup intermediate files
                        if clean_unnecessary:
                            try:
                                os.remove(nonhost_bam)
                                os.remove(nonvirus_bam)
                                os.remove(nonhost_se)
                                os.remove(f"{output_dir}/data/{se}_trim.fq")
                                os.remove(f"{output_dir}/data/{se}_trim_Unp.fq")
                            except Exception as e:
                                pass
                    else:shell("cp {se} {nonvirus_se}")
                
                shell("touch {output_dir}/log/reads_QC_SE.done")
    
        rule fastp_SE:
            input:
                "{output_dir}/log/downsample_SE.done"
            output:
                "{output_dir}/log/fastp_SE.done"
            log:
                "{output_dir}/log/fastp_SE.log"
            benchmark:
                "{output_dir}/benchmark/fastp_SE.benchmark"
            threads: 1
            run:
                for i,sample in enumerate(sample_list):
                    se = se_fq_list[i]
                    if se == '':
                        continue
                    se_output = f"{output_dir}/data/{sample}.fp.se.fastq"
                    json = f"{output_dir}/qc/{sample}.fastp.json"
                    html = f"{output_dir}/qc/{sample}.fastp.html"
                    if downsample_reads:
                        se_downsample = f"{output_dir}/data/{sample}.ds.se.fastq"
                        shell("{fastp} -i {se_downsample} -j {json} \
                               -h {html} -o {se_output}  > {log} 2>&1")
                    else:
                        shell("{fastp} -i {se} -j {json} -h {html} -o {se_output} \
                               > {log} 2>&1")
                shell("touch {output_dir}/log/fastp_SE.done")
                
                
        rule downsample_SE:
            input:
                [_ for _ in se_fq_list if _ != '']
            output:
                "{output_dir}/log/downsample_SE.done"
            log:
                "{output_dir}/log/downsample_SE.log"
            benchmark:
                "{output_dir}/benchmark/downsample_SE.benchmark"
            threads: 1
            run:
                os.makedirs(f"{output_dir}/data", exist_ok=True)
                os.makedirs(f"{output_dir}/qc", exist_ok=True)
                os.makedirs(f"{output_dir}/log", exist_ok=True)
                
                for i,sample in enumerate(sample_list):
                    se = se_fq_list[i]
                    if se == '':
                        continue
                    # downsample reads if too large
                    if downsample_reads:
                        se_downsample = f"{output_dir}/data/{sample}.ds.se.fastq"
                        shell("{seqkit} head -n {downsample_reads} {se} > {se_downsample}")
                shell("touch {output_dir}/log/downsample_SE.done")




















