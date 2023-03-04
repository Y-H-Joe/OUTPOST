# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:05:54 2023

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

================================ description ==================================
concatenate genes, contigs, samples and taxa.

generate and save intermediate tables.
visulize using two-dim cluster heatmap.
=================================== input =====================================
#FILE	SEQUENCE	START	END	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION
GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010003464.1	186	633	Streptomyces_cinnamoneus_EF-Tu_mutants_conferring_resistance_to_elfamycin	713-1163/1194	......../======	11-Apr	37.19	75.82	card	X98831.1:362-1556
GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010003464.1	731	831	Streptomyces_cinnamoneus_EF-Tu_mutants_conferring_resistance_to_elfamycin	2-102/1194	==.............	0/0	8.46	85.15	card	X98831.1:362-1556
GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010006077.1	1	437	emrK	285-721/1057	....=======....	0/0	41.34	97.25	card	D78168:537-1593

JAIZPE010000001.1	1262830	Bacteria	Firmicutes	Clostridia	Eubacteriales	Clostridiaceae	Clostridium	Clostridium sp. CAG:632	6	0	9	12	0	1	0	0	0	0	0	0	0	11	11	0
JAIZPE010000002.1	1262792	Bacteria	Firmicutes	Clostridia	Eubacteriales	Clostridiaceae	Clostridium	Clostridium sp. CAG:299	174	78	1420	353	170	153	6	38	114	21	6	24	63	51	284	10

=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================
"""
import sys
from os import path
from itertools import chain
import pandas as pd
import heapq
import numpy as np
import string
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


def visualize_abricate(gene_tb_dp, counts_tb_dp, output_db, 
                   group1_index, group1_name, group2_index, group2_name, 
                   taxa_level, plasmid):
    genes_tb = pd.read_csv(gene_tb_dp, sep = '\t')
    counts_tb = pd.read_csv(counts_tb_dp, sep = '\t', header = None)
    
    counts_tb.columns = ['contigID'] + taxa_level + list(range(counts_tb.shape[1] - len(taxa_level) -1))
    
    contigs_set = list(set(genes_tb['SEQUENCE']))
    
    contigs_taxa_counts = counts_tb.loc[counts_tb['contigID'].isin(contigs_set)]
    contigs_taxa_counts['group1'] = contigs_taxa_counts[group1_index].sum(axis = 1)
    contigs_taxa_counts['group2'] = contigs_taxa_counts[group2_index].sum(axis = 1)
    contigs_taxa_counts.index = contigs_taxa_counts['contigID']
    contigs_taxa_counts = contigs_taxa_counts[taxa_level + ['group1','group2']]
    #contigs_taxa_counts.drop(labels = group1_index + group2_index + ['contigID'], axis = 1, inplace = True)

    # nimA gene is NIMA gene
    # tet(W)_6 gene is tet(W) gene
    # wbaP/rfbP is wbaP and rfbP
    # b(pilb should be deleted
    # tet37 is tet-37 is tet(37)
    def has_numbers(inputString):
        return any(char.isdigit() for char in inputString)
    def clean_genes(genes: list, plasmid):
        genes = [gene.lower().strip("'") for gene in genes]
        genes = list(chain(*[gene.split('/') for gene in genes]))
        genes = [gene.split('_')[0] if gene.count('_') == 1 else gene for gene in genes ]
        genes = [gene for gene in genes  if sum([('(' in gene), (')' in gene)]) in [0,2] ]
        if not plasmid:
            #genes = list(chain(*[gene.split('(') if has_numbers(gene) else [gene] for gene in genes]))
            #genes = list(chain(*[gene.split(')') if has_numbers(gene) else [gene] for gene in genes]))
            #genes = list(chain(*[gene.replace(')','').replace('(','') if len(gene) <= 6 else [gene] for gene in genes]))
            genes = list(chain(*[gene.split('-') if gene.count('_') <= 3 else [gene] for gene in genes]))
            genes = [gene for gene in genes if not gene.isnumeric()]
            genes = [gene for gene in genes if 'iii' not in gene]
            genes = [gene for gene in genes if "''" not in gene]
            genes = [gene.rstrip(string.digits) for gene in genes]
            genes = [gene for gene in genes if len(gene) >= 3]
        return list(set(genes))
        
    
    contigs_genes_dict = {}
    for contig in contigs_set:
        contigs_genes_dict[contig] = clean_genes(genes_tb.loc[genes_tb['SEQUENCE'].isin([contig])]['GENE'], plasmid)
    
    # generate the intermediate table
    tb_dict = {}
    idx = 0
    for contig,genes in contigs_genes_dict.items():
        for gene in genes:
            try:
                tb_dict[idx] = [gene] + list(contigs_taxa_counts.loc[contig])
                print(tb_dict[idx])
                idx += 1
            except: pass
    tb_df = pd.DataFrame(tb_dict).T
    tb_df.columns = ['gene'] + taxa_level + ['group1','group2']
    tb_df['total'] = tb_df[['group1','group2']].sum(axis = 1)
    tb_df.to_csv(path.join(output_db, f'genes_taxa_counts_{group1_name}_vs_{group2_name}.tsv'), sep = '\t', header = True, index = None)
    
    # genes vs taxa heatmap
    for taxa in taxa_level:
        temp_df = tb_df[['gene',taxa,'total']]
        temp_df = temp_df[temp_df[taxa].notna()]
        temp_df_groupby = temp_df.groupby(['gene']).agg({taxa:'first','total':'sum'})
        #temp_df_log10 =temp_df_groupby.apply(np.log)
        #temp_df = temp_df_log10.reset_index()
        # only select top20 genes
        #temp_df_top20 = temp_df.groupby(['gene']).agg({taxa:'first','total':'sum'})
        gene_total_dict = {}
        for gene in temp_df_groupby.index:
            try:
                gene_total_dict[gene] = sum(temp_df_groupby.loc[gene]['total'])
            except:
                gene_total_dict[gene] = sum([temp_df_groupby.loc[gene]['total']])
        if len(gene_total_dict) >= 20:
            gene_total_dict = sorted(gene_total_dict.items(), key=lambda x: x[1], reverse=True)[:20]
        else:
            gene_total_dict = sorted(gene_total_dict.items(), key=lambda x: x[1], reverse=True)
        gene_top20 = [x[0] for x in gene_total_dict]
        temp_df_top20 = temp_df.loc[temp_df['gene'].isin(gene_top20)]
        temp_df_top20.index = temp_df_top20['gene']
        temp_df_top20.drop(labels = ['gene'], inplace = True, axis = 1)
        taxa_list = list(set(temp_df_top20[taxa]))
        # create gene vs taxa heatmap
        gene_taxa_dict = {}
        for gene in gene_top20:
            temp_df = temp_df_top20.loc[gene]
            gene_taxa_dict[gene] = []
            if type(temp_df[taxa]) is not pd.core.series.Series:
                hit_taxa = [temp_df[taxa]]
            else:
                hit_taxa = list(temp_df[taxa])

            for t in taxa_list:
                if t in hit_taxa:
                    if type(temp_df) is pd.core.series.Series:
                        gene_taxa_dict[gene] += [temp_df['total']]
                    else:
                        gene_taxa_dict[gene] += [sum(list(temp_df.loc[temp_df[taxa] == t]['total']))]
                else:
                    gene_taxa_dict[gene] += [0]
        gene_taxa_df = pd.DataFrame.from_dict(gene_taxa_dict, orient = 'index', columns = taxa_list).T
        second_min=heapq.nsmallest(2,set(gene_taxa_df.to_numpy().flatten()))[1]
        gene_taxa_df.replace(0,second_min/100,inplace=True)
        gene_taxa_df = gene_taxa_df.apply(np.log)
        
        try:
            plt.figure(figsize=(gene_taxa_df.shape[1] / 3, gene_taxa_df.shape[0] * 2))
            sns.clustermap(data=gene_taxa_df, xticklabels=True, yticklabels=True)
            plt.tight_layout()
            plt.ioff()
            plt.savefig(path.join(output_db, f'genes_{taxa}_{group1_name}_vs_{group2_name}_heatmap.pdf'))
            plt.close()
        except: plt.close()
    
    
    # genes vs group1/group2 distrution plot
    tb_df_gene = tb_df.groupby(['gene']).agg({'group1':'sum','group2':'sum','total':'sum'})
    tb_df_sort = tb_df_gene.sort_values(by=['total'], ascending = False).reset_index(drop = False)
    if tb_df_sort.shape[0] >= 20:
        tb_df_top20 = tb_df_sort.loc[:19]
    else:
        tb_df_top20 = tb_df_sort
    second_min=heapq.nsmallest(2,set(tb_df_top20[['group1','group2','total']].to_numpy().flatten()))[1]
    tb_df_top20.replace(0,second_min/100,inplace=True)
    tb_df_top20.index = tb_df_top20['gene']
    tb_df_top20.drop(['gene'], axis = 1, inplace = True)
    tb_df_top20 = tb_df_top20.astype(float).apply(np.log10)
    
    x = np.arange(tb_df_top20.shape[0])
    y1 = list(tb_df_top20['group1'])
    y2 = list(tb_df_top20['group2'])
    width = 0.4
    plt.figure(figsize=(tb_df_top20.shape[0] / 4, 5))
    plt.bar(x-0.2, y1, width, align='edge')
    plt.bar(x+0.2, y2, width, align='edge')
    plt.xticks(x, list(tb_df_top20.index), rotation='vertical')
    plt.xlabel("genes")
    plt.ylabel("log10")
    plt.legend([group1_name, group2_name])
    plt.tight_layout()
    plt.ioff()
    plt.savefig(path.join(output_db, f'genes_{group1_name}_vs_{group2_name}_distrplot.pdf'))
    plt.close()


if __name__=='__main__':
    gene_tb_dp = sys.argv[1]
    counts_tb_dp = sys.argv[2]
    output_db = sys.argv[3]
    group1_index = [int(_) for _ in sys.argv[4].split(',')]
    group1_name = sys.argv[5]
    group2_index = [int(_) for _ in sys.argv[6].split(',')]
    group2_name = sys.argv[7]
    taxa_level = sys.argv[8].split(',')
    try:
        plasmid = eval(sys.argv[9])
    except:
        plasmid = False
    """
    #gene_tb_dp = r"D:\Projects\GEMINI\cat\antibiotic_analysis\cat.antibiotic.tsv"
    gene_tb_dp = r"D:\Projects\GEMINI\cat\plasmid_analysis\cat.plasmidfinder.tsv"
    #gene_tb_dp = r"D:\Projects\GEMINI\cat\virulence_analysis\cat.virulence.tsv"
    counts_tb_dp = r"D:\Projects\GEMINI\cat\assembly_analysis\cat.taxa_counts.tsv"
    #output_db = r"D:\Projects\GEMINI\cat\antibiotic_analysis"
    output_db = r"D:\Projects\GEMINI\cat\plasmid_analysis"
    #output_db = r"D:\Projects\GEMINI\cat\virulence_analysis"
    group1_index = [0,1,2,3,4,13,14,15]
    group1_name = 'normal'
    group2_index = [5,6,7,8,9,10,11,12]
    group2_name = 'obese'
    taxa_level = ['taxaID','superkingdom','phylum','class','order','family','genus','species']
    plasmid = True
    """

    try:
        visualize_abricate(gene_tb_dp, counts_tb_dp, output_db, 
                           group1_index, group1_name, group2_index, group2_name, 
                           taxa_level, plasmid)
    except Exception as e:
        print(e)
        print("GEMINI: Error. visualize_abricate: ",gene_tb_dp," has problem.skip.")
