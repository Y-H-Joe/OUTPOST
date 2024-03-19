# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 00:00:16 2023

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

================================ description ==================================
#mpa_vOct22_CHOCOPhlAnSGB_202212		
clade_name	profiled_metagenome	dogd_metagenome
k__Bacteria	100	100
k__Bacteria|p__Bacteroidetes	63.82247	48.29516
k__Bacteria|p__Firmicutes	21.04578	10.91573
k__Bacteria|p__Proteobacteria	14.56321	32.25881
k__Bacteria|p__Actinobacteria	0.56854	0
k__Bacteria|p__Bacteroidetes|c__Bacteroidia	59.19866	40.06818
k__Bacteria|p__Firmicutes|c__Clostridia	16.9041	7.81597
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria	9.28481	30.86117
k__Bacteria|p__Proteobacteria|c__Betaproteobacteria	5.2784	0
k__Bacteria|p__Bacteroidetes|c__Flavobacteriia	3.92409	0
k__Bacteria|p__Firmicutes|c__CFGB75522	2.31427	0
k__Bacteria|p__Firmicutes|c__CFGB49531	1.82741	0
k__Bacteria|p__Bacteroidetes|c__CFGB76132	0.69972	4.22433
k__Bacteria|p__Actinobacteria|c__Actinomycetia	0.56854	0
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales	59.19866	40.06818
k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales	16.9041	7.81597

=================================== input =====================================

=================================== output ====================================
	Bacteroidetes	Firmicutes	Proteobacteria	Actinobacteria	Fusobacteria
profiled_metagenome	63.82247	21.04578	14.56321	0.56854	0.0
dogd_metagenome	48.29516	10.91573	32.25881	0.0	8.5303

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================
"""
import pandas as pd
import sys
import os

def taxa_annotation2tsv(taxa_annotation, sample_list, output_dir):
    # 读取数据
    df = pd.read_csv(taxa_annotation, sep='\t', comment='#')
    df.columns[1:] = sample_list
    # 分类等级
    levels = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    real_levels = ['superkingdom','phylum','class','order','family','genus','species']
    # 存储每个等级的DataFrame
    taxa_dfs = {level: pd.DataFrame(index = sample_list) for level in levels}
    
    # 遍历原始数据
    for index, row in df.iterrows():
        clade_names = row['clade_name'].split('|')
        
        # 对每个级别进行操作
        for level in levels:
            # 查找与当前级别匹配的分类名称
            assigned_taxo =  clade_names[-1]
            if assigned_taxo.startswith(f"{level}__"):
                taxonomy = assigned_taxo.split('__')[1]  # 提取当前级别的分类名称
                taxa_dfs[level][taxonomy] = row[1:]
                break
    # 保存为TSV文件
    for level, df_level in taxa_dfs.items():
        df_level.to_csv(os.path.join(output_dir, f"taxonomy_{real_levels[levels.index(level)]}.tsv"), sep='\t')

if __name__ == '__main__':
    taxa_annotation = sys.argv[1]
    sample_list = sys.argv[2].split(',')
    output_dir = sys.argv[3]
    # TODO type = metaphly4
    taxa_annotation2tsv(taxa_annotation, sample_list, output_dir)

