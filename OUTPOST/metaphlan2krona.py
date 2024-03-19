# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 17:59:04 2023

@author: Yihang (Ethan) Zhou


Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

================================ description ==================================

=================================== input =====================================
#mpa_vOct22_CHOCOPhlAnSGB_202212
#/home/z/zhouyh/anaconda3/envs/OUTPOST/bin/metaphlan /home/z/zhouyh/softwares/OUTPOST/results/data/cat_1_noncatvirus.fq --offline --bowtie2out results/metaphlan_analysis/bowtie2/metaphlan_cat_1.bowtie2.bz2 --nproc 2 --input_type fastq -o results/metaphlan_analysis/taxonomy/metaphlan_cat_1_taxa.txt
#16218 reads processed
#SampleID	Metaphlan_Analysis
#clade_name	NCBI_tax_id	relative_abundance	additional_species
k__Bacteria	2	100.0	
k__Bacteria|p__Bacteroidetes	2|976	100.0	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia	2|976|200643	100.0	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales	2|976|200643|171549	100.0	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae	2|976|200643|171549|171552	100.0	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella	2|976|200643|171549|171552|838	100.0	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A	2|976|200643|171549|171552|838|165179	100.0	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626	2|976|200643|171549|171552|838|165179|	100.0	


=================================== output ====================================
2	Fats	Saturated fat
3	Fats	Unsaturated fat	Monounsaturated fat
3	Fats	Unsaturated fat	Polyunsaturated fat
13	Carbohydrates	Sugars
4	Carbohydrates	Dietary fiber
21	Carbohydrates
5	Protein
================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================

"""
import pandas as pd
import os
import sys

def metaphlan2krona(dp, output):
    df = pd.read_csv(dp, header = None, sep = '\t', comment = '#')
    df = df.dropna(how='all')
    # 分类等级
    levels = ['s', 'g', 'f', 'o', 'c', 'p', 'k']
    # real_levels = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
    # 遍历原始数据
    krona_df = []
    for level in levels:
        hit = 0
        for index, row in df.iterrows():
            clade_names = row[0].split('|')
            # 对每个级别进行操作
            # 查找与当前级别匹配的分类名称
            assigned_taxo =  clade_names[-1]
            if assigned_taxo.startswith(f"{level}__"):
                hit += 1
                krona_df.append([row[2]] + [x.split('__')[-1] for x in clade_names])
        if hit > 0:
            break
    # 保存为TSV文件
    krona_df = pd.DataFrame(krona_df)
    krona_df.to_csv( output, sep='\t', header=None, index=None)

if __name__ == '__main__':
    metaphlan = sys.argv[1]
    output = sys.argv[2]
    metaphlan2krona(metaphlan, output)
    

