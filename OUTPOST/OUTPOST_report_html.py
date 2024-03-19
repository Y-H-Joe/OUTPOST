# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:03:14 2024

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

================================ description ==================================

=================================== input =====================================

=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================
"""
import sys
import os
import glob
import shutil
import traceback
from datetime import datetime
import subprocess
from jinja2 import Template

results_dir = sys.argv[1]
snakefile_config = sys.argv[2]
outpost_config = sys.argv[3]
python3 = sys.argv[4]
html_template = sys.argv[5]

assert os.path.exists(results_dir), "error, {results_dir} does not exist!"
assert os.path.exists(snakefile_config), "error, {snakefile_config} does not exist!"
assert os.path.exists(outpost_config), "error, {outpost_config} does not exist!"

report_dir = os.path.join(results_dir, 'report')
materials_dir = os.path.join(report_dir, 'materials')
html_template_basename = os.path.basename(html_template)
report_1 = os.path.join('materials', '1_configuration')
report_2 = os.path.join('materials', '2_quality_control')
report_3 = os.path.join('materials', '3_taxonomy')
report_4 = os.path.join('materials', '4_diversity')
report_5 = os.path.join('materials', '5_metaphlan')
report_6 = os.path.join('materials', '6_function')
report_7 = os.path.join('materials', '7_lda')
report_8 = os.path.join('materials', '8_antibiotic')
report_9 = os.path.join('materials', '9_virulence_factors')
report_10 = os.path.join('materials', '10_plasmids')
report_11 = os.path.join('materials', '11_assembly')
report_12 = os.path.join('materials', '12_biomarkers')

report_1_abs = os.path.join(materials_dir, '1_configuration')
report_2_abs = os.path.join(materials_dir, '2_quality_control')
report_3_abs = os.path.join(materials_dir, '3_taxonomy')
report_4_abs = os.path.join(materials_dir, '4_diversity')
report_5_abs = os.path.join(materials_dir, '5_metaphlan')
report_6_abs = os.path.join(materials_dir, '6_function')
report_7_abs = os.path.join(materials_dir, '7_lda')
report_8_abs = os.path.join(materials_dir, '8_antibiotic')
report_9_abs = os.path.join(materials_dir, '9_virulence_factors')
report_10_abs = os.path.join(materials_dir, '10_plasmids')
report_11_abs = os.path.join(materials_dir, '11_assembly')
report_12_abs = os.path.join(materials_dir, '12_biomarkers')

report_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

# clean
try:
    shutil.rmtree(report_dir)
except Exception as e:
    print(e)
for _ in [report_1_abs, report_2_abs, report_3_abs, report_4_abs, report_5_abs, 
         report_6_abs, report_7_abs, report_8_abs, report_9_abs, report_10_abs, 
         report_11_abs, report_12_abs]:
    os.makedirs(_, exist_ok= True)

#%% template
report_dict = {
    'c_1': os.path.join(results_dir, 'report','materials', '1_configuration'),
    'c_2': os.path.join(results_dir, 'report','materials', '2_quality_control'),
    'c_3': os.path.join(results_dir, 'report','materials', '3_taxonomy'),
    'c_4': os.path.join(results_dir, 'report','materials', '4_diversity'),
    'c_5': os.path.join(results_dir, 'report','materials', '5_metaphlan'),
    'c_6': os.path.join(results_dir, 'report','materials', '6_function'),
    'c_7': os.path.join(results_dir, 'report','materials', '7_lda'),
    'c_8': os.path.join(results_dir, 'report','materials', '8_antibiotic'),
    'c_9': os.path.join(results_dir, 'report','materials', '9_virulence_factors'),
    'c_10': os.path.join(results_dir, 'report','materials', '10_plasmids'),
    'c_11': os.path.join(results_dir, 'report','materials', '11_assembly'),
    'c_12': os.path.join(results_dir, 'report','materials', '12_biomarkers'),
    }

source_dict = {
    'f_11':'',
    'f_12':outpost_config,
    'f_13':os.path.join(results_dir, 'log'),
    'f_14':os.path.join(results_dir, 'benchmark'),
    'f_211':'',
    'f_22':'',
    'f_212': os.path.join(results_dir, 'qc'),
    'f_311':'',
    'f_312':'',
    'f_321':'',
    'f_322':'',
    'f_331':os.path.join(results_dir, 'taxonomy_analysis', 'temp'),
    'f_332':os.path.join(results_dir, 'taxonomy_analysis', 'kaiju'),
    'f_333':os.path.join(results_dir, 'taxonomy_analysis', 'counts_tables'),
    'f_334':os.path.join(results_dir, 'taxonomy_analysis', 'top_taxa_*'),
    'f_335':os.path.join(results_dir, 'taxonomy_analysis', 'utest_*'),
    'f_41':os.path.join(results_dir, 'taxonomy_analysis', 'alpha_beta_*_vs_*'),
    'f_411':'',
    'f_412':os.path.join(results_dir, 'taxonomy_analysis', 'alpha_beta_*_vs_*'),
    'f_421':os.path.join(results_dir, 'taxonomy_analysis', 'alpha_beta_*_vs_*'),
    'f_422':os.path.join(results_dir, 'taxonomy_analysis', 'alpha_beta_*_vs_*'),
    'f_51':os.path.join(results_dir, 'metaphlan_analysis', 'figs'),
    'f_511':'',
    'f_512':'',
    'f_513':'',
    'f_52':os.path.join(results_dir, 'metaphlan_analysis', 'diversity'),
    'f_531':os.path.join(results_dir, 'metaphlan_analysis', 'taxonomy'),
    'f_532':os.path.join(results_dir, 'metaphlan_analysis', 'krona'),
    'f_54':os.path.join(results_dir, 'metaphlan_analysis', 'bowtie2'),
    'f_611':'',
    'f_612':'',
    'f_621':os.path.join(results_dir, 'function_analysis', 'humann3', 'ori_results'),
    'f_622':os.path.join(results_dir, 'function_analysis', 'humann3', 'top_humann_*_vs_*'),
    'f_623':os.path.join(results_dir, 'function_analysis', 'humann3', 'utest_*_vs_*'),
    'f_711':'',
    'f_712':os.path.join(results_dir, 'LDA_analysis', 'taxa_*_vs_*'),
    'f_721':'',
    'f_722':os.path.join(results_dir, 'LDA_analysis', 'humann_*_vs_*'),
    'f_81':'',
    'f_82':'',
    'f_831':os.path.join(results_dir, 'antibiotic_genes_analysis', 'assembly_based', 'genes_taxa_counts_*_vs_*.tsv'),
    'f_832':os.path.join(results_dir, 'antibiotic_genes_analysis', 'assembly_based', '*.antibiotic.tsv'),
    'f_91':'',
    'f_92':'',
    'f_931':os.path.join(results_dir, 'virulence_factors_analysis', 'assembly_based', 'genes_taxa_counts_*_vs_*.tsv'),
    'f_932':os.path.join(results_dir, 'virulence_factors_analysis', 'assembly_based', '*.virulence.tsv'),
    'f_101':'',
    'f_102':'',
    'f_1031':os.path.join(results_dir, 'plasmids_analysis', 'assembly_based', 'genes_taxa_counts_*_vs_*.tsv'),
    'f_1032':os.path.join(results_dir, 'plasmids_analysis', 'assembly_based', '*.plasmidfinder.tsv'),
    'f_111':os.path.join(results_dir, 'assembly_analysis', 'contigs'),
    'f_112':os.path.join(results_dir, 'assembly_analysis', 'outpost_contigs'),
    'f_113':'',
    'f_114':os.path.join(results_dir, 'assembly_analysis', 'annotation'),
    'f_1141':os.path.join(results_dir, 'assembly_analysis', 'annotation', 'metagenemark'),
    'f_1142':os.path.join(results_dir, 'assembly_analysis', 'annotation', 'prodigal'),
    'f_115':os.path.join(results_dir, 'assembly_analysis', 'annotation'),
    'f_1151':os.path.join(results_dir, 'assembly_analysis', 'annotation', 'eggmapper'),
    'f_1152':os.path.join(results_dir, 'assembly_analysis', 'annotation', 'rgi'),
    'f_1153':os.path.join(results_dir, 'assembly_analysis', 'annotation', 'gtdbtk'),
    'f_116':os.path.join(results_dir, 'assembly_analysis', 'quantify'),
    'f_117':os.path.join(results_dir, 'assembly_analysis', 'taxa_counts'),
    'f_121':os.path.join(results_dir, 'biomarkers_analysis', 'OUTPOST_biomarker_scores*'),
    'f_122':os.path.join(results_dir, 'biomarkers_analysis', 'ANCOM_identification'),
    'f_1221':'',
    'f_1222':'',
    'f_123':os.path.join(results_dir, 'biomarkers_analysis', 'ANCOM_identification'),
    }
target_dict = {
    'f_11':'',
    'f_211':'',
    'f_22':'',
    'f_311':'',
    'f_312':'',
    'f_321':'',
    'f_322':'',
    'f_411':'',
    'f_421':'',
    'f_511':'',
    'f_512':'',
    'f_513':'',
    'f_611':'',
    'f_612':'',
    'f_711':'',
    'f_721':'',
    'f_81':'',
    'f_82':'',
    'f_91':'',
    'f_92':'',
    'f_101':'',
    'f_102':'',
    'f_113':'',
    'f_1221':'',
    'f_1222':''
    }


#%% put materials to folder
def find_and_copy(pattern, original_dir, target_dir):
    # 使用os.walk遍历original_dir及其所有子目录
    for root, dirs, files in os.walk(original_dir):
        for name in files:
            # 检查文件名是否匹配给定的模式
            if glob.fnmatch.fnmatch(name, pattern):
                # 构建完整的文件路径
                file_path = os.path.join(root, name)
                # 确保目标目录存在
                if not os.path.exists(target_dir):
                    try:
                        os.makedirs(target_dir)
                    except OSError as error:
                        print(f"Error creating directory {target_dir}: {error}")
                        return ''
                # 尝试复制找到的文件到目标目录
                try:
                    shutil.copy(file_path, target_dir)
                    print(f"File {file_path} copied to {target_dir}")
                    return os.path.join(target_dir, name)
                except Exception as e:
                    print(f"Error copying file: {e}")
                    traceback.print_exc()
                    return ''
    # 如果遍历完毕都没有找到匹配的文件
    print("No files matched the given pattern.")
    return ''


def run_command(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stderr:
        print("Error:", result.stderr)

#%%%% copy files from results to materials
# 1
shutil.copy(snakefile_config, report_1_abs)
source_dict['f_11'] = snakefile_config
target_dict['f_11'] = os.path.join(report_1, os.path.basename(snakefile_config))

# 2
target_dp1 = find_and_copy("*.fastp.html", os.path.join(results_dir, 'qc'), report_2_abs)
source_dict['f_211'] = os.path.join(results_dir, 'qc', "*.fastp.html")
target_dict['f_211'] = os.path.join(report_2, os.path.basename(target_dp1))
target_dp2 = find_and_copy("all_samples.taxa_counts.rel_abun.*.rmU*pdf", os.path.join(results_dir, 'batch_effect'), report_2_abs)
source_dict['f_22'] = os.path.join(results_dir, 'batch_effect', "all_samples.taxa_counts.rel_abun.*.rmU*pdf")
target_dict['f_22'] = os.path.join(report_2, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_2_abs])

# 3
target_dp1 = find_and_copy("all_samples.rel_abun.*_vs_*.at_*.rel_abun.*.*.boxplot.pdf", os.path.join(results_dir, 'taxonomy_analysis'), report_3_abs)
source_dict['f_311'] = os.path.join(results_dir, 'taxonomy_analysis', "all_samples.rel_abun.*_vs_*.at_*.rel_abun.*.*.boxplot.pdf")
target_dict['f_311'] = os.path.join(report_3, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("all_samples.rel_abun.*_vs_*.at_*.rel_abun.*.top*.fillmin.scaled.heatmap.pdf", os.path.join(results_dir, 'taxonomy_analysis'), report_3_abs)
source_dict['f_312'] = os.path.join(results_dir, 'taxonomy_analysis', "all_samples.rel_abun.*_vs_*.at_*.rel_abun.*.top*.fillmin.scaled.heatmap.pdf")
target_dict['f_312'] = os.path.join(report_3, os.path.basename(target_dp2.replace('.pdf','.png')))
target_dp3 = find_and_copy("all_samples.taxa_counts.rel_abun.*.rmU.top*.barplot.pdf", os.path.join(results_dir, 'taxonomy_analysis'), report_3_abs)
source_dict['f_321'] = os.path.join(results_dir, 'taxonomy_analysis', "all_samples.taxa_counts.rel_abun.*.rmU.top*.barplot.pdf")
target_dict['f_321'] = os.path.join(report_3, os.path.basename(target_dp3.replace('.pdf','.png')))
target_dp4 = find_and_copy("all_samples.taxa_counts.rel_abun.*.rmU.top*.fillmin.scaled.heatmap.pdf", os.path.join(results_dir, 'taxonomy_analysis'), report_3_abs)
source_dict['f_322'] = os.path.join(results_dir, 'taxonomy_analysis', "all_samples.taxa_counts.rel_abun.*.rmU.top*.fillmin.scaled.heatmap.pdf")
target_dict['f_322'] = os.path.join(report_3, os.path.basename(target_dp4.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_3_abs])

#4
target_dp1 = find_and_copy("*alpha_diversity.at_*.pdf", os.path.join(results_dir, 'diversity_analysis'), report_4_abs)
source_dict['f_411'] = os.path.join(results_dir, 'diversity_analysis', "*alpha_diversity.at_*.pdf")
target_dict['f_411'] = os.path.join(report_4, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("PCoA*.*.at_*.pdf", os.path.join(results_dir, 'diversity_analysis'), report_4_abs)
source_dict['f_421'] = os.path.join(results_dir, 'diversity_analysis', "PCoA*.*.at_*.pdf")
target_dict['f_421'] = os.path.join(report_4, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_4_abs])

#5
target_dp1 = find_and_copy("graphlan.pdf", os.path.join(results_dir, 'metaphlan_analysis', 'figs'), report_5_abs)
source_dict['f_511'] = os.path.join(results_dir, 'metaphlan_analysis', 'figs', "graphlan.pdf")
target_dict['f_511'] = os.path.join(report_5, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("metaphlan_merge_taxa.heatmap.pdf",  os.path.join(results_dir, 'metaphlan_analysis', 'figs'), report_5_abs)
source_dict['f_512'] = os.path.join(results_dir, 'metaphlan_analysis', 'figs', "metaphlan_merge_taxa.heatmap.pdf")
target_dict['f_512'] = os.path.join(report_5, os.path.basename(target_dp2.replace('.pdf','.png')))
target_dp3 = find_and_copy("krona.html",  os.path.join(results_dir, 'metaphlan_analysis', 'figs'), report_5_abs)
source_dict['f_513'] = os.path.join(results_dir, 'metaphlan_analysis', 'figs', "krona.html")
target_dict['f_513'] = os.path.join(report_5, os.path.basename(target_dp3))
run_command([python3, "OUTPOST/pdf2png.py", report_5_abs])

#6
target_dp1 = find_and_copy("allSamples_genefamilies_uniref90names_relab_*_unstratified.named.rel_abun_format.top*.fillmin.scaled.heatmap.pdf", os.path.join(results_dir, 'function_analysis', 'figs'), report_6_abs)
source_dict['f_611'] = os.path.join(results_dir, 'function_analysis', 'figs', "allSamples_genefamilies_uniref90names_relab_*_unstratified.named.rel_abun_format.top*.fillmin.scaled.heatmap.pdf")
target_dict['f_611'] = os.path.join(report_6, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("allSamples_genefamilies_uniref90names_relab_*_unstratified.named.rel_abun_format.*_vs_*.rel_abun.*.top*.fillmin.scaled.heatmap.pdf", os.path.join(results_dir, 'function_analysis', 'figs'), report_6_abs)
source_dict['f_612'] = os.path.join(results_dir, 'function_analysis', 'figs', "allSamples_genefamilies_uniref90names_relab_*_unstratified.named.rel_abun_format.*_vs_*.rel_abun.*.top*.fillmin.scaled.heatmap.pdf")
target_dict['f_612'] = os.path.join(report_6, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_6_abs])

#7
target_dp1 = find_and_copy("all_samples.rel_abun.*_vs_*.at_*.rel_abun.*.lefse.pdf", os.path.join(results_dir, 'LDA_analysis', 'figs_taxa'), report_7_abs)
source_dict['f_711'] = os.path.join(results_dir, 'LDA_analysis', 'figs_taxa', "all_samples.rel_abun.*_vs_*.at_*.rel_abun.*.lefse.pdf")
target_dict['f_711'] = os.path.join(report_7, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("allSamples_genefamilies_uniref90names_relab_*unstratified.named.rel_abun_format.*_vs_*.rel_abun.*.lefse.pdf", os.path.join(results_dir, 'LDA_analysis', 'figs_humann'), report_7_abs)
source_dict['f_721'] = os.path.join(results_dir, 'LDA_analysis', 'figs_humann', "allSamples_genefamilies_uniref90names_relab_*unstratified.named.rel_abun_format.*_vs_*.rel_abun.*.lefse.pdf")
target_dict['f_721'] = os.path.join(report_7, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_7_abs])

#8
target_dp1 = find_and_copy("genes_*_*_vs_*_heatmap.pdf", os.path.join(results_dir, 'antibiotic_genes_analysis', 'assembly_based'), report_8_abs)
source_dict['f_81'] = os.path.join(results_dir, 'antibiotic_genes_analysis', 'assembly_based', "genes_*_*_vs_*_heatmap.pdf")
target_dict['f_81'] = os.path.join(report_8, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("genes_*_vs_*_distrplot.pdf", os.path.join(results_dir, 'antibiotic_genes_analysis', 'assembly_based'), report_8_abs)
source_dict['f_82'] = os.path.join(results_dir, 'antibiotic_genes_analysis', 'assembly_based', "genes_*_vs_*_distrplot.png")
target_dict['f_82'] = os.path.join(report_8, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_8_abs])


#9
target_dp1 = find_and_copy("genes_*_*_vs_*_heatmap.pdf", os.path.join(results_dir, 'virulence_factors_analysis', 'assembly_based'), report_9_abs)
source_dict['f_91'] = os.path.join(results_dir, 'virulence_factors_analysis', 'assembly_based', "genes_*_*_vs_*_heatmap.pdf")
target_dict['f_91'] = os.path.join(report_9, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("genes_*_vs_*_distrplot.pdf", os.path.join(results_dir, 'virulence_factors_analysis', 'assembly_based'), report_9_abs)
source_dict['f_92'] = os.path.join(results_dir, 'virulence_factors_analysis', 'assembly_based', "genes_*_vs_*_distrplot.png")
target_dict['f_92'] = os.path.join(report_9, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_9_abs])

#10
target_dp1 = find_and_copy("genes_*_*_vs_*_heatmap.pdf", os.path.join(results_dir, 'plasmids_analysis', 'assembly_based'), report_10_abs)
source_dict['f_101'] = os.path.join(results_dir, 'plasmids_analysis', 'assembly_based', "genes_*_*_vs_*_heatmap.pdf")
target_dict['f_101'] = os.path.join(report_10, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("genes_*_vs_*_distrplot.pdf", os.path.join(results_dir, 'plasmids_analysis', 'assembly_based'), report_10_abs)
source_dict['f_102'] = os.path.join(results_dir, 'plasmids_analysis', 'assembly_based', "genes_*_vs_*_distrplot.png")
target_dict['f_102'] = os.path.join(report_10, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_10_abs])

#11
os.system(f"cp -r {os.path.join(results_dir, 'assembly_analysis', 'quast')} {report_11_abs}")
source_dict['f_113'] = os.path.join(results_dir, 'assembly_analysis', 'quast')
target_dict['f_113'] = os.path.join(report_11, 'quast', 'report.html')

#12
target_dp1 = find_and_copy("ancom_biomarkers.dotplot.*_vs_*.at_*.pdf", os.path.join(results_dir, 'biomarkers_analysis', 'ANCOM_identification'), report_12_abs)
source_dict['f_1221'] = os.path.join(results_dir, 'biomarkers_analysis', 'ANCOM_identification', "ancom_biomarkers.dotplot.*_vs_*.at_*.pdf")
target_dict['f_1221'] = os.path.join(report_12, os.path.basename(target_dp1.replace('.pdf','.png')))
target_dp2 = find_and_copy("ancom_biomarkers.volcano.*_vs_*.at_*.pdf", os.path.join(results_dir, 'biomarkers_analysis', 'ANCOM_identification'), report_12_abs)
source_dict['f_1222'] = os.path.join(results_dir, 'biomarkers_analysis', 'ANCOM_identification', "ancom_biomarkers.volcano.*_vs_*.at_*.pdf")
target_dict['f_1222'] = os.path.join(report_12, os.path.basename(target_dp2.replace('.pdf','.png')))
run_command([python3, "OUTPOST/pdf2png.py", report_12_abs])

#%% jinja2 render
# 直接从文件读取模板内容
with open(html_template, 'r', encoding='utf-8') as file:
    template_string = file.read()
# 创建Template实例
jinja_template = Template(template_string)

# 准备字典
jinjia_dict = {
    'report_time': report_time,
    'report_dict': report_dict,
    'source_dict': source_dict,
    'target_dict': target_dict,
    'results_dir': results_dir,
    }

# print("#################")
# print(jinjia_dict)
# print("#################")
# 渲染
rendered_html = jinja_template.render(jinjia_dict)
# 保存渲染后的HTML到文件
with open(os.path.join(report_dir, "OUTPOST_report.html"), 'w', encoding='utf-8') as f:
    f.write(rendered_html)
