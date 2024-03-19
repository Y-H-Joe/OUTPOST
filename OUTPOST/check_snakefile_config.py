# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 13:30:52 2024

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
import yaml
import sys
import os
from subprocess import Popen, PIPE

try:
    snakefile_config = sys.argv[1]
    with open(snakefile_config, 'r') as f:
        config = yaml.safe_load(f)
except:
    script_dir = os.path.dirname(__file__)  # 获取脚本所在目录
    snakefile_config = os.path.join(script_dir, 'Snakefile_config.yml')
    with open(snakefile_config, 'r') as f:
        config = yaml.safe_load(f)
        
error = 0

def check_command(command):
    """检查命令是否可执行"""
    global error
    try:
        process = Popen(['which', command], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f"OUTPOST error: command '{command}' cannot be executed. Need to fix.")
            error += 1
        return command
    except Exception as e:
        print(e)
        print(f"OUTPOST error: {command} has error.")
        error += 1

def check_path(file_path):
    """检查文件是否存在"""
    global error
    try:
        if not os.path.exists(file_path):
            print(f"OUTPOST error: path '{file_path}' does not exist. Need to fix.")
            error += 1
        return file_path
    except Exception as e:
        print(e)
        print(f"OUTPOST error: {file_path} has error.")
        error += 1

def check_int(parameter):
    global error
    try:
        if isinstance(parameter, int) and( parameter > 0):
            return parameter
        else:
            print(f"OUTPOST error: parameter '{parameter}' is not integer. Need to fix.")
            error += 1
    except Exception as e:
        print(e)
        print(f"OUTPOST error: {parameter} has error.")
        error += 1

def check_list(parameter, list_):
    """检查parameter是否在list_中。如果parameter是列表，则检查每个元素是否在list_中。"""
    global error
    try:
        if isinstance(parameter, list):
            # 检查parameter列表中的每个元素是否都在list_中
            if all(item in list_ for item in parameter):
                return parameter
            else:
                print(f"OUTPOST error: parameter '{parameter}' is not in list {list_}. Need to fix.")
                error += 1
        elif parameter in list_:
            return parameter
        else:
            print(f"OUTPOST error: parameter '{parameter}' is not in list {list_}. Need to fix.")
            error += 1
    except Exception as e:
        print(e)
        print(f"OUTPOST error: {parameter} has error.")
        error += 1
            
            
def check_range(value, range_):
    global error
    try:
        """检查value是否在range中"""
        if min(range_) <= value <= max(range_):
            return value
        else:
            print(f"OUTPOST error: value '{value}' is not in range {range_}. Need to fix.")
            error += 1
    except Exception as e:
        print(e)
        print(f"OUTPOST error: {value} has error.")
        error += 1

# %% command parameters
abricate                      = check_command(config['abricate'])
adapters                      = check_path(config['adapters'])
bamToFastq                    = check_command(config['bamToFastq'])
bwa                           = check_command(config['bwa'])
cdhit                         = check_command(config['cdhit'])
cdhit_cutoff                  = check_range(config['cdhit_cutoff'],[0,1])
diamond                       = check_command(config['diamond'])
emapper                       = check_command(config['emapper'])
emapper_db                    = check_path(config['emapper_db'])
export2graphlan               = check_path(config['export2graphlan'])
fastp                         = check_command(config['fastp'])
graphlan                      = check_command(config['graphlan'])
graphlan_annotate             = check_command(config['graphlan_annotate'])
humann                        = check_command(config['humann'])
humann_join_tables            = check_command(config['humann_join_tables'])
humann_regroup_table          = check_command(config['humann_regroup_table'])
humann_rename_table           = check_command(config['humann_rename_table'])
humann_renorm_table           = check_command(config['humann_renorm_table'])
humann_split_stratified_table = check_command(config['humann_split_stratified_table'])
kaiju                         = check_command(config['kaiju'])
kaiju_addTaxonNames           = check_command(config['kaiju_addTaxonNames'])
kaiju_fmi                     = check_path(config['kaiju_fmi'])
kaiju_names                   = check_path(config['kaiju_names'])
kaiju_nodes                   = check_path(config['kaiju_nodes'])
ktimporttext                  = check_command(config['ktimporttext'])
lefse_format_input            = check_command(config['lefse_format_input'])
lefse_run                     = check_command(config['lefse_run'])
megahit                       = check_command(config['megahit'])
merge_metaphlan_tables        = check_command(config['merge_metaphlan_tables'])
metaphlan                     = check_command(config['metaphlan'])
metaspades                    = check_command(config['metaspades'])
mgm                           = check_command(config['mgm'])
mod_file                      = check_path(config['mod_file'])
OUTPOST_config                = check_path(config["config"])
prodigal                      = check_command(config['prodigal'])
python3                       = check_command(config['python3'])
quast                         = check_command(config['quast'])
rgi                           = check_command(config['rgi'])
Rscript                       = check_command(config['Rscript'])
salmon                        = check_command(config['salmon'])
samtools                      = check_command(config['samtools'])
seqkit                        = check_command(config['seqkit'])
trimmomatic                   = check_command(config['trimmomatic'])
virus_genome                  = check_path(config['virus_genome'])

# %% settings
assemble_contigs              = check_list(config['assemble_contigs'],[True, False]) 
assembly_method               = check_list(config['assembly_method'], ['megahit', 'metaspades']) 
biomarker_num                 = check_int(config['biomarker_num'])
clean_unnecessary             = check_list(config['clean_unnecessary'], [True, False])
cores                         = check_int(config['cores'])
databases                     = check_list(config['databases'], ['rxn','eggnog','ko','level4ec','pfam']) 
downsample_reads              = check_int(int(config['downsample_reads']) - (int(config['downsample_reads']) % 4))
LDA_cutoff                    = check_range(config['LDA_cutoff'], [0,99])
memory_use                    = check_list(config['memory_use'],['minimum', 'maximum'])
output_dir                    = check_path(config['output_dir'])
paired                        = check_list(config['paired'],[True, False])
process_batch_size            = check_int(config['process_batch_size'])
qvalue                        = check_list(config["qvalue"], ['Bonferroni', 'Bonferroni-Holm', 'Benjamini-Hochberg']) 
rm_batch_effect               = check_list(config['rm_batch_effect'],[True, False]) 
skip_assembly_qtest           = check_list(config['skip_assembly_qtest'],[True, False])
skip_humann_init              = check_list(config['skip_humann_init'],[True, False]) 
skip_kaiju                    = check_list(config['skip_kaiju'],[True, False])
taxa_level                    = check_list(config['taxa_level'],['taxaID','superkingdom','phylum','class','order','family','genus','species'])
top                           = check_int(config['top'])
two_sided                     = check_list(config['two_sided'],[True, False])

if error > 0:
    sys.exit()
print(f"Snakefile {snakefile_config} is good! You now can run your OUTPOST as long as your OUTPOST configuratin file is set.")
