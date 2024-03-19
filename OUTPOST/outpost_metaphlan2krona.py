# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 18:46:55 2024

@author: Yihang (Ethan) Zhou


Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

================================ description ==================================

=================================== input =====================================
#mpa_vOct22_CHOCOPhlAnSGB_202212
#/home/z/zhouyh/anaconda3/envs/OUTPOST/bin/metaphlan /home/z/zhouyh/softwares/OUTPOST/results/data/dog_1_nondogvirus.fq --offline --bowtie2out results/metaphlan_analysis/bowtie2/metaphlan_dog_1.bowtie2.bz2 --nproc 2 --input_type fastq -o results/metaphlan_analysis/taxonomy/metaphlan_dog_1_taxa.txt
#14537 reads processed
#SampleID	Metaphlan_Analysis
#clade_name	NCBI_tax_id	relative_abundance	additional_species
k__Bacteria	2	100.0	
k__Bacteria|p__Proteobacteria	2|1224	75.96972	
k__Bacteria|p__Bacteroidetes	2|976	24.03028	
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria	2|1224|1236	75.96972	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia	2|976|200643	24.03028	
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales	2|1224|1236|91347	75.96972	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales	2|976|200643|171549	24.03028	
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae	2|1224|1236|91347|543	75.96972	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae	2|976|200643|171549|171552	24.03028	
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__GGB1354	2|1224|1236|91347|543|	75.96972	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella	2|976|200643|171549|171552|838	24.03028	
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__GGB1354|s__GGB1354_SGB1818	2|1224|1236|91347|543||	75.96972	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A	2|976|200643|171549|171552|838|165179	24.03028	
k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__GGB1354|s__GGB1354_SGB1818|t__SGB1818	2|1224|1236|91347|543|||	75.96972	
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri_clade_A|t__SGB1626	2|976|200643|171549|171552|838|165179|	24.03028	
=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================

"""
#!/usr/bin/env python

# ============================================================================== 
# Conversion script: from MetaPhlAn output to Krona text input file
# Author: Daniel Brami (daniel.brami@gmail.com)
# ==============================================================================

import sys
import optparse
import re

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-p', '--profile', dest='profile', default='', action='store', help='The input file is the MetaPhlAn standard result file' )
    parser.add_option( '-k', '--krona', dest='krona', default='krona.out', action='store', help='the Krons output file name' )
    ( options, spillover ) = parser.parse_args()

    if not options.profile or not options.krona:
        parser.print_help()
        sys.exit()

    re_candidates = re.compile(r"s__")
    re_replace = re.compile(r"\w__")
    re_bar = re.compile(r"\|")

    metaPhLan = list()
    with open(options.profile,'r') as f:
        metaPhLan = f.readlines()
    f.close()

    krona_tmp = options.krona 
    metaPhLan_FH = open(krona_tmp, 'w')

    for aline in (metaPhLan):
        if(re.search(re_candidates, aline)):
            x=re.sub(re_replace, '\t', aline)
            x=re.sub(re_bar, '', x)
            x_cells = x.split('\t')
            lineage = '\t'.join(x_cells[:-3])# remove the NCBI ID column
            try:
                abundance = float(x_cells[-2].rstrip('\n')) 
                metaPhLan_FH.write('%s\n'%(str(abundance) + '\t' + lineage))
            except Exception as e:
                pass

    metaPhLan_FH.close()

if __name__ == '__main__':
    main()