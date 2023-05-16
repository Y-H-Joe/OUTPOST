# GEMINI
An multi-group comparative, flexible, light and fast, downstream pipeline for metagenomics whole genome sequencing analysis.

# GEMINI installation
0. GEMINI only supports Illumina reads, not for Ion Torrent yet.
1. Download `GEMINI.yml`
2. Create the conda environment `conda env create --name GEMINI --file GEMINI.yml `
3. Activate the environment `conda activate GEMINI`
4. Check you're using right humann `which humann`
5. Download 1st humann databases `humann_databases --download chocophlan full /path/to/databases --update-config yes`
6. Download 2nd humann databases `humann_databases --download uniref uniref90_diamond /path/to/databases --update-config yes`
7. Download 3rd humann databases `humann_databases --download utility_mapping full /path/to/databases --update-config yes`
8. Check your humann databases `cd /path/to/databases` then run `ll chocophlan/ | wc -l` you get a number >= 11289. run `ll uniref` you should see a `uniref90_201901b_full.dmnd` (or newer) file with >= 34G size. Run `ll utility_mapping` you should see >= 21 files and all of them have > 3M size (or some of them truncated during download).
9. Check your humann's misc folder. Located at `/path/to/your/anaconda/envs/GEMINI/lib/python3.9/site-packages/humann/data/misc/`. Due to conda's unknown error, usually the files are missing. To know full file list, check README.txt in the misc folder. To download the files, check https://github.com/biobakery/humann/tree/master/humann/data/misc . Or unzip the `misc.zip` in the `utils` folder.
10. Check your humann by running `humann -i sample_reads.fastq -o sample_results` (prepare sample_reads.fastq by yourself)
11. Check you're using right kaiju `which kaiju`
12. Download kaiju databases `mkdir /path/to/kaijudb` then `cd /path/to/kaijudb` then `kaiju-makedb -s nr_euk` (this takes a long time and space and memory). Or you can download and unzip the annotation files from [kaiju servier](https://kaiju.binf.ku.dk/server).
13. Check abricate databases `abricate --list`, you should see 9 databases (argannot,card,ecoh,ecoli_vf,megares,ncbi,plasmidfinder,resfinder,vfdb). If you didn't, go to download abricate [databases](https://github.com/tseemann/abricate/tree/master/db), or use the `db.zip` file in `utils` folder, unzip them under your `db` folder. The location of `db` folder can be seen by running `abricate --help`, see the `--datadir` line. Then `cd` to the `db` folder, run `abricate --setupdb`.
14. Install `iPaper` R package. type `R` in command line, then in the R concle, type in `install.packages("remotes")`, then type in `remotes::install_github("kongdd/Ipaper")`. When R concle asks you whether to update other packages, choose `none`. After installation, type in `library(Ipaper)`, if no error occurs, then you're good to contine.
15. Open `Snakefile.py`, modify the bwa,kaiju,python3,Rscript,...lefse_run parameters to the executable command lines in your environment. To make sure all command line works, please test the command line one by one in your linux shell.
16. Test GEMINI. `cd parent/folder/of/GEMINI`. Prepare some example data. Open and modify the `GEMINI/GEMINI_contig.tsv` to make sure the data_dir is right. then run `snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete`. The test will last 3 hours.
17. If you occured any errors. check the printed log to debug. Or check the log file in `name_of_your_assembly/log` folder. You can use time stamps to refer which rule is error, or to understand the error information. After debugging, delete the `name_of_your_assembly/log/name_of_the_error_rule.done`. and rerun the `snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete`. GEMINI will auto-resume.
18. You should see the outputs in `horsedonkey` folder. If no errors occured. then you're good to go.
19. In a summary, GEMINI itself is just a list of scripts, which is easy to use. The above installation guide actually is helping you to install other tools, such as humann3 and kaiju. On the other hand, if you have installed these tools somewhere else, you can just modify the `Snakemake.py` file to skip the above installation procedure.

# format of GEMINI_config.tsv
| samples                      | fq_dir                                                                      | bam_dir                                                                         | assembly              | assembly_dir                                         | group                | batch |
|------------------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------------|-----------------------|------------------------------------------------------|----------------------|-------|
| old_healthy_male_asian_1     | /home/yh/GEMINI/data/old_healthy_male_asian_1_onlyPE_nonhuman64virus.fq     | /home/yh/GEMINI/bam_sam/old_healthy_male_asian_1_contigs_sorted.human63.bam     | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | asian,healthy,male   | 1     |
| old_healthy_male_asian_2     | /home/yh/GEMINI/data/old_healthy_male_asian_2_onlyPE_nonhuman64virus.fq     | /home/yh/GEMINI/bam_sam/old_healthy_male_asian_2_contigs_sorted.human63.bam     | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | asian,healthy,male   | 1     |
| old_healthy_male_asian_3     | /home/yh/GEMINI/data/old_healthy_male_asian_3_onlyPE_nonhuman64virus.fq     | /home/yh/GEMINI/bam_sam/old_healthy_male_asian_3_contigs_sorted.human63.bam     | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | asian,healthy,male   | 1     |
| old_healthy_male_asian_4     | /home/yh/GEMINI/data/old_healthy_male_asian_4_onlyPE_nonhuman64virus.fq     | /home/yh/GEMINI/bam_sam/old_healthy_male_asian_4_contigs_sorted.human63.bam     | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | asian,healthy,male   | 1     |
| young_ill_female_euro_1      | /home/yh/GEMINI/data/young_ill_female_euro_1_onlyPE_nonhuman64virus.fq      | /home/yh/GEMINI/bam_sam/young_ill_female_euro_1_contigs_sorted.human63.bam      | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | euro,ill,female      | 5     |
| young_ill_female_euro_2      | /home/yh/GEMINI/data/young_ill_female_euro_2_onlyPE_nonhuman64virus.fq      | /home/yh/GEMINI/bam_sam/young_ill_female_euro_2_contigs_sorted.human63.bam      | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | euro,ill,female      | 6     |
| young_ill_female_euro_3      | /home/yh/GEMINI/data/young_ill_female_euro_3_onlyPE_nonhuman64virus.fq      | /home/yh/GEMINI/bam_sam/young_ill_female_euro_3_contigs_sorted.human63.bam      | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | euro,ill,female      | 5     |
| young_ill_female_euro_4      | /home/yh/GEMINI/data/young_ill_female_euro_4_onlyPE_nonhuman64virus.fq      | /home/yh/GEMINI/bam_sam/young_ill_female_euro_4_contigs_sorted.human63.bam      | human62_batch_effect2 | /home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa | euro,ill,female      | 6     |

The ***samples*** column contains the ids of each sample. Each sample id should be unique. The id can contain special characters like '_', but no ',' allowed.

The ***fq_dir*** column contains the absolute location of each fastq file. If you use PE end sequencing, `cat` them to one file, because humann3 asked so. 

The ***bam_dir*** column contains the absolute directory of each ***bam***. The bam file should be generated by users, via aligning reads (both SE and PE mapping is acceptable) to assembly. 

The ***assembly*** column contains the name of the assembly. 

The ***assembly_dir*** column contains the absolute directory of ***assembly***.

The ***group*** column contains the group name of the each sample. If one sample belong to multiple groups, seperate by ','. Each group should at least have 3 samples for statistical power concern.

The ***batch*** is the batch of the data. If `rm_batch_effect` is True, then you should have multiple batches, otherwise set all to be 1 or blank. 

***do not change the name of columns***

In the above case. We have one reference, human62_batch_effect2. For groups, we have 'aisan','euro', 'female', 'male', 'healthy','ill'.

The pipeline will pair-wise compare all of them. Which are, 
['asian_vs_euro','male_vs_female','healthy_vs_ill'...].

# GEMINI usuage
The usuage of GEMINI is very easy and light. If you successfully installed GEMINI and got all the expected outputs files in `human62_batch_effect2` folder (same name as your assembly), you just need to replace the data locations with your own, and modify the `GEMINI/GEMINI_config.tsv`.

1. The working directory must be the parent directory of `GEMINI`, becasue some scripts in Snakemake.py use relative path.
2. Snakemake (the framework GEMINI relied on) will lock working directory during running. So you should prepare two working directories if you're running two GEMINI pipelines (or any other Snakemake based softwares). To be more specific, `mkdir folder1/` and `mkdir folder2/`, copy the entire GEMINI folder to `folder1/` and `folder2/`. Then `cd folder1`, run GEMINI. Then `cd folder2`, run GEMINI.
3. One GEMINI process only take one assembly. If you have multi assemblies to analyze, run GEMINI multiple times (in parallel) (in different working directories).
4. Humann analysis is very time/computation consuption. I personally prefer to use computer cluster to distributedly run Humann. The Snakemake based GEMINI can also be deployed on cluster, but maybe not [easy](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). So GEMINI provide `skip_humann_init` option. Set `skip_humann_init = True` in `GEMINI_config.tsv`, then GEMINI will not run human_init rule, but to check the human results under folder `name_of_the_assembly/metabolism_analysis/humann3/ori_results/`, so you need to put the humann output *genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the folder. If `skip_humann_init = True`, GEMINI will check the outputs existence first then skip the humann step. Make sure these humann outptus are from the same fastq you offered to GEMINI.
5. GEMINI offers `skip_assembly_analysis`. You can skip the assembly analysis if you have a large assembly and a number of groups which will save a lot of time about MAG tables generation. 

# GEMINI outpus
Each GEMINI run accepts one assembly, all outputs are categorized in the folder named by the assembly.
```
(base) yh@superServer:human62_batch_effect2$ l
antibiotic_analysis/
assembly_analysis/
batch_effect/
benchmark/
diversity_analysis/
LDA_analysis/
log/
metabolism_analysis/
plasmid_analysis/
taxa_analysis/
virulence_analysis/
```
## assembly analysis
```
(base) yh@superServer:human62_batch_effect2$ l assembly_analysis/
human62_batch_effect2.taxa_counts.tsv
```
I set `skip_assembly_analysis` to True, so I only have one table. The MAG analysis table can refer to Table S1 of our paper.


## taxa analysis
```
(base) yh@superServer:human62_batch_effect2$ l taxa_analysis/
boxplot_asian_vs_euro/
boxplot_asian_vs_female/
boxplot_asian_vs_healthy/
...
boxplot_healthy_vs_male/
boxplot_male_vs_euro/
boxplot_male_vs_female/
boxplot_male_vs_ill/
figs/ 
human62_batch_effect2.taxa_counts.rel_abun.class.csv
human62_batch_effect2.taxa_counts.rel_abun.class.rmU.csv
human62_batch_effect2.taxa_counts.rel_abun.class.rmU.top20.addOthers.csv
human62_batch_effect2.taxa_counts.rel_abun.class.rmU.top20.csv
human62_batch_effect2.taxa_counts.rel_abun.class.rmU.top20.fillmin.scaled.csv
human62_batch_effect2.taxa_counts.rel_abun.class.rmU.top20.fillmin.scaled.index
...
kaiju/
temp/
top_taxa_asian_vs_euro/
top_taxa_asian_vs_female/
top_taxa_asian_vs_healthy/
...
top_taxa_male_vs_ill/
utest_asian_vs_euro/
utest_asian_vs_female/
...
utest_male_vs_female/
utest_male_vs_ill/
```
The boxplot folder contains the boxplots for all significant taxa crossing all taxonomy levels for all group-pair comparisons.
For example, boxplot_healthy_vs_ill/human62_batch_effect2.rel_abun.healthy_vs_ill.at_species.rel_abun.unequal.Alistipes.finegoldii.CAG.68.boxplot.pdf . We designed all file names to make them easy to understand. From this name, we know it is relative abundance of a significant species Alistipes.finegoldii.CAG.68 comparing between healthy and ill. Also, all plots generated by GEMINI are vector PDF format for the convenience of publication.

 ![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.rel_abun.healthy_vs_ill.at_species.rel_abun.unequal.Alistipes.finegoldii.CAG.68.boxplot.png)
 
 The figs fold contains the heatmap and barplots:
 ```
 (base) yh@superServer:taxa_analysis$ l figs/
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf
human62_batch_effect2.rel_abun.asian_vs_euro.at_family.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
....
human62_batch_effect2.taxa_counts.rel_abun.taxaID.rmU.top20.barplot.pdf
human62_batch_effect2.taxa_counts.rel_abun.taxaID.rmU.top20.fillmin.scaled.heatmap.pdf
```
 For example, figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_family.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf . GEMINI not only draws heatmap for unequal taxa, but also for equal taxa.
![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_family.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf)
 
 GEMINI produces barplot for all taxonomy levels. Here we have top 20 taxa, the 20 here is an adjustable parameter.
 ![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.phylum.rmU.top20.barplot.pdf)
 
 The csv files here are intermediate tables. GEMINI include these tables for user's convenience.
 
 The kaiju folder contains the taxonomy annotated contigs table:
 ```
 (base) yh@superServer:taxa_analysis$ l kaiju/
human62_batch_effect2_kaiju.ref  human62_batch_effect2_kaiju.ref.nm  human62_batch_effect2_kaiju.ref.nm.tsv
 ```
 
 The top_taxa folders contain the intermediate tables for top abundant taxa crossing all taxonomy levels for every group pair.
 ```
 (base) yh@superServer:taxa_analysis$ l top_taxa_asian_vs_euro/
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.rel_abun.equal.top20.csv
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.rel_abun.equal.top20.fillmin.scaled.csv
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.rel_abun.equal.top30.fillmin.scaled.index
...
human62_batch_effect2.rel_abun.asian_vs_euro.at_taxaID.rel_abun.unequal.top20.csv
human62_batch_effect2.rel_abun.asian_vs_euro.at_taxaID.rel_abun.unequal.top20.fillmin.scaled.csv
human62_batch_effect2.rel_abun.asian_vs_euro.at_taxaID.rel_abun.unequal.top20.fillmin.scaled.index
 ```
 
 The utest folders contain the statistical results.
 ```
 (base) yh@superServer:taxa_analysis$ l utest_asian_vs_euro/
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.ave_change.equal.csv     human62_batch_effect2.rel_abun.asian_vs_euro.at_phylum.ave_change.equal.csv
...
human62_batch_effect2.rel_abun.asian_vs_euro.at_order.u-test.two_sided.csv     human62_batch_effect2.rel_abun.asian_vs_euro.at_taxaID.u-test.two_sided.csv
 ```
 
## metabolism analysis
## alpha/beta diversity analysis
## lefse analysis


# common questions
### When running humann3, get `No MetaPhlAn BowTie2 database found (--index option)!` error.
Basically, the cause is the wrong installation of humann3.
1. google it. Try `metaphlan —install`. Then re-run snakemake.
2. if 1. not work. check the error log, try to find this sentence `Expecting location /the/expection/location`; `cd` to the location, examine any file truncation/loss. If there be, remove all files in the `/the/expection/location`, re-run `metaphlan —install`. Then re-run snakemake.

### How to deploy GEMINI on Slurm/PBS system?
Refer to this [answer](https://stackoverflow.com/questions/53545690/how-to-activate-a-specific-python-environment-as-part-of-my-submission-to-slurm). I tried and succeeded. 

### Can't locate File/Slurp.pm in @INC (@INC contains: /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at /home/yzz0191/anaconda3/envs/GEMINI/bin/abricate line 9.
I occured this when submitting GEMINI to Slurm system. The reason is `~/anaconda3/envs/GEMINI/bin/abricate`, the `abricate` was written in Perl, and the first line of `abricate` is `#!/usr/bin/env perl`. So, comment this line out, add a new head line `#! /path/to/your/anaconda3/bin/perl` (modify /path/to/your !), will solve it.
### When use GEMINI.yml create conda environment, occurred Bioconductor related issues.
`conda env remove GEMINI` to clean the failed GEMINI environment. Then use `GEMINI_without_bioconductor.yml` to create a new GEMINI environment. Then `conda activate GEMINI`, then type `R` to open the R command line, then install the R library in R scripts manually one by one. Then follow the left normal GEMINI install instruction. The R libraries include `mixOmics`,`gridExtra`,`sva`,`ggplot2`,`limma`,`grid`,`edgeR`,`DESeq2`,`pheatmap`,`wesanderson`,`ggpubr`,`Ipaper`,`reshape2`

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('mixOmics')
BiocManager::install("sva")
BiocManager::install("DESeq2")
install.packages('remotes')

### Humann error, `Please update your version of MetaPhlAn2 to v3.0`
run `conda install -c biobakery humann=3.1.1`. This is because humann internal error. Need to update to newer version of humann. After the installation, the original chocophlan database will be out of use (you will get `Please install the latest version of the database: v201901_v31` error). So run `humann_databases --download chocophlan full /path/to/your/databases/ --update-config yes` to update your chocophlan database for humann. Then, since you just updated the humann, you need to update the nucleotide and protein database record in humann_config too so that huamann can find the databases. Run `humann_config --update database_folders protein /path/to/your/uniref/database/` and `humann_config --update database_folders protein /path/to/your/utility_mapping/database/`.
### The assembly looks normal, but no plasmids/antibiotic/virus table generated
Reformat your assembly file.
from
`>human63_000000000001
TTTCCTTCGATGAGTTCTATGCCGTATATAATAAAAAGCATTCCGCTATTGAACAGCGTCTCGCAGAAAAAGGATTGCCGGAACATCTGCTTCATCGTAAGGAACGCAGACAGGAAAAACTGAATCATCCTGCTGTAAAAACGACAAAGCCCCACAGAAAGAAGAAAAAGAAACAGGTGTTCGAGCCGCTCTTGGAACAGAATGATGATTTCTTCTTTATTGCTGGTTATACTTCTGGCGGTGCCCCTTATGGTGTCACATGGGAAGAAATGGGACTAGAGCCTTGGGAAGAACTTGTATAAATATTATTGCCATTGCCGATTGCCAAAAGCA`
to 
`>human63_000000000001
TTTCCTTCGATGAGTTCTATGCCGTATATAATAAAAAGCATTCCGCTATTGAACAGCGTC
TCGCAGAAAAAGGATTGCCGGAACATCTGCTTCATCGTAAGGAACGCAGACAGGAAAAAC
TGAATCATCCTGCTGTAAAAACGACAAAGCCCCACAGAAAGAAGAAAAAGAAACAGGTGT
TCGAGCCGCTCTTGGAACAGAATGATGATTTCTTCTTTATTGCTGGTTATACTTCTGGCG
GTGCCCCTTATGGTGTCACATGGGAAGAAATGGGACTAGAGCCTTGGGAAGAACTTGTAT
AAATATTATTGCCATTGCCGATTGCCAAAAGCA`

