# GEMINI
GEMINI: a comprehensive, reliable, and user-friendly downstream analysis pipeline for whole-metagenome shotgun sequence.

# GEMINI installation
0. GEMINI only supports Illumina reads, not for Ion Torrent yet.
1. Download `GEMINI.yml`
2. Create the conda environment `conda env create --name GEMINI --file GEMINI.yml `
3. Activate the environment `conda activate GEMINI`.
4. Check you're using right humann `which humann`. Read step 19.
5. Download 1st humann databases `humann_databases --download chocophlan full /path/to/databases --update-config yes`
6. Download 2nd humann databases `humann_databases --download uniref uniref90_diamond /path/to/databases --update-config yes`
7. Download 3rd humann databases `humann_databases --download utility_mapping full /path/to/databases --update-config yes`
8. Check your humann databases `cd /path/to/databases` then run `ll chocophlan/ | wc -l` you get a number >= 11289. run `ll uniref` you should see a `uniref90_201901b_full.dmnd` (or newer) file with >= 34G size. Run `ll utility_mapping` you should see >= 21 files and all of them have > 3M size (or some of them truncated during download).
9. Check your humann's misc folder. Located at `/path/to/your/anaconda/envs/GEMINI/lib/python3.9/site-packages/humann/data/misc/`. Due to conda's unknown error, usually the files are missing. To know full file list, check README.txt in the misc folder. To download the files, check https://github.com/biobakery/humann/tree/master/humann/data/misc . Or unzip the `misc.zip` in the `utils` folder.
10. Check your humann by running `humann -i sample_reads.fastq -o sample_results` (prepare sample_reads.fastq by yourself)
11. Check you're using right kaiju `which kaiju`
12. Download kaiju databases `mkdir /path/to/kaijudb` then `cd /path/to/kaijudb` then `kaiju-makedb -s nr_euk` (this takes a long time and space and memory). Or you can download and unzip the annotation files from [kaiju server](https://kaiju.binf.ku.dk/server).
13. Check abricate databases `abricate --list`, you should see 9 databases (argannot,card,ecoh,ecoli_vf,megares,ncbi,plasmidfinder,resfinder,vfdb). If you didn't, go to download abricate [databases](https://github.com/tseemann/abricate/tree/master/db), or use the `db.zip` file in `utils` folder, unzip them under your `db` folder. The location of `db` folder can be seen by running `abricate --help`, see the `--datadir` line. Then `cd` to the `db` folder, run `abricate --setupdb`.
14. Install `iPaper` R package. type `R` in command line, then in the R concle, type in `install.packages("remotes")`, then type in `remotes::install_github("kongdd/Ipaper")`. When R concle asks you whether to update other packages, choose `none`. After installation, type in `library(Ipaper)`, if no error occurs, then you're good to contine.
15. Open `Snakefileconfig.yml` under `GEMINI` folder, modify the bwa,kaiju,python3,Rscript,...lefse_run parameters to the executable command lines in your environment. To make sure all command line works, you could test the command line one by one in your linux shell.
16. Test GEMINI. `cd parent/folder/of/GEMINI`. Prepare some example data. Open and modify the `GEMINI/GEMINI_contig.tsv` to make sure the data_dir is right.Then run `nohup snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete &`. This command will run GEMINI in backend with a log file `nohup.out`.
17. If you occured any errors. check the log to debug. Or check the log file in `name_of_your_assembly/log` folder. You can use time stamps to refer which rule is error, or to understand the error information. After debugging, delete the `name_of_your_assembly/log/name_of_the_error_rule.done`. and rerun the `snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete`. GEMINI will auto-resume.
18. You should see all the outputs in `name_of_your_assembly` folder. If no errors occured. then you're good to go.
19. In a summary, GEMINI itself is just a list of scripts, which is easy to use. The above installation guide actually is helping you to install other tools, such as humann3 and kaiju. On the other hand, if you have installed these tools somewhere else, you can just modify the `Snakefile_config.yml` file to skip the above installation procedure.

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

The ***fq_dir*** column contains the absolute location of each fastq file. If you use PE end sequencing, `cat` them to one file, because humann3 asked so. For example `cat sample1_seq_r1.fq sample1_seq_r2.fq > sample1_seq.fq`. You should fill the absolute location of sample1_seq.fq in this column, do not fill multiple files.

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
The usuage of GEMINI is very easy and light. The first step is to install GEMINI, the second step is to prepare the `GEMINI/GEMINI_config.tsv` (for experiment information), the third step is to modify the `GEMINI/Snakemake_config.yml` (for software parameters). The last step is to run `nohup snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete &`.

Here are some suggestions for beginners:
1. The working directory must be the parent directory of `GEMINI`, becasue some scripts in Snakemake.py use relative path. To be more specific, you have a folder `test`, you have all scripts of GEMINI in `test` folder, you have `test\Snakefile.py` and `test\GEMINI`. You `cd test`, then run  `nohup snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete &`, you will get all outputs in `test\name_of_your_assembly` folder.
2. Snakemake (the framework GEMINI relied on) will lock working directory during running. So you should prepare two working directories if you're running two GEMINI pipelines (or any other Snakemake based softwares). To be more specific, `mkdir folder1/` and `mkdir folder2/`, copy the entire GEMINI folder to `folder1/` and `folder2/`. Then `cd folder1`, run GEMINI. Then `cd folder2`, run GEMINI.
3. One GEMINI process only take one assembly. If you have multi assemblies to analyze, run GEMINI multiple times (in parallel) (in different working directories).
4. Humann analysis is very time/computation consuption. I personally prefer to use computer cluster to distributedly run Humann. The Snakemake based GEMINI can also be deployed on cluster, but maybe not [easy](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). So GEMINI provide `skip_humann_init` option. Set `skip_humann_init = True` in `GEMINI_config.tsv`, then GEMINI will not run human_init rule, but to check the human results under folder `name_of_the_assembly/metabolism_analysis/humann3/ori_results/`, so you need to put the humann output with suffix as genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the folder. If `skip_humann_init = True`, GEMINI will check the outputs existence first then skip the humann step. Make sure these humann outptus are from the same fastq you offered to GEMINI.
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


 ![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.rel_abun.healthy_vs_ill.at_species.rel_abun.unequal.Alistipes.finegoldii.CAG.68.boxplot.pdf) (If image failed to show, click it to view.)
 
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
 (If image failed to show, click it to view.)
 GEMINI produces barplot for all taxonomy levels. Here we have top 20 taxa, the 20 here is an adjustable parameter.
 
 ![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.phylum.rmU.top20.barplot.pdf)
 (If image failed to show, click it to view.)
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
 
## batch effect
If rm_batch_effect is True, GEMINI will visualize the principal components (PCA) as well as the variance distribution. To check the batch effect, users can compare the plots before and after batch effect removal. 
For example,

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.family.rmU.batch_effect_PCA.pdf)
(If image failed to show, click it to view.)
```
(base) yh@superServer:human62_batch_effect2$ l batch_effect/
human62_batch_effect2.taxa_counts.rel_abun.class.rmU.batch_effect_PCA.pdf   human62_batch_effect2.taxa_counts.rel_abun.phylum.rmU.batch_effect_PCA.pdf
human62_batch_effect2.taxa_counts.rel_abun.family.rmU.batch_effect_PCA.pdf  human62_batch_effect2.taxa_counts.rel_abun.species.rmU.batch_effect_PCA.pdf
human62_batch_effect2.taxa_counts.rel_abun.genus.rmU.batch_effect_PCA.pdf   human62_batch_effect2.taxa_counts.rel_abun.superkingdom.rmU.batch_effect_PCA.pdf
human62_batch_effect2.taxa_counts.rel_abun.order.rmU.batch_effect_PCA.pdf   human62_batch_effect2.taxa_counts.rel_abun.taxaID.rmU.batch_effect_PCA.pdf
```

## benchmark
This folder contains the log information of each computation module in GEMINI. Users can check the benchmark files for diagnosis.
For example,
```
(base) yh@superServer:human62_batch_effect2$ l benchmark/
alpha_beta_diversity.benchmark   heatmap_taxa.benchmark     humann_utest.benchmark                 rel_abun2lefse_taxa.benchmark          taxa_boxplot.benchmark
antibiotic_analysis.benchmark    humann2rel_abun.benchmark  lefse_taxa.benchmark                   rel_abun_utest.benchmark               virulence_analysis.benchmark
counts_table2rel_abun.benchmark  humann_annotate.benchmark  merge_counts_pure_and_kaiju.benchmark  scale_humann_rel_abun_table.benchmark  visualize_batch_effect.benchmark
extract_top_humann.benchmark     humann_group.benchmark     paste_counts_table.benchmark           scale_taxa_rel_abun_table.benchmark
extract_top_taxa.benchmark       humann_init.benchmark      plasmid_analysis.benchmark             skip_kaiju.benchmark
heatmap_humann.benchmark         humann_output.benchmark    rel_abun2lefse_humann.benchmark        taxa_barplots.benchmark
(base) yh@superServer:human62_batch_effect2$ cat benchmark/heatmap_taxa.benchmark
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time
1353.9893	0:22:33	697.34	4668.55	656.18	666.59	25.41	34.90	1.33	10.34
```

## diversity_analysis
This folder contains all the alpha and beta diversity analysis crossing every group-pair.
```
(base) yh@superServer:human62_batch_effect2$ l diversity_analysis/
alpha_beta_asian_vs_euro/     alpha_beta_asian_vs_ill/   alpha_beta_female_vs_euro/   alpha_beta_healthy_vs_female/  alpha_beta_male_vs_euro/
alpha_beta_asian_vs_female/   alpha_beta_asian_vs_male/  alpha_beta_female_vs_ill/    alpha_beta_healthy_vs_ill/     alpha_beta_male_vs_female/
alpha_beta_asian_vs_healthy/  alpha_beta_euro_vs_ill/    alpha_beta_healthy_vs_euro/  alpha_beta_healthy_vs_male/    alpha_beta_male_vs_ill/
```
Each sub-folder contains alpha and beta diversity analysis intermediate tables and graphs.
```
(base) yh@superServer:diversity_analysis$ l alpha_beta_healthy_vs_ill
alpha_diversity.at_genus.csv                                               index.txt                                    PCoA23.bray.at_genus.pdf
alpha_diversity.at_species.csv                                             Inv_Simpson.alpha_diveristy.at_.genus.pdf    PCoA23.bray.at_species.pdf
beta_pcoa_P-value.bray.at_genus.csv                                        Inv_Simpson.alpha_diveristy.at_.species.pdf  PCoA23.jaccard.at_genus.pdf
beta_pcoa_P-value.bray.at_species.csv                                      jaccard.at_genus.csv                         PCoA23.jaccard.at_species.pdf
beta_pcoa_P-value.jaccard.at_genus.csv                                     jaccard.at_species.csv                       Pielou_evenness.alpha_diveristy.at_.genus.pdf
beta_pcoa_P-value.jaccard.at_species.csv                                   PCoA12.bray.at_genus.pdf                     Pielou_evenness.alpha_diveristy.at_.species.pdf
bray.at_genus.csv                                                          PCoA12.bray.at_species.pdf                   Shannon.alpha_diveristy.at_.genus.pdf
bray.at_species.csv                                                        PCoA12.jaccard.at_genus.pdf                  Shannon.alpha_diveristy.at_.species.pdf
bray_ln.at_genus.csv                                                       PCoA12.jaccard.at_species.pdf                Simpson.alpha_diveristy.at_.genus.pdf
bray_ln.at_species.csv                                                     PCoA13.bray.at_genus.pdf                     Simpson.alpha_diveristy.at_.species.pdf
group.csv                                                                  PCoA13.bray.at_species.pdf                   Simpson_evenness.alpha_diveristy.at_.genus.pdf
human62_batch_effect2.taxa_counts.rel_abun.genus.rmU.healthy_vs_ill.csv    PCoA13.jaccard.at_genus.pdf                  Simpson_evenness.alpha_diveristy.at_.species.pdf
human62_batch_effect2.taxa_counts.rel_abun.species.rmU.healthy_vs_ill.csv  PCoA13.jaccard.at_species.pdf
```
For example, 

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/PCoA12.bray.at_species.pdf)

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/Shannon.alpha_diveristy.at_.species.pdf)
(If image failed to show, click it to view.)
## LDA analysis
This folder contains the metabolism and taxonomy LDA results for every group-pair.
```
(base) yh@superServer:human62_batch_effect2$ l LDA_analysis/
figs_humann/              humann_asian_vs_ill/    humann_healthy_vs_euro/    humann_male_vs_female/  taxa_asian_vs_ill/    taxa_healthy_vs_euro/    taxa_male_vs_female/
figs_taxa/                humann_asian_vs_male/   humann_healthy_vs_female/  humann_male_vs_ill/     taxa_asian_vs_male/   taxa_healthy_vs_female/  taxa_male_vs_ill/
humann_asian_vs_euro/     humann_euro_vs_ill/     humann_healthy_vs_ill/     taxa_asian_vs_euro/     taxa_euro_vs_ill/     taxa_healthy_vs_ill/
humann_asian_vs_female/   humann_female_vs_euro/  humann_healthy_vs_male/    taxa_asian_vs_female/   taxa_female_vs_euro/  taxa_healthy_vs_male/
humann_asian_vs_healthy/  humann_female_vs_ill/   humann_male_vs_euro/       taxa_asian_vs_healthy/  taxa_female_vs_ill/   taxa_male_vs_euro/
```
The figs_humann and figs_taxa contain the lefSE results.
```
(base) yh@superServer:LDA_analysis$ l figs_humann/
allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.asian_vs_euro.rel_abun.unequal.lefse.pdf
...
allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.male_vs_female.rel_abun.unequal.lefse.pdf
allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.male_vs_ill.rel_abun.unequal.lefse.pdf
(base) yh@superServer:LDA_analysis$ l figs_taxa/
human62_batch_effect2.rel_abun.asian_vs_euro.at_genus.rel_abun.unequal.lefse.pdf
...
human62_batch_effect2.rel_abun.female_vs_ill.at_taxaID.rel_abun.unequal.lefse.pdf 
```
For example, 

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.asian_vs_euro.rel_abun.unequal.lefse.pdf)

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_genus.rel_abun.unequal.lefse.pdf)
(If image failed to show, click it to view.)
Other sub-folders contain the intermediate tables for user's convenience.
```
(base) yh@superServer:LDA_analysis$ l taxa_asian_vs_euro
human62_batch_effect2.rel_abun.asian_vs_euro.at_class.rel_abun.equal.lefse.format     human62_batch_effect2.rel_abun.asian_vs_euro.at_phylum.rel_abun.equal.lefse.format
...
human62_batch_effect2.rel_abun.asian_vs_euro.at_order.rel_abun.unequal.lefse.tsv      human62_batch_effect2.rel_abun.asian_vs_euro.at_taxaID.rel_abun.unequal.lefse.tsv
```

## metabolism analysis
The folder contains the metabolism analysis results, mainly based on the 5 independent databases (which are selectable in `Snakemake_config.yml` file)
```
(base) yh@superServer:human62_batch_effect2$ l metabolism_analysis
figs/  humann3/
```
The `humann3` folder contains all intermediate tables for our humann3 process. If users want to skip the humann_init (as former suggested), they should `mkdir -p path_to_metabolism_analysis/metabolism_analysis/humann3/ori_results` and prepare all the humann3 generated tables with suffix as genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the `ori_results` folder, then set skip_humann_init to be True.
```
(base) yh@superServer:metabolism_analysis$ l humann3/
allSamples/  male/                         top_humann_asian_vs_ill/     top_humann_healthy_vs_female/  utest_asian_vs_euro/     utest_female_vs_euro/     utest_male_vs_euro/
asian/       ori_results/                  top_humann_asian_vs_male/    top_humann_healthy_vs_ill/     utest_asian_vs_female/   utest_female_vs_ill/      utest_male_vs_female/
euro/        output/                       top_humann_euro_vs_ill/      top_humann_healthy_vs_male/    utest_asian_vs_healthy/  utest_healthy_vs_euro/    utest_male_vs_ill/
female/      top_humann_asian_vs_euro/     top_humann_female_vs_euro/   top_humann_male_vs_euro/       utest_asian_vs_ill/      utest_healthy_vs_female/
healthy/     top_humann_asian_vs_female/   top_humann_female_vs_ill/    top_humann_male_vs_female/     utest_asian_vs_male/     utest_healthy_vs_ill/
ill/         top_humann_asian_vs_healthy/  top_humann_healthy_vs_euro/  top_humann_male_vs_ill/        utest_euro_vs_ill/       utest_healthy_vs_male/
```
The figs folder contains all the heatmap plots for every independent database crossing every group-pair.
```
(base) yh@superServer:metabolism_analysis$ l figs/
allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.asian_vs_euro.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.asian_vs_euro.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf
...
allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.male_vs_ill.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.male_vs_ill.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf
allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.top20.fillmin.scaled.heatmap.pdf
```
For example,

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.healthy_vs_ill.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf)
(If image failed to show, click it to view.)
## antibiotic_analysis/plasmid_analysis/virulence_analysis
These folders contain the dist plot and heatmap plot for all features with responding taxonomy, as well as intermediate tables.
```
(base) yh@superServer:human62_batch_effect2$ l antibiotic_analysis
genes_asian_vs_euro_distrplot.pdf           genes_female_vs_euro_distrplot.pdf         genes_order_male_vs_euro_heatmap.pdf         genes_taxa_counts_asian_vs_ill.tsv
...
genes_family_male_vs_female_heatmap.pdf     genes_order_healthy_vs_ill_heatmap.pdf     genes_taxa_counts_asian_vs_female.tsv
genes_family_male_vs_ill_heatmap.pdf        genes_order_healthy_vs_male_heatmap.pdf    genes_taxa_counts_asian_vs_healthy.tsv
```
For example,

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/genes_asian_vs_euro_distrplot.pdf)

![image](https://github.com/Y-H-Joe/GEMINI/blob/main/figs/genes_class_healthy_vs_ill_heatmap.pdf)
(If image failed to show, click it to view.)
GEMINI also provide all the annotation information in `utils/abricate.annoatations.txt`, which can assist users to determine the antibiotic genes/plasmids/virulence factors. Also, the intermediate tables are helpful. The `*.antibiotic.tsv` `*.virulence.tsv` `*.plasmidfinder.tsv` are summary tables.
```
(base) yh@superServer:antibiotic_analysis$ head human62_batch_effect2.antibiotic.tsv
#FILE	SEQUENCE	START	END	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION
/home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa	human63_000000001502	165	274	erm(G)_1	1-110/735	===............	0/0	14.97	100.00	resfinder	M15332
/home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa	human63_000000002481	1	126	tet(Q)_1	916-1041/1926	.......==......	0/0	6.54	99.21	resfinder	L33696
/home/yh/GEMINI/assembly/human63.contigs.nonrd.rn.fa	human63_000000002481	115	350	tet(Q)_1	868-1103/1926	......===......	0/0	12.25	96.61	resfinder	L336
```

## log
This folder contains the log for all processing modules of GEMINI. Notably, GEMINI use these log files with suffix of `.done` to judge the status of each rule. Users can check the `Snakemake.py`, if the `done` files for certain rule exists, GEMINI will regart it as successfully finished. If users want to re-run certain rules, they should remove the corressponding `done` files in the first place.
```
(base) yh@superServer:human62_batch_effect2$ cd log/
(base) yh@superServer:log$ l
alpha_beta_diversity.done                         extract_top_taxa.done                           heatmap_taxa.log_taxa               plasmid_analysis.done
alpha_beta_diversity.log                          extract_top_taxa.log                            humann2rel_abun.done                plasmid_analysis.log
...
extract_top_humann.log_male_vs_ill_equal          heatmap_taxa.done                               paste_counts_table.done             visualize_batch_effect.done
extract_top_humann.log_male_vs_ill_unequal        heatmap_taxa.log                                paste_counts_table.log              visualize_batch_effect.log
```

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


```
>human63_000000000001
TTTCCTTCGATGAGTTCTATGCCGTATATAATAAAAAGCA
TTCCGCTATTGAACAGCGTCTCGCAGAAAAAGGATTGCCG
GAACATCTGCTTCATCGTAAGGAACGCAGACAGGAAAAAC
TGAATCATCCTGCTGTAAAAACGACAAAGCCCCACAGAAA
GAAGAAAAAGAAACAGGTGTTCGAGCCGCTCTTGGAACAG
AATGATGATTTCTTCTTTATTGCTGGTTATACTTCTGGCG
GTGCCCCTTATGGTGTCACATGGGAAGAAATGGGACTAGA
GCCTTGGGAAGAACTTGTATAAATATTATTGCCATTGCCG
ATTGCCAAAAGCA
```


