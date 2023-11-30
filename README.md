# OUTPOST: a comprehensive and useful downstream analysis pipeline for whole-metagenome shotgun sequencing.
---
 

# OUTPOST installation
---
```
# 1. Download `OUTPOST.yml` in `install` folder.
git clone https://github.com/Y-H-Joe/OUTPOST.git
cd OUTPOST

# 2. Create the conda environment `conda env create --name OUTPOST --file OUTPOST.yml `, the bioconductor has unknown errors.
conda env create --name OUTPOST --file install/OUTPOST_without_bioconductor.yml

# 3. Activate the environment `conda activate OUTPOST`.
conda activate OUTPOST

# 4. install. Databases will also be donwnloaded, require ~500GB. The downloading time is based on your bandwith.
# make sure your working directory is OUTPOST folder
bash install/install_guideline.sh
```
---
 

  
# OUTPOST usuage
---
```
# 1. prepare the `OUTPOST/OUTPOST_config.tsv` (for experiment, data, group information)
# 2. modify the `OUTPOST/Snakemake_config.yml` (for OUTPOST parameters)
# 3. run OUTPOST
nohup snakemake --cores 32 --verbose -s ./OUTPOST_run.py --rerun-incomplete &
```
We prepared the OUTPOST results of cat microbiome dataset (described in our article) in [figshare](https://figshare.com/articles/figure/OUTPOST_results_of_cat_microbiome_dataset/24082542), which can be an example for users to testify their OUTPOST installation or usage. 

To Be Remind:
```
1. The working directory must be the parent directory of `OUTPOST`, becasue some scripts in Snakemake.py use relative path. To be more specific, you have a folder `test`, you have all scripts of OUTPOST in `test` folder, you have `test\Snakefile.py` and `test\OUTPOST`. You `cd test`, then run  `nohup snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete &`, you will get all outputs in `test\name_of_your_assembly` folder.
2. Snakemake (the framework OUTPOST relied on) will lock working directory during running. So you should prepare two working directories if you're running two OUTPOST pipelines (or any other Snakemake based softwares). To be more specific, `mkdir folder1/` and `mkdir folder2/`, copy the entire OUTPOST folder to `folder1/` and `folder2/`. Then `cd folder1`, run OUTPOST. Then `cd folder2`, run OUTPOST.
3. One OUTPOST process only take one assembly. If you have multi assemblies to analyze, run OUTPOST multiple times (in parallel) (in different working directories).
4. Humann analysis is very time/computation consuption. I personally prefer to use computer cluster to distributedly run Humann. The Snakemake based OUTPOST can also be deployed on cluster, but maybe not [easy](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). So OUTPOST provide `skip_humann_init` option. Set `skip_humann_init = True` in `Snakefile_config.yml`, then OUTPOST will not run human_init rule, but to check the human results under folder `name_of_the_assembly/metabolism_analysis/humann3/ori_results/`, so you need to put the humann output with suffix as genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the folder. If `skip_humann_init = True`, OUTPOST will check the outputs existence first then skip the humann step. Make sure these humann outptus are from the same fastq you offered to OUTPOST.
5. OUTPOST offers `skip_assembly_analysis`. You can skip the assembly analysis if you have a large assembly and a number of groups which will save a lot of time about MAG tables generation.
6. Do not change the column names of OUTPOST_config.tsv.
7. Refer to OUTPOST/OUTPOST_config.explanation.txt for more details.
```


---
 


  
# OUTPOST outputs
---
We prepared the OUTPOST results of cat microbiome dataset (described in our article) in [figshare](https://figshare.com/articles/figure/OUTPOST_results_of_cat_microbiome_dataset/24082542), which can be an example for users to testify their OUTPOST installation or usage. 


Each OUTPOST run accepts one assembly, all outputs are categorized in the folder named by the assembly.
```
total 60K
drwxrwxr-x 15 yihang yihang 4.0K Sep  4 12:03 ./
drwxrwxr-x  5 yihang yihang 4.0K Sep  4 12:58 ../
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 09:37 antibiotic_genes_analysis/
drwxrwxr-x  3 yihang yihang 4.0K Sep  4 10:08 assembly_analysis/
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 10:08 batch_effect/
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 09:37 benchmark/
drwxrwxr-x  3 yihang yihang 4.0K Sep  4 12:06 biomarkers_analysis/
drwxrwxr-x  3 yihang yihang 4.0K Sep  4 10:06 diversity_analysis/
drwxrwxr-x  4 yihang yihang 4.0K Sep  4 06:55 function_analysis/
drwxrwxr-x  6 yihang yihang 4.0K Sep  4 10:10 LDA_analysis/
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 12:32 log/
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 08:22 plasmids_analysis/
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 12:31 report/
drwxrwxr-x  9 yihang yihang 4.0K Sep  4 12:06 taxonomy_analysis/
drwxrwxr-x  2 yihang yihang 4.0K Sep  4 08:52 virulence_factors_analysis/
```


### OUTPOST outputs: assembly analysis
---
```
assembly_analysis/
├── [175M]  cat.normal_vs_obese.contig_table.processed.tsv
├── [ 61M]  cat.normal_vs_obese.contig_table.tsv
├── [ 954]  cat.taxa_counts.summarized.tsv
├── [ 53M]  cat.taxa_counts.tsv
└── [4.0K]  metagenemark
    ├── [1.3G]  cat.gtf
    ├── [939M]  cat.nucl.fa
    ├── [920M]  cat.nucl.nonrd.fa
    ├── [ 56M]  cat.nucl.nonrd.fa.clstr
    └── [417M]  cat.prot.fa
```
I set `skip_assembly_analysis` to True, so I only have one table. The MAG analysis table can refer to Table S1 of our paper.


### OUTPOST outputs: taxonomy analysis
---
```
taxonomy_analysis/
├── [628K]  boxplot_normal_vs_obese
│   ├── [4.0K]  class
│   ├── [ 28K]  family
│   ├── [ 92K]  genus
│   ├── [ 12K]  order
│   ├── [4.0K]  phylum
│   ├── [ 92K]  species
│   ├── [4.0K]  superkingdom
│   └── [444K]  taxaID
├── [4.0K]  counts_tables
│   ├── [ 31K]  cat.taxa_counts.rel_abun.class.csv
│   ├── [ 31K]  cat.taxa_counts.rel_abun.class.rmU.csv
│   ├── [7.0K]  cat.taxa_counts.rel_abun.class.rmU.top20.addOthers.csv
...
│   ├── [6.2K]  cat.taxa_counts.rel_abun.taxaID.rmU.top20.csv
│   ├── [6.4K]  cat.taxa_counts.rel_abun.taxaID.rmU.top20.fillmin.scaled.csv
│   └── [ 193]  cat.taxa_counts.rel_abun.taxaID.rmU.top20.fillmin.scaled.index
├── [4.0K]  figs
│   ├── [7.9K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
│   ├── [7.9K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf
...
│   ├── [8.1K]  cat.taxa_counts.rel_abun.taxaID.rmU.top20.barplot.pdf
│   ├── [7.7K]  cat.taxa_counts.rel_abun.taxaID.rmU.top20.fillmin.scaled.heatmap.pdf
│   └── [8.1K]  Rplots.pdf
├── [4.0K]  kaiju
│   ├── [ 66M]  cat_kaiju.ref
│   ├── [ 96M]  cat_kaiju.ref.nm
│   └── [ 36M]  cat_kaiju.ref.nm.tsv
├── [4.0K]  temp
│   ├── [ 19M]  cat.counts.pure
│   ├── [9.1M]  D001_obese_contigs_sorted.cat.bam.idxstats
│   ├── [808K]  D001_obese_contigs_sorted.cat.bam.idxstats.pure
...
│   ├── [9.1M]  Z116_normal_contigs_sorted.cat.bam.idxstats
│   └── [783K]  Z116_normal_contigs_sorted.cat.bam.idxstats.pure
├── [ 12K]  top_taxa_normal_vs_obese
│   ├── [6.8K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.equal.top20.csv
│   ├── [6.5K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.equal.top20.fillmin.scaled.csv
...
│   └── [ 297]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.unequal.top20.fillmin.scaled.index
└── [4.0K]  utest_normal_vs_obese
    ├── [2.3K]  cat.rel_abun.normal_vs_obese.at_class.ave_change.equal.csv
...
    ├── [833K]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.unequal.csv
    └── [630K]  cat.rel_abun.normal_vs_obese.at_taxaID.u-test.two_sided.csv
```
The `boxplot` folder contains the boxplots for all significant taxa crossing all taxonomy levels for all group-pair comparisons.
For example, boxplot_healthy_vs_ill/human62_batch_effect2.rel_abun.healthy_vs_ill.at_species.rel_abun.unequal.Alistipes.finegoldii.CAG.68.boxplot.pdf . We designed all file names to make them easy to understand. From this name, we know it is relative abundance of a significant species Alistipes.finegoldii.CAG.68 comparing between healthy and ill. Also, all plots generated by OUTPOST are vector PDF format for the convenience of publication.


 ![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.rel_abun.healthy_vs_ill.at_species.rel_abun.unequal.Alistipes.finegoldii.CAG.68.boxplot.pdf) (If image failed to show, click it to view.)
 
 The `figs` folder contains the heatmap and barplots:
 OUTPOST not only draws heatmap for unequal taxa, but also for equal taxa.
 
![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_family.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf)
 (If image failed to show, click it to view.)
 OUTPOST produces barplot for all taxonomy levels. Here we have top 20 taxa, the 20 here is an adjustable parameter.
 
 ![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.phylum.rmU.top20.barplot.pdf)
 (If image failed to show, click it to view.)

 The `csv` files here are intermediate tables. OUTPOST include these tables for user's convenience.
 The `kaiju` folder contains the taxonomy annotated contigs table:
 The `top_taxa` folders contain the intermediate tables for top abundant taxa crossing all taxonomy levels for every group pair.
 The `utes`t folders contain the statistical results.


### OUTPOST outputs: batch effect
---
If `rm_batch_effect` is `True`, OUTPOST will visualize the principal components (PCA) as well as the variance distribution. To check the batch effect, users can compare the plots before and after batch effect removal. 
For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.family.rmU.batch_effect_PCA.pdf)
(If image failed to show, click it to view.)
```
batch_effect/
├── [ 17K]  cat.taxa_counts.rel_abun.class.rmU.not_rm_batch_effect_PCA.pdf
├── [ 17K]  cat.taxa_counts.rel_abun.family.rmU.not_rm_batch_effect_PCA.pdf
...
├── [ 17K]  cat.taxa_counts.rel_abun.species.rmU.not_rm_batch_effect_PCA.pdf
├── [ 16K]  cat.taxa_counts.rel_abun.superkingdom.rmU.not_rm_batch_effect_PCA.pdf
└── [ 17K]  cat.taxa_counts.rel_abun.taxaID.rmU.not_rm_batch_effect_PCA.pdf
```


### OUTPOST outputs: diversity analysis
---
This folder contains all the alpha and beta diversity analysis crossing every group-pair.
```
(base) yh@superServer:human62_batch_effect2$ l diversity_analysis/
alpha_beta_asian_vs_euro/     alpha_beta_asian_vs_ill/   alpha_beta_female_vs_euro/   alpha_beta_healthy_vs_female/  alpha_beta_male_vs_euro/
alpha_beta_asian_vs_female/   alpha_beta_asian_vs_male/  alpha_beta_female_vs_ill/    alpha_beta_healthy_vs_ill/     alpha_beta_male_vs_female/
alpha_beta_asian_vs_healthy/  alpha_beta_euro_vs_ill/    alpha_beta_healthy_vs_euro/  alpha_beta_healthy_vs_male/    alpha_beta_male_vs_ill/
```
Each sub-folder contains alpha and beta diversity analysis intermediate tables and graphs.
```
diversity_analysis/
└── [4.0K]  alpha_beta_normal_vs_obese
    ├── [1.8K]  alpha_diversity.at_genus.csv
    ├── [1.8K]  alpha_diversity.at_species.csv
    ├── [  32]  beta_pcoa_P-value.bray.at_genus.csv
    ├── [  31]  beta_pcoa_P-value.bray.at_species.csv
...
    ├── [4.9K]  bray_ln.at_species.csv
    ├── [497K]  cat.taxa_counts.rel_abun.genus.rmU.normal_vs_obese.csv
    ├── [2.0M]  cat.taxa_counts.rel_abun.species.rmU.normal_vs_obese.csv
    ├── [ 453]  group.csv
    ├── [ 297]  index.txt
    ├── [5.7K]  Inv_Simpson.alpha_diveristy.at_genus.pdf
    ├── [5.7K]  Inv_Simpson.alpha_diveristy.at_species.pdf
    ├── [4.9K]  jaccard.at_genus.csv
...
    ├── [6.9K]  PCoA12.jaccard.at_genus.pdf
    ├── [6.9K]  PCoA12.jaccard.at_species.pdf
...
    ├── [5.7K]  Shannon.alpha_diveristy.at_genus.pdf
    ├── [5.7K]  Shannon.alpha_diveristy.at_species.pdf
    ├── [5.7K]  Simpson.alpha_diveristy.at_genus.pdf
    ├── [5.7K]  Simpson.alpha_diveristy.at_species.pdf
    ├── [5.8K]  Simpson_evenness.alpha_diveristy.at_genus.pdf
    └── [5.8K]  Simpson_evenness.alpha_diveristy.at_species.pdf
```
For example, 

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/PCoA12.bray.at_species.pdf)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/Shannon.alpha_diveristy.at_.species.pdf)
(If image failed to show, click it to view.)


### OUTPOST outputs: LDA analysis
---
This folder contains the metabolism and taxonomy LDA results for every group-pair.
```
LDA_analysis/
├── [4.0K]  figs_humann
│   ├── [ 70K]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.lefse.pdf
│   ├── [ 27K]  allSamples_genefamilies_uniref90names_relab_ko_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.lefse.pdf
...
│   └── [ 39K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.lefse.pdf
├── [4.0K]  figs_taxa
│   ├── [ 11K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.unequal.lefse.pdf
│   ├── [ 14K]  cat.rel_abun.normal_vs_obese.at_family.rel_abun.unequal.lefse.pdf
...
│   ├── [ 13K]  cat.rel_abun.normal_vs_obese.at_species.rel_abun.equal.lefse.pdf
│   ├── [ 36K]  cat.rel_abun.normal_vs_obese.at_species.rel_abun.unequal.lefse.pdf
│   └── [ 26K]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.unequal.lefse.pdf
├── [ 12K]  humann_normal_vs_obese
│   ├── [3.9M]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.lefse.format
│   ├── [959K]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.lefse.res
│   ├── [3.8M]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.lefse.tsv
...
│   ├── [181K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.lefse.format
│   ├── [ 93K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.lefse.res
│   └── [253K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.lefse.tsv
└── [ 12K]  taxa_normal_vs_obese
    ├── [ 13K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.equal.lefse.format
    ├── [2.3K]  cat.rel_abun.normal_vs_obese.at_class.rel_abun.equal.lefse.res
...
    ├── [312K]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.lefse.res
    ├── [548K]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.unequal.lefse.format
    ├── [156K]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.unequal.lefse.res
    └── [833K]  cat.rel_abun.normal_vs_obese.at_taxaID.rel_abun.unequal.lefse.tsv
```

The `figs_humann` and `figs_taxa` folders contain the lefSE results.

For example, 

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.asian_vs_euro.rel_abun.unequal.lefse.pdf)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_genus.rel_abun.unequal.lefse.pdf)
(If image failed to show, click it to view.)


Other sub-folders contain the intermediate tables for user's convenience.


### OUTPOST outputs: function analysis
---
The folder contains the metabolism analysis results, mainly based on the 5 independent databases (which are selectable in `Snakemake_config.yml` file)
```
function_analysis/
├── [4.0K]  figs
│   ├── [7.9K]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
...
│   ├── [8.0K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.top20.fillmin.scaled.heatmap.pdf
│   ├── [8.3K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf
│   ├── [8.1K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.top20.fillmin.scaled.heatmap.pdf
│   └── [8.2K]  Rplots.pdf
└── [4.0K]  humann3
    ├── [ 32K]  allSamples
    │   ├── [3.1M]  111030_normal_genefamilies_cpm_eggnog.tsv
    │   ├── [451K]  111030_normal_genefamilies_cpm_ko.tsv
...
    │   ├── [163K]  Z116_normal_pathcoverage_cpm.tsv
    │   ├── [166K]  Z116_normal_pathcoverage_relab.tsv
    │   └── [175K]  Z116_normal_pathcoverage.tsv
    ├── [ 16K]  normal
    │   ├── [3.1M]  111030_normal_genefamilies_cpm_eggnog.tsv
    │   ├── [451K]  111030_normal_genefamilies_cpm_ko.tsv
...
    │   ├── [166K]  Z116_normal_pathcoverage_relab.tsv
    │   └── [175K]  Z116_normal_pathcoverage.tsv
    ├── [ 16K]  obese
    │   ├── [2.3M]  D001_obese_genefamilies_cpm_eggnog.tsv
    │   ├── [379K]  D001_obese_genefamilies_cpm_ko.tsv
    │   ├── [808K]  D001_obese_genefamilies_cpm_level4ec.tsv
...
    │   ├── [207K]  L001_obese_pathcoverage_cpm.tsv
    │   ├── [212K]  L001_obese_pathcoverage_relab.tsv
    │   └── [221K]  L001_obese_pathcoverage.tsv
    ├── [ 28K]  ori_results
    │   ├── [3.1M]  111030_normal_genefamilies_cpm_eggnog.tsv
    │   ├── [451K]  111030_normal_genefamilies_cpm_ko.tsv
...
    │   ├── [179K]  Z116_normal_pathabundance.tsv
    │   ├── [163K]  Z116_normal_pathcoverage_cpm.tsv
    │   ├── [166K]  Z116_normal_pathcoverage_relab.tsv
    │   └── [175K]  Z116_normal_pathcoverage.tsv
    ├── [ 20K]  output
    │   ├── [ 16M]  allSamples_genefamilies_uniref90names_cpm_eggnog_stratified.tsv
...
    ├── [ 12K]  top_humann_normal_vs_obese
    │   ├── [5.8K]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.equal.top20.csv
...
    │   ├── [7.1K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.top20.csv
    │   ├── [7.6K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.top20.fillmin.scaled.csv
    │   └── [ 297]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.top20.fillmin.scaled.index
    └── [4.0K]  utest_normal_vs_obese
        ├── [957K]  allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.normal_vs_obese.ave_change.equal.csv
...
        ├── [253K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.rel_abun.unequal.csv
        └── [420K]  allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.normal_vs_obese.u-test.two_sided.csv
```

The `humann3` folder contains all intermediate tables for our humann3 process. If users want to skip the humann_init (as former suggested), they should `mkdir -p path_to_metabolism_analysis/metabolism_analysis/humann3/ori_results` and prepare all the humann3 generated tables with suffix as genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the `ori_results` folder, then set skip_humann_init to be True.

The figs folder contains all the heatmap plots for every independent database crossing every group-pair.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/allSamples_genefamilies_uniref90names_relab_pfam_unstratified.named.rel_abun_format.healthy_vs_ill.rel_abun.unequal.top20.fillmin.scaled.heatmap.pdf)
(If image failed to show, click it to view.)


### OUTPOST outputs: antibiotic genes analysis/plasmids analysis/virulence factors analysis
---
These folders contain the dist plot and heatmap plot for all features with responding taxonomy, as well as intermediate tables.
```
antibiotic_genes_analysis/
├── [250K]  cat.antibiotic.tsv
├── [ 32K]  cat.argannot.tsv
├── [ 93K]  cat.card.tsv
├── [ 614]  cat.ecoh.tsv
├── [ 63K]  cat.megares.tsv
...
├── [105K]  genes_taxa_counts_normal_vs_obese.tsv
└── [ 49K]  genes_taxaID_normal_vs_obese_heatmap.pdf
```

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/genes_asian_vs_euro_distrplot.pdf)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/genes_class_healthy_vs_ill_heatmap.pdf)
(If image failed to show, click it to view.)


OUTPOST also provide all the annotation information in `utils/abricate.annoatations.txt`, which can assist users to determine the antibiotic genes/plasmids/virulence factors. Also, the intermediate tables are helpful. The `*.antibiotic.tsv` `*.virulence.tsv` `*.plasmidfinder.tsv` are summary tables.
```
$ head antibiotic_genes_analysis/cat.antibiotic.tsv
#FILE	SEQUENCE	START	END	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION
/home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010000403.1	1526	2052	nimA_1	1-525/531	========/======	2/2	98.87	98.10	resfinder	X71444
/home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010006356.1	26667	26812	lsa(B)_1	1-146/1479	==.............	0/0	9.87	76.03	resfinder	AJ579365
/home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010006955.1	1336	1476	VanR-A_1	556-696/696	...........====	0/0	20.26	78.72	resfinder	FJ866609
/home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010007078.1	1635	1688	lnu(C)_1	442-495/495	.............==	0/0	10.91	100.00	resfinder	AY928180
...0
/home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010022089.1	8462	8940	VanY-G_1	340-818/840	......=========	0/0	57.02	75.99	resfinder	DQ212986
```


### OUTPOST outputs: biomarkers analysis
---
The `biomarkers_analysis` folder contains the biomarkers inference across all taxonomy levels for all group-pairs. The ANCOM (Analysis of compositions of microbiomes) figures are in `ANCOM_identification`.

```
biomarkers_analysis/
├── [4.0K]  ANCOM_identification
│   ├── [ 11K]  ancom_biomarker.normal_vs_obese.at_class.tsv
│   ├── [ 49K]  ancom_biomarker.normal_vs_obese.at_family.tsv
...
│   ├── [ 33K]  ancom_biomarkers.volcano.normal_vs_obese.at_genus.pdf
│   ├── [ 12K]  ancom_biomarkers.volcano.normal_vs_obese.at_order.pdf
│   ├── [9.0K]  ancom_biomarkers.volcano.normal_vs_obese.at_phylum.pdf
│   ├── [323K]  ancom_biomarkers.volcano.normal_vs_obese.at_species.pdf
│   ├── [4.8K]  ancom_biomarkers.volcano.normal_vs_obese.at_superkingdom.pdf
│   └── [193K]  ancom_biomarkers.volcano.normal_vs_obese.at_taxaID.pdf
├── [3.5K]  OUTPOST_biomarker_scores.normal_vs_obese.class.tsv
├── [ 16K]  OUTPOST_biomarker_scores.normal_vs_obese.family.tsv
...
├── [ 241]  OUTPOST_biomarker_scores.normal_vs_obese.superkingdom.tsv
└── [432K]  OUTPOST_biomarker_scores.normal_vs_obese.taxaID.tsv
```

The higher OUTPOST_biomarker_score, the more likely the taxonomy is the biomarker. Full score is 7. The following example shows that `Escherichia coli` is most likely to be biomarker species.
```
$ head biomarkers_analysis/OUTPOST_biomarker_scores.normal_vs_obese.species.tsv
	taxa_abun_top_100	taxa_abun_signicant	LDA_larger_than_2	virulence_top_100	plasmid_top_100	antibiotic_top_100	ANCOM_identified	OUTPOST_biomarker_score
Escherichia coli	1.0	1.0	1.0	1.0	1.0	1.0	1.0	7.0
Holdemanella biformis	1.0	1.0	1.0	1.0	0.0	1.0	1.0	6.0
...
uncultured bacterium	1.0	1.0	1.0	1.0	0.0	1.0	1.0	6.0
Bifidobacterium pseudolongum	1.0	1.0	1.0	1.0	0.0	1.0	0.0	5.0
Helicobacter canis	1.0	1.0	1.0	1.0	0.0	0.0	1.0	5.0
```

### OUTPOST outputs: OUTPOST report
---
The `report` contains the automatic generated reports for all group-pairs.
```
report/
└── [5.2M]  OUTPOST_report_normal_vs_obese.pdf
```
![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/OUTPOST_report_normal_vs_obese.pdf)
(If image failed to show, click it to view.)

## log
The `log` folder contains the log for all processing modules of OUTPOST. Notably, OUTPOST use these log files with suffix of `.done` to judge the status of each rule. Users can check the `Snakemake.py`, if the `done` files for certain rule exists, OUTPOST will regart it as successfully finished. If users want to re-run certain rules, they should remove the corressponding `done` files in the first place.
```
log/
├── [   0]  alpha_beta_diversity.done
...
├── [   0]  virulence_factors_analysis.done
├── [   0]  virulence_factors_analysis.log
├── [   0]  visualize_batch_effect.done
└── [   0]  visualize_batch_effect.log
```


### OUTPOST outputs: benchmark
---
The `benchmark` folder contains the log information of each computation module in OUTPOST. Users can check the benchmark files for diagnosis.
For example,
```
benchmark/
├── [ 141]  alpha_beta_diversity.benchmark
...
├── [ 143]  virulence_analysis.benchmark
├── [ 145]  virulence_factors_analysis.benchmark
└── [ 140]  visualize_batch_effect.benchmark

$cat benchmark/alpha_beta_diversity.benchmark 
s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time
14.7261	0:00:14	487.06	19615.18	446.61	461.81	5.42	3.06	186.02	29.58
```

---
# Common questions
---
### 1. When running humann3, get `No MetaPhlAn BowTie2 database found (--index option)!` error.
Basically, the cause is the wrong installation of humann3.
1. google it. Try `metaphlan —install`. Then re-run snakemake.
2. if 1. not work. check the error log, try to find this sentence `Expecting location /the/expection/location`; `cd` to the location, examine any file truncation/loss. If there be, remove all files in the `/the/expection/location`, re-run `metaphlan —install`. Then re-run snakemake.

### 2. How to deploy OUTPOST on Slurm/PBS system?
Refer to this [answer](https://stackoverflow.com/questions/53545690/how-to-activate-a-specific-python-environment-as-part-of-my-submission-to-slurm). I tried and succeeded. 

### 3. Can't locate File/Slurp.pm in @INC (@INC contains: /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor_perl /usr/share/perl5/vendor_perl /usr/lib64/perl5 /usr/share/perl5 .) at /home/yzz0191/anaconda3/envs/OUTPOST/bin/abricate line 9.
I occured this when submitting OUTPOST to Slurm system. The reason is `~/anaconda3/envs/OUTPOST/bin/abricate`, the `abricate` was written in Perl, and the first line of `abricate` is `#!/usr/bin/env perl`. So, comment this line out, add a new head line `#! /path/to/your/anaconda3/bin/perl` (modify /path/to/your !), will solve it.

### 4. When use OUTPOST.yml create conda environment, occurred Bioconductor related issues.
`conda env remove -n OUTPOST` to clean the failed OUTPOST environment. Then use `OUTPOST_without_bioconductor.yml` to create a new OUTPOST environment.

### 5. Humann error, `Please update your version of MetaPhlAn2 to v3.0`
run `conda install -c biobakery humann=3.1.1`. This is because humann internal error. Need to update to newer version of humann. After the installation, the original chocophlan database will be out of use (you will get `Please install the latest version of the database: v201901_v31` error). So run `humann_databases --download chocophlan full /path/to/your/databases/ --update-config yes` to update your chocophlan database for humann. Then, since you just updated the humann, you need to update the nucleotide and protein database record in humann_config too so that huamann can find the databases. Run `humann_config --update database_folders protein /path/to/your/uniref/database/` and `humann_config --update database_folders protein /path/to/your/utility_mapping/database/`.

### 6. The assembly looks normal, but no plasmids/antibiotic/virus table generated
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


