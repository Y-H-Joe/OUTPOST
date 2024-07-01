# OUTPOST: a comprehensive analysis software for whole-metagenome shotgun sequencing incorporating group stratification.

***

# Installation

***

    # 1. Download `OUTPOST.yml` in `install` folder.
    git clone https://github.com/Y-H-Joe/OUTPOST.git
    cd OUTPOST

    # 2. Create the conda environment `conda env create --name OUTPOST --file OUTPOST.yml `, the bioconductor has unknown errors.
    export PIP_DEFAULT_TIMEOUT=100
    mkdir -p ~/.pip
    echo -e "[global]\ntimeout = 100\n" > ~/.pip/pip.conf
    conda env create --name OUTPOST --file install/OUTPOST.yml

    # 3. Activate the environment `conda activate OUTPOST`.
    conda activate OUTPOST

    # 4. install. Databases will also be donwnloaded, require ~500GB. The downloading time is based on your bandwith.
    # make sure your working directory is OUTPOST folder
    bash install/OUTPOST_install.sh

***

# Usage

***

    # 1. prepare the `OUTPOST/OUTPOST_config.tsv` (for experiment, data, group information)
    # 2. modify the `OUTPOST/Snakemake_config.yml` (for OUTPOST parameters)
    # 3. run `python path/to/OUTPOST/check_snakefile_config.py path/to/Snakefile_config.yml` . If no errors, then you are good to go.
    # 4. run OUTPOST
    `nohup snakemake --cores 32 --verbose -s ./OUTPOST_run.py --rerun-incomplete &`

We prepared the example dataset and corresponding results in [figshare](https://figshare.com/articles/dataset/The_OUTPOST_results_for_example_dataset_and_report_/25434553), which can be used for testing your installation.

We also uploaded our independant analysis of a cat icrobiome dataset (described in our article) in [figshare](https://figshare.com/articles/figure/OUTPOST_results_of_cat_microbiome_dataset/24082542).

To Be Remind:

1.  The working directory (where you run OUTPOST) must be the parent directory of `OUTPOST`, becasue some scripts in Snakemake.py use relative path. To be more specific, you have a folder `test`, you have all scripts of OUTPOST in `test` folder, you have `test\Snakefile.py` and `test\OUTPOST`. You `cd test`, then run  `nohup snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete &`, you will get all outputs in `test\name_of_your_assembly` folder.

2.  Snakemake (the framework OUTPOST relied on) will lock working directory during running. So you should prepare two working directories if you're running two OUTPOST pipelines (or any other Snakemake based softwares). To be more specific, `mkdir folder1/` and `mkdir folder2/`, copy the entire OUTPOST folder to `folder1/` and `folder2/`, respectively. Then `cd folder1`, run OUTPOST. Then `cd folder2`, run OUTPOST.

3.  The Snakemake based OUTPOST can also be deployed on [cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) .

4.  Do not change the column names of OUTPOST\_config.tsv or the parameter names of Snakefile\_config.yml.

5.  Refer to OUTPOST/OUTPOST\_config.explanation.txt for more details.

***

# Outputs

***

We prepared the example dataset and corresponding results in [figshare](https://figshare.com/articles/dataset/The_OUTPOST_results_for_example_dataset_and_report_/25434553), which can be used for testing your installation.

We also uploaded our independant analysis of a cat icrobiome dataset (described in our article) in [figshare](https://figshare.com/articles/figure/OUTPOST_results_of_cat_microbiome_dataset/24082542) using OUTPOST version 1.0 .

We strongly suggest to read the OUTPOST\_report.html in `OUTPOST_report_example.zip`.

Each OUTPOST run accepts one assembly, all outputs are categorized in the folder named by the assembly.

    antibiotic_genes_analysis
    assembly_analysis
    batch_effect
    benchmark
    biomarkers_analysis
    data
    diversity_analysis
    function_analysis
    LDA_analysis
    log
    metaphlan_analysis
    plasmids_analysis
    qc
    report
    taxonomy_analysis
    virulence_factors_analysis

### 1. OUTPOST reports

***

    [   4]  report/
    ├── [  14]  materials
    └── [ 38K]  OUTPOST_report.html

The [OUTPOST report](https://github.com/Y-H-Joe/OUTPOST/blob/main/OUTPOST_report_example.zip) is a HTML file with materials, which contains the general results description and annoatation.

For more annotation, please refer to OUTPOST report html.

### 2. log

***

    log/
    ├── [   0]  alpha_beta_diversity.done
    ...
    ├── [   0]  virulence_factors_analysis.done
    ├── [   0]  virulence_factors_analysis.log
    ├── [   0]  visualize_batch_effect.done
    └── [   0]  visualize_batch_effect.log

The `log` folder contains the log for all processing modules of OUTPOST. Notably, OUTPOST use these log files with suffix of `.done` to judge the status of each rule. Users can check the `Snakemake.py`, if the `done` files for certain rule exists, OUTPOST will regart it as successfully finished. If users want to re-run certain rules, they should remove the corressponding `done` files in the first place.

### 3. benchmark

***

    benchmark/
    ├── [ 141]  alpha_beta_diversity.benchmark
    ...
    ├── [ 143]  virulence_analysis.benchmark
    ├── [ 145]  virulence_factors_analysis.benchmark
    └── [ 140]  visualize_batch_effect.benchmark

    $cat benchmark/alpha_beta_diversity.benchmark 
    s	h:m:s	max_rss	max_vms	max_uss	max_pss	io_in	io_out	mean_load	cpu_time
    14.7261	0:00:14	487.06	19615.18	446.61	461.81	5.42	3.06	186.02	29.58

The `benchmark` folder contains the log information of each computation module in OUTPOST. Users can check the benchmark files for diagnosis. For example,

For more annotation, please refer to OUTPOST report html.

### 4. quality control

***

    [  56]  qc/
    ├── [234K]  buffalo_1.fastp.html
    ...
    ├── [ 51K]  pig_3.fastp.json
    ├── [   0]  TrimSummmry_buffalo_1.txt
    ...
    ├── [   0]  TrimSummmry_pig_1.txt
    ├── [   0]  TrimSummmry_pig_2.txt
    └── [   0]  TrimSummmry_pig_3.txt

The `qc` folder contains quality control results for reads. For example,

### 5. batch effect

***

    [  10]  batch_effect/
    ├── [ 32K]  all_samples.taxa_counts.rel_abun.class.rmU.batch_effect_PCA.pdf
    ...
    ├── [ 25K]  all_samples.taxa_counts.rel_abun.species.rmU.batch_effect_PCA.pdf
    ├── [ 32K]  all_samples.taxa_counts.rel_abun.superkingdom.rmU.batch_effect_PCA.pdf
    └── [ 25K]  all_samples.taxa_counts.rel_abun.taxaID.rmU.batch_effect_PCA.pdf

If `rm_batch_effect` is `True`, OUTPOST will visualize the principal components (PCA) as well as the variance distribution. To check the batch effect, users can compare the plots before and after batch effect removal. For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.family.rmU.batch_effect_PCA.png) (If image failed to show, click it to view.)

For more annotation, please refer to OUTPOST report html.

### 6. taxonomy analysis

***

    taxonomy_analysis/
    ├── [628K]  boxplot_normal_vs_obese
    │   ├── [4.0K]  class
    │   ├── [ 28K]  family
    ...
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

The `boxplot` folder contains the boxplots for all significant taxa crossing all taxonomy levels for all group-pair comparisons. Also, all plots generated by OUTPOST are vector PDF format for the convenience of publication.

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.rel_abun.healthy_vs_ill.at_species.rel_abun.unequal.Alistipes.finegoldii.CAG.68.boxplot.png) (If image failed to show, click it to view.)

The `figs` folder contains the heatmap and barplots: OUTPOST not only draws heatmap for unequal taxa, but also for equal taxa.

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_family.rel_abun.equal.top20.fillmin.scaled.heatmap.png) (If image failed to show, click it to view.) OUTPOST produces barplot for all taxonomy levels. Here we have top 20 taxa, the 20 here is an adjustable parameter.

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.taxa_counts.rel_abun.phylum.rmU.top20.barplot.png) (If image failed to show, click it to view.)

The `csv` files here are intermediate tables. OUTPOST include these tables for user's convenience. The `kaiju` folder contains the taxonomy annotated contigs table: The `top_taxa` folders contain the intermediate tables for top abundant taxa crossing all taxonomy levels for every group pair. The `utes`t folders contain the statistical results.

For more annotation, please refer to OUTPOST report html.

### 7. diversity analysis

***

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

This folder contains all the alpha and beta diversity analysis crossing every group-pair.

    (base) yh@superServer:human62_batch_effect2$ l diversity_analysis/
    alpha_beta_asian_vs_euro/     alpha_beta_asian_vs_ill/   alpha_beta_female_vs_euro/   alpha_beta_healthy_vs_female/  alpha_beta_male_vs_euro/
    alpha_beta_asian_vs_female/   alpha_beta_asian_vs_male/  alpha_beta_female_vs_ill/    alpha_beta_healthy_vs_ill/     alpha_beta_male_vs_female/
    alpha_beta_asian_vs_healthy/  alpha_beta_euro_vs_ill/    alpha_beta_healthy_vs_euro/  alpha_beta_healthy_vs_male/    alpha_beta_male_vs_ill/

Each sub-folder contains alpha and beta diversity analysis intermediate tables and graphs.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/PCoA12.bray.at_species.png)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/Shannon.alpha_diveristy.at_.species.png) (If image failed to show, click it to view.)

For more annotation, please refer to OUTPOST report html.

### 8. MetaPhlAn4 analysis

***

    [   7]  metaphlan_analysis/
    ├── [  20]  bowtie2
    ├── [  11]  diversity
    ├── [   5]  figs
    ├── [  20]  krona
    └── [  24]  taxonomy

This folder contains all the analysis results using MetaPhlan4, including some results specially generated by OUTPOST. This folder include diversity analysis, taxonomy analysis, and more.

OUTPOST generated GraPhlAn plot, heatmap and Krona Plot for each Sample. The tables for taxonomy, Krona and Bowtie2 intermediates are also stored.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/graphlan.png) (If image failed to show, click it to view.)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/metaphlan_merge_taxa.heatmap.png) (If image failed to show, click it to view.)

For more annotation, please refer to OUTPOST report html.

### 9. function analysis

***

The folder contains the metabolism analysis results, mainly based on the 5 independent databases (which are selectable in `Snakemake_config.yml` file)

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

The `humann3` folder contains all intermediate tables for our humann3 process. If users want to skip the humann\_init (as former suggested), they should `mkdir -p path_to_metabolism_analysis/metabolism_analysis/humann3/ori_results` and prepare all the humann3 generated tables with suffix as genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the `ori_results` folder, then set skip\_humann\_init to be True.

The figs folder contains all the heatmap plots for every independent database crossing every group-pair.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/allSamples_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.omnivorous_vs_dog.rel_abun.equal.top20.fillmin.scaled.heatmap.png) (If image failed to show, click it to view.)

For more annotation, please refer to OUTPOST report html.

### 10. LDA analysis

***

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

This folder contains the metabolism and taxonomy LDA results for every group-pair.

The `figs_humann` and `figs_taxa` folders contain the lefSE results.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/allSamples_genefamilies_uniref90names_relab_eggnog_unstratified.named.rel_abun_format.asian_vs_euro.rel_abun.unequal.lefse.png)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/human62_batch_effect2.rel_abun.asian_vs_euro.at_genus.rel_abun.unequal.lefse.png) (If image failed to show, click it to view.)

Other sub-folders contain the intermediate tables for user's convenience.

For more annotation, please refer to OUTPOST report html.

### 11. antibiotic genes analysis/plasmids analysis/virulence factors analysis

***

    antibiotic_genes_analysis/
    ├── [250K]  cat.antibiotic.tsv
    ├── [ 32K]  cat.argannot.tsv
    ├── [ 93K]  cat.card.tsv
    ├── [ 614]  cat.ecoh.tsv
    ├── [ 63K]  cat.megares.tsv
    ...
    ├── [105K]  genes_taxa_counts_normal_vs_obese.tsv
    └── [ 49K]  genes_taxaID_normal_vs_obese_heatmap.pdf

These folders contain the dist plot and heatmap plot for all features with responding taxonomy, as well as intermediate tables.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/genes_asian_vs_euro_distrplot.png)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/genes_class_healthy_vs_ill_heatmap.png) (If image failed to show, click it to view.)

OUTPOST also provide all the annotation information in `utils/abricate.annoatations.txt`, which can assist users to determine the antibiotic genes/plasmids/virulence factors. Also, the intermediate tables are helpful. The `*.antibiotic.tsv` `*.virulence.tsv` `*.plasmidfinder.tsv` are summary tables.

    $ head antibiotic_genes_analysis/cat.antibiotic.tsv
    #FILE	SEQUENCE	START	END	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION
    /home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010000403.1	1526	2052	nimA_1	1-525/531	========/======	2/2	98.87	98.10	resfinder	X71444
    /home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010006356.1	26667	26812	lsa(B)_1	1-146/1479	==.............	0/0	9.87	76.03	resfinder	AJ579365
    /home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010006955.1	1336	1476	VanR-A_1	556-696/696	...........====	0/0	20.26	78.72	resfinder	FJ866609
    /home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010007078.1	1635	1688	lnu(C)_1	442-495/495	.............==	0/0	10.91	100.00	resfinder	AY928180
    ...0
    /home/yihang/refgenomes/cat_microbiome/GCA_022675345.1_ASM2267534v1_genomic.fna	JAIZPE010022089.1	8462	8940	VanY-G_1	340-818/840	......=========	0/0	57.02	75.99	resfinder	DQ212986

For more annotation, please refer to OUTPOST report html.

### 12. assembly analysis

***

        [   8]  assembly_analysis/
    ├── [   7]  annotation
    │   ├── [   4]  eggmapper
    │   ├── [  10]  gtdbtk
    │   ├── [   7]  metagenemark
    │   ├── [   9]  prodigal
    │   └── [   6]  rgi
    ├── [   8]  contigs
    │   ├── [ 262]  checkpoints.txt
    │   ├── [   0]  done
    │   ├── [ 15M]  final.contigs.fa
    │   ├── [  80]  intermediate_contigs
    │   ├── [184K]  log
    │   └── [2.9K]  options.json
    ├── [  11]  outpost_contigs
    │   ├── [ 23K]  outpost_assembly_split
    │   ├── [ 15M]  outpost_nonrd_contigs.fasta
    │   ├── [ 28M]  outpost_nonrd_contigs.fasta.0123
    │   ├── [  17]  outpost_nonrd_contigs.fasta.amb
    │   ├── [1.2M]  outpost_nonrd_contigs.fasta.ann
    │   ├── [ 45M]  outpost_nonrd_contigs.fasta.bwt.2bit.64
    │   ├── [940K]  outpost_nonrd_contigs.fasta.clstr
    │   ├── [3.5M]  outpost_nonrd_contigs.fasta.pac
    │   └── [   2]  stdin.split
    ├── [   4]  quantify
    │   ├── [   4]  index
    │   └── [   4]  quantify
    ├── [  15]  quast
    │   ├── [   6]  basic_stats
    │   ├── [ 52K]  icarus.html
    │   ├── [   3]  icarus_viewers
    │   ├── [   6]  predicted_genes
    │   ├── [4.4K]  quast.log
    │   ├── [418K]  report.html
    │   ├── [ 31K]  report.pdf
    │   ├── [1.5K]  report.tex
    │   ├── [ 778]  report.tsv
    │   ├── [1.7K]  report.txt
    │   ├── [1.3K]  transposed_report.tex
    │   ├── [ 778]  transposed_report.tsv
    │   └── [1.3K]  transposed_report.txt
    └── [   4]  taxa_counts
        ├── [1.0K]  all_samples.taxa_counts.summarized.tsv
        └── [4.4M]  all_samples.taxa_counts.tsv

This folder contains assemble contigs, final assembly, assembly quality control, assembly genes prediction, assembly genes annotation, assembly genes quantification, and assembly genes taxonomy counts.

For more annotation, please refer to OUTPOST report html.

### 13. biomarkers analysis

***

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

The higher OUTPOST\_biomarker\_score, the more likely the taxonomy is the biomarker. Full score is 7. The following example shows that `Escherichia coli` is most likely to be biomarker species.

    $ head biomarkers_analysis/OUTPOST_biomarker_scores.normal_vs_obese.species.tsv
    	taxa_abun_top_100	taxa_abun_signicant	LDA_larger_than_2	virulence_top_100	plasmid_top_100	antibiotic_top_100	ANCOM_identified	OUTPOST_biomarker_score
    Escherichia coli	1.0	1.0	1.0	1.0	1.0	1.0	1.0	7.0
    Holdemanella biformis	1.0	1.0	1.0	1.0	0.0	1.0	1.0	6.0
    ...
    uncultured bacterium	1.0	1.0	1.0	1.0	0.0	1.0	1.0	6.0
    Bifidobacterium pseudolongum	1.0	1.0	1.0	1.0	0.0	1.0	0.0	5.0
    Helicobacter canis	1.0	1.0	1.0	1.0	0.0	0.0	1.0	5.0

The `biomarkers_analysis` folder contains the biomarkers inference across all taxonomy levels for all group-pairs. The ANCOM (Analysis of compositions of microbiomes) figures are in `ANCOM_identification`.

For example,

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/ancom_biomarkers.dotplot.pig_vs_horse.at_order.png) (If image failed to show, click it to view.)

![image](https://github.com/Y-H-Joe/OUTPOST/blob/main/figs/ancom_biomarkers.volcano.omnivorous_vs_pig.at_phylum.png) (If image failed to show, click it to view.)

For more annotation, please refer to OUTPOST report html.

***

# Common questions

***

### 1. When running humann3, get `No MetaPhlAn BowTie2 database found (--index option)!` error.

Basically, the cause is the wrong installation of humann3.

1.  google it. Try `metaphlan —install`. Then re-run snakemake.

2.  if 1. not work. check the error log, try to find this sentence `Expecting location /the/expection/location`; `cd` to the location, examine any file truncation/loss. If there be, remove all files in the `/the/expection/location`, re-run `metaphlan —install`. Then re-run snakemake.

### 2. How to deploy OUTPOST on Slurm/PBS system?

Refer to this [answer](https://stackoverflow.com/questions/53545690/how-to-activate-a-specific-python-environment-as-part-of-my-submission-to-slurm). I tried and succeeded.

### 3. Can't locate File/Slurp.pm in @INC (@INC contains: /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor\_perl /usr/share/perl5/vendor\_perl /usr/lib64/perl5 /usr/share/perl5 .) at /home/yzz0191/anaconda3/envs/OUTPOST/bin/abricate line 9.

I occured this when submitting OUTPOST to Slurm system. The reason is `~/anaconda3/envs/OUTPOST/bin/abricate`, the `abricate` was written in Perl, and the first line of `abricate` is `#!/usr/bin/env perl`. So, comment this line out, add a new head line `#! /path/to/your/anaconda3/bin/perl` (modify /path/to/your !), will solve it.

### 4. When use OUTPOST.yml create conda environment, occurred Bioconductor related issues.

`conda env remove -n OUTPOST` to clean the failed OUTPOST environment. Then use `OUTPOST_without_bioconductor.yml` to create a new OUTPOST environment.

### 5. Humann error, `Please update your version of MetaPhlAn2 to v3.0`

run `conda install -c biobakery humann=3.1.1`. This is because humann internal error. Need to update to newer version of humann. After the installation, the original chocophlan database will be out of use (you will get `Please install the latest version of the database: v201901_v31` error). So run `humann_databases --download chocophlan full /path/to/your/databases/ --update-config yes` to update your chocophlan database for humann. Then, since you just updated the humann, you need to update the nucleotide and protein database record in humann\_config too so that huamann can find the databases. Run `humann_config --update database_folders protein /path/to/your/uniref/database/` and `humann_config --update database_folders protein /path/to/your/utility_mapping/database/`.

### 6. The assembly looks normal, but no plasmids/antibiotic/virus table generated

Reformat your assembly file. from `>human63_000000000001 TTTCCTTCGATGAGTTCTATGCCGTATATAATAAAAAGCATTCCGCTATTGAACAGCGTCTCGCAGAAAAAGGATTGCCGGAACATCTGCTTCATCGTAAGGAACGCAGACAGGAAAAACTGAATCATCCTGCTGTAAAAACGACAAAGCCCCACAGAAAGAAGAAAAAGAAACAGGTGTTCGAGCCGCTCTTGGAACAGAATGATGATTTCTTCTTTATTGCTGGTTATACTTCTGGCGGTGCCCCTTATGGTGTCACATGGGAAGAAATGGGACTAGAGCCTTGGGAAGAACTTGTATAAATATTATTGCCATTGCCGATTGCCAAAAGCA` to

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

### 7. Skip humann and assembly

Humann analysis is very time/computation consuption. You can use the `downsample_reads` parameters, or you can set `skip_humann_init` if you already have the init resutls. Set `skip_humann_init = True` in `Snakefile_config.yml`, then OUTPOST will not run human\_init rule, but to check the human results under folder `name_of_the_assembly/metabolism_analysis/humann3/ori_results/`, so you need to put the humann output with suffix as genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the folder.

OUTPOST offers `skip_assembly_analysis`. You can skip the assembly analysis if you have a large assembly and a number of groups which will save a lot of time about MAG tables generation.
