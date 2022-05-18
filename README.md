# GEMINI
An auto-comparative, flexible, light and fast, downstream pipeline for metagenomics whole genome sequencing

# GEMINI installation
0. GEMINI only supports Illumina reads, not for Ion Torrent yet.
1. download `GEMINI.yml`
2. create the conda environment `conda env create --name GEMINI --file GEMINI.yml `
3. activate the environment `conda activate GEMINI`
4. check you're using right humann `which humann`
5. download 1st humann databases `humann_databases --download chocophlan full /path/to/databases --update-config yes`
6. download 2nd humann databases `humann_databases --download uniref uniref90_diamond /path/to/databases --update-config yes`
7. download 3rd humann databases `humann_databases --download utility_mapping full /path/to/databases --update-config yes`
8. check your humann databases `cd /path/to/databases` then run `ll chocophlan/ | wc -l` you get a number >= 11289. run `ll uniref` you should see a `uniref90_201901b_full.dmnd` (or newer) file with >= 34G size. run `ll utility_mapping` you should see >= 21 files and all of them have > 3M size (or some of them truncated during download).
9. check your humann's misc folder. Located at `/path/to/your/anaconda/envs/GEMINI/lib/python3.9/site-packages/humann/data/misc/`. Due to conda's unknown error, usually the files are missing. To know full file list, check README.txt in the misc folder. To download the files, check https://github.com/biobakery/humann/tree/master/humann/data/misc . 
10. check your humann by running `humann -i sample_reads.fastq -o sample_results` (prepare sample_reads.fastq by yourself)
11. check you're using right kaiju `which kaiju`
12. download kaiju databases `mkdir /path/to/kaijudb` then `cd /path/to/kaijudb` then `kaiju-makedb -s nr_euk` (this takes a long time and space and memory)
13. open `Snakefile.py`, modify the bwa,kaiju,python3,Rscript,...lefse_run parameters to the executable command lines in your environment. To make sure all command line works, please test the command line one by one in your linux shell.
14. test GEMINI. `cd parent/folder/of/GEMINI`. open and modify the `GEMINI/contig.tsv` to make sure the data_dir is right. then run `snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete`. The test will last 3 hours.
15. if you occured any errors. check the printed log to debug. or check the log file in `name_of_your_assembly/log` folder. You can use time stamps to refer which rule is error, or to understand the error information. After debugging, delete the `name_of_your_assembly/log/name_of_the_error_rule.done`. and rerun the `snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete`. GEMINI will auto-resume.
16. You should see the outputs in `horsedonkey` folder. If no errors occured. then you're good to go.

# format of GEMINI_config.tsv
| samples | fq_dir                                                       | bam_dir                                                  | assembly    | assembly_dir                                                 | group               |   |   |
|---------|--------------------------------------------------------------|----------------------------------------------------------|-------------|--------------------------------------------------------------|---------------------|---|---|
| donkey1 | /analysis1/yihang_analysis/pipeline/data/reads/donkey1.fq.gz | /analysis1/yihang_analysis/pipeline/data/bam/donkey1.bam | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | donkey,horsedonkey  |   |   |
| hinny2  | /analysis1/yihang_analysis/pipeline/data/reads/hinny2.fq.gz  | /analysis1/yihang_analysis/pipeline/data/bam/hinny2.bam  | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | hinny               |   |   |
| hinny3  | /analysis1/yihang_analysis/pipeline/data/reads/hinny3.fq.gz  | /analysis1/yihang_analysis/pipeline/data/bam/hinny3.bam  | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | hinny               |   |   |
| horse1  | /analysis1/yihang_analysis/pipeline/data/reads/horse1.fq.gz  | /analysis1/yihang_analysis/pipeline/data/bam/horse1.bam  | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | horse,horsedonkey   |   |   |
| horse3  | /analysis1/yihang_analysis/pipeline/data/reads/horse3.fq.gz  | /analysis1/yihang_analysis/pipeline/data/bam/horse3.bam  | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | horse,horsedonkey   |   |   |
| donkey2 | /analysis1/yihang_analysis/pipeline/data/reads/donkey2.fq.gz | /analysis1/yihang_analysis/pipeline/data/bam/donkey2.bam | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | donkey,horsedonkey  |   |   |
| donkey3 | /analysis1/yihang_analysis/pipeline/data/reads/donkey3.fq.gz | /analysis1/yihang_analysis/pipeline/data/bam/donkey3.bam | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | donkey,horsedonkey  |   |   |
| horse2  | /analysis1/yihang_analysis/pipeline/data/reads/horse2.fq.gz  | /analysis1/yihang_analysis/pipeline/data/bam/horse2.bam  | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | horse,horsedonkey   |   |   |
| hinny1  | /analysis1/yihang_analysis/pipeline/data/reads/hinny1.fq.gz  | /analysis1/yihang_analysis/pipeline/data/bam/hinny1.bam  | horsedonkey | /analysis1/yihang_analysis/pipeline/data/assembly/sample1.fa | hinny               |   |   |


The ***samples*** column contains the ids of each sample. Each sample id should be unique. The id can contain special characters like '_', but no ',' allowed.

The ***fq_dir*** column contains the absolute location of each fastq file. If you use PE end sequencing, `cat` them to one file, because humann3 asked so. 

The ***bam_dir*** column contains the absolute directory of each ***bam***. The bam file should be generated by users, via aligning reads (both SE and PE mapping is acceptable) to assembly. 

The ***assembly*** column contains the name of 1st assembly. 

The ***assembly_dir*** column contains the absolute directory of ***assembly***.

The ***group*** column contains the group name of the each sample. If one sample belong to multiple groups, seperate by ','. Each group should at least have 3 samples for statistical power concern.

***do not change the name of columns***

In the above case. We have one reference, hmdh. For groups, we have 'donkey', 'horse', 'horsedonkey', 'hinny'.

The pipeline will pair-wise compare all of them. Which are, 
['donkey_vs_horse','donkey_vs_horsedonkey','donkey_vs_hinny','horse_vs_horsedonky','horse_vs_hinny','horsedonkey_vs_hinny'].

# GEMINI usuage
The usuage of GEMINI is very easy and light. If you successfully installed GEMINI and got all the expected outputs files in `horsedonkey` folder, you just need to replace the sample data with your own data, and modify the `GEMINI/config.tsv`.

I just list some warnings.

1. The working directory must be the parent directory of `GEMINI`, becasue the command line in Snakemake.py contains many sentences like `python3 GEMINI/python_script.py`.
2. Snakemake (the framework GEMINI relied on) will lock working directory during running. So you should prepare two working directories if you're running two GEMINI pipelines (or any other Snakemake based softwares). To be more specific, `mkdir folder1/GEMINI` and `mkdir folder2/GEMINI`, copy the `Snakemake.py` to `folder1/Snakemake.py` and `folder2/Snakemake.py`. Then `cd folder1`, run GEMINI. Then `cd folder2`, run GEMINI.
3. One GEMINI process only take one assembly. If you have multi assemblies to analyze, run GEMINI multiple times (fake parallel) on multiple computing nodes (in different working directories)
4. Humann analysis is very time/computation consuption. I personally prefer to use computer cluster to distributedly run Humann. The Snakemake based GEMINI can also be deployed on cluster, but maybe not [easy](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). So GEMINI provide `skip_humann_init` option. Set `skip_humann_init = True` in `GEMINI_config.tsv`, then GEMINI will not run human_init rule, but to check the human results under folder `name_of_the_assembly/metabolism_analysis/humann3/ori_results/`, so you need to put the humann output *genefamilies.tsv/pathabundance.tsv/pathcoverage.tsv under the folder. If `skip_humann_init = True`, GEMINI will check the outputs existence first. Make sure these humann outptus are from the same fastq you offered to GEMINI.
5. GEMINI offers `skip_assembly_analysis` in `GEMINI_config.tsv`. 

# GEMINI outpus
## assembly analysis
1. `horsedonkey.taxa_counts.tsv`: contig_ID  NCBI_taxa_ID  superkingdom  phylum  class order family  genus species number_of_sample1_reads_aligned_to_this_contig  number_of_sample2_reads_aligned_to_this_contig...
2. `horsedonkey.horsedonkey_vs_hinny.contig_table.tsv`: similar as `horsedonkey.taxa_counts.tsv` besides q value and p value (group horsedonkey v.s. group hinny)
3. `horsedonkey.horsedonkey_vs_hinny.contig_table.processed.tsv`: similar as `horsedonkey.horsedonkey_vs_hinny.contig_table.tsv` besides relative abundance and counts per million reads for each sample in each contig.

## taxa analysis
`horsedonkey.taxa_counts.rel_abun.phylum.rmU.top30.barplot.pdf`: the barplot of all samples, only top 30 abundant phlyum were labeled.
![image](https://user-images.githubusercontent.com/28485220/166333192-718409d2-3fce-47cf-b0ae-682bccdbf086.png)

`horsedonkey.taxa_counts.rel_abun.phylum.rmU.top30.fillmin.log10.heatmap.pdf`: the top 30 abundant phlyum 
![image](https://user-images.githubusercontent.com/28485220/166333440-0a4163a3-baac-4f97-832a-2c9dd0bbf6d7.png)

`horsedonkey.rel_abun.horsedonkey_vs_hinny.at_species.rel_abun.equal.top30.fillmin.log10.heatmap.pdf`: the top 30 equal species comparing horsedonkey and hinny. Missing values (relative abundance) were filled using global minimum value. All values were log10 transfered.

![image](https://user-images.githubusercontent.com/28485220/166330581-d9906ab6-332f-423d-9cb9-49fdb2edb00f.png)
`horsedonkey.rel_abun.horsedonkey_vs_hinny.at_species.rel_abun.unequal.top30.fillmin.log10.heatmap.pdf`: the top 30 unequal species comparing horsedonkey and hinny. Missing values (relative abundance) were filled using global minimum value. All values were log10 transfered.

![image](https://user-images.githubusercontent.com/28485220/166333038-9a1153f0-0aea-4344-9206-ceced0e9d47e.png)

`horsedonkey.rel_abun.horsedonkey_vs_hinny.at_species.rel_abun.unequal.Bacteroides.sp..CAG.144.boxplot.pdf`: the boxplot of Bacteroides.sp..CAG.144 showing statistic significant different 
![image](https://user-images.githubusercontent.com/28485220/166333863-731a9747-b775-4fb3-a585-afdbab300883.png)

## metabolism analysis
## alpha/beta diversity analysis
## lefse analysis


# common questions
## When running humann3, get `No MetaPhlAn BowTie2 database found (--index option)!` error.
Basically, the cause is the wrong installation of humann3.
1. google it. Try `metaphlan —install`. Then re-run snakemake.
2. if 1. not work. check the error log, try to find this sentence `Expecting location /the/expection/location`; `cd` to the location, examine any file truncation/loss. If there be, remove all files in the `/the/expection/location`, re-run `metaphlan —install`. Then re-run snakemake.
