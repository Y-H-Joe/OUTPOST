# 0. OUTPOST only supports Illumina reads, not for Ion Torrent yet.

# 1. Download `OUTPOST.yml` in `install` folder.
# git https://github.com/Y-H-Joe/OUTPOST.git
# cd OUTPOST

# 2. Create the conda environment `conda env create --name OUTPOST --file OUTPOST.yml `, the bioconductor has unknown errors.
# conda env create --name OUTPOST --file install/OUTPOST_without_bioconductor.yml

# 3. Activate the environment `conda activate OUTPOST`.
# conda activate OUTPOST
mkdir -p install/log/
# 4. Check you're using right humann `which humann`.
humann_path=$(which humann)
if [ -z "$humann_path" ]; then
    echo "humann not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# 5. Download 1st humann databases `humann_databases --download chocophlan full /path/to/databases --update-config yes`
# if failed, run following:
# wget -c http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz --no-check-certificate
# humann_databases --download chocophlan full  $(pwd)/databases --database-location ./full_chocophlan.v201901_v31.tar.gz  --update-config yes
# touch install/log/humann_databases.utility_mapping.done
DONE_FILE="install/log/humann_databases.chocophlan.done"
if [ ! -f "$DONE_FILE" ]; then
	mkdir -p databases
	humann_databases --download chocophlan full  $(pwd)/databases --update-config yes
	touch "$DONE_FILE"
else
    echo "humann_databases.chocophlan already installed. Skipping this step."
fi

# 6. Download 2nd humann databases `humann_databases --download uniref uniref90_diamond /path/to/databases --update-config yes`
# if failed, run following:
# wget -c http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz --no-check-certificate
# 
# touch install/log/humann_databases.utility_mapping.done
DONE_FILE="install/log/humann_databases.uniref.done"
if [ ! -f "$DONE_FILE" ]; then
	humann_databases --download uniref uniref90_diamond  $(pwd)/databases --update-config yes
	touch "$DONE_FILE"
else
    echo "humann_databases.uniref already installed. Skipping this step."
fi
# 7. Download 3rd humann databases `humann_databases --download utility_mapping full /path/to/databases --update-config yes`
# if failed, run following:
# wget -c http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz --no-check-certificate
# humann_databases --download utility_mapping full $(pwd)/databases --database-location ./full_mapping_v201901b.tar.gz --update-config yes
# touch install/log/humann_databases.utility_mapping.done
DONE_FILE="install/log/humann_databases.utility_mapping.done"
if [ ! -f "$DONE_FILE" ]; then
	humann_databases --download utility_mapping full $(pwd)/databases --update-config yes
	touch "$DONE_FILE"
else
    echo "humann_databases.utility_mapping already installed. Skipping this step."
fi
# 8. Check your humann databases `cd /path/to/databases` then run `ll chocophlan/ | wc -l` you get a number >= 11289. run `ll uniref` you should see a `uniref90_201901b_full.dmnd` (or newer) file with >= 34G size. Run `ll utility_mapping` you should see >= 21 files and all of them have > 3M size (or some of them truncated during download).
ls -l $(pwd)/databases/chocophlan | wc -l
ls -l $(pwd)/databases/uniref
ls -l $(pwd)/databases/utility_mapping

# 9. Check your humann's misc folder. Located at `/path/to/your/anaconda/envs/OUTPOST/lib/python3.9/site-packages/humann/data/misc/`. Due to conda's unknown error, usually the files are missing. To know full file list, check README.txt in the [`misc`](https://github.com/biobakery/humann/tree/master/humann/data/misc) folder. Or unzip the `misc.zip` in the `utils` folder.
DONE_FILE="install/log/misc.done"
if [ ! -f "$DONE_FILE" ]; then
	python_version=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
	misc_path=$(dirname $(dirname $humann_path))/lib/python${python_version}/site-packages/humann/data/misc/
	cp utils/misc/* $misc_path/
	touch "$DONE_FILE"
else
    echo "misc already updated. Skipping this step."
fi

# 10. Check your humann by running `humann -i sample_reads.fastq -o sample_results` (prepare sample_reads.fastq by yourself)
DONE_FILE="install/log/humann_test.done"
if [ ! -f "$DONE_FILE" ]; then
	humann_test
	touch "$DONE_FILE"
else
    echo "humann_test already finished. Skipping this step."
fi

# 11. Check you're using right kaiju `which kaiju`
kaiju_path=$(which kaiju)
if [ -z "$kaiju_path" ]; then
    echo "kaiju not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# 12. Download kaiju databases `mkdir /path/to/kaijudb` then `cd /path/to/kaijudb` then `kaiju-makedb -s nr_euk` (this takes a long time and space and memory). Or you can download and unzip the annotation files from [kaiju server](https://kaiju.binf.ku.dk/server).
DONE_FILE="install/log/kaiju.done"
if [ ! -f "$DONE_FILE" ]; then
	mkdir -p databases/kaijudb
	wget -c -O databases/kaijudb/kaiju_db_nr_euk.tgz https://kaiju-idx.s3.eu-central-1.amazonaws.com/2022/kaiju_db_nr_euk_2022-03-10.tgz databases/kaijudb/
	tar -xzf databases/kaijudb/kaiju_db_nr_euk.tgz -C databases/kaijudb/
	touch "$DONE_FILE"
else
    echo "humann_test already finished. Skipping this step."
fi

# 13. Check abricate databases `abricate --list`, you should see 9 databases (argannot,card,ecoh,ecoli_vf,megares,ncbi,plasmidfinder,resfinder,vfdb). If you didn't, go to download abricate [databases](https://github.com/tseemann/abricate/tree/master/db), or use the `db.zip` file in `utils` folder, unzip them under your `db` folder. The location of `db` folder can be seen by running `abricate --help`, see the `--datadir` line. Then `cd` to the `db` folder, run `abricate --setupdb`.
abricate_path=$(which abricate)
if [ -z "$abricate_path" ]; then
    echo "abricate not found. Please ensure it is installed and in your PATH."
    exit 1
fi
DONE_FILE="install/log/abricate.done"
if [ ! -f "$DONE_FILE" ]; then
	abricate_db_path=$(dirname $(dirname $abricate_path))/db
	unzip -o utils/db.zip -d $abricate_db_path
	current_path=$(pwd)
	cd $abricate_db_path
	abricate --setupdb
	abricate --help
	cd $current_path
	touch "$DONE_FILE"
else
    echo "abricate already updated. Skipping this step."
fi

# 14. Install [`Metagenemark`](http://exon.gatech.edu/license_download.cgi). The `utils` folder also provides a copy. Note to put the key file at `~/.gm_key`.
DONE_FILE="install/log/metagenemark.done"
if [ ! -f "$DONE_FILE" ]; then
	tar -xzf utils/MetaGeneMark_linux_64.tar.gz -C utils/
	gunzip -c utils/gm_key_64.gz > ~/.gm_key
	touch "$DONE_FILE"
else
    echo "MetaGeneMark already installed. Skipping this step."
fi

# 15. Manually install R packages. Very easy. 
# 15.1 Type: R 
# 15.2 Type: options(repos = c(CRAN = "https://cloud.r-project.org"))
# 15.3 Type: if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 15.4 Type: BiocManager::install(c('mixOmics', 'sva','DESeq2'), update = FALSE)
# 15.5 Type: install.packages("remotes")
# 15.6 Type: remotes::install_github("kongdd/Ipaper", upgrade = "never"))
# 15.7 When R concle asks you whether to update other packages, choose `none`. After installation, type in `library(Ipaper)`, if no error occurs, then you're good to contine. Install [`fastANCOM`](https://github.com/ZRChao/fastANCOM), which is easy by local installation.
DONE_FILE="install/log/R.done"
if [ ! -f "$DONE_FILE" ]; then
	Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("mixOmics", "sva", "DESeq2"), update = FALSE); install.packages("remotes"); remotes::install_github("kongdd/Ipaper", upgrade = "never"); install.packages("utils/fastANCOM-main", repos = NULL, type = "source")'
	touch "$DONE_FILE"
else
    echo "R packages already installed. Skipping this step."
fi

# 16. Open `Snakefile_config.yml` under `OUTPOST` folder, modify the bwa,kaiju,python3,Rscript,cd-hit-est...lefse_run parameters to the executable command lines in your environment. To make sure all command line works, you could test the command line one by one in your linux shell.

# 17. Test OUTPOST. `cd parent/folder/of/OUTPOST`. Prepare some example data. Open and modify the `OUTPOST/OUTPOST_contig.tsv` to make sure the data_dir is right.Then run `nohup snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete &`. This command will run OUTPOST in backend with a log file `nohup.out`.

# 18. If you occured any errors. check the log to debug. Or check the log file in `name_of_your_assembly/log` folder. You can use time stamps to refer which rule is error, or to understand the error information. After debugging, delete the `name_of_your_assembly/log/name_of_the_error_rule.done`. and rerun the `snakemake --cores 32 --verbose -s ./Snakefile.py --rerun-incomplete`. OUTPOST will auto-resume.

# 19. You should see all the outputs in `name_of_your_assembly` folder. If no errors occured. then you're good to go.

# 20. In a summary, OUTPOST itself is just a list of scripts, which is easy to use. The above installation guide actually is helping you to install other tools, such as humann3 and kaiju. On the other hand, if you have installed these tools somewhere else, you can just modify the `Snakefile_config.yml` file to skip the above installation procedure.

# 21. fastp
DONE_FILE="install/log/fastp.done"
if [ ! -f "$DONE_FILE" ]; then
	chmod a+x utils/fastp
	touch "$DONE_FILE"
else
	echo "fastp already installed. Skipping this step."
fi

# 22. bwa-mem2
DONE_FILE="install/log/bwa_mem2.done"
if [ ! -f "$DONE_FILE" ]; then
	tar -jxf utils/bwa-mem2-2.2.1_x64-linux.tar.bz2 -C utils/
	chmod a+x utils/bwa-mem2-2.2.1_x64-linux/bwa-mem2
	touch "$DONE_FILE"
else
	echo "bwa-mem2 already installed. Skipping this step."
fi

# 23. gtdbtk
DONE_FILE="install/log/gtdbtk.done"
if [ ! -f "$DONE_FILE" ]; then
	python -m pip install gtdbtk
	gtdbtk_db=$(pwd)/databases/gtdbdk
	mkdir -p $gtdbtk_db
	conda env config vars set GTDBTK_DATA_PATH=$gtdbtk_db
	conda deactivate
	conda activate OUTPOST
	bash utils/gtdk_db_download.sh
	touch "$DONE_FILE"
else
	echo "gtdbtk already installed. Skipping this step."
fi

# 24. salmon
DONE_FILE="install/log/salmon.done"
if [ ! -f "$DONE_FILE" ]; then
	tar -xzf utils/v1.10.0.tar.gz -C utils/
	chmod a+x utils/salmon-latest_linux_x86_64/bin/salmon
	touch "$DONE_FILE"
else
	echo "salmon already installed. Skipping this step."
fi

# 25. eggnog-mapper
DONE_FILE="install/log/eggnog_mapper.done"
if [ ! -f "$DONE_FILE" ]; then
	emapper_db=$(pwd)/databases/emapper_db
	mkdir -p $emapper_db
	python utils/outpost_download_eggnog_data.py -y --data_dir $emapper_db
	touch "$DONE_FILE"
else
	echo "eggnog-mapper already installed. Skipping this step."
fi

# 26. CARD rgi
DONE_FILE="install/log/rgi.done"
if [ ! -f "$DONE_FILE" ]; then
	rgi_db=$(pwd)/databases/rgi
	mkdir -p $rgi_db
	tar -jxf utils/card-data.tar.bz2 -C $rgi_db
	rgi load --local -i $rgi_db/card.json
	rgi card_annotation -i $rgi_db/card.json
	touch "$DONE_FILE"
else
	echo "eggnog-mapper already installed. Skipping this step."
fi

# 27. viruses.fa unzip
DONE_FILE="install/log/viruses.done"
if [ ! -f "$DONE_FILE" ]; then
	unzip utils/viruses.fa.zip -d utils/
	touch "$DONE_FILE"
else
	echo "viruses.fa already unzipped. Skipping this step."
fi


echo "OUTPOST installed successfully!"

