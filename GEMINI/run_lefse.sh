#!/bin/bash
## conda activate biobakery_env

namefile=namefile.lefse.txt

:<<"Cm"
## version 1
draw_lefse(){
	input=$1

	lefse-format_input.py ${input} ${input}.format -c 1 -s 2 -u 3 -o 1000000
	run_lefse.py ${input}.format  ${input}.format.res
	lefse-plot_res.py ${input}.format.res  ${input}.format.res.pdf --format pdf 
	lefse-plot_cladogram.py ${input}.format.res  ${input}.format.res.cladogram.pdf --format pdf

}
Cm

## version 2
draw_lefse(){
	input=$1

	lefse_format_input.py ${input} ${input}.format -c 1 -s 2 -u 3 -o 1000000
	lefse_run.py ${input}.format  ${input}.format.res
	lefse_plot_res.py ${input}.format.res  ${input}.format.res.pdf --format pdf 
	lefse_plot_cladogram.py ${input}.format.res  ${input}.format.res.cladogram.pdf --format pdf

}


cat ${namefile} | while read id 
do
	for level in {1..8}
		do
			input=sample16_rel_abun.${level}.rmU.euk.${id}.lefse.tsv
			echo "${input} is processing"
			draw_lefse ${input} > ${input}.log 2>&1 &
		done
done



#draw_lefse horsemuledonkeyhinny_genefamilies_uniref90names_cpm_rxn_unstratifiednamedtsvrela_abun_style.hmdh.lefse.tsv
