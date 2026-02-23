#!/bin/bash
#$ -cwd
#$ -l mem_free=16G,h_vmem=16G,h_fsize=100G
#$ -N integra_gtex
#$ -o code/02_integra_normaliz/logs/01_integra_normaliz.txt
#$ -e code/02_integra_normaliz/logs/01_integra_normaliz.txt
#$ -m e

echo "**** Job starts ****"
date

# Ejecutar script de R
Rscript code/02_integra_normaliz/01_integra_normaliz.R

echo "**** Job ends ****"
date