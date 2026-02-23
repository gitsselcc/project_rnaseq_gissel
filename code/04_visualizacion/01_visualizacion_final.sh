#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -N viz_gtex
#$ -o code/04_visualizacion/logs/01_visualizacion_final.txt
#$ -e code/04_visualizacion/logs/01_visualizacion_final.txt
#$ -m e

echo "**** Job starts ****"
date

# Ejecutar script de R
Rscript code/04_visualizacion/01_visualizacion_final.R

echo "**** Job ends ****"
date