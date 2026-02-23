#!/bin/bash
#$ -cwd
#$ -l mem_free=16G,h_vmem=16G,h_fsize=100G
#$ -N analisis_gtex
#$ -o code/03_analisis_expresion/logs/01_analisis_expresion.txt
#$ -e code/03_analisis_expresion/logs/01_analisis_expresion.txt
#$ -m e

echo "**** Job starts ****"
date

# Ejecutar script de R
Rscript code/03_analisis_expresion/01_analisis_expresion.R

echo "**** Job ends ****"
date