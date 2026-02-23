#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -N descarga_gtex
#$ -o logs/01_descarga_gtex.txt
#$ -e logs/01_descarga_gtex.txt
#$ -m e

echo "**** Job starts ****"
date

# Configurar entorno
module load conda_R/4.4.x

# Ejecutar script
Rscript 01_descarga_gtex.R

echo "**** Job ends ****"
date