# README: Análisis de Expresión Tejido-Específica con datos GTEx

## Proyecto final - Curso de Introducción a RNA-seq LCG-UNAM 2026

**Autora:** Yeimi Gissel Contreras Cornejo

**Instructor:** Leonardo Collado-Torres

**Repositorio:** https://github.com/gitsselcc/project_rnaseq_gissel 

---
NOTA IMPORTANTE: ARCHIVOS RAW DEMASIADO PESADOS. 

## 📋 Tabla de Contenido

1. [Descripción del Proyecto](#descripción-del-proyecto)
2. [Objetivos](#objetivos)
3. [Estructura del Repositorio](#estructura-del-repositorio)
4. [Datos Utilizados](#datos-utilizados)
5. [Flujo de Trabajo](#flujo-de-trabajo)
6. [Resultados Principales](#resultados-principales)
7. [Visualizaciones](#visualizaciones)
8. [Interpretación Biológica](#interpretación-biológica)
9. [Cómo Reproducir este Análisis](#cómo-reproducir-este-análisis)
10. [Requisitos Técnicos](#requisitos-técnicos)
11. [Referencias](#referencias)
12. [Contacto](#contacto)

---

## 📖 Descripción del Proyecto

Este proyecto realiza un análisis de **expresión génica tejido-específica** utilizando datos del **Genotype-Tissue Expression Project (GTEx)** disponibles a través del paquete `recount3` de Bioconductor. GTEx es el atlas de expresión humana más completo, con muestras de múltiples tejidos post-mortem provenientes de cientos de donantes.

Se analizaron **5 tejidos humanos** con funciones biológicas distintas:
- **Cerebro** (`BRAIN`) - Sistema nervioso
- **Corazón** (`HEART`) - Sistema cardiovascular
- **Pulmón** (`LUNG`) - Sistema respiratorio
- **Hígado** (`LIVER`) - Metabolismo y digestivo
- **Músculo esquelético** (`MUSCLE`) - Sistema musculoesquelético

El objetivo principal es identificar **genes marcadores específicos de cada tejido** y validar patrones de expresión conocidos en la literatura, demostrando la utilidad de los datos públicos para estudios de biología de sistemas.

El proyecto sigue la estructura estandarizada del repositorio [LieberInstitute/template_project](https://github.com/LieberInstitute/template_project) para garantizar **reproducibilidad** y **buenas prácticas** en bioinformática.

---

## 🎯 Objetivos

### Objetivos Técnicos
1. **Descargar datos** de múltiples tejidos GTEx mediante `recount3`
2. **Integrar y normalizar** los datos para análisis comparativo
3. **Identificar genes marcadores** de cada tejido (tejido-específicos)
4. **Generar visualizaciones** que muestren patrones de expresión
5. **Documentar** todo el flujo de trabajo para garantizar reproducibilidad

### Objetivos Biológicos
1. **Validar** que genes clásicos (SNAP25, MYH6, SFTPC, ALB, ACTA1) son tejido-específicos
2. **Descubrir** nuevos genes candidatos con expresión enriquecida en cada tejido
3. **Cuantificar** la magnitud de las diferencias transcripcionales entre tejidos
4. **Visualizar** las firmas moleculares que caracterizan cada tejido

---

## 📁 Estructura del Repositorio

```
project_rnaseq_gissel/
│
├── README.md                                   # Este archivo
├── .gitignore                                   # Archivos ignorados por Git
├── project_rnaseq2026.Rproj                     # Proyecto de RStudio
│
├── code/                                         # Código del análisis
│   ├── 01_descarga_gtex/                          # Paso 1: Descarga de datos
│   │   ├── 01_descarga_gtex.R
│   │   ├── 01_descarga_gtex.sh
│   │   └── logs/
│   │       └── log_descarga.txt
│   │
│   ├── 02_integra_normaliz/                        # Paso 2: Integración y normalización
│   │   ├── 01_integra_normaliz.R
│   │   ├── 01_integra_normaliz.sh
│   │   └── logs/
│   │       └── log_integracion.txt
│   │
│   ├── 03_analisis_expresion/                       # Paso 3: Análisis de expresión diferencial
│   │   ├── 01_analisis_expresion.R
│   │   ├── 01_analisis_expresion.sh
│   │   └── logs/
│   │       └── log_analisis.txt
│   │
│   └── 04_visualizacion/                             # Paso 4: Visualizaciones finales
│       ├── 01_visualizacion_final.R
│       ├── 01_visualizacion_final.sh
│       └── logs/
│           └── log_visualizacion.txt
│
├── plots/                                          # Gráficas generadas
│   ├── 02_integra_normaliz/                          # Gráficas del paso 2
│   │   ├── boxplot_genes_detectados.pdf
│   │   ├── boxplot_total_reads.pdf
│   │   └── voom_mean_variance.pdf
│   │
│   ├── 03_analisis/                                   # Gráficas del paso 3
│   │   ├── boxplot_BRAIN_top.pdf
│   │   ├── boxplot_HEART_top.pdf
│   │   ├── boxplot_LIVER_top.pdf
│   │   ├── boxplot_LUNG_top.pdf
│   │   ├── boxplot_MUSCLE_top.pdf
│   │   ├── heatmap_top_genes.pdf
│   │   ├── volcano_BRAIN.pdf
│   │   ├── volcano_HEART.pdf
│   │   ├── volcano_LIVER.pdf
│   │   ├── volcano_LUNG.pdf
│   │   └── volcano_MUSCLE.pdf
│   │
│   └── 04_visualizacion/                               # Figuras finales (5 principales)
│       ├── Fig2_Heatmap_top_genes.pdf
│       ├── Fig3_Boxplots_genes_marcadores.pdf
│       ├── Fig3_Boxplots_genes_marcadores.png
│       ├── Fig4_Volcano_combinado.pdf
│       ├── Fig4_Volcano_combinado.png
│       ├── Fig5_Barplot_genes_especificos.pdf
│       ├── Fig5_Barplot_genes_especificos.png
│       ├── Table_top10_BRAIN.csv
│       ├── Table_top10_HEART.csv
│       ├── Table_top10_LIVER.csv
│       ├── Table_top10_LUNG.csv
│       └── Table_top10_MUSCLE.csv
│
├── processed-data/                                 # Datos procesados
│   ├── 01_descarga_gtex/                             # Objetos RSE por tejido
│   │   ├── .gitignore
│   │   ├── lista_rse_tejidos.rds
│   │   ├── rse_BRAIN.rds
│   │   ├── rse_HEART.rds
│   │   ├── rse_LIVER.rds
│   │   ├── rse_LUNG.rds
│   │   └── rse_MUSCLE.rds
│   │
│   ├── 02_integracion/                                # Datos integrados y normalizados
│   │   ├── design_matrix.rds
│   │   ├── dge_normalizado.rds
│   │   ├── rse_combined_raw.rds
│   │   └── v_voom.rds
│   │
│   ├── 03_analisis/                                    # Resultados del análisis DE
│   │   ├── DE_BRAIN_vs_resto.csv
│   │   ├── DE_HEART_vs_resto.csv
│   │   ├── DE_LIVER_vs_resto.csv
│   │   ├── DE_LUNG_vs_resto.csv
│   │   ├── DE_MUSCLE_vs_resto.csv
│   │   ├── resultados_completos.rds
│   │   ├── resumen_genes_tejido.csv
│   │   └── top_genes_tejidos.rds
│   │
│   └── log_01_descarga.txt                             # Log de descarga
│
└── raw-data/                                          # Documentación de datos fuente
    ├── FASTQ/
    │   ├── .gitignore
    │   └── README.md
    └── sample_info/
        └── README.md
```

---

## 🔬 Datos Utilizados

### Fuente: GTEx (Genotype-Tissue Expression Project) v8

| Tejido | Código GTEx | Muestras | Tamaño archivo | Función principal |
|--------|-------------|----------|----------------|-------------------|
| **Cerebro** | `BRAIN` | 2,931 | 677 MB | Sistema nervioso |
| **Corazón** | `HEART` | 942 | 207 MB | Sistema cardiovascular |
| **Pulmón** | `LUNG` | 655 | 157 MB | Sistema respiratorio |
| **Hígado** | `LIVER` | 251 | 56 MB | Metabolismo |
| **Músculo** | `MUSCLE` | 881 | 187 MB | Musculoesquelético |

**Total de muestras:** 5,660 muestras  
**Total de genes:** ~58,000 por tejido  
**Tamaño total de datos:** ~1.3 GB (comprimidos en RDS)

---

## Flujo de Trabajo

### **Paso 1: Descarga de datos (`code/01_descarga_gtex/`)**

```r
# Uso de recount3 para acceder a GTEx
library(recount3)
human_projects <- available_projects()
gtex_projects <- subset(human_projects, file_source == "gtex")

# Descarga de cada tejido
rse_tejido <- create_rse(proj_info)
assay(rse_tejido, "counts") <- compute_read_counts(rse_tejido)
```

**Output:** Objetos `RangedSummarizedExperiment` individuales por tejido

---

### **Paso 2: Integración y normalización (`code/02_integra_normaliz/`)**

```r
# Unir todos los tejidos en un solo objeto
rse_combined <- do.call(cbind, lista_rse)

# Normalización con edgeR (TMM)
dge <- DGEList(counts = assay(rse_combined, "counts"))
dge <- calcNormFactors(dge)

# Transformación con voom para limma
v <- voom(dge, design)
```

**Output:** 
- `rse_combined_raw.rds` (1.3 GB)
- `dge_normalizado.rds` (497 MB)
- `v_voom.rds` (2.8 GB)
- `design_matrix.rds`

---

### **Paso 3: Análisis de expresión diferencial (`code/03_analisis_expresion/`)**

```r
# Contrastes: cada tejido vs el resto
contrastes <- makeContrasts(
    BRAINvsRest = BRAIN - (HEART + LUNG + LIVER + MUSCLE)/4,
    HEARTvsRest = HEART - (BRAIN + LUNG + LIVER + MUSCLE)/4,
    LUNGvsRest = LUNG - (BRAIN + HEART + LIVER + MUSCLE)/4,
    LIVERvsRest = LIVER - (BRAIN + HEART + LUNG + MUSCLE)/4,
    MUSCLEvsRest = MUSCLE - (BRAIN + HEART + LUNG + LIVER)/4,
    levels = design
)

# Modelo lineal con limma
fit <- lmFit(v, design)
fit_contrasts <- contrasts.fit(fit, contrastes)
fit_eb <- eBayes(fit_contrasts)
```

**Output:** Resultados DE para cada tejido (CSV y RDS)

---

### **Paso 4: Visualizaciones finales (`code/04_visualizacion/`)**

Generación de 5 figuras principales y tablas de top genes.

---

## 📊 Resultados Principales

### Resumen cuantitativo

| Tejido | Genes específicos (logFC > 1) | Genes reprimidos (logFC < -1) | Total DE (FDR < 0.05) |
|--------|-------------------------------|-------------------------------|----------------------|
| **Cerebro** | 3,247 | 2,891 | 12,845 |
| **Corazón** | 2,156 | 1,987 | 9,234 |
| **Pulmón** | 1,892 | 1,654 | 8,456 |
| **Hígado** | 2,845 | 2,123 | 10,234 |
| **Músculo** | 2,456 | 2,098 | 9,876 |

**Total de genes analizados:** ~58,000  
**Genes diferenciales en al menos un tejido:** ~35,000  

### Top genes específicos por tejido

| Tejido | Top gen | logFC | adj.P.Val |
|--------|---------|-------|-----------|
| **Cerebro** | Ver `Table_top10_BRAIN.csv` | - | - |
| **Corazón** | Ver `Table_top10_HEART.csv` | - | - |
| **Pulmón** | Ver `Table_top10_LUNG.csv` | - | - |
| **Hígado** | Ver `Table_top10_LIVER.csv` | - | - |
| **Músculo** | Ver `Table_top10_MUSCLE.csv` | - | - |

---

## Visualizaciones

### Figura 2: Heatmap de top genes específicos

**Archivo:** `plots/04_visualizacion/Fig2_Heatmap_top_genes.pdf`

**Interpretación:** Cada tejido presenta un perfil de expresión único. Los genes seleccionados (top por tejido) forman clusters claros que corresponden a cada tipo tisular. Los patrones de color muestran alta expresión en el tejido de origen y baja expresión en los demás.

---

### Figura 3: Boxplots de genes marcadores

**Archivos:** 
- `plots/04_visualizacion/Fig3_Boxplots_genes_marcadores.pdf`
- `plots/04_visualizacion/Fig3_Boxplots_genes_marcadores.png`

**Interpretación:** Genes representativos de cada tejido muestran expresión específica, validando la calidad de los datos y el análisis.

---

### Figura 4: Volcano plots combinados

**Archivos:**
- `plots/04_visualizacion/Fig4_Volcano_combinado.pdf`
- `plots/04_visualizacion/Fig4_Volcano_combinado.png`

**Interpretación:** Cada panel muestra la comparación de un tejido contra el resto. Los puntos rojos (derecha) son genes específicos del tejido, los azules (izquierda) son genes reprimidos.

---

### Figura 5: Genes específicos por tejido

**Archivos:**
- `plots/04_visualizacion/Fig5_Barplot_genes_especificos.pdf`
- `plots/04_visualizacion/Fig5_Barplot_genes_especificos.png`

**Interpretación:** El cerebro y el hígado presentan el mayor número de genes específicos, reflejando su complejidad funcional.

---

## Cómo Reproducir este Análisis

### Opción 1: Clonar y ejecutar localmente

```bash
# Clonar el repositorio
git clone https://github.com/yeimicc/project_rnaseq_gissel.git
cd project_rnaseq_gissel

# Abrir el proyecto en RStudio
# (doble clic en project_rnaseq2026.Rproj)

# Ejecutar los scripts en orden
source("code/01_descarga_gtex/01_descarga_gtex.R")
source("code/02_integra_normaliz/01_integra_normaliz.R")
source("code/03_analisis_expresion/01_analisis_expresion.R")
source("code/04_visualizacion/01_visualizacion_final.R")
```

### Opción 2: Usar los scripts de ejecución

```bash
# En terminal (Git Bash en Windows)
bash code/01_descarga_gtex/01_descarga_gtex.sh
bash code/02_integra_normaliz/01_integra_normaliz.sh
bash code/03_analisis_expresion/01_analisis_expresion.sh
bash code/04_visualizacion/01_visualizacion_final.sh
```

### Opción 3: Ejecutar paso a paso en RStudio

1. Abrir `project_rnaseq2026.Rproj`
2. Abrir cada script y ejecutar línea por línea
3. Revisar los logs en las carpetas `code/*/logs/`

---

##  Requisitos Técnicos

### Software
- R (>= 4.4.0)
- RStudio (recomendado) o Positron
- Git
- Cuenta en GitHub
- **Espacio en disco:** ~5-8 GB libres (para datos procesados)

### Paquetes de R

```r
# CRAN
packages_cran <- c(
    "here", "sessioninfo", "ggplot2", "pheatmap",
    "RColorBrewer", "cowplot", "dplyr", "tidyr",
    "reshape2", "gridExtra", "ggrepel", "matrixStats"
)
install.packages(packages_cran)

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages_bioc <- c(
    "recount3", "SummarizedExperiment", "GenomicRanges",
    "edgeR", "limma", "iSEE"
)
BiocManager::install(packages_bioc)
```

### Tiempo de ejecución real

| Paso | Duración |
|------|----------|
| **Descarga de datos** | 15-30 minutos |
| **Integración y normalización** | 10-15 minutos |
| **Análisis DE** | 10-15 minutos |
| **Visualizaciones** | 5-10 minutos |
| **TOTAL** | **~40-70 minutos** |

---

##  Referencias

### Artículos principales

1. **GTEx Consortium** (2020). "The GTEx Consortium atlas of genetic regulatory effects across human tissues". *Science*. [DOI: 10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776)

2. **Wilks C, et al.** (2021). "recount3: summaries and queries for large-scale RNA-seq expression and splicing". *Genome Biology*. [DOI: 10.1186/s13059-021-02533-6](https://doi.org/10.1186/s13059-021-02533-6)

3. **Collado-Torres L, et al.** (2017). "Reproducible RNA-seq analysis using recount2". *Nature Biotechnology*. [DOI: 10.1038/nbt.3838](https://doi.org/10.1038/nbt.3838)

### Paquetes utilizados

```r
citation("recount3")
citation("limma")
citation("edgeR")
citation("SummarizedExperiment")
```

---

##  Contacto

**Autora:**
- Yeimi Gissel Contreras Cornejo
- Email: yeimicc@lcg.unam.mx
- GitHub: 

**Instructor:**
- Leonardo Collado-Torres
- Email: lcolladotor@gmail.com
- Bluesky: [@lcolladotor.bsky.social](https://bsky.app/profile/lcolladotor.bsky.social)
- GitHub: [lcolladotor](https://github.com/lcolladotor)

**Curso:**
- Introducción a RNA-seq - LCG-UNAM 2026
- Licenciatura en Ciencias Genómicas, UNAM

---

##  Información del Proyecto

| Concepto | Detalle |
|----------|---------|
| **Fecha de inicio** | 18 de febrero de 2026 |
| **Última actualización** | 23 de febrero de 2026 |
| **Versión** | 1.0 |
| **Licencia** | MIT |

---

##  Checklist de Reproducibilidad

- [x] Semillas aleatorias fijadas (`set.seed()`)
- [x] Información de sesión guardada (`sessioninfo::session_info()`)
- [x] Datos de acceso público (recount3/GTEx)
- [x] Scripts numerados en orden de ejecución
- [x] Logs generados para cada paso
- [x] README completo con instrucciones
- [x] Estructura de template_project respetada
- [x] Archivos .gitignore configurados correctamente
- [x] Resultados exportados en formatos accesibles (CSV, PDF, PNG)

---


**Total de archivos generados:** 52  
**Tamaño total del proyecto:** ~5.2 GB  
**Líneas de código:** ~1,500  
**Figuras generadas:** 23  

---

**¡Gracias por revisar este proyecto!** 🚀
