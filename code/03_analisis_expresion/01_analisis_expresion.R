#!/usr/bin/env Rscript

# ==================================================
# 01_analisis_expresion.R
# Análisis de expresión diferencial entre 5 tejidos GTEx
# Identificación de genes marcadores tejido-específicos
# Proyecto: Análisis de expresión tejido-específica
# ==================================================

# Cargar librerías
library("limma")
library("edgeR")
library("SummarizedExperiment")
library("sessioninfo")
library("here")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")

# Configurar directorios
dir_processed <- here::here("processed-data")
dir_02 <- file.path(dir_processed, "02_integracion")
dir_03 <- file.path(dir_processed, "03_analisis")
dir_plots <- here::here("plots", "03_analisis")
dir_logs <- here::here("code", "03_analisis_expresion", "logs")

dir.create(dir_03, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)

# Iniciar log
log_file <- file.path(dir_logs, "log_analisis.txt")
sink(log_file, split = TRUE)
cat("========================================\n")
cat("Inicio de análisis:", date(), "\n")
cat("========================================\n\n")

# 1. Cargar datos normalizados
cat("PASO 1: Cargando datos normalizados...\n")

# Cargar objetos del paso anterior
v <- readRDS(file.path(dir_02, "v_voom.rds"))
dge <- readRDS(file.path(dir_02, "dge_normalizado.rds"))
design <- readRDS(file.path(dir_02, "design_matrix.rds"))

# Cargar también el objeto combinado para metadatos
rse_combined <- readRDS(file.path(dir_02, "rse_combined_raw.rds"))

cat("  Datos cargados correctamente\n")
cat("  Dimensiones de v$E:", paste(dim(v$E), collapse = " x "), "\n")
cat("  Tejidos:", paste(colnames(design), collapse = ", "), "\n\n")

# 2. Definir contrastes (comparaciones)
cat("PASO 2: Definiendo contrastes de interés...\n")

# Matriz de contrastes: cada tejido vs el promedio de los demás
# Esto identifica genes específicos de cada tejido
contrastes <- makeContrasts(
  BRAINvsRest = BRAIN - (HEART + LUNG + LIVER + MUSCLE) / 4,
  HEARTvsRest = HEART - (BRAIN + LUNG + LIVER + MUSCLE) / 4,
  LUNGvsRest = LUNG - (BRAIN + HEART + LIVER + MUSCLE) / 4,
  LIVERvsRest = LIVER - (BRAIN + HEART + LUNG + MUSCLE) / 4,
  MUSCLEvsRest = MUSCLE - (BRAIN + HEART + LUNG + LIVER) / 4,
  levels = design
)

cat("  Contrastes definidos:\n")
cat("  - BRAIN vs Resto\n")
cat("  - HEART vs Resto\n")
cat("  - LUNG vs Resto\n")
cat("  - LIVER vs Resto\n")
cat("  - MUSCLE vs Resto\n\n")

# 3. Ajustar modelo lineal
cat("PASO 3: Ajustando modelo lineal con limma...\n")

# Ajustar modelo
fit <- lmFit(v, design)

# Aplicar contrastes
fit_contrasts <- contrasts.fit(fit, contrastes)

# Aplicar eBayes para estadísticas moderadas
fit_eb <- eBayes(fit_contrasts)

cat("  Modelo ajustado correctamente\n")
cat("  Grados de libertad:", fit_eb$df.total[1], "\n\n")

# 4. Extraer resultados para cada tejido
cat("PASO 4: Extrayendo genes diferenciales por tejido...\n")

# Lista para guardar resultados
resultados_tejidos <- list()
top_genes_tejidos <- list()

for (tejido in c("BRAIN", "HEART", "LUNG", "LIVER", "MUSCLE")) {
  contraste <- paste0(tejido, "vsRest")

  # Obtener todos los resultados
  tt <- topTable(fit_eb, coef = contraste, number = Inf, sort.by = "p")

  # Agregar información del gen
  tt$gene_id <- rownames(tt)

  # Clasificar dirección de cambio
  tt$direction <- ifelse(tt$logFC > 0, "up", "down")

  # Guardar
  resultados_tejidos[[tejido]] <- tt

  # Top 50 genes (los más significativos con logFC positivo = específicos del tejido)
  top50 <- head(tt[tt$logFC > 0, ], 50)
  top_genes_tejidos[[tejido]] <- top50

  cat("\n  ---", tejido, "---\n")
  cat("    Genes totales:", nrow(tt), "\n")
  cat(
    "    Genes up (específicos):",
    sum(tt$logFC > 0 & tt$adj.P.Val < 0.05),
    "\n"
  )
  cat("    Genes down:", sum(tt$logFC < 0 & tt$adj.P.Val < 0.05), "\n")
  cat("    Top 3 genes específicos:\n")
  for (i in 1:3) {
    if (i <= nrow(top50)) {
      cat(
        "      ",
        i,
        "-",
        top50$gene_id[i],
        "(logFC =",
        round(top50$logFC[i], 2),
        ", adj.P =",
        format(top50$adj.P.Val[i], scientific = TRUE),
        ")\n"
      )
    }
  }
}

# 5. Guardar resultados completos
cat("\nPASO 5: Guardando resultados...\n")

saveRDS(
  resultados_tejidos,
  file = file.path(dir_03, "resultados_completos.rds")
)
saveRDS(top_genes_tejidos, file = file.path(dir_03, "top_genes_tejidos.rds"))

# Guardar como CSV para fácil consulta
for (tejido in names(resultados_tejidos)) {
  write.csv(
    resultados_tejidos[[tejido]],
    file = file.path(dir_03, paste0("DE_", tejido, "_vs_resto.csv")),
    row.names = FALSE
  )
}

cat("  Resultados guardados en:", dir_03, "\n\n")

# 6. Validación con genes marcadores conocidos
cat("PASO 6: Validando con genes marcadores conocidos...\n")

# Lista de genes clásicos tejido-específicos
genes_conocidos <- list(
  BRAIN = c("SNAP25", "GFAP", "SYN1", "MAP2"),
  HEART = c("MYH6", "TNNT2", "NPPA", "MYL2"),
  LUNG = c("SFTPC", "SFTPA1", "SCGB1A1", "AGER"),
  LIVER = c("ALB", "CYP3A4", "SERPINA1", "APOB"),
  MUSCLE = c("ACTA1", "MYL1", "CKMT2", "DES")
)

# Verificar si estos genes están en nuestros datos
# Necesitamos mapear símbolos a IDs (esto es aproximado)
# En GTEx, los IDs son ENSG, pero podemos buscar por símbolo en rowData
row_data <- as.data.frame(rowData(rse_combined))

# Buscar coincidencias aproximadas (esto puede no encontrar todos)
cat("\nGenes marcadores encontrados:\n")
for (tejido in names(genes_conocidos)) {
  cat("\n", tejido, ":\n")
  for (gen in genes_conocidos[[tejido]]) {
    # Buscar en rowData (esto depende de cómo se llame la columna)
    if ("gene_name" %in% colnames(row_data)) {
      idx <- grep(gen, row_data$gene_name, ignore.case = TRUE)
    } else {
      idx <- grep(gen, rownames(rse_combined), ignore.case = TRUE)
    }

    if (length(idx) > 0) {
      # Ver su logFC en el contraste correspondiente
      res <- resultados_tejidos[[tejido]]
      gene_id <- rownames(rse_combined)[idx[1]]
      if (gene_id %in% rownames(res)) {
        logfc <- res[gene_id, "logFC"]
        adjp <- res[gene_id, "adj.P.Val"]
        cat(
          gen,
          ": logFC =",
          round(logfc, 2),
          ", adj.P =",
          format(adjp, scientific = TRUE),
          "\n"
        )
      } else {
        cat(gen, ": encontrado pero no en resultados\n")
      }
    } else {
      cat(gen, ": no encontrado\n")
    }
  }
}

# 7. Visualizaciones principales
cat("\nPASO 7: Generando visualizaciones...\n")

### 7.1 Heatmap de top genes por tejido
cat("  Generando heatmap de top genes...\n")

# Seleccionar top 20 genes de cada tejido (los más específicos)
genes_seleccionados <- c()
for (tejido in names(top_genes_tejidos)) {
  top20 <- head(top_genes_tejidos[[tejido]], 20)
  genes_seleccionados <- c(genes_seleccionados, rownames(top20))
}
genes_seleccionados <- unique(genes_seleccionados)

# Extraer expresión de estos genes
exprs_heatmap <- v$E[genes_seleccionados, ]

# Crear anotaciones para las muestras
annotation_col <- data.frame(
  Tejido = factor(rse_combined$tejido_nombre)
)
rownames(annotation_col) <- colnames(exprs_heatmap)

# Colores para tejidos
colores <- brewer.pal(5, "Set1")
names(colores) <- levels(annotation_col$Tejido)
annotation_colors <- list(Tejido = colores)

# Heatmap
pdf(file.path(dir_plots, "heatmap_top_genes.pdf"), width = 12, height = 14)
pheatmap(
  exprs_heatmap,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  main = "Top genes específicos por tejido"
)
dev.off()

### 7.2 Boxplots de genes marcadores conocidos
cat("  Generando boxplots de genes marcadores...\n")

# Función para hacer boxplot de un gen
plot_gene_boxplot <- function(gene_id, gene_name) {
  if (!gene_id %in% rownames(v$E)) {
    return(NULL)
  }

  df <- data.frame(
    expresion = v$E[gene_id, ],
    tejido = rse_combined$tejido_nombre
  )

  p <- ggplot(df, aes(x = tejido, y = expresion, fill = tejido)) +
    geom_boxplot() +
    theme_bw() +
    labs(
      title = paste("Expresión de", gene_name),
      x = "Tejido",
      y = "Expresión normalizada (log2)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(p)
}

# Buscar un gen representativo de cada tejido
genes_ejemplo <- list(
  BRAIN = list(
    id = rownames(top_genes_tejidos[["BRAIN"]])[1],
    name = "Top BRAIN"
  ),
  HEART = list(
    id = rownames(top_genes_tejidos[["HEART"]])[1],
    name = "Top HEART"
  ),
  LUNG = list(id = rownames(top_genes_tejidos[["LUNG"]])[1], name = "Top LUNG"),
  LIVER = list(
    id = rownames(top_genes_tejidos[["LIVER"]])[1],
    name = "Top LIVER"
  ),
  MUSCLE = list(
    id = rownames(top_genes_tejidos[["MUSCLE"]])[1],
    name = "Top MUSCLE"
  )
)

for (tejido in names(genes_ejemplo)) {
  p <- plot_gene_boxplot(
    genes_ejemplo[[tejido]]$id,
    genes_ejemplo[[tejido]]$name
  )
  if (!is.null(p)) {
    ggsave(
      file.path(dir_plots, paste0("boxplot_", tejido, "_top.pdf")),
      p,
      width = 8,
      height = 6
    )
  }
}

### 7.3 Volcano plots
cat("  Generando volcano plots...\n")

for (tejido in names(resultados_tejidos)) {
  res <- resultados_tejidos[[tejido]]

  # Crear columna para colorear
  res$color <- "No significativo"
  res$color[res$logFC > 1 & res$adj.P.Val < 0.05] <- "Específico"
  res$color[res$logFC < -1 & res$adj.P.Val < 0.05] <- "Bajo en tejido"

  # Top 5 genes para etiquetar
  top5 <- head(res[res$logFC > 1 & res$adj.P.Val < 0.05, ], 5)

  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(
      values = c(
        "Específico" = "red",
        "Bajo en tejido" = "blue",
        "No significativo" = "gray"
      )
    ) +
    theme_bw() +
    labs(
      title = paste("Volcano plot:", tejido, "vs Resto"),
      x = "log2 Fold Change",
      y = "-log10(adj.P.Val)"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    theme(legend.position = "bottom")

  # Añadir etiquetas de top genes
  if (nrow(top5) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = top5,
        aes(label = rownames(top5)),
        size = 2,
        max.overlaps = 5
      )
  }

  ggsave(
    file.path(dir_plots, paste0("volcano_", tejido, ".pdf")),
    p,
    width = 8,
    height = 7
  )
}

# 8. Resumen de resultados
cat("\n========================================\n")
cat("RESUMEN DE GENES ESPECÍFICOS POR TEJIDO:\n")
cat("========================================\n")
for (tejido in names(resultados_tejidos)) {
  res <- resultados_tejidos[[tejido]]
  n_up <- sum(res$logFC > 1 & res$adj.P.Val < 0.05)
  n_down <- sum(res$logFC < -1 & res$adj.P.Val < 0.05)
  cat(sprintf(
    "%-10s: %4d genes específicos, %4d genes bajos\n",
    tejido,
    n_up,
    n_down
  ))
}

# 9. Guardar resumen en CSV
resumen <- data.frame(
  Tejido = names(resultados_tejidos),
  Genes_especificos = sapply(resultados_tejidos, function(x) {
    sum(x$logFC > 1 & x$adj.P.Val < 0.05)
  }),
  Genes_bajos = sapply(resultados_tejidos, function(x) {
    sum(x$logFC < -1 & x$adj.P.Val < 0.05)
  }),
  Total_genes = sapply(resultados_tejidos, nrow)
)
write.csv(
  resumen,
  file = file.path(dir_03, "resumen_genes_tejido.csv"),
  row.names = FALSE
)

# 10. Información de sesión
cat("\n========================================\n")
cat("Información de sesión:\n")
sessioninfo::session_info()

sink()

cat("\n Análisis completado.\n")
cat(" Resultados guardados en:", dir_03, "\n")
cat(" Gráficos guardados en:", dir_plots, "\n")
cat(" Log guardado en:", log_file, "\n")
