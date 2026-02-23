#!/usr/bin/env Rscript

# ==================================================
# 01_integra_normaliz.R
# Integración y normalización de 5 tejidos GTEx
# Proyecto: Análisis de expresión tejido-específica
# ==================================================

# Cargar librerías
library("recount3")
library("SummarizedExperiment")
library("edgeR")
library("limma")
library("sessioninfo")
library("here")
library("ggplot2")

# Configurar directorios
dir_processed <- here::here("processed-data")
dir_01 <- file.path(dir_processed, "01_descarga_gtex")
dir_02 <- file.path(dir_processed, "02_integracion")
dir_plots <- here::here("plots", "02_integra_normaliz")
dir_logs <- here::here("code", "02_integra_normaliz", "logs")

dir.create(dir_02, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)

# Iniciar log
log_file <- file.path(dir_logs, "log_integracion.txt")
sink(log_file, split = TRUE)
cat("========================================\n")
cat("Inicio de integración:", date(), "\n")
cat("========================================\n\n")

# 1. Cargar los datos descargados
cat("PASO 1: Cargando datos de los 5 tejidos...\n")

# Lista de tejidos que descargamos
tejidos <- c("BRAIN", "HEART", "LUNG", "LIVER", "MUSCLE")

# Cargar cada tejido
lista_rse <- list()
for (tejido in tejidos) {
    archivo <- file.path(dir_01, paste0("rse_", tejido, ".rds"))
    if (file.exists(archivo)) {
        cat("  Cargando", tejido, "...\n")
        lista_rse[[tejido]] <- readRDS(archivo)
    } else {
        cat(" No se encontró", archivo, "\n")
    }
}

# Verificar que todos se cargaron
cat("\nTejidos cargados:", paste(names(lista_rse), collapse = ", "), "\n")
cat("Número de tejidos:", length(lista_rse), "\n\n")

# 2. Unir todos los tejidos en un solo objeto
cat("PASO 2: Uniendo todos los tejidos...\n")

# Primero, necesitamos que todos tengan los mismos genes
# Tomamos los nombres de genes del primer tejido
genes_comunes <- rownames(lista_rse[[1]])

# Verificamos que los demás tengan los mismos genes
for (tejido in names(lista_rse)[-1]) {
    genes_actuales <- rownames(lista_rse[[tejido]])
    if (!all(genes_comunes == genes_actuales)) {
        cat(tejido, "tiene diferentes genes. Usando intersección...\n")
        genes_comunes <- intersect(genes_comunes, genes_actuales)
    }
}

cat("  Genes comunes a todos los tejidos:", length(genes_comunes), "\n")

# Subsetear cada tejido a los genes comunes
for (tejido in names(lista_rse)) {
    lista_rse[[tejido]] <- lista_rse[[tejido]][genes_comunes, ]
}

# Ahora unimos por columnas (cbind)
# Pero antes, agregamos una columna con el nombre del tejido en colData
for (tejido in names(lista_rse)) {
    lista_rse[[tejido]]$tejido <- tejido
    lista_rse[[tejido]]$tejido_nombre <- switch(
        tejido,
        "BRAIN" = "Cerebro",
        "HEART" = "Corazón",
        "LUNG" = "Pulmón",
        "LIVER" = "Hígado",
        "MUSCLE" = "Músculo",
        tejido
    )
}

# Unir todos los objetos
rse_combined <- do.call(cbind, lista_rse)

cat("  Objeto combinado creado:\n")
cat("  - Genes:", nrow(rse_combined), "\n")
cat("  - Muestras totales:", ncol(rse_combined), "\n")
cat(
    "  - Tejidos incluidos:",
    paste(unique(rse_combined$tejido), collapse = ", "),
    "\n\n"
)

# 3. Control de calidad básico
cat("PASO 3: Control de calidad...\n")

# Verificar valores missing
if (any(is.na(assay(rse_combined, "counts")))) {
    cat(" Hay valores NA en las cuentas. Eliminando...\n")
    assay(rse_combined, "counts")[is.na(assay(rse_combined, "counts"))] <- 0
}

# Calcular estadísticas por muestra
rse_combined$total_reads <- colSums(assay(rse_combined, "counts"))
rse_combined$detected_genes <- colSums(assay(rse_combined, "counts") > 0)

cat("\nResumen de lecturas totales por tejido:\n")
for (tejido in unique(rse_combined$tejido)) {
    idx <- which(rse_combined$tejido == tejido)
    cat(
        "  -",
        tejido,
        ": media =",
        round(mean(rse_combined$total_reads[idx]) / 1e6, 1),
        "M lecturas, rango = [",
        round(min(rse_combined$total_reads[idx]) / 1e6, 1),
        "-",
        round(max(rse_combined$total_reads[idx]) / 1e6, 1),
        "M]\n"
    )
}

# 4. Normalización con edgeR
cat("\nPASO 4: Normalización con edgeR...\n")

# Crear objeto DGEList
dge <- DGEList(
    counts = assay(rse_combined, "counts"),
    genes = rowData(rse_combined),
    samples = colData(rse_combined)
)

# Calcular factores de normalización (TMM)
dge <- calcNormFactors(dge)

cat("  Factores de normalización calculados\n")
cat("  Rango de factores:", round(range(dge$samples$norm.factors), 3), "\n")

# 5. Transformación con voom para limma
# 5. Transformación con voom para limma
cat("\nPASO 5: Transformación voom...\n")

# Diseño experimental (modelo)
tejido_factor <- factor(dge$samples$tejido)
design <- model.matrix(~ 0 + tejido_factor)
colnames(design) <- levels(tejido_factor)

# Voom - SIN plot=TRUE para evitar conflictos
v <- voom(dge, design, plot = FALSE)

# Ahora SÍ generamos nuestro propio plot
pdf(file.path(dir_plots, "voom_mean_variance.pdf"), width = 8, height = 6)

# Calcular medias y varianzas
mean_exp <- rowMeans(v$E)
var_exp <- rowVars(v$E)

# Plot
plot(
    mean_exp,
    var_exp,
    pch = 16,
    cex = 0.5,
    col = rgb(0, 0, 0, 0.3),
    xlab = "log2(mean count)",
    ylab = "log2(variance)",
    main = "Voom: Mean-variance trend"
)

# Añadir línea de tendencia (lowess)
lines(lowess(mean_exp, var_exp), col = "red", lwd = 2)

# Añadir leyenda
legend("topright", legend = "Tendencia", col = "red", lwd = 2, bty = "n")

dev.off()

cat("  Transformación voom completada\n")
cat(
    "  Gráfico guardado en:",
    file.path(dir_plots, "voom_mean_variance.pdf"),
    "\n"
)
# 6. Guardar objetos procesados
cat("\nPASO 6: Guardando objetos procesados...\n")

# Guardar objeto combinado
saveRDS(rse_combined, file = file.path(dir_02, "rse_combined_raw.rds"))

# Guardar DGEList
saveRDS(dge, file = file.path(dir_02, "dge_normalizado.rds"))

# Guardar datos voom
saveRDS(v, file = file.path(dir_02, "v_voom.rds"))

# Guardar matriz de diseño
saveRDS(design, file = file.path(dir_02, "design_matrix.rds"))

# 7. Visualización exploratoria rápida
cat("\nPASO 7: Generando visualizaciones exploratorias...\n")

# Boxplot de lecturas totales por tejido
df_reads <- data.frame(
    tejido = rse_combined$tejido_nombre,
    total_reads = rse_combined$total_reads / 1e6,
    detected_genes = rse_combined$detected_genes
)

p1 <- ggplot(df_reads, aes(x = tejido, y = total_reads, fill = tejido)) +
    geom_boxplot() +
    theme_bw() +
    labs(
        title = "Lecturas totales por tejido",
        x = "Tejido",
        y = "Total reads (millones)"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave(
    file.path(dir_plots, "boxplot_total_reads.pdf"),
    p1,
    width = 8,
    height = 6
)

p2 <- ggplot(df_reads, aes(x = tejido, y = detected_genes, fill = tejido)) +
    geom_boxplot() +
    theme_bw() +
    labs(
        title = "Genes detectados por tejido",
        x = "Tejido",
        y = "Número de genes (counts > 0)"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave(
    file.path(dir_plots, "boxplot_genes_detectados.pdf"),
    p2,
    width = 8,
    height = 6
)

cat("  Gráficos guardados en:", dir_plots, "\n")

# 8. Resumen final
cat("\n========================================\n")
cat("RESUMEN FINAL:\n")
cat("========================================\n")
cat(
    "Tejidos integrados:",
    paste(unique(rse_combined$tejido), collapse = ", "),
    "\n"
)
cat("Total genes:", nrow(rse_combined), "\n")
cat("Total muestras:", ncol(rse_combined), "\n")
cat("Distribución de muestras por tejido:\n")
print(table(rse_combined$tejido_nombre))

# 9. Información de sesión
cat("\n========================================\n")
cat("Información de sesión:\n")
sessioninfo::session_info()

sink()

cat("\n Proceso completado. Revisa el log en:", log_file, "\n")
cat(" Objetos guardados en:", dir_02, "\n")
cat(" Gráficos guardados en:", dir_plots, "\n")
