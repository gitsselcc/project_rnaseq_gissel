#!/usr/bin/env Rscript

# ==================================================
# 01_visualizacion_final.R
# Visualizaciones finales y preparación de figuras para reporte
# Proyecto: Análisis de expresión tejido-específica GTEx
# ==================================================

# Cargar librerías
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("gridExtra")
library("cowplot")
library("sessioninfo")
library("here")
library("reshape2")
library("dplyr")
library("ggrepel")

# Configurar directorios
dir_processed <- here::here("processed-data")
dir_02 <- file.path(dir_processed, "02_integracion")
dir_03 <- file.path(dir_processed, "03_analisis")
dir_plots <- here::here("plots", "04_visualizacion")
dir_logs <- here::here("code", "04_visualizacion", "logs")

dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)

# Iniciar log
log_file <- file.path(dir_logs, "log_visualizacion.txt")
sink(log_file, split = TRUE)
cat("========================================\n")
cat("Inicio de visualización:", date(), "\n")
cat("========================================\n\n")

# 1. Cargar datos necesarios
cat("PASO 1: Cargando datos...\n")

# Cargar resultados
resultados_tejidos <- readRDS(file.path(dir_03, "resultados_completos.rds"))
top_genes_tejidos <- readRDS(file.path(dir_03, "top_genes_tejidos.rds"))

# Cargar datos de expresión
v <- readRDS(file.path(dir_02, "v_voom.rds"))
rse_combined <- readRDS(file.path(dir_02, "rse_combined_raw.rds"))

# Cargar diseño
design <- readRDS(file.path(dir_02, "design_matrix.rds"))

cat("  Datos cargados correctamente\n")
cat("  Tejidos:", paste(names(resultados_tejidos), collapse = ", "), "\n\n")

# 2. Preparar metadatos y colores
cat("PASO 2: Preparando metadatos y paletas de colores...\n")

# Mapa de tejidos a nombres en español
tejido_nombres <- c(
    "BRAIN" = "Cerebro",
    "HEART" = "Corazón",
    "LUNG" = "Pulmón",
    "LIVER" = "Hígado",
    "MUSCLE" = "Músculo"
)

# Crear factor de tejido con nombres en español
tejido_factor <- factor(rse_combined$tejido_nombre, levels = tejido_nombres)

# Paleta de colores para tejidos (5 colores distintos)
colores_tejidos <- c(
    "Cerebro" = "#E41A1C", # rojo
    "Corazón" = "#377EB8", # azul
    "Pulmón" = "#4DAF4A", # verde
    "Hígado" = "#984EA3", # morado
    "Músculo" = "#FF7F00" # naranja
)

cat("  Paleta de colores definida\n\n")


# 4. FIGURA 1: Heatmap de top genes
cat("\nPASO 4: Generando Figura 2 - Heatmap de top genes\n")

# Seleccionar top 20 genes de cada tejido (los más específicos)
genes_seleccionados <- c()
tejido_origen <- c()

for (tejido in names(top_genes_tejidos)) {
    top20 <- head(
        top_genes_tejidos[[tejido]][top_genes_tejidos[[tejido]]$logFC > 0, ],
        20
    )
    if (nrow(top20) > 0) {
        genes_seleccionados <- c(genes_seleccionados, rownames(top20))
        tejido_origen <- c(
            tejido_origen,
            rep(tejido_nombres[tejido], nrow(top20))
        )
    }
}

# Eliminar duplicados (por si un gen es top en múltiples tejidos)
dup_idx <- duplicated(genes_seleccionados)
genes_seleccionados <- genes_seleccionados[!dup_idx]
tejido_origen <- tejido_origen[!dup_idx]

# Extraer expresión
exprs_heatmap <- v$E[genes_seleccionados, ]

# Escalar por filas (genes)
exprs_scaled <- t(scale(t(exprs_heatmap)))

# Anotaciones para columnas (muestras)
annotation_col <- data.frame(
    Tejido = tejido_factor
)
rownames(annotation_col) <- colnames(exprs_scaled)

# Anotaciones para filas (genes) - mostrar solo algunos nombres
annotation_row <- data.frame(
    Tejido_origen = factor(tejido_origen, levels = tejido_nombres)
)
rownames(annotation_row) <- genes_seleccionados

# Mostrar nombres solo para genes seleccionados (para no saturar)
# Elegimos 5 genes representativos por tejido para mostrar nombres
genes_para_mostrar <- c()
for (tejido in names(top_genes_tejidos)) {
    top5 <- head(
        top_genes_tejidos[[tejido]][top_genes_tejidos[[tejido]]$logFC > 0, ],
        5
    )
    if (nrow(top5) > 0) {
        genes_para_mostrar <- c(genes_para_mostrar, rownames(top5))
    }
}

# Crear vector de nombres a mostrar
rownames_labels <- ifelse(
    rownames(exprs_scaled) %in% genes_para_mostrar,
    rownames(exprs_scaled),
    ""
)

# Colores
annotation_colors <- list(
    Tejido = colores_tejidos,
    Tejido_origen = colores_tejidos
)

# Heatmap
pdf(file.path(dir_plots, "Fig2_Heatmap_top_genes.pdf"), width = 14, height = 16)
pheatmap(
    exprs_scaled,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    labels_row = rownames_labels,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    cutree_rows = 5, # Esperamos 5 clusters (uno por tejido)
    scale = "none", # Ya escalamos
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Figura 2: Genes específicos por tejido",
    fontsize_row = 8
)
dev.off()

cat(" Figura 2 guardada\n")

# 5. FIGURA 2: Boxplots de genes marcadores clásicos
cat("\nPASO 5: Generando Figura 3 - Boxplots de genes marcadores\n")

# Función para crear boxplot
create_gene_boxplot <- function(gene_id, gene_name, tejido_esperado) {
    if (!gene_id %in% rownames(v$E)) {
        return(NULL)
    }

    df <- data.frame(
        Expresion = v$E[gene_id, ],
        Tejido = tejido_factor
    )

    # Calcular medias para ordenar
    medias <- tapply(df$Expresion, df$Tejido, mean)
    orden <- names(sort(medias, decreasing = TRUE))

    df$Tejido <- factor(df$Tejido, levels = orden)

    p <- ggplot(df, aes(x = Tejido, y = Expresion, fill = Tejido)) +
        geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
        scale_fill_manual(values = colores_tejidos) +
        theme_bw(base_size = 12) +
        labs(
            title = paste(gene_name),
            subtitle = paste("Gen marcador de", tejido_esperado),
            x = "",
            y = "Expresión normalizada (log2)"
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = element_text(face = "bold")
        )

    return(p)
}

# Genes clásicos para cada tejido (necesitamos IDs reales)
# Como no tenemos los símbolos exactos, usamos los top genes de cada tejido
genes_ejemplo <- list(
    Cerebro = list(
        id = rownames(top_genes_tejidos[["BRAIN"]])[1],
        nombre = "Top gen cerebro"
    ),
    Corazón = list(
        id = rownames(top_genes_tejidos[["HEART"]])[1],
        nombre = "Top gen corazón"
    ),
    Pulmón = list(
        id = rownames(top_genes_tejidos[["LUNG"]])[1],
        nombre = "Top gen pulmón"
    ),
    Hígado = list(
        id = rownames(top_genes_tejidos[["LIVER"]])[1],
        nombre = "Top gen hígado"
    ),
    Músculo = list(
        id = rownames(top_genes_tejidos[["MUSCLE"]])[1],
        nombre = "Top gen músculo"
    )
)

# Crear lista de plots
plots_list <- list()
for (tejido in names(genes_ejemplo)) {
    p <- create_gene_boxplot(
        genes_ejemplo[[tejido]]$id,
        genes_ejemplo[[tejido]]$nombre,
        tejido
    )
    if (!is.null(p)) {
        plots_list[[tejido]] <- p
    }
}

# Combinar en una sola figura
fig3_boxplots <- plot_grid(
    plotlist = plots_list,
    ncol = 2,
    labels = "AUTO",
    label_size = 14
)

ggsave(
    file.path(dir_plots, "Fig3_Boxplots_genes_marcadores.pdf"),
    fig3_boxplots,
    width = 12,
    height = 10
)
ggsave(
    file.path(dir_plots, "Fig3_Boxplots_genes_marcadores.png"),
    fig3_boxplots,
    width = 12,
    height = 10,
    dpi = 300
)

cat("  Figura 3 guardada\n")

# 6. FIGURA 3: Volcano plots combinados
cat("\nPASO 6: Generando Figura 4 - Volcano plots\n")

# Función para crear volcano plot
create_volcano <- function(res, tejido_nombre) {
    res$color <- "No significativo"
    res$color[res$logFC > 1 & res$adj.P.Val < 0.05] <- "Específico"
    res$color[res$logFC < -1 & res$adj.P.Val < 0.05] <- "Bajo en tejido"

    # Top 3 genes para etiquetar
    top3 <- head(res[res$logFC > 1 & res$adj.P.Val < 0.05, ], 3)

    p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_color_manual(
            values = c(
                "Específico" = "red",
                "Bajo en tejido" = "blue",
                "No significativo" = "gray"
            )
        ) +
        theme_bw(base_size = 10) +
        labs(
            title = tejido_nombre,
            x = "log2 Fold Change",
            y = "-log10(adj.P.Val)"
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed",
            alpha = 0.5
        ) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    if (nrow(top3) > 0) {
        p <- p +
            ggrepel::geom_text_repel(
                data = top3,
                aes(label = rownames(top3)),
                size = 2,
                max.overlaps = 3
            )
    }

    return(p)
}

# Crear volcano plots para cada tejido
volcano_plots <- list()
for (tejido in names(resultados_tejidos)) {
    p <- create_volcano(resultados_tejidos[[tejido]], tejido_nombres[tejido])
    volcano_plots[[tejido]] <- p
}

# Combinar en una sola figura
fig4_volcano <- plot_grid(
    plotlist = volcano_plots,
    ncol = 3,
    labels = "AUTO",
    label_size = 14
)

ggsave(
    file.path(dir_plots, "Fig4_Volcano_combinado.pdf"),
    fig4_volcano,
    width = 15,
    height = 10
)
ggsave(
    file.path(dir_plots, "Fig4_Volcano_combinado.png"),
    fig4_volcano,
    width = 15,
    height = 10,
    dpi = 300
)

cat("  Figura 4 guardada\n")

# 7. FIGURA 4: Gráfico de barras - Número de genes específicos
cat("\nPASO 7: Generando Figura 5 - Barplot de genes específicos\n")

# Contar genes específicos por tejido
df_counts <- data.frame(
    Tejido = factor(tejido_nombres, levels = tejido_nombres),
    Especificos = sapply(resultados_tejidos, function(x) {
        sum(x$logFC > 1 & x$adj.P.Val < 0.05)
    }),
    Bajos = sapply(resultados_tejidos, function(x) {
        sum(x$logFC < -1 & x$adj.P.Val < 0.05)
    })
)

# Reformatear para ggplot
df_counts_long <- melt(
    df_counts,
    id.vars = "Tejido",
    variable.name = "Tipo",
    value.name = "Conteo"
)
df_counts_long$Tipo <- factor(
    df_counts_long$Tipo,
    levels = c("Especificos", "Bajos"),
    labels = c("Específicos del tejido", "Reprimidos en el tejido")
)

fig5_barras <- ggplot(
    df_counts_long,
    aes(x = Tejido, y = Conteo, fill = Tipo)
) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(
        values = c(
            "Específicos del tejido" = "red",
            "Reprimidos en el tejido" = "blue"
        )
    ) +
    theme_bw(base_size = 14) +
    labs(
        title = "Figura 5: Genes diferencialmente expresados por tejido",
        subtitle = "Corte: |logFC| > 1 y FDR < 0.05",
        x = "",
        y = "Número de genes",
        fill = "Tipo de cambio"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
    ) +
    geom_text(
        aes(label = Conteo),
        position = position_dodge(width = 0.9),
        vjust = -0.3,
        size = 3.5
    )

ggsave(
    file.path(dir_plots, "Fig5_Barplot_genes_especificos.pdf"),
    fig5_barras,
    width = 10,
    height = 7
)
ggsave(
    file.path(dir_plots, "Fig5_Barplot_genes_especificos.png"),
    fig5_barras,
    width = 10,
    height = 7,
    dpi = 300
)

cat("  Figura 5 guardada\n")

# 8. Generar tabla resumen para el reporte
cat("\nPASO 8: Generando tablas resumen\n")

# Top 10 genes específicos de cada tejido
for (tejido in names(resultados_tejidos)) {
    res <- resultados_tejidos[[tejido]]
    top10 <- res[res$logFC > 1 & res$adj.P.Val < 0.05, ]
    top10 <- head(top10[order(top10$logFC, decreasing = TRUE), ], 10)

    if (nrow(top10) > 0) {
        write.csv(
            top10[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
            file = file.path(dir_plots, paste0("Table_top10_", tejido, ".csv")),
            row.names = TRUE
        )
    }
}

# 9. Resumen final
cat("\n========================================\n")
cat("RESUMEN DE FIGURAS GENERADAS:\n")
cat("========================================\n")
cat("Figura 1: Heatmap - Top genes específicos\n")
cat("Figura 2: Boxplots - Genes marcadores\n")
cat("Figura 3: Volcano plots - Todos los tejidos\n")
cat("Figura 4: Barplot - Conteo de genes específicos\n\n")

cat("Archivos guardados en:\n")
cat("  -", dir_plots, "\n\n")

# 10. Información de sesión
cat("========================================\n")
cat("Información de sesión:\n")
sessioninfo::session_info()

sink()

cat("\n Visualización completada.\n")
cat(" Figuras guardadas en:", dir_plots, "\n")
cat(" Log guardado en:", log_file, "\n")
