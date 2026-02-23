#!/usr/bin/env Rscript

# ==================================================
# 01_descarga_gtex.R
# Descarga de datos GTEX desde recount3
# Tejidos seleccionados: BRAIN, HEART, LUNG, LIVER, MUSCLE
# Proyecto: Análisis de expresión tejido-específica
# ==================================================

# Cargar librerías
library("recount3")
library("sessioninfo")
library("here")

# Configurar directorios
dir_processed <- here::here("processed-data", "01_descarga_gtex")
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)
dir_logs <- here::here("code", "01_descarga_gtex", "logs")
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)

# Iniciar log
log_file <- file.path(dir_logs, "log_descarga.txt")
sink(log_file, split = TRUE)
cat("========================================\n")
cat("Inicio de descarga:", date(), "\n")
cat("========================================\n\n")

# 1. Obtener proyectos disponibles
cat("PASO 1: Obteniendo proyectos humanos de recount3...\n")
human_projects <- available_projects()
cat("Total proyectos humanos:", nrow(human_projects), "\n")

# 2. Obtener proyectos GTEx
cat("\nPASO 2: Obteniendo proyectos GTEx...\n")
gtex_projects <- subset(human_projects, file_source == "gtex")
cat("Total proyectos GTEx:", nrow(gtex_projects), "\n")

# 3. Definir los 5 tejidos de interés (OPCIÓN 1 - Diversidad biológica)
tejidos_interes <- c(
  "BRAIN", # Cerebro - Sistema nervioso
  "HEART", # Corazón - Cardiovascular
  "LUNG", # Pulmón - Respiratorio
  "LIVER", # Hígado - Digestivo/Metabolismo
  "MUSCLE" # Músculo - Musculoesquelético
)

# Verificar que existen en GTEx
tejidos_disponibles <- tejidos_interes[
  tejidos_interes %in% gtex_projects$project
]

cat("\nPASO 3: Tejidos seleccionados:\n")
for (t in tejidos_disponibles) {
  n_samples <- gtex_projects[gtex_projects$project == t, "n_samples"]
  cat("  -", t, ":", n_samples, "muestras\n")
}

if (length(tejidos_disponibles) < length(tejidos_interes)) {
  faltantes <- tejidos_interes[!tejidos_interes %in% gtex_projects$project]
  cat("\n Tejidos no encontrados:", paste(faltantes, collapse = ", "), "\n")
  cat("Usando solo los disponibles.\n")
}
cat("\n")

# 4. Descargar datos para cada tejido
cat("PASO 4: Descargando datos por tejido...\n")
lista_rse <- list()

for (tejido in tejidos_disponibles) {
  cat("\n--- Procesando:", tejido, "---\n")

  # Seleccionar proyecto
  proj_info <- subset(gtex_projects, project == tejido)

  # Descargar datos
  cat("  Descargando datos...\n")
  rse_temp <- tryCatch(
    {
      create_rse(proj_info)
    },
    error = function(e) {
      cat(" ERROR:", e$message, "\n")
      return(NULL)
    }
  )

  if (!is.null(rse_temp)) {
    # Convertir cuentas por nucleótido a cuentas por lectura
    cat("  Convirtiendo cuentas...\n")
    assay(rse_temp, "counts") <- compute_read_counts(rse_temp)

    # Expandir atributos (para GTEx)
    cat("  Expandiendo metadatos...\n")
    rse_temp <- expand_sra_attributes(rse_temp)

    # Guardar en lista
    lista_rse[[tejido]] <- rse_temp

    # Guardar copia individual
    archivo_salida <- file.path(dir_processed, paste0("rse_", tejido, ".rds"))
    saveRDS(rse_temp, file = archivo_salida)

    cat(
      " Completado. Dimensiones:",
      paste(dim(rse_temp), collapse = " x "),
      "\n"
    )
    cat("  Archivo guardado:", archivo_salida, "\n")
  }
}

# 5. Guardar lista completa
cat("PASO 5: Guardando lista completa de tejidos...")
if (length(lista_rse) > 0) {
  archivo_lista <- file.path(dir_processed, "lista_rse_tejidos.rds")
  saveRDS(lista_rse, file = archivo_lista)
  cat(" Archivo guardado en:", archivo_lista, "\n")

  # Mostrar resumen
  cat("RESUMEN DE TEJIDOS DESCARGADOS:")
  for (nombre in names(lista_rse)) {
    dims <- dim(lista_rse[[nombre]])
    cat("  -", nombre, ":", dims[1], "genes,", dims[2], "muestras\n")
  }
} else {
  cat(" No se pudo descargar ningún tejido\n")
}

# 6. Información de sesión
cat("\n========================================\n")
cat("Información de sesión:\n")
sessioninfo::session_info()

sink()

cat("\n Proceso completado. Revisa el log en:", log_file, "\n")
