## SCRIPT 2. Obtención de las dianas de miRNAs con más relevancia científica de miRTarBase.
# Autor: Manuel Remesal González
# Version: 2/12/2023
# Creación propia.

if (!require(stringr)) {
  install.packages("stringr")
}


library(dplyr)
library(openxlsx)
library(ggplot2)
library(stringr)
library(openxlsx)


# Ruta al archivo de Excel
ruta_archivo <- "hsa_MTI_renamed.xlsx"

# Leer el archivo de Excel
datos <- read.xlsx(ruta_archivo, sheet = 1)

# Crear las nuevas columnas
new_columns <- c("Western blot", "qPCR", "Reporter assay", "Microarray", "NGS", "CLIP-SEQ", "pSILAC", "Other")

# Inicializar las nuevas columnas con valores de cero
datos[, new_columns] <- 0




# Contar la frecuencia de términos por columna para cada fila
for (i in 1:nrow(datos)) {
  experiments <- strsplit(datos$Experiments[i], ";")[[1]]
  term_counts <- table(factor(experiments, levels = new_columns))
  datos[i, names(term_counts)] <- as.numeric(term_counts)
  datos[i, "Total_evidence"] <- sum(term_counts)
}

# Calcular la suma de las columnas "Western blot", "qPCR" y "Reporter assay" en una nueva columna "Strong_evidence"
datos$Strong_evidence <- datos$`Western blot` + datos$qPCR + datos$`Reporter assay`

datos <- datos %>% 
  relocate("Total_evidence", .after = "Strong_evidence")


write.xlsx(datos, "hsa_MTI_cuantificado.xlsx", rowNames = FALSE)

#datos <- read.xlsx("hsa_MTI_cuantificado.xlsx")

# Crear un nuevo dataframe sin los valores NA en la columna "support.type"
datos <- na.omit(datos)

# Filtrar las columnas 10 a 19
columnas_suma <- datos[, 10:19]

# Calcular la suma de las columnas para cada fila con el mismo ID
suma_por_id <- aggregate(columnas_suma, by = list(datos$miRTarBase.ID), FUN = sum)

# Renombrar la columna del ID
colnames(suma_por_id)[1] <- "ID"

# Calcular la suma de las columnas de evidencias para cada fila
suma_por_id$Suma <- rowSums(suma_por_id[, -1], na.rm = TRUE)

# Reordenar las columnas
suma_por_id <- suma_por_id[, c(1, ncol(suma_por_id), 2:(ncol(suma_por_id)-1))]

# Filtrar la tabla original por ID y seleccionar las columnas 1 a 7
tabla_filtrada <- datos[, c("miRTarBase.ID", "miRNA", "Species.(miRNA)", "Target.Gene", "Target.Gene.(Entrez.ID)", "Species.(Target.Gene)", "Experiments")]

# Remover duplicados por ID
tabla_filtrada <- tabla_filtrada[!duplicated(tabla_filtrada$miRTarBase.ID), ]

# Unir las columnas de la tabla suma_por_id a tabla_filtrada
tabla_final <- cbind(tabla_filtrada, suma_por_id)

# Eliminar colunmnas imnecesarias
tabla_final[9:19] <- lapply(tabla_final[9:19], as.numeric)
tabla_final <- tabla_final[, -c(8, 9)]

#Guardamos el progreso como 'tabla_final_filtrado'
write.table(tabla_final, "tabla_final_filtrado.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#tabla_final <- read.table("tabla_final_filtrado.txt", sep = "\t", header = TRUE, quote = "")

# Comprobar valores duplicados
duplicados <- duplicated(tabla_final[,1])
num_duplicados <- sum(duplicados)
duplicados_valores <- tabla_final[,1][duplicados]
print(duplicados_valores)

# Vemos como se distribuyen las evidencias para aplicar un filtrado de las mejores evidencias
summary(tabla_final[,10:17])
boxplot(tabla_final$Strong_evidence, main = "Gráfico de Caja")

# Crear gráfica de barras
frecuencias <- table(tabla_final$Strong_evidence)
barplot(frecuencias, xlab = "Evidencia Umbral", ylab = "Frecuencia",
        main = "Entradas según número de evidencias sin filtrar", col = "darkblue")

# Calcular los percentiles
percentiles <- quantile(tabla_final$Strong_evidence, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999))

# Calcular los permiles
permiles <- quantile(tabla_final$Strong_evidence, probs = seq(0.001, 0.999, by = 0.001))

# Imprimir los resultados
print(permiles)

# Imprimir los resultados
print(percentiles)



# FILTRADO
# Filtrar filas por Strong_evidence >= 5
tabla_filtrada <- tabla_final[tabla_final$Strong_evidence >= 5, ]

#Creamos algunos graficos y calculamos resumen estadistico
summary(tabla_filtrada[,10:17])
boxplot(tabla_filtrada$Strong_evidence, main = "Gráfico de Caja")
hist(tabla_filtrada$Strong_evidence, breaks = 50, col = "lightblue", main = "Histograma de datos tras filtrar")

# Calcular los percentiles
percentiles <- quantile(tabla_filtrada$Strong_evidence, probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999))

# Calcular los permiles
permiles <- quantile(tabla_filtrada$Strong_evidence, probs = seq(0.001, 0.999, by = 0.001))

# Imprimir los resultados
print(permiles)

# Imprimir los resultados
print(percentiles)

# Guardar tabla_final en un archivo CSV
write.table(tabla_filtrada, "tabla_filtrada_mayor_igual5.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#tabla_filtrada <- read.table("tabla_filtrada_mayor_igual5.txt", sep = "\t", header = TRUE, quote = "")

# Obtenemos la lista de genes en formato Entrez ID
tabla_filtrada_genes_unicos <- tabla_filtrada[!duplicated(tabla_filtrada[,5]), ]

# Obtenemos la frecuencia de la columna 'Strong_Evidence'
frecuencias <- table(tabla_filtrada_genes_unicos$Strong_evidence)

# Crear gráfica de barras
barplot(frecuencias, xlab = "Evidencia Umbral", ylab = "Frecuencia",
        main = "Entradas con evidencias mayores o iguales a 5", col = "darkblue")

# Convertir los genes a un único string separado por comas
genes_coma_separados <- paste(tabla_filtrada_genes_unicos, collapse = ",")

# Guardar la cadena de genes en un archivo
writeLines(genes_coma_separados, "lista_genes.txt")

# Creamos la tabla que vamos a usar finalmente: gene_target_list
gene_target_list <- tabla_filtrada_genes_unicos[,c(4, 5)]

# Añadimos una tercera columna que nos servirá más adelante con el paquete DREXmL
gene_target_list$is_target <- TRUE

# Cambiar los nombres de las primeras dos columnas
names(gene_target_list)[1:2] <- c("symbol_id","entrez_id")

# Guardar gene_target_list en un archivo TSV
write.table(gene_target_list, "gene_target_list.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

