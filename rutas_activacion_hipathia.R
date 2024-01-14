#setwd("~/OneDrive - Universidad Europea/Doctorado/TFM")
setwd("C:/Users/manur/OneDrive - Universidad Europea/Doctorado/TFM")

# Instalar BiocManager si no está presente
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Función para instalar paquetes de Bioconductor si no están ya instalados
install_bioc_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

# Instalar paquetes de Bioconductor
bioc_packages <- c("org.Hs.eg.db", "biomaRt", "AnnotationDbi", "edgeR", "hipathia")
lapply(bioc_packages, install_bioc_if_missing)

# Instalar paquetes de CRAN si no están presentes
cran_packages <- c("dplyr", "data.table")
lapply(cran_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
})


library(org.Hs.eg.db)
library(data.table)
library(AnnotationDbi)
library(hipathia)
library(dplyr)
library(edgeR)



# Leemos un archivo de texto que contiene los datos de lectura de genes de RNA-Seq y lo asignamos a 'hipathia_data'.
# El argumento 'skip=2' omite las primeras 2 líneas del archivo, que pueden ser líneas de encabezado o metadatos.
# 'header=TRUE' indica que la primera línea que se lee (la tercera en el archivo) contiene los nombres de las columnas.
hipathia_data <- read.table("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip = 2, header = TRUE)

# Creamos un nuevo dataframe 'datos' que es una copia de 'hipathia_data'. Lo hacemos para mantener los datos 
# originales sin cambios.
datos <- hipathia_data

# Limpiamos los identificadores de las isoformas en la columna 'Name' del dataframe 'datos'. La función 'sub' se usa para 
# reemplazar cualquier cosa que siga a un punto en los nombres de los genes con una cadena vacía, para que puedan
# ser traducidos correctamente
datos$Name <- sub("\\..*", "", datos$Name)

# Mapeamos los identificadores Ensembl de los genes a identificadores ENTREZ utilizando la base de datos 'org.Hs.eg.db'.
# 'keys' especifica la columna de 'datos' que contiene los identificadores Ensembl,
# 'column' especifica que queremos obtener los identificadores ENTREZ,
# 'keytype' asegura que los identificadores que estamos proporcionando son del tipo Ensembl,
# y 'multiVals' indica que si hay múltiples valores, tomamos el primero.
datos$Name = mapIds(org.Hs.eg.db,
                    keys=datos$Name,
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

# Filtramos 'datos' para mantener solo las filas donde 'Name' no es NA (es decir, tiene un identificador ENTREZ válido).
datos <- subset(datos, !is.na(Name))

# Elimino la columna de Gene Symbol porque ya no nos hace falta
datos <- datos[,-2]




## ELIMINACIÓN DE DUPLICADOS ##
# Eliminamos duplicados sumando los valores de los
# IDs solo los duplicados para hacerlo más eficiente

# Identificamos los 'Name' duplicados
nombres_duplicados <- unique(datos$Name[duplicated(datos$Name)])

# Filtramos solo las filas con 'Name' duplicados
datos_duplicados <- datos[datos$Name %in% nombres_duplicados, ]

# Transformamos la tabla a DT
setDT(datos_duplicados)

# Sumamos/agrupamos por 'Name' solo para los duplicados
datos_agrupados_duplicados <- datos_duplicados[, lapply(.SD, sum), by = Name]

# Creamos un índice lógico para las filas que no tienen nombres duplicados
indices_no_duplicados <- !(datos$Name %in% nombres_duplicados)

# Filtramos las filas del dataframe original que no tienen nombres duplicados
datos_no_duplicados <- datos[indices_no_duplicados, ]

# Unir los datos no duplicados con los datos agrupados
datos <- rbind(datos_no_duplicados, datos_agrupados_duplicados)

#  Establecemos los nombres de las filas
setDF(datos)  # Convertimos de nuevo a dataframe

# Ahora podemos establecer los nombres de las filas sin problemas
rownames(datos) <- datos$Name

# Eliminamos la primera columna de la matriz que contiene los identificadored de ENTREZ
datos <- datos[,-1]


#NORMALIZACIÓN
# 'Datos' es nuestra matriz de recuentos donde las filas representan genes
# y las columnas representan muestras individuales
dge <- DGEList(counts=datos)


# Calculamos los factores de normalización con el método TMM
dge <- calcNormFactors(dge, method = 'TMM')

# Ahora el objeto 'dge' contiene los factores de normalización
# Podemos revisar los factores de normalización con
dge$samples$norm.factors

# Para obtener los recuentos normalizados, usamos los factores de esta manera
normalized_counts <- cpm(dge, log = TRUE, prior.count = 3)
write.table(normalized_counts, file = "normalized_counts.tsv", sep = "\t", 
             row.names = TRUE, col.names = TRUE, quote = FALSE)
#normalized_counts <- read.table(file = "normalized_counts.tsv", header = TRUE, sep = "\t", row.names = 1)


# Sacamos la tabla de expresión de genes de GTEX, de normalized counts, tabulado y transpuesto
normalized_counts_transpuesto <- as.data.frame(t(normalized_counts))
write.table(normalized_counts_transpuesto, file = "gtex.tsv", sep = "\t", 
            row.names = TRUE, col.names = TRUE, quote = FALSE)






# CALCULO DE LAS RUTAS DE ACTIVACIÓN CON HIPATHIA
# Transformamos nuestra tabla de conteos normalizados en una matriz de datos para 
#procesar con hipathia.
matrix_data <- as.matrix(normalized_counts)

# Utilizamos la función translate_data de hipathia para convertir los datos de la matriz 
# a un formato que pueda ser interpretado por hipathia. 'hsa' indica que los datos se 
# refieren al genoma humano. 'verbose=TRUE' permite imprimir mensajes adicionales durante 
# la ejecución.
trans_data <- translate_data(matrix_data, "hsa", verbose=TRUE)

# Guardamos la tabla de conteos normalizados y los datos traducidos en archivos de texto 
# para su posterior recuperación si es necesario. 'sep = "\t"' indica que el separador 
# es un tabulador, 'row.names = TRUE' y 'col.names = TRUE' indican que se incluyan los 
# nombres de las filas y columnas, y 'quote = FALSE' evita que los textos se guarden 
# entre comillas.
#write.table(trans_data, file = "transdata.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# La siguiente línea está comentada y serviría para leer los datos traducidos desde un archivo si fuese necesario.
# trans_data <- read.table("transdata1.txt", sep = "\t", header = TRUE)

# Cargamos un conjunto específico de rutas de señalización relacionadas con la inflamación en 
# hipathia utilizando sus identificadores KEGG.
pathways_infl <- load_pathways(species = "hsa", 
                               pathways_list = c("hsa04630", "hsa04010", "hsa04064", "hsa04620", 
                                                 "hsa04668", "hsa04621", "hsa04750", "hsa04151"))

# Ejecutamos el análisis de hipathia en los datos traducidos y las rutas de inflamación 
# cargadas sin descomponer las rutas en subrutas, y con mensajes de consola desactivados 
# (verbose=FALSE).
results <- hipathia(trans_data, pathways_infl, decompose=FALSE, verbose=FALSE)

# Imprimimos los resultados obtenidos para verificar el output.
results

# Extraemos los valores de las rutas obtenidos del análisis
# Esto facilita su manipulación y visualización en R.
path_vals<-get_paths_data(results)

# Imprimimos las primeras 4 líneas de la matriz de resultados de las rutas
hhead(path_vals, 4)

# Extraemos los datos
path_vals_matrix <- assay(path_vals)

# Ahora podemos guardar el mapa de la enfermedad, path_vals_matrix en un archivo
path_vals_transpuesto <- as.data.frame(t(path_vals_matrix))
write.table(path_vals_transpuesto, file = "pathvals.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Guardamos los nombres de las vías que vamos usar con DREXML
circuits <- data.frame(
  index = rownames(path_vals),
  in_disease = rep(1, nrow(path_vals)))
write.table(circuits, file = "circuits.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Ahora, obtenemos los nombres comprensibles de las subvías para reemplazar los IDs.
path_names <- get_path_names(pathways_infl, rownames(path_vals))
path_names

# Ahora podemos agregar los nombres legibles a la matriz de valores de las rutas.
rownames(path_vals_matrix) <- path_names

# Ahora transponemos el archivo
path_vals_transpuesto_legible <- as.data.frame(t(path_vals_matrix))
write.table(path_vals_transpuesto_legible, file = "pathvals_nombres_legibles.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)