## SCRIPT 5. Análisis de los resultados de DRExM³L y para la obtención de gráficos relacionados con estos datos
# Autor: Manuel Remesal González
# Version: 20/12/2023
# Creación propia.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("enrichplot")

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(dplyr)
library(tidyr)
library(tidyverse)
library(scales)



# Cargamos las tablas de DRExM3L
#shap_selection <- read.table(file = "shap_selection.tsv", header = TRUE, sep = "\t", row.names = 1)
#shap_summary <- read.table(file = "shap_summary.tsv", header = TRUE, sep = "\t", row.names = 1)
#stability_results <- read.table(file = "stability_results.tsv", header = TRUE, sep = "\t", row.names = 1)

shap_selection_symbol <- read.table(file = "shap_selection_symbol.tsv", header = TRUE, sep = "\t", row.names = 1)
shap_summary_symbol <- read.table(file = "shap_summary_symbol.tsv", header = TRUE, sep = "\t", row.names = 1)
stability_results_symbol <- read.table(file = "stability_results_symbol.tsv", header = TRUE, sep = "\t", row.names = 1)

# Creamos un vector con los nombres de las filas donde 'stability' >= 0.75
selected_rows <- rownames(stability_results_symbol[stability_results_symbol$stability >= 0.75, ])

# Filtrar las filas en todos los archivos usando los nombres seleccionados y eliminar filas con NA
filtered_shap_selection_symbol <- na.omit(shap_selection_symbol[selected_rows, , drop = FALSE])
filtered_shap_summary_symbol <- na.omit(shap_summary_symbol[selected_rows, , drop = FALSE])
filtered_stability_results_symbol <- na.omit(stability_results_symbol[selected_rows, , drop = FALSE])


# Guardar los data frames filtrados en archivos TSV
write.table(filtered_shap_selection_symbol, file = "filtered_shap_selection_symbol.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(filtered_shap_summary_symbol, file = "filtered_shap_summary_symbol.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(filtered_stability_results_symbol, file = "filtered_stability_results_symbol.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Nos aseguramos de que los valores 'True' y 'False' se conviertan a valores lógicos
filtered_shap_selection_symbol_True <- sapply(filtered_shap_selection_symbol, function(x) x == "True")

# Obtener los nombres de las columnas que tienen al menos un valor TRUE
genes_with_true <- colnames(filtered_shap_selection_symbol_True)[colSums(filtered_shap_selection_symbol_True, na.rm = TRUE) > 0]

# Ver los resultados
genes_with_true

# Crear un data frame con los genes
genes_data_frame <- data.frame(miRTar = genes_with_true)

# Guardar el data frame en un archivo TSV
write.table(genes_data_frame, file = "miRTar_list.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




# CONTEO DE RUTAS POR DIANA
# Limpiar los nombres de las rutas en los nombres de las filas
filtered_shap_selection_symbol <- filtered_shap_selection_symbol[, genes_with_true]
lista_nombres_rutas <- row.names(filtered_shap_selection_symbol)  # Extraer los nombres actuales de las filas
nombres_rutas <- sub(":.*$", "", lista_nombres_rutas)  # Eliminar todo después de los dos puntos
filtered_shap_selection_symbol <- data.frame(nombres_rutas, filtered_shap_selection_symbol)

# Asegurarse de que 'df' es tu dataframe y que 'nuevos_nombres' es la columna con los nombres de las rutas
# Convertir todos los valores 'True'/'False' de tipo caracter a valores lógicos
filtered_shap_selection_symbol[, -1] <- lapply(filtered_shap_selection_symbol[, -1], function(x) ifelse(x == 'True', TRUE, FALSE))

# Realizar el conteo
conteo <- filtered_shap_selection_symbol %>% 
  group_by(nombres_rutas) %>% 
  summarise(across(everything(), sum, na.rm = TRUE))


write.table(conteo, file = "conteo.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# Sumamos las interacciones para cada ruta
conteo_aggregated <- conteo %>%
  rowwise() %>%
  mutate(TotalInteractions = sum(c_across(-nombres_rutas))) %>%
  ungroup() %>%
  select(nombres_rutas, TotalInteractions)

# Ordenamos los datos por interar
conteo_aggregated <- conteo_aggregated %>%
  arrange(desc(TotalInteractions)) %>%
  mutate(nombres_rutas = factor(nombres_rutas, levels = nombres_rutas))

ggplot(conteo_aggregated, aes(x = reorder(nombres_rutas, -TotalInteractions), y = TotalInteractions)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  labs(y = "Número total de interacciones", x = "Rutas", title = "Interacciones de Dianas por ruta") +
  theme(legend.position = "none") 


# Sumar las interacciones para cada ruta
conteo_aggregated <- conteo %>%
  rowwise() %>%
  mutate(TotalInteractions = sum(c_across(-nombres_rutas))) %>%
  ungroup()

# Calcular el total de interacciones en todas las rutas
total_interacciones <- sum(conteo_aggregated$TotalInteractions)

# Calcular el porcentaje de interacciones para cada ruta
conteo_aggregated <- conteo_aggregated %>%
  mutate(PercentageOfInteractions = (TotalInteractions / total_interacciones) * 100) %>%
  arrange(desc(PercentageOfInteractions)) %>%
  mutate(nombres_rutas = factor(nombres_rutas, levels = nombres_rutas))

# Crear un gráfico de barras con porcentajes
ggplot(conteo_aggregated, aes(x = nombres_rutas, y = PercentageOfInteractions)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  labs(y = "Porcentaje del Total de Interacciones", x = "Rutas", title = "Porcentaje del Total de Dianas por Ruta") +
  theme(legend.position = "none")




# Convert the data from wide to long format, excluding 'TotalInteractions'
conteo_long <- conteo %>%
  pivot_longer(cols = -c(nombres_rutas, TotalInteractions), 
               names_to = "Target", 
               values_to = "Count")

# Calculate the percentage
conteo_long <- conteo_long %>%
  mutate(Percentage = (Count / TotalInteractions) * 100)

# Create the heatmap
ggplot(conteo_long, aes(x = Target, y = nombres_rutas, fill = Percentage)) +
  geom_tile() + # Use geom_tile for heatmap squares
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "white", 
                      name = "Percentage", labels = percent_format()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(face = "bold")) +
  labs(fill = "Interaction\nPercentage")






# Transformar el dataframe a formato largo
conteo_largo <- conteo %>%
  gather(key = "Diana", value = "Interacciones", -nombres_rutas)

# Calcular el total de interacciones por diana
totales_dianas <- conteo_largo %>%
  group_by(Diana) %>%
  summarise(Total = sum(Interacciones))

# Unir los totales con el conteo largo para calcular los porcentajes
conteo_largo <- conteo_largo %>%
  left_join(totales_dianas, by = "Diana") %>%
  mutate(Porcentaje = (Interacciones / Total) * 100)

# Ahora el dataframe contiene los porcentajes correctos, procedemos a graficar
ggplot(conteo_largo, aes(x = nombres_rutas, y = Porcentaje, fill = Diana)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Rutas", y = "Porcentaje de Interacciones", fill = "Diana") +
  coord_flip() # Girar el gráfico para mejor visualización




#HEATMAP PORCENTAJE DE INTERACCIONES DE DIANAS RUTA

# Transformar el dataframe a formato largo
conteo_largo <- conteo %>%
  gather(key = "Diana", value = "Interacciones", -nombres_rutas)

# Calcular el total de interacciones por diana
totales_dianas <- conteo_largo %>%
  group_by(Diana) %>%
  summarise(Total = sum(Interacciones))

# Unir los totales con el conteo largo para calcular los porcentajes
conteo_largo <- conteo_largo %>%
  left_join(totales_dianas, by = "Diana") %>%
  mutate(Porcentaje = (Interacciones / Total) * 100)

# Ordenar las dianas por el total de interacciones
# Crear un factor ordenado para Diana basado en la suma de interacciones
dianas_ordenadas <- totales_dianas %>%
  arrange(desc(Total)) %>%
  .$Diana

conteo_largo$Diana <- factor(conteo_largo$Diana, levels = dianas_ordenadas)

# Crear un gráfico de calor con las dianas ordenadas
ggplot(conteo_largo, aes(x = Diana, y = nombres_rutas, fill = Porcentaje)) +
  geom_tile() + # Usar geom_tile para crear el gráfico de calor
  scale_fill_gradient2(low = "white", high = "blue", mid = "white", 
                       midpoint = 0, limit = c(0,100), space = "Lab", 
                       name="Porcentaje de\nInteracciones") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Diana", y = "Rutas")


# Calculando las interacciones totales por diana
interacciones_totales_diana <- conteo_largo %>%
  group_by(Diana) %>%
  summarise(TotalInteracciones = sum(Interacciones)) %>%
  ungroup() %>%
  arrange(desc(TotalInteracciones))

# Gráfico de barras para las interacciones totales por diana, ordenadas y en color gris
ggplot(interacciones_totales_diana, aes(x = reorder(Diana, -TotalInteracciones), y = TotalInteracciones)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  labs(x = "Diana", y = "Total de Interacciones", title = "Interacciones Totales por Diana")



# Calcular el total de interacciones
total_interacciones <- sum(interacciones_totales_diana$TotalInteracciones)

# Calcular porcentaje de interacciones por diana respecto al total
interacciones_totales_diana <- interacciones_totales_diana %>%
  mutate(Porcentaje = (TotalInteracciones / total_interacciones) * 100) %>%
  arrange(desc(Porcentaje))

# Gráfico de barras para las interacciones por diana en porcentaje
ggplot(interacciones_totales_diana, aes(x = reorder(Diana, -Porcentaje), y = Porcentaje, fill = Diana)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  labs(x = "Diana", y = "Porcentaje del Total de Interacciones", title = "Interacciones por Diana en Porcentaje del Total")


# Calcular la suma total de interacciones para cada diana
suma_interacciones_diana <- interacciones_diana_ruta %>%
  group_by(Diana) %>%
  summarise(SumaInteracciones = sum(TotalInteracciones)) %>%
  arrange(desc(SumaInteracciones))

# Crear un factor ordenado para Diana basado en la suma de interacciones
interacciones_diana_ruta$Diana <- factor(interacciones_diana_ruta$Diana, levels = suma_interacciones_diana$Diana)
# Creando el heatmap con las dianas en el eje x y las rutas en el eje y
ggplot(interacciones_diana_ruta, aes(x = Diana, y = nombres_rutas, fill = TotalInteracciones)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "Dianas", y = "Rutas", fill = "Total de Interacciones")

# Transformar el dataframe a formato largo
conteo_largo <- conteo %>%
  gather(key = "Diana", value = "Interacciones", -nombres_rutas) %>%
  # Filtrar para incluir solo las rutas con al menos una interacción
  filter(Interacciones > 0)

# Calcular el número de rutas por diana y ordenar por ese número
rutas_por_diana <- conteo_largo %>%
  group_by(Diana) %>%
  summarise(NumeroRutas = n_distinct(nombres_rutas)) %>%
  arrange(desc(NumeroRutas))

# Crear un gráfico de barras
ggplot(rutas_por_diana, aes(x = reorder(Diana, -NumeroRutas), y = NumeroRutas)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  theme(legend.position = "none") +
  labs(x = "Diana", y = "Número de Rutas", title = "Dianas por ruta") +
  guides(fill=FALSE)  # Eliminar la leyenda para 'fill' ya que no es necesaria




# Transformar el dataframe a formato largo
conteo_largo <- conteo %>%
  gather(key = "Diana", value = "Interacciones", -nombres_rutas) %>%
  # Filtrar para incluir solo las rutas con al menos una interacción
  filter(Interacciones > 0)

# Calcular el número de rutas por diana
rutas_por_diana <- conteo_largo %>%
  group_by(Diana) %>%
  summarise(NumeroRutas = n_distinct(nombres_rutas))

# Calcular el porcentaje de rutas para cada diana
total_rutas <- n_distinct(conteo$nombres_rutas)
rutas_por_diana <- rutas_por_diana %>%
  mutate(PorcentajeRutas = (NumeroRutas / total_rutas) * 100)

# Crear un gráfico de barras con porcentajes
ggplot(rutas_por_diana, aes(x = reorder(Diana, -PorcentajeRutas), y = PorcentajeRutas)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  theme(legend.position = "none") +
  labs(x = "Diana", y = "Porcentaje de Rutas", title = "Porcentaje de Dianas Totales por Ruta") +
  guides(fill=FALSE)




#HEATMAP de miRNAs por ruta

# Primero, uniremos las tablas basándonos en la columna común que es el Target.Gene con Diana
conteo_agregado <- conteo_largo %>%
  inner_join(miRNAs_relevantes, by = c("Diana" = "Target.Gene"))

# Ahora, sumamos el número de dianas relevantes por miRNA
conteo_miRNA <- conteo_agregado %>%
  group_by(miRNA) %>%
  summarise(TotalDianas = n_distinct(Diana))

# También, sumamos el número de circuitos relevantes por ruta
conteo_ruta <- conteo_agregado %>%
  group_by(nombres_rutas) %>%
  summarise(TotalCircuitos = n_distinct(miRNA))

# Preparar los datos para el heatmap
heatmap_data <- conteo_agregado %>%
  group_by(miRNA, nombres_rutas) %>%
  summarise(TotalInteracciones = sum(Interacciones)) %>%
  ungroup() %>%
  complete(miRNA, nombres_rutas, fill = list(TotalInteracciones = 0))

# Crear el heatmap
ggplot(heatmap_data, aes(x = miRNA, y = nombres_rutas, fill = TotalInteracciones)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "miRNA", y = "Ruta", fill = "Total de Interacciones")



#HEATMAP de miRNAs por ruta pero ordenado
# Realizamos el inner join
conteo_agregado <- conteo_largo %>%
  inner_join(miRNAs_relevantes, by = c("Diana" = "Target.Gene"))

# Suma las interacciones para obtener TotalInteracciones
conteo_agregado <- conteo_agregado %>%
  group_by(miRNA, nombres_rutas) %>%
  summarise(TotalInteracciones = sum(Interacciones), .groups = 'drop')

# Calculamos la suma total de interacciones para cada miRNA
suma_interacciones_miRNA <- conteo_agregado %>%
  group_by(miRNA) %>%
  summarise(SumaInteracciones = sum(TotalInteracciones)) %>%
  arrange(desc(SumaInteracciones))

# Nos aseguramos de que miRNA esté en el orden correcto para el heatmap
conteo_agregado$miRNA <- factor(conteo_agregado$miRNA, levels = suma_interacciones_miRNA$miRNA)

# Genera el heatmap con los miRNAs ordenados y sin malla de fondo
ggplot(conteo_agregado, aes(x = miRNA, y = nombres_rutas, fill = TotalInteracciones)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(), # Elimina la malla principal
    panel.grid.minor = element_blank(), # Elimina la malla secundaria
    panel.background = element_rect(fill = "white", colour = "white"), # Fondo blanco
    axis.line = element_line(colour = "black")
  ) +
  labs(x = "miRNA", y = "Ruta", fill = "Total de Interacciones")




# miRNAs POR RUTA
# Sumar las interacciones para cada combinación única de miRNA y ruta
interacciones_miRNA_ruta <- conteo_agregado %>%
  group_by(miRNA, nombres_rutas) %>%
  summarise(TotalInteracciones = sum(TotalInteracciones), .groups = 'drop')  # Sumar y quitar los grupos

# Crear un gráfico de barras para el número de interacciones por combinación de miRNA y ruta
ggplot(interacciones_miRNA_ruta, aes(x = nombres_rutas, y = TotalInteracciones, fill = miRNA)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
    legend.position = "right",
    legend.text = element_text(size = 7) 
  ) +
  labs(x = "Ruta", y = "Total de Interacciones", title = "Interacciones Totales de miRNA por Ruta") +
  scale_fill_discrete(name = "miRNA")


# Agrupar por miRNA y ruta, y luego sumar las interacciones
interacciones_miRNA_ruta <- conteo_agregado %>%
  group_by(miRNA, nombres_rutas) %>%
  summarise(TotalInteracciones = sum(TotalInteracciones), .groups = 'drop')  # Sumar y quitar los grupos

# Preparar los datos para el heatmap
heatmap_data <- interacciones_miRNA_ruta %>%
  complete(miRNA, nombres_rutas, fill = list(TotalInteracciones = 0))

# Crear el heatmap
ggplot(heatmap_data, aes(x = nombres_rutas, y = miRNA, fill = TotalInteracciones)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(x = "Ruta", y = "miRNA", fill = "Total de Interacciones")


# Sumar las interacciones para cada miRNA
total_interacciones_miRNA <- conteo_agregado %>%
  group_by(miRNA) %>%
  summarise(TotalInteracciones = sum(TotalInteracciones), .groups = 'drop') %>%
  arrange(desc(TotalInteracciones))  # Ordenar los miRNAs por el total de interacciones

# Crear un gráfico de barras para el número total de interacciones por miRNA
ggplot(total_interacciones_miRNA, aes(x = miRNA, y = TotalInteracciones, fill = miRNA)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none"
  ) +
  labs(x = "miRNA", y = "Total de Interacciones", title = "Interacciones Totales por miRNA")




# Agrupar los totales por ruta y sumar las interacciones
interacciones_por_ruta <- conteo_agregado %>%
  group_by(nombres_rutas) %>%
  summarise(TotalInteracciones = sum(TotalInteracciones)) %>%
  ungroup() %>%
  arrange(desc(TotalInteracciones))

# Reordenar las rutas para el gráfico
interacciones_por_ruta$nombres_rutas <- factor(interacciones_por_ruta$nombres_rutas, 
                                               levels = interacciones_por_ruta$nombres_rutas)

# Crear un gráfico de barras para las interacciones totales por ruta
ggplot(interacciones_por_ruta, aes(x = nombres_rutas, y = TotalInteracciones)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  theme(legend.position = "none") +
  labs(x = "Ruta", y = "Total de Interacciones", title = "Interacciones Totales de miRNA por Ruta")




# Calcular el percentil 90
percentil90 <- quantile(total_interacciones_miRNA$TotalInteracciones, 0.9)

# Filtrar miRNAs que superan el percentil 90 en interacciones totales
miRNAs_superan_percentil90 <- total_interacciones_miRNA %>%
  filter(TotalInteracciones > percentil90)

# Ver los resultados
print(miRNAs_superan_percentil90$miRNA)

# Filtrar miRNAs_relevantes para obtener los registros correspondientes
dianas_miRNAs_filtradas <- miRNAs_relevantes %>%
  filter(miRNA %in% miRNAs_superan_percentil90$miRNA)

# Unir con totales_diana
dianas_unidas <- dianas_miRNAs_filtradas %>%
  inner_join(totales_dianas, by = c("Target.Gene" = "Diana"))

# Ordenar las dianas según los totales
dianas_ordenadas <- dianas_unidas %>%
  arrange(desc(Total))

# Mostrar las dianas ordenadas
print(unique(dianas_ordenadas$Target.Gene))



# ANÁLISIS FUNCIONAL DE LOS GENES DIANAS CON LOS miRNAs QUE AFECTAN A MÁS RUTAS

# Lista de genes para el análisis
genes <- c("PTBP1", "AHR", "NFKB1", "HCN2", "SPI1", "PAX6", "SERPINB5", 
           "IKBKE", "TAGLN2", "CDH1", "SNCA", "PPP1R13L", "E2F2", "INPP5D", 
           "RELA", "KLF4")

# Mapeo de genes a identificadores de Entrez
gene_list <- bitr(genes, fromType = "SYMBOL", 
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Análisis de enriquecimiento de GO para Procesos Biológicos
ego_BP <- enrichGO(gene = gene_list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Visualización del resultado para Procesos Biológicos
bp_plot <- barplot(ego_BP, showCategory = 10) + 
  ggtitle("Análisis de Enriquecimiento de GO: Procesos Biológicos") +
  theme(plot.title = element_text(hjust = 0.5))
bp_plot




# Análisis de enriquecimiento de GO para Funciones Moleculares
ego_MF <- enrichGO(gene = gene_list$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Visualización del resultado para Funciones Moleculares
mf_plot <- barplot(ego_MF, showCategory = 10) + 
  ggtitle("Análisis de Enriquecimiento de GO: Funciones Moleculares") +
  theme(plot.title = element_text(hjust = 0.5))
mf_plot

ego_result <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE)
barplot(ego_result, showCategory = 10)

dotplot(ego_result) + ggtitle("Análisis de enriquecimiento GO")






















---------------------------------------------------------------------------
# ENRIQUECIMIENTO FUNCIONAL DE LAS DIANAS OBTENIDAS CON DREXML3
---------------------------------------------------------------------------
  
# Realizar el enriquecimiento funcional de GO
ego <- enrichGO(gene = genes_with_true,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Análisis de enriquecimiento de GO para Procesos Biológicos
ego_BP <- enrichGO(gene = genes_with_true,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Visualización del resultado para Procesos Biológicos
bp_plot <- barplot(ego_BP, showCategory = 15) + 
  ggtitle("Análisis de Enriquecimiento de GO: Procesos Biológicos") +
  theme(plot.title = element_text(hjust = 0.5))
bp_plot


# Análisis de enriquecimiento de GO para Funciones Moleculares
ego_MF <- enrichGO(gene = genes_with_true,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Visualización del resultado para Funciones Moleculares
mf_plot <- barplot(ego_MF, showCategory = 15) + 
  ggtitle("Análisis de Enriquecimiento de GO: Funciones Moleculares") +
  theme(plot.title = element_text(hjust = 0.5))
mf_plot



# Crear un dotplot con los resultados
dotplot(ego) + ggtitle("GO Enrichment Analysis")

# Crear un cnetplot con los resultados
cnetplot(ego, foldChange = genes_with_true)

barplot(ego, showCategory=15) # Muestra las 20 principales categorías

---------------------------------------------------------------------------





















# ESTO REALIZAR CUANDO TENGA UN VALOR ESTADISTICO (como log fold change, p-valores, etc.)
# gsea_result <- gseGO(geneList   = genes_with_true,
#                      OrgDb      = org.Hs.eg.db,
#                      ont        = "BP",
#                      nPerm      = 1000,
#                      minGSSize  = 10,
#                      maxGSSize  = 500,
#                      pAdjustMethod = "BH",
#                      pvalueCutoff = 0.05,
#                      verbose    = TRUE)

#runningplot(gsea_result)



# SELECCION DE miRNA con dianas relevantes
tabla_mti <- read.table("tabla_mti_filtrada.tsv", sep = "\t", header = TRUE, quote = "")

# Crear una nueva tabla seleccionando solo las columnas 'miRNA' y 'Target.Gene'
tabla_mti <- tabla_mti %>% 
  select(miRNA, Target.Gene)

# Filtrar para obtener solo aquellas filas donde 'Target.Gene' está en 'genes_with_true'
miRNAs_relevantes <- tabla_mti %>% 
  filter(Target.Gene %in% genes_with_true)

# Seleccionamos los miRNAs unico en una nueva lista
miRNA_list <- unique(miRNAs_relevantes$miRNA)



# Calculando la frecuencia de cada miRNA
frecuencia_miRNA <- miRNAs_relevantes %>%
  count(miRNA) %>%
  arrange(desc(n))

# Creando un factor con los niveles en el orden de frecuencia
miRNAs_relevantes$miRNA <- factor(miRNAs_relevantes$miRNA, levels = frecuencia_miRNA$miRNA)

# Gráfico de barras para la frecuencia de miRNAs, ordenado por frecuencia
ggplot(miRNAs_relevantes, aes(x = miRNA)) +
  geom_bar(colour = "white") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) +
  labs(x = "miRNA", y = "Frecuencia", title = "Dianas por miRNA")

# Calcula los conteos de cada miRNA
miRNA_counts <- as.data.frame(table(miRNAs_relevantes$miRNA))

# Crea el histograma usando ggplot2
ggplot(miRNA_counts, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") +
  labs(title = "Histogram of Target Genes per miRNA", 
       x = "Number of Target Genes", 
       y = "Frequency of miRNAs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank())

# Calculando la frecuencia de cada gen diana en miRNAs_relevantes
frecuencia_genes <- miRNAs_relevantes %>%
  count(Target.Gene) %>%
  arrange(desc(n))

# Creando un factor con los niveles en el orden de frecuencia
miRNAs_relevantes$Target.Gene <- factor(miRNAs_relevantes$Target.Gene, levels = frecuencia_genes$Target.Gene)

# Frecuencia de los genes dianas relevantes
ggplot(miRNAs_relevantes, aes(x = Target.Gene)) +
  geom_bar(colour = "white") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank()) + 
  labs(x = "Target Gene", y = "Frequency", title = "miRNA por Diana")


# Calcula los conteos de cada miRNA
target_gene_counts <- as.data.frame(table(miRNAs_relevantes$Target.Gene))

# Crea el histograma usando ggplot2
ggplot(target_gene_counts, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "gray", color = "black") +
  labs(title = "Histogram of miRNAs per Target Gene", 
       x = "Number of miRNAs", 
       y = "Frequency of Target Gene") +
  theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank())

# Crear un boxplot para visualizar la distribución de las interacciones totales por miRNA
ggplot(total_interacciones_miRNA, aes(x = "", y = TotalInteracciones)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribución de Interacciones Totales por miRNA", 
       y = "Total de Interacciones", 
       x = "")




# Usar summary() para obtener un resumen estadístico básico
resumen_estadistico <- summary(total_interacciones_miRNA$TotalInteracciones)

estadisticas_detalladas <- total_interacciones_miRNA %>%
  summarise(
    Media = mean(TotalInteracciones),
    Mediana = median(TotalInteracciones),
    DesviacionEstandar = sd(TotalInteracciones),
    Minimo = min(TotalInteracciones),
    Maximo = max(TotalInteracciones),
    Cuartil1 = quantile(TotalInteracciones, 0.25),
    Cuartil3 = quantile(TotalInteracciones, 0.75),
    Percentil90 = quantile(TotalInteracciones, 0.9),
    Percentil95 = quantile(TotalInteracciones, 0.95),
    Percentil99 = quantile(TotalInteracciones, 0.99),
    Percentil99.9 = quantile(TotalInteracciones, 0.999)
  )

