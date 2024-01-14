setwd("C:/Users/manur/OneDrive - Universidad Europea/Doctorado/TFM")

 
# url <- "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/9.0/miRTarBase_MTI.xlsx"  # URL del archivo con todas las interacciones miRNA-Target publicadas en formato excel

# Cargar el paquete splitstackshape
library(splitstackshape)

# Cargar el paquete openxlsx
library(openxlsx)

# Ruta al archivo de Excel
ruta_archivo <- "hsa_MTI.xlsx"

# Leer el archivo de Excel
nueva_tabla <- read.xlsx(ruta_archivo, sheet = 1)

# Explorar los datos
head(nueva_tabla)  # Muestra las primeras filas del DataFrame

# Realizar el reemplazo de ';' y '/' por '//' en la columna 'Experiments'
nueva_tabla$Experiments <- gsub("Immunoprecipitaion銊uciferase reporter assay", "Immunoprecipitaion;luciferase reporter assay", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("Luciferase reporter assay and Western blot", "Luciferase reporter assay;Western blot", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("luciferasereporterassay", "Luciferase reporter assay", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("qPCRWestern blot", "qPCR;Western blot", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("Western blot, luciferase assay", "Western blot;Luciferase assay", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("Reporter Assay and Western blot", "Reporter Assay;Western blot", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("Reporter AssayyOtherss", "Reporter Assay;Otherss", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("QRTPCRWestern blot ", "QRTPCR;Western blot ", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("QRTPCRWestern blott blot", "QRTPCR;Western blot", nueva_tabla$Experiments)


nueva_tabla$Experiments <- gsub("//", ";", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("/", ";", nueva_tabla$Experiments)


Western_blot <- c("/Western blot",
                  "Weastern blot",
                  "Weastern blot",
                  "Wstern blot",
                  "Western",
                  "Western blo",
                  "Western blog",
                  "western blot",
                  "Western blot",
                  "Western Blot",
                  "Western blot",
                  "Western blot blo",
                  "Western blot blog",
                  "Western blot Blot",
                  "Western blotting",
                  "Western Blotting",
                  "westernblotting",
                  "Westren blot",
                  "Westren Blot",
                  "Wetsern blot",
                  "Wstern blot",
                  "Western blot",
                  "Western blott",
                  "Western blott blo",
                  "Western blott blog",
                  "Western blott blot",
                  "Western blott Blot",
                  "Western blott blot",
                  "Western blott blotting",
                  "Western blott Blotting",
                  "Western blott blot",
                  "Western blot",
                  "Western blot blo",
                  "Western blot blog",
                  "Western blot blot",
                  "Western blot Blot",
                  "Western blot blot",
                  "Western blot blotting",
                  "Western blot Blotting",
                  "Western blott",
                  "Western blotg",
                  "Western blottt",
                  "Western blotting",
                  "Western blottting",
                  "Western blot ",
                  " Western blot",
                  "Western blot "
                  )


qPCR <- c("AGO2 binding RNA immunoprecipitation qRT-PCR",
          "q-PCR",
          "qRT_PCR",
          "QRTPCR",
          "qrt-pcr",
          "qRT-PCR",
          "Quantitative proteomic approach",
          "Real time PCR",
          "real time Qpcr",
          "real time qRT-PCR",
          "Real time RT-PCR",
          "real-time RT-PCR",
          "RTPCR",
          "RT-PCR",
          "semi-qPCR",
          "semi-qRT-PCR",
          "TaqMan miRNA assay",
          "TaqMan miRNA assay",
          "TaqMan miRNA assay/RT-PCR",
          "qqPCR",
          "qRT-PCR",
          "real time qPCR")

Reporter_assay <- c("B-globin reporter assay",
                    "dual-luciferase reporter assay",
                    "EGFP reporter assay",
                    "EGFR reporter assay",
                    "GFP reporter assay",
                    "GFP reporter",
                    "Gluc reporter assay",
                    "GUS reporter assay",
                    "Immunoflourescence",
                    "Immunofluorescence",
                    "Immunofluorescence analysis",
                    "immunofluorescence assays",
                    "Immunofluorescence microscopy",
                    "Immunofluorescence staining",
                    "Immunofluorescent Assay",
                    "In situ hybridization",
                    "LacZ reporter assay",
                    "LacZ assay",
                    "Luciferase assay",
                    "Luciferase reporter assay",
                    "Luciferase reporter assay",
                    "Luciferase reporter assay and western blot",
                    "Luciferase reporter assay",
                    "Luciferase reporter assayMTT",
                    "luciferase reporter assays",
                    "Luciferasereporterassay",
                    "Reporter assay",
                    "Reporter Assay",
                    "Reporter Assay analysis",
                    "Reporter Assay microscopy",
                    "Reporter Assay staining",
                    "Reporter AssayOther",
                    "Reporter Assays",
                    "Reporter Assayy",
                    "TOP茂卢鈥歛sh/FOP茂卢鈥歛sh reporter assay",
                    "TOP铿俛sh/FOP铿俛sh reporter assay",
                    "YFP expresión",
                    "YFP expression",
                    "Luciferase reporter assa",
                    "luciferase reporter assay",
                    " Reporter assay",
                    " Reporter assay analysis",
                   " Reporter assay and Western blot",
                   " Reporter assay microscopy",
                   " Reporter assay staining",
                    " Reporter assayOtherss",
                   "Reporter assay ",
                   "Reporter assay analysis ",
                   "Reporter assay and Western blot ",
                   "Reporter assay microscopy ",
                   "Reporter assay staining ",
                   "Reporter assayOtherss ",
                   "Reporter assayanalysis", 
                   "Reporter assayandWestern blot", 
                   "Reporter assaymicroscopy", 
                   "Reporter assayOtherss",
                   "Reporter assaystaining"
                   )

Microarray <- c("Micorarray",
              "microarray",
              "PCR array",
              "miR Microarray system")

NGS <- c("Chip",
         "ChIP",
         "ChIP-PCR",
         "chip-seq",
         "ChIP-seq",
         "Genotyping",
         "IlluminaExpressionArrays",
         "NGS immunoprecipitation",
         "NGS-PCR",
         "NGS-seq",
         "rna-seq",
         "RNA-seq",
         "Sequencing",
         "Next Generation NGS (NGS)",
         "Next Generation Sequencing (NGS)",
         "Next Generation NGS")

CLIP_SEQ <- c("CLIP-seq",
              "CLIP-Seq dataset analysis",
              "CLIP-SEQ",
              "HITS-CLIP",
              "PAR-CLIP")

pSILAC <- c("pSILAC",
              "SILAC (Stable Isotope Labeling of Amino acids in Culture)",
              "Psilac")


Others <- c("2-D Gel Electrophoresis (2DGE)",
            "2DGE",
            "3'LIFE",
            "5''RACE",
            "5\"RACE",
            "5\\'RACE",
            "5RACE",
            "5'RACE",
            "7 assay",
            "AGO2 Immunoprecipitation",
            "AGO2 Other",
            "Ago2-IP",
            "Ago2-IP/IgG-IP",
            "Alizarin red S staining",
            "Annexin V-FITC",
            "anoikis assay",
            "apoptosis",
            "ASO assay",
            "CAM assay",
            "Caspase-Glo庐 3/7 assay",
            "CC tissues and cells (C33A, HeLa, CaSki, SiHa, and ME-180)",
            "cell cycle assays",
            "Cell proliferation",
            "Cell proliferation assay",
            "ChIP immunoprecipitation",
            "Chromatin immunoprecipitation",
            "Chromatin Immunoprecipitation",
            "Chromatin Other",
            "Chromogenic in situ hybridization",
            "CLASH",
            "Coimmunoprecipitation",
            "Co-immunoprecipitation",
            "colony formation",
            "Colony formation assay",
            "Communoprecipitaion",
            "DNA methylation analysis",
            "ELISA",
            "EMSA",
            "enzyme-linked immunosorbent assay",
            "enzyme-linked Other assay",
            "FACS",
            "flow",
            "Flow",
            "Flow cytometry",
            "FOP luciferase assay",
            "FOP茂卢鈥歛sh reporter assay",
            "FOP铿俛sh reporter assay",
            "GAGs contents assay",
            "Gluc assay",
            "IgG-IP",
            "immunoblot",
            "Immunoblot",
            "Immunoblot analysis",
            "Immunoblotting",
            "Immunocytochemistry",
            "Immunohistochemical (IHC) staining",
            "Immunohistochemical analysis",
            "Immunohistochemistry",
            "Immunohistochemistry (IHC)",
            "Immunohistochemistry analysis",
            "Immunoprecipitaion",
            "immunoprecipitation",
            "Immuno-precipitation",
            "immunosorbent",
            "immunostaining",
            "Immunostaining",
            "Immuohistochemistry",
            "in vitro cullelar assays",
            "in vivo gain",
            "in vivo gain/loss-of-function experiments",
            "intrarenal expression",
            "LC-MS",
            "LC-MS/MS",
            "loss-of-function experiments",
            "Mass spectrometry",
            "mice xenograft",
            "miR PCR array system",
            "miRNA-masking antisense ODN (miR-Mask) assay",
            "Motility assay",
            "mRNA decay",
            "MS",
            "mtt",
            "MTT",
            "MTT assay",
            "Northern blot",
            "Northern blotting",
            "Other (IHC)",
            "Other (IP)",
            "Other analysis",
            "Other assay",
            "Others (ICC)",
            "Others (IHC)",
            "Others (IP)",
            "Others analysis",
            "Others assay",
            "Otherss",
            "Otherss (ICC)",
            "Otherss (IHC)",
            "Otherss (IP)",
            "Otherss analysis",
            "Otherss assay",
            "Otherss assays",
            "Otherss cytometry",
            "Otherssting",
            "Othersting",
            "Others茂卢鈥歛sh",
            "Others铿俛sh",
            "Otherting",
            "pMIR-REPORT",
            "Protein Immunoblot Analyses",
            "Protein Otherss Analyses",
            "Proteomics",
            "Proteomics",
            "proteomics analysis",
            "Reverse-phase protein array",
            "RISC analysis",
            "RNA-binding protein immunoprecipitation",
            "RISC-IP",
            "RNA immunopercipitation",
            "RNA immunoprecipitation assay (RIP)",
            "RNA immunoprecipitation assay (RIP)",
            "RNA immunoprecipitation assay (RIP)",
            "safranin o staining",
            "safranin o staining/GAGs contents assay",
            "SILAC (Stable Isotope Labeling of Amino acids in Culture)",
            "SILAC (Stable Isotope Labeling of Amino acids in Culture)",
            "silico analysis",
            "TargetScan",
            "To test if miR-141 directly targets the PR transcript, we analyzed four predicted miR-141-binding sites (Figure 4c)",
            "TOP",
            "TOP/FOP luciferase assay",
            "TOP茂卢鈥歛sh",
            "TOP铿俛sh",
            "Translational profiling",
            "transwell insert",
            "TRAP",
            "wound healing assays",
            "Zymography",
            "nan",
            "Others",
            "Other")


library(stringr)

for (elemento in Western_blot) {
    nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "Western blot")
}

for (elemento in qPCR) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "qPCR")
}

for (elemento in Microarray) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "Microarray")
}

for (elemento in NGS) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "NGS")
}

for (elemento in CLIP_SEQ) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "CLIP-SEQ")
}

for (elemento in pSILAC) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "pSILAC")
}

for (elemento in Others) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "Others")
}

for (elemento in Reporter_assay) {
  nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, elemento, "Reporter assay")
}

nueva_tabla$Experiments <- gsub("qPCRWestern blot", "qPCR;Western blott", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub(fixed("Next Generation Sequencing (NGS)"), "NGS", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("Next Generation NGS \\(NGS\\)", "NGS", nueva_tabla$Experiments)
nueva_tabla$Experiments <- gsub("Next Generation Sequencing \\(NGS\\)", "NGS", nueva_tabla$Experiments)
nueva_tabla$Experiments <- str_replace(nueva_tabla$Experiments, "NGS \\(NGS\\)", "NGS")

valid_options <- c("Western blot", "qPCR", "Reporter assay", "Microarray", "NGS", "CLIP-SEQ", "pSILAC")

nueva_tabla$Experiments <- str_extract_all(nueva_tabla$Experiments, paste(valid_options, collapse = "|")) %>%
  lapply(function(x) ifelse(length(x) > 0, paste(x, collapse = ";"), "Other")) %>%
  unlist()

write.xlsx(nueva_tabla, "hsa_MTI_mejorado.xlsx", rowNames = FALSE)

### Apartado de unificación de nombres de experimentos

# Convertir la columna "Experiments" en un vector de caracteres
nueva_tabla$Experiments <- as.character(nueva_tabla$Experiments)

# Dividir los experimentos en una lista separados por ";"
experimentos_lista <- strsplit(nueva_tabla$Experiments, ";", fixed = TRUE)

# Unir todos los experimentos en un vector
Experimentos <- unlist(experimentos_lista)

# Obtener una lista de experimentos diferentes y su frecuencia
experimentos_diferentes <- table(Experimentos)

# Abre el archivo en modo escritura
archivo <- file("numero_experimentos.txt", "w")
sink(archivo)

print(experimentos_diferentes)
sink()
# Cierra el archivo y restaura la salida estándar
close(archivo)


