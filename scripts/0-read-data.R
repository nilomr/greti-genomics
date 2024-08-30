# Import configuration file
config <- config::get()

# parameters
percent.train <- 0.7
num.threads <- 44
num.eigenvectors <- 20
set.seed(555)

# Input files
immigrant.info <- file.path(config$path$data, "600K_immigrants.fam")
immigrant.vcf <- file.path(config$path$data, "600K_immigrants.vcf")
resident.vcf <- file.path(config$path$data, "600K_residents.vcf")

immigrant.gds <- file.path(config$path$data, "immigrant_recode.gds")
resident.gds <- file.path(config$path$data, "resident_recode.gds")
all.gds <- file.path(config$path$data, "all_recode.gds")

# remove any existing gds files before creating new ones
file.remove(list.files(path = config$path$data, pattern = "*.gds$", full.names = TRUE))
SNPRelate::snpgdsVCF2GDS(immigrant.vcf, immigrant.gds, method = "biallelic.only")
SNPRelate::snpgdsVCF2GDS(resident.vcf, resident.gds, method = "biallelic.only")
gdsfmt::showfile.gds(closeall = TRUE)
# Append the immigrant and resident data
SNPRelate::snpgdsCombineGeno(c(immigrant.gds, resident.gds), all.gds)
