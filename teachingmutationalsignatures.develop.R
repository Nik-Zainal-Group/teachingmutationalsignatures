#file with commands to set up the teachingmutationalsignatures R package
# 
# devtools::create("teachingmutationalsignatures")

outdir <- "~/Google Drive/My Drive/NikZainal/projects/TeachingMaterial/results/simulatedData/"
SD001 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD001/"))
SD002 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD002/"))
SD003 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD003/"))

dataObj <- list()
dataObj[["SD001"]] <- SD001
dataObj[["SD002"]] <- SD002
dataObj[["SD003"]] <- SD003

usethis::use_data(dataObj,
                  internal = T,
                  overwrite = T)

devtools::document()
devtools::install()




