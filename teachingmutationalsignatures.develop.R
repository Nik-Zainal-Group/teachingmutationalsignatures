#file with commands to set up the teachingmutationalsignatures R package
# 
# devtools::create("teachingmutationalsignatures")

outdir <- "~/Google Drive/My Drive/NikZainal/projects/TeachingMaterial/results/simulatedData/"
SD001 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD001/"))
SD002 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD002/"))
SD003 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD003/"))
SD004 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD004/"))
SD005 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD005/"))
SD006 <- signature.tools.lib::loadSimulatedData(paste0(outdir,"SD006/"))

dataObj <- list()
dataObj[["SD001"]] <- SD001
dataObj[["SD002"]] <- SD002
dataObj[["SD003"]] <- SD003
dataObj[["SD004"]] <- SD004
dataObj[["SD005"]] <- SD005
dataObj[["SD006"]] <- SD006

usethis::use_data(dataObj,
                  internal = T,
                  overwrite = T)

devtools::document()
devtools::install()




