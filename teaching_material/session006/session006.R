# Session 006
# 
# In this session we perform a signature fit analysis on the simulated
# dataset SD006. This dataset contains only one catalogue with
# common breast cancer signatures and two rare signature

# load the required R packages
library(signature.tools.lib)
library(teachingmutationalsignatures)

# set and create an output directory for the analysis
outdir <- "~/TeachingMutationalSignatures/Session006/"
dir.create(outdir,showWarnings = F,recursive = T)

# fetch the mutational catalogue that we will work with in this session
catalogues_SD006 <- fetchData(datasetname = "SD006")

# we know this is a biliary cancer sample so let's fetch the biliary cancer
# mutational signatures, both common and rare
SBS_Biliary_sigs <- getSignaturesForFitting(organ = "Biliary")

# this gives us 7 common signatures and 54 possible rare signatures
# what would happen if we try to use them all?
objFit <- Fit(catalogues = catalogues_SD006,
              signatures = cbind(SBS_Biliary_sigs$common,SBS_Biliary_sigs$rare),
              useBootstrap = TRUE,
              nboot = 200)
# plot the results
plotFitResults(fitObj = objFit,outdir = paste0(outdir,"fitAll/"))

# let's try FitMS instead, which will only add rare signatures when they 
# significantly improve the fit
objFitMS <- FitMS(catalogues = catalogues_SD006,
                  organ = "Biliary",
                  useBootstrap = TRUE,
                  nboot = 200)
# plot the results
plotFitResults(fitObj = objFitMS,outdir = paste0(outdir,"fitMS/"))

# it appears that there are multiple solutions, let's see if therre are 
# multiple rare signatures at once
objFitMS_n2 <- FitMS(catalogues = catalogues_SD006,
                     organ = "Biliary",
                     maxRareSigsPerSample = 2,
                     useBootstrap = TRUE,
                     nboot = 200)
# plot the results
plotFitResults(fitObj = objFitMS_n2,outdir = paste0(outdir,"fitMS_n2/"))

# let's check and compare the exposure performance
expPerf <- checkPerformanceExposuresList(estimated_exposures_list = list(fitAll=objFit$exposures,
                                                                         fitMS=objFitMS$exposures,
                                                                         fitMS_n2=objFitMS_n2$exposures),
                                         datasetname = "SD006",
                                         outfile = paste0(outdir,"performance/SD006exposures.pdf"))

