# Session 004
# 
# In this session we perform a signature fit analysis on the simulated
# dataset SD004. This dataset contains only one catalogue with only
# common breast cancer signatures

# load the required R packages
library(signature.tools.lib)
library(teachingmutationalsignatures)

# set and create an output directory for the analysis
outdir <- "~/TeachingMutationalSignatures/Session004/"
dir.create(outdir,showWarnings = F,recursive = T)

# fetch the mutational catalogue that we will work with in this session
catalogues_SD004 <- fetchData(datasetname = "SD004")

# we know this is a breast cancer sample so let's fetch the breast cancer
# mutational signatures, both common and rare
SBS_Breast_sigs <- getSignaturesForFitting(organ = "Breast")

# this gives us 9 common signatures and 54 possible rare signatures
# what would happen if we try to use them all?
objFit <- Fit(catalogues = catalogues_SD004,
              signatures = cbind(SBS_Breast_sigs$common,SBS_Breast_sigs$rare),
              useBootstrap = TRUE,
              nboot = 200)
# plot the results
plotFitResults(fitObj = objFit,outdir = paste0(outdir,"fitAll/"))

# let's try FitMS instead, which will only add rare signatures when they 
# significantly improve the fit
objFitMS <- FitMS(catalogues = catalogues_SD004,
                  commonSignatures = SBS_Breast_sigs$common,
                  rareSignatures = SBS_Breast_sigs$rare,
                  useBootstrap = TRUE,
                  nboot = 200)
# plot the results
plotFitResults(fitObj = objFitMS,outdir = paste0(outdir,"fitMS/"))

# let's check and compare the exposure perfomance
expPerf <- checkPerformanceExposuresList(estimated_exposures_list = list(fitAll=objFit$exposures,
                                                                         fitMS=objFitMS$exposures),
                                         datasetname = "SD004",
                                         outfile = paste0(outdir,"performance/SD004exposures.pdf"))

