# Session 001
# 
# In this session we perform a signature extraction analysis on the simulated
# dataset SD001. This dataset contains only common breast cancer signatures.

# load the required R packages
library(signature.tools.lib)
library(teachingmutationalsignatures)

# set and create an output directory for the analysis
outdir <- "~/TeachingMutationalSignatures/Session001/"
dir.create(outdir,showWarnings = F,recursive = T)

# fetch the mutational catalogues that we will work with in this session
catalogues_SD001 <- fetchData(datasetname = "SD001")

# we begin by exploring the data using clustering
resCl_SD001 <- cataloguesClustering(catalogues_SD001,nclusters = 1:10,
                                    outdir = paste0(outdir,"cataloguesClustering/"))

# let's attempt an extraction with all samples
SignatureExtraction(cat = catalogues_SD001,
                    outFilePath = paste0(outdir,"Extraction/"),
                    nrepeats = 50,nboots = 4,filterBestOfEachBootstrap = T,
                    nparallel = 4,nsig = 8:10,plotResultsFromAllClusteringMethods = F,
                    parallel = T)
# load the optimal signatures
estimated_signatures_SD001 <- readTable(paste0(outdir,"Extraction/round_1/sig_9/Sigs_plot_extraction_ns9_nboots4.tsv"))

# is there any sample that was not fully explained by the signatures we found?
unexplSamples_SD001 <- unexplainedSamples(outfileRoot = paste0(outdir,"unexplained/SD001"),
                                          catalogues = catalogues_SD001,
                                          sigs = estimated_signatures_SD001,
                                          nmuts_threshold = 300,
                                          pvalue_threshold = 0.05)

# now fit the signatures to the samples
resFit_SD001 <- Fit(catalogues = catalogues_SD001,
                    signatures = estimated_signatures_SD001,
                    useBootstrap = T,
                    nboot = 50,
                    threshold_percent = 5,
                    threshold_p.value = 0.05,
                    nparallel = 4)
plotFitResults(fitObj = resFit_SD001,outdir = paste0(outdir,"sigfit/"))
saveFitToFile(fitObj = resFit_SD001,
              filename = paste0(outdir,"sigfit/SD001_fit.rData"))

# check whether the signatures have been extracted correctly
perfSigs_SD001 <- checkPerformanceSignatures(estimated_signatures = estimated_signatures_SD001,
                                             datasetname = "SD001",
                                             outfile = paste0(outdir,"performance/SD001_extraction_signatures.pdf"))
# check if the exposures and signatures assignments to samples are correct
estimated_exposures_SD001 <- resFit_SD001$exposures
# we need to correct the names
colnames(estimated_exposures_SD001) <- updateSigNamesWithMatchTable(signames = colnames(estimated_exposures_SD001),
                                                                    matchTable = perfSigs_SD001$matchTable)
# now the performance
perfExp_SD001 <- checkPerformanceExposures(estimated_exposures = estimated_exposures_SD001,
                                           datasetname = "SD001",
                                           outfile = paste0(outdir,"performance/SD001_extraction_exposures.pdf"))


