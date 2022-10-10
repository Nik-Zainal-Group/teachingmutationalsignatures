# Session 002
# 
# In this session we perform a signature extraction analysis on the simulated
# dataset SD002. This dataset contains both common and rare breast cancer signatures.

# load the required R packages
library(signature.tools.lib)
library(teachingmutationalsignatures)

# set and create an output directory for the analysis
outdir <- "~/TeachingMutationalSignatures/Session002/"
dir.create(outdir,showWarnings = F,recursive = T)

# fetch the mutational catalogues that we will work with in this session
catalogues_SD002 <- fetchData(datasetname = "SD002")

# we begin by exploring the data using clustering
resCl_SD002 <- cataloguesClustering(catalogues_SD002,nclusters = 1:15,
                                    outdir = paste0(outdir,"cataloguesClustering/"))

# let's attempt an extraction with all samples
outdir_all <- paste0(outdir,"with_all_samples/")
SignatureExtraction(cat = catalogues_SD002,
                    outFilePath = paste0(outdir_all,"Extraction/"),
                    nrepeats = 50,nboots = 4,filterBestOfEachBootstrap = T,
                    nparallel = 4,nsig = 12:16,plotResultsFromAllClusteringMethods = F,
                    parallel = T)
# load the optimal signatures (14)
estimated_signatures_SD002_all <- readTable(paste0(outdir_all,"Extraction/round_1/sig_14/Sigs_plot_extraction_ns14_nboots4.tsv"))

# is there any sample that was not fully explained by the signatures we found?
unexplSamples_SD002_all <- unexplainedSamples(outfileRoot = paste0(outdir_all,"unexplained/SD002"),
                                              catalogues = catalogues_SD002,
                                              sigs = estimated_signatures_SD002_all,
                                              nmuts_threshold = 300,
                                              pvalue_threshold = 0.05)

# now fit the signatures to the samples
resFit_SD002_all <- Fit(catalogues = catalogues_SD002,
                        signatures = estimated_signatures_SD002_all,
                        useBootstrap = T,
                        nboot = 50,
                        threshold_percent = 5,
                        threshold_p.value = 0.05,
                        nparallel = 4)
plotFitResults(fitObj = resFit_SD002_all,outdir = paste0(outdir_all,"sigfit/"))
saveFitToFile(fitObj = resFit_SD002_all,
              filename = paste0(outdir_all,"sigfit/SD001_fit.rData"))

# Alternatively, let's try to remove samples that may have are signatures first
outdir_sel <- paste0(outdir,"with_selected_samples/")
nclustersSelection <- "9"
clusters_to_keep <- c(1:5)
samplesSelected <- rownames(resCl_SD002$clusters_table)[resCl_SD002$clusters_table[,nclustersSelection] %in% clusters_to_keep]
SignatureExtraction(cat = catalogues_SD002[,samplesSelected],
                    outFilePath = paste0(outdir_sel,"Extraction/"),
                    nrepeats = 25,nboots = 4,filterBestOfEachBootstrap = T,
                    nparallel = 4,nsig = 8:10,plotResultsFromAllClusteringMethods = F,
                    parallel = T)

# load the optimal signatures (9)
estimated_signatures_SD002_sel <- readTable(paste0(outdir_sel,"Extraction/round_1/sig_9/Sigs_plot_extraction_ns9_nboots4.tsv"))

# is there any sample that was not fully explained by the signatures we found?
unexplSamples_SD002_sel <- unexplainedSamples(outfileRoot = paste0(outdir_sel,"unexplained/SD002"),
                                                                   catalogues = catalogues_SD002,
                                                                   sigs = estimated_signatures_SD002_sel,
                                                                   nmuts_threshold = 300,
                                                                   pvalue_threshold = 0.05)
# yes, 5 samples, get their residual
significant_residuals <- unexplSamples_SD002_sel$unexplSamples_residuals
# now we cluster the residuals
resCl_residuals_SD002_sel <- cataloguesClustering(significant_residuals,
                                                  nclusters = 1:5,
                                                  outdir = paste0(outdir_sel,"residualClustering/"))
# we decide that there are 5 clusters, each having a single sample in it
# now we need to extract the rare signatures
resRareSigs_SD002_sel <- rareSignatureExtraction(outfileRoot = paste0(outdir_sel,"ExtractionRare/SD002"),
                                                 catalogues = catalogues_SD002,
                                                 residuals = unexplSamples_SD002_sel$all_residuals,
                                                 unexpl_samples = unexplSamples_SD002_sel$unexplSamples,
                                                 clusters = resCl_residuals_SD002_sel$clusters_table[,"5"],
                                                 useclusters = list(c(1),c(2),c(3),c(4),c(5)),
                                                 commonSignatures = estimated_signatures_SD002_sel,
                                                 commonSigsToIgnore = list(NA,"S5",NA,NA,NA),
                                                 commonExposures = unexplSamples_SD002_sel$exposures)

resFinalExpo_SD002_sel <- finaliseCommonRareSignatureExposures(outfileRoot = paste0(outdir_sel,"ExtractionRare/SD002"),
                                                               catalogues = catalogues_SD002,
                                                               commonSigs = estimated_signatures_SD002_sel,
                                                               listofsignatures = resRareSigs_SD002_sel$listofsignatures,
                                                               listofsamples = resRareSigs_SD002_sel$listofsamples,
                                                               nboot = 50,nparallel = 4)



# check whether the signatures have been extracted correctly
perfSigs_SD002 <- checkPerformanceSignaturesList(estimated_signatures_list = list(allsamples = estimated_signatures_SD002_all,
                                                                                  selectedsamples = resRareSigs_SD002_sel$commonAndRareSignatures),
                                                 datasetname = "SD002",
                                                 outfile = paste0(outdir,"performance/SD002_extraction_signatures.pdf"))
# check if the exposures and signatures assignments to samples are correct
estimated_exposures_SD002_all <- resFit_SD002_all$exposures
estimated_exposures_SD002_sel <- t(resFinalExpo_SD002_sel$fitWithRare$exposures)
# we need to correct the names
colnames(estimated_exposures_SD002_all) <- updateSigNamesWithMatchTable(signames = colnames(estimated_exposures_SD002_all),
                                                                        matchTable = perfSigs_SD002$perfList$allsamples$matchTable)
colnames(estimated_exposures_SD002_sel) <- updateSigNamesWithMatchTable(signames = colnames(estimated_exposures_SD002_sel),
                                                                        matchTable = perfSigs_SD002$perfList$selectedsamples$matchTable)
# now the performance
perfExp_SD002 <- checkPerformanceExposuresList(estimated_exposures_list = list(allsamples = estimated_exposures_SD002_all,
                                                                               selectedsamples = estimated_exposures_SD002_sel),
                                               datasetname = "SD002",
                                               outfile = paste0(outdir,"performance/SD002_extraction_exposures.pdf"))


