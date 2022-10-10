# Session 003
# 
# In this session we perform a signature extraction analysis on the simulated
# dataset SD003. This dataset contains both common and rare breast cancer signatures.

# load the required R packages
library(signature.tools.lib)
library(teachingmutationalsignatures)

# set and create an output directory for the analysis
outdir <- "~/TeachingMutationalSignatures/Session003/"
dir.create(outdir,showWarnings = F,recursive = T)

# fetch the mutational catalogues that we will work with in this session
catalogues_SD003 <- fetchData(datasetname = "SD003")

# we begin by exploring the data using clustering
resCl_SD003 <- cataloguesClustering(catalogues_SD003,nclusters = 1:15,
                                    outdir = paste0(outdir,"cataloguesClustering/"))

# Alternatively, let's try to remove samples that may have are signatures first
outdir_sel <- paste0(outdir,"with_selected_samples/")
nclustersSelection <- "9"
clusters_to_keep <- c(1:4,6,7,9)
samplesSelected <- rownames(resCl_SD003$clusters_table)[resCl_SD003$clusters_table[,nclustersSelection] %in% clusters_to_keep]
SignatureExtraction(cat = catalogues_SD003[,samplesSelected],
                    outFilePath = paste0(outdir_sel,"Extraction/"),
                    nrepeats = 25,nboots = 4,filterBestOfEachBootstrap = T,
                    nparallel = 4,nsig = 8:12,plotResultsFromAllClusteringMethods = F,
                    parallel = T)

# load the optimal signatures (9)
estimated_signatures_SD003_sel <- readTable(paste0(outdir_sel,"Extraction/round_1/sig_9/Sigs_plot_extraction_ns9_nboots4.tsv"))

# is there any sample that was not fully explained by the signatures we found?
unexplSamples_SD003_sel <- unexplainedSamples(outfileRoot = paste0(outdir_sel,"unexplained/SD003"),
                                                                   catalogues = catalogues_SD003,
                                                                   sigs = estimated_signatures_SD003_sel,
                                                                   nmuts_threshold = 300,
                                                                   pvalue_threshold = 0.01)
# yes, 4 samples, get their residual
significant_residuals <- unexplSamples_SD003_sel$unexplSamples_residuals
# now we cluster the residuals
resCl_residuals_SD003_sel <- cataloguesClustering(significant_residuals,
                                                  nclusters = 1:3,
                                                  outdir = paste0(outdir_sel,"residualClustering/"))
# we decide that there are 2 clusters, each having 2 samples in it
# now we need to extract the rare signatures
resRareSigs_SD003_sel <- rareSignatureExtraction(outfileRoot = paste0(outdir_sel,"ExtractionRare/SD003"),
                                                 catalogues = catalogues_SD003,
                                                 residuals = unexplSamples_SD003_sel$all_residuals,
                                                 unexpl_samples = unexplSamples_SD003_sel$unexplSamples,
                                                 clusters = resCl_residuals_SD003_sel$clusters_table[,"2"],
                                                 useclusters = list(c(1),c(2)),
                                                 commonSignatures = estimated_signatures_SD003_sel,
                                                 commonExposures = unexplSamples_SD003_sel$exposures)
# let's compare with a limiting maxiter for cluster 2
resRareSigs_SD003_sel_maxiter <- rareSignatureExtraction(outfileRoot = paste0(outdir_sel,"ExtractionRare_maxiter/SD003"),
                                                         catalogues = catalogues_SD003,
                                                         residuals = unexplSamples_SD003_sel$all_residuals,
                                                         unexpl_samples = unexplSamples_SD003_sel$unexplSamples,
                                                         clusters = resCl_residuals_SD003_sel$clusters_table[,"2"],
                                                         useclusters = list(c(1),c(2)),
                                                         maxiter = c(1000,100),
                                                         commonSignatures = estimated_signatures_SD003_sel,
                                                         commonExposures = unexplSamples_SD003_sel$exposures)

# check whether the signatures have been extracted correctly
perfSigs_SD003 <- checkPerformanceSignaturesList(estimated_signatures_list = list(selectedsamples = resRareSigs_SD003_sel$commonAndRareSignatures,
                                                                                  selectedsamplesMaxiter = resRareSigs_SD003_sel_maxiter$commonAndRareSignatures),
                                                 datasetname = "SD003",
                                                 outfile = paste0(outdir,"performance/SD003_extraction_signatures.pdf"))

# now let's get the exposures
resFinalExpo_SD003_sel <- finaliseCommonRareSignatureExposures(outfileRoot = paste0(outdir_sel,"ExtractionRare/SD003"),
                                                               catalogues = catalogues_SD003,
                                                               commonSigs = estimated_signatures_SD003_sel,
                                                               listofsignatures = resRareSigs_SD003_sel_maxiter$listofsignatures,
                                                               listofsamples = resRareSigs_SD003_sel_maxiter$listofsamples,
                                                               nboot = 50,nparallel = 4)

# check if the exposures and signatures assignments to samples are correct
estimated_exposures_SD003_sel <- t(resFinalExpo_SD003_sel$fitWithRare$exposures)
# we need to correct the names
colnames(estimated_exposures_SD003_sel) <- updateSigNamesWithMatchTable(signames = colnames(estimated_exposures_SD003_sel),
                                                                        matchTable = perfSigs_SD003$perfList$selectedsamplesMaxiter$matchTable)
# now the performance
perfExp_SD003 <- checkPerformanceExposuresList(estimated_exposures_list = list(selectedsamplesMaxiter = estimated_exposures_SD003_sel),
                                               datasetname = "SD003",
                                               outfile = paste0(outdir,"performance/SD003_extraction_exposures.pdf"))


