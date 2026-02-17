# Functions to process keju
#' @import basilisk

## usethis namespace: start
#' @importFrom basilisk basiliskRun
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
## usethis namespace: end
NULL

#' Create keju object from count data. Counts and metadata should be matched by barcode and batch.
#'
#' @param barcode a length N vector of barcodes
#' @param R a length N vector of RNA counts
#' @param D a length N vector of DNA counts
#' @param batch a length N vector denoting RNA batches
#' @param dna_batch a length N vector denoting DNA batches
#' @param architecture a length N vector of enhancer architectures or names. Transcription rate and effect size estimates will be fit on these.
#' @param treatment a length N vector detailing the treatments.
#' @param is_control_treatment a length N boolean vector denoting the control treatment.
#' @param is_control_architecture a length N boolean vector denoting architectures that are negative controls.
#' @param covariates a length N vector denoting covariates to be used in effect size correction and potentially slope/intercept fitting. If none, will default to a vector of ones of length N, to perform global effect size correction.
#' @param motif a length N vector denoting motifs shared across architectures. Used to shrink architecture-level estimates to a shared motif-level estimate. Necessary for motif_shrinkage and covariate_motif_slope_intercept models.
#' @returns an initialized keju object with a counts object
#' @export
keju_from_counts <- function(barcode,
                      R,
                      D,
                      batch,
                      dna_batch,
                      architecture,
                      treatment,
                      is_control_treatment = FALSE,
                      is_control_architecture = FALSE,
                      covariates = 'correction',
                      motif = NULL
) {
    counts <- data.frame(
                barcode = barcode,
                R = R,
                D = D,
                batch = batch,
                dna_batch = dna_batch,
                architecture = architecture,
                treatment = treatment,
                is_control_treatment = is_control_treatment,
                is_control_architecture = is_control_architecture,
                covariates = covariates
    )

    if (!is.null(motif)) {
        counts$motif <- motif
    }

    keju <- list(counts = counts)
    return(keju)
}

#' Filters architectures without at least an average of barcode_threshold barcodes per batch with DNA counts > dna_threshold and RNA counts > rna_threshold.
#'
#' @param keju a keju object containing counts
#' @param dna_threshold minimum DNA counts per barcode
#' @param rna_threshold minimum RNA counts per barcode
#' @param barcode_threshold minimum average barcodes per batch
#' @returns a keju object with the filtered_counts object
#' @export
filter <- function(keju, dna_threshold=5, rna_threshold=5, barcode_threshold=10) {
    counts <- keju$counts

    message(paste0('Converting ', sum(is.na(counts)), ' NA values to zero'))
    counts[is.na(counts)] = 0
    
    n_batch <- length(unique(counts$batch))
    threshold_number_of_barcodes <- n_batch * barcode_threshold # assessing average per-batch

    counts$thresholded <- (counts$D < dna_threshold) | (counts$R < rna_threshold)
    counts$to_keep <- !counts$thresholded

    architecture_threshold_map <- counts %>% 
                                    group_by(architecture) %>% 
                                    summarize(barcode_passed_threshold = sum(to_keep),
                                      thresholded = sum(thresholded))
    architecture_threshold_map$to_keep <- architecture_threshold_map$barcode_passed_threshold > threshold_number_of_barcodes

    message(paste0('Filtering ', sum(!architecture_threshold_map$to_keep), ' architectures with insufficient barcodes.'))
    to_keep <- architecture_threshold_map[architecture_threshold_map$to_keep, ]$architecture
    keju$filtered_counts <- counts[counts$architecture %in% to_keep, ]
    return(keju)
}

#' Process data to fit keju model without motif-level shrinkage. Prefers to run on filtered_counts, but will run on raw counts.
#'
#' @param keju a keju object containing counts
#' @param G number of enhancers per dispersion parameter
#' @param infer_differential_activity boolean to use effect size estimation with a control and alternate treatment
#' @returns a keju object processed for the no_motif model
#' @export
use_no_motif <- function(keju, 
                    G=50,
                    infer_differential_activity=FALSE
) {
    # assert keju has counts
    if (is.null(keju$counts) & is.null(keju$filtered_counts)) {
        stop('No counts object or filtered_counts object detected. Please run keju_from_counts to create a counts first, and optionally filter counts using the built-in filtering function.')
    }

    # py_require("numpy")
    # py_require("pandas")
    # py_require("formulaic")
    proc <- basiliskStart(processing_env)
    on.exit(basiliskStop(proc))

    keju <- basiliskRun(proc, fun=function(keju, G, infer_differential_activity) {
        source_python(system.file('python', 'process.py', package='keju'))
        return(py_process(keju, G, infer_differential_activity)) 
    }, keju=keju, G=G, infer_differential_activity=infer_differential_activity)

    return(keju)
}

#' Process counts to fit keju model with motif-level shrinkage. 
#'
#' @param keju a keju object containing counts or filtered_counts object.
#' @param G number of enhancers per dispersion parameter
#' @param infer_differential_activity boolean to use effect size estimation with a control and alternate treatment
#' @returns a keju object processed for the motif_shrinkage model
#' @export
use_motif_shrinkage <- function(keju,
                                G=50,
                                infer_differential_activity=FALSE

) {
    # py_require("numpy")
    # py_require("pandas")
    # py_require("formulaic")

    proc <- basiliskStart(processing_env)
    on.exit(basiliskStop(proc))

    keju <- basiliskRun(proc, fun=function(keju, G, infer_differential_activity) {
        source_python(system.file('python', 'process.py', package='keju'))
        keju <- use_no_motif(keju, G, infer_differential_activity)

        if(is.null(keju$df$motif)) {
            stop('No motif information provided. Please set the motif parameter in keju_from_counts.')
        }
        keju <- py_use_motif_shrinkage(keju, infer_differential_activity)

        return(keju) 
    }, keju=keju, G=G, infer_differential_activity=infer_differential_activity)
    
    return(keju)
}

#' Process data to fit keju model with motif-level shrinkage and covariate-level effects on transcription rate.
#'
#' @param keju a keju object containing counts or filtered_counts object.
#' @param G number of enhancers per dispersion parameter
#' @param infer_differential_activity boolean to use effect size estimation with a control and alternate treatment
#' @returns a keju object processed for the covariate_motif_slope_intercept model
#' @export
use_covariate_slope_intercept <- function(keju,
                                          G=50,
                                          infer_differential_activity=FALSE
) {
    proc <- basiliskStart(processing_env)
    on.exit(basiliskStop(proc))

    keju <- basiliskRun(proc, fun=function(keju, G, infer_differential_activity) {
        source_python(system.file('python', 'process.py', package='keju'))
            keju <- use_motif_shrinkage(keju, G, infer_differential_activity)
            keju <- py_use_covariate_slope_intercept(keju)
        return(keju) 
    }, keju=keju, G=G, infer_differential_activity=infer_differential_activity)
    
    return(keju)
}
