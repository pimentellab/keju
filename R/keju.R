# Functions to run keju
#' @import cmdstanr

## usethis namespace: start
#' @importFrom cmdstanr cmdstan_model
## usethis namespace: end
NULL


#' Fits keju model to data.
#'
#' @param keju a filtered and processed keju object
#' @param output_dir directory for output files
#' @param model model type. Must be one of 'no_motif', 'motif_shrinkage', or 'covariate_motif_slope_intercept.
#' @param infer_differential_activity boolean to use effect size estimation with a control and alternate treatment
#' @param seed random seed for MCMC sampling. See https://mc-stan.org/cmdstanr/reference/model-method-sample.html
#' @param chains number of Markov chains to use.  See https://mc-stan.org/cmdstanr/reference/model-method-sample.html
#' @param parallel_chains maximum number of MCMC chains to run in parallel. See https://mc-stan.org/cmdstanr/reference/model-method-sample.html
#' @returns a fitted keju object
#' @export
keju_fit <- function(
    keju,
    output_dir,
    model = 'no_motif',
    infer_differential_activity=FALSE,
    seed = 1,
    chains = 4,
    parallel_chains = 4
) {
    # check model
    if (model %in% c('no_motif', 'motif_shrinkage', 'covariate_motif_slope_intercept')) {
        if (infer_differential_activity) {
            stan_path <- system.file('stan', paste0('keju_', model, '.stan'), package='keju')
        } else {
            stan_path <- system.file('stan', paste0('keju_', model, '_', 'no_differential_activity', '.stan'), package='keju')
        }
    } else {
        stop(paste0('Model ', model, ' not recognized. Must be one of no_motif, motif_shrinkage, or covariate_motif_slope_intercept.'))
    }

    # create output directory and folder with outputs
    if (!dir.exists(output_dir)) {
        warning(paste0("Output directory ", output_dir, " does not exist. Creating...\n"))
        dir.create(output_dir, recursive=T)
    }
    output_folder <- gsub("//", "/", file.path(output_dir, "keju"))
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }

    mod <- cmdstan_model(stan_path)

    fit <- mod$sample(data = keju$data,
                    seed = seed,
                    chains = chains,
                    parallel_chains = parallel_chains,
                    # refresh = 250,
                    show_messages=TRUE,
                    show_exceptions=FALSE
    )

    keju$diagnostics <- fit$diagnostic_summary()
    print(keju$diagnostics)

    fit_summary = fit$summary()
    rownames(fit_summary) = fit_summary$variable
    keju$fit = fit_summary

    write.csv(fit_summary, gsub("//", "/", file.path(output_folder, 'fit_summary.csv')))
    keju <- pretty_summarize(keju, model, output_folder, infer_differential_activity)
    return(keju)
}

#' Processes parameter estimates using existing metadata from the keju object. Calls significance for effect sizes and writes to output folder.
#'
#' @param keju a fitted keju object
#' @param model the type of model fit by keju
#' @param output_folder folder for output estimate csv files
#' @param infer_differential_activity boolean to use effect size estimation with a control and alternate treatment
pretty_summarize <- function(
                             keju,
                             model,
                             output_folder,
                             infer_differential_activity
) {
    fit <- keju$fit
    tres <- keju$tres
    effects <- keju$effects
    correction_effects <- keju$correction_effects
    
    alphas <- fit[startsWith(rownames(fit), 'alpha['),]
    alphas$architecture <- keju$tres
    rownames(alphas) <- keju$tres
    keju$alphas_estimate <- alphas
    write.csv(alphas, gsub("//", "/", file.path(output_folder, 'alphas.csv')))

    if (infer_differential_activity) {
        betas <- fit[startsWith(rownames(fit), 'beta['),]
        betas$is_significant = (betas$q5 > 0) | (betas$q95 < 0)
        betas$architecture <- keju$effects
        rownames(betas) <- keju$effects
        keju$betas_estimate <- betas
        write.csv(betas, gsub("//", "/", file.path(output_folder, 'betas.csv')))   

        covariates <- fit[startsWith(rownames(fit), 'beta_correction['),]  
        covariates$name <- keju$covariates
        rownames(covariates) <- keju$covariates
        keju$covariates_estimate <- covariates
        write.csv(covariates, gsub("//", "/", file.path(output_folder, 'covariate_effects.csv')))
    }

    

    if (model != 'no_motif') {
        alpha_motifs <- fit[startsWith(rownames(fit), 'motif_transcription_rate['),]
        alpha_motifs$motif <- keju$alpha_motifs
        rownames(alpha_motifs) <- keju$alpha_motifs
        keju$alpha_motifs_estimate <- alpha_motifs
        write.csv(alpha_motifs, gsub("//", "/", file.path(output_folder, 'motifs_alphas.csv')))

        if (infer_differential_activity) {
            beta_motifs <- fit[startsWith(rownames(fit), 'motif_effect['),]
            beta_motifs$is_significant = (beta_motifs$q5 > 0) | (beta_motifs$q95 < 0)
            beta_motifs$motif <- keju$beta_motifs
            rownames(beta_motifs) <- keju$beta_motifs
            keju$beta_motifs_estimate <- beta_motifs
            write.csv(beta_motifs, gsub("//", "/", file.path(output_folder, 'motif_betas.csv')))
        }

        if (model == 'covariate_motif_slope_intercept') {
            slope <- fit[startsWith(rownames(fit), 'slope['),]
            intercept <- fit[startsWith(rownames(fit), 'intercept['),]

            slope$covariate <- keju$covariates
            intercept$covariate <- keju$covariates

            rownames(slope) <- keju$covariates
            rownames(intercept) <- keju$covariates

            keju$covariate_slope_estimate <- slope
            keju$covariate_intercept_estimate <- intercept

            write.csv(slope, gsub("//", "/", file.path(output_folder, 'covariate_slopes.csv')))
            write.csv(intercept, gsub("//", "/", file.path(output_folder, 'covariate_intercepts.csv')))
        }
    }
    return(keju)
}