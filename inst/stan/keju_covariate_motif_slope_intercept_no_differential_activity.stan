data {
  int<lower=0> P; // number of minimal correction effects
  int<lower=0> E; // number of TREs (including controls)
  int<lower=0> N; // total number of paired DNA - RNA readouts

  int<lower=0> N_ALPHA_MOTIF;
  int<lower=1> N_DISPERSION_GROUPS;
  
  
  array[E] int<lower=1, upper=P> tre_to_correction_factor;
  array[E] int<lower=1, upper=N_ALPHA_MOTIF> tre_to_motif;
  array[N] int<lower=1, upper=E> to_tre; 
  array[N] int<lower=1, upper=N_DISPERSION_GROUPS> to_dispersion_index;

  array[N] int<lower=0> R;
  vector[N] d;
  vector[N] S_r;
}

parameters {
  vector[P] intercept;
  vector<lower=-1>[P] slope;
  vector<lower=0>[N_DISPERSION_GROUPS] disp_r;

  vector[E] alpha;
  vector<lower=0>[N_ALPHA_MOTIF] sigma2_alpha;
  vector[N_ALPHA_MOTIF] motif_transcription_rate;
}

model {
  sigma2_alpha ~ inv_gamma(1, 1);
  motif_transcription_rate ~ normal(0, 1);

  slope ~ normal(0, 1);
  intercept ~ normal(0, 1);
  alpha ~ normal((1 + slope[tre_to_correction_factor]) .* motif_transcription_rate[tre_to_motif] + intercept[tre_to_correction_factor],
                 sqrt(sigma2_alpha[tre_to_motif]));

  disp_r ~ exponential(1);

  R ~ neg_binomial_2_log(
    log(S_r) +
    log(d) + 
    alpha[to_tre],
    disp_r[to_dispersion_index]
  );
}
