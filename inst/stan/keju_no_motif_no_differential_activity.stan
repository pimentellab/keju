data {
  int<lower=0> E; // number of TREs (including controls)
  int<lower=0> N; // total number of paired DNA - RNA readouts
  int<lower=1> N_DISPERSION_GROUPS;
  
  array[N] int<lower=1, upper=E> to_tre; 
  array[N] int<lower=1, upper=N_DISPERSION_GROUPS> to_dispersion_index;

  array[N] int<lower=0> R;
  vector[N] d;
  vector[N] S_r;
}

parameters {
  vector<lower=0>[N_DISPERSION_GROUPS] disp_r;
  vector[E] alpha;
}

model {
  alpha ~ normal(0, 1);
  disp_r ~ exponential(1);

  R ~ neg_binomial_2_log(
    log(S_r) +
    log(d) + 
    alpha[to_tre],
    disp_r[to_dispersion_index]
  );
}
