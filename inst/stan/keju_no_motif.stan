data {
  int<lower=0> P; // number of minimal correction effects
  int<lower=0> E; // number of TREs (including controls)
  int<lower=0> N; // total number of paired DNA - RNA readouts

  int<lower=1> N_EFFECTS;
  int<lower=1> N_DISPERSION_GROUPS;
  
  array[N] int<lower=1, upper=E> to_tre; 
  array[N] int<lower=1, upper=N_DISPERSION_GROUPS> to_dispersion_index;

  array[N] int<lower=0> R;
  vector[N] d;
  vector[N] S_r;

  int n_w;
  int n_v;
  int n_u;
  vector[n_w] Xw;
  array[n_v] int Xv;
  array[n_u] int Xu;
  int Xm;
  int Xn;

  int y_n_w;
  int y_n_v;
  int y_n_u;
  vector[y_n_w] Yw;
  array[y_n_v] int Yv;
  array[y_n_u] int Yu;
  int Ym;
  int Yn;
}

parameters {
  vector[P] beta_correction;
  vector<lower=0>[N_DISPERSION_GROUPS] disp_r;
  vector[E] alpha;
  vector[N_EFFECTS] beta;
}

model {
  alpha ~ normal(0, 1);
  beta_correction ~ normal(0, 1);
  beta ~ normal(0, 1);
  disp_r ~ exponential(1);

  R ~ neg_binomial_2_log(
    log(S_r) +
    log(d) + 
    alpha[to_tre] + 
    csr_matrix_times_vector(Xm, Xn, Xw, Xv, Xu, beta) +
    csr_matrix_times_vector(Ym, Yn, Yw, Yv, Yu, beta_correction), 
    disp_r[to_dispersion_index]
  );
}
