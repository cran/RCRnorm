

fitWithPosCtrl = function(y, x)
{
  mod1 = stats::lm(y ~ x)
  coefs = stats::coef(mod1)
  unname(coefs)
}


get_range = function(x, mm = 5)
{
  c(mean(x)-mm*stats::sd(x), mean(x)+mm*stats::sd(x))
}

get_residual = function(log_dat, RNA_conc, coefs)
  #in matrix format
{
  log_dat - sweep(sweep(RNA_conc, 2, coefs[2, ], '*'), 2, coefs[1, ], '+')
}



#'@title An Integrated Regression Model for Normalizing 'NanoString nCounter' Data
#'@description 'NanoString nCounter' is a medium-throughput platform that measures gene or microRNA expression levels.
#'Here is a publication that introduces this platform: Malkov (2009) <doi:10.1186/1756-0500-2-80>. Here is the webpage of NanoString
#'nCounter where you can find detailed information about this platform <https://www.nanostring.com/scientific-content/technology-overview/ncounter-technology>.
#'It has great clinical application, such as diagnosis and prognosis of cancer. This function implements an integrated
#'system of random-coefficient hierarchical regression model for normalizing 'NanoString nCounter' data.
#'It removes noise from the data so that expression levels of genes can be compared
#'across patients.
#'@param dat A list containing data for the 4 probe types: positive control, negative control, housekeeping gene and regular gene.
#'The names for the 4 elements in the list should exactly be: pos_dat, neg_dat, hk_dat and reg_dat, respectively. For an example
#'of the input data format, please refer to the FFPE_dat included in the dataset.
#'The data for each probe type should be a dataframe with rows being genes and column being patients. The number of columns (patients)
#'should be the same for data of all four probe types. The rows of positive control data should have the same order as the postive control
#'RNA amount vector supplied to the function.
#'@param pos_conc A vector of log10 RNA amount of the positive controls. The order of these controls should be the same as the rows of positive control
#'data in dat. The defaut is: log10(c(128, 32, 8, 2, 0.5, 0.125)).
#'@param iter Total number of iterations for Monte Carlo simulation. Default is 1000.
#'@param warmup Number of burnin cycles for Monte Carlo simulation. Default is 500.
#'@param seed Seed for the MCMC sampling for reproducibility. Default is 1.
#'@param random_init Whether to estimate the starting point from data
#'@param all_dat Whether should all data be used to update a_i and b_i.
#'@param mm Number of standard deviations for the prior uniform range.
#'@param m_ab Number of variance for the prior distribution of mu_a and mu_b.
#'@details 'NanoString nCounter' platform includes several internal controls (Positive control; Negative control; Housekeeping genes) to remove noise and normalize data to enable inter-patient
#'gene expression comparasion: 1. removing lane-by-lane experimental variation with positive controls;
#'2. removing background noise introduced by none specific binding with negative controls;
#'3. removing sample loading amount variation or difference in RNA degradation level with housekeeping genes.
#'Our IBMnorm model integrates information from these 3 types of internal controls and get the normalized expression levels of genes we are interested in.
#'Detailed models are in the publication.
#'@return The function returns a list of elements including: summary statistics of key parameters in the model and a list of MCMC samples. The number of MCMC samples
#'equals iter-warmup.
#'@export
#'@examples
#'data(FFPE_dat)
#'result = RCRnorm(FFPE_dat, iter = 20, warmup = 0)
#'@import truncnorm


RCRnorm = function(dat, pos_conc = log10(c(128, 32, 8, 2, 0.5, 0.125)),
                   iter = 8000, warmup = 5000, random_init = F, all_dat = T, seed = 1, mm = 3, m_ab = 9)
{
  ptm <- proc.time()

  set.seed(seed*3723)
  # cat(seed)
  # cat('\n')

  #log10 transform the original count data
  pos_dat = log10(dat$pos_dat + 1)
  neg_dat = log10(dat$neg_dat + 1)
  hk_dat = log10(dat$hk_dat + 1)
  reg_dat = log10(dat$reg_dat + 1)

  #inverse gamma parameter
  u = v = .01
  #number of each class of genes.
  n_hk = dim(hk_dat)[1]
  n_reg = dim(reg_dat)[1]
  n_neg = dim(neg_dat)[1]
  n_pos = dim(pos_dat)[1]
  #number of patients or samples
  n_patient = dim(pos_dat)[2]
  #number of MCMC iteration used to calculate posterior after convergence
  iter_keep = iter - warmup
  #calculate the coefficient for each patient; note: with positive controls.
  # all_coef: First row: a+5se; second row: b+5se; third row: a-5se; fourth row: b-5se.
  all_coef = apply(pos_dat, 2, fitWithPosCtrl, pos_conc)

  mu_a_itm = mean(all_coef[1,])
  mu_b_itm = mean(all_coef[2,])

  #Jacknife to estimate mean and variance of mu_a and mu_b
  mu_a = numeric()
  mu_b = numeric()
  for (i in 1:500)
  {
    mu_a[i] = mean(sample(all_coef[1, ], n_patient - 2))
    mu_b[i] = mean(sample(all_coef[2, ], n_patient - 2))
  }

  mu_a_mu = mean(mu_a)
  cat(mu_a_mu, '\n')
  mu_b_mu = mean(mu_b)
  cat(mu_b_mu, '\n')
  sigma2_mu_a = m_ab * stats::var(mu_a)
  cat(sigma2_mu_a, '\n')
  sigma2_mu_b = m_ab * stats::var(mu_b)
  cat(sigma2_mu_b, '\n')

  # a_L = mean(all_coef[1,]) - mm*stats::sd(all_coef[1,])
  # a_U = mean(all_coef[1,]) + mm*stats::sd(all_coef[1,])
  # b_L = mean(all_coef[2,]) - mm*stats::sd(all_coef[2,])
  # b_U = mean(all_coef[2,]) + mm*stats::sd(all_coef[2,])

  hk_RNA = sweep(sweep(hk_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')
  reg_RNA = sweep(sweep(reg_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')

  ##estimate genes' mean expression level range
  lambda_hk_range = apply(hk_RNA, 1, get_range, mm = mm)
  lambda_reg_range = apply(reg_RNA, 1, get_range, mm = mm)

  #estimate patient effect range by two way ANOVA with patient's regular gene expression level.
  gene = factor(rep(1:n_reg, n_patient))
  patient = factor(rep(1:n_patient, each = n_reg))
  mod = stats::lm(unlist(reg_RNA) ~ patient + gene, contrasts = list(patient = 'contr.sum', gene = 'contr.sum'))

  phi = numeric(n_patient)
  phi[1:(n_patient - 1)] = summary(mod)$coefficients[2:n_patient, 1]
  phi[n_patient] = -sum(phi)

  phi_L = phi - mm * summary(mod)$coefficients[2, 2]
  phi_U = phi + mm * summary(mod)$coefficients[2, 2]


  #initialize all the parameters
  aa = matrix(NA, ncol = n_patient, nrow = iter_keep)
  bb = matrix(NA, ncol = n_patient, nrow = iter_keep)
  cc = numeric(iter_keep)
  phi_return = matrix(NA, ncol = n_patient, nrow = iter_keep) #patient effect
  kappa_hk = matrix(NA, ncol = n_patient * n_hk, nrow = iter_keep)
  kappa_reg = matrix(NA, ncol = n_patient * n_reg, nrow = iter_keep)
  lambda_hk = matrix(NA, ncol = n_hk, nrow = iter_keep)
  lambda_reg = matrix(NA, ncol = n_reg, nrow = iter_keep)
  d_neg = matrix(NA, ncol = n_neg, nrow = iter_keep)
  d_pos = matrix(NA, ncol = n_pos, nrow = iter_keep)
  d_hk = matrix(NA, ncol = n_hk, nrow = iter_keep)
  d_reg = matrix(NA, ncol = n_reg, nrow = iter_keep)
  mu_a = numeric(iter_keep)
  mu_b = numeric(iter_keep)
  sigma2e_neg = numeric(iter_keep)
  sigma2e_phr = numeric(iter_keep)
  sigma2a = numeric(iter_keep)
  sigma2b = numeric(iter_keep)
  sigma2kappa_hk = numeric(iter_keep)
  sigma2kappa_reg = numeric(iter_keep)
  sigma2d_neg = numeric(iter_keep)
  sigma2d_phr = numeric(iter_keep)

  if (random_init == T)
  {
    # mu_a_itm = stats::runif(1, 1, 4)
    # mu_b_itm = stats::runif(1, 0, 2)
    sigma2a_itm = stats::runif(1, 0, .01)
    sigma2b_itm = stats::runif(1, 0, .01)
    a_itm = stats::rnorm(n_patient, 2.5, .1)
    b_itm = stats::rnorm(n_patient, .9, .1)
    cc_itm = stats::runif(1, -6, -1)

    phi_itm = stats::rnorm(n_patient, 0, 2)
    phi_itm[n_patient] = -sum(phi_itm[1:(n_patient - 1)])

    kappa_hk_itm = stats::rnorm(n_hk * n_patient, 0, 1)
    sigma2kappa_hk_itm = stats::runif(1, 0, 1)
    kappa_reg_itm = stats::rnorm(n_reg * n_patient, 0, 1)
    sigma2kappa_reg_itm = stats::runif(1, 0, 1)


    lambda_hk_itm = stats::rnorm(n_hk, 0, 1)
    lambda_reg_itm = stats::rnorm(n_reg, 0, 1)

    d_neg_itm = stats::rnorm(n_neg, 0, .01)
    sigma2d_neg_itm = stats::runif(1, 0, .1)
    d_pos_itm = stats::rnorm(n_pos, 0, .01)
    sigma2d_phr_itm = stats::runif(1, 0, .1)
    d_hk_itm = stats::rnorm(n_hk, 0, .01)
    d_reg_itm = stats::rnorm(n_reg, 0, .01)


    sigma2e_neg_itm = stats::runif(1, 0, .1)
    sigma2e_phr_itm = stats::runif(1, 0, .1)
  }

  #get initial values; itm: intermediate
  if (random_init == F)
  {
    sigma2a_itm = stats::var(all_coef[1,])
    sigma2b_itm = stats::var(all_coef[2,])
    a_itm = all_coef[1,]
    b_itm = all_coef[2,]
    cc_itm = mean(unlist(sweep(sweep(neg_dat, 2, all_coef[1,], '-'), 2, all_coef[2,], '/')))

    phi_itm = phi

    lambda_hk_itm = apply(lambda_hk_range, 2, mean)
    lambda_reg_itm = apply(lambda_reg_range, 2, mean)

    estimate_kappa = sweep(rbind(hk_RNA, reg_RNA), 2, phi_itm, '-')
    estimate_kappa_var = sweep(sweep(rbind(hk_RNA, reg_RNA), 2, phi_itm, '-'), 1, c(lambda_hk_itm, lambda_reg_itm), '-')

    kappa_hk_itm = as.vector(unlist(estimate_kappa[1:n_hk,]))
    sigma2kappa_hk_itm = stats::var(as.vector(unlist(estimate_kappa_var[1:n_hk,])))     #0.02497646
    kappa_reg_itm = as.vector(unlist(estimate_kappa[(1+n_hk):(n_hk+n_reg),]))
    sigma2kappa_reg_itm = stats::var(as.vector(unlist(estimate_kappa_var[(1+n_hk):(n_hk+n_reg),])) )   #0.1220459

    pos_RNA = matrix(rep(pos_conc, n_patient), ncol = n_patient)
    neg_RNA = matrix(rep(cc_itm, n_neg * n_patient), ncol = n_patient)

    d_neg_itm = apply(get_residual(neg_dat, neg_RNA, all_coef), 1, mean)
    sigma2d_neg_itm = stats::var(d_neg_itm)
    d_pos_itm = apply(get_residual(pos_dat, pos_RNA, all_coef), 1, mean)
    sigma2d_phr_itm = stats::var(d_pos_itm)
    d_hk_itm = rep(0, n_hk) #stats::rnorm(n_hk, 0, sqrt(sigma2d_phr_itm))
    d_reg_itm = rep(0, n_reg) #stats::rnorm(n_reg, 0, sqrt(sigma2d_phr_itm))


    sigma2e_neg_itm = stats::var(unlist(sweep(get_residual(neg_dat, neg_RNA, all_coef), 1, d_neg_itm, '-')))
    sigma2e_phr_itm = stats::var(unlist(sweep(get_residual(pos_dat, pos_RNA, all_coef), 1, d_pos_itm, '-')))
  }

  for (i in 1:iter)
  {
    A2 = colSums(sweep(pos_dat - pos_conc %o% b_itm, 1, d_pos_itm, '-'))

    #
    B2 = colSums(sweep(sweep(sweep(pos_dat, 2, a_itm, '-'), 1, d_pos_itm, '-'), 1, pos_conc, '*'))

    # b_itm = stats::rnorm(rep(1, n_patient), ((B1/sigma2e_neg_itm) + (B2)/sigma2e_phr_itm + mu_b_itm/sigma2b_itm)/
    #                 (n_neg*(cc_itm^2)/sigma2e_neg_itm + (sum(pos_conc^2))/sigma2e_phr_itm + 1/sigma2b_itm),
    #               sqrt(1/(n_neg*(cc_itm^2)/sigma2e_neg_itm + (sum(pos_conc^2))/sigma2e_phr_itm + 1/sigma2b_itm)))


    if (all_dat == F)
    {
      a_itm = stats::rnorm(rep(1, n_patient), ((A2)/sigma2e_phr_itm + mu_a_itm/sigma2a_itm)/
                      ((n_pos)/sigma2e_phr_itm + 1/sigma2a_itm),
                    sqrt(1/((n_pos)/sigma2e_phr_itm + 1/sigma2a_itm)))
      b_itm = stats::rnorm(rep(1, n_patient), ((B2)/sigma2e_phr_itm + mu_b_itm/sigma2b_itm)/
                      ((sum(pos_conc^2))/sigma2e_phr_itm + 1/sigma2b_itm),
                    sqrt(1/((sum(pos_conc^2))/sigma2e_phr_itm + 1/sigma2b_itm)))
    } else
    {
      A1 = colSums(sweep(sweep(neg_dat, 2, b_itm*cc_itm, '-'), 1, d_neg_itm, '-'))
      A3 = colSums(sweep(hk_dat - matrix(rep(b_itm, each = n_hk) * (rep(phi_itm, each = n_hk) + kappa_hk_itm), nrow = n_hk), 1, d_hk_itm, '-'))
      A4 = colSums(sweep(reg_dat - matrix(rep(b_itm, each = n_reg) * (rep(phi_itm, each = n_reg) + kappa_reg_itm), nrow = n_reg), 1, d_reg_itm, '-'))
      a_itm = stats::rnorm(rep(1, n_patient), ((A1/sigma2e_neg_itm) + (A2+A3+A4)/sigma2e_phr_itm + mu_a_itm/sigma2a_itm)/
                      (n_neg/sigma2e_neg_itm + (n_pos+n_hk+n_reg)/sigma2e_phr_itm + 1/sigma2a_itm),
                    sqrt(1/(n_neg/sigma2e_neg_itm + (n_pos+n_hk+n_reg)/sigma2e_phr_itm + 1/sigma2a_itm)))

      B1 = cc_itm * colSums(sweep(sweep(neg_dat, 2, a_itm, '-'), 1, d_neg_itm, '-'))
      B3 = colSums(matrix(rep(phi_itm, each = n_hk) + kappa_hk_itm, nrow = n_hk)*sweep(sweep(hk_dat, 2, a_itm, '-'), 1, d_hk_itm, '-'))
      B4 = colSums(matrix(rep(phi_itm, each = n_reg) + kappa_reg_itm, nrow = n_reg)*sweep(sweep(reg_dat, 2, a_itm, '-'), 1, d_reg_itm, '-'))
      b_itm = stats::rnorm(rep(1, n_patient), ((B1/sigma2e_neg_itm) + (B2+B3+B4)/sigma2e_phr_itm + mu_b_itm/sigma2b_itm)/
                      (n_neg*(cc_itm^2)/sigma2e_neg_itm + (sum(pos_conc^2)+colSums(matrix((rep(phi_itm, each = n_hk) + kappa_hk_itm)^2, nrow = n_hk))+
                                                             colSums(matrix((rep(phi_itm, each = n_reg) + kappa_reg_itm)^2, nrow = n_reg)))/sigma2e_phr_itm + 1/sigma2b_itm),
                    sqrt(1/(n_neg*(cc_itm^2)/sigma2e_neg_itm + (sum(pos_conc^2)+colSums(matrix((rep(phi_itm, each = n_hk) + kappa_hk_itm)^2, nrow = n_hk))+
                                                                  colSums(matrix((rep(phi_itm, each = n_reg) + kappa_reg_itm)^2, nrow = n_reg)))/sigma2e_phr_itm + 1/sigma2b_itm)))
    }

    #print(b_itm)
    cc_itm = truncnorm::rtruncnorm(1, -6, -1, sum(sweep(sweep(sweep(neg_dat, 2, a_itm, '-'), 1, d_neg_itm, '-'), 2, b_itm, '*'))/
                          (n_neg * sum(b_itm^2)), sqrt(sigma2e_neg_itm/(n_neg * sum(b_itm^2))))
    #print(cc_itm)

    phi_itm = numeric(n_patient)
    phi_itm[1:(n_patient - 1)] = truncnorm::rtruncnorm(rep(1, n_patient), phi_L, phi_U,
                                            (colSums(hk_dat-matrix(rep(a_itm, each = n_hk) + rep(b_itm, each = n_hk) * kappa_hk_itm + rep(d_hk_itm, n_patient), nrow = n_hk))+
                                               colSums(reg_dat-matrix(rep(a_itm, each = n_reg) + rep(b_itm, each = n_reg) * kappa_reg_itm + rep(d_reg_itm, n_patient), nrow = n_reg)))/
                                              ((n_hk+n_reg)*b_itm), sqrt(sigma2e_phr_itm/((b_itm^2) * (n_hk+n_reg))))[-n_patient]
    phi_itm[n_patient] = -sum(phi_itm)

    #print(phi_itm)
    kappa_hk_itm = stats::rnorm(rep(1, n_patient * n_hk),
                         ((rep(b_itm, each = n_hk)*(unlist(hk_dat) - rep(a_itm, each = n_hk) - rep(b_itm*phi_itm, each = n_hk) - rep(d_hk_itm, n_patient)))/
                            sigma2e_phr_itm + rep(lambda_hk_itm, n_patient)/sigma2kappa_hk_itm)/rep((b_itm^2)/sigma2e_phr_itm + 1/sigma2kappa_hk_itm, each = n_hk),
                         sqrt(rep(1/((b_itm^2)/sigma2e_phr_itm + 1/sigma2kappa_hk_itm), each = n_hk)))
    #print(kappa_hk_itm)
    kappa_reg_itm = stats::rnorm(rep(1, n_patient * n_reg),
                          ((rep(b_itm, each = n_reg)*(unlist(reg_dat) - rep(a_itm, each = n_reg) - rep(b_itm*phi_itm, each = n_reg) - rep(d_reg_itm, n_patient)))/
                             sigma2e_phr_itm + rep(lambda_reg_itm, n_patient)/sigma2kappa_reg_itm)/rep((b_itm^2)/sigma2e_phr_itm + 1/sigma2kappa_reg_itm, each = n_reg),
                          sqrt(rep(1/((b_itm^2)/sigma2e_phr_itm + 1/sigma2kappa_reg_itm), each = n_reg)))

    lambda_hk_itm = truncnorm::rtruncnorm(rep(1, n_hk), lambda_hk_range[1,], lambda_hk_range[2,],
                               rowMeans(matrix(kappa_hk_itm, nrow = n_hk)), sqrt(sigma2kappa_hk_itm/n_patient))
    lambda_reg_itm = truncnorm::rtruncnorm(rep(1, n_reg), lambda_reg_range[1,], lambda_reg_range[2,],
                                rowMeans(matrix(kappa_reg_itm, nrow = n_reg)), sqrt(sigma2kappa_reg_itm/n_patient))
    #print(kappa_reg_itm)
    d_neg_itm = stats::rnorm(n_neg, rowSums(sweep(neg_dat, 2, a_itm + cc_itm*b_itm, '-')/sigma2e_neg_itm)/(n_patient/sigma2e_neg_itm + 1/sigma2d_neg_itm), sqrt(1/(n_patient/sigma2e_neg_itm + 1/sigma2d_neg_itm)))
    #print(d_neg_itm)
    d_pos_itm = stats::rnorm(n_pos, rowSums(sweep(pos_dat-pos_conc %o% b_itm, 2, a_itm, '-')/sigma2e_phr_itm)/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm), sqrt(1/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm)))
    #print(d_pos_itm)
    d_hk_itm = stats::rnorm(n_hk, rowSums(sweep(hk_dat-matrix(rep(b_itm, each = n_hk)*(rep(phi_itm, each = n_hk)+kappa_hk_itm),nrow = n_hk), 2, a_itm, '-')/sigma2e_phr_itm)/
                       (n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm), sqrt(1/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm)))
    #print(d_hk_itm)
    d_reg_itm = stats::rnorm(n_reg, rowSums(sweep(reg_dat-matrix(rep(b_itm, each = n_reg)*(rep(phi_itm, each = n_reg)+kappa_reg_itm),nrow = n_reg), 2, a_itm, '-')/sigma2e_phr_itm)/
                        (n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm), sqrt(1/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm)))
    #print(d_reg_itm)
    #mu_a_itm = rtruncnorm(1, a_L, a_U, mean(a_itm), sqrt(sigma2a_itm/n_patient))
    mu_a_itm = stats::rnorm(1, (sum(a_itm)/sigma2a_itm + mu_a_mu/sigma2_mu_a)/(n_patient/sigma2a_itm + 1/sigma2_mu_a), sqrt(1/(n_patient/sigma2a_itm + 1/sigma2_mu_a)))
    # #print(mu_a_itm)
    #mu_b_itm = rtruncnorm(1, b_L, b_U, mean(b_itm), sqrt(sigma2b_itm/n_patient))
    mu_b_itm = stats::rnorm(1, (sum(b_itm)/sigma2b_itm + mu_b_mu/sigma2_mu_b)/(n_patient/sigma2b_itm + 1/sigma2_mu_b), sqrt(1/(n_patient/sigma2b_itm + 1/sigma2_mu_b)))
    #print(mu_b_itm)

    sigma2e_neg_itm = 1/stats::rgamma(1, u + n_patient*n_neg/2,  v+sum(sweep(sweep(neg_dat, 2, a_itm + b_itm*cc_itm, '-'), 1, d_neg_itm, '-')^2)/2)
    # sigma2e_phr_itm = 1/stats::rgamma(1, u+n_patient*(n_pos)/2,  v+
    #                               (sum(sweep(sweep(pos_dat - pos_conc %o% b_itm, 2, a_itm, '-'), 1, d_pos_itm, '-')^2))/2)

    sigma2e_phr_itm = 1/stats::rgamma(1, u+n_patient*(n_pos+n_hk+n_reg)/2, v+
                                 (sum(sweep(sweep(pos_dat - pos_conc %o% b_itm, 2, a_itm, '-'), 1, d_pos_itm, '-')^2)+
                                    sum(sweep(hk_dat - matrix(rep(a_itm, each = n_hk) + rep(b_itm, each = n_hk) *
                                                                (rep(phi_itm, each = n_hk) + kappa_hk_itm), nrow = n_hk), 1, d_hk_itm, '-')^2)+
                                    sum(sweep(reg_dat - matrix(rep(a_itm, each = n_reg) + rep(b_itm, each = n_reg) *
                                                                 (rep(phi_itm, each = n_reg) + kappa_reg_itm), nrow = n_reg), 1, d_reg_itm, '-')^2))/2)
    #print(sigma2e_phr_itm)
    sigma2a_itm = 1/stats::rgamma(1, u + n_patient/2,  v + sum((a_itm - mu_a_itm)^2)/2)
    #print(sigma2a_itm)
    sigma2b_itm = 1/stats::rgamma(1, u + n_patient/2,  v + sum((b_itm - mu_b_itm)^2)/2)
    #print(sigma2b_itm)
    sigma2kappa_hk_itm = 1/stats::rgamma(1, u + n_patient*n_hk/2,  v + sum((kappa_hk_itm - rep(lambda_hk_itm, n_patient))^2)/2)
    #print(sigma2kappa_hk_itm)
    sigma2kappa_reg_itm = 1/stats::rgamma(1, u + n_patient*n_reg/2,  v + sum((kappa_reg_itm - rep(lambda_reg_itm, n_patient))^2)/2)
    cat(paste(round(sigma2kappa_reg_itm, 3),''))
    #print(sigma2kappa_reg_itm)
    sigma2d_neg_itm = 1/stats::rgamma(1, u + n_neg/2,  v + sum(d_neg_itm^2)/2)
    #print(sigma2d_neg_itm)
    #sigma2d_phr_itm = 1/stats::rgamma(1, u + n_pos/2,  v + sum(d_pos_itm^2, d_reg_itm^2, d_hk_itm^2)/2)
    sigma2d_phr_itm = 1/stats::rgamma(1, u + n_pos/2,  v + sum(d_pos_itm^2)/2)
    #print(sigma2d_phr_itm)

    if(i > warmup){

      j = i - warmup

      aa[j,] = a_itm
      bb[j,] = b_itm
      cc[j] = cc_itm
      phi_return[j,] = phi_itm
      kappa_hk[j,] = kappa_hk_itm
      kappa_reg[j,] = kappa_reg_itm
      lambda_hk[j,] = lambda_hk_itm
      lambda_reg[j,] = lambda_reg_itm
      d_neg[j,] = d_neg_itm
      d_pos[j,] = d_pos_itm
      d_hk[j,] = d_hk_itm
      d_reg[j,] = d_reg_itm
      mu_a[j] = mu_a_itm
      mu_b[j] = mu_b_itm
      sigma2e_neg[j] = sigma2e_neg_itm
      sigma2e_phr[j] = sigma2e_phr_itm
      sigma2a[j] = sigma2a_itm
      sigma2b[j] = sigma2b_itm
      sigma2kappa_hk[j] = sigma2kappa_hk_itm
      sigma2kappa_reg[j] = sigma2kappa_reg_itm
      sigma2d_neg[j] = sigma2d_neg_itm
      sigma2d_phr[j] = sigma2d_phr_itm
    }
  }

  #get the mcmc samples of key parameters.
  #mcmc.samples = list(mu_a = mu_a, mu_b = mu_b, coef_a = aa, coef_b = bb, patient_eff = phi, kappa_hk = kappa_hk, kappa_reg = kappa_reg,
  #                    sigma2kappa_reg = sigma2kappa_reg, sigma2kappa_hk = sigma2kappa_hk, sigma2e_neg = sigma2e_neg,
  #                    sigma2e_phr = sigma2e_phr, sigma2d_neg = sigma2d_neg, sigma2d_phr = sigma2d_phr)
  mcmc.samples = list(aa = aa,
                      bb = bb,
                      d_pos = d_pos,
                      d_neg = d_neg,
                      cc = cc,
                      mu_a = mu_a,
                      mu_b = mu_b,
                      phi = phi_return,
                      kappa_reg = kappa_reg,
                      d_hk = d_hk,
                      sigma2a = sigma2a,
                      sigma2b = sigma2b,
                      sigma2kappa_reg = sigma2kappa_reg,
                      sigma2kappa_hk = sigma2kappa_hk,
                      sigma2e_neg = sigma2e_neg,
                      sigma2e_phr = sigma2e_phr,
                      sigma2d_neg = sigma2d_neg,
                      sigma2d_phr = sigma2d_phr)


  #calculate summary statistics for the parameters.
  kappa_reg = matrix(colMeans(kappa_reg), nrow = n_reg)
  kappa_hk = matrix(colMeans(kappa_hk), nrow = n_hk)
  cc = mean(cc)
  lambda = c(colMeans(lambda_hk), colMeans(lambda_reg))
  mu_a = mean(mu_a)
  mu_b = mean(mu_b)
  phi = colMeans(phi_return)
  sigma2a = stats::median(sigma2a)
  sigma2b = stats::median(sigma2b)
  sigma2kappa_hk = stats::median(sigma2kappa_hk)
  sigma2kappa_reg = stats::median(sigma2kappa_reg)
  sigma2e_neg = stats::median(sigma2e_neg)
  sigma2e_phr = stats::median(sigma2e_phr)
  sigma2d_neg = stats::median(sigma2d_neg)
  sigma2d_phr = stats::median(sigma2d_phr)

  print(proc.time() - ptm)

  #return(list(mu_a = mu_a, mu_b = mu_b, kappa_hk = kappa_hk, kappa_hk_summary = kappa_hk_summary, kappa_reg = kappa_reg,
  #            kappa_reg_summary = kappa_reg_summary, coef_aa = aa, coef_bb = bb, patient_eff = patient_eff,
  #            sigma2kappa_reg = sigma2kappa_reg, sigma2kappa_hk = sigma2kappa_hk, sigma2e_phr = sigma2e_phr,
  #            sigma2e_neg = sigma2e_neg, sigma2d_neg = sigma2d_neg, sigma2d_phr = sigma2d_phr, mcmc.samples = mcmc.samples))
  #return(list(kappa_reg = kappa_reg, mcmc.samples = mcmc.samples))
  return(list(mu_a = mu_a, mu_b = mu_b, cc = cc, lambda = lambda, phi = phi, kappa_hk = kappa_hk, kappa_reg = kappa_reg,
              sigma2a = sigma2a, sigma2b = sigma2b, sigma2kappa_hk = sigma2kappa_hk, sigma2kappa_reg = sigma2kappa_reg,
              sigma2d_neg = sigma2d_neg, sigma2d_phr = sigma2d_phr, sigma2e_neg = sigma2e_neg, sigma2e_phr = sigma2e_phr))
}


