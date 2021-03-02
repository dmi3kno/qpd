functions {
    vector metalog_coefficients(data vector p, data vector quantiles) {
    int n = rows(p);
    matrix[n, n] Y;
    if (n == 0) reject("p cannot be of size zero");
    if (rows(quantiles) == n) {
      vector[n] log_odds = n > 1 ? logit(p) : rep_vector(not_a_number(), n);
      vector[n] pmhalf =   n > 2 ? p - 0.5  : rep_vector(not_a_number(), n);
      int odd = 1;
      Y[ , 1] = rep_vector(1, n);
      if (n > 1) Y[ , 2] = log_odds;
      if (n > 2) Y[ , 3] = pmhalf .* log_odds;
      if (n > 3) Y[ , 4] = pmhalf;
      for (m in 5:n) {
        if (odd) {
          pmhalf .*= pmhalf;
          Y[ , m]  = pmhalf;
        } else Y[ , m] = pmhalf .* log_odds;
        odd = odd == 0;
      }
    } else reject("p and quantiles must be of the same size");
    return Y \ quantiles;
  }

  vector metalog_lpdfm(data vector p, data vector a) {
    int n = rows(a);
    vector[n] log_odds = n > 1 ? logit(p) : rep_vector(not_a_number(), n);
    vector[n] p_bin =  n > 1 ? p .* (1 - p): rep_vector(not_a_number(), n);
    vector[n] pmhalf = n > 2 ? p - 0.5 : rep_vector(not_a_number(), n);
    int odd = 1;
    vector[n] mm1_half = rep_vector(1, n);
    vector[n] m_halfm1 = rep_vector(1, n);
    vector[n] pmhalf_o = rep_vector(1, n);
    vector[n] pmhalf_e = rep_vector(1, n);
    vector[n] res;

    if (n>1) res  = a[2] ./ p_bin;
    if (n>2) res = res + a[3] * (p_bin ./ pmhalf + log_odds);
    if (n>3) res = res + a[4];
    for (m in 5:n) {
      if (odd) {
        pmhalf_o .*= pmhalf;
        mm1_half += rep_vector(1, n);
        res += a[m] * mm1_half .* pmhalf_o;
      } else {
        pmhalf_e .*= pmhalf;
        m_halfm1 += rep_vector(1, n);
        res += a[m] * (pmhalf_e .* pmhalf ./ p_bin + m_halfm1 .* pmhalf_e .* log_odds);
        }
      odd = odd == 0;
    }

    return 1 ./ res;
  }

}
