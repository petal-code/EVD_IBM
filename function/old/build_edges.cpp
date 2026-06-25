#include <Rcpp.h>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// ==============================================================================
// Global state — initialized once per territory via init_edge_builder()
// ==============================================================================

// Bucket kernel weights: 4 buckets
// 0: d=0.5km (same cell), 1: 1-10km, 2: 10-100km, 3: 100-1000km
static double g_w_bucket[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

// Bucket index matrix: n_cells × n_cells (uint8, row-major)
// Value 0-3 = bucket index, 255 = beyond max radius (unused with full matrix)
static std::vector<uint8_t> g_bucket_mat;
static int g_n_cells = 0;

// cell_age_members[cell_idx * 16 + age_group] = vector of person indices
// Flattened for cache efficiency
static std::vector<std::vector<int>> g_cell_age_members;  // size: n_cells * 16

// Prem multinomial probabilities: q[age_i * 16 + age_j] = M_ij / m_i
static std::vector<double> g_prem_q;   // size: 16 * 16
static std::vector<double> g_prem_mu;  // size: 16 (row sums = m_i)

// NB dispersion parameter
static double g_nb_size = 0.1;

// [[Rcpp::export]]
void init_edge_builder(
    NumericMatrix   bucket_mat_r,     // n_cells × n_cells bucket indices (0-3)
    List            cell_age_members, // list of length n_cells*16, each = int vector
    NumericMatrix   prem_matrix,      // 16 × 16 Prem contact matrix
    NumericVector   bucket_weights,   // length 4: kernel weight per bucket
    double          nb_size) {

  int n_cells = bucket_mat_r.nrow();
  g_n_cells   = n_cells;
  g_nb_size   = nb_size;

  // Store bucket weights
  for (int b = 0; b < 5; ++b) g_w_bucket[b] = bucket_weights[b];

  // Store bucket matrix as uint8 vector (row-major)
  g_bucket_mat.resize((long long)n_cells * n_cells);
  for (int i = 0; i < n_cells; ++i)
    for (int j = 0; j < n_cells; ++j)
      g_bucket_mat[(long long)i * n_cells + j] = (uint8_t)bucket_mat_r(i, j);

  // Store cell-age members
  g_cell_age_members.resize(n_cells * 16);
  for (int idx = 0; idx < n_cells * 16; ++idx) {
    IntegerVector v = cell_age_members[idx];
    g_cell_age_members[idx] = std::vector<int>(v.begin(), v.end());
  }

  // Store Prem matrix as probabilities q[i*16+j] = M_ij / m_i
  g_prem_q.resize(16 * 16);
  g_prem_mu.resize(16);
  for (int i = 0; i < 16; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < 16; ++j) row_sum += prem_matrix(i, j);
    g_prem_mu[i] = row_sum;
    for (int j = 0; j < 16; ++j)
      g_prem_q[i * 16 + j] = (row_sum > 0) ? prem_matrix(i, j) / row_sum : 0.0;
  }
}

// [[Rcpp::export]]
DataFrame build_edges_cpp(
    IntegerVector active_ids,    // 1-based person indices
    IntegerVector cell_ids,      // cell_id per person (1-based)
    IntegerVector hh_ids,        // household id per person
    IntegerVector age_groups,    // age group per person (1-based, 1-16)
    int           seed) {

  std::mt19937 rng(seed);

  std::vector<int> out_from, out_to;
  out_from.reserve(active_ids.size() * 8);
  out_to.reserve(active_ids.size() * 8);

  int n_active = active_ids.size();

  // Pre-compute cumulative Prem probabilities for Multinomial sampling
  // cum_q[i][j] = sum_{k<=j} q[i*16+k]
  std::vector<std::vector<double>> cum_q(16, std::vector<double>(16));
  for (int i = 0; i < 16; ++i) {
    cum_q[i][0] = g_prem_q[i * 16];
    for (int j = 1; j < 16; ++j)
      cum_q[i][j] = cum_q[i][j-1] + g_prem_q[i * 16 + j];
  }

  for (int k = 0; k < n_active; ++k) {
    if (k % 10000 == 0)
      Rcpp::Rcout << "    Layer 2: " << k << " / " << n_active << " persons\n" << std::flush;

    int id_k     = active_ids[k];           // 1-based
    int c_k      = cell_ids[id_k - 1] - 1; // 0-based cell index
    int hh_k     = hh_ids[id_k - 1];
    int age_i    = age_groups[id_k - 1] - 1; // 0-based age group (0-15)

    // ── Step 1: Draw total contacts D_k ~ NB(m_i, theta) ──────────────────
    double mu_i = g_prem_mu[age_i];
    if (mu_i <= 0.0) continue;

    // NB parameterization: size=theta, mu=mu_i
    // Use gamma-Poisson mixture: X ~ Poisson(lambda), lambda ~ Gamma(size, size/mu)
    std::gamma_distribution<double> gamma_dist(g_nb_size, mu_i / g_nb_size);
    double lambda = gamma_dist(rng);
    std::poisson_distribution<int> pois_dist(lambda);
    int D_k = pois_dist(rng);
    if (D_k == 0) continue;

    // ── Step 2: Multinomial(D_k, q_i) → d_j per age group ────────────────
    int d_j[16] = {0};
    int remaining = D_k;
    for (int j = 0; j < 15 && remaining > 0; ++j) {
      // Draw from Binomial(remaining, p_j_given_remaining)
      double p_remaining = (1.0 - (j > 0 ? cum_q[age_i][j-1] : 0.0));
      if (p_remaining <= 0.0) break;
      double p_j = g_prem_q[age_i * 16 + j] / p_remaining;
      p_j = std::min(std::max(p_j, 0.0), 1.0);
      std::binomial_distribution<int> binom(remaining, p_j);
      d_j[j] = binom(rng);
      remaining -= d_j[j];
    }
    d_j[15] = remaining;  // Last group gets remainder

    // ── Step 3: For each age group j, sample d_j contacts ─────────────────
    for (int j = 0; j < 16; ++j) {
      if (d_j[j] == 0) continue;

      // Build weighted cell candidate list.
      // Weight = w_bucket[b] * (n_age_j[c] / bucket_total_pop[b])
      // This ensures inter-bucket allocation follows kernel proportions,
      // while within each bucket, cells are selected proportional to population.

      // First pass: compute per-bucket total population for normalization
      double bucket_pop[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (int c = 0; c < g_n_cells; ++c) {
        int n_members = (int)g_cell_age_members[c * 16 + j].size();
        if (n_members == 0) continue;
        uint8_t bucket = g_bucket_mat[(long long)c_k * g_n_cells + c];
        bucket_pop[bucket] += n_members;
      }

      // Second pass: normalized weights
      std::vector<int>    cand_cells;
      std::vector<double> cand_weights;
      cand_cells.reserve(g_n_cells);
      cand_weights.reserve(g_n_cells);

      double w_total = 0.0;
      for (int c = 0; c < g_n_cells; ++c) {
        int n_members = (int)g_cell_age_members[c * 16 + j].size();
        if (n_members == 0) continue;
        uint8_t bucket = g_bucket_mat[(long long)c_k * g_n_cells + c];
        if (bucket_pop[bucket] <= 0.0) continue;
        // w_bucket controls inter-bucket ratio; n_members/bucket_pop within-bucket ratio
        double w = g_w_bucket[bucket] * (double)n_members / bucket_pop[bucket];
        if (w <= 0.0) continue;
        cand_cells.push_back(c);
        cand_weights.push_back(w);
        w_total += w;
      }

      if (cand_cells.empty() || w_total <= 0.0) continue;

      // Sample d_j contacts from weighted cells (with replacement)
      std::uniform_real_distribution<double> udist(0.0, w_total);

      // Pre-compute cumulative weights
      std::vector<double> cum_w(cand_cells.size());
      cum_w[0] = cand_weights[0];
      for (int ci = 1; ci < (int)cand_cells.size(); ++ci)
        cum_w[ci] = cum_w[ci-1] + cand_weights[ci];

      for (int s = 0; s < d_j[j]; ++s) {
        // Sample a cell proportional to kernel × population
        double u = udist(rng);
        int chosen_c = (int)(std::lower_bound(cum_w.begin(), cum_w.end(), u)
                               - cum_w.begin());
        chosen_c = std::min(chosen_c, (int)cand_cells.size() - 1);
        int c_j = cand_cells[chosen_c];

        // Sample a person uniformly from age_group j in cell c_j
        const std::vector<int>& members = g_cell_age_members[c_j * 16 + j];
        std::uniform_int_distribution<int> uid(0, (int)members.size() - 1);
        int cand_id = members[uid(rng)];  // 1-based person id

        // Exclude self and same household
        if (cand_id == id_k) continue;
        if (hh_ids[cand_id - 1] == hh_k) continue;

        // Store as undirected edge (i < j convention)
        if (cand_id > id_k) {
          out_from.push_back(id_k);
          out_to.push_back(cand_id);
        } else {
          out_from.push_back(cand_id);
          out_to.push_back(id_k);
        }
      }
    }
  }

  return DataFrame::create(
    Named("from") = IntegerVector(out_from.begin(), out_from.end()),
    Named("to")   = IntegerVector(out_to.begin(),   out_to.end())
  );
}
