#include <Rcpp.h>
#include <random>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <numeric>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// ==============================================================================
// Global state — initialized once per territory via init_edge_builder()
//
// Design (v2 — physical / non-physical drawn independently):
//   - Non-physical degree per participant: D_np ~ Poisson(mu_nonphys_i)
//   - Physical    degree per participant: D_p  ~ NegBin(mu_phys_i, size=phys_nb_size)
//   - Each degree is split across target age groups via Multinomial using
//     its OWN row-normalized probability vector (q_nonphys / q_phys), since
//     physical and non-physical contacts can have different age patterns.
//   - Stratum (daily/weekly/monthly) is still allocated via
//     Multinomial(1, [p_daily, p_weekly, p_monthly]) per participant age,
//     independent of physical status.
//   - Collision handling: for a given participant k, the PHYSICAL pass is
//     sampled first; any partner already claimed as physical is recorded in
//     a per-person used-target set. The NON-PHYSICAL pass then skips any
//     candidate already in that set — i.e. "physical wins" on collision.
//     This guarantees zero physical/non-physical duplicate edges within a
//     single participant's own draws.
// ==============================================================================

static double g_w_bucket[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

static std::vector<uint8_t> g_bucket_mat;
static int g_n_cells = 0;

// cell_age_members[cell_idx * 16 + age_group] = vector of person indices
static std::vector<std::vector<int>> g_cell_age_members;

// Non-physical unique-contact matrix -> row-normalized probs + row sums (mean degree)
static std::vector<double> g_nonphys_q;   // [age_i*16 + age_j]
static std::vector<double> g_nonphys_mu;  // [age_i]

// Physical unique-contact matrix -> row-normalized probs + row sums (mean degree)
static std::vector<double> g_phys_q;      // [age_i*16 + age_j]
static std::vector<double> g_phys_mu;     // [age_i]

// Stratum allocation probs: g_stratum_p[age_i * 3 + s]
// s=0: daily, s=1: weekly, s=2: monthly
static std::vector<double> g_stratum_p;

static double g_phys_nb_size = 0.25;

// [[Rcpp::export]]
void init_edge_builder(
    NumericMatrix   bucket_mat_r,       // n_cells x n_cells bucket indices (0-4)
    List            cell_age_members,   // list length n_cells*16, each = int vector
    NumericMatrix   nonphys_unique,     // 16x16 unique NON-PHYSICAL contact matrix
    NumericMatrix   phys_unique,        // 16x16 unique PHYSICAL contact matrix
    NumericMatrix   stratum_prob_mat,   // 16x3 stratum allocation probs [p_daily, p_weekly, p_monthly]
    NumericVector   bucket_weights,     // length 5: kernel weight per bucket
    double          phys_nb_size) {     // NB dispersion for physical draws (e.g. 0.25)

  int n_cells    = bucket_mat_r.nrow();
  g_n_cells      = n_cells;
  g_phys_nb_size = phys_nb_size;

  for (int b = 0; b < 5; ++b) g_w_bucket[b] = bucket_weights[b];

  // Bucket matrix as uint8 (row-major)
  g_bucket_mat.resize((long long)n_cells * n_cells);
  for (int i = 0; i < n_cells; ++i)
    for (int j = 0; j < n_cells; ++j)
      g_bucket_mat[(long long)i * n_cells + j] = (uint8_t)bucket_mat_r(i, j);

  // Cell-age member index
  g_cell_age_members.resize(n_cells * 16);
  for (int idx = 0; idx < n_cells * 16; ++idx) {
    IntegerVector v = cell_age_members[idx];
    g_cell_age_members[idx] = std::vector<int>(v.begin(), v.end());
  }

  // Helper to row-normalize a 16x16 matrix into q (probs) + mu (row sums)
  auto build_q_mu = [&](NumericMatrix mat,
                        std::vector<double>& q_out,
                        std::vector<double>& mu_out) {
    q_out.resize(16 * 16);
    mu_out.resize(16);
    for (int i = 0; i < 16; ++i) {
      double row_sum = 0.0;
      for (int j = 0; j < 16; ++j) row_sum += mat(i, j);
      mu_out[i] = row_sum;
      for (int j = 0; j < 16; ++j)
        q_out[i * 16 + j] = (row_sum > 0) ? mat(i, j) / row_sum : 0.0;
    }
  };

  build_q_mu(nonphys_unique, g_nonphys_q, g_nonphys_mu);
  build_q_mu(phys_unique,    g_phys_q,    g_phys_mu);

  // Stratum allocation probabilities: 16 x 3 (flat, row-major)
  g_stratum_p.resize(16 * 3);
  for (int i = 0; i < 16; ++i)
    for (int s = 0; s < 3; ++s)
      g_stratum_p[i * 3 + s] = stratum_prob_mat(i, s);
}

// ==============================================================================
// [Helper] Multinomial split of a total draw D across 16 age-group buckets,
// using a row-normalized cumulative probability vector cum_q[age_i][*]
// ==============================================================================
static void multinomial_split(
    int D, int age_i,
    const std::vector<std::vector<double>>& cum_q,
    const std::vector<double>& q_flat,
    std::mt19937& rng,
    int d_j_out[16]) {

  for (int j = 0; j < 16; ++j) d_j_out[j] = 0;
  if (D <= 0) return;

  int remaining = D;
  for (int j = 0; j < 15 && remaining > 0; ++j) {
    double p_remaining = 1.0 - (j > 0 ? cum_q[age_i][j-1] : 0.0);
    if (p_remaining <= 0.0) break;
    double p_j = q_flat[age_i * 16 + j] / p_remaining;
    p_j = std::min(std::max(p_j, 0.0), 1.0);
    std::binomial_distribution<int> binom(remaining, p_j);
    d_j_out[j] = binom(rng);
    remaining -= d_j_out[j];
  }
  d_j_out[15] = remaining;
}

// [[Rcpp::export]]
DataFrame build_edges_cpp(
    IntegerVector active_ids,
    IntegerVector cell_ids,
    IntegerVector hh_ids,
    IntegerVector age_groups,
    int           seed) {

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif01(0.0, 1.0);

  std::vector<int> out_from, out_to;
  std::vector<int> out_stratum;     // 0=daily, 1=weekly, 2=monthly
  std::vector<int> out_is_physical; // 0/1

  out_from.reserve(active_ids.size() * 8);
  out_to.reserve(active_ids.size() * 8);
  out_stratum.reserve(active_ids.size() * 8);
  out_is_physical.reserve(active_ids.size() * 8);

  int n_active = active_ids.size();

  // Pre-compute cumulative probabilities for Multinomial age-partner sampling
  std::vector<std::vector<double>> cum_q_np(16, std::vector<double>(16));
  std::vector<std::vector<double>> cum_q_p (16, std::vector<double>(16));
  for (int i = 0; i < 16; ++i) {
    cum_q_np[i][0] = g_nonphys_q[i * 16];
    cum_q_p[i][0]  = g_phys_q[i * 16];
    for (int j = 1; j < 16; ++j) {
      cum_q_np[i][j] = cum_q_np[i][j-1] + g_nonphys_q[i * 16 + j];
      cum_q_p[i][j]  = cum_q_p[i][j-1]  + g_phys_q[i * 16 + j];
    }
  }

  // Pre-compute cumulative stratum probs per age group
  std::vector<std::vector<double>> cum_s(16, std::vector<double>(3));
  for (int i = 0; i < 16; ++i) {
    cum_s[i][0] = g_stratum_p[i * 3];
    cum_s[i][1] = cum_s[i][0] + g_stratum_p[i * 3 + 1];
    cum_s[i][2] = 1.0;  // ensures last bucket catches remainder
  }

  for (int k = 0; k < n_active; ++k) {
    if (k % 10000 == 0)
      Rcpp::Rcout << "    edges: " << k << " / " << n_active << "\n" << std::flush;

    int id_k  = active_ids[k];
    int c_k   = cell_ids[id_k - 1] - 1;    // 0-based cell index
    int hh_k  = hh_ids[id_k - 1];
    int age_i = age_groups[id_k - 1] - 1;  // 0-based (0-15)

    // ── Step 1: independent degree draws ──────────────────────────────────
    // Non-physical: Poisson(mu_nonphys_i)
    double mu_np = g_nonphys_mu[age_i];
    int D_np = 0;
    if (mu_np > 0.0) {
      std::poisson_distribution<int> pois_np(mu_np);
      D_np = pois_np(rng);
    }

    // Physical: NegBin(mu_phys_i, size = phys_nb_size), via Gamma-Poisson mixture
    double mu_p = g_phys_mu[age_i];
    int D_p = 0;
    if (mu_p > 0.0) {
      std::gamma_distribution<double> gamma_dist(g_phys_nb_size, mu_p / g_phys_nb_size);
      double lambda = gamma_dist(rng);
      std::poisson_distribution<int> pois_p(lambda);
      D_p = pois_p(rng);
    }

    if (D_np == 0 && D_p == 0) continue;

    // ── Step 2: Multinomial split across target age groups (own q for each) ──
    int d_j_np[16], d_j_p[16];
    multinomial_split(D_np, age_i, cum_q_np, g_nonphys_q, rng, d_j_np);
    multinomial_split(D_p,  age_i, cum_q_p,  g_phys_q,    rng, d_j_p);

    // Per-participant used-target set: physical pass fills it first, then
    // non-physical pass skips anything already present ("physical wins").
    std::unordered_set<int> used_targets;
    used_targets.reserve(32);

    // ── Step 3: for each target age group j, sample partners (physical first) ──
    for (int j = 0; j < 16; ++j) {
      if (d_j_np[j] == 0 && d_j_p[j] == 0) continue;

      // Build kernel-weighted cell candidate list once, shared by both passes
      double bucket_pop[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (int c = 0; c < g_n_cells; ++c) {
        int n_members = (int)g_cell_age_members[c * 16 + j].size();
        if (n_members == 0) continue;
        uint8_t bucket = g_bucket_mat[(long long)c_k * g_n_cells + c];
        bucket_pop[bucket] += n_members;
      }

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
        double w = g_w_bucket[bucket] * (double)n_members / bucket_pop[bucket];
        if (w <= 0.0) continue;
        cand_cells.push_back(c);
        cand_weights.push_back(w);
        w_total += w;
      }

      if (cand_cells.empty() || w_total <= 0.0) continue;

      std::vector<double> cum_w(cand_cells.size());
      cum_w[0] = cand_weights[0];
      for (int ci = 1; ci < (int)cand_cells.size(); ++ci)
        cum_w[ci] = cum_w[ci-1] + cand_weights[ci];

      std::uniform_real_distribution<double> udist(0.0, w_total);

      auto sample_one_partner = [&]() -> int {
        double u = udist(rng);
        int chosen_c = (int)(std::lower_bound(cum_w.begin(), cum_w.end(), u)
                               - cum_w.begin());
        chosen_c = std::min(chosen_c, (int)cand_cells.size() - 1);
        int c_j = cand_cells[chosen_c];

        const std::vector<int>& members = g_cell_age_members[c_j * 16 + j];
        std::uniform_int_distribution<int> uid(0, (int)members.size() - 1);
        return members[uid(rng)];
      };

      // ── Physical pass first ──────────────────────────────────────────────
      for (int s = 0; s < d_j_p[j]; ++s) {
        int cand_id = sample_one_partner();
        if (cand_id == id_k) continue;
        if (hh_ids[cand_id - 1] == hh_k) continue;
        if (used_targets.count(cand_id)) continue;  // dedup within physical pass
        used_targets.insert(cand_id);

        double u_s   = unif01(rng);
        int stratum  = 2;
        if      (u_s < cum_s[age_i][0]) stratum = 0;
        else if (u_s < cum_s[age_i][1]) stratum = 1;

        if (cand_id > id_k) { out_from.push_back(id_k);     out_to.push_back(cand_id); }
        else                { out_from.push_back(cand_id);  out_to.push_back(id_k);    }
        out_stratum.push_back(stratum);
        out_is_physical.push_back(1);
      }

      // ── Non-physical pass — skip anything already claimed as physical ────
      for (int s = 0; s < d_j_np[j]; ++s) {
        int cand_id = sample_one_partner();
        if (cand_id == id_k) continue;
        if (hh_ids[cand_id - 1] == hh_k) continue;
        if (used_targets.count(cand_id)) continue;  // physical wins on collision
        used_targets.insert(cand_id);

        double u_s   = unif01(rng);
        int stratum  = 2;
        if      (u_s < cum_s[age_i][0]) stratum = 0;
        else if (u_s < cum_s[age_i][1]) stratum = 1;

        if (cand_id > id_k) { out_from.push_back(id_k);     out_to.push_back(cand_id); }
        else                { out_from.push_back(cand_id);  out_to.push_back(id_k);    }
        out_stratum.push_back(stratum);
        out_is_physical.push_back(0);
      }
    }
  }

  return DataFrame::create(
    Named("from")        = IntegerVector(out_from.begin(),        out_from.end()),
    Named("to")          = IntegerVector(out_to.begin(),          out_to.end()),
    Named("stratum")     = IntegerVector(out_stratum.begin(),     out_stratum.end()),
    Named("is_physical") = IntegerVector(out_is_physical.begin(), out_is_physical.end())
  );
}
