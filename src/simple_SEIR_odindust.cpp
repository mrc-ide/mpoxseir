// Generated by dust (version 0.15.1) - do not edit
#include <cpp11.hpp>

[[cpp11::register]]
cpp11::sexp dust_simple_SEIR_odindust_gpu_info();
[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                         cpp11::sexp r_n_particles, int n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp gpu_config, cpp11::sexp ode_control);

[[cpp11::register]]
cpp11::sexp dust_cpu_simple_SEIR_odindust_capabilities();

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_run(SEXP ptr, cpp11::sexp r_time_end);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_simulate(SEXP ptr, cpp11::sexp time_end);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_run_adjoint(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state,
                                           SEXP index, SEXP reset_step_size);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_time(SEXP ptr);

[[cpp11::register]]
void dust_cpu_simple_SEIR_odindust_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_rng_state(SEXP ptr, bool first_only, bool last_only);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_set_data(SEXP ptr, cpp11::list data, bool shared);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_compare_data(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood);

[[cpp11::register]]
void dust_cpu_simple_SEIR_odindust_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_cpu_simple_SEIR_odindust_n_state(SEXP ptr);

[[cpp11::register]]
void dust_cpu_simple_SEIR_odindust_set_stochastic_schedule(SEXP ptr, SEXP time);

[[cpp11::register]]
SEXP dust_cpu_simple_SEIR_odindust_ode_statistics(SEXP ptr);
#include <dust/r/dust.hpp>

// Generated by odin.dust (version 0.3.9) - do not edit
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
// [[dust::class(simple_SEIR_odindust)]]
// [[dust::time_type(discrete)]]
// [[dust::param(beta, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta_zoonotic, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(CFR, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(D0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dt, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(E0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(E02, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_E, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_I, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_Id, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_Ir, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Id0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Ir0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(m, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_age, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(S0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
class simple_SEIR_odindust {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    real_type beta;
    real_type beta_zoonotic;
    std::vector<real_type> CFR;
    std::vector<real_type> D0;
    int dim_CFR;
    int dim_D;
    int dim_D0;
    int dim_delta_D;
    int dim_delta_E1;
    int dim_delta_E2;
    int dim_delta_Id;
    int dim_delta_Ir;
    int dim_delta_R;
    int dim_E0;
    int dim_E02;
    int dim_E1;
    int dim_E2;
    int dim_I;
    int dim_Id;
    int dim_Id0;
    int dim_Ir;
    int dim_Ir0;
    int dim_lambda;
    int dim_m;
    int dim_m_1;
    int dim_m_2;
    int dim_n_E1E2;
    int dim_n_E2I;
    int dim_n_E2Id;
    int dim_n_E2Ir;
    int dim_n_IdD;
    int dim_n_IrR;
    int dim_n_SE1;
    int dim_p_SE;
    int dim_R;
    int dim_R0;
    int dim_S;
    int dim_s_ij;
    int dim_s_ij_1;
    int dim_s_ij_2;
    int dim_S0;
    real_type dt;
    std::vector<real_type> E0;
    std::vector<real_type> E02;
    real_type gamma_E;
    real_type gamma_I;
    real_type gamma_Id;
    real_type gamma_Ir;
    std::vector<real_type> Id0;
    std::vector<real_type> initial_D;
    std::vector<real_type> initial_E1;
    std::vector<real_type> initial_E2;
    std::vector<real_type> initial_Id;
    std::vector<real_type> initial_Ir;
    std::vector<real_type> initial_R;
    std::vector<real_type> initial_S;
    real_type initial_time;
    std::vector<real_type> Ir0;
    std::vector<real_type> m;
    int N_age;
    int offset_variable_D;
    int offset_variable_E1;
    int offset_variable_E2;
    int offset_variable_Id;
    int offset_variable_Ir;
    int offset_variable_R;
    real_type p_EE;
    real_type p_EI;
    real_type p_IdD;
    real_type p_IrR;
    std::vector<real_type> R0;
    std::vector<real_type> S0;
  };
  struct internal_type {
    std::vector<real_type> delta_D;
    std::vector<real_type> delta_E1;
    std::vector<real_type> delta_E2;
    std::vector<real_type> delta_Id;
    std::vector<real_type> delta_Ir;
    std::vector<real_type> delta_R;
    std::vector<real_type> I;
    std::vector<real_type> lambda;
    std::vector<real_type> n_E1E2;
    std::vector<real_type> n_E2I;
    std::vector<real_type> n_E2Id;
    std::vector<real_type> n_E2Ir;
    std::vector<real_type> n_IdD;
    std::vector<real_type> n_IrR;
    std::vector<real_type> n_SE1;
    std::vector<real_type> p_SE;
    std::vector<real_type> s_ij;
  };
  simple_SEIR_odindust(const dust::pars_type<simple_SEIR_odindust>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() const {
    return shared->dim_D + shared->dim_E1 + shared->dim_E2 + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 1;
  }
  std::vector<real_type> initial(size_t step, rng_state_type& rng_state) {
    std::vector<real_type> state(shared->dim_D + shared->dim_E1 + shared->dim_E2 + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 1);
    state[0] = shared->initial_time;
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + 1);
    std::copy(shared->initial_E1.begin(), shared->initial_E1.end(), state.begin() + shared->offset_variable_E1);
    std::copy(shared->initial_E2.begin(), shared->initial_E2.end(), state.begin() + shared->offset_variable_E2);
    std::copy(shared->initial_Ir.begin(), shared->initial_Ir.end(), state.begin() + shared->offset_variable_Ir);
    std::copy(shared->initial_Id.begin(), shared->initial_Id.end(), state.begin() + shared->offset_variable_Id);
    std::copy(shared->initial_R.begin(), shared->initial_R.end(), state.begin() + shared->offset_variable_R);
    std::copy(shared->initial_D.begin(), shared->initial_D.end(), state.begin() + shared->offset_variable_D);
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type * S = state + 1;
    const real_type * E1 = state + shared->offset_variable_E1;
    const real_type * E2 = state + shared->offset_variable_E2;
    const real_type * Ir = state + shared->offset_variable_Ir;
    const real_type * Id = state + shared->offset_variable_Id;
    const real_type * R = state + shared->offset_variable_R;
    const real_type * D = state + shared->offset_variable_D;
    state_next[0] = (step + 1) * shared->dt;
    for (int i = 1; i <= shared->dim_I; ++i) {
      internal.I[i - 1] = Ir[i - 1] + Id[i - 1];
    }
    for (int i = 1; i <= shared->dim_n_E1E2; ++i) {
      internal.n_E1E2[i - 1] = dust::random::binomial<real_type>(rng_state, E1[i - 1], shared->p_EE);
    }
    for (int i = 1; i <= shared->dim_n_E2I; ++i) {
      internal.n_E2I[i - 1] = dust::random::binomial<real_type>(rng_state, E2[i - 1], shared->p_EI);
    }
    for (int i = 1; i <= shared->dim_n_IdD; ++i) {
      internal.n_IdD[i - 1] = dust::random::binomial<real_type>(rng_state, Id[i - 1], shared->p_IdD);
    }
    for (int i = 1; i <= shared->dim_n_IrR; ++i) {
      internal.n_IrR[i - 1] = dust::random::binomial<real_type>(rng_state, Ir[i - 1], shared->p_IrR);
    }
    for (int i = 1; i <= shared->dim_delta_D; ++i) {
      internal.delta_D[i - 1] = internal.n_IdD[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_R; ++i) {
      internal.delta_R[i - 1] = internal.n_IrR[i - 1];
    }
    for (int i = 1; i <= shared->dim_n_E2Id; ++i) {
      internal.n_E2Id[i - 1] = dust::random::binomial<real_type>(rng_state, internal.n_E2I[i - 1], shared->CFR[i - 1]);
    }
    for (int i = 1; i <= shared->dim_s_ij_1; ++i) {
      for (int j = 1; j <= shared->dim_s_ij_2; ++j) {
        internal.s_ij[i - 1 + shared->dim_s_ij_1 * (j - 1)] = shared->m[shared->dim_m_1 * (j - 1) + i - 1] * internal.I[j - 1];
      }
    }
    for (int i = 1; i <= shared->dim_delta_Id; ++i) {
      internal.delta_Id[i - 1] = internal.n_E2Id[i - 1] - internal.n_IdD[i - 1];
    }
    for (int i = 1; i <= shared->dim_lambda; ++i) {
      internal.lambda[i - 1] = (shared->beta + shared->beta_zoonotic) * odin_sum2<real_type>(internal.s_ij.data(), i - 1, i, 0, shared->dim_s_ij_2, shared->dim_s_ij_1);
    }
    for (int i = 1; i <= shared->dim_n_E2Ir; ++i) {
      internal.n_E2Ir[i - 1] = internal.n_E2I[i - 1] - internal.n_E2Id[i - 1];
    }
    for (int i = 1; i <= shared->dim_D; ++i) {
      state_next[shared->offset_variable_D + i - 1] = D[i - 1] + internal.delta_D[i - 1];
    }
    for (int i = 1; i <= shared->dim_R; ++i) {
      state_next[shared->offset_variable_R + i - 1] = R[i - 1] + internal.delta_R[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_E2; ++i) {
      internal.delta_E2[i - 1] = internal.n_E1E2[i - 1] - internal.n_E2Ir[i - 1] - internal.n_E2Id[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_Ir; ++i) {
      internal.delta_Ir[i - 1] = internal.n_E2Ir[i - 1] - internal.n_IrR[i - 1];
    }
    for (int i = 1; i <= shared->dim_p_SE; ++i) {
      internal.p_SE[i - 1] = 1 - dust::math::exp(- internal.lambda[i - 1] * shared->dt);
    }
    for (int i = 1; i <= shared->dim_Id; ++i) {
      state_next[shared->offset_variable_Id + i - 1] = Id[i - 1] + internal.delta_Id[i - 1];
    }
    for (int i = 1; i <= shared->dim_n_SE1; ++i) {
      internal.n_SE1[i - 1] = dust::random::binomial<real_type>(rng_state, S[i - 1], internal.p_SE[i - 1]);
    }
    for (int i = 1; i <= shared->dim_E2; ++i) {
      state_next[shared->offset_variable_E2 + i - 1] = E2[i - 1] + internal.delta_E2[i - 1];
    }
    for (int i = 1; i <= shared->dim_Ir; ++i) {
      state_next[shared->offset_variable_Ir + i - 1] = Ir[i - 1] + internal.delta_Ir[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_E1; ++i) {
      internal.delta_E1[i - 1] = internal.n_SE1[i - 1] - internal.n_E1E2[i - 1];
    }
    for (int i = 1; i <= shared->dim_S; ++i) {
      state_next[1 + i - 1] = S[i - 1] - internal.n_SE1[i - 1];
    }
    for (int i = 1; i <= shared->dim_E1; ++i) {
      state_next[shared->offset_variable_E1 + i - 1] = E1[i - 1] + internal.delta_E1[i - 1];
    }
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1) {
  real_type tot = 0.0;
  for (int j = from_j; j < to_j; ++j) {
    int jj = j * dim_x_1;
    for (int i = from_i; i < to_i; ++i) {
      tot += x[i + jj];
    }
  }
  return tot;
}
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_type, typename container>
__host__ __device__
real_type odin_sum1(const container x, size_t from, size_t to) {
  real_type tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
namespace dust {
template<>
dust::pars_type<simple_SEIR_odindust> dust_pars<simple_SEIR_odindust>(cpp11::list user) {
  using real_type = typename simple_SEIR_odindust::real_type;
  auto shared = std::make_shared<simple_SEIR_odindust::shared_type>();
  simple_SEIR_odindust::internal_type internal;
  shared->initial_time = 0;
  shared->beta = NA_REAL;
  shared->beta_zoonotic = NA_REAL;
  shared->dt = NA_REAL;
  shared->gamma_E = NA_REAL;
  shared->gamma_I = NA_REAL;
  shared->gamma_Id = NA_REAL;
  shared->gamma_Ir = NA_REAL;
  shared->N_age = NA_INTEGER;
  shared->beta = user_get_scalar<real_type>(user, "beta", shared->beta, NA_REAL, NA_REAL);
  shared->beta_zoonotic = user_get_scalar<real_type>(user, "beta_zoonotic", shared->beta_zoonotic, NA_REAL, NA_REAL);
  shared->dt = user_get_scalar<real_type>(user, "dt", shared->dt, NA_REAL, NA_REAL);
  shared->gamma_E = user_get_scalar<real_type>(user, "gamma_E", shared->gamma_E, NA_REAL, NA_REAL);
  shared->gamma_I = user_get_scalar<real_type>(user, "gamma_I", shared->gamma_I, NA_REAL, NA_REAL);
  shared->gamma_Id = user_get_scalar<real_type>(user, "gamma_Id", shared->gamma_Id, NA_REAL, NA_REAL);
  shared->gamma_Ir = user_get_scalar<real_type>(user, "gamma_Ir", shared->gamma_Ir, NA_REAL, NA_REAL);
  shared->N_age = user_get_scalar<int>(user, "N_age", shared->N_age, NA_INTEGER, NA_INTEGER);
  shared->dim_CFR = shared->N_age;
  shared->dim_D = shared->N_age;
  shared->dim_D0 = shared->N_age;
  shared->dim_delta_D = shared->N_age;
  shared->dim_delta_E1 = shared->N_age;
  shared->dim_delta_E2 = shared->N_age;
  shared->dim_delta_Id = shared->N_age;
  shared->dim_delta_Ir = shared->N_age;
  shared->dim_delta_R = shared->N_age;
  shared->dim_E0 = shared->N_age;
  shared->dim_E02 = shared->N_age;
  shared->dim_E1 = shared->N_age;
  shared->dim_E2 = shared->N_age;
  shared->dim_I = shared->N_age;
  shared->dim_Id = shared->N_age;
  shared->dim_Id0 = shared->N_age;
  shared->dim_Ir = shared->N_age;
  shared->dim_Ir0 = shared->N_age;
  shared->dim_lambda = shared->N_age;
  shared->dim_m_1 = shared->N_age;
  shared->dim_m_2 = shared->N_age;
  shared->dim_n_E1E2 = shared->N_age;
  shared->dim_n_E2I = shared->N_age;
  shared->dim_n_E2Id = shared->N_age;
  shared->dim_n_E2Ir = shared->N_age;
  shared->dim_n_IdD = shared->N_age;
  shared->dim_n_IrR = shared->N_age;
  shared->dim_n_SE1 = shared->N_age;
  shared->dim_p_SE = shared->N_age;
  shared->dim_R = shared->N_age;
  shared->dim_R0 = shared->N_age;
  shared->dim_S = shared->N_age;
  shared->dim_s_ij_1 = shared->N_age;
  shared->dim_s_ij_2 = shared->N_age;
  shared->dim_S0 = shared->N_age;
  shared->p_EE = 1 - dust::math::exp(- shared->gamma_E * shared->dt);
  shared->p_EI = 1 - dust::math::exp(- shared->gamma_I * shared->dt);
  shared->p_IdD = 1 - dust::math::exp(- shared->gamma_Id * shared->dt);
  shared->p_IrR = 1 - dust::math::exp(- shared->gamma_Ir * shared->dt);
  internal.delta_D = std::vector<real_type>(shared->dim_delta_D);
  internal.delta_E1 = std::vector<real_type>(shared->dim_delta_E1);
  internal.delta_E2 = std::vector<real_type>(shared->dim_delta_E2);
  internal.delta_Id = std::vector<real_type>(shared->dim_delta_Id);
  internal.delta_Ir = std::vector<real_type>(shared->dim_delta_Ir);
  internal.delta_R = std::vector<real_type>(shared->dim_delta_R);
  internal.I = std::vector<real_type>(shared->dim_I);
  shared->initial_D = std::vector<real_type>(shared->dim_D);
  shared->initial_E1 = std::vector<real_type>(shared->dim_E1);
  shared->initial_E2 = std::vector<real_type>(shared->dim_E2);
  shared->initial_Id = std::vector<real_type>(shared->dim_Id);
  shared->initial_Ir = std::vector<real_type>(shared->dim_Ir);
  shared->initial_R = std::vector<real_type>(shared->dim_R);
  shared->initial_S = std::vector<real_type>(shared->dim_S);
  internal.lambda = std::vector<real_type>(shared->dim_lambda);
  internal.n_E1E2 = std::vector<real_type>(shared->dim_n_E1E2);
  internal.n_E2I = std::vector<real_type>(shared->dim_n_E2I);
  internal.n_E2Id = std::vector<real_type>(shared->dim_n_E2Id);
  internal.n_E2Ir = std::vector<real_type>(shared->dim_n_E2Ir);
  internal.n_IdD = std::vector<real_type>(shared->dim_n_IdD);
  internal.n_IrR = std::vector<real_type>(shared->dim_n_IrR);
  internal.n_SE1 = std::vector<real_type>(shared->dim_n_SE1);
  internal.p_SE = std::vector<real_type>(shared->dim_p_SE);
  shared->CFR = user_get_array_fixed<real_type, 1>(user, "CFR", shared->CFR, {shared->dim_CFR}, NA_REAL, NA_REAL);
  shared->D0 = user_get_array_fixed<real_type, 1>(user, "D0", shared->D0, {shared->dim_D0}, NA_REAL, NA_REAL);
  shared->dim_m = shared->dim_m_1 * shared->dim_m_2;
  shared->dim_s_ij = shared->dim_s_ij_1 * shared->dim_s_ij_2;
  shared->E0 = user_get_array_fixed<real_type, 1>(user, "E0", shared->E0, {shared->dim_E0}, NA_REAL, NA_REAL);
  shared->E02 = user_get_array_fixed<real_type, 1>(user, "E02", shared->E02, {shared->dim_E02}, NA_REAL, NA_REAL);
  shared->Id0 = user_get_array_fixed<real_type, 1>(user, "Id0", shared->Id0, {shared->dim_Id0}, NA_REAL, NA_REAL);
  shared->Ir0 = user_get_array_fixed<real_type, 1>(user, "Ir0", shared->Ir0, {shared->dim_Ir0}, NA_REAL, NA_REAL);
  shared->offset_variable_D = shared->dim_E1 + shared->dim_E2 + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 1;
  shared->offset_variable_E1 = shared->dim_S + 1;
  shared->offset_variable_E2 = shared->dim_E1 + shared->dim_S + 1;
  shared->offset_variable_Id = shared->dim_E1 + shared->dim_E2 + shared->dim_Ir + shared->dim_S + 1;
  shared->offset_variable_Ir = shared->dim_E1 + shared->dim_E2 + shared->dim_S + 1;
  shared->offset_variable_R = shared->dim_E1 + shared->dim_E2 + shared->dim_Id + shared->dim_Ir + shared->dim_S + 1;
  shared->R0 = user_get_array_fixed<real_type, 1>(user, "R0", shared->R0, {shared->dim_R0}, NA_REAL, NA_REAL);
  shared->S0 = user_get_array_fixed<real_type, 1>(user, "S0", shared->S0, {shared->dim_S0}, NA_REAL, NA_REAL);
  internal.s_ij = std::vector<real_type>(shared->dim_s_ij);
  for (int i = 1; i <= shared->dim_D; ++i) {
    shared->initial_D[i - 1] = shared->D0[i - 1];
  }
  for (int i = 1; i <= shared->dim_E1; ++i) {
    shared->initial_E1[i - 1] = shared->E0[i - 1];
  }
  for (int i = 1; i <= shared->dim_E2; ++i) {
    shared->initial_E2[i - 1] = shared->E02[i - 1];
  }
  for (int i = 1; i <= shared->dim_Id; ++i) {
    shared->initial_Id[i - 1] = shared->Id0[i - 1];
  }
  for (int i = 1; i <= shared->dim_Ir; ++i) {
    shared->initial_Ir[i - 1] = shared->Ir0[i - 1];
  }
  for (int i = 1; i <= shared->dim_R; ++i) {
    shared->initial_R[i - 1] = shared->R0[i - 1];
  }
  for (int i = 1; i <= shared->dim_S; ++i) {
    shared->initial_S[i - 1] = shared->S0[i - 1];
  }
  shared->m = user_get_array_fixed<real_type, 2>(user, "m", shared->m, {shared->dim_m_1, shared->dim_m_2}, NA_REAL, NA_REAL);
  return dust::pars_type<simple_SEIR_odindust>(shared, internal);
}
template <>
cpp11::sexp dust_info<simple_SEIR_odindust>(const dust::pars_type<simple_SEIR_odindust>& pars) {
  const std::shared_ptr<const simple_SEIR_odindust::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "S", "E1", "E2", "Ir", "Id", "R", "D"});
  cpp11::writable::list dim(8);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({shared->dim_S});
  dim[2] = cpp11::writable::integers({shared->dim_E1});
  dim[3] = cpp11::writable::integers({shared->dim_E2});
  dim[4] = cpp11::writable::integers({shared->dim_Ir});
  dim[5] = cpp11::writable::integers({shared->dim_Id});
  dim[6] = cpp11::writable::integers({shared->dim_R});
  dim[7] = cpp11::writable::integers({shared->dim_D});
  dim.names() = nms;
  cpp11::writable::list index(8);
  index[0] = cpp11::writable::integers({1});
  index[1] = integer_sequence(2, shared->dim_S);
  index[2] = integer_sequence(shared->offset_variable_E1 + 1, shared->dim_E1);
  index[3] = integer_sequence(shared->offset_variable_E2 + 1, shared->dim_E2);
  index[4] = integer_sequence(shared->offset_variable_Ir + 1, shared->dim_Ir);
  index[5] = integer_sequence(shared->offset_variable_Id + 1, shared->dim_Id);
  index[6] = integer_sequence(shared->offset_variable_R + 1, shared->dim_R);
  index[7] = integer_sequence(shared->offset_variable_D + 1, shared->dim_D);
  index.names() = nms;
  size_t len = shared->offset_variable_D + shared->dim_D;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

cpp11::sexp dust_simple_SEIR_odindust_gpu_info() {
  return dust::gpu::r::gpu_info();
}
using model_cpu = dust::dust_cpu<simple_SEIR_odindust>;

cpp11::sexp dust_cpu_simple_SEIR_odindust_capabilities() {
  return dust::r::dust_capabilities<model_cpu>();
}

SEXP dust_cpu_simple_SEIR_odindust_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                             cpp11::sexp r_n_particles, int n_threads,
                             cpp11::sexp r_seed, bool deterministic,
                             cpp11::sexp gpu_config, cpp11::sexp ode_control) {
  return dust::r::dust_cpu_alloc<simple_SEIR_odindust>(r_pars, pars_multi, r_time, r_n_particles,
                                        n_threads, r_seed, deterministic,
                                        gpu_config, ode_control);
}

SEXP dust_cpu_simple_SEIR_odindust_run(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_run<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_simple_SEIR_odindust_simulate(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_simulate<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_simple_SEIR_odindust_run_adjoint(SEXP ptr) {
  return dust::r::dust_run_adjoint<model_cpu>(ptr);
}

SEXP dust_cpu_simple_SEIR_odindust_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust::r::dust_set_index<model_cpu>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_cpu_simple_SEIR_odindust_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size) {
  return dust::r::dust_update_state<model_cpu>(ptr, r_pars, r_state, r_time,
                                                      r_set_initial_state, index, reset_step_size);
}

SEXP dust_cpu_simple_SEIR_odindust_state(SEXP ptr, SEXP r_index) {
  return dust::r::dust_state<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_simple_SEIR_odindust_time(SEXP ptr) {
  return dust::r::dust_time<model_cpu>(ptr);
}

void dust_cpu_simple_SEIR_odindust_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust::r::dust_reorder<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_simple_SEIR_odindust_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust::r::dust_resample<model_cpu>(ptr, r_weights);
}

SEXP dust_cpu_simple_SEIR_odindust_rng_state(SEXP ptr, bool first_only, bool last_only) {
  return dust::r::dust_rng_state<model_cpu>(ptr, first_only, last_only);
}

SEXP dust_cpu_simple_SEIR_odindust_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust::r::dust_set_rng_state<model_cpu>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_cpu_simple_SEIR_odindust_set_data(SEXP ptr, cpp11::list data,
                                       bool shared) {
  dust::r::dust_set_data<model_cpu>(ptr, data, shared);
  return R_NilValue;
}

SEXP dust_cpu_simple_SEIR_odindust_compare_data(SEXP ptr) {
  return dust::r::dust_compare_data<model_cpu>(ptr);
}

SEXP dust_cpu_simple_SEIR_odindust_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood) {
  return dust::r::dust_filter<model_cpu>(ptr, time_end,
                                                save_trajectories,
                                                time_snapshot,
                                                min_log_likelihood);
}

void dust_cpu_simple_SEIR_odindust_set_n_threads(SEXP ptr, int n_threads) {
  return dust::r::dust_set_n_threads<model_cpu>(ptr, n_threads);
}

int dust_cpu_simple_SEIR_odindust_n_state(SEXP ptr) {
  return dust::r::dust_n_state<model_cpu>(ptr);
}

void dust_cpu_simple_SEIR_odindust_set_stochastic_schedule(SEXP ptr, SEXP time) {
  dust::r::dust_set_stochastic_schedule<model_cpu>(ptr, time);
}

SEXP dust_cpu_simple_SEIR_odindust_ode_statistics(SEXP ptr) {
  return dust::r::dust_ode_statistics<model_cpu>(ptr);
}
