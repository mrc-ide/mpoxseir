// Generated by odin.dust (version 0.3.9) - do not edit
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
template <typename real_type, typename T, typename U>
__host__ __device__ real_type fmodr(T x, U y) {
  real_type tmp = std::fmod(static_cast<real_type>(x),
                            static_cast<real_type>(y));
  if (tmp * y < 0) {
    tmp += y;
  }
  return tmp;
}

// These exist to support the model on the gpu, as in C++14 std::min
// and std::max are constexpr and error without --expt-relaxed-constexpr
template <typename T>
__host__ __device__ T odin_min(T x, T y) {
  return x < y ? x : y;
}

template <typename T>
__host__ __device__ T odin_max(T x, T y) {
  return x > y ? x : y;
}

template <typename T>
__host__ __device__ T odin_sign(T x) {
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}
// [[dust::class(model)]]
// [[dust::time_type(discrete)]]
// [[dust::param(beta_h, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(beta_z, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(CFR, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(D0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Ea0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Eb0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_E, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_I, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_Id, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(gamma_Ir, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Id0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Ir0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(m, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(n_group, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(S0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dt, has_default = TRUE, default_value = 1L, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class model {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    real_type beta_h;
    std::vector<real_type> beta_z;
    std::vector<real_type> CFR;
    std::vector<real_type> D0;
    int dim_beta_z;
    int dim_CFR;
    int dim_D;
    int dim_D0;
    int dim_delta_D;
    int dim_delta_Ea;
    int dim_delta_Eb;
    int dim_delta_Id;
    int dim_delta_Ir;
    int dim_delta_R;
    int dim_E;
    int dim_Ea;
    int dim_Ea0;
    int dim_Eb;
    int dim_Eb0;
    int dim_I;
    int dim_Id;
    int dim_Id0;
    int dim_Ir;
    int dim_Ir0;
    int dim_lambda;
    int dim_m;
    int dim_m_1;
    int dim_m_2;
    int dim_N;
    int dim_n_EaEb;
    int dim_n_EbI;
    int dim_n_EbId;
    int dim_n_EbIr;
    int dim_n_IdD;
    int dim_n_IrR;
    int dim_n_SEa;
    int dim_p_SE;
    int dim_R;
    int dim_R0;
    int dim_S;
    int dim_s_ij;
    int dim_s_ij_1;
    int dim_s_ij_2;
    int dim_S0;
    real_type dt;
    std::vector<real_type> Ea0;
    std::vector<real_type> Eb0;
    real_type gamma_E;
    real_type gamma_I;
    real_type gamma_Id;
    real_type gamma_Ir;
    std::vector<real_type> Id0;
    std::vector<real_type> initial_D;
    real_type initial_D_tot;
    std::vector<real_type> initial_E;
    real_type initial_E_tot;
    std::vector<real_type> initial_Ea;
    std::vector<real_type> initial_Eb;
    std::vector<real_type> initial_I;
    real_type initial_I_tot;
    std::vector<real_type> initial_Id;
    std::vector<real_type> initial_Ir;
    std::vector<real_type> initial_N;
    real_type initial_N_tot;
    std::vector<real_type> initial_R;
    real_type initial_R_tot;
    std::vector<real_type> initial_S;
    real_type initial_S_tot;
    real_type initial_weekly_cases;
    real_type initial_weekly_cases_00_04;
    real_type initial_weekly_cases_05_14;
    real_type initial_weekly_cases_15;
    real_type initial_weekly_deaths;
    real_type initial_weekly_deaths_00_04;
    real_type initial_weekly_deaths_05_14;
    real_type initial_weekly_deaths_15;
    std::vector<real_type> Ir0;
    std::vector<real_type> m;
    int n_group;
    int offset_variable_D;
    int offset_variable_E;
    int offset_variable_Ea;
    int offset_variable_Eb;
    int offset_variable_I;
    int offset_variable_Id;
    int offset_variable_Ir;
    int offset_variable_N;
    int offset_variable_R;
    real_type p_EE;
    real_type p_EI;
    real_type p_IdD;
    real_type p_IrR;
    std::vector<real_type> R0;
    std::vector<real_type> S0;
    real_type steps_per_week;
  };
  struct internal_type {
    std::vector<real_type> delta_D;
    std::vector<real_type> delta_Ea;
    std::vector<real_type> delta_Eb;
    std::vector<real_type> delta_Id;
    std::vector<real_type> delta_Ir;
    std::vector<real_type> delta_R;
    real_type initial_time;
    std::vector<real_type> lambda;
    std::vector<real_type> n_EaEb;
    std::vector<real_type> n_EbI;
    std::vector<real_type> n_EbId;
    std::vector<real_type> n_EbIr;
    std::vector<real_type> n_IdD;
    std::vector<real_type> n_IrR;
    std::vector<real_type> n_SEa;
    std::vector<real_type> p_SE;
    std::vector<real_type> s_ij;
  };
  model(const dust::pars_type<model>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() const {
    return shared->dim_D + shared->dim_E + shared->dim_Ea + shared->dim_Eb + shared->dim_I + shared->dim_Id + shared->dim_Ir + shared->dim_N + shared->dim_R + shared->dim_S + 15;
  }
  std::vector<real_type> initial(size_t step, rng_state_type& rng_state) {
    std::vector<real_type> state(shared->dim_D + shared->dim_E + shared->dim_Ea + shared->dim_Eb + shared->dim_I + shared->dim_Id + shared->dim_Ir + shared->dim_N + shared->dim_R + shared->dim_S + 15);
    internal.initial_time = step;
    state[0] = internal.initial_time;
    state[1] = shared->initial_S_tot;
    state[2] = shared->initial_E_tot;
    state[3] = shared->initial_I_tot;
    state[4] = shared->initial_R_tot;
    state[5] = shared->initial_D_tot;
    state[6] = shared->initial_N_tot;
    state[7] = shared->initial_weekly_cases;
    state[8] = shared->initial_weekly_cases_00_04;
    state[9] = shared->initial_weekly_cases_05_14;
    state[10] = shared->initial_weekly_cases_15;
    state[11] = shared->initial_weekly_deaths;
    state[12] = shared->initial_weekly_deaths_00_04;
    state[13] = shared->initial_weekly_deaths_05_14;
    state[14] = shared->initial_weekly_deaths_15;
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + 15);
    std::copy(shared->initial_Ea.begin(), shared->initial_Ea.end(), state.begin() + shared->offset_variable_Ea);
    std::copy(shared->initial_Eb.begin(), shared->initial_Eb.end(), state.begin() + shared->offset_variable_Eb);
    std::copy(shared->initial_Ir.begin(), shared->initial_Ir.end(), state.begin() + shared->offset_variable_Ir);
    std::copy(shared->initial_Id.begin(), shared->initial_Id.end(), state.begin() + shared->offset_variable_Id);
    std::copy(shared->initial_R.begin(), shared->initial_R.end(), state.begin() + shared->offset_variable_R);
    std::copy(shared->initial_D.begin(), shared->initial_D.end(), state.begin() + shared->offset_variable_D);
    std::copy(shared->initial_E.begin(), shared->initial_E.end(), state.begin() + shared->offset_variable_E);
    std::copy(shared->initial_I.begin(), shared->initial_I.end(), state.begin() + shared->offset_variable_I);
    std::copy(shared->initial_N.begin(), shared->initial_N.end(), state.begin() + shared->offset_variable_N);
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type * S = state + 15;
    const real_type * Ea = state + shared->offset_variable_Ea;
    const real_type * Eb = state + shared->offset_variable_Eb;
    const real_type * Ir = state + shared->offset_variable_Ir;
    const real_type * Id = state + shared->offset_variable_Id;
    const real_type * R = state + shared->offset_variable_R;
    const real_type * D = state + shared->offset_variable_D;
    const real_type * E = state + shared->offset_variable_E;
    const real_type * I = state + shared->offset_variable_I;
    const real_type * N = state + shared->offset_variable_N;
    const real_type weekly_cases = state[7];
    const real_type weekly_cases_00_04 = state[8];
    const real_type weekly_cases_05_14 = state[9];
    const real_type weekly_deaths = state[11];
    const real_type weekly_deaths_00_04 = state[12];
    const real_type weekly_deaths_05_14 = state[13];
    state_next[0] = (step + 1) * shared->dt;
    real_type is_same_week = fmodr<real_type>(step, shared->steps_per_week) > 0;
    state_next[5] = odin_sum1<real_type>(D, 0, shared->dim_D);
    for (int i = 1; i <= shared->dim_E; ++i) {
      state_next[shared->offset_variable_E + i - 1] = Ea[i - 1] + Eb[i - 1];
    }
    state_next[2] = odin_sum1<real_type>(E, 0, shared->dim_E);
    for (int i = 1; i <= shared->dim_I; ++i) {
      state_next[shared->offset_variable_I + i - 1] = Ir[i - 1] + Id[i - 1];
    }
    state_next[3] = odin_sum1<real_type>(I, 0, shared->dim_I);
    for (int i = 1; i <= shared->dim_N; ++i) {
      state_next[shared->offset_variable_N + i - 1] = S[i - 1] + Ea[i - 1] + Eb[i - 1] + Ir[i - 1] + Id[i - 1] + R[i - 1] + D[i - 1];
    }
    state_next[6] = odin_sum1<real_type>(N, 0, shared->dim_N);
    state_next[4] = odin_sum1<real_type>(R, 0, shared->dim_R);
    state_next[1] = odin_sum1<real_type>(S, 0, shared->dim_S);
    for (int i = 1; i <= shared->dim_n_EaEb; ++i) {
      internal.n_EaEb[i - 1] = dust::random::binomial<real_type>(rng_state, Ea[i - 1], shared->p_EE);
    }
    for (int i = 1; i <= shared->dim_n_EbI; ++i) {
      internal.n_EbI[i - 1] = dust::random::binomial<real_type>(rng_state, Eb[i - 1], shared->p_EI);
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
    for (int i = 1; i <= shared->dim_delta_Eb; ++i) {
      internal.delta_Eb[i - 1] = internal.n_EaEb[i - 1] - internal.n_EbI[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_R; ++i) {
      internal.delta_R[i - 1] = internal.n_IrR[i - 1];
    }
    for (int i = 1; i <= shared->dim_n_EbId; ++i) {
      internal.n_EbId[i - 1] = dust::random::binomial<real_type>(rng_state, internal.n_EbI[i - 1], shared->CFR[i - 1]);
    }
    for (int i = 1; i <= shared->dim_s_ij_1; ++i) {
      for (int j = 1; j <= shared->dim_s_ij_2; ++j) {
        internal.s_ij[i - 1 + shared->dim_s_ij_1 * (j - 1)] = shared->m[shared->dim_m_1 * (j - 1) + i - 1] * I[j - 1];
      }
    }
    state_next[11] = weekly_deaths * is_same_week + odin_sum1<real_type>(internal.n_IdD.data(), 0, shared->dim_n_IdD);
    state_next[12] = weekly_deaths_00_04 * is_same_week + internal.n_IdD[0];
    state_next[13] = weekly_deaths_05_14 * is_same_week + odin_sum1<real_type>(internal.n_IdD.data(), 1, 3);
    state_next[14] = weekly_deaths * is_same_week + odin_sum1<real_type>(internal.n_IdD.data(), 3, shared->n_group);
    for (int i = 1; i <= shared->dim_delta_Id; ++i) {
      internal.delta_Id[i - 1] = internal.n_EbId[i - 1] - internal.n_IdD[i - 1];
    }
    for (int i = 1; i <= shared->dim_lambda; ++i) {
      internal.lambda[i - 1] = shared->beta_h * odin_sum2<real_type>(internal.s_ij.data(), i - 1, i, 0, shared->dim_s_ij_2, shared->dim_s_ij_1) + shared->beta_z[i - 1];
    }
    for (int i = 1; i <= shared->dim_n_EbIr; ++i) {
      internal.n_EbIr[i - 1] = internal.n_EbI[i - 1] - internal.n_EbId[i - 1];
    }
    for (int i = 1; i <= shared->dim_D; ++i) {
      state_next[shared->offset_variable_D + i - 1] = D[i - 1] + internal.delta_D[i - 1];
    }
    for (int i = 1; i <= shared->dim_Eb; ++i) {
      state_next[shared->offset_variable_Eb + i - 1] = Eb[i - 1] + internal.delta_Eb[i - 1];
    }
    for (int i = 1; i <= shared->dim_R; ++i) {
      state_next[shared->offset_variable_R + i - 1] = R[i - 1] + internal.delta_R[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_Ir; ++i) {
      internal.delta_Ir[i - 1] = internal.n_EbIr[i - 1] - internal.n_IrR[i - 1];
    }
    for (int i = 1; i <= shared->dim_p_SE; ++i) {
      internal.p_SE[i - 1] = 1 - dust::math::exp(- internal.lambda[i - 1] * shared->dt);
    }
    for (int i = 1; i <= shared->dim_Id; ++i) {
      state_next[shared->offset_variable_Id + i - 1] = Id[i - 1] + internal.delta_Id[i - 1];
    }
    for (int i = 1; i <= shared->dim_n_SEa; ++i) {
      internal.n_SEa[i - 1] = dust::random::binomial<real_type>(rng_state, S[i - 1], internal.p_SE[i - 1]);
    }
    for (int i = 1; i <= shared->dim_Ir; ++i) {
      state_next[shared->offset_variable_Ir + i - 1] = Ir[i - 1] + internal.delta_Ir[i - 1];
    }
    for (int i = 1; i <= shared->dim_delta_Ea; ++i) {
      internal.delta_Ea[i - 1] = internal.n_SEa[i - 1] - internal.n_EaEb[i - 1];
    }
    for (int i = 1; i <= shared->dim_S; ++i) {
      state_next[15 + i - 1] = S[i - 1] - internal.n_SEa[i - 1];
    }
    state_next[7] = weekly_cases * is_same_week + odin_sum1<real_type>(internal.n_SEa.data(), 0, shared->dim_n_SEa);
    state_next[8] = weekly_cases_00_04 * is_same_week + internal.n_SEa[0];
    state_next[9] = weekly_cases_05_14 * is_same_week + odin_sum1<real_type>(internal.n_SEa.data(), 1, 3);
    state_next[10] = weekly_cases * is_same_week + odin_sum1<real_type>(internal.n_SEa.data(), 3, shared->n_group);
    for (int i = 1; i <= shared->dim_Ea; ++i) {
      state_next[shared->offset_variable_Ea + i - 1] = Ea[i - 1] + internal.delta_Ea[i - 1];
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
dust::pars_type<model> dust_pars<model>(cpp11::list user) {
  using real_type = typename model::real_type;
  auto shared = std::make_shared<model::shared_type>();
  model::internal_type internal;
  shared->initial_weekly_cases = 0;
  shared->initial_weekly_cases_00_04 = 0;
  shared->initial_weekly_cases_05_14 = 0;
  shared->initial_weekly_cases_15 = 0;
  shared->initial_weekly_deaths = 0;
  shared->initial_weekly_deaths_00_04 = 0;
  shared->initial_weekly_deaths_05_14 = 0;
  shared->initial_weekly_deaths_15 = 0;
  shared->beta_h = NA_REAL;
  shared->gamma_E = NA_REAL;
  shared->gamma_I = NA_REAL;
  shared->gamma_Id = NA_REAL;
  shared->gamma_Ir = NA_REAL;
  shared->n_group = NA_INTEGER;
  shared->dt = 1;
  internal.initial_time = 0;
  shared->beta_h = user_get_scalar<real_type>(user, "beta_h", shared->beta_h, NA_REAL, NA_REAL);
  shared->dt = user_get_scalar<real_type>(user, "dt", shared->dt, NA_REAL, NA_REAL);
  shared->gamma_E = user_get_scalar<real_type>(user, "gamma_E", shared->gamma_E, NA_REAL, NA_REAL);
  shared->gamma_I = user_get_scalar<real_type>(user, "gamma_I", shared->gamma_I, NA_REAL, NA_REAL);
  shared->gamma_Id = user_get_scalar<real_type>(user, "gamma_Id", shared->gamma_Id, NA_REAL, NA_REAL);
  shared->gamma_Ir = user_get_scalar<real_type>(user, "gamma_Ir", shared->gamma_Ir, NA_REAL, NA_REAL);
  shared->n_group = user_get_scalar<int>(user, "n_group", shared->n_group, NA_INTEGER, NA_INTEGER);
  shared->dim_beta_z = shared->n_group;
  shared->dim_CFR = shared->n_group;
  shared->dim_D = shared->n_group;
  shared->dim_D0 = shared->n_group;
  shared->dim_delta_D = shared->n_group;
  shared->dim_delta_Ea = shared->n_group;
  shared->dim_delta_Eb = shared->n_group;
  shared->dim_delta_Id = shared->n_group;
  shared->dim_delta_Ir = shared->n_group;
  shared->dim_delta_R = shared->n_group;
  shared->dim_E = shared->n_group;
  shared->dim_Ea = shared->n_group;
  shared->dim_Ea0 = shared->n_group;
  shared->dim_Eb = shared->n_group;
  shared->dim_Eb0 = shared->n_group;
  shared->dim_I = shared->n_group;
  shared->dim_Id = shared->n_group;
  shared->dim_Id0 = shared->n_group;
  shared->dim_Ir = shared->n_group;
  shared->dim_Ir0 = shared->n_group;
  shared->dim_lambda = shared->n_group;
  shared->dim_m_1 = shared->n_group;
  shared->dim_m_2 = shared->n_group;
  shared->dim_N = shared->n_group;
  shared->dim_n_EaEb = shared->n_group;
  shared->dim_n_EbI = shared->n_group;
  shared->dim_n_EbId = shared->n_group;
  shared->dim_n_EbIr = shared->n_group;
  shared->dim_n_IdD = shared->n_group;
  shared->dim_n_IrR = shared->n_group;
  shared->dim_n_SEa = shared->n_group;
  shared->dim_p_SE = shared->n_group;
  shared->dim_R = shared->n_group;
  shared->dim_R0 = shared->n_group;
  shared->dim_S = shared->n_group;
  shared->dim_s_ij_1 = shared->n_group;
  shared->dim_s_ij_2 = shared->n_group;
  shared->dim_S0 = shared->n_group;
  shared->p_EE = 1 - dust::math::exp(- shared->gamma_E * shared->dt);
  shared->p_EI = 1 - dust::math::exp(- shared->gamma_I * shared->dt);
  shared->p_IdD = 1 - dust::math::exp(- shared->gamma_Id * shared->dt);
  shared->p_IrR = 1 - dust::math::exp(- shared->gamma_Ir * shared->dt);
  shared->steps_per_week = 7 / (real_type) shared->dt;
  internal.delta_D = std::vector<real_type>(shared->dim_delta_D);
  internal.delta_Ea = std::vector<real_type>(shared->dim_delta_Ea);
  internal.delta_Eb = std::vector<real_type>(shared->dim_delta_Eb);
  internal.delta_Id = std::vector<real_type>(shared->dim_delta_Id);
  internal.delta_Ir = std::vector<real_type>(shared->dim_delta_Ir);
  internal.delta_R = std::vector<real_type>(shared->dim_delta_R);
  shared->initial_D = std::vector<real_type>(shared->dim_D);
  shared->initial_E = std::vector<real_type>(shared->dim_E);
  shared->initial_Ea = std::vector<real_type>(shared->dim_Ea);
  shared->initial_Eb = std::vector<real_type>(shared->dim_Eb);
  shared->initial_I = std::vector<real_type>(shared->dim_I);
  shared->initial_Id = std::vector<real_type>(shared->dim_Id);
  shared->initial_Ir = std::vector<real_type>(shared->dim_Ir);
  shared->initial_N = std::vector<real_type>(shared->dim_N);
  shared->initial_R = std::vector<real_type>(shared->dim_R);
  shared->initial_S = std::vector<real_type>(shared->dim_S);
  internal.lambda = std::vector<real_type>(shared->dim_lambda);
  internal.n_EaEb = std::vector<real_type>(shared->dim_n_EaEb);
  internal.n_EbI = std::vector<real_type>(shared->dim_n_EbI);
  internal.n_EbId = std::vector<real_type>(shared->dim_n_EbId);
  internal.n_EbIr = std::vector<real_type>(shared->dim_n_EbIr);
  internal.n_IdD = std::vector<real_type>(shared->dim_n_IdD);
  internal.n_IrR = std::vector<real_type>(shared->dim_n_IrR);
  internal.n_SEa = std::vector<real_type>(shared->dim_n_SEa);
  internal.p_SE = std::vector<real_type>(shared->dim_p_SE);
  shared->beta_z = user_get_array_fixed<real_type, 1>(user, "beta_z", shared->beta_z, {shared->dim_beta_z}, NA_REAL, NA_REAL);
  shared->CFR = user_get_array_fixed<real_type, 1>(user, "CFR", shared->CFR, {shared->dim_CFR}, NA_REAL, NA_REAL);
  shared->D0 = user_get_array_fixed<real_type, 1>(user, "D0", shared->D0, {shared->dim_D0}, NA_REAL, NA_REAL);
  shared->dim_m = shared->dim_m_1 * shared->dim_m_2;
  shared->dim_s_ij = shared->dim_s_ij_1 * shared->dim_s_ij_2;
  shared->Ea0 = user_get_array_fixed<real_type, 1>(user, "Ea0", shared->Ea0, {shared->dim_Ea0}, NA_REAL, NA_REAL);
  shared->Eb0 = user_get_array_fixed<real_type, 1>(user, "Eb0", shared->Eb0, {shared->dim_Eb0}, NA_REAL, NA_REAL);
  shared->Id0 = user_get_array_fixed<real_type, 1>(user, "Id0", shared->Id0, {shared->dim_Id0}, NA_REAL, NA_REAL);
  shared->Ir0 = user_get_array_fixed<real_type, 1>(user, "Ir0", shared->Ir0, {shared->dim_Ir0}, NA_REAL, NA_REAL);
  shared->offset_variable_D = shared->dim_Ea + shared->dim_Eb + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 15;
  shared->offset_variable_E = shared->dim_D + shared->dim_Ea + shared->dim_Eb + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 15;
  shared->offset_variable_Ea = shared->dim_S + 15;
  shared->offset_variable_Eb = shared->dim_Ea + shared->dim_S + 15;
  shared->offset_variable_I = shared->dim_D + shared->dim_E + shared->dim_Ea + shared->dim_Eb + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 15;
  shared->offset_variable_Id = shared->dim_Ea + shared->dim_Eb + shared->dim_Ir + shared->dim_S + 15;
  shared->offset_variable_Ir = shared->dim_Ea + shared->dim_Eb + shared->dim_S + 15;
  shared->offset_variable_N = shared->dim_D + shared->dim_E + shared->dim_Ea + shared->dim_Eb + shared->dim_I + shared->dim_Id + shared->dim_Ir + shared->dim_R + shared->dim_S + 15;
  shared->offset_variable_R = shared->dim_Ea + shared->dim_Eb + shared->dim_Id + shared->dim_Ir + shared->dim_S + 15;
  shared->R0 = user_get_array_fixed<real_type, 1>(user, "R0", shared->R0, {shared->dim_R0}, NA_REAL, NA_REAL);
  shared->S0 = user_get_array_fixed<real_type, 1>(user, "S0", shared->S0, {shared->dim_S0}, NA_REAL, NA_REAL);
  internal.s_ij = std::vector<real_type>(shared->dim_s_ij);
  for (int i = 1; i <= shared->dim_D; ++i) {
    shared->initial_D[i - 1] = shared->D0[i - 1];
  }
  shared->initial_D_tot = odin_sum1<real_type>(shared->D0.data(), 0, shared->dim_D0);
  for (int i = 1; i <= shared->dim_E; ++i) {
    shared->initial_E[i - 1] = shared->Ea0[i - 1] + shared->Eb0[i - 1];
  }
  shared->initial_E_tot = odin_sum1<real_type>(shared->Ea0.data(), 0, shared->dim_Ea0) + odin_sum1<real_type>(shared->Eb0.data(), 0, shared->dim_Eb0);
  for (int i = 1; i <= shared->dim_Ea; ++i) {
    shared->initial_Ea[i - 1] = shared->Ea0[i - 1];
  }
  for (int i = 1; i <= shared->dim_Eb; ++i) {
    shared->initial_Eb[i - 1] = shared->Eb0[i - 1];
  }
  for (int i = 1; i <= shared->dim_I; ++i) {
    shared->initial_I[i - 1] = shared->Ir0[i - 1] + shared->Id0[i - 1];
  }
  shared->initial_I_tot = odin_sum1<real_type>(shared->Ir0.data(), 0, shared->dim_Ir0) + odin_sum1<real_type>(shared->Id0.data(), 0, shared->dim_Id0);
  for (int i = 1; i <= shared->dim_Id; ++i) {
    shared->initial_Id[i - 1] = shared->Id0[i - 1];
  }
  for (int i = 1; i <= shared->dim_Ir; ++i) {
    shared->initial_Ir[i - 1] = shared->Ir0[i - 1];
  }
  for (int i = 1; i <= shared->dim_N; ++i) {
    shared->initial_N[i - 1] = shared->S0[i - 1] + shared->Ea0[i - 1] + shared->Eb0[i - 1] + shared->Ir0[i - 1] + shared->Id0[i - 1] + shared->R0[i - 1] + shared->D0[i - 1];
  }
  shared->initial_N_tot = odin_sum1<real_type>(shared->S0.data(), 0, shared->dim_S0) + odin_sum1<real_type>(shared->Ea0.data(), 0, shared->dim_Ea0) + odin_sum1<real_type>(shared->Eb0.data(), 0, shared->dim_Eb0) + odin_sum1<real_type>(shared->Ir0.data(), 0, shared->dim_Ir0) + odin_sum1<real_type>(shared->Id0.data(), 0, shared->dim_Id0) + odin_sum1<real_type>(shared->R0.data(), 0, shared->dim_R0) + odin_sum1<real_type>(shared->D0.data(), 0, shared->dim_D0);
  for (int i = 1; i <= shared->dim_R; ++i) {
    shared->initial_R[i - 1] = shared->R0[i - 1];
  }
  shared->initial_R_tot = odin_sum1<real_type>(shared->R0.data(), 0, shared->dim_R0);
  for (int i = 1; i <= shared->dim_S; ++i) {
    shared->initial_S[i - 1] = shared->S0[i - 1];
  }
  shared->initial_S_tot = odin_sum1<real_type>(shared->S0.data(), 0, shared->dim_S0);
  shared->m = user_get_array_fixed<real_type, 2>(user, "m", shared->m, {shared->dim_m_1, shared->dim_m_2}, NA_REAL, NA_REAL);
  return dust::pars_type<model>(shared, internal);
}
template <>
cpp11::sexp dust_info<model>(const dust::pars_type<model>& pars) {
  const std::shared_ptr<const model::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "S_tot", "E_tot", "I_tot", "R_tot", "D_tot", "N_tot", "weekly_cases", "weekly_cases_00_04", "weekly_cases_05_14", "weekly_cases_15", "weekly_deaths", "weekly_deaths_00_04", "weekly_deaths_05_14", "weekly_deaths_15", "S", "Ea", "Eb", "Ir", "Id", "R", "D", "E", "I", "N"});
  cpp11::writable::list dim(25);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({1});
  dim[6] = cpp11::writable::integers({1});
  dim[7] = cpp11::writable::integers({1});
  dim[8] = cpp11::writable::integers({1});
  dim[9] = cpp11::writable::integers({1});
  dim[10] = cpp11::writable::integers({1});
  dim[11] = cpp11::writable::integers({1});
  dim[12] = cpp11::writable::integers({1});
  dim[13] = cpp11::writable::integers({1});
  dim[14] = cpp11::writable::integers({1});
  dim[15] = cpp11::writable::integers({shared->dim_S});
  dim[16] = cpp11::writable::integers({shared->dim_Ea});
  dim[17] = cpp11::writable::integers({shared->dim_Eb});
  dim[18] = cpp11::writable::integers({shared->dim_Ir});
  dim[19] = cpp11::writable::integers({shared->dim_Id});
  dim[20] = cpp11::writable::integers({shared->dim_R});
  dim[21] = cpp11::writable::integers({shared->dim_D});
  dim[22] = cpp11::writable::integers({shared->dim_E});
  dim[23] = cpp11::writable::integers({shared->dim_I});
  dim[24] = cpp11::writable::integers({shared->dim_N});
  dim.names() = nms;
  cpp11::writable::list index(25);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = cpp11::writable::integers({6});
  index[6] = cpp11::writable::integers({7});
  index[7] = cpp11::writable::integers({8});
  index[8] = cpp11::writable::integers({9});
  index[9] = cpp11::writable::integers({10});
  index[10] = cpp11::writable::integers({11});
  index[11] = cpp11::writable::integers({12});
  index[12] = cpp11::writable::integers({13});
  index[13] = cpp11::writable::integers({14});
  index[14] = cpp11::writable::integers({15});
  index[15] = integer_sequence(16, shared->dim_S);
  index[16] = integer_sequence(shared->offset_variable_Ea + 1, shared->dim_Ea);
  index[17] = integer_sequence(shared->offset_variable_Eb + 1, shared->dim_Eb);
  index[18] = integer_sequence(shared->offset_variable_Ir + 1, shared->dim_Ir);
  index[19] = integer_sequence(shared->offset_variable_Id + 1, shared->dim_Id);
  index[20] = integer_sequence(shared->offset_variable_R + 1, shared->dim_R);
  index[21] = integer_sequence(shared->offset_variable_D + 1, shared->dim_D);
  index[22] = integer_sequence(shared->offset_variable_E + 1, shared->dim_E);
  index[23] = integer_sequence(shared->offset_variable_I + 1, shared->dim_I);
  index[24] = integer_sequence(shared->offset_variable_N + 1, shared->dim_N);
  index.names() = nms;
  size_t len = shared->offset_variable_N + shared->dim_N;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}
