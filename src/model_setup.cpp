#include <Rcpp.h>

template <typename T>
T* r_data(SEXP x) {
  static_assert(sizeof(T) == 0, "Only specializations of r_data can be used");
}

template <>
double* r_data(SEXP x) {
  return REAL(x);
}

template <>
int * r_data(SEXP x) {
  return INTEGER(x);
}

template<typename T, typename... Args>
auto parse_data(const Rcpp::List data, const std::string& key, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::array<int, rank> dimensions{ static_cast<int>(dims)... };

  int length = std::accumulate(dimensions.begin(), dimensions.end(), 1, std::multiplies<int>());
  SEXP array_data = data[key];
  // In cases where the input data has project years we might not use all of it model fit
  // So we can take create a Map over a smaller slice of the data
  // As long as this is true we can be confident we're not referencing invalid memory
  if (LENGTH(array_data) < length) {
    Rcpp::stop("Invalid size of data for '%s', expected %d got %d",
               key,
               length,
               LENGTH(array_data));
  }

  return Eigen::TensorMap<Eigen::Tensor<T, rank>>(r_data<T>(array_data), static_cast<int>(dims)...);
}

template<typename real_type, leapfrog::HivAgeStratification S>
leapfrog::Parameters <real_type> setup_model(const Rcpp::List data,
                                             const leapfrog::Options<double, S> opts,
                                             const int proj_years) {
  constexpr auto ss = StateSpace<S>();
  const leapfrog::TensorMap2<double> base_pop = parse_data<double>(data, "basepop",
                                                                   ss.age_groups_pop, ss.num_genders);
  const leapfrog::TensorMap3<double> survival = parse_data<double>(data, "Sx",
                                                                   ss.age_groups_pop + 1, ss.num_genders, proj_years);
  const leapfrog::TensorMap3<double> net_migration = parse_data<double>(data, "netmigr_adj",
                                                                        ss.age_groups_pop, ss.num_genders, proj_years);
  const leapfrog::TensorMap2<double> age_sex_fertility_ratio = parse_data<double>(data, "asfr",
                                                                                  ss.age_groups_fert, proj_years);
  const leapfrog::TensorMap2<double> births_sex_prop = parse_data<double>(data, "births_sex_prop",
                                                                          ss.num_genders, proj_years);
  const leapfrog::TensorMap1<double> incidence_rate = parse_data<double>(data, "incidinput",
                                                                         proj_years);
  const leapfrog::TensorMap3<double> incidence_relative_risk_age = parse_data<double>(
      data, "incrr_age", ss.age_groups_pop - ss.hiv_adult_first_age_group, ss.num_genders, proj_years);
  const leapfrog::TensorMap1<double> incidence_relative_risk_sex = parse_data<double>(data,
                                                                                      "incrr_sex", proj_years);
  const leapfrog::TensorMap3<double> cd4_mortality = parse_data<double>(data, "cd4_mort",
                                                                        ss.disease_stages, ss.age_groups_hiv,
                                                                        ss.num_genders);
  const leapfrog::TensorMap3<double> cd4_progression = parse_data<double>(data, "cd4_prog",
                                                                          ss.disease_stages - 1, ss.age_groups_hiv,
                                                                          ss.num_genders);

  leapfrog::Tensor1<int> artcd4elig_idx = parse_data<int>(data, "artcd4elig_idx", proj_years + 1);
  for (int i = 0; i <= proj_years; ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    artcd4elig_idx(i) = artcd4elig_idx(i) - 1;
  }

  const leapfrog::TensorMap3<double> cd4_initdist = parse_data<double>(data, "cd4_initdist",
                                                                       ss.disease_stages, ss.age_groups_hiv,
                                                                       ss.num_genders);
  const leapfrog::TensorMap4<double> art_mortality = parse_data<double>(data, "art_mort",
                                                                        ss.treatment_stages, ss.disease_stages,
                                                                        ss.age_groups_hiv, ss.num_genders);
  const leapfrog::TensorMap2<double> artmx_timerr = parse_data<double>(data, "artmx_timerr",
                                                                       ss.treatment_stages, proj_years);
  leapfrog::Tensor1<double> h_art_stage_dur(ss.treatment_stages - 1);
  h_art_stage_dur.setConstant(0.5);

  const leapfrog::TensorMap1<double> art_dropout = parse_data<double>(data, "art_dropout", proj_years);
  const leapfrog::TensorMap2<double> art15plus_num = parse_data<double>(data, "art15plus_num",
                                                                        ss.num_genders, proj_years);
  const leapfrog::TensorMap2<int> art15plus_isperc = parse_data<int>(data, "art15plus_isperc",
                                                                     ss.num_genders, proj_years);

  const leapfrog::Parameters<double> params = {incidence_rate,
                                               base_pop,
                                               survival,
                                               net_migration,
                                               age_sex_fertility_ratio,
                                               births_sex_prop,
                                               incidence_relative_risk_age,
                                               incidence_relative_risk_sex,
                                               cd4_mortality,
                                               cd4_progression,
                                               artcd4elig_idx,
                                               cd4_initdist,
                                               art_mortality,
                                               artmx_timerr,
                                               h_art_stage_dur,
                                               art_dropout,
                                               art15plus_num,
                                               art15plus_isperc};
  return params;
}
