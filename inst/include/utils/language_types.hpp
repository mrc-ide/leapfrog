#pragma once

#ifdef CALLER_R
#include <Rcpp.h>

using InputData = Rcpp::List;

struct OutputData {
  Rcpp::List data;
  Rcpp::CharacterVector names;
};

#else

#include <filesystem>

using InputData = std::filesystem::path;

using OutputData = std::filesystem::path;

#endif
