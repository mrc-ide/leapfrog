#pragma once

#ifdef CALLER_R
#include <Rcpp.h>

using InputData = Rcpp::List;

struct OutputData {

  Rcpp::List data;
  Rcpp::CharacterVector names;

};

OutputData initialize_output(int output_size) {
  Rcpp::List data(output_size);
  Rcpp::CharacterVector names(output_size);
  data.attr("names") = names;
  return OutputData{
    data,
    names
  };
}

#else

#include <filesystem>

using InputData = std::filesystem::path;

#endif
