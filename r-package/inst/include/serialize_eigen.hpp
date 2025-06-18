#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <unsupported/Eigen/CXX11/Tensor>

namespace serialize {

namespace internal {

template<typename T>
const char *data_name();

template<>
inline const char *data_name<int>() {
  return "int";
}

template<>
inline const char *data_name<double>() {
  return "double";
}

const std::string WHITESPACE = " \n\r\t\f\v";

std::string ltrim(const std::string &s) {
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string &s) {
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string &s) {
  return rtrim(ltrim(s));
}

template<typename T>
struct csv_contents {
  std::vector<size_t> dim;
  std::vector<T> data;
};

template<typename T>
std::vector<T> strsplit(const std::string &s) {
  std::vector<T> ret;
  std::stringstream ss(s);

  while (ss.good()) {
    std::string substr;
    getline(ss, substr, ',');
    ret.push_back(static_cast<T>(std::stod(substr)));
  }

  return ret;
}

template<typename T>
csv_contents<T> parse_csv(const std::string &path) {
  std::ifstream src(path);
  if (!src.good()) {
    std::stringstream ss;
    ss << "File at path '" << path << "' does not exist." << std::endl;
    throw std::runtime_error(ss.str());
  }
  std::string line;
  getline(src, line);
  line = trim(line);
  if (line != data_name<T>()) {
    std::stringstream ss;
    ss << "Data at path '" << path << "' is of wrong type. Trying to read as: '" << data_name<T>()
       << "', data saved as: '" << line << "'." << std::endl;
    throw std::runtime_error(ss.str());
  }
  std::string line2;
  getline(src, line2);
  line2 = trim(line2);
  auto dim = strsplit<size_t>(line2);

  std::string line3;
  getline(src, line3);
  line2 = trim(line3);
  std::vector<T> data = strsplit<T>(line3);
  // TODO: validate data.size() is prod(dim)

  return csv_contents<T>{dim, data};
}

}

template<typename T, size_t rank>
Eigen::Tensor<T, rank> deserialize_tensor(const std::string &path) {
  const auto contents = internal::parse_csv<T>(path);
  if (contents.dim.size() != rank) {
    std::stringstream ss;
    ss << "Data at path '" << path << "' is of wrong rank. Trying to read with rank: '" << rank
       << "', data saved as rank: '" << contents.dim.size() << "'." << std::endl;
    throw std::runtime_error(ss.str());
  }

  Eigen::array <Eigen::Index, rank> dim;
  for (size_t i = 0; i < rank; ++i) {
    dim[i] = static_cast<Eigen::Index>(contents.dim[i]);
  }

  Eigen::Tensor <T, rank> ret(dim);
  for (size_t i = 0; i < contents.data.size(); ++i) {
    ret(i) = contents.data[i];
  }

  return ret;
}

template<typename T>
T deserialize_scalar(const std::string &path) {
  const auto contents = internal::parse_csv<T>(path);
  if (contents.data.size() != 1) {
    std::stringstream ss;
    ss << "Data at path '" << path << "' is of wrong size. Trying to read as scalar '"
       << "but data has size: '" << contents.data.size() << "'." << std::endl;
    throw std::runtime_error(ss.str());
  }

  return contents.data[0];
}

// Basic seralization of tensor
// Writes out a file with 3 lines
// First line is the type, only double or int supported
// 2nd line is the dimensions comma separated e.g. 2,3,2
// 3rd line is the comma separated data
template<typename T, int rank>
void serialize_tensor(const Eigen::Tensor <T, rank> &data, const std::string &path) {
  std::ofstream dest(path);
  dest << internal::data_name<T>() << std::endl;
  const auto d = data.dimensions();
  for (int i = 0; i < rank; ++i) {
    dest << d[i];
    if (i == rank - 1) {
      dest << std::endl;
    } else {
      dest << ",";
    }
  }
  for (int i = 0; i < data.size(); ++i) {
    dest << data.data()[i];
    if (i == data.size() - 1) {
      dest << std::endl;
    } else {
      dest << ",";
    }
  }
  dest.close();
}

template<typename T, int rank>
void serialize_tensor_map(const Eigen::TensorMap <Eigen::Tensor<T, rank>> &data, const std::string &path) {
  Eigen::Tensor <T, rank> d = data;
  serialize_tensor(d, path);
}

}
