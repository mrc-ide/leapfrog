#include <nanobind/nanobind.h>

int return_20() { return 20; }

// namespace nb = nanobind;

NB_MODULE(_core, m) {
  m.doc() = "nanobind example module";

  m.def("return_20", &return_20, R"pbdoc(
      A function that returns 20.
  )pbdoc");
}
