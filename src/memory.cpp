// memory.cpp: physical memory for prioritize_lr's memory budget (memory_os.cpp asks the system)
#include <Rcpp.h>

double physical_memory_bytes();

// Bytes of physical memory, NA when the system does not tell
// [[Rcpp::export]]
double physical_memory_cpp(){
  const double bytes = physical_memory_bytes();

  return bytes > 0 ? bytes : NA_REAL;
}
