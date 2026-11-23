// memory_os.cpp: the computer's physical memory, in a file of its own so that windows.h never meets R's headers (both define ERROR
// and other names)
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif

// Bytes of physical memory, or -1 when the system does not tell
double physical_memory_bytes(){
#ifdef _WIN32
  MEMORYSTATUSEX m;
  m.dwLength = sizeof m;
  if(!GlobalMemoryStatusEx(&m)) return -1;
  return (double) m.ullTotalPhys;
#else
  const long pages = sysconf(_SC_PHYS_PAGES), size = sysconf(_SC_PAGESIZE);
  if(pages <= 0 || size <= 0) return -1;
  return (double) pages * (double) size;
#endif
}
