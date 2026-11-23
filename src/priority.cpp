// priority.cpp: lowers the scheduling priority of a kernel worker thread, so a run on all cores leaves the computer responsive:
// other programs get the CPU first and the workers take what is left (on an idle computer they run as fast as before). Outputs do
// not depend on it. Best effort: a call the system refuses changes nothing. Kept out of the kernel files so that windows.h does not
// meet R's headers.
#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <cerrno>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>
#elif defined(__APPLE__) && defined(__has_include)
#if __has_include(<pthread/qos.h>)
#include <pthread.h>
#include <pthread/qos.h>
#define HYPERCOM_QOS 1
#endif
#endif

namespace hc {

const int NICE_STEP = 5;                // Linux: how much nicer than the R process a worker runs (19 is the nicest)
const int QOS_STEP = -5;                // macOS: priority below the default class, which keeps the performance cores (the
                                        // utility class would move the workers to the efficiency cores)

// Called by each worker thread once, at its start; R's own thread keeps its priority
void lower_thread_priority(){
#if defined(_WIN32)
  SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_BELOW_NORMAL);
#elif defined(__linux__)
  const id_t tid = (id_t) syscall(SYS_gettid);
  errno = 0;
  const int nice = getpriority(PRIO_PROCESS, tid);
  if(errno == 0) setpriority(PRIO_PROCESS, tid, nice + NICE_STEP > 19 ? 19 : nice + NICE_STEP);
#elif defined(HYPERCOM_QOS)
  pthread_set_qos_class_self_np(QOS_CLASS_DEFAULT, QOS_STEP);
#endif
}

}  // namespace hc
