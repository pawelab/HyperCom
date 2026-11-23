// threads.h: the worker threads of the kernels (score_grids.cpp, screen_grids.cpp, influence.cpp)
//
// The calling thread does no kernel work: it starts the workers, then waits for them and asks R every POLL_MS whether the user has
// interrupted (Ctrl-C, Esc), which only R's main thread may ask. A worker that throws, a thread that cannot start, or an interrupt
// sets Stop::failed, which the workers check at their barriers and between work units, and return. After the join the calling thread
// raises the first error as an R error, or the interrupt as R's own interrupt (Rcpp's export wrapper signals it again), so nothing
// a worker does can end the R session.
// Workers run at a lower scheduling priority than other programs (priority.cpp), so a run on all cores leaves the computer responsive.
#ifndef HYPERCOM_THREADS_H
#define HYPERCOM_THREADS_H

#include <Rcpp.h>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <mutex>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace hc {

const int POLL_MS = 100;                // how often the waiting thread asks R for an interrupt

// Lowers the calling worker thread's scheduling priority, so other programs stay responsive (priority.cpp)
void lower_thread_priority();

// Whether the user has interrupted R; main thread only. R_CheckUserInterrupt longjmps on an interrupt, so Rcpp runs it inside
// R_ToplevelExec, which turns the jump into an exception here
inline bool user_interrupt(){
  try {
    Rcpp::checkUserInterrupt();
  } catch(Rcpp::internal::InterruptedException &){
    return true;
  }
  return false;
}

// user_interrupt at most every POLL_MS, for work that runs on the calling thread between its units
class Poll {
 public:
  bool operator()(){
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    if(now < next) return false;
    next = now + std::chrono::milliseconds(POLL_MS);
    return user_interrupt();
  }

 private:
  std::chrono::steady_clock::time_point next = std::chrono::steady_clock::now() + std::chrono::milliseconds(POLL_MS);
};

// Why a run ended early. failed is set once: by the first worker error (its message kept), by a thread that could not start, or by
// an interrupt
class Stop {
 public:
  std::atomic<bool> failed{false};

  // In a worker's catch(...): records the exception being handled if it is the first failure
  void caught() noexcept {
    try {
      throw;
    } catch(const std::bad_alloc &){
      set("not enough memory");
    } catch(const std::exception &e){
      set(e.what());
    } catch(...){
      set("unknown error");
    }
  }

  void set(const char *msg) noexcept {
    if(failed.exchange(true)) return;
    try {
      message = msg;
    } catch(...){}
  }

  // Calling thread only
  void interrupt() noexcept {
    interrupted = true;
    failed.store(true);
  }
  bool was_interrupted() const {
    return interrupted;
  }

  // Calling thread, after the join: R's interrupt, or the R error "what: message"
  void raise(const char *what) const {
    if(interrupted) throw Rcpp::internal::InterruptedException();
    if(failed.load()) Rcpp::stop(std::string(what) + ": " + message);
  }

 private:
  std::string message;                  // written by the thread that set failed, read after the join
  bool interrupted = false;
};

const unsigned SPIN = 100;              // busy-wait rounds before a thread waiting at a barrier yields

// The members of a team wait here for each other: they spin, then yield after SPIN rounds. wait returns false once the run has
// failed (Stop::failed), which also ends the wait of members whose team will never be complete; sync throws instead, for code that
// unwinds to run_threads.
class Barrier {
 public:
  explicit Barrier(size_t n) : n(n) {}

  bool wait(const std::atomic<bool> &failed){
    const size_t ph = phase.load(std::memory_order_acquire);
    if(count.fetch_add(1, std::memory_order_acq_rel) + 1 == n){
      count.store(0, std::memory_order_relaxed);
      phase.store(ph + 1, std::memory_order_release);
      return !failed.load();
    }
    for(unsigned spin = 0; phase.load(std::memory_order_acquire) == ph; spin += spin <= SPIN){
      if(failed.load(std::memory_order_relaxed)) return false;
      if(spin > SPIN) std::this_thread::yield();
    }
    return !failed.load();
  }

  void sync(const std::atomic<bool> &failed){
    if(!wait(failed)) throw std::runtime_error("stopped: another thread failed");
  }

 private:
  const size_t n;
  std::atomic<size_t> count{0}, phase{0};
};

// fn(k) on n threads, k = 0..n-1, while the calling thread waits and polls for an interrupt. A thread that cannot start fails the run
// (the members of a team wait for each other at barriers, so a team cannot run short); the threads that did start see failed and
// return.
template <class Fn> void run_threads(size_t n, Stop &stop, Fn fn){
  std::mutex mu;
  std::condition_variable cv;
  size_t done = 0;
  auto body = [&](size_t k){
    try {
      lower_thread_priority();
      fn(k);
    } catch(...){
      stop.caught();
    }
    std::lock_guard<std::mutex> lock(mu);
    ++done;
    cv.notify_one();
  };
  std::vector<std::thread> pool;
  try {
    pool.reserve(n);
    for(size_t k = 0; k < n; ++k) pool.emplace_back(body, k);
  } catch(...){
    char msg[64];
    std::snprintf(msg, sizeof msg, "could not start %d threads", (int) n);
    stop.set(msg);
  }
  {
    std::unique_lock<std::mutex> lock(mu);
    while(!cv.wait_for(lock, std::chrono::milliseconds(POLL_MS), [&]{ return done == pool.size(); })){
      lock.unlock();
      if(!stop.was_interrupted() && user_interrupt()) stop.interrupt();
      lock.lock();
    }
  }
  for(std::thread &th : pool) th.join();
}

}  // namespace hc

#endif
