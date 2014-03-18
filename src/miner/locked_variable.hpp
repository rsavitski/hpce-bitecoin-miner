#ifndef locked_variable_hpp
#define locked_variable_hpp
#include <mutex>
#include <condition_variable>
#include <functional>

template <typename T> class LockedVariable {
  std::mutex m;
  std::condition_variable cv;

public:
  T data; // DO NOT ACCESS THIS WITHOUT OBTAINING A LOCK!
  LockedVariable() {}
  LockedVariable(T data) : data(data) {}
  LockedVariable(const LockedVariable &obj) = delete;
  void operator=(const LockedVariable &obj) = delete;

  void lockAndDo(std::function<void(T &)> action) {
    std::unique_lock<std::mutex> lk = lock();
    action(data);
    lk.unlock();
  }

  std::unique_lock<std::mutex> waitFor(std::function<bool(T &)> condition) {
  	std::unique_lock<std::mutex> lk = lock();
  	cv.wait(lk, condition);
  	return lk;
  }

  std::unique_lock<std::mutex> lock() {
    return std::unique_lock<std::mutex>(m);
  }

  void unlock(std::unique_lock<std::mutex> &&lk) { lk.unlock(); }

  void notify() { cv.notify_all(); }
};

#endif
