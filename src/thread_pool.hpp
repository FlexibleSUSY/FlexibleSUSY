// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include "logger.hpp"

#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

namespace flexiblesusy {

/**
 * @class Thread_pool
 * @brief A pool of threads
 *
 * Thread_pool represents a collection of threads.  Tasks (callables)
 * can be added to an internal queue.  The tasks will be executed as
 * soon as there is an idle thread.  The destructor of the Thread_pool
 * will wait until all tasks are finished and the queue is empty.
 *
 * @param pool_size number of threads in the pool
 */
class Thread_pool {
public:
   explicit Thread_pool(std::size_t pool_size = std::thread::hardware_concurrency())
   {
      VERBOSE_MSG("launching " << pool_size << " threads ...");
      for (std::size_t i = 0; i < pool_size; ++i)
         threads.emplace_back(
            [this] () {
               for (;;) {
                  std::function<void()> task;

                  {
                     std::unique_lock<std::mutex> lock(mutex);
                     condition.wait(lock, [this]{ return stop || !tasks.empty(); });
                     if (stop && tasks.empty())
                        return;
                     task = std::move(tasks.front());
                     tasks.pop();
                  }

                  task();
               }
            });
   }

   Thread_pool(const Thread_pool&) = delete;
   Thread_pool(Thread_pool&&) = delete;
   Thread_pool& operator=(const Thread_pool&) = delete;
   Thread_pool& operator=(Thread_pool&&) = delete;

   /// waits for all tasks to finish and closes threads
   ~Thread_pool()
   {
      {
         std::unique_lock<std::mutex> lock(mutex);
         stop = true;
      }

      condition.notify_all();

      try {
         for (auto& t: threads)
            t.join();
      } catch (const std::exception& e) {
         ERROR(e.what());
      }
   }

   /// runs task and returns future
   template <typename Task>
   auto run_packaged_task(Task&& task) -> std::future<decltype(task())>
   {
      using return_t = decltype(task());

      auto ptask = std::make_shared<std::packaged_task<return_t()>>([task](){ return task(); });

      std::future<return_t> fut = ptask->get_future();

      if (threads.empty()) {
         (*ptask)();
      } else {
         {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.emplace([ptask](){ (*ptask)(); });
         }
         condition.notify_one();
      }

      return fut;
   }

   /// runs task
   template <typename Task>
   void run_task(Task&& task)
   {
      if (threads.empty()) {
         task();
      } else {
         {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.emplace(std::forward<Task>(task));
         }
         condition.notify_one();
      }
   }

   std::size_t size() const { return threads.size(); }

private:
   std::vector<std::thread> threads{};
   std::queue<std::function<void()>> tasks{};
   std::mutex mutex{};
   std::condition_variable condition{};
   bool stop{false};
};

} // namespace flexiblesusy

#endif
