/* Copyright 2014 - Mario Kostelac (mario.kostelac@gmail.com) */
#ifndef SRC_TIMER_TIMER_H_
#define SRC_TIMER_TIMER_H_
#include <ctime>

class Timer {
 public:
     explicit Timer(const char * name);
     ~Timer();
     Timer* end(bool print = true);

 private:
     char * name;
     std::clock_t start_time;
     std::clock_t end_time;
};
#endif  // SRC_TIMER_TIMER_H_
