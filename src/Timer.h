#ifndef TIMER_H
#define TIMER_H

#include <stdexcept>
#include <sys/time.h>
#include <time.h>

class Timer
{
  private:
    bool running;
    double accumulated_cpu_time;
    double accumulated_wall_time;
    double start_cpu_time;
    double start_wall_time;

    double get_cpu_time() const;
    double get_wall_time() const;

  public:
    Timer(const bool = false);
    tm *current_time() const;
    double elapsed_cpu_time() const;
    double elapsed_wall_time() const;
    void reset();
    void restart();
    void start();
    void stop();
};

#endif