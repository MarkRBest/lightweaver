#pragma once

#include <iostream>
#include <stdint.h>
#include <sys/time.h>

#define CLOCK_MODE CLOCK_REALTIME

namespace lw {

#define ASM_CLOCK

#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}
#elif defined(__x86_64__)

static __inline__ uint64_t rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (uint64_t)lo)|( ((uint64_t)hi)<<32 );
}

#endif

// Used for nanosecond level timing of functions
class Timer {

private:
    timespec tsStart;
    int64_t time1;

public:
    Timer() {
        reset();
    }

    inline void reset() {
#ifndef ASM_CLOCK
        clock_gettime(CLOCK_MODE, &tsStart);
#else
        time1 = rdtsc();
#endif
    }

    inline int64_t elapsedTime() {
#ifndef ASM_CLOCK
        timespec tsEnd;
        clock_gettime(CLOCK_MODE, &tsEnd);
        return (tsEnd.tv_sec - tsStart.tv_sec) * 1E9L
                + (tsEnd.tv_nsec - tsStart.tv_nsec);
#else
        return (rdtsc() - time1)/2;
#endif
    }

    friend std::ostream& operator <<(std::ostream& stream, Timer& timer) {
        int64_t et = timer.elapsedTime();
        stream << "Timer {";
        stream << et << " ns, ";
        stream << et / 1000 << " mics, ";
        stream << et / 1000000 << " ms";
        stream << "}";
        return stream;
    }
};

} // namespace lw
