#ifndef TIME_WINDOWS_H
#define TIME_WINDOWS_H


// Thanks to this stackoverflow post for providing a windows drop-in for gettimeofday
// https://stackoverflow.com/questions/34831145/using-gettimeofday-equivalents-on-windows

#include <stdio.h>
#include <Windows.h>
#include <stdint.h>
#include <sys/timeb.h>
#include <time.h>

inline int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    time = ((uint64_t)file_time.dwLowDateTime);
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec = (long)((time - EPOCH) / 10000000L);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
    return 0;
}

#endif
