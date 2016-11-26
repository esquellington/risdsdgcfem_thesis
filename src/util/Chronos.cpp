#include <util/Chronos.h>

#ifdef WIN32
#include <windows.h> // Windows timer related methods
#else
#include <sys/time.h>
#endif

namespace util
{

void Chronos::ResetTime()
{
    m_BaseTime = GetSysTime_s();
    m_CurrentTime = m_BaseTime;
}

Chronos::time_type Chronos::GetTime()
{
    m_CurrentTime = GetSysTime_s();
    return (m_CurrentTime - m_BaseTime);
}

Chronos::time_type Chronos::GetDt()
{
    Chronos::time_type last_time = m_CurrentTime;
    m_CurrentTime = GetSysTime_s();
    m_CurrentDt = m_CurrentTime - last_time;
    return m_CurrentDt;
}

Chronos::time_type Chronos::GetSysTime_s()
{
#ifdef WIN32
    /*
      __declspec(naked) unsigned long GetCounter()
      {
      __asm { rtdsc }
      }
    */
    /* Nebula2
       Chronos::time_type
       nTimeServer::GetTime()
       {
       if (this->lock_delta_t > 0.0) return this->lock_time;
       else {
       #       ifdef __WIN32__
       LONGLONG time,freq;
       QueryPerformanceCounter((LARGE_INTEGER *)&time);
       if (this->stopped) time = this->time_stop;
       QueryPerformanceFrequency((LARGE_INTEGER *)&freq);
       LONGLONG td = time - this->time_diff;
       Chronos::time_type d_time = ((Chronos::time_type)td) / ((Chronos::time_type)freq);
       return d_time;
       return 0.0;
       }
       #   elif defined(__LINUX__) || defined(__MACOSX__)
       long long int time;
       long long int td;
       Chronos::time_type d_time;
       if (this->stopped) time = this->time_stop;
       else {
       struct timeval tv;
       gettimeofday(&tv,NULL);
       time = tv2micro(tv);
       }
       td = time - this->time_diff;
       d_time = ((Chronos::time_type)td) / N_MICROSEC_FLOAT;
       return d_time;
       }
       #   else
       #error "Method not implemented!"
       #   endif
       }
    */

    return ((Chronos::time_type)GetTickCount())*1000.0f;
#else
    // gettimeofday has a precision of 1 microsecond (tested in test_Util)
    Chronos::time_type v;
    struct timeval tv;
    gettimeofday(&tv, 0);
    v = (Chronos::time_type) tv.tv_sec + (Chronos::time_type) tv.tv_usec / 1000000.0;
    return v;
#endif
}

Chronos::time_type Chronos::GetSysTime_ms() { return GetSysTime_s()*1000.0; }
Chronos::time_type Chronos::GetSysTime_us() { return GetSysTime_s()*1000000.0; }

}; // namespace util
