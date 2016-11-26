#ifndef UTIL_CHRONOS_H
#define UTIL_CHRONOS_H

namespace util
{

//! God of Time
/*! This class keeps track of the aplication time.
  IMPORTANT: All times are given in seconds.

  Methods are provided to:
  - Get the current time.
  - Get the time interval (dt) between two moments
  - Get the last time interval.

  The tipical sequence will be:
  1 - Reset at initialization, new level, etc... before starting to compute frames.
  2 - GetDt at the beginning of the current frame
  3 - GetCurrentDt when the dt value is needed during the computation of the current frame.
  4 - Render
  5 - Goto 3

*/
class Chronos
{
public:
    typedef double time_type;
    /*! Returns the system time from the beginning of the app execution.
      The implementation of this method is platform dependent (win32 != Linux). */
    static time_type GetSysTime_s();
    static time_type GetSysTime_ms();
    static time_type GetSysTime_us();

public:
    Chronos()
    {
        ResetTime();
    }

    /*! Resets current time to 0. */
    void ResetTime();

    /*! Returns current time. */
    time_type GetTime();

    /*! Returns the time between the last call to GetTime or GetDt. */
    time_type GetDt();

    /*! Returns the current value of dt (the time interval) without changing it.
      The dt value will remain constant while GetDt() is not called. */
    time_type GetCurrentDt() const { return m_CurrentDt; }

private:

    time_type m_CurrentDt;
    time_type m_BaseTime;
    time_type m_CurrentTime;
};

} // namespace util

#endif // UTIL_CHRONOS_H
