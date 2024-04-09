#include <cstddef>
#include <sys/time.h>

class Timer {
    private:
        struct timeval m_StartTime, m_EndTime;
    public:
        void Start() {
            gettimeofday(&m_StartTime, NULL);
        }
        long Stop() {
            long useconds, seconds, mseconds;
            gettimeofday(&m_EndTime, NULL);
            useconds = m_EndTime.tv_usec - m_StartTime.tv_usec;
            seconds = m_EndTime.tv_sec - m_StartTime.tv_sec;
            mseconds = ((seconds) * 1000 + useconds/1000.0 + 0.5);
            return mseconds;
        }
        
};

