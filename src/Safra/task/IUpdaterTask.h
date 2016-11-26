#ifndef SFR_TASK_IUPDATER_TASK_H
#define SFR_TASK_IUPDATER_TASK_H

#include <Safra/core/ITask.h>
#include <util/Chronos.h>

namespace sfr
{

//! Anything receiving timed updates
class IUpdateTarget
{
public:
    virtual bool Update( double dt ) = 0;
};

//! Task that updates an UpdateTarget with a given policy
class IUpdaterTask: public ITask
{
public:
    IUpdaterTask( IUpdateTarget *p_ut ) : m_pUpdateTarget(p_ut) {}
    ~IUpdaterTask() {}
protected:
    IUpdateTarget *m_pUpdateTarget;
};

/*! Fixed timestep update at every tick. */
class IdleUT: public IUpdaterTask
{
public:
    IdleUT( IUpdateTarget *p_ut, float timestep )
    : IUpdaterTask( p_ut )
    , m_TimeStep(timestep)
    {}
    ~IdleUT() {}

    void Run( unsigned int tick )
    {
        m_pUpdateTarget->Update( m_TimeStep );
    }
    
private:
    float m_TimeStep;
};

/*! Realtime updates at fixed, minimum or variable
  - timestep > 0 => Fixed dt (ensure constant frequency)
  - timestep < 0 => Minimum dt (ensure maximum frequency)
  - timestep = 0 => dt = realtime step (variable frequency)
*/
class RealtimeUT: public IUpdaterTask
{
public:
    RealtimeUT( IUpdateTarget *p_ut, float timestep = 0.0f )
    : IUpdaterTask( p_ut )
    , m_TimeStep(timestep)
    , m_TotalTime(0)
    , m_AccTime(0)
    , m_TotalSteps(0)        
    {}

    void Run( unsigned int tick )
    {
        double dt = m_Chronos.GetDt();
        m_TotalTime += dt;
    
        // Check if fixed or variable timestep
        if( m_TimeStep > 0.0 )
        {
            m_AccTime += dt;
            while( m_AccTime >= m_TimeStep )
            {
                m_pUpdateTarget->Update( m_TimeStep );
                m_AccTime -= m_TimeStep;
                m_TotalSteps++;
            }                  
            m_TotalTime -= m_AccTime; // Restem part del temps que NO calculem ara
        }
        else if( m_TimeStep < 0.0f )
        {
            m_AccTime += dt;
            if( m_AccTime >= -m_TimeStep )            
                m_pUpdateTarget->Update( -m_TimeStep );            
            m_TotalSteps++;            
            while( m_AccTime >= -m_TimeStep )
                m_AccTime += m_TimeStep;            
            m_TotalTime -= m_AccTime; // Restem part del temps que NO calculem ara
        }
        else // m_TimeStep == 0.0
        {
            m_pUpdateTarget->Update( dt );
            m_TotalSteps++;
        }
    }        
    
private:
    float m_TimeStep;
    double m_TotalTime;
    double m_AccTime;
    unsigned int m_TotalSteps;
    util::Chronos m_Chronos;
};

} // namespace sfr

#endif // SFR_TASK_IUPDATER_TASK_H
