#ifndef SFR_CORE_ITASK_H
#define SFR_CORE_ITASK_H

namespace sfr
{

class ITask
{
protected:
    ITask() : m_LastTick(0) {}
    virtual ~ITask() {}    
    friend class IKernel;    
public:
    virtual void Run( unsigned int tick ) { m_LastTick = tick; }
    unsigned int GetLastTick() const { return m_LastTick; }
    
private:
    unsigned int m_LastTick;
};

} // namespace sfr

#endif // SFR_CORE_ITASK_H

