/*!
  \todo Move code to .cpp
*/

#ifndef SFR_CORE_KERNEL_CONSOLE_KERNEL_H
#define SFR_CORE_KERNEL_CONSOLE_KERNEL_H

#include <Safra/core/IKernel.h>
#include <list>
#ifdef WIN32
#include <conio.h>
#else
#include <Safra/core/platform/posix/conio_nonstandard.h>
#endif

namespace sfr
{

class ConsoleKernel: public IKernel
{
public:
    ConsoleKernel()
    : m_LastTick(0)
    {}
    ~ConsoleKernel() {}

    bool AddTask( ITask *p_task )
    {
        m_listTasks.push_back( p_task );
        return true;
    }

    IView *CreateView( const char *name, Flags32 flags, int width, int height ) { return NULL; }
    ITextView *CreateTextView( const char *name, Flags32 flags,
                               int width, int height,
                               int text_rows, int text_columns ) { return NULL; }

    sfr::gfx::IRenderer *GetDefaultRenderer() const { return NULL; }

    void Run()
    {
        while( !(kbhit() && getch() == 'q') )
        {
            for(std::list<ITask *>::iterator it_task = m_listTasks.begin();
                it_task != m_listTasks.end();
                ++it_task)
                (*it_task)->Run(m_LastTick);
            m_LastTick++;
        }
    }

    void Kill()
    {
        exit(0);
    }

private:
    unsigned int m_LastTick;
    std::list<ITask *> m_listTasks;
};

} // namespace sfr

#endif // SFR_CORE_KERNEL_CONSOLE_KERNEL_H
