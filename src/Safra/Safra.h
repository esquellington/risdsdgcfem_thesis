#ifndef SFR_SAFRA_H
#define SFR_SAFRA_H

#include <Safra/core/IKernel.h>

namespace util { class LogStream; }

class Safra
{
public:
    enum EKernelType {
        eKernelNone = 0,
        eKernelConsole,
        eKernelGLUT
    };

public:    
    static bool Init( EKernelType kt );
    static bool Run();
    static bool ShutDown();

    static bool AddTask( sfr::ITask *p_task );
    
    static sfr::IView *CreateView( const char *name, unsigned int flags, int width, int height );
    static sfr::ITextView *CreateTextView( const char *name, unsigned int flags, 
                                           int width, int height,
                                           int text_rows, int text_columns );

    static sfr::gfx::IRenderer *GetDefaultRenderer();
    
    static void SetLogStream( util::LogStream *p_log_stream );
    static util::LogStream &GetLogStream();

private:
    Safra() {}
    ~Safra() {}
    
private:
    static EKernelType m_KernelType;
    static sfr::IKernel *m_pKernel;
    static util::LogStream *m_pLogStream;
};

#endif // SFR_SAFRA_H
