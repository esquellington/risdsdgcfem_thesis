#include <Safra/Safra.h>
#include <Safra/core/kernel/ConsoleKernel.h>
#include <Safra/core/kernel/GLUTKernel.h>
#include <util/LogStream.h>

#include <iostream>

Safra::EKernelType Safra::m_KernelType;
sfr::IKernel *Safra::m_pKernel;
util::LogStream *Safra::m_pLogStream = NULL;
    
bool Safra::Init( EKernelType kt )
{    
    std::cout << "Initializing Kernel: ";
    m_KernelType = kt;
    switch( m_KernelType )
    {
    case eKernelConsole:
        m_pKernel = new sfr::ConsoleKernel();
        std::cout << "Console" << std::endl;
        break;
    case eKernelGLUT:
        m_pKernel = new sfr::GLUTKernel();
        std::cout << "GLUT" << std::endl;
        break;
    default:
        std::cout << "Unknown!" << std::endl;
        return false;
    }
    return true;
}

bool Safra::Run()
{
    /*\todo logstream should be an ItemStream
    if( m_pLogStream )
        (*m_pLogStream) << "Run Lola Run!!" << util::Endl();
    */
    if( !m_pKernel )
        return false;
    m_pKernel->Run();
    return true;
}

bool Safra::ShutDown()
{    
    if( m_pKernel )
    {
        m_pKernel->Kill();
        delete m_pKernel;
    }
    return true;
}

bool Safra::AddTask( sfr::ITask *p_task )
{
    if( !m_pKernel )
        return false;
    return m_pKernel->AddTask( p_task );
}

sfr::IView *Safra::CreateView( const char *name, unsigned int flags, int width, int height )
{
    if( !m_pKernel )
        return NULL;
    return m_pKernel->CreateView(name,flags,width,height);
}

sfr::ITextView *Safra::CreateTextView( const char *name, unsigned int flags,
                                       int width, int height,
                                       int text_rows, int text_columns )
{
    if( !m_pKernel )
        return NULL;
    return m_pKernel->CreateTextView(name,flags,width,height,text_rows,text_columns);
}


sfr::gfx::IRenderer *Safra::GetDefaultRenderer()
{
    if( !m_pKernel )
        return NULL;
    return m_pKernel->GetDefaultRenderer();    
}

void Safra::SetLogStream( util::LogStream *p_log_stream )
{
    m_pLogStream = p_log_stream;
}

util::LogStream &Safra::GetLogStream()
{
    return *m_pLogStream;
}
