#ifndef SFR_CORE_IKERNEL_H
#define SFR_CORE_IKERNEL_H

#include <Safra/core/ITask.h>
#include <Safra/core/IView.h>
#include <Safra/core/gfx/IRenderer.h>

namespace sfr
{

class IKernel
{
public:

    virtual ~IKernel() {}

    // Task Management
    virtual bool AddTask( ITask *p_task ) = 0;

    // View Management
    virtual IView *CreateView( const char *name, Flags32 flags, int width, int height ) = 0;
    virtual ITextView *CreateTextView( const char *name, Flags32 flags,
                                       int width, int height,
                                       int text_rows, int text_columns ) = 0;

    // Renderer Management
    virtual sfr::gfx::IRenderer *GetDefaultRenderer() const = 0;

    // Execution
    virtual void Run() = 0; //!< May not end
    virtual void Kill() = 0;
};

} // namespace sfr

#endif // SFR_CORE_IKERNEL_H
