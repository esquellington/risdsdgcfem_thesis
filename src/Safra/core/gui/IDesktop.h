#ifndef SFR_CORE_GUI_IDESKTOP_H
#define SFR_CORE_GUI_IDESKTOP_H

#include <Safra/Config.h>
#include "IMouseListener.h"
#include "IKeyboardListener.h"
#include "IWidget.h"

namespace sfr { namespace gui
{

typedef boost::function<void (PropertyIt)> D_PropertyTreeSyncCall;

class IDesktop: public IMouseListener, public IKeyboardListener
{
public:
    IDesktop() {}
    virtual ~IDesktop() {}

    virtual void Message( const char *message, float duration ) = 0;
    virtual IWidget *CreateWidget( const char *name, Flags32 flags, Vec2f pos ) = 0;
    virtual IWidget *CreateWidget_Tree( const char *name, Flags32 flags, Vec2 pos ) = 0;
    virtual IWidget *CreateWidget_PropertyTree( const char *name, Flags32 flags, Vec2f pos,
                                                PropertyIt property_tree_it,
                                                D_PropertyTreeSyncCall sync_call = D_PropertyTreeSyncCall() ) = 0;

    virtual void Resize( int w, int h ) = 0;
    virtual Vec2 GetSizes() const = 0;
    virtual bool Draw() = 0;
    virtual bool Update( float dt ) = 0;

    virtual void MinimizeAll() = 0;
    virtual void MaximizeAll() = 0;
    virtual void Clear() = 0;

    virtual bool IsEnabled() = 0;
    virtual void SetEnabled( bool b_enabled ) = 0;

    virtual void SetBackgroundColor( float r, float g, float b, float a ) = 0;
};

} } // namespace sfr::gui

#endif // SFR_CORE_GUI_IDESKTOP_H
