#ifndef SFR_GUI_DESKTOP_H
#define SFR_GUI_DESKTOP_H

#include <Safra/core/gui/IDesktop.h>
#include "BaseControl.h"
#include "WidgetControl.h"
#include <string>

namespace sfr { class IView; }

namespace sfr { namespace gui
{

class MinimizedControl;

class Desktop: public IDesktop, public BaseControl
{
public:
    Desktop( IView *p_view );
    ~Desktop();

    //!\name IDesktop implementation
    //@{
    void Message( const char *message, float duration );

    IWidget *CreateWidget( const char *name, Flags32 flags, Vec2 pos );
    IWidget *CreateWidget_Tree( const char *name, Flags32 flags, Vec2 pos );
    IWidget *CreateWidget_PropertyTree( const char *name, Flags32 flags, Vec2 pos,
                                        PropertyIt property_tree_it,
                                        D_PropertyTreeSyncCall sync_call = D_PropertyTreeSyncCall() );

    void Resize( int w, int h );
    Vec2 GetSizes() const { return BaseControl::GetSizes(); }
    bool Draw();
    bool Update( float dt );

    void MinimizeAll() {}
    void MaximizeAll() {}
    void Clear() {}

    bool IsEnabled() { return m_bIsEnabled; }
    void SetEnabled( bool b_enabled ) { m_bIsEnabled = b_enabled; }

    void SetBackgroundColor( float r, float g, float b, float a );
    //@}

public:
    bool OnMouseButtonPressed( EMouseButton mb, int x, int y );
    bool OnMouseButtonReleased( EMouseButton mb, int x, int y );
    bool OnMouseMotion( int x, int y );
    bool OnKeyPressed( EKey key, int x, int y );
    bool OnKeyReleased( EKey key, int x, int y );

    void Minimize( WidgetControl *p_widget );
    void Restore( WidgetControl *p_widget );

private:
    void OnButtonCall_CreateWidget( void *p_user_data );
    void OnButtonCall_CreateTree( void *p_user_data );
    void OnButtonCall_CreatePropertyTree( void *p_user_data );

private:
    IView *m_pView;
    bool m_bIsEnabled;
    bool m_bUpdatedSinceLastDraw;

    double m_MessageTimeOut;
    std::string m_MessageString;
    MinimizedControl *m_pMinimizedControl;
};

} } // namespace sfr::gui

#endif // SFR_GUI_DESKTOP_H
