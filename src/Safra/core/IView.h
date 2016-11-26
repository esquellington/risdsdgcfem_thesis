#ifndef SFR_CORE_IVIEW_H
#define SFR_CORE_IVIEW_H

#include <Safra/task/IUpdaterTask.h>
#include <Safra/core/gui/IDesktop.h>
#include <Safra/core/gui/IMouseListener.h>
#include <Safra/core/gui/IKeyboardListener.h>

#include <Safra/core/gfx/IRenderer.h>
#include <Safra/core/gfx/ICamera.h>

namespace sfr
{

/*! View Interface, represents a dimensionless window that accepts
  Updates and Mouse/Keyboard inputs.
*/
class IView: public IUpdateTarget
{
public:
    typedef enum {
        eFlags_None = 0,
        eFlags_EnableDesktop = 1<<0,
        eFlags_EnableFPS     = 1<<1,
        eFlags_Default = eFlags_None
    } EViewFlags;

public:
    IView() {}
    virtual ~IView() {}

    //! \name Visualization Messages
    //@{
    virtual void BeginDraw() {}
    virtual bool Draw() = 0;
    virtual void EndDraw() {}

    virtual void Reshape( int w, int h ) = 0;
    //@}

    //! \name Consultors
    //@{
    virtual void GetViewportSizes( int &w, int &h ) const = 0;
    //@}

    //!\name Camera Management
    //@{
    virtual gfx::ICamera *GetCamera() = 0;
    virtual void SetCamera( gfx::ICamera *p_camera ) = 0;
    //@}

    //!\name Renderer Management
    //@{
    virtual void SetRenderer( gfx::IRenderer *p_renderer ) = 0;
    virtual gfx::IRenderer *GetRenderer() const = 0;
    //@}

    //!\name GUI related methods
    //@{
    virtual Flags32 GetFlags() const = 0;
    virtual void SetFlags( Flags32 flags ) = 0;

    virtual void SetMouseListener( gui::IMouseListener *p_mouse_listener ) = 0;
    virtual void SetKeyboardListener( gui::IKeyboardListener *p_keyboard_listener ) = 0;

    virtual gui::IDesktop *GetDesktop() const { return 0; }
    virtual gui::IMouseListener *GetMouseListener() const = 0;
    virtual gui::IKeyboardListener *GetKeyboardListener() const = 0;
    //}@

    //!\name Mouse events
    //@{
    virtual bool OnMouseButtonPressed( gui::EMouseButton mb, int x, int y )
    {
        //\todo While !pListener->OnMouseButtonPressed() , next...
        if( GetDesktop() && GetDesktop()->IsEnabled() && GetDesktop()->OnMouseButtonPressed( mb, x, y ) ) return true;
        else if( GetMouseListener() && GetMouseListener()->OnMouseButtonPressed( mb, x, y ) ) return true;
        else return false;
    }
    virtual bool OnMouseButtonReleased( gui::EMouseButton mb, int x, int y )
    {
        //\todo While !pListener->OnMouseButtonReleased() , next...
        if( GetDesktop() && GetDesktop()->IsEnabled() && GetDesktop()->OnMouseButtonReleased( mb, x, y ) ) return true;
        else if( GetMouseListener() && GetMouseListener()->OnMouseButtonReleased( mb, x, y ) ) return true;
        else return false;
    }
    virtual bool OnMouseMotion( int x, int y )
    {
        //\todo While !pListener->OnMouseMotion() , next...
        if( GetDesktop() && GetDesktop()->IsEnabled() && GetDesktop()->OnMouseMotion( x, y ) ) return true;
        else if( GetMouseListener() && GetMouseListener()->OnMouseMotion( x, y ) ) return true;
        else return false;
    }
    //@}

    //!\name Keyboard events
    //@{
    virtual bool OnKeyPressed( gui::EKey key, int x, int y )
    {
        //\todo While !pListener->OnKeyPressed() , next...
        if( GetDesktop() && GetDesktop()->IsEnabled() && GetDesktop()->OnKeyPressed( key, x, y ) ) return true;
        else if( GetKeyboardListener() && GetKeyboardListener()->OnKeyPressed( key, x, y ) ) return true;
        else return false;
    }
    virtual bool OnKeyReleased( gui::EKey key, int x, int y )
    {
        //\todo While !pListener->OnKeyReleased() , next...
        if( GetDesktop() && GetDesktop()->IsEnabled() && GetDesktop()->OnKeyReleased( key, x, y ) ) return true;
        else if( GetKeyboardListener() && GetKeyboardListener()->OnKeyReleased( key, x, y ) ) return true;
        else return false;
    }
    //@}
};

class ITextView: public IView
{
public:
    ITextView() {}
    virtual ~ITextView() {}

    virtual void SetFontSize( int size_in_pixels ) = 0;
    virtual void SetDefaultColor( const gfx::Color &color ) = 0;

    virtual void GetCursor( int &row, int &col ) const = 0;
    virtual void SetCursor( int row, int col ) = 0;
    virtual void MoveCursor( int rows, int cols ) = 0;
    virtual void SetCursorVisible( bool b_visible ) = 0;

    virtual char Read() const = 0;

    virtual int Write( const char *text, int len ) = 0;
    virtual int Write( const char *text ) = 0;
    virtual int Write( char c ) = 0;

    virtual int Space( int count ) = 0;
    virtual int Tab() = 0;

    virtual int BackSpace() = 0;
    virtual int Del() = 0;

    virtual void Endl() = 0;

    virtual void Clear() = 0;

    virtual const char *GetText() const = 0;
    virtual const char *GetLine() const = 0;
};


} // namespace sfr

#endif // SFR_CORE_IVIEW_H
