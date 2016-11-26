#ifndef SFR_CORE_GUI_IMOUSE_LISTENER_H
#define SFR_CORE_GUI_IMOUSE_LISTENER_H

namespace sfr { namespace gui
{

enum EMouseButton {    
    eMB_Left,
    eMB_Middle,
    eMB_Right    
};

class IMouseListener
{
    
public:
    IMouseListener() {}
    virtual ~IMouseListener() {}
    virtual bool OnMouseButtonPressed( EMouseButton mb, int x, int y ) { return false; }
    virtual bool OnMouseButtonReleased( EMouseButton mb, int x, int y ) { return false; }
    virtual bool OnMouseMotion( int x, int y ) { return false; }
};

} } // namespace sfr::gui

#endif // SFR_CORE_GUI_IMOUSE_LISTENER_H
