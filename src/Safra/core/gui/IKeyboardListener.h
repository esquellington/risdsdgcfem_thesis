#ifndef SFR_CORE_GUI_IKEYBOARD_LISTENER_H
#define SFR_CORE_GUI_IKEYBOARD_LISTENER_H

namespace sfr { namespace gui
{

enum EKey {
    
    // Normal (=ASCII code)
    eKey_BackSpace = 8,
    eKey_Tab = 9,
    eKey_Return = 13,
    eKey_Escape = 27,
    eKey_Space = 32,
    eKey_Del = 127,
    
    // Modifiers (>255)
    eKey_Control = 256,
    eKey_Alt,
    eKey_Shift,

    // Specials (>255)    
    eKey_Up,
    eKey_Down,
    eKey_Left,
    eKey_Right,

    eKey_F1,
    eKey_F2,
    eKey_F3,
    eKey_F4,
    eKey_F5,
    eKey_F6,
    eKey_F7,
    eKey_F8,
    eKey_F9,
    eKey_F10,
    eKey_F11,
    eKey_F12,

    // Invalid
    eKey_Invalid = 0xFFFFFFFF,    
};

class IKeyboardListener
{

public:
    IKeyboardListener() {}
    virtual ~IKeyboardListener() {}

    virtual bool OnKeyPressed( EKey key, int x, int y ) { return false; }
    virtual bool OnKeyReleased( EKey key, int x, int y ) { return false; }
};

} } // namespace sfr::gui

#endif // SFR_CORE_GUI_IKEYBOARD_LISTENER_H
