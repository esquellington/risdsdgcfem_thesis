#ifndef SFR_GUI_BUTTON_CONTROL_H
#define SFR_GUI_BUTTON_CONTROL_H

#include "BaseControl.h"
#include "LabelControl.h"
#include <boost/function.hpp>

namespace sfr { namespace gui
{

/*! Simple text button that calls a callback when mouse press-release inside
  m_bIsElastic: Controls wether the button automatically returns to unpressed when released.
*/
class ButtonControl: public LabelControl
{
public:    
    typedef boost::function<void (void *)> D_OnButtonCall;
    
public:
    ButtonControl( const char *name = 0, void *p_user_data = 0 )
    : m_pUserData(p_user_data)
    , m_FontColorReleased(0.5,0.5,0.5,1)
    , m_FontColorPressed(1,1,1,1)
    , m_bIsElastic(true)
    , m_bIsPressed(false)
    {
        if( name ) SetLabel( name );
        else SetLabel( "unnamed" );
        LabelControl::SetFontColor( m_FontColorReleased );
    }
    ~ButtonControl() {}
    
    inline void SetOnButtonCall( D_OnButtonCall on_button_call ) { m_OnButtonCall = on_button_call; }
    
    inline void SetUserData( void *p_user_data ) { m_pUserData = p_user_data; }
    inline void *GetUserData() const { return m_pUserData; }
    
    inline void SetButtonColors( const gfx::Color &color_released, const gfx::Color &color_pressed )
    {
        m_FontColorReleased = color_released;
        m_FontColorPressed = color_pressed;
        SetFontColor( (m_bIsPressed) ? m_FontColorPressed : m_FontColorReleased );
    }
    
    inline void SetFontColor( const gfx::Color &color )
    {
        m_FontColorReleased = color;
        m_FontColorPressed = color;
        LabelControl::SetFontColor(color);
    }

    inline void SetElastic( bool b_elastic ) { m_bIsElastic = b_elastic; }
    inline void SetPressed( bool b_pressed ) //! Mostly intended for non-elastic buttons
        {
            m_bIsPressed = b_pressed;
            LabelControl::SetFontColor( (m_bIsPressed) ? m_FontColorPressed : m_FontColorReleased );
        }
    
    //! Toggle color when pressed
    bool OnMouseButtonPressed( EMouseButton mb, Vec2 p )
    {
        if( !IsPointInside(p) )
            return false;
        else
        {
            if( m_bIsElastic ) m_bIsPressed = true;
            //else m_bIsPressed = !m_bIsPressed;
            LabelControl::SetFontColor( (m_bIsPressed) ? m_FontColorPressed : m_FontColorReleased );
            return true;
        }
    }

    //! Call delegate and toggle color when released
    bool OnMouseButtonReleased( EMouseButton mb, Vec2 p )
    {
        if( !IsPointInside(p) )
        {
            if( HasFocus() ) OnFocusLost(); //Explicitly lose focus when released outside
            return false;
        }
        else
        {
            if( m_bIsElastic ) m_bIsPressed = false;
            else m_bIsPressed = !m_bIsPressed;
            LabelControl::SetFontColor( (m_bIsPressed) ? m_FontColorPressed : m_FontColorReleased );
            if( !m_OnButtonCall.empty() ) m_OnButtonCall(m_pUserData);
            else SFR_LOG( "ButtonControl(%s): Invalid D_OnButtonCall", GetLabel() );
            return true;
        }
    }
    
    //! Toggle color if focus lost
    void OnFocusLost()
    {
        if( m_bIsElastic ) m_bIsPressed = false;
        LabelControl::SetFontColor( (m_bIsPressed) ? m_FontColorPressed : m_FontColorReleased );
        LabelControl::OnFocusLost();
    }
                       
private:
    void *m_pUserData;
    D_OnButtonCall m_OnButtonCall;
    gfx::Color m_FontColorReleased;
    gfx::Color m_FontColorPressed;
    bool m_bIsElastic;
    bool m_bIsPressed;
};

}} //namespace sfr::gui

#endif //SFR_GUI_BUTTON_CONTROL_H
