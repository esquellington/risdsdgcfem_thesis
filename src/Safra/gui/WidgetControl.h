#ifndef SFR_GUI_WIDGET_CONTROL_H
#define SFR_GUI_WIDGET_CONTROL_H

#include <Safra/core/gui/IWidget.h>
#include "BaseControl.h"
#include "LabelControl.h"
#include "ButtonControl.h"
#include "LayoutControl.h"
#include <vector>

namespace sfr { namespace gui
{

class Desktop;
class WidgetControl;

/*! Static text caption with click-and-drag control */
class Caption: public LabelControl
{
public:
    Caption( WidgetControl *p_widget );
    ~Caption();
    
    //!\name Specializations
    //@{
    bool OnMouseButtonPressed( EMouseButton mb, Vec2 p );
    bool OnMouseButtonReleased( EMouseButton mb, Vec2 p );
    bool OnMouseMotion( Vec2 p );
    //@}
        
private:
    WidgetControl *m_pWidget;
    
    //!\name Interaction status
    //@{
    enum EInteractionState { eIS_Standby, eIS_Moving } m_InteractionState;
    Vec2 m_InteractionPos0;
    //@}
};

/*! Static text headbar that allows:
  - Dragging
  - Roll/Unroll
  - Minimize
  - Close
*/    
class HeadBar: public HorizontalLayout
{
public:
    HeadBar( WidgetControl *p_widget );
    ~HeadBar();
    
    const Caption &GetCaption() const { return m_Caption; }
    Caption &GetCaption() { return m_Caption; }

    //!\name Action Button Calls
    //@{
    void ToggleRoll( void *p_user_data );
    void ToggleMinimize( void *p_user_data );
    void Close( void *p_user_data );
    //@}

    void NotifyRolledChange();
    void NotifyMinimizedChange();
        
private:
    WidgetControl *m_pWidget;
        
    //!\name Sub-controls
    //@{
    Caption m_Caption;
    ButtonControl m_buttonToggleRoll;
    ButtonControl m_buttonToggleMinimize;
    ButtonControl m_buttonClose;
    //@}
};

/*! FootBar area that accepts externally defined buttons.
    \note It's ONLY drawn if it contains N>0 buttons, otherwise the
    layout has empty sizes and is not drawn.
*/
class FootBar: public HorizontalLayout
{
public:
    FootBar( WidgetControl *p_widget );
    ~FootBar();
    bool AddButton( const char *name, ButtonControl::D_OnButtonCall on_button_call, void *p_user_data );
    bool RemoveButton( const char *name );
        
private:
    WidgetControl *m_pWidget;
    std::vector<ButtonControl *> m_vecButtons;
};    

/*! Gui Widget
  Encapsulates a movable panel with:
  - HeadBar
  - Externally defined body
  - FootBar
*/
class WidgetControl: public VerticalLayout, public IWidget
{
    typedef boost::function<void (WidgetControl *)> D_OnCloseCall;
    
public:
    WidgetControl( gfx::Color color, Flags32 flags, BaseControl *p_body = 0 );
    ~WidgetControl();

    void SetColor( gfx::Color color ); //from BasicControl
    //void SetFontSizes( const Vec2 &aFontSizes );
    
    inline void SetCaption( const char *caption ) { m_HeadBar.GetCaption().SetLabel(caption); }
    inline const char *GetCaption() const { return m_HeadBar.GetCaption().GetLabel(); }
    
    Flags32 GetWidgetFlags() const { return m_WidgetFlags; }
    
    //! \name IWidget BasicAPI implementation
    //@{
    void Close();
    bool AddFootButton( const char *name, D_NotifyCall notify_call, void *p_user_data );
    void SetRolled( bool b_rolled );
    bool IsRolled() const { return m_IsRolled; }
    void SetMinimized( bool b_minimized );
    bool IsMinimized() const { return m_IsMinimized; } 
    //@}    

    //! \name Widget keybindings
    //@{
    bool OnKeyReleased( EKey key, Vec2 p );
    //@}

    //! \name Parent and Subparts
    //@{
    Desktop *GetDesktop() const;
    inline HeadBar *GetHeadBar() { return &m_HeadBar; }
    inline BaseControl *GetBody() const { return m_pBody; }
    inline FootBar *GetFootBar() { return &m_FootBar; }
    //@}

    //! notifications
    //@{
    void SetOnCloseCall( D_OnCloseCall on_close_call ) { m_OnCloseCall = on_close_call; }
    //{}
    
protected:

    // Allow full access to subparts
    friend class HeadBar;
    friend class FootBar;

    Flags32 m_WidgetFlags;
    
    //!\name State
    //@{
    bool m_IsRolled;
    bool m_IsMinimized;
    //@}

    //\name Sub-controls
    //@{
    HeadBar m_HeadBar;
    FrameLayout m_BodyFrame;
    FootBar m_FootBar;    
    BaseControl *m_pBody;
    //@}
    
    D_OnCloseCall m_OnCloseCall;

public:
    //\name IWidget APIs implementation
    //@{
    class ImplBasicAPI: public IWidget::IBasicAPI
    {
    public:
        ImplBasicAPI( WidgetControl *p_wc ) : m_pWidgetControl(p_wc) {}
        void Close() { m_pWidgetControl->Close(); }
        void Minimize() { m_pWidgetControl->SetMinimized(true); }
        void Restore() { m_pWidgetControl->SetMinimized(false); }
        void Roll() { m_pWidgetControl->SetRolled(true); }
        void Unroll() { m_pWidgetControl->SetRolled(false); }
        bool AddFootButton( const char *name, D_NotifyCall notify_call, void *p_user_data )
            { return m_pWidgetControl->AddFootButton(name,notify_call,p_user_data); }
        void SetColor( float r, float g, float b, float a ) { m_pWidgetControl->SetColor( gfx::Color(r,g,b,a) ); }
    private:
        WidgetControl *m_pWidgetControl;
    };
    ImplBasicAPI m_BasicAPI;
    IWidget::IBasicAPI *GetBasicAPI() { return &m_BasicAPI; }
    IWidget::IPropertyTreeAPI *GetPropertyTreeAPI() { return m_pBody->GetAPI_PropertyTree(); }
    //@}
};

}} //namespace krm::krt::dbg::dsk

#endif //SFR_GUI_WIDGET_CONTROL_H
