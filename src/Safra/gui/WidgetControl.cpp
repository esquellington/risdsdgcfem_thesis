#include "WidgetControl.h"
#include "Desktop.h"

namespace sfr { namespace gui
{

//------------------------------------------------------------------
//---- WidgetControl Implementation
//------------------------------------------------------------------

WidgetControl::WidgetControl( gfx::Color color, Flags32 flags, BaseControl *p_body )
: m_WidgetFlags(flags)
, m_IsRolled(false)
, m_IsMinimized(false)
, m_HeadBar(this)
, m_FootBar(this)
, m_pBody(p_body)
, m_BasicAPI(this)
{
    SFR_ASSERT( 0 != p_body );

    // Add children in order Head,Body,Foot
    AddChild( &m_HeadBar );
    AddChild( &m_BodyFrame );
    m_BodyFrame.AddChild( m_pBody );
    AddChild( &m_FootBar );
    m_FootBar.SetSizes(Vec2(10,10));
    m_FootBar.SetFlags( m_FootBar.GetFlags().Enable( eDCF_Corner_BR ) );
    SetColor(color);
    SetMinSeparation(0);
}

WidgetControl::~WidgetControl()
{
    delete m_pBody;
}

void WidgetControl::SetColor( gfx::Color color )
{
    // Widget layout and Body are transparent
    BaseControl::SetColor( gfx::Color(1,1,1,0) );
    m_pBody->SetColor( gfx::Color(1,1,1,0) );
    // Body Frame has the given color
    m_BodyFrame.SetColor( color );
    // Head,Foot have body-derived colors
    gfx::Color lHeadColor( 0.5f*color[0],
                           0.5f*color[1],
                           0.5f*color[2],
                           0.5f );
    gfx::Color lFootColor( 0.75f*color[0],
                           0.75f*color[1],
                           0.75f*color[2],
                           0.5f );
    m_HeadBar.SetColor( lHeadColor  );
    m_FootBar.SetColor( lFootColor );
}

/*
void WidgetControl::SetFontSizes( const TVec2 &aFontSizes )
{
    m_HeadBar.GetCaption().SetFontSizes(aFontSizes);
    //m_pBody->SetFontSizes(aFontSizes);
    //m_FootBar.SetFontSizes(aFontSizes);
}
*/

Desktop *WidgetControl::GetDesktop() const
{
    /* Direct Parent is NOT always desktop directly, as the Widget may
       be embedded in some other control (eg: MinimizedControl), thus
       we search for the root control, which will always be te
       desktop.
     */
    BaseControl *pParent( GetParent() );
    while( 0 != pParent->GetParent() ) pParent = pParent->GetParent();
    return static_cast<Desktop*>( pParent );
}

//---- Widget actions
void WidgetControl::SetRolled( bool b_rolled )
{
    if( m_IsRolled != b_rolled )
    {
        //\todo MUST NOTIFY HEADBAR TO UPDATE +/- LABEL!!
        m_pBody->SetActive(!b_rolled);
        if( b_rolled )
            RemoveChild(&m_BodyFrame);
            //RemoveChild(m_pBody);
        else
        {
            // Preserve order Head->Body->Foot in layout
            RemoveChild(&m_FootBar);
            //AddChild(m_pBody);
            AddChild(&m_BodyFrame);
            AddChild(&m_FootBar);
        }
        m_IsRolled = b_rolled;
        m_HeadBar.NotifyRolledChange();
    }
}

void WidgetControl::SetMinimized( bool b_minimized )
{
    if( m_IsMinimized != b_minimized )
    {
        if( b_minimized ) { SetRolled(true); GetDesktop()->Minimize( this ); }
        else { GetDesktop()->Restore( this ); SetRolled(false); }
        m_IsMinimized = b_minimized;
        m_HeadBar.NotifyMinimizedChange();
    }
}

//---- IWidget implementation
void WidgetControl::Close()
{
    if( !m_OnCloseCall.empty() ) m_OnCloseCall(this);
    if( 0 != GetParent() ) GetParent()->RemoveChild(this);
    else SFR_LOG_ERROR("WidgetControl::Close() on a widget with no parent. This may only happen if you call Destroy() on a Minimized widget." );
    delete this;
}

bool WidgetControl::AddFootButton( const char *name, D_NotifyCall notify_call, void *p_user_data )
{
    if( !m_WidgetFlags.Test(eDWF_CanAddButton) ) return false;
    m_FootBar.AddButton( name, notify_call, p_user_data );
    return true;
}

bool WidgetControl::OnKeyReleased( EKey key, Vec2 p )
{
    if( key == gui::eKey_Escape && m_WidgetFlags.Test(eDWF_CanClose) )
    {
        Close();
        return true;
    }
    else
        return VerticalLayout::OnKeyReleased(key,p);
}

//------------------------------------------------------------------
//---- HeadBar Implementation
//------------------------------------------------------------------

HeadBar::HeadBar( WidgetControl *p_widget )
: m_pWidget(p_widget)
, m_Caption(p_widget)
, m_buttonToggleRoll("-")
, m_buttonToggleMinimize("_")
, m_buttonClose("X")
{
    if( m_pWidget->GetWidgetFlags().Test(eDWF_CanRoll) )
    {
        AddChild(&m_buttonToggleRoll);
        m_buttonToggleRoll.SetOnButtonCall( std::bind1st( std::mem_fun(&HeadBar::ToggleRoll), this) );
    }

    AddChild(&m_Caption);
    m_Caption.SetFontColor(gfx::Color(1,1,1,1));
    m_Caption.SetFontSizes(Vec2(8,13));
    m_Caption.SetLabel("Unnamed");

    if( m_pWidget->GetWidgetFlags().Test(eDWF_CanMinimize) )
    {
        AddChild(&m_buttonToggleMinimize);
        m_buttonToggleMinimize.SetOnButtonCall( std::bind1st( std::mem_fun(&HeadBar::ToggleMinimize), this) );
    }

    if( m_pWidget->GetWidgetFlags().Test(eDWF_CanClose) )
    {
        AddChild(&m_buttonClose);
        m_buttonClose.SetOnButtonCall( std::bind1st( std::mem_fun(&HeadBar::Close), this) );
    }
}

HeadBar::~HeadBar()
{
}

//---- Action button delegates
void HeadBar::ToggleRoll( void* /*aUserData*/ )
{
    m_pWidget->SetRolled(!m_pWidget->IsRolled());
}
void HeadBar::ToggleMinimize( void* /*aUserData*/ )
{
    m_pWidget->SetMinimized( !m_pWidget->IsMinimized() );
}
void HeadBar::Close( void* /*aUserData*/ ) { m_pWidget->GetBasicAPI()->Close(); }

void HeadBar::NotifyRolledChange()
{
    m_buttonToggleRoll.SetLabel( m_pWidget->IsRolled() ? "+" : "-" );
}
void HeadBar::NotifyMinimizedChange()
{
    m_buttonToggleMinimize.SetLabel( m_pWidget->IsMinimized() ? "^" : "_" );
}

//------------------------------------------------------------------
//---- WidgetControl::Caption Implementation
//------------------------------------------------------------------

Caption::Caption( WidgetControl *p_widget )
: m_pWidget(p_widget)
, m_InteractionState(eIS_Standby)
{
}

Caption::~Caption()
{
}

bool Caption::OnMouseButtonPressed( EMouseButton mb, Vec2 p )
{
    //\todo Ignore if minimized... could Restore also...
    if( m_pWidget->IsMinimized() ) return false;

    if( m_pWidget->GetWidgetFlags().Test(eDWF_CanMove) )
    {
        if( !IsPointInside(p) )
            return false;

        switch(m_InteractionState)
        {
        case eIS_Standby: m_InteractionState = eIS_Moving; m_InteractionPos0=p; return true; break;
        case eIS_Moving: SFR_LOG("WidgetControl::OnMouseButtonPressed() ignored while eMoving"); return true; break;
        default: return false; break;
        }
    }
    else
        return false;
}

bool Caption::OnMouseButtonReleased( EMouseButton mb, Vec2 p )
{
    if( !HasFocus() ) return false;

    switch(m_InteractionState)
    {
    case eIS_Standby: return false; break;
    case eIS_Moving: m_InteractionState = eIS_Standby; return true; break;
    default: return false; break;
    }
}

/*! Caption movement is clamped to Desktop, but widget body can extend outside of it. */
bool Caption::OnMouseMotion( Vec2 p )
{
    // \todo Could highlight color OnMouseEnter()/OnMouseLeave()
    switch(m_InteractionState)
    {
    case eIS_Standby: return false; break;
    case eIS_Moving:
        {
            Vec2 new_pos_rel( m_pWidget->GetPosRel() + (p-m_InteractionPos0) );
            new_pos_rel = mal::Clamp( new_pos_rel,
                                      Vec2(0,0),
                                      m_pWidget->GetDesktop()->GetSizes() - GetParent()->GetSizes() );
            m_pWidget->SetPosRel( new_pos_rel );
            m_InteractionPos0 = p;
            return true;
        }
        break;
    default: return false; break;
    };
}

//------------------------------------------------------------------
//---- FootBar Implementation
//------------------------------------------------------------------

FootBar::FootBar( WidgetControl *p_widget )
: m_pWidget(p_widget)
{
}

FootBar::~FootBar()
{
    for(std::vector<ButtonControl *>::iterator it_button = m_vecButtons.begin(); it_button != m_vecButtons.end(); ++it_button)
    {
        delete *it_button;
    }
    m_vecButtons.clear();
}

bool FootBar::AddButton( const char *name, ButtonControl::D_OnButtonCall on_button_call, void *p_user_data )
{
    if( m_pWidget->GetWidgetFlags().Test(eDWF_CanAddButton) )
    {
        ButtonControl *pButton = new ButtonControl( name, p_user_data );
        pButton->SetOnButtonCall( on_button_call );
        AddChild(pButton);
        m_vecButtons.push_back(pButton);
        return true;
    }
    else
        return false;
}

bool FootBar::RemoveButton( const char* name )
{
    if( m_pWidget->GetWidgetFlags().Test(eDWF_CanAddButton) )
    {
        SFR_LOG_ERROR("FootBar::RemoveButton() NOT YET IMPLEMENTED");
        return false;
    }
    else
        return false;
}

}} //namespace sfr::gui
