#include "BaseControl.h"

namespace sfr { namespace gui
{

BaseControl::BaseControl()
: m_Flags(eDCF_Default)
, m_PosRel(0,0)
, m_Sizes(10,10)
, m_Color( 0.5, 0.5, 0.5, 1 )
, m_Time(0)
, m_StatusFlags(eStatus_Default)
, m_TimeStamp(0)
, m_pParent(0)
, m_pNextSibling(0)
, m_pFirstChild(0)
, m_pFocusedChild(0)
//, m_PDS(0)
, m_ParentTimeStamp(0)
{}

BaseControl::~BaseControl()
{
    if( m_Flags.Test( eDCF_AutoRemoveAndDeleteChildren ) )
        RemoveAndDeleteChildren();
}

//---- Control Nesting/Hierarchy
bool BaseControl::AddChild( BaseControl *p_ctrl )
{
    //SFR_LOG( "BaseControl::AddChild( %llx )", (machine_uint_type)p_ctrl );
    SFR_ASSERT( 0 != p_ctrl && 0 == p_ctrl->GetParent() );
    m_TimeStamp++;
    if( 0 == m_pFirstChild )
        m_pFirstChild = p_ctrl;
    else
    {
        BaseControl *pChild(m_pFirstChild);
        while( 0 != pChild->GetNextSibling() )
            pChild = pChild->GetNextSibling();
        pChild->SetNextSibling(p_ctrl);
    }
    p_ctrl->SetParent(this);
    return true;
}

bool BaseControl::RemoveChild( BaseControl *p_ctrl )
{
    //SFR_LOG( "BaseControl::RemoveChild( %llx )", (machine_uint_type)p_ctrl );
    SFR_ASSERT( 0 != p_ctrl && this == p_ctrl->GetParent() );
    m_TimeStamp++;

    // Invalidate focused child if removed
    if( m_pFocusedChild == p_ctrl )
        m_pFocusedChild = 0;

    // Remove child (first or interior)
    if( m_pFirstChild == p_ctrl )
    {
        m_pFirstChild = m_pFirstChild->GetNextSibling();
        p_ctrl->UnlinkParentAndSibling();
        return true;
    }
    else
    {
        BaseControl *pChild0( m_pFirstChild );
        while( 0 != pChild0->GetNextSibling() )
        {
            BaseControl *pChild1( pChild0->GetNextSibling() );
            if( 0 != pChild1 && pChild1 == p_ctrl )
            {
                pChild0->SetNextSibling( pChild1->GetNextSibling() );
                p_ctrl->UnlinkParentAndSibling();
                return true;
            }
            else
                pChild0 = pChild1;
        }
    }
    return false;
}

void BaseControl::RemoveAndDeleteChildren()
{
    m_TimeStamp++;
    BaseControl *pChild(m_pFirstChild);
    while( pChild != 0 )
    {
        BaseControl *pTmp = pChild;
        pChild = pChild->GetNextSibling();
        delete pTmp;
    }
    m_pFirstChild = 0;
}

//---- Input Handling

//! Send to focused child if exists. Process globally if nonexistent or rejected.
bool BaseControl::OnKeyPressed( EKey key, Vec2 p )
{
    if( !IsActive() )
        return false;
    else if( m_pFocusedChild && m_pFirstChild->OnKeyPressed(key,p) )
        return true;
    else if (m_Flags.Test( eDCF_AlwaysConsumeKeyboardEvents ) && IsPointInside(p) )
        return true;
    else
        return false;
}

bool BaseControl::OnKeyReleased( EKey key, Vec2 p )
{
    if( !IsActive() )
        return false;
    else if( m_pFocusedChild && m_pFirstChild->OnKeyReleased(key,p) )
        return true;
    else if (m_Flags.Test( eDCF_AlwaysConsumeKeyboardEvents ) && IsPointInside(p) )
        return true;
    else
        return false;
}

bool BaseControl::OnMouseButtonPressed( EMouseButton mb, Vec2 p )
{
    if( !IsActive() ) return false;
    // Send mouse event to focused child
    if( m_pFocusedChild && m_pFocusedChild->OnMouseButtonPressed(mb,p) )
        return true;
    else
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            if( pChild->IsPointInside(p) )
            {
                // Check if focus must change
                if( pChild != m_pFocusedChild )
                {
                    if( m_pFocusedChild ) m_pFocusedChild->OnFocusLost();
                    if( pChild->OnFocusGained() ) m_pFocusedChild = pChild;
                }
                // Send event to control, regardless of focus
                if( pChild->OnMouseButtonPressed(mb,p) )
                    return true;
            }
    // Process here if not accepted by any child
    if( m_Flags.Test(eDCF_AlwaysConsumeMouseEvents) && IsPointInside(p) )
        return true;
    return false;
}

bool BaseControl::OnMouseButtonReleased( EMouseButton mb, Vec2 p )
{
    if( !IsActive() ) return false;
    // Send mouse event to focused child or, if none, to all children hit
    if( m_pFocusedChild && m_pFocusedChild->OnMouseButtonReleased(mb,p) )
        return true;
    else
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            if( pChild->IsPointInside(p) && pChild->OnMouseButtonReleased(mb,p) )
                return true;
    // Process here if not accepted by any child
    if( m_Flags.Test(eDCF_AlwaysConsumeMouseEvents) && IsPointInside(p) )
        return true;
    return false;
}

bool BaseControl::OnMouseMotion( Vec2 p )
{
    if( !IsActive() ) return false;
    // Send mouse event to focused child or, if none, to all children hit
    if( m_pFocusedChild && m_pFocusedChild->OnMouseMotion(p) )
        return true;
    else
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            if( pChild->IsPointInside(p) && pChild->OnMouseMotion(p) )
                return true;
    // Process here if not accepted by any child
    if( m_Flags.Test(eDCF_AlwaysConsumeMouseEvents) && IsPointInside(p) )
        return true;
    return false;
}

//---- Update-cycle methods
bool BaseControl::Draw( gfx::IRenderer *p_renderer )
{
    if( !IsVisible() ) return false;

    //---- Draw myself
    Vec2 pos_abs( GetPosAbs() );

    // Draw if not transparent
    if( m_Color[3] != 0 )
    {
        //Simple AABB p_renderer->DrawAABB2( pos_abs, pos_abs+GetSizes(), gfx::Style( m_Color, 1.0f, gfx::Style::eSolid ) );
        // Optional 45º corners
        Vec2 sizes( GetSizes() );
        Vec2 pos_ul( pos_abs
                     + Vec2( m_Flags.Test(eDCF_Corner_UL) ? sizes[1] : 0, 0 ) );
        Vec2 pos_ur( pos_abs + Vec2(sizes[0],0)
                     - Vec2( m_Flags.Test(eDCF_Corner_UR) ? sizes[1] : 0, 0 ) );
        Vec2 pos_bl( pos_abs + Vec2(0,sizes[1])
                     + Vec2( m_Flags.Test(eDCF_Corner_BL) ? sizes[1] : 0, 0 ) );
        Vec2 pos_br( pos_abs + sizes
                     - Vec2( m_Flags.Test(eDCF_Corner_BR) ? sizes[1] : 0, 0 ) );
        p_renderer->DrawQuad2( pos_ul, pos_bl, pos_br, pos_ur,
                               gfx::Style( m_Color, 1.0f, gfx::Style::eSolid ) );
    }

    // Draw focus
    if( HasFocus() && m_Flags.Test(eDCF_HighlightOnFocus) )
    {
        Vec2 half_offset(2,2);
        p_renderer->DrawAABB2( pos_abs-half_offset, pos_abs+GetSizes()+half_offset,
                               gfx::Style( gfx::Color(1,1,1,1), 1.0f, gfx::Style::eWire ) );
    }

    //---- Draw children
    for( BaseControl *pChild = m_pFirstChild; 0 != pChild; pChild = pChild->GetNextSibling() )
        pChild->Draw(p_renderer);

    return true;
}

bool BaseControl::Update( float dt )
{
    if( !IsActive() ) return false;
    // Update myself
    m_Time += dt;
    // Update children
    for( BaseControl *pChild = m_pFirstChild; 0 != pChild; pChild = pChild->GetNextSibling() )
        pChild->Update(dt);
    return true;
}

} } //namespace sfr::gui
