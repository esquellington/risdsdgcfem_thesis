#include "Desktop.h"
#include "BaseControl.h"
#include <Safra/core/IView.h>

//#include "LabelControl.h"
#include "ButtonControl.h"
#include "LayoutControl.h"
#include "WidgetControl.h"
#include "TreeControl.h"
#include "PropertyTreeControl.h"
#include "MinimizedControl.h"

#include <Mal/GRandom.h> //TEMPORAL: to Restore at random pos...

#include <GL/gl.h> //TEMPORAL: to disable zbuffer

namespace sfr { namespace gui
{

Desktop::Desktop( IView *p_view )
: m_pView(p_view)
, m_bIsEnabled(true)
, m_bUpdatedSinceLastDraw(false)
, m_MessageTimeOut(0)
, m_MessageString("message")
{
    SetFlags( GetFlags().Disable( eDCF_AlwaysConsumeMouseEvents ) );
    SetColor( gfx::Color(1,1,1,0.1f) ); //\todo Draw perimeter only
    SetPosRel( Vec2(10,10) );
    Resize(300,300);
    //TEMPORAL: Test
    //CreateWidget( "Test", eDWF_Default, Vec2(80,80) );

    // Create Minimized widget
    m_pMinimizedControl = new MinimizedControl();
    WidgetControl *pWidget = new WidgetControl( gfx::Color(0.33,0.33,0.33,0.25),
                                                /*eDWF_CanMove |*/ eDWF_CanRoll | eDWF_CanAddButton,
                                                m_pMinimizedControl );
    pWidget->SetCaption( "Minimized" );
    AddChild( pWidget );
    pWidget->SetPosRel( Vec2(0,0) );
}

Desktop::~Desktop()
{
}

void Desktop::Message( const char *message, float duration )
{
    m_MessageTimeOut = duration;
    m_MessageString = message;
    //Update(0);
}

IWidget *Desktop::CreateWidget( const char *name, Flags32 flags, Vec2 pos )
{
    // Test Body
    ILayoutControl *pLayout( new VerticalLayout );
    {
        LabelControl *pLabel = new LabelControl( "Label" );
        pLayout->AddChild( pLabel );
        ButtonControl *pButton = new ButtonControl( "Button" );
        pLayout->AddChild( pButton );
    }
    pLayout->SetPosRel( pos );
    WidgetControl *pWidget( new WidgetControl( gfx::Color(1,0,0,0.5), flags, pLayout ) );
    pWidget->SetCaption( name );
    pWidget->GetFootBar()->AddButton( "(W)", std::bind1st( std::mem_fun(&Desktop::OnButtonCall_CreateWidget), this), 0 );
    pWidget->GetFootBar()->AddButton( "(T)", std::bind1st( std::mem_fun(&Desktop::OnButtonCall_CreateTree), this), 0 );
    pWidget->GetFootBar()->AddButton( "(PT)", std::bind1st( std::mem_fun(&Desktop::OnButtonCall_CreatePropertyTree), this), 0 );
    AddChild( pWidget );
    pWidget->SetPosRel( pos );
    //Update(0);
    return pWidget;
}

IWidget *Desktop::CreateWidget_Tree( const char *name, Flags32 flags, Vec2 pos )
{
    TreeControl *pTree( new TreeControl() );
    TreeItemId id000 = pTree->Add( "child000", 0 );
    TreeItemId id010 = pTree->Add( "child010", 0, id000 );
    pTree->Add( "child011", 0, id010 );
    pTree->Add( "child012_FDFDFDfdfdsfdsfdfdsfdsddfds", 0, id010 );
    TreeItemId id020 = pTree->Add( "child020", 0 );
    TreeItemId id030 = pTree->Add( "child030", 0 );
    pTree->Add( "child021", 0, id020 );
    pTree->Add( "child022", 0, id020 );
    pTree->Add( "child031", 0, id030 );
    pTree->Add( "child032", 0, id030 );
    pTree->Add( "child033", 0, id030 );
    WidgetControl *pWidget( new WidgetControl( gfx::Color(0,1,0,0.5), flags, pTree ) );
    pWidget->SetCaption( name );
    AddChild( pWidget );
    pWidget->SetPosRel( pos );
    //Update(0);
    return pWidget;
}

IWidget *Desktop::CreateWidget_PropertyTree( const char *name, Flags32 flags, Vec2 pos,
                                             PropertyIt property_tree_it,
                                             D_PropertyTreeSyncCall sync_call )
{
    SFR_ASSERT( property_tree_it.IsComplex() );
    SFR_ASSERT( eType_Property_Object == property_tree_it.GetType()
                || eType_Property_Group == property_tree_it.GetType() );
    PropertyTreeControl *pPT = new PropertyTreeControl( property_tree_it, sync_call );
    WidgetControl *pWidget = new WidgetControl( gfx::Color(0,0,1,0.5), flags, pPT );
    pWidget->SetCaption( (std::string(name) + " " + property_tree_it.GetName()).c_str() );
    AddChild( pWidget );
    pWidget->SetPosRel( pos );
    //Update(0);
    return pWidget;
}

void Desktop::Resize( int w, int h )
{
    SetSizes( Vec2(w,h) - Vec2(20,20) );
    //Update(0);
}

bool Desktop::Draw()
{
    /*TEMPORAL: MUST ensure that Update() is called at least once each
      frame to guarantee that Layout will be up-to-date... CONSIDER
      splitting Update() into Step(dt) and PreDraw(), where PreDraw()
      would update the layout, etc, and Step(dt) would do strictly
      time-related stuff.
    */
    if( !m_bUpdatedSinceLastDraw ) Update(0);

    //\todo We should set camera-specific state or kha::Phase for
    // the desktop, with no Zbuffer, transparency, etc... and
    // sned draw primitives here but draw according to kha
    // target phase ordering...
    m_pView->GetCamera()->PrepareToShootViewport();
    glDisable(GL_DEPTH_TEST);
    BaseControl::Draw( m_pView->GetRenderer() );
    //m_pView->GetRenderer()->DrawLabel( Vec3(GetPosAbs().x(), GetPosAbs().y(), 0), "Desktop!!", gfx::Style(gfx::Color(0,1,0,1)) );
    if( m_MessageTimeOut > 0 )
        m_pView->GetRenderer()->DrawLabel( Vec3(GetPosAbs().x(), GetPosAbs().y() + GetSizes().y(), 0),
                                           m_MessageString.c_str(), gfx::Style(gfx::Color(0.5,0.5,0,1)) );
    //TEMPORAL: Restore old state
    glEnable(GL_DEPTH_TEST);
    m_pView->GetCamera()->PrepareToShoot();

    m_bUpdatedSinceLastDraw = false;

    return true;
}

bool Desktop::Update( float dt )
{
    m_bUpdatedSinceLastDraw = true;
    m_MessageTimeOut -= dt;
    return BaseControl::Update(dt);
}

void Desktop::SetBackgroundColor( float r, float g, float b, float a )
{
    SetColor( gfx::Color(r,g,b,a) );
}

void Desktop::Minimize( WidgetControl *p_widget )
{
    /*
    SFR_LOG( "Minimizing widget %s", p_widget->GetCaption() );
    {
        SFR_LOG( "Desktop children:" );
        BaseControl *p_child( m_pFirstChild );
        while( 0 != p_child  )
        {
            SFR_LOG( "%llx", (machine_uint_type)p_child );
            p_child = p_child->GetNextSibling();
        }
    }
    */
    if( RemoveChild( p_widget ) )
    {
        /*
        SFR_LOG( "Removed widget %s from Desktop", p_widget->GetCaption() );
        {
            SFR_LOG( "Desktop children:" );
            BaseControl *p_child( m_pFirstChild );
            while( 0 != p_child  )
            {
                SFR_LOG( "%llx", (machine_uint_type)p_child );
                p_child = p_child->GetNextSibling();
            }
        }
        */
        m_pMinimizedControl->AddChild( p_widget );
        /*
        SFR_LOG( "Added %s to MinimizedControl", p_widget->GetCaption() );
        {
            SFR_LOG( "Minimized children:" );
            BaseControl *p_child( m_pMinimizedControl->GetFirstChild() );
            while( 0 != p_child  )
            {
                SFR_LOG( "%llx", (machine_uint_type)p_child );
                p_child = p_child->GetNextSibling();
            }
        }
        */
    }
    else SFR_LOG_ERROR( "Cannot minimize non-child widget %s", p_widget->GetCaption() );
}

void Desktop::Restore( WidgetControl *p_widget )
{
    //SFR_LOG( "Restoring widget %s", p_widget->GetCaption() );
    if( m_pMinimizedControl->RemoveChild( p_widget ) )
    {
        AddChild( p_widget );
        /*\todo We restore it in a random pos. Ideally we should
          either save previous pos or find a non-overlapping one
          automatically
        */
        p_widget->SetPosRel( mal::RandomV( 0.1f*GetSizes(), 0.9f*p_widget->GetSizes() ) );
    }
    else SFR_LOG_ERROR( "Cannot restore non-child widget %s", p_widget->GetCaption() );
}

bool Desktop::OnMouseButtonPressed( EMouseButton mb, int x, int y )
{
    //SFR_LOG( "MB at (%d,%d)", x, y );
    if( mb != eMB_Left ) return false;
    else return BaseControl::OnMouseButtonPressed(mb,Vec2(x,y));
}
bool Desktop::OnMouseButtonReleased( EMouseButton mb, int x, int y )
{
    if( mb != eMB_Left ) return false;
    else return BaseControl::OnMouseButtonReleased(mb,Vec2(x,y));
}
bool Desktop::OnMouseMotion( int x, int y ) { return BaseControl::OnMouseMotion(Vec2(x,y)); }
bool Desktop::OnKeyPressed( EKey key, int x, int y ) { return BaseControl::OnKeyPressed(key,Vec2(x,y)); }
bool Desktop::OnKeyReleased( EKey key, int x, int y ) { return BaseControl::OnKeyReleased(key,Vec2(x,y)); }


//TEMPORAL: TESTING---------------------------------------
void Desktop::OnButtonCall_CreateWidget( void *p_user_data )
{
    IWidget *pWidget = CreateWidget( "I'm Alive", eDWF_Default, Vec2(110,110) );
    pWidget->GetBasicAPI()->Close();
}
void Desktop::OnButtonCall_CreateTree( void *p_user_data )
{
    CreateWidget_Tree( "I'm Alive", eDWF_Default, Vec2(110,110) );
}

void Desktop::OnButtonCall_CreatePropertyTree( void *p_user_data )
{
    util::ItemStream *pIS( new util::ItemStream(1024,1024) );
    pIS->BeginComplex( "[PropertyObject]", eType_Property_Object );
    {
        //pIS->Write( "int32", int32(-1) );
        pIS->Write( "NIR_int32", Property_NIR_int32(10,-5,5) );
        pIS->Write( "NIR_uint32", Property_NIR_uint32(10,0,20) );
        pIS->Write( "NIR_float32", Property_NIR_float32(1.0f,0.0f,2.0f) );
        pIS->Write( "NIR_float64", Property_NIR_float64(1.0,0.0,2.0) );
        pIS->BeginComplex( "<PropertyGroup>", eType_Property_Group );
        {
            pIS->Write( "is_closed", true );
        }
        pIS->EndComplex();
    }
    pIS->EndComplex();
    CreateWidget_PropertyTree( "I'm Alive", eDWF_Default, Vec2(110,110), pIS->BeginRW() );
}


} } // namespace sfr::gui
