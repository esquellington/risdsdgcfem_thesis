#include "TreeControl.h"
#include "LayoutControl.h"
#include <stdio.h> //TEMPORAL: for sprintf...

namespace sfr { namespace gui {

/*! Collapsable treeview item
  - Expand/Collapse
  - User-defined buttons
*/
class TreeItemControl: public FrameLayout
{
public:
    TreeItemControl( const char *name, void *p_user_data, int margin_h, int margin_v )
    : m_pUserData(p_user_data)
    , m_bIsRolled(true)
    , m_buttonToggleRoll(name)
    , m_NumSubItems(0)
        {
            m_NameLength = strlen(name);
            SFR_ASSERT( m_NameLength < LabelControl::cMaxLength );
            strcpy( &m_Name[0], name );

            AddChild(&m_Layout);
            //---- Add buttons row
            if( m_NameLength > 0 )
                m_Layout.AddChild( &m_ButtonsLayout );
            else
            {
                //Set permanently to unrolled, there'll be no way to roll it
                m_bIsRolled = false;
                m_Layout.AddChild( &m_SubItemsLayout );
            }
            m_Layout.SetMinSeparation( margin_v );
            // Rolled but with no childs on init
            m_buttonToggleRoll.SetFieldLength( 16 ); //\todo This should be adapted globaly to avoid label cropping
            m_buttonToggleRoll.SetColor( gfx::Color(0,0,0,0) );
            m_buttonToggleRoll.SetFontColor( gfx::Color(1,1,1,1) );
            m_buttonToggleRoll.SetOnButtonCall( std::bind1st( std::mem_fun(&TreeItemControl::OBC_ToggleRoll), this) );
            m_ButtonsLayout.AddChild( &m_buttonToggleRoll );
            //---- Var layout, initially rolled
            m_SubItemsLayout.SetMargins( Vec2(margin_h,0), Vec2(0,0) );
            m_SubItemsLayout.SetMinSeparation( margin_v );
        }
    ~TreeItemControl()
        {
            m_ButtonsLayout.RemoveChild( &m_buttonToggleRoll );
            m_ButtonsLayout.RemoveAndDeleteChildren();
            m_SubItemsLayout.RemoveAndDeleteChildren();
        }

    //TEMPORAL: For lazy recompute monitoring
    //bool Update( float dt ) { UpdateSubItemCountText(); FrameLayout::Update(dt); }

    TreeItemId Add( const char *name, void *p_user_data, int margin_h, int margin_v )
        {
            TreeItemControl *pItem = new TreeItemControl( name, p_user_data, margin_h, margin_v );
            // If first child, enable rolling
            if( 0 == m_SubItemsLayout.GetFirstChild() )
            {
                m_buttonToggleRoll.SetColor( gfx::Color(1,1,1,1) );
                m_buttonToggleRoll.SetFontColor( gfx::Color(0,0,0,1) );
            }
            m_SubItemsLayout.AddChild( pItem );

            m_NumSubItems++;
            UpdateSubItemCountText();

            return reinterpret_cast<TreeItemId>(pItem);
        }

    TreeItemId FindByName( const char *name ) const
        {
            SFR_LOG_ERROR("TreeItemId::FindByName() NOT YET IMPLEMENTED");
            return TreeItemId(TreeControl::cInvalidItemId);
        }

    TreeItemId FindByUserData( void *p_user_data ) const
        {
            // DFS search
            for( BaseControl *it_subitem = m_SubItemsLayout.GetFirstChild();
                 it_subitem != 0;
                 it_subitem = it_subitem->GetNextSibling() )
            {
                TreeItemControl *pItem( static_cast<TreeItemControl*>(it_subitem) );
                if( pItem->m_pUserData == p_user_data )
                    return TreeItemId(pItem);
                else
                {
                    TreeItemId subitem_id = pItem->FindByUserData(p_user_data);
                    if( TreeItemId(TreeControl::cInvalidItemId) != subitem_id )
                        return subitem_id;
                }
            }
            return TreeItemId(TreeControl::cInvalidItemId);
        }

    void Remove( TreeItemId id )
        {
            SFR_LOG_ERROR("TreeItemId::Remove() NOT YET IMPLEMENTED");
            /*
              m_NumSubItems--;
              UpdateSubItemCountText();
              // If no more children, disable rolling
              if( 0 == m_SubItemsLayout.GetFirstChild() )
              {
              m_buttonToggleRoll.SetColor( gfx::Color(0,0,0,0) );
              m_buttonToggleRoll.SetFontColor( gfx::Color(1,1,1,1) );
              }
            */
        }

    void Clear()
        {
            m_NumSubItems = 0;
            m_SubItemsLayout.RemoveAndDeleteChildren();
            m_buttonToggleRoll.SetColor( gfx::Color(1,1,1,1) );
            m_buttonToggleRoll.SetFontColor( gfx::Color(0,0,0,1) );
            UpdateSubItemCountText();
        }

    void AddButton( const char *button_name, ButtonControl::D_OnButtonCall on_button_call )
        {
            ButtonControl *pButton = new ButtonControl( button_name, m_pUserData );
            pButton->SetOnButtonCall(on_button_call);
            m_ButtonsLayout.AddChild(pButton);
        }
    void AddControl( BaseControl *p_control )
        {
            m_ButtonsLayout.AddChild(p_control);
        }

    void OBC_ToggleRoll( void* p_user_data )
        {
            // Only if it has sub-items
            if( 0 != m_SubItemsLayout.GetFirstChild() )
            {
                if( m_bIsRolled )
                {
                    m_Layout.AddChild( &m_SubItemsLayout );
                    m_buttonToggleRoll.SetColor( gfx::Color(0,0,0,1) );
                    m_buttonToggleRoll.SetFontColor( gfx::Color(1,1,1,1) );
                }
                else
                {
                    m_Layout.RemoveChild( &m_SubItemsLayout );
                    m_buttonToggleRoll.SetColor( gfx::Color(1,1,1,1) );
                    m_buttonToggleRoll.SetFontColor( gfx::Color(0,0,0,1) );
                }
                m_bIsRolled = !m_bIsRolled;
            }
        }

    void OBC_RecomputeLayout( void* p_user_data )
        {
            m_ButtonsLayout.RecomputeLayout();
            RecomputeLayout();
        }

    inline bool IsRolled() const { return m_bIsRolled; }

    inline TreeItemId GetParentItem() const { return TreeItemId( static_cast<TreeItemControl*>( GetParent()->GetParent() ) ); } //item.this -> parent.m_SubItemsLayout -> parent
    inline TreeItemId GetFirstChildItem() const { return TreeItemId( static_cast<TreeItemControl*>( m_SubItemsLayout.GetFirstChild() ) ); }
    inline TreeItemId GetNextSiblingItem() const { return TreeItemId( static_cast<TreeItemControl*>( GetNextSibling() ) ); }

private:
    void UpdateSubItemCountText()
        {
            if( false )//m_NumSubItems > 0 ) //TEMPORAL: This should be optional... but by default it's ugly
            {
                char str[128];
                sprintf( str,"%s (%2d)", m_Name, m_NumSubItems );
                //sprintf( str,"%s#%10d", m_Name, m_CountRL ); //TEMPORAL: Testing lazy update
                m_buttonToggleRoll.SetLabel(str);
            }
            else
                m_buttonToggleRoll.SetLabel(m_Name);
        }

protected:
    void *m_pUserData;
    bool m_bIsRolled;
    VerticalLayout m_Layout;
    HorizontalLayout m_ButtonsLayout;
    VerticalLayout m_SubItemsLayout;
    ButtonControl m_buttonToggleRoll;

    int m_NumSubItems;
    int m_NameLength;
    char m_Name[LabelControl::cMaxLength];
};

//--------------------------------------------------------
//---- TreeControl Implementation
//--------------------------------------------------------

TreeControl::TreeControl( unsigned int hint_max_items, int margin_h, int margin_v )
: m_pRoot(0)
, m_MarginH(margin_h)
, m_MarginV(margin_v)
{
    // Create root rolled, by default
    //m_pRoot = new TreeItemControl( "<root>", 0, m_MarginH, m_MarginV );
    m_pRoot = new TreeItemControl( "", 0, m_MarginH, m_MarginV );
    AddChild(m_pRoot);
}

TreeControl::~TreeControl()
{
    RemoveAndDeleteChildren();
}

TreeItemId TreeControl::Add( const char *name, void *p_user_data, TreeItemId parent_id )
{
    if( TreeItemId(TreeControl::cRootItemId) == parent_id ) return m_pRoot->Add( name, p_user_data, m_MarginH, m_MarginV );
    else return reinterpret_cast<TreeItemControl*>(parent_id)->Add( name, p_user_data, m_MarginH, m_MarginV );
}

TreeItemId TreeControl::FindByName( const char *name, TreeItemId parent_id ) const
{
    if( TreeItemId(TreeControl::cRootItemId) == parent_id ) return m_pRoot->FindByName( name );
    else return reinterpret_cast<TreeItemControl*>(parent_id)->FindByName( name );
}

TreeItemId TreeControl::FindByUserData( void *p_user_data, TreeItemId parent_id ) const
{
    if( TreeItemId(TreeControl::cRootItemId) == parent_id ) return m_pRoot->FindByUserData( p_user_data );
    else return reinterpret_cast<TreeItemControl*>(parent_id)->FindByUserData( p_user_data );
}

void TreeControl::Remove( TreeItemId id )
{
    m_pRoot->Remove(id);
}

void TreeControl::Clear( TreeItemId id )
{
    if( TreeItemId(TreeControl::cRootItemId) == id ) m_pRoot->Clear();
    else reinterpret_cast<TreeItemControl*>(id)->Clear();
}

void TreeControl::SetRolled( TreeItemId id, bool b_rolled )
{
    TreeItemControl *pItem(0);
    if( TreeItemId(TreeControl::cRootItemId) == id ) pItem = m_pRoot;
    else pItem = reinterpret_cast<TreeItemControl*>(id);
    // Call existing toggle_roll if required (we avoid a SetRolled() method with this...)
    if( pItem->IsRolled() != b_rolled ) pItem->OBC_ToggleRoll(0);
}

void TreeControl::AddButton( TreeItemId id, const char *button_name, ButtonControl::D_OnButtonCall on_button_call )
{
    SFR_ASSERT( TreeItemId(TreeControl::cRootItemId) != id ); //NO buttons in Root
    reinterpret_cast<TreeItemControl*>(id)->AddButton( button_name, on_button_call );
}

void TreeControl::AddControl( TreeItemId id, BaseControl *p_control )
{
    SFR_ASSERT( TreeItemId(TreeControl::cRootItemId) != id ); //NO controls in Root
    reinterpret_cast<TreeItemControl*>(id)->AddControl( p_control );
}

TreeItemId TreeControl::GetParentItem( TreeItemId id ) const
{
    if( TreeItemId(TreeControl::cRootItemId) == id ) return TreeItemId(0);
    else return reinterpret_cast<TreeItemControl*>(id)->GetParentItem();
}
TreeItemId TreeControl::GetFirstChildItem( TreeItemId id ) const
{
    if( TreeItemId(TreeControl::cRootItemId) == id ) return m_pRoot->GetFirstChildItem();
    else return reinterpret_cast<TreeItemControl*>(id)->GetFirstChildItem();
}
TreeItemId TreeControl::GetNextSiblingItem( TreeItemId id ) const
{
    if( TreeItemId(TreeControl::cRootItemId) == id ) return TreeItemId(0);
    else return reinterpret_cast<TreeItemControl*>(id)->GetNextSiblingItem();
}

}} //namespace sfr::gui
