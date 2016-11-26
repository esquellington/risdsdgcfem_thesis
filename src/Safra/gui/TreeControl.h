#ifndef SFR_GUI_TREE_CONTROL_H
#define SFR_GUI_TREE_CONTROL_H

#include "LayoutControl.h"
#include "ButtonControl.h"

namespace sfr { namespace gui {

class TreeItemControl;
typedef void* TreeItemId;

class TreeControl: public FrameLayout
{
public:
    enum EConstants { cRootItemId = 0, cInvalidItemId = 0xFFFFFFFF };

public:
    TreeControl( unsigned int hint_max_items = 10, int margin_h = 10, int margin_v = 5 );
    ~TreeControl();

    TreeItemId Add( const char *name, void *p_user_data, TreeItemId parent_id = TreeItemId(cRootItemId) );
    TreeItemId FindByName( const char *name, TreeItemId parent_id = TreeItemId(cRootItemId) ) const;
    TreeItemId FindByUserData( void *p_user_data, TreeItemId parent_id = TreeItemId(cRootItemId) ) const;

    void Remove( TreeItemId id );
    void Clear( TreeItemId id );
    void SetRolled( TreeItemId id, bool b_rolled );

    void AddButton( TreeItemId id, const char *button_name, ButtonControl::D_OnButtonCall on_button_call );
    void AddControl( TreeItemId id, BaseControl *p_control ); //!< Adds any ctrl to a treeitem, transferring ownership

    //static TreeControl *GetTree( TreeItemId id ) { return 0; } //\todo If required, to retrieve TI's parent Tree

    TreeItemId GetParentItem( TreeItemId id ) const;
    TreeItemId GetFirstChildItem( TreeItemId id ) const;
    TreeItemId GetNextSiblingItem( TreeItemId id ) const;

private:
    TreeItemControl *m_pRoot;
    int m_MarginH;
    int m_MarginV;
};

}} //namespace sfr::gui

#endif //SFR_GUI_TREE_CONTROL_H
