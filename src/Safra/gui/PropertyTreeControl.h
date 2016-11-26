#ifndef SFR_GUI_PROPERTY_TREE_CONTROL_H
#define SFR_GUI_PROPERTY_TREE_CONTROL_H

#include "LayoutControl.h"
#include "TreeControl.h"
#include <util/ItemStream.h>
#include <boost/function.hpp>

namespace sfr { namespace gui {

typedef util::ItemStream::ItemItRW PropertyIt;
typedef boost::function<void (PropertyIt)> D_SyncCall;

class PropertyTreeControl: public FrameLayout
{
public:
    PropertyTreeControl( PropertyIt property_tree_it, D_SyncCall sync_call = D_SyncCall() );
    ~PropertyTreeControl();

    void Clear();
    void Rebuild( PropertyIt property_tree_it );    
    bool Update( float dt );
    
private:
    void AddComplexProperty( PropertyIt pit, TreeItemId tid );
    void AddComplexProperty_NIR( PropertyIt pit, TreeItemId tid );
    void AddComplexProperty_Enum( PropertyIt pit, TreeItemId tid );
    void AddSimpleProperty( PropertyIt pit, TreeItemId tid );
    
private:
    PropertyIt m_PropertyTreeIT;
    TreeControl m_TreeControl;
    D_SyncCall m_SyncCall;

private:
    //!\name Implemetation of IPropertyTreeAPI
    //@{
    class Impl_PropertyTreeAPI: public IWidget::IPropertyTreeAPI
    {
    public:
        Impl_PropertyTreeAPI( PropertyTreeControl *p_ptc ) : m_pPropertyTreeControl(p_ptc) {}
        void Clear() { m_pPropertyTreeControl->Clear(); }
        virtual void Rebuild( PropertyIt property_tree_it ) { m_pPropertyTreeControl->Rebuild(property_tree_it); }
    private:
        PropertyTreeControl *m_pPropertyTreeControl;
    };
    Impl_PropertyTreeAPI m_PropertyTreeAPI;
    virtual IWidget::IPropertyTreeAPI *GetAPI_PropertyTree() { return &m_PropertyTreeAPI; }
    //@}
};

}} //namespace sfr::gui

#endif //SFR_GUI_PROPERTY_TREE_CONTROL_H
