#ifndef SFR_CORE_GUI_IWIDGET_H
#define SFR_CORE_GUI_IWIDGET_H

#include <Safra/Config.h>
#include <boost/function.hpp>
#include <util/ItemStream.h>

namespace sfr { namespace gui
{

enum EDesktopWidgetFlags {
    eDWF_None         = 0,
    eDWF_CanMove      = (1<<0),
    eDWF_CanRoll      = (1<<1),
    eDWF_CanMinimize  = (1<<2),
    eDWF_CanClose     = (1<<3),
    eDWF_CanAddButton = (1<<4),
    eDWF_Default      = ( eDWF_CanMove
                          | eDWF_CanRoll
                          | eDWF_CanMinimize
                          | eDWF_CanClose
                          | eDWF_CanAddButton )
};

typedef boost::function<void (void *)> D_NotifyCall;


// \todo Define public interfaces and implement them, ensure they are PURE VIRTUAL!
class IWidget
{
public:
    //\name Specific widget APIs
    //@{
    class IBasicAPI;
    class ITreeAPI;
    class IPropertyTreeAPI;
    class IPlotAPI;
    //@}
    
public:
    IWidget() {}
    virtual ~IWidget() {}

    virtual IBasicAPI *GetBasicAPI() { return 0; }
    virtual ITreeAPI *GetTreeAPI() { return 0; }
    virtual IPropertyTreeAPI *GetPropertyTreeAPI() { return 0; }
    virtual IPlotAPI *GetPlotAPI() { return 0; }
};

/*\todo Widget types CANNOT be subclasses of IWidget because they are
  all implemented by WidgetControl.h with different Body
  controls... in krtDD this worked because krtDW was NOT a parent
  class of CWidgetCtrl, but a PIMPL... maybe it'd be best to use the
  same approach here to... or retrieve specific widget APIs from
  IWidget, implemented by the Body, not the Widget, which could then
  derive from a proper ITreeAPI, IPropertyTreeAPI interface...  */
class IWidget::IBasicAPI
{
public:
    virtual void Close() = 0;    
    virtual void Minimize() = 0;
    virtual void Restore() = 0;
    virtual void Roll() = 0;
    virtual void Unroll() = 0;
    virtual bool AddFootButton( const char *name, D_NotifyCall notify_call, void *p_user_data = 0 ) = 0;

    virtual void SetColor( float r, float g, float b, float a ) = 0;
};

typedef util::ItemStream::ItemItRW PropertyIt;

class IWidget::IPropertyTreeAPI
{
public:
    virtual void Clear() = 0;
    virtual void Rebuild( PropertyIt property_tree_it ) = 0;
};

} } // namespace sfr::gui

#endif // SFR_CORE_GUI_IWIDGET_H
