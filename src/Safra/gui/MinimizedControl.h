#ifndef SFR_GUI_MINIMIZED_CONTROL_H
#define SFR_GUI_MINIMIZED_CONTROL_H

#include "LayoutControl.h"

namespace sfr { namespace gui
{

class MinimizedControl: public VerticalLayout
{
public:
    MinimizedControl()
    {
        SetMinSeparation( 0 );
        SetColor( gfx::Color(0,0,0,0) );
    }
    virtual ~MinimizedControl() {}

    /* TEMP: Only for debug
    bool AddChild( BaseControl *p_ctrl )
        {
            SFR_LOG("MinimizedControl::AddChild( %llx )", (machine_uint_type)p_ctrl );
            return VerticalLayout::AddChild(p_ctrl);
        }
    bool RemoveChild( BaseControl *p_ctrl )
        {
            SFR_LOG("MinimizedControl::RemoveChild( %llx )", (machine_uint_type)p_ctrl );
            return VerticalLayout::RemoveChild(p_ctrl);
        }
    */
};

}} //namespace sfr::gui

#endif //SFR_GUI_MINIMIZED_CONTROL_H
