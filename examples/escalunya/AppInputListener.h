#ifndef ESCALUNYA_APP_INPUT_LISTENER_H
#define ESCALUNYA_APP_INPUT_LISTENER_H

#include "AppTask.h"

/*! Centralized input listener for app main window
  - If Ctrl is pressed:
    - Dispatch to camera controller
  - Otherwise:
    ???
*/
class AppInputListener: public sfr::gui::IKeyboardListener, public sfr::gui::IMouseListener
{
public:
    enum EMode { eMode_Camera
#ifdef __ENABLE_BVH_TEST
                 , eMode_TestBVH
#endif
    };
public:
    AppInputListener( AppTask &app_task, sfr::gfx::ICamera *p_camera, sfr::gui::ECameraControllerMode cm )
    : m_rAppTask( app_task )
    , m_CC( p_camera, cm )
    , m_Mode( eMode_Camera )
        {}

    //! \name Keyboard Event Processing
    //@{
    bool OnKeyPressed( sfr::gui::EKey key, int x, int y )
        {
            return false;
        }
    bool OnKeyReleased( sfr::gui::EKey key, int x, int y )
        {
            if( sfr::gui::eKey_Tab == key )
            {
                m_Mode = eMode_Camera;
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Camera", 2 );
                return true;
            }
#ifdef __ENABLE_BVH_TEST
            else if( 'b' == key )
            {
                m_Mode = eMode_TestBVH;
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Test BVH", 2 );
                return true;
            }
#endif
            return false;
        }
    //@}

    //! \name Mouse Event Processing
    //@{
    bool OnMouseButtonPressed( sfr::gui::EMouseButton mb, int x, int y )
        {
            if( eMode_Camera == m_Mode ) return m_CC.OnMouseButtonPressed(mb,x,y);
#ifdef __ENABLE_BVH_TEST
            else if( eMode_TestBVH == m_Mode ) { return m_rAppTask.TestBVH( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) ); } //\todo 3D RAYCAST
#endif
            else return false;
        }
    bool OnMouseButtonReleased( sfr::gui::EMouseButton mb, int x, int y )
        {
            if( eMode_Camera == m_Mode ) return m_CC.OnMouseButtonReleased(mb,x,y);
            else return false;
        }
    bool OnMouseMotion( int x, int y )
        {
            if( eMode_Camera == m_Mode ) return m_CC.OnMouseMotion(x,y);
            else return false;
        }
    //@}

private:
    AppTask &m_rAppTask;
    sfr::gui::CameraControllerML m_CC;
    EMode m_Mode;
};

#endif //ESCALUNYA_APP_INPUT_LISTENER_H
