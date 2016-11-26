#ifndef TEST_FE_APPINPUTLISTENER_H
#define TEST_FE_APPINPUTLISTENER_H

#include "Config.h"
#include "AppTask.h"

#include <Safra/Safra.h>
#include <Safra/gui/CameraControllerML.h>

/*! Centralized input listener for app main window
  - If Ctrl is pressed:
    - Dispatch to camera controller
  - Otherwise:
    \todo
*/
class AppInputListener: public sfr::gui::IKeyboardListener, public sfr::gui::IMouseListener
{
public:
    enum EMode {
        eMode_Camera         = 0,
        eMode_Drag           = 1,
        eNumModes,
        eMode_Default        = eMode_Camera
    };
public:
    AppInputListener( AppTask &app_task, sfr::gfx::ICamera *p_camera, sfr::gui::ECameraControllerMode cm )
    : m_rAppTask( app_task )
    , m_CC( p_camera, cm )
    , m_Mode( eMode_Default )
        {}

    //! \name Keyboard Event Processing
    //@{
    bool OnKeyPressed( sfr::gui::EKey key, int x, int y )
        {
            return false;
        }
    bool OnKeyReleased( sfr::gui::EKey key, int x, int y )
        {
            if( 'r' == key )
            {
                m_rAppTask.SetActive( !m_rAppTask.IsActive() );
                if( m_rAppTask.IsActive() ) m_rAppTask.GetAppView()->GetDesktop()->Message( "Run", 2 );
                else m_rAppTask.GetAppView()->GetDesktop()->Message( "Paused", 2 );
                return true;
            }
            else if( 's' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Step", 2 );
                m_rAppTask.SetActive(true);
                m_rAppTask.Update(cDefaultTimeStep);
                m_rAppTask.SetActive(false);
                return true;
            }
            else if( 'd' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Drag", 2 );
                m_Mode = eMode_Drag;
                return true;
            }
            else if( sfr::gui::eKey_Tab == key )
            {
                m_Mode = eMode_Camera;
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Camera", 2 );
                return true;
            }
            else if( sfr::gui::eKey_F2 == key )
            {
                m_rAppTask.SetElementDimension(2);
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Element: 2D", 2 );
                return true;
            }
            else if( sfr::gui::eKey_F3 == key )
            {
                m_rAppTask.SetElementDimension(3);
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Element: 3D", 2 );
                return true;
            }
            else if( '2' == key )
            {
                m_CC.SetMode( sfr::gui::eCameraController_Planar );
                int w,h;
                m_rAppTask.GetAppView()->GetCamera()->GetViewportSizes(w,h);
                m_rAppTask.GetAppView()->GetCamera()->InitOrtho( sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                                                 10.0, float(w)/h,
                                                                 1.0, 500.0,
                                                                 w, h ); //\note Must adapt to current res...
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: Planar", 2 );
                return true;
            }
            else if( '3' == key )
            {
                m_CC.SetMode( sfr::gui::eCameraController_Spherical );
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: Spherical", 2 );
                return true;
            }
            return false;
        }
    //@}

    //! \name Mouse Event Processing
    //@{
    bool OnMouseButtonPressed( sfr::gui::EMouseButton mb, int x, int y )
        {
            Vec2f p2d( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
            if( eMode_Camera == m_Mode )
                return m_CC.OnMouseButtonPressed(mb,x,y);
            else if( sfr::gui::eMB_Left == mb )
            {
                if( eMode_Drag == m_Mode ) m_rAppTask.BeginDrag( p2d );
                return true;
            }
            else
                return false;
        }
    bool OnMouseButtonReleased( sfr::gui::EMouseButton mb, int x, int y )
        {
            if( eMode_Camera == m_Mode )
                return m_CC.OnMouseButtonReleased(mb,x,y);
            else if( sfr::gui::eMB_Left == mb )
            {
                if( eMode_Drag == m_Mode ) m_rAppTask.EndDrag();
                return true;
            }
            else
                return false;
        }
    bool OnMouseMotion( int x, int y )
        {
            Vec2f p2d( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
            if( eMode_Camera == m_Mode )
                return m_CC.OnMouseMotion(x,y);
            else if( eMode_Drag == m_Mode )
            {
                m_rAppTask.Drag( p2d );
                return true;
            }
            else
                return false;
        }
    //@}

private:
    AppTask &m_rAppTask;
    sfr::gui::CameraControllerML m_CC;
    EMode m_Mode;
};

#endif //TEST_FE_APPINPUTLISTENER_H
