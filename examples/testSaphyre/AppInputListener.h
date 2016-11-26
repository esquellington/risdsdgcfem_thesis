#ifndef TEST_SAPHYRE_APP_INPUT_LISTENER_H
#define TEST_SAPHYRE_APP_INPUT_LISTENER_H

#include "AppTask.h"

#include <Safra/Safra.h>
#include <Safra/gui/CameraControllerML.h>
#include <GL/freeglut.h> //TEMP: Required for glutFullScreenToggle, NOT in standard glut.h

/*! Centralized input listener for app main window */
class AppInputListener: public sfr::gui::IKeyboardListener, public sfr::gui::IMouseListener
{
public:
    enum EMode {
        eMode_Camera         = 0,
        eMode_CreateParticle = 1,
        eMode_Pressure       = 2,
        eMode_RayCast        = 3,
        eMode_Drag           = 4,
        eMode_Nail           = 5,
        eMode_MoveOrRotate   = 6,
        eMode_Contact        = 7,
        eMode_FixBoundaryPlane = 8, //TEMP for DCLFEM 3D
#ifdef __ENABLE_WIND
        eMode_Wind           = 9,
#endif
        eNumModes,
        eMode_Default        = eMode_Camera
    };

    enum EAxisFlags {
        eAxis_None = 0,
        eAxis_X    = 1,
        eAxis_Y    = 2,
        eAxis_All  = 0xFFFFFFFF
    };

public:
    AppInputListener( AppTask &app_task, sfr::gfx::ICamera *p_camera, sfr::gui::ECameraControllerMode cm )
    : m_rAppTask( app_task )
    , m_CC( p_camera, cm )
    , m_Mode( eMode_Default )
    , m_AxisFlags( eAxis_All )
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
            else if( 'q' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Raycast", 2 );
                m_Mode = eMode_RayCast;
                return true;
            }
            else if( 'p' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Pressure", 2 );
                m_Mode = eMode_Pressure;
                return true;
            }
            else if( 'd' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Drag", 2 );
                m_Mode = eMode_Drag;
                return true;
            }
            else if( 'n' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Nail", 2 );
                m_Mode = eMode_Nail;
                return true;
            }
            else if( 'm' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Move | Rotate", 2 );
                m_Mode = eMode_MoveOrRotate;
                return true;
            }
            else if( 'c' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Contact", 2 );
                m_Mode = eMode_Contact;
                return true;
            }
            else if( sfr::gui::eKey_Tab == key )
            {
                m_Mode = eMode_Camera;
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Camera", 2 );
                return true;
            }
            else if( sfr::gui::eKey_Space == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Particle", 2 );
                m_Mode = eMode_CreateParticle;
                return true;
            }
            else if( 'x' == key )
            {
                if( m_AxisFlags.Test( eAxis_X ) ) m_rAppTask.GetAppView()->GetDesktop()->Message( "Axis: X OFF", 2 );
                else m_rAppTask.GetAppView()->GetDesktop()->Message( "Axis: X ON", 2 );
                m_AxisFlags.Toggle( eAxis_X );
                return true;
            }
            else if( 'y' == key )
            {
                if( m_AxisFlags.Test( eAxis_Y ) ) m_rAppTask.GetAppView()->GetDesktop()->Message( "Axis: Y OFF", 2 );
                else m_rAppTask.GetAppView()->GetDesktop()->Message( "Axis: Y ON", 2 );
                m_AxisFlags.Toggle( eAxis_Y );
                return true;
            }
            else if( '2' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Dimension: 2, Mode = Camera", 2 );
                m_Mode = eMode_Camera;
                m_rAppTask.SetDimension( 2 );
                m_rAppTask.GetAppView()->GetCamera()->InitOrtho( sfr::Vec3(-1.0,-3.0,10.0), sfr::Vec3(-1,-3,0), sfr::Vec3(0,1,0), //sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                                                 10.0, 1.0,
                                                                 -100.0, 500.0,
                                                                 600, 600 );
            }
            else if( '3' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Dimension: 3, Mode = Camera", 2 );
                m_Mode = eMode_Camera;
                m_rAppTask.SetDimension( 3 );
                m_rAppTask.GetAppView()->GetCamera()->InitPerspective( sfr::Vec3(0,1.5,3), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0), //sfr::Vec3(5,7,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                                                       45.0, 1.0,
                                                                       0.1, 1000.0,
                                                                       600, 600 );
            }
            else if( sfr::gui::eKey_F1 == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: FirstPerson", 2 );
                m_Mode = eMode_Camera;
                m_CC.SetMode( sfr::gui::eCameraController_FirstPerson );
                return true;
            }
            else if( sfr::gui::eKey_F2 == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: Planar", 2 );
                m_Mode = eMode_Camera;
                m_CC.SetMode( sfr::gui::eCameraController_Planar );
                return true;
            }
            else if( sfr::gui::eKey_F3 == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: Spherical", 2 );
                m_Mode = eMode_Camera;
                m_CC.SetMode( sfr::gui::eCameraController_Spherical );
                return true;
            }
            else if( sfr::gui::eKey_F12 == key )
            {
                if( m_CC.GetSpeedScale() == 1.0f )
                {
                    m_CC.SetSpeedScale(1.0f/4);
                    m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: 1/4x Speed", 2 );
                }
                else if( m_CC.GetSpeedScale() == 1.0f/4 )
                {
                    m_CC.SetSpeedScale(1.0f/16);
                    m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: 1/16x Speed", 2 );
                }
                else
                {
                    m_CC.SetSpeedScale(1.0f);
                    m_rAppTask.GetAppView()->GetDesktop()->Message( "Camera: 1x Speed", 2 );
                }
                return true;
            }
            else if( 'b' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Fix Boundary Plane", 2 );
                m_Mode = eMode_FixBoundaryPlane;
                return true;
            }
#ifdef __ENABLE_WIND
            else if( 'w' == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->Message( "Mode: Wind", 2 );
                m_Mode = eMode_Wind;
                return true;
            }
#endif
            else if( sfr::gui::eKey_F4 == key )
            {
                m_rAppTask.GetAppView()->GetDesktop()->SetEnabled( !m_rAppTask.GetAppView()->GetDesktop()->IsEnabled() );
                return true;
            }
            else if( sfr::gui::eKey_F5 == key )
            {
                glutFullScreenToggle();
                return true;
            }
            else if( sfr::gui::eKey_F6 == key )
            {
                m_rAppTask.GetAppView()->SetFlags( m_rAppTask.GetAppView()->GetFlags().Toggle(sfr::IView::eFlags_EnableFPS) );
                return true;
            }
            else if( sfr::gui::eKey_F7 == key )
            {
                char command_str[1024];
                snprintf( command_str, 1024, "gnome-screenshot -w -B -f screenshot-%2.2f.png", m_rAppTask.GetScene().GetTime() );
                system( command_str );
                return true;
            }
            return false;
        }
    //@}

    //\todo Apply BeginXXXX/EndXXX/XXX protocol to all input commands, calling AppTask for actual execution!!

    //! \name Mouse Event Processing
    //@{
    bool OnMouseButtonPressed( sfr::gui::EMouseButton mb, int x, int y )
        {
            if( eMode_Camera == m_Mode )
                return m_CC.OnMouseButtonPressed(mb,x,y);
            else if( sfr::gui::eMB_Left == mb )
            {
                if( 2 == m_rAppTask.GetDimension() )
                {
                    Vec2f p2d( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
                    if( !m_AxisFlags.Test(eAxis_X) ) p2d.x() = 0;
                    if( !m_AxisFlags.Test(eAxis_Y) ) p2d.y() = 0;

                    if( eMode_CreateParticle == m_Mode ) m_rAppTask.CreateParticle2D( p2d );
                    else if( eMode_RayCast == m_Mode )
                    {
                        m_rAppTask.m_bIsEnabledRayCast = true;
                        m_rAppTask.m_RayCastPoint = p2d;
                    }
                    else if( eMode_Pressure == m_Mode )
                    {
                        m_rAppTask.m_bIsEnabledPressure = true;
                        m_rAppTask.m_PressurePoint = p2d;
                    }
                    else if( eMode_Drag == m_Mode ) m_rAppTask.BeginDrag2D( p2d, false );
                    else if( eMode_Nail == m_Mode ) m_rAppTask.BeginDrag2D( p2d, true );
                    else if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.BeginMove2D( p2d );
                    else if( eMode_Contact == m_Mode ) m_rAppTask.SelectContactEntities2D( p2d );
#ifdef __ENABLE_WIND
                    else if( eMode_Wind == m_Mode ) m_rAppTask.BeginWind( p2d );
#endif
                    return true;
                }
                else // 3 == m_rAppTask.GetDimension()
                {
                    Vec3f p0 = m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0) );
                    Vec3f p1 = m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0.1) );
                    Vec3f dir = mal::Normalized(p1-p0);
                    if( eMode_Drag == m_Mode ) m_rAppTask.BeginDrag3D( p0, dir, false );
                    else if( eMode_Nail == m_Mode ) m_rAppTask.BeginDrag3D( p0, dir, true );
                    else if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.BeginMove3D( p0, dir );
                    else if( eMode_FixBoundaryPlane == m_Mode ) m_rAppTask.FixBoundaryPlane3D( p0, dir );
                    else if( eMode_Contact == m_Mode ) m_rAppTask.SelectContactEntities3D( p0, dir );
                    return true;
                }
            }
            else if( sfr::gui::eMB_Right == mb )
            {
                if( 2 == m_rAppTask.GetDimension() )
                {
                    Vec2f p2d( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
                    if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.BeginRotate2D( p2d );
                }
                else // 3 == m_rAppTask.GetDimension()
                {
                    Vec3f p0 = m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0) );
                    Vec3f p1 = m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0.1) );
                    Vec3f dir = mal::Normalized(p1-p0);
                    if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.BeginRotate3D( p0, dir );
                }
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
                if( 2 == m_rAppTask.GetDimension() )
                {
                    if( eMode_RayCast == m_Mode ) m_rAppTask.m_bIsEnabledRayCast = false;
                    else if( eMode_Pressure == m_Mode ) m_rAppTask.m_bIsEnabledPressure = false;
                    else if( eMode_Drag == m_Mode ) m_rAppTask.EndDrag2D( false );
                    else if( eMode_Nail == m_Mode ) m_rAppTask.EndDrag2D( true );
                    else if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.EndMove2D( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
#ifdef __ENABLE_WIND
                    else if( eMode_Wind == m_Mode ) m_rAppTask.EndWind( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
#endif
                    return true;
                }
                else // 3 == m_rAppTask.GetDimension()
                {
                    if( eMode_Drag == m_Mode ) m_rAppTask.EndDrag3D( false );
                    else if( eMode_Nail == m_Mode ) m_rAppTask.EndDrag3D( true );
                    else if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.EndMove3D( m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0) ) );
                    return true;
                }
            }
            else if( sfr::gui::eMB_Right == mb )
            {
                if( 2 == m_rAppTask.GetDimension() )
                {
                    if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.EndRotate2D( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
                }
                else // 3 == m_rAppTask.GetDimension()
                {
                    if( eMode_MoveOrRotate == m_Mode ) m_rAppTask.EndRotate3D( m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0) ) );
                }
                return true;
            }
            else
                return false;
        }
    bool OnMouseMotion( int x, int y )
        {
            if( eMode_Camera == m_Mode )
                return m_CC.OnMouseMotion(x,y);
            else if( 2 == m_rAppTask.GetDimension() )
            {
                Vec2f p2d( m_CC.GetCamera()->UnProject2( sfr::Vec3(x,y,0) ) );
                if( !m_AxisFlags.Test(eAxis_X) ) p2d.x() = 0;
                if( !m_AxisFlags.Test(eAxis_Y) ) p2d.y() = 0;
                if( eMode_Pressure == m_Mode )
                {
                    m_rAppTask.m_PressurePoint = p2d;
                    return true;
                }
                else if( eMode_RayCast == m_Mode )
                {
                    m_rAppTask.m_RayCastPoint = p2d;
                    return true;
                }
                else if( eMode_Drag == m_Mode )
                {
                    m_rAppTask.Drag2D( p2d, false );
                    return true;
                }
                else if( eMode_Nail == m_Mode )
                {
                    m_rAppTask.Drag2D( p2d, true );
                    return true;
                }
                else if( eMode_MoveOrRotate == m_Mode )
                {
                    m_rAppTask.MoveOrRotate2D( p2d );
                    return true;
                }
#ifdef __ENABLE_WIND
                else if( eMode_Wind == m_Mode )
                {
                    m_rAppTask.Wind( p2d );
                    return true;
                }
#endif
                else
                    return false;
            }
            else // 3 == m_rAppTask.GetDimension()
            {
                Vec3f p0 = m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0) );
                Vec3f p1 = m_CC.GetCamera()->UnProject( sfr::Vec3(x,y,0.1) );
                Vec3f dir = mal::Normalized(p1-p0);
                if( eMode_Drag == m_Mode )
                {
                    m_rAppTask.Drag3D( p0, dir, false );
                    return true;
                }
                else if( eMode_Nail == m_Mode )
                {
                    m_rAppTask.Drag3D( p0, dir, true );
                    return true;
                }
                else if( eMode_MoveOrRotate == m_Mode )
                {
                    m_rAppTask.MoveOrRotate3D( p0, dir );
                    return true;
                }
                return false;
            }
        }
    //@}

private:
    AppTask &m_rAppTask;
    sfr::gui::CameraControllerML m_CC;
    EMode m_Mode;
    Flags32 m_AxisFlags;
};

#endif //TEST_SAPHYRE_APP_INPUT_LISTENER_H
