#ifndef SFR_GUI_CAMERA_CONTROLLER_ML_H
#define SFR_GUI_CAMERA_CONTROLLER_ML_H

#include <Safra/Config.h>
#include <Safra/core/gui/IMouseListener.h>
#include <Safra/core/gfx/ICamera.h>
#include <Safra/util/IObserver.h>

#include <Mal/GConversion.h>

//
//#include "../../../dep/GLUI/glui_v2_2/Arcball.h"

namespace sfr { namespace gui
{

enum ECameraControllerMode {
    eCameraController_Planar,     //Translate on camera X,Y and zoom on cameraZ
    eCameraController_Spherical,
    eCameraController_ArcBall,
    eCameraController_FirstPerson
};

class CameraControllerML: public IMouseListener, public Util::IObserver
{

public:
    inline CameraControllerML( gfx::ICamera *p_camera, ECameraControllerMode cm )
    : m_pCamera( p_camera )
    , m_Mode( cm )
    , m_SpeedScale(1)
    , m_LastX(0), m_LastY(0)
    , m_bIsLeftPressed( false )
    , m_bIsRightPressed( false )
    , m_bIsMiddlePressed( false )
    {
        AddObservable( m_pCamera );
        m_pCamera->GetViewportSizes( m_Width, m_Height );
    }
    ~CameraControllerML() {}

    //!\name Configuration
    //@{
    inline void SetMode( ECameraControllerMode cm ) { m_Mode = cm; }
    inline void SetSpeedScale( Real speed ) { m_SpeedScale = speed; }
    inline ECameraControllerMode GetMode() const { return m_Mode; }
    inline Real GetSpeedScale() const { return m_SpeedScale; }
    //@}

    const gfx::ICamera *GetCamera() const { return m_pCamera; }

    //! \name Mouse Event Processing
    //@{
    inline bool OnMouseButtonPressed( EMouseButton mb, int x, int y )
    {
        m_bIsLeftPressed = (eMB_Left == mb);
        m_bIsRightPressed = (eMB_Right == mb);
        m_bIsMiddlePressed = (eMB_Middle == mb);
        if( m_bIsLeftPressed || m_bIsRightPressed || m_bIsMiddlePressed )
        {
            m_LastX = x;
            m_LastY = y;
        }
        return true;
    }

    inline bool OnMouseButtonReleased( EMouseButton mb, int x, int y )
    {
        if (eMB_Left == mb) m_bIsLeftPressed = false;
        if (eMB_Right == mb) m_bIsRightPressed = false;
        if (eMB_Middle == mb) m_bIsMiddlePressed = false;
        return true;
    }

    inline bool OnMouseMotion( int x, int y )
    {
        switch( m_Mode )
        {
        case eCameraController_Planar:
            if( m_bIsRightPressed )
            {
                // Drag in viewport coords, compute 1:1 movement in world coords
                // \todo Maybe move it to a DragOrtho() camera method...
                Vec3 displ_viewport( -Real(x-m_LastX), Real(y-m_LastY), 0 );
                Vec3 displ_ruf( displ_viewport * m_pCamera->GetOrtho_ViewWidth() / (m_Width * m_pCamera->GetZoom()) );
                m_pCamera->StrafeRUF( displ_ruf );

                /* Simpler alternative
                Vec3 last_point( m_pCamera->UnProject( Vec3(m_LastX,m_LastY,0) ) );
                Vec3 current_point( m_pCamera->UnProject( Vec3(x,y,0) ) );
                m_pCamera->StrafeWorld( last_point - current_point );
                */

                m_LastX = x;
                m_LastY = y;
            }
            if( m_bIsMiddlePressed )
            {
                // Zoom in/out uniformly (must divide increment if scale < 0)
                Real zoom_inc( (Real(y-m_LastY)/m_Height) );
                Real zoom_factor = m_pCamera->GetZoom() + ( (m_pCamera->GetZoom() > 0) ? zoom_inc : 1.0f/zoom_inc );
                m_pCamera->SetZoom( mal::Max( 0.001f, zoom_factor ) );
                m_LastX = x;
                m_LastY = y;
            }
            break;
        case eCameraController_Spherical:
            if( m_bIsLeftPressed )
            {
                Real angle_up( -Real(x-m_LastX)/m_Width*Real(MAL_CONSTANT_PI) );
                Real angle_right( -Real(y-m_LastY)/m_Height*Real(MAL_CONSTANT_PI) );
                m_pCamera->OrbitRUF( Vec3(angle_right,angle_up,0) );
                m_LastX = x;
                m_LastY = y;
            }
            if( m_bIsRightPressed )
            {
                /* The following code does not not work well with
                   perspective-cameras due to distortion... should use
                   a scene raycast to determine proper camera
                   depth-plane on which the world strafe displacement
                   must be computed
                Vec3 last_point( m_pCamera->UnProject( Vec3(m_LastX,m_LastY,0) ) );
                Vec3 current_point( m_pCamera->UnProject( Vec3(x,y,0) ) );

                Vec3 last_point( last_ray_cast_hit_point );
                Vec3 current_point( m_pCamera->UnProject( Vec3(x,y,same_viewport_z_coord_as_lrhp) ) );

                m_pCamera->StrafeWorld( last_point - current_point );
                */
                Real scale(m_SpeedScale*25.0f);
                m_pCamera->StrafeRUF( scale*Vec3( -Real(x-m_LastX)/m_Width, Real(y-m_LastY)/m_Height, 0) );

                m_LastX = x;
                m_LastY = y;
            }
            if( m_bIsMiddlePressed )
            {
                //This didn't work properly... seems to be forgotten old code...
                Real scale(m_SpeedScale*50.0f);
                m_pCamera->MoveObs( scale * Real(y-m_LastY)/m_Height * m_pCamera->GetViewDir() );
                m_LastX = x;
                m_LastY = y;
                //

                /*\todo THIS DOES NOTHING IN PERSPECTIVE!! BUT WORKS IN ORTHO... Zoom in/out uniformly (must divide increment if scale < 0)
                Real zoom_inc( (Real(y-m_LastY)/m_Height) );
                Real zoom_factor = m_pCamera->GetZoom() + ( (m_pCamera->GetZoom() > 0) ? zoom_inc : 1.0f/zoom_inc );
                m_pCamera->SetZoom( mal::Max( 0.001f, zoom_factor ) );
                m_LastX = x;
                m_LastY = y;
                */
            }
            break;
        case eCameraController_ArcBall:
            break;
        case eCameraController_FirstPerson:
            if( m_bIsLeftPressed )
            {
                Real angle_up( -Real(x-m_LastX)/m_Width*Real(MAL_CONSTANT_PI) );
                Real angle_right( -Real(y-m_LastY)/m_Height*Real(MAL_CONSTANT_PI) );
                m_pCamera->RotateRUF( Vec3(angle_right,angle_up,0) );
                m_LastX = x;
                m_LastY = y;
            }
            if( m_bIsRightPressed )
            {
                Real scale(m_SpeedScale*25.0f);
                m_pCamera->StrafeRUF( scale*Vec3( Real(x-m_LastX)/m_Width, -Real(y-m_LastY)/m_Height, 0) );
                m_LastX = x;
                m_LastY = y;
            }
            if( m_bIsMiddlePressed )
            {
                Real scale(m_SpeedScale*50.0f);
                m_pCamera->StrafeRUF( Vec3( 0, 0, scale * Real(y-m_LastY)/m_Height ) );
                m_LastX = x;
                m_LastY = y;
            }
            break;
        default: break;
        };
        return true;
    }
    //@}


    /*! Observable change notification
      - Camera has changed
    */
    inline void NotifyChanged( const Util::Observable *p_observable )
    {
        if( p_observable == (Util::Observable*)m_pCamera )
        {
            m_pCamera->GetViewportSizes( m_Width, m_Height );
            m_bIsLeftPressed = false;
            m_bIsRightPressed = false;
            m_bIsMiddlePressed = false;
        }
    }

private:
    gfx::ICamera *m_pCamera;
    int m_Width;
    int m_Height;
    ECameraControllerMode m_Mode;
    Real m_SpeedScale;
    //Arcball m_Arcball;

    //! \name Mouse-Controller State
    //@{
    int m_LastX, m_LastY;
    bool m_bIsLeftPressed, m_bIsRightPressed, m_bIsMiddlePressed;
    //@}
};

} } // namespace sfr::gui

#endif // SFR_GUI_CAMERA_CONTROLLER_ML_H
