#ifndef SFR_CORE_GFX_ICAMERA_H
#define SFR_CORE_GFX_ICAMERA_H

#include <Safra/Config.h>
#include <Safra/util/Observable.h>

namespace sfr { namespace gfx
{

enum ECameraProjectionType {
    eCameraProjection_Perspective,
    eCameraProjection_Ortho
    //Could have other camera-projection types here...
};

//! 3D Camera Interface
/*! Observable 3D camera interface. */
class ICamera: public Util::Observable
{
public:
    ICamera() {}
    virtual ~ICamera() {}

    //!\name Camera definition and modification
    //@{
    virtual void InitPerspective ( const Vec3 &from, const Vec3 &at, const Vec3 &up_vector,
                                   Real angle, Real aspect_ratio, //ratio = w/h
                                   Real near_plane, Real far_plane,
                                   int viewport_width, int viewport_height ) {}
    virtual void InitOrtho ( const Vec3 &from, const Vec3 &at, const Vec3 &up_vector,
                             Real view_width, Real aspect_ratio, //ratio = w/h
                             Real near_plane, Real far_plane,
                             int viewport_width, int viewport_height ) {}

    virtual void SetZoom( Real zoom_factor ) = 0;    

    // Global X,Y,Z camera manipulation
    virtual void LookAt( const Vec3 &at ) = 0; //!< Change Vrp
    virtual void MoveObs( const Vec3 &displacement_global ) = 0; //!< Change Obs
    virtual void StrafeWorld( const Vec3 &displacement_world ) = 0; //!< Translate Obs & Vrp locally to camera Right,Up,Front
    /* Not needed by now...
    virtual void RotateWorld( const Vec3 &angles_world ) = 0;   //!< Rotate around Obs, rotation in world axis
    virtual void OrbitWorld( const Vec3 &angles_world ) = 0;    //!< Orbit around Vrp, rotation in world axis
    */           

    // Local Right,Up,Front camera manipulation
    virtual void StrafeRUF( const Vec3 &displacement_ruf ) = 0; //!< Translate Obs & Vrp locally to camera Right,Up,Front
    virtual void RotateRUF( const Vec3 &angles_ruf ) = 0;   //!< Rotate around Obs, rotation in camera axis
    virtual void OrbitRUF( const Vec3 &angles_ruf ) = 0;    //!< Orbit around Vrp, rotation in camera axis
    //@}
    
    //! \name Viewport related methods
    //@{
    virtual void SetAspect( Real aspect_ratio ) = 0;
    virtual Real GetAspect() const = 0;
    virtual Real GetOrtho_ViewWidth() const = 0;
    virtual Real GetPerspective_Angle() const = 0;
    virtual void SetViewportSizes( int w, int h ) = 0;
    virtual void GetViewportSizes( int &w, int &h ) = 0;
    //@}

    //! \name Consultors
    //@{
    virtual Vec3 GetObs() const = 0;
    virtual Vec3 GetVrp() const = 0;
    virtual Vec3 GetViewDir() const = 0;
    virtual Vec3 GetUpVect() const = 0;
    virtual Vec3 GetRightVect() const = 0;
    virtual ECameraProjectionType GetProjectionType() const = 0;
    virtual Real GetZoom() const = 0;
    //@}

    //! \name Camera preparation before rendering methods
    //@{
    virtual void PrepareToShoot() const = 0;
    virtual void PrepareToSelect( int px, int py ) const = 0;
    
    virtual void PrepareToShootViewport() const = 0;    
    //@}

    //! \name Utility Methods
    //@{
    virtual Vec3 UnProject( const Vec3 &p3d_viewport ) const = 0;
    virtual Vec3 Project( const Vec3 &p3d_world ) const = 0;
    finline Vec2 UnProject2( const Vec3 &p3d_viewport ) const { Vec3 p3d( UnProject(p3d_viewport) ); return Vec2(p3d.x(),p3d.y()); }
    //@}
};

} } // namespace sfr::gfx
  
#endif //SFR_CORE_GFX_ICAMERA_H
