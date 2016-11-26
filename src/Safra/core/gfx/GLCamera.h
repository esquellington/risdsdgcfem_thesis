#ifndef SFR_CORE_GFX_GLCAMERA_H
#define SFR_CORE_GFX_GLCAMERA_H

#include <Safra/core/gfx/ICamera.h>

namespace sfr { namespace gfx
{

class GLCamera: public ICamera
{
public:
    GLCamera();
    ~GLCamera() {}

    //!\name Camera definition and modification
    //@{
    void InitPerspective ( const Vec3 &from, const Vec3 &at, const Vec3 &up_vector,
                           Real angle, Real aspect_ratio, //ratio = w/h
                           Real near_plane, Real far_plane,
                           int viewport_width, int viewport_height );
    void InitOrtho ( const Vec3 &from, const Vec3 &at, const Vec3 &up_vector,
                     Real view_width, Real aspect_ratio, //ratio = w/h
                     Real near_plane, Real far_plane,
                     int viewport_width, int viewport_height );
    
    void SetZoom( Real zoom_factor ) { m_ZoomFactor = zoom_factor; }

    // Global X,Y,Z camera manipulation
    void LookAt( const Vec3 &at );
    void MoveObs( const Vec3 &displacement_world );
    void StrafeWorld( const Vec3 &displacement_world );
    
    // Local Right,Up,Front camera manipulation    
    void StrafeRUF( const Vec3 &displacement_ruf );
    void RotateRUF( const Vec3 &angles_ruf );
    void OrbitRUF( const Vec3 &angles_ruf );
    //@}

    //! \name Viewport related methods
    //@{
    void SetAspect( Real aspect_ratio );
    Real GetAspect() const { return m_Aspect; }
    Real GetOrtho_ViewWidth() const { return m_Ortho_ViewWidth; }
    Real GetPerspective_Angle() const { return m_Perspective_Angle; }
    void SetViewportSizes( int w, int h );
    void GetViewportSizes( int &w, int &h ) { w = m_ViewportWidth; h = m_ViewportHeight; }
    //@}
    
    //! \name Consultors
    //@{
    Vec3 GetObs() const { return m_Obs; }
    Vec3 GetVrp() const { return m_Vrp; }
    Vec3 GetViewDir() const { return (m_Vrp-m_Obs).Normalized(); }
    Vec3 GetUpVect() const { return m_UpVec; }
    Vec3 GetRightVect() const { return -m_UpVec % GetViewDir(); }
    ECameraProjectionType GetProjectionType() const { return m_ProjectionType; }
    Real GetZoom() const { return m_ZoomFactor; }
    //@}
    
    //! \name Camera preparation before rendering methods (TEMPORAL!!)
    //@{
    void PrepareToShoot() const;
    void PrepareToSelect( int px, int py ) const;

    void PrepareToShootPerspective() const;
    void PrepareToSelectPerspective( int px, int py ) const;
    
    void PrepareToShootOrtho() const;
    void PrepareToSelectOrtho( int px, int py ) const;

    void PrepareToShootViewport() const;    
    //@}

    //! \name Utility Methods
    //@{
    Vec3 UnProject( const Vec3 &p3d_window ) const;
    Vec3 Project( const Vec3 &p3d_world ) const;
    //@}

private:
    void SyncUpToViewDir(); //!< Called whenever Obs or Vrp change
        
protected:
    ECameraProjectionType m_ProjectionType;
    
    Vec3 m_Obs;
    Vec3 m_Vrp;
    Vec3 m_UpVec;

    Real m_ZoomFactor;
    Real m_Perspective_Angle;
    Real m_Ortho_ViewWidth;
    Real m_Aspect;
    Real m_NearClip, m_FarClip;

    //! \name View related attribs
    //@{
    int m_ViewportWidth;
    int m_ViewportHeight;
    //@}
};

} } // namespace sfr::gfx
  
#endif //SFR_CORE_GFX_GLCAMERA_H
