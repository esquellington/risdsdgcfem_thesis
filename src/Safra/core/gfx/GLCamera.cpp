#include <Safra/core/gfx/GLCamera.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <Mal/GConversion.h>

namespace sfr { namespace gfx
{

GLCamera::GLCamera()
{
    InitPerspective( Vec3(0,0,1), Vec3(0,0,0), Vec3(0,1,0), 45, 1.0, 1.0, 1000.0, 200, 200 );
}

void GLCamera::InitPerspective( const Vec3 &from, const Vec3 &at, const Vec3 &up_vector,
                                Real angle, Real aspect_ratio,
                                Real near_plane, Real far_plane,
                                int viewport_width, int viewport_height )
{
    m_ProjectionType = eCameraProjection_Perspective;
    m_Obs = from;
    m_Vrp = at;
    m_UpVec = up_vector;
    SyncUpToViewDir();

    m_ZoomFactor = 1.0f;
    m_Perspective_Angle = angle;
    m_Aspect = aspect_ratio;
    m_NearClip = near_plane;
    m_FarClip = far_plane;
    m_ViewportWidth = viewport_width;
    m_ViewportHeight = viewport_height;
}

void GLCamera::InitOrtho( const Vec3 &from, const Vec3 &at, const Vec3 &up_vector,
                          Real view_width, Real aspect_ratio,
                          Real near_plane, Real far_plane,
                          int viewport_width, int viewport_height )
{
    m_ProjectionType = eCameraProjection_Ortho;
    m_Obs = from;
    m_Vrp = at;
    m_UpVec = up_vector;
    SyncUpToViewDir();

    m_ZoomFactor = 1.0f;
    m_Ortho_ViewWidth = view_width;
    m_Aspect = aspect_ratio;
    m_NearClip = near_plane;
    m_FarClip = far_plane;
    m_ViewportWidth = viewport_width;
    m_ViewportHeight = viewport_height;
}

void GLCamera::LookAt( const Vec3 &at )
{
    m_Vrp = at;
    SyncUpToViewDir();
}

void GLCamera::MoveObs( const Vec3 &displacement_world )
{
    m_Obs += displacement_world;
    SyncUpToViewDir();
}

void GLCamera::StrafeWorld( const Vec3 &displacement_world )
{
    m_Obs += displacement_world;
    m_Vrp += displacement_world;
}

void GLCamera::StrafeRUF( const Vec3 &displacement_ruf )
{
    // Compute global vector
    Vec3 displacement_world( displacement_ruf.x() * GetRightVect()
                             + displacement_ruf.y() * GetUpVect()
                             + displacement_ruf.z() * GetViewDir() );
    StrafeWorld( displacement_world );
}

void GLCamera::RotateRUF( const Vec3 &angles_ruf )
{
    // Compute combined rotation in local camera axis
    Mat3x3 rot_up,rot_right,rot_front;
    rot_up = mal::GRotation3x3_From( GetUpVect(), angles_ruf.y() );
    rot_right = mal::GRotation3x3_From( GetRightVect(), angles_ruf.x() );
    Mat3x3 rot = rot_right * rot_up;
    // Rotate Vrp around Obs
    m_Vrp = m_Obs + rot*(m_Vrp-m_Obs); //\todo MAY accumulate numeric drift and change actual distance obs->vrp
    // Rotate camera vectors
    m_UpVec = rot*m_UpVec;
    m_UpVec.Normalize();
}

void GLCamera::OrbitRUF( const Vec3 &angles_ruf )
{
    // Compute combined rotation in local camera axis
    Mat3x3 rot_up,rot_right,rot_front;
    rot_up = mal::GRotation3x3_From( GetUpVect(), angles_ruf.y() );
    rot_right = mal::GRotation3x3_From( GetRightVect(), angles_ruf.x() );
    Mat3x3 rot = rot_right * rot_up;
    // Rotate Obs around Vrp
    m_Obs = m_Vrp + rot*(m_Obs-m_Vrp); //\todo MAY accumulate numeric drift and change actual distance obs->vrp
    // Rotate camera vectors
    m_UpVec = rot*m_UpVec;
    m_UpVec.Normalize();
}

void GLCamera::SetAspect( Real aspect_ratio )
{
    m_Aspect = aspect_ratio;
    Util::Observable::NotifyChanged();
}

void GLCamera::SetViewportSizes( int w, int h )
{
    m_ViewportWidth = w;
    m_ViewportHeight = h;
    Util::Observable::NotifyChanged();
}

void GLCamera::PrepareToShoot() const
{
    switch( m_ProjectionType )
    {
    case eCameraProjection_Perspective: PrepareToShootPerspective(); break;
    case eCameraProjection_Ortho: PrepareToShootOrtho(); break;
    default: break;
    }
}

void GLCamera::PrepareToSelect( int px, int py ) const
{
    switch( m_ProjectionType )
    {
    case eCameraProjection_Perspective: PrepareToSelectPerspective(px,py); break;
    case eCameraProjection_Ortho: PrepareToSelectOrtho(px,py); break;
    default: break;
    }
}

void GLCamera::PrepareToShootPerspective() const
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective(m_Perspective_Angle,m_Aspect,m_NearClip,m_FarClip);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt( m_Obs.x(),m_Obs.y(),m_Obs.z(),
               m_Vrp.x(),m_Vrp.y(),m_Vrp.z(),
               m_UpVec.x(),m_UpVec.y(),m_UpVec.z() );
}

void GLCamera::PrepareToSelectPerspective( int px, int py ) const
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPickMatrix(px,viewport[3]-py,5,5,viewport);
    gluPerspective(m_Perspective_Angle,m_Aspect,m_NearClip,m_FarClip);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt( m_Obs.x(),m_Obs.y(),m_Obs.z(),
               m_Vrp.x(),m_Vrp.y(),m_Vrp.z(),
               m_UpVec.x(),m_UpVec.y(),m_UpVec.z() );
}

void GLCamera::PrepareToShootOrtho() const
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

    Real half_ortho_view_width = (m_Ortho_ViewWidth / 2.0f) / m_ZoomFactor;
    Real half_ortho_view_height = half_ortho_view_width / m_Aspect;
    glOrtho( -half_ortho_view_width, half_ortho_view_width,
             -half_ortho_view_height, half_ortho_view_height,
             m_NearClip, m_FarClip );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt( m_Obs.x(),m_Obs.y(),m_Obs.z(),
               m_Vrp.x(),m_Vrp.y(),m_Vrp.z(),
               m_UpVec.x(),m_UpVec.y(),m_UpVec.z() );
}

void GLCamera::PrepareToSelectOrtho(int px,int py) const
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    Real half_ortho_view_width = (m_Ortho_ViewWidth / 2.0f) / m_ZoomFactor;
    Real half_ortho_view_height = half_ortho_view_width / m_Aspect;
    glOrtho( -half_ortho_view_width, half_ortho_view_width,
             -half_ortho_view_height, half_ortho_view_height,
             m_NearClip, m_FarClip );

    gluPickMatrix(px,viewport[3]-py,5,5,viewport);

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt( m_Obs.x(),m_Obs.y(),m_Obs.z(),
               m_Vrp.x(),m_Vrp.y(),m_Vrp.z(),
               m_UpVec.x(),m_UpVec.y(),m_UpVec.z() );
}

void GLCamera::PrepareToShootViewport() const
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    // (0,0) -> (width,height)
    glOrtho( 0, m_ViewportWidth,
             0, m_ViewportHeight,
             0, 1 ); //Clipping??
    glMatrixMode( GL_MODELVIEW );
    //glLoadIdentity();
    // Invert Y axis, and move center to the top of the viewport
    double matrix_mirror_y[] = { 1,  0, 0, 0,
                                 0, -1, 0, 0,
                                 0,  0, 1, 0,
                                 0,  double(m_ViewportHeight-1), 0, 1 };
    glLoadMatrixd( matrix_mirror_y );
}

//---- Utility Methods

/*! Computes World pos of a given window point.

\note OpenGL viewport pos starts at lower-left window corner (0,0) and ends at
upper-right corner (resolution.x,resolution.y)
Window-system coords, on the contrary, typically increase downwards,
and therefore we correct the projected viewport position Y coord to match this.
*/
Vec3 GLCamera::UnProject( const Vec3 &p3d_window ) const
{
    // IMPORTANT to set matrices before UnProject
    PrepareToShoot();

    GLint viewport[4];
    GLdouble mvmatrix[16], projmatrix[16];

    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);

    GLdouble res[3];
    gluUnProject( p3d_window.x(), GLdouble(viewport[3]-1) - p3d_window.y(), p3d_window.z(),
                  mvmatrix, projmatrix, viewport,
                  &res[0], &res[1], &res[2] );

    return Vec3( res[0], res[1], res[2] );
}

/*! Computes Window pos of a given World point.

\note OpenGL viewport pos starts at lower-left window corner (0,0) and ends at
upper-right corner (resolution.x,resolution.y)
Window-system coords, on the contrary, typically increase downwards,
and therefore we correct the projected viewport position Y coord to match this.
*/
Vec3 GLCamera::Project( const Vec3 &p3d_world ) const
{
    // IMPORTANT to set matrices before Project
    PrepareToShoot();

    GLint viewport[4];
    GLdouble mvmatrix[16], projmatrix[16];

    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);

    GLdouble p3d_viewport[3];
    gluProject( p3d_world.x(), p3d_world.y(), p3d_world.z(),
                mvmatrix, projmatrix, viewport,
                &p3d_viewport[0], &p3d_viewport[1], &p3d_viewport[2] );

    return Vec3( p3d_viewport[0], (GLdouble(viewport[3]+1) - p3d_viewport[1]), p3d_viewport[2] );
}

void GLCamera::SyncUpToViewDir()
{
    // Project current Up vector to current view plane
    Vec3 view_dir( GetViewDir() );
    m_UpVec -= mal::Dot(m_UpVec,view_dir)*view_dir;
    m_UpVec.Normalize();
}

} } // namespace sfr::gfx
