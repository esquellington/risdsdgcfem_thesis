#include <Safra/core/gfx/GLRenderer.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <Mal/GConversion.h>
#include <string.h> //TEMPORAL: For strlen()

namespace sfr { namespace gfx
{

//---- Utility methods
inline void ComputeGLMatrixf( float gl_matrix[16], const Vec3 &p, const Mat3x3 &m )
{
    gl_matrix[0] = m(0,0); gl_matrix[4] = m(0,1); gl_matrix[8]  = m(0,2); gl_matrix[12] = p[0];
    gl_matrix[1] = m(1,0); gl_matrix[5] = m(1,1); gl_matrix[9]  = m(1,2); gl_matrix[13] = p[1];
    gl_matrix[2] = m(2,0); gl_matrix[6] = m(2,1); gl_matrix[10] = m(2,2); gl_matrix[14] = p[2];
    gl_matrix[3] = 0;      gl_matrix[7] = 0;      gl_matrix[11] = 0;      gl_matrix[15] = 1;
}
inline void ComputeGLMatrixf( float gl_matrix[16], const Vec3 &p, const Quat &q )
{
    Mat3x3 m = mal::GRotation3x3_From(q);
    ComputeGLMatrixf( gl_matrix, p, m );
}
inline void ComputeGLMatrixf( float gl_matrix[16], const Vec2 &p, const Mat2x2 &m )
{
    gl_matrix[0] = m(0,0); gl_matrix[4] = m(0,1); gl_matrix[8]  = 0; gl_matrix[12] = p[0];
    gl_matrix[1] = m(1,0); gl_matrix[5] = m(1,1); gl_matrix[9]  = 0; gl_matrix[13] = p[1];
    gl_matrix[2] = 0;      gl_matrix[6] = 0;      gl_matrix[10] = 1; gl_matrix[14] = 0;
    gl_matrix[3] = 0;      gl_matrix[7] = 0;      gl_matrix[11] = 0; gl_matrix[15] = 1;
}


//---- Immediate Mode rendering
//-- Linear shapes
void GLRenderer::DrawPoint( const Vec3 &p, const Style &style )
{
    glPointSize( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    glBegin( GL_POINTS );
      float array3f[3];
      p.ToArray(array3f);
      glVertex3fv( array3f );
    glEnd();
}
void GLRenderer::DrawSegment( const Vec3 &p0, const Vec3 &p1, const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    glBegin( GL_LINES );
      glVertex3f( p0[0], p0[1], p0[2] );
      glVertex3f( p1[0], p1[1], p1[2] );
    glEnd();
}
void GLRenderer::DrawVector( const Vec3 &p0, const Vec3 &v, const Style &style )
{
    DrawSegment( p0, p0+v, style );
    // Arrow head ENDS at vector length
    Vec3f v0( mal::SafeNormalized( v ) );
    Vec3f v1( mal::SafeNormalized( mal::Cross(v0,Vec3(0,0,1)) ) );
    float length( mal::Norm(v) );
    Real size( mal::Clamp(0.25f*length,0.003f,0.01f) );
    DrawTriangle( p0 + v,
                  p0 + v - size*v0 + 0.5f*size*v1,
                  p0 + v - size*v0 - 0.5f*size*v1,
                  Style(style.m_Color,style.m_PenSize,Style::eSolid) );
}
void GLRenderer::DrawRefSys( const Vec3 &p, const Mat3x3 &m, Real scale, const Style &style )
{
    glLineWidth( style.m_PenSize );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
      glMultMatrixf(transf);
      glBegin( GL_LINES );
        glColor3f(1.0,0.0,0.0);
        glVertex3f(0.0,0.0,0.0);
        glVertex3f(scale,0.0,0.0);

        glColor3f(0.0,1.0,0.0);
        glVertex3f(0.0,0.0,0.0);
        glVertex3f(0.0,scale,0.0);

        glColor3f(0.0,0.0,1.0);
        glVertex3f(0.0,0.0,0.0);
        glVertex3f(0.0,0.0,scale);
      glEnd();
    glPopMatrix();
}

//-- Planar shapes
void GLRenderer::DrawTriangle( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Style &style )
{
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glBegin( GL_TRIANGLES );
          glVertex3f( p0[0], p0[1], p0[2] );
          glVertex3f( p1[0], p1[1], p1[2] );
          glVertex3f( p2[0], p2[1], p2[2] );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
          glVertex3f( p0[0], p0[1], p0[2] );
          glVertex3f( p1[0], p1[1], p1[2] );
          glVertex3f( p2[0], p2[1], p2[2] );
          glVertex3f( p0[0], p0[1], p0[2] );
        glEnd();
    }
}

void GLRenderer::DrawQuad( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Vec3 &p3, const Style &style )
{
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glBegin( GL_QUADS );
          glVertex3f( p0[0], p0[1], p0[2] );
          glVertex3f( p1[0], p1[1], p1[2] );
          glVertex3f( p2[0], p2[1], p2[2] );
          glVertex3f( p3[0], p3[1], p3[2] );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
          glVertex3f( p0[0], p0[1], p0[2] );
          glVertex3f( p1[0], p1[1], p1[2] );
          glVertex3f( p2[0], p2[1], p2[2] );
          glVertex3f( p3[0], p3[1], p3[2] );
          glVertex3f( p0[0], p0[1], p0[2] );
        glEnd();
    }
}

void GLRenderer::DrawRectangle( const Vec3 &p, const Mat3x3 &m, const Vec2 &half_sizes, const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
      glMultMatrixf(transf);
      if( style.m_Flags.Test(Style::eSolid) )
      {
          glBegin( GL_TRIANGLE_FAN );
            glVertex3f(0,0,0);
            glVertex3f( +half_sizes.x(), +half_sizes.y(), 0 );
            glVertex3f( -half_sizes.x(), +half_sizes.y(), 0 );
            glVertex3f( -half_sizes.x(), -half_sizes.y(), 0 );
            glVertex3f( +half_sizes.x(), -half_sizes.y(), 0 );
            glVertex3f( +half_sizes.x(), +half_sizes.y(), 0 );
          glEnd();
      }
      if( style.m_Flags.Test(Style::eWire) )
      {
          glBegin( GL_LINE_STRIP );
            glVertex3f( +half_sizes.x(), +half_sizes.y(), 0 );
            glVertex3f( -half_sizes.x(), +half_sizes.y(), 0 );
            glVertex3f( -half_sizes.x(), -half_sizes.y(), 0 );
            glVertex3f( +half_sizes.x(), -half_sizes.y(), 0 );
            glVertex3f( +half_sizes.x(), +half_sizes.y(), 0 );
          glEnd();
      }
    glPopMatrix();
}

void GLRenderer::DrawDisk( const Vec3 &p, const Mat3x3 &m, Real radius, const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
      glMultMatrixf(transf);
      int num_segments = 10;
      Real angle_step( mal::TwoPi<Real>()/num_segments );
      if( style.m_Flags.Test(Style::eSolid) )
      {
          glBegin( GL_TRIANGLE_FAN );
          glVertex3f(0,0,0);
          for( int i=0; i<num_segments+1; i++ )
              glVertex3f(radius*mal::Cos(angle_step*i),radius*mal::Sin(angle_step*i),0);
          glEnd();
      }
      if( style.m_Flags.Test(Style::eWire) )
      {
          glBegin( GL_LINE_STRIP );
          for( int i=0; i<num_segments+1; i++ )
              glVertex3f(radius*mal::Cos(angle_step*i),radius*mal::Sin(angle_step*i),0);
          glEnd();
      }
    glPopMatrix();
}

//! Disk discretized into num_segments segments, starting at the topmost vertex!
void GLRenderer::DrawDiskAligned( const Vec3 &p, Real radius, const Style &style, unsigned int num_segments )
{
    glColor4fv( style.m_Color.AsArray() );
    Real angle_step( mal::TwoPi<Real>()/num_segments );
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glBegin( GL_TRIANGLE_FAN );
        glVertex3f(p[0],p[1],p[2]);
        for( unsigned int i=0; i<num_segments+1; i++ )
            glVertex3f( p[0]+radius*mal::Sin(angle_step*i),
                        p[1]+radius*mal::Cos(angle_step*i),
                        p[2] );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
        for( unsigned int i=0; i<num_segments+1; i++ )
            glVertex3f( p[0]+radius*mal::Sin(angle_step*i),
                        p[1]+radius*mal::Cos(angle_step*i),
                        p[2] );
        glEnd();
    }
}

void GLRenderer::DrawVibratingDiskAligned( const Vec3 &p, Real radius, Real amplitude, Real freq, Real time,
                                           const Style &style, unsigned int num_segments )
{
    glColor4fv( style.m_Color.AsArray() );
    Real angle_step( mal::TwoPi<Real>()/num_segments );
    Real integer_freq( mal::Ceil(freq) ); //Integer-freq guarantees periodic wave in domain [0,2Pi]
    //float sin_freq_x_time = mal::Sin(freq*time);
    Real sin_freq_x_time( mal::Sin(freq*time) );
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glBegin( GL_TRIANGLE_FAN );
        glVertex3f(p[0],p[1],p[2]);
        for( unsigned int i=0; i<num_segments+1; i++ )
        {
            Real angle( angle_step*i );
            Real r( radius + amplitude * Real(0.5f) * ( Real(1) + mal::Sin(integer_freq*angle) * sin_freq_x_time ) );
            glVertex3f( p[0]+r*mal::Cos(angle),
                        p[1]+r*mal::Sin(angle),
                        p[2] );
        }
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
        for( unsigned int i=0; i<num_segments+1; i++ )
        {
            Real angle( angle_step*i );
            Real r( radius + amplitude * Real(0.5f) * ( Real(1) + mal::Sin(integer_freq*angle) * sin_freq_x_time ) );
            glVertex3f( p[0]+r*mal::Cos(angle),
                        p[1]+r*mal::Sin(angle),
                        p[2] );
        }
        glEnd();
    }
}

//-- Volumetric shapes
void GLRenderer::DrawBox( const Vec3 &p, const Mat3x3 &m, const Vec3 &half_sizes, const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
      glMultMatrixf(transf);
      glScalef( 2.0f*half_sizes[0], 2.0f*half_sizes[1], 2.0f*half_sizes[2] );
      if( style.m_Flags.Test(Style::eSolid) )
          glutSolidCube( 1.0f );
      if( style.m_Flags.Test(Style::eWire) )
          glutWireCube( 1.0f );
    glPopMatrix();
}
void GLRenderer::DrawSphere( const Vec3 &p, const Mat3x3 &m, Real radius, const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
      glMultMatrixf(transf);
      if( style.m_Flags.Test(Style::eSolid) )
          glutSolidSphere( radius, 10, 10 );
      if( style.m_Flags.Test(Style::eWire) )
          glutWireSphere( radius, 10, 10 );
    glPopMatrix();
}

// Misc stuff
void GLRenderer::DrawLabel( const Vec3 &p, const char *label, const Style &style )
{
    void *glut_bitmap_type;
    if( style.m_PenSize <= 8.0f ) glut_bitmap_type = GLUT_BITMAP_8_BY_13;
    else if( style.m_PenSize <= 9.0f ) glut_bitmap_type = GLUT_BITMAP_9_BY_15;
    else if( style.m_PenSize <= 10.0f ) glut_bitmap_type = GLUT_BITMAP_HELVETICA_10;
    else if( style.m_PenSize <= 12.0f ) glut_bitmap_type = GLUT_BITMAP_HELVETICA_12;
    else if( style.m_PenSize <= 18.0f ) glut_bitmap_type = GLUT_BITMAP_HELVETICA_18;
    else glut_bitmap_type = GLUT_BITMAP_TIMES_ROMAN_24;

    glColor4fv( style.m_Color.AsArray() );
    glRasterPos3f( p[0], p[1], p[2] );
    for( unsigned int i=0; i < strlen(label); i++ )
        glutBitmapCharacter(glut_bitmap_type, label[i]);

    /* Available fonts (Freeglut)
      GLUT_STROKE_ROMAN
      GLUT_STROKE_MONO_ROMAN
      GLUT_BITMAP_9_BY_15
      GLUT_BITMAP_8_BY_13
      GLUT_BITMAP_TIMES_ROMAN_10
      GLUT_BITMAP_TIMES_ROMAN_24
      GLUT_BITMAP_HELVETICA_10
      GLUT_BITMAP_HELVETICA_12
      GLUT_BITMAP_HELVETICA_18
    */

    /* This is howto draw Viewport coords returned by gluProject(), which are (0,0) at the lower-left corner and
       (resolution_x,resolution_y) at the upper-right corner

    GLint viewport[4];
    GLdouble mvmatrix[16], projmatrix[16];

    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);

    GLdouble pos_viewport[3];
    gluProject( p[0], p[1], p[2],
                mvmatrix, projmatrix, viewport,
                &pos_viewport[0], &pos_viewport[1], &pos_viewport[2] );

    char str[128];
    sprintf(str,"(%f,%f,%f)", pos_viewport[0], pos_viewport[1], pos_viewport[2] );
    glRasterPos3f( p[0], p[1], p[2] );
    for( int i=0; i < strlen(str); i++ )
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, str[i]);
    */
}

void GLRenderer::DrawGrid( const Vec3 &p,
                           const Vec3 &dir1, const Vec3 &dir2,
                           Real length1, Real length2,
                           int size1, int size2,
                           const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    glBegin( GL_LINES );
      float array3f[3];
      for( int i=0; i<size1+1; i++ )
      {
          Vec3 p0 = p + Real(length1)*Real(float(i)/size1)*dir1;
          Vec3 p1 = p0 + length2*dir2;
          p0.ToArray( array3f );
          glVertex3fv( array3f );
          p1.ToArray( array3f );
          glVertex3fv( array3f );
      }
      for( int i=0; i<size2+1; i++ )
      {
          Vec3 p0 = p + Real(length2)*Real(float(i)/size2)*dir2;
          Vec3 p1 = p0 + length1*dir1;
          p0.ToArray( array3f );
          glVertex3fv( array3f );
          p1.ToArray( array3f );
          glVertex3fv( array3f );
      }
    glEnd();
}

//---- 2D Stuff
void GLRenderer::DrawPoint2( const Vec2 &p, const Style &style )
{
    glPointSize( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    glBegin( GL_POINTS );
      glVertex2f( p[0], p[1] );
    glEnd();
}
void GLRenderer::DrawSegment2( const Vec2 &p0, const Vec2 &p1, const Style &style )
{
    glLineWidth( style.m_PenSize );
    glColor4fv( style.m_Color.AsArray() );
    glBegin( GL_LINES );
      glVertex2f( p0[0], p0[1] );
      glVertex2f( p1[0], p1[1] );
    glEnd();
}

void GLRenderer::DrawVector2( const Vec2 &p0, const Vec2 &v, const Style &style )
{
    DrawSegment2(p0,p0+v,style);
    DrawPoint2( p0+v, Style(style.m_Color,3.0f*style.m_PenSize) );
}

void GLRenderer::DrawTriangle2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Style &style )
{
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glBegin( GL_TRIANGLES );
          glVertex2f( p0[0], p0[1] );
          glVertex2f( p1[0], p1[1] );
          glVertex2f( p2[0], p2[1] );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
          glVertex2f( p0[0], p0[1] );
          glVertex2f( p1[0], p1[1] );
          glVertex2f( p2[0], p2[1] );
          glVertex2f( p0[0], p0[1] );
        glEnd();
    }
}

void GLRenderer::DrawBox2( const Vec2 &p, const Mat2x2 &m, const Vec2 &half_sizes, const Style &style )
{
    glLineWidth( style.m_PenSize );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
    {
        glMultMatrixf(transf);
        if( style.m_Flags.Test(Style::eSolid) )
        {
            glColor4fv( style.m_Color.AsArray() );
            glBegin( GL_QUADS );
            glVertex2f( -half_sizes[0], -half_sizes[1] );
            glVertex2f( +half_sizes[0], -half_sizes[1] );
            glVertex2f( +half_sizes[0], +half_sizes[1] );
            glVertex2f( -half_sizes[0], +half_sizes[1] );
            glEnd();
        }
        if( style.m_Flags.Test(Style::eWire) )
        {
            glColor4fv( style.m_Color.AsArray() );
            glLineWidth( style.m_PenSize );
            glBegin( GL_LINE_STRIP );
            glVertex2f( -half_sizes[0], -half_sizes[1] );
            glVertex2f( +half_sizes[0], -half_sizes[1] );
            glVertex2f( +half_sizes[0], +half_sizes[1] );
            glVertex2f( -half_sizes[0], +half_sizes[1] );
            glVertex2f( -half_sizes[0], -half_sizes[1] );
            glEnd();
        }
    }
    glPopMatrix();
}

void GLRenderer::DrawAABB2( const Vec2 &pos_min, const Vec2 &pos_max, const Style &style )
{
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glBegin( GL_QUADS );
          glVertex2f( pos_min[0], pos_min[1] );
          glVertex2f( pos_max[0], pos_min[1] );
          glVertex2f( pos_max[0], pos_max[1] );
          glVertex2f( pos_min[0], pos_max[1] );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
          glVertex2f( pos_min[0], pos_min[1] );
          glVertex2f( pos_max[0], pos_min[1] );
          glVertex2f( pos_max[0], pos_max[1] );
          glVertex2f( pos_min[0], pos_max[1] );
          glVertex2f( pos_min[0], pos_min[1] );
        glEnd();
    }
}

void GLRenderer::DrawQuad2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Vec2 &p3, const Style &style )
{
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glBegin( GL_QUADS );
          glVertex2f( p0[0], p0[1] );
          glVertex2f( p1[0], p1[1] );
          glVertex2f( p2[0], p2[1] );
          glVertex2f( p3[0], p3[1] );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glColor4fv( style.m_Color.AsArray() );
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
          glVertex2f( p0[0], p0[1] );
          glVertex2f( p1[0], p1[1] );
          glVertex2f( p2[0], p2[1] );
          glVertex2f( p3[0], p3[1] );
          glVertex2f( p0[0], p0[1] );
        glEnd();
    }
}

void GLRenderer::DrawDisk2( const Vec2 &p, Real radius, const Style &style, unsigned int num_segments )
{
    glColor4fv( style.m_Color.AsArray() );
    Real angle_step( mal::TwoPi<Real>()/num_segments );
    if( style.m_Flags.Test(Style::eSolid) )
    {
        glBegin( GL_TRIANGLE_FAN );
        glVertex2f(p[0],p[1]);
        for( unsigned int i=0; i<num_segments+1; i++ )
            glVertex2f( p[0]+radius*mal::Sin(angle_step*i),
                        p[1]+radius*mal::Cos(angle_step*i) );
        glEnd();
    }
    if( style.m_Flags.Test(Style::eWire) )
    {
        glLineWidth( style.m_PenSize );
        glBegin( GL_LINE_STRIP );
        for( unsigned int i=0; i<num_segments+1; i++ )
            glVertex2f( p[0]+radius*mal::Sin(angle_step*i),
                        p[1]+radius*mal::Cos(angle_step*i) );
        glEnd();
    }
}

void GLRenderer::DrawRefSys2( const Vec2 &p, const Mat2x2 &m, Real scale, const Style &style )
{
    glLineWidth( style.m_PenSize );
    float transf[16];
    ComputeGLMatrixf( transf, p, m );
    glPushMatrix();
      glMultMatrixf(transf);
      glBegin( GL_LINES );
        glColor3f(1.0,0.0,0.0);
        glVertex2f(0.0,0.0);
        glVertex2f(scale,0.0);
        glColor3f(0.0,1.0,0.0);
        glVertex2f(0.0,0.0);
        glVertex2f(0.0,scale);
      glEnd();
    glPopMatrix();
}

void GLRenderer::DrawLabel2( const Vec2 &p, const char *label, const Style &style )
{
    DrawLabel( Vec3(p[0],p[1],0), label, style );
}

}} // namespace sfr::gfx
