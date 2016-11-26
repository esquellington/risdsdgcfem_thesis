#ifndef SFR_CORE_GFX_GLRENDERER_H
#define SFR_CORE_GFX_GLRENDERER_H

#include <Safra/core/gfx/IRenderer.h>

namespace sfr { namespace gfx
{

//! Immediate Mode GL Renderer for basic geometric primitives
class GLRenderer: public IRenderer
{
public:
    GLRenderer() {}
    ~GLRenderer() {}

    //! \name Immediate Mode rendering
    //@{
    // Linear shapes
    void DrawPoint( const Vec3 &p, const Style &style );
    void DrawSegment( const Vec3 &p0, const Vec3 &p1, const Style &style );
    void DrawVector( const Vec3 &p0, const Vec3 &v, const Style &style );
    void DrawRefSys( const Vec3 &p, const Mat3x3 &m, Real scale, const Style &style );

    // Planar shapes
    void DrawTriangle( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Style &style );
    void DrawQuad( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Vec3 &p3, const Style &style );
    void DrawRectangle( const Vec3 &p, const Mat3x3 &m, const Vec2 &half_sizes, const Style &style );
    void DrawDisk( const Vec3 &p, const Mat3x3 &m, Real radius, const Style &style );
    void DrawDiskAligned( const Vec3 &p, Real radius, const Style &style, unsigned int num_segments );
    void DrawVibratingDiskAligned( const Vec3 &p, Real radius, Real amplitude, Real freq, Real time,
                             const Style &style, unsigned int num_segments );

    // Volumetric shapes
    void DrawBox( const Vec3 &p, const Mat3x3 &m, const Vec3 &half_sizes, const Style &style );
    void DrawSphere( const Vec3 &p, const Mat3x3 &m, Real radius, const Style &style );

    // Misc stuff
    void DrawLabel( const Vec3 &p, const char *label, const Style &style );
    void DrawGrid( const Vec3 &p,
                   const Vec3 &dir1, const Vec3 &dir2,
                   Real length1, Real length2,
                   int size1, int size2,
                   const Style &style );

    // 2D Stuff
    void DrawPoint2( const Vec2 &p, const Style &style );
    void DrawSegment2( const Vec2 &p0, const Vec2 &p1, const Style &style );
    void DrawVector2( const Vec2 &p0, const Vec2 &v, const Style &style );
    void DrawTriangle2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Style &style );
    void DrawBox2( const Vec2 &p, const Mat2x2 &m, const Vec2 &half_sizes, const Style &style );
    void DrawAABB2( const Vec2 &p, const Vec2 &half_sizes, const Style &style );
    void DrawQuad2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Vec2 &p3, const Style &style );
    void DrawDisk2( const Vec2 &p, Real radius, const Style &style, unsigned int num_segments );
    void DrawRefSys2( const Vec2 &p, const Mat2x2 &m, Real scale, const Style &style );
    void DrawLabel2( const Vec2 &p, const char *label, const Style &style );
    //@}
};

}} // namespace sfr::gfx

#endif // SFR_CORE_GFX_GLRENDERER_H
