#ifndef SFR_CORE_GFX_IRENDERER_H
#define SFR_CORE_GFX_IRENDERER_H

#include <Safra/Config.h>
#include "types.h"

namespace sfr {
namespace gfx {

//! Immediate Mode Renderer for basic geometric primitives
class IRenderer
{
public:
    IRenderer() {}
    virtual ~IRenderer() {}

    //! \name Retained Mode Rendering
    //@{
    virtual bool Render()
    {
        //DrawRefSys( Vec3(0,0,0), Mat3x3::Identity(), 25.0f, Style() );
        return true;
    }
    //@}

    //! \name Immediate Mode rendering
    //@{
    // Linear shapes
    virtual void DrawPoint( const Vec3 &p, const Style &style ) {}
    virtual void DrawSegment( const Vec3 &p0, const Vec3 &p1, const Style &style ) {}
    virtual void DrawVector( const Vec3 &p0, const Vec3 &v, const Style &style ) {}
    virtual void DrawRefSys( const Vec3 &p, const Mat3x3 &m, Real scale, const Style &style ) {}

    // Planar shapes
    virtual void DrawTriangle( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Style &style ) {}
    virtual void DrawQuad( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Vec3 &p3, const Style &style ) {}
    virtual void DrawRectangle( const Vec3 &p, const Mat3x3 &m, const Vec2 &half_sizes, const Style &style ) {}
    virtual void DrawDisk( const Vec3 &p, const Mat3x3 &m, Real radius, const Style &style ) {}
    virtual void DrawDiskAligned( const Vec3 &p, Real radius, const Style &style, unsigned int num_segments ) {}
    virtual void DrawVibratingDiskAligned( const Vec3 &p, Real radius,
                                           Real amplitude, Real freq, Real time,
                                           const Style &style, unsigned int num_segments ) {}

    // Volumetric shapes
    virtual void DrawBox( const Vec3 &p, const Mat3x3 &m, const Vec3 &half_sizes, const Style &style ) {}
    virtual void DrawSphere( const Vec3 &p, const Mat3x3 &m, Real radius, const Style &style ) {}
    virtual void DrawCapsule() {} //todo
    virtual void DrawCylinder() {} //todo

    // Misc
    virtual void DrawLabel( const Vec3 &p, const char *label, const Style &style ) {}
    virtual void DrawGrid( const Vec3 &p,
                           const Vec3 &dir1, const Vec3 &dir2,
                           Real length1, Real length2,
                           int size1, int size2,
                           const Style &style ) {}

    // 2D stuff
    virtual void DrawPoint2( const Vec2 &p, const Style &style ) {}
    virtual void DrawSegment2( const Vec2 &p0, const Vec2 &p1, const Style &style ) {}
    virtual void DrawVector2( const Vec2 &p0, const Vec2 &v, const Style &style ) {}
    virtual void DrawTriangle2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Style &style ) {}
    virtual void DrawBox2( const Vec2 &p, const Mat2x2 &m, const Vec2 &half_sizes, const Style &style ) {}
    virtual void DrawAABB2( const Vec2 &pos_min, const Vec2 &pos_max, const Style &style ) {}
    virtual void DrawQuad2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Vec2 &p3, const Style &style ) {}
    virtual void DrawDisk2( const Vec2 &p, Real radius, const Style &style, unsigned int num_segments ) {}
    virtual void DrawRefSys2( const Vec2 &p, const Mat2x2 &m, Real scale, const Style &style ) {}
    virtual void DrawLabel2( const Vec2 &p, const char *label, const Style &style ) {}
    //@}
};

//! Propagates all render calls to the IRenderer received on construction
/*! Utility class that propagates all render calls by deriving from it.
*/
class RendererPropagated: public IRenderer
{
public:
    RendererPropagated( IRenderer *p_renderer )
    : m_pR( p_renderer )
    {}
    ~RendererPropagated() {}

    IRenderer *GetRenderer() const { return m_pR; }

    //! \name IRenderer Implementation by Propagation
    //@{
    virtual bool Render() { return m_pR->Render(); }

    // Linear shapes
    virtual void DrawPoint( const Vec3 &p, const Style &style ) { m_pR->DrawPoint(p,style); }
    virtual void DrawSegment( const Vec3 &p0, const Vec3 &p1, const Style &style ) { m_pR->DrawSegment(p0,p1,style); }
    virtual void DrawVector( const Vec3 &p0, const Vec3 &v, const Style &style ) { m_pR->DrawVector(p0,v,style); }
    virtual void DrawRefSys( const Vec3 &p, const Mat3x3 &m, Real scale,
                             const Style &style ) { m_pR->DrawRefSys(p,m,scale,style); }

    // Planar shapes
    virtual void DrawTriangle( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2,
                               const Style &style ) { m_pR->DrawTriangle(p0,p1,p2,style); }
    virtual void DrawQuad( const Vec3 &p0, const Vec3 &p1, const Vec3 &p2, const Vec3 &p3,
                           const Style &style ) { m_pR->DrawQuad(p0,p1,p2,p3,style); }
    virtual void DrawRectangle( const Vec3 &p, const Mat3x3 &m, const Vec2 &half_sizes,
                                const Style &style ) { m_pR->DrawRectangle(p,m,half_sizes,style); }
    virtual void DrawDisk( const Vec3 &p, const Mat3x3 &m, Real radius,
                           const Style &style ) { m_pR->DrawDisk(p,m,radius,style); }
    virtual void DrawDiskAligned( const Vec3 &p, Real radius,
                            const Style &style, unsigned int num_segments ) { m_pR->DrawDiskAligned(p,radius,style,num_segments); }
    virtual void DrawVibratingDiskAligned( const Vec3 &p, Real radius,
                                     Real amplitude, Real freq, Real time,
                                     const Style &style,
                                     unsigned int num_segments ) { m_pR->DrawVibratingDiskAligned(p,radius,amplitude,freq,time, style,num_segments); }
    // Volumetric shapes
    virtual void DrawBox( const Vec3 &p, const Mat3x3 &m, const Vec3 &half_sizes,
                          const Style &style ) { m_pR->DrawBox(p,m,half_sizes,style); }
    virtual void DrawSphere( const Vec3 &p, const Mat3x3 &m, Real radius,
                             const Style &style ) { m_pR->DrawSphere(p,m,radius,style); }

    // Misc stuff
    virtual void DrawLabel( const Vec3 &p, const char *label, const Style &style ) { m_pR->DrawLabel(p,label,style); }
    virtual void DrawGrid( const Vec3 &p,
                           const Vec3 &dir1, const Vec3 &dir2,
                           Real length1, Real length2,
                           int size1, int size2,
                           const Style &style ) { m_pR->DrawGrid(p,dir1,dir2,length1,length2,size1,size2,style); }
    // 2D Stuff
    virtual void DrawPoint2( const Vec2 &p, const Style &style ) { m_pR->DrawPoint2(p,style); }
    virtual void DrawSegment2( const Vec2 &p0, const Vec2 &p1, const Style &style ) { m_pR->DrawSegment2(p0,p1,style); }
    virtual void DrawVector2( const Vec2 &p0, const Vec2 &v, const Style &style ) { m_pR->DrawVector2(p0,v,style); }
    virtual void DrawTriangle2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Style &style ) { m_pR->DrawTriangle2(p0,p1,p2,style); }
    virtual void DrawBox2( const Vec2 &p, const Mat2x2 &m, const Vec2 &half_sizes, const Style &style ) { m_pR->DrawBox2(p,m,half_sizes,style); }
    virtual void DrawAABB2( const Vec2 &pos_min, const Vec2 &pos_max, const Style &style ) { m_pR->DrawAABB2(pos_min,pos_max,style); }
    virtual void DrawQuad2( const Vec2 &p0, const Vec2 &p1, const Vec2 &p2, const Vec2 &p3, const Style &s ) { m_pR->DrawQuad2(p0,p1,p2,p3,s); }
    virtual void DrawDisk2( const Vec2 &p, Real radius, const Style &style, unsigned int num_segments ) { m_pR->DrawDisk2(p,radius,style,num_segments); }
    virtual void DrawRefSys2( const Vec2 &p, const Mat2x2 &m, Real scale, const Style &style ) { m_pR->DrawRefSys2(p,m,scale,style); }
    virtual void DrawLabel2( const Vec2 &p, const char *label, const Style &style ) { m_pR->DrawLabel2(p,label,style); }
    //@}

private:
    IRenderer *m_pR;
};

} } // namespace sfr::gfx

#endif // SFR_CORE_GFX_IRENDERER_H
