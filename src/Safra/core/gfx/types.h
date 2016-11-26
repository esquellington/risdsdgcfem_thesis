#ifndef SFR_CORE_GFX_TYPES_H
#define SFR_CORE_GFX_TYPES_H

#include <Safra/Config.h>

namespace sfr {
namespace gfx {

//! Simple (R,G,B,A) color class
// struct Color
// {
//     float m_R, m_G, m_B, m_A;
//     Color() : m_R(1), m_G(1), m_B(1), m_A(1) {}
//     Color( float r, float g, float b, float a = 1.0f ) : m_R(r), m_G(g), m_B(b), m_A(a) {}
//     Color( const Vec3f &vec, float a = 1.0f ) { m_R = vec[0]; m_G = vec[1]; m_B = vec[2]; m_A = a; }
//     Color( const Vec4f &vec ) { m_R = vec[0]; m_G = vec[1]; m_B = vec[2]; m_A = vec[3]; }
//     operator const float *() const { return &m_R; }
//     operator Vec4f() const { return Vec4f(m_R,m_G,m_B,m_A); }
// };

typedef Vec4f Color;

struct Style
{
    enum EFlags { eNone = 0, eSolid = 1, eWire = 2, eDefault = eWire };
    Color m_Color;
    float m_PenSize; //!< aka line_width or point_size...
    Flags32 m_Flags;

    Style() : m_Color(1,1,1,1), m_PenSize(1), m_Flags(eDefault) {}
    Style( const Color &color ) : m_Color(color), m_PenSize(1), m_Flags(eDefault) {}
    Style( const Color &color, float pen_size ) : m_Color(color), m_PenSize(pen_size), m_Flags(eDefault) {}
    Style( const Color &color, float pen_size, Flags32 flags ) : m_Color(color), m_PenSize(pen_size), m_Flags(flags) {}
    Style( const Color &color, Flags32 flags ) : m_Color(color), m_PenSize(1), m_Flags(flags) {}
    operator const Color &() const { return m_Color; }
};

} } // namespace sfr::gfx

#endif // SFR_CORE_GFX_TYPES_H
