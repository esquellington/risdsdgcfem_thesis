#ifndef SAFRA_CONFIG_H
#define SAFRA_CONFIG_H

#include <boost/static_assert.hpp>

#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GMat.h>
#include <Mal/GTransform.h>
#include <Mal/GMatUtils.h>
#include <Mal/GQuatUtils.h>

#include <pla_types.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define SFR_ASSERT( x )
#  define SFR_LOG( x, ... )
#  define SFR_LOG_WARNING( x, ... )
#  define SFR_LOG_ERROR( x, ... )
#  define SFR_LOG_ASSERT( c, x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define SFR_ASSERT( x ) assert(x)
#  define SFR_LOG( x, ... ) { printf("<SFR LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define SFR_LOG_WARNING( x, ... ) { printf("<SFR WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define SFR_LOG_ERROR( x, ... ) { printf("<SFR ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define SFR_LOG_ASSERT( c, x, ... ) { if(!(c)) { printf("<SFR ASSERT FAILED> "); printf( x, ##__VA_ARGS__ ); printf("\n"); assert(c); } }
#endif

namespace sfr
{

typedef float Real;
typedef mal::GVec<Real,3> Vec3;
typedef mal::GVec<Real,2> Vec2;
typedef mal::GQuat<Real> Quat;
typedef mal::GTransform<Real,3> Transform2;
typedef mal::GTransform<Real,3> Transform3;
typedef mal::GMat<Real,2,2> Mat2x2;
typedef mal::GMat<Real,3,3> Mat3x3;
typedef mal::GMat<Real,4,4> Mat4x4;

}

#endif //SAFRA_CONFIG_H
