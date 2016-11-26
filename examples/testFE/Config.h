#ifndef TEST_FE_CONFIG_H
#define TEST_FE_CONFIG_H

// Some useful macros
#ifdef PROFILE_FINAL
#  define APP_ASSERT( x )
#  define APP_LOG( x, ... )
#  define APP_LOG_WARNING( x, ... )
#  define APP_LOG_ERROR( x, ... )
#  define APP_LOG_ASSERT( x, y, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define APP_ASSERT( x ) assert(x)
#  define APP_LOG( x, ... ) { printf("<APP LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_WARNING( x, ... ) { printf("<APP WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_ERROR( x, ... ) { printf("<APP ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_ASSERT( c, x, ... )  { if(!(c)) { printf("<APP ASSERT FAILED> "); printf( x, ##__VA_ARGS__ ); printf("\n"); assert(c); } }
#endif

#include <Mal/GVec.h>
#include <Mal/GMat.h>
#include <Mal/GMatUtils.h>
#include <Mal/GConversion.h>

//#define TEST_FE_REAL_IS_DOUBLE
#ifdef TEST_FE_REAL_IS_DOUBLE
  typedef double Real;
#else
  typedef float Real;
#endif

//Element2D
typedef mal::GVec<Real,2> Vec2r;
typedef mal::GVec<Real,6> Vec6r;
typedef mal::GMat<Real,2,2> Mat2x2r;
typedef mal::GMat<Real,2,6> Mat2x6r;
typedef mal::GMat<Real,3,3> Mat3x3r;
typedef mal::GMat<Real,3,6> Mat3x6r;
typedef mal::GMat<Real,6,6> Mat6x6r;

//Element3D
typedef mal::GVec<Real,3> Vec3r;
typedef mal::GVec<Real,12> Vec12r;
typedef mal::GMat<Real,3,3> Mat3x3r;
typedef mal::GMat<Real,4,4> Mat4x4r;
typedef mal::GMat<Real,6,6> Mat6x6r;
typedef mal::GMat<Real,6,12> Mat6x12r;
typedef mal::GMat<Real,12,12> Mat12x12r;

#define __USE_NOC

#define __USE_GSL
#ifdef __USE_GSL
#  define __USE_GSL_SVD //\see http://www.gnu.org/software/gsl/manual/html_node/Singular-Value-Decomposition.html
#  define __USE_GSL_QR //\see http://www.gnu.org/software/gsl/manual/html_node/QR-Decomposition.html
//#  define __USE_GSL_QRPT //\see http://www.gnu.org/software/gsl/manual/html_node/QR-Decomposition-with-Column-Pivoting.html
#endif

#endif //TEST_FE_CONFIG_H
