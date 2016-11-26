#ifndef S2_CONFIG_H
#define S2_CONFIG_H

#include <Mal/GInterval.h>
#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GMat.h>
#include <Mal/GTransform.h>
#include <Mal/GQuatUtils.h>
#include <Mal/GMatUtils.h>
#include <Mal/GSRV.h>

namespace S2
{
//#define S2_REAL_IS_DOUBLE
#ifdef S2_REAL_IS_DOUBLE
  typedef double Real;
#else
  typedef float Real;
#endif

// Standard Real types
typedef mal::GInterval<Real> Interval;

typedef mal::GVec<Real,3> Vec3;
typedef mal::GVec<Real,3> Point3;
typedef mal::GQuat<Real> Quat;
typedef mal::GMat<Real,3,3> Mat3x3;
typedef mal::GTransform<Real,3> Transform3;

typedef mal::GVec<Real,2> Vec2;
typedef mal::GVec<Real,2> Point2;
typedef mal::GMat<Real,2,2> Mat2x2;
typedef mal::GTransform<Real,2> Transform2;

typedef mal::GSRV<Real> SRV;
//typedef GSRM<Real> SRM;

}

#endif // S2_CONFIG_H
