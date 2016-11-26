#ifndef S2_MS_FEM_ROTATION_H
#define S2_MS_FEM_ROTATION_H

#include "../Config.h"
#include <Mal/GMatDecomposition.h>
#include "TetrahedronElement3.h"

namespace S2 {
namespace ms {
namespace fem {

////////////////////////////////////////////////////////////////
// Rotation
////////////////////////////////////////////////////////////////

inline Mat2x2 Compute_R_QR( const Mat2x2& F ) { return mal::GRotation2x2_GramSchmidtOrthonormalization_YX( F, Mat2x2::Identity() ); }
inline Mat3x3 Compute_R_QR( const Mat3x3& F ) { return mal::GRotation3x3_GramSchmidtOrthonormalization_XYZ( F, Mat3x3::Identity() ); }

inline Mat2x2 Compute_R_PD( const Mat2x2& F, Real detF ) { return mal::GRotation2x2_PolarDecomposition( F, detF ); }
inline Mat3x3 Compute_R_PD( const Mat3x3& F, Real detF ) { return mal::GRotation3x3_PolarDecomposition( F, detF ); }

Mat2x2 Compute_R_SVD( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F );
Mat3x3 Compute_R_SVD( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F );

Mat2x2 Compute_R_PDP_Degenerate( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F,
                                 int32 noc,
                                 const Vec2& p0, const Vec2& p1, const Vec2& p2,
                                 Real area, const Mat2x2& inv_Dm );
Mat3x3 Compute_R_PDP_Degenerate( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F,
                                 TetrahedronElement3::DoC doc,
                                 const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                                 Real volume, const Mat3x3& inv_Dm );

Mat2x2 Compute_R_PDR_Degenerate( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F,
                                 int32 noc,
                                 const Vec2& p0, const Vec2& p1, const Vec2& p2 );
Mat3x3 Compute_R_PDR_Degenerate( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F,
                                 TetrahedronElement3::DoC doc,
                                 const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3 );

Mat2x2 Compute_R_DAPD_Degenerate( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F,
                                  int32 noc,
                                  const Vec2& p0, const Vec2& p1, const Vec2& p2,
                                  Real area, const Mat2x2& inv_Dm,
                                  Real factor_L, Real factor_NL, Real exponent_NL );
Mat3x3 Compute_R_DAPD_Degenerate( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F,
                                  TetrahedronElement3::DoC doc,
                                  const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                                  Real volume, const Mat3x3& inv_Dm,
                                  Real factor_L, Real factor_NL, Real exponent_NL );

Mat2x2 Compute_R_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                       int32 noc,
                       const Vec2& y0, const Vec2& y1, const Vec2& y2,
                       const Mat2x2& invDm, Real element_area );
Mat3x3 Compute_R_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                       TetrahedronElement3::DoC doc,
                       const Vec3& y0, const Vec3& y1, const Vec3& y2, const Vec3& y3,
                       const Mat3x3& invDm, Real element_volume );

////////////////////////////////////////////////////////////////
// Rotation Differential
// \note All methods returns false if dR is 0 (optimization)
//       or infinity (avoid div-by-zero)
////////////////////////////////////////////////////////////////

bool Compute_dR_QR_YX( const Vec2& dy0, const Vec2& dy1, const Vec2& dy2,
                       const Mat2x2& invDm, const Mat2x2& F,
                       Mat2x2& dR );

bool Compute_dR( const Vec2& dx0, const Vec2& dx1, const Vec2& dx2,
                 const Mat2x2& invDm, const Mat2x2& F, const Mat2x2& R, const Mat2x2& Rt,
                 Mat2x2& dR );
bool Compute_dR( const Vec3& dx0, const Vec3& dx1, const Vec3& dx2, const Vec3& dx3,
                 const Mat3x3& invDm, const Mat3x3& F, const Mat3x3& R, const Mat3x3& Rt,
                 Mat3x3& dR );

bool Compute_dR_Numerical( const Vec2& x0, const Vec2& x1, const Vec2& x2,
                           const Vec2& dx0, const Vec2& dx1, const Vec2& dx2,
                           const Mat2x2& invDm, Real degenerate_threshold_det_F,
                           Real h,
                           Mat2x2& dR );
bool Compute_dR_Numerical( const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3,
                           const Vec3& dx0, const Vec3& dx1, const Vec3& dx2, const Vec3& dx3,
                           const Mat3x3& invDm, Real degenerate_threshold_det_F,
                           Real h,
                           Mat3x3& dR );

bool Compute_dR_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                      int32 noc,
                      const Vec2& y0, const Vec2& y1, const Vec2& y2,
                      const Vec2& dy0, const Vec2& dy1, const Vec2& dy2,
                      const Mat2x2& invDm, const Mat2x2& F, const Mat2x2& R, const Mat2x2& Rt, Real element_area,
                      Mat2x2& dR );
bool Compute_dR_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                      TetrahedronElement3::DoC doc,
                      const Vec3& y0, const Vec3& y1, const Vec3& y2, const Vec3& y3,
                      const Vec3& dy0, const Vec3& dy1, const Vec3& dy2, const Vec3& dy3,
                      const Mat3x3& invDm, const Mat3x3& F, const Mat3x3& R, const Mat3x3& Rt, Real element_area,
                      Mat3x3& dR );

bool Compute_dR_DAPD_Numerical( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                                int32 noc,
                                const Vec2& x0, const Vec2& x1, const Vec2& x2,
                                const Vec2& dx0, const Vec2& dx1, const Vec2& dx2,
                                const Mat2x2& invDm, Real element_area,
                                Real h,
                                Mat2x2& dR );
bool Compute_dR_DAPD_Numerical( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                                TetrahedronElement3::DoC doc,
                                const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3,
                                const Vec3& dx0, const Vec3& dx1, const Vec3& dx2, const Vec3& dx3,
                                const Mat3x3& invDm, Real element_volume,
                                Real h,
                                Mat3x3& dR );

}}} //namespace S2::ms::fem

#endif //S2_MS_FEM_ROTATION_H
