#ifndef GEO_NP_CLOSEST_POINTS_H
#define GEO_NP_CLOSEST_POINTS_H

#include <Geo/Config.h>
#include "Context.h"

namespace geo {
namespace np {

/*! Computes point on the segment closest to point. Supports degenerate 0-length segments */
template< unsigned D >
mal::GVec<Real,D> GClosestPoint_Point_Segment( const mal::GVec<Real,D> &point,
                                               const mal::GVec<Real,D> &p0, const mal::GVec<Real,D> &p1 )
{
    mal::GVec<Real,D> d( p1 - p0 );
    mal::GVec<Real,D> r( point - p0 );
    Real dot_d_r( mal::Dot( d, r ) );
    // Determine if point is "behind", "after" or "inbetween" (p0,p1)
    if( dot_d_r < g_pDefaultContext->m_Epsilon_Length ) return p0; //Includes 0-lenght segment case
    else if( mal::Dot( point - p1, d ) > g_pDefaultContext->m_Epsilon_Length ) return p1;
    else return p0 + ( dot_d_r / mal::NormSq(d) ) * d;
}

/*! Taken from RTCD 5.1.9 pag 149, but ignoring the degenerate
  (zero-length segments) case.

  Computes closest points P' and Q' of SP(s)=P0+s*(P1-P0) and
  SQ(t)=Q0+t*(Q1-Q0), returning s and t.
*/
template< unsigned D >
void GClosestPoints_Segment_Segment( const mal::GVec<Real,D> &p0, const mal::GVec<Real,D> &p1,
                                     const mal::GVec<Real,D> &q0, const mal::GVec<Real,D> &q1,
                                     Real &lambda_p, Real &lambda_q )
{
    mal::GVec<Real,D> dp(p1-p0);
    mal::GVec<Real,D> dq(q1-q0);
    
    mal::GVec<Real,D> r(p0 - q0);
    Real sq_norm_dp( dp.NormSq() ); //(a) Squared length of segment SP, always nonnegative 
    Real sq_norm_dq( dq.NormSq() ); //(e) Squared length of segment SQ, always nonnegative
    Real dot_dp_r( mal::Dot(dp,r) ); //(c)
    Real dot_dq_r( mal::Dot(dq,r) ); //(f)
    
    // The general nondegenerate case starts here
    Real dot_dp_dq( mal::Dot(dp,dq) ); //(b)
    Real denom( sq_norm_dp*sq_norm_dq - dot_dp_dq*dot_dp_dq ); // Always nonnegative
    
    // If segments not parallel, compute closest point on SP to SQ, and
    // clamp to segment SP. Else pick arbitrary s (here 0)
    if( denom != Real(0) )
        lambda_p = mal::Clamp01( Real(dot_dp_dq*dot_dq_r - dot_dp_r*sq_norm_dq) / denom );
    else
        lambda_p = Real(0);

    // Compute point on SQ closest to SP(s) using
    // t = Dot((P0+DP*s)-Q0,DQ) / Dot(DQ,DQ) = (b*s + f) / e
    lambda_q = (dot_dp_dq*lambda_p + dot_dq_r) / sq_norm_dq;

    // If t in [0,1] done. Else clamp t, recompute s for the new value
    // of t using s = Dot((Q0+DQ*t)-P0,DQ) / Dot(DP,DP) = (t*b - c) / a
    // and clamp s to [0, 1]
    if ( lambda_q < Real(0) )
    {
        lambda_q = Real(0);
        lambda_p = mal::Clamp01( -dot_dp_r / sq_norm_dp );
    }
    else if ( lambda_q > Real(1) )
    {
        lambda_q = Real(1);
        lambda_p = mal::Clamp01( (dot_dp_dq - dot_dp_r) / sq_norm_dp );
    }
}

}} //namespace geo::np

#endif // GEO_NP_CLOSEST_POINTS_H
