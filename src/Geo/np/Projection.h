#ifndef GEO_NP_PROJECTION_H
#define GEO_NP_PROJECTION_H

#include <Geo/Config.h>

namespace geo {
namespace np {

//----------------------------------------------------------------------------------------
//---- Arbitrary Axis Shape Projections
//----------------------------------------------------------------------------------------

// Compute projection interval between ray and sphere
template< unsigned D >
inline void GProjection_Sphere( const mal::GVec<Real,D> &sphere_pos, Real sphere_radius,
                                const mal::GVec<Real,D> &axis,
                                Interval &interval )
{
    interval.SetCenterHalfSizes( mal::Dot(axis,sphere_pos), sphere_radius );
}

/*\todo
  template< unsigned D >
  inline void GProjection_Box( const mal::GTransform<Real,D> &box_transform, const mal::GVec<Real,D> &box_half_sizes,
  const mal::GVec<Real,D> &axis,
  Interval &interval )
  {
  interval.SetCenterHalfSizes( mal::Dot(axis,box_tra), sphere_radius );
  }
*/

template< unsigned D >
inline void GProjection_LSS( const mal::GVec<Real,D> &lss_p0, const mal::GVec<Real,D> &lss_p1, Real lss_radius,
                             const mal::GVec<Real,D> &axis,
                             Interval &interval )
{
    interval.SetCenterHalfSizes( Real(0.5)*mal::Dot( axis, lss_p0+lss_p1 ),
                                 Real(0.5)*mal::Abs( mal::Dot( axis, lss_p1-lss_p0 ) + lss_radius ) );
}

//\todo AABB....


//----------------------------------------------------------------------------------------
//---- Coordinate Axis Shape Projections
//----------------------------------------------------------------------------------------

template< unsigned D, unsigned C >
inline void GProjection_Sphere( const mal::GVec<Real,D> &sphere_pos, Real sphere_radius,
                                Interval &interval )
{
    GEO_STATIC_ASSERT( C < D );
    interval.SetCenterHalfSizes( sphere_pos[C], sphere_radius );
}

template< unsigned D, unsigned C >
inline void GProjection_LSS( const mal::GVec<Real,D> &lss_p0, const mal::GVec<Real,D> &lss_p1, Real lss_radius,
                             Interval &interval )
{
    GEO_STATIC_ASSERT( C < D );
    if( lss_p0[C] < lss_p1[C] )
        interval.SetMinMax( lss_p0[C]-lss_radius, lss_p1[C]+lss_radius );
    else
        interval.SetMinMax( lss_p1[C]-lss_radius, lss_p0[C]+lss_radius );
}

}} //namespace geo::np

#endif // GEO_NP_PROJECTION_H
