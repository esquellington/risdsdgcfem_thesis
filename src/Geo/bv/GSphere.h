#ifndef GEO_BV_GSPHERE_H
#define GEO_BV_GSPHERE_H

#include <Geo/bv/BoundingVolume.h>
#include "stats.h"

namespace geo {
namespace bv {

template <unsigned D>
class GSphere: public GBoundingVolumeD<D>
{
public:
    enum EConstants { cDimension = D };
    typedef mal::GVec<Real,D> vec_type;

public:
    finline GSphere() : GBoundingVolumeD<D>((D==2)?eBV_Sphere2:eBV_Sphere3) {}
    finline explicit GSphere( const vec_type& pos ) : GBoundingVolumeD<D>((D==2)?eBV_Sphere2:eBV_Sphere3)
                                                    , m_Pos(pos), m_Radius(0) {}
    finline GSphere( const vec_type& pos, Real radius ) : GBoundingVolumeD<D>((D==2)?eBV_Sphere2:eBV_Sphere3)
                                                        , m_Pos(pos), m_Radius(radius) {}
    finline ~GSphere() {}

    finline GSphere& SetPos( const vec_type& pos ) { m_Pos = pos; IBoundingVolume::Touch(); return *this; }
    finline GSphere& SetRadius( Real radius ) { m_Radius = radius; IBoundingVolume::Touch(); return *this; }
    finline GSphere& SetPosRadius( const vec_type& pos, Real radius ) { m_Pos = pos; m_Radius = radius; IBoundingVolume::Touch(); return *this; }

    finline const vec_type& GetPos() const { return m_Pos; }
    finline Real GetRadius() const { return m_Radius; }

    finline GSphere& Extend( Real thickness ) { m_Radius += thickness; IBoundingVolume::Touch(); return *this; }

    finline GSphere& Merge( const vec_type& point, Real radius )
        {
            //\note Code adapted from RTCD pg.268
            vec_type d = point - m_Pos;
            Real dist_sq = mal::NormSq( d );
            if( mal::Sq(m_Radius - radius) >= dist_sq )
            {
                // a sphere completely inside the other
                if( m_Radius < radius )
                {
                    // this sphere completely inside the other
                    m_Pos = point;
                    m_Radius = radius;
                    IBoundingVolume::Touch();
                }
            }
            else
            {
                // spheres partially overlapping or disjoint
                Real dist = mal::Sqrt(dist_sq);
                Real r0 = m_Radius;
                m_Radius = Real(0.5f)*(dist + r0 + radius);
                if( dist > mal::Epsilon<Real>() ) m_Pos += ((m_Radius-r0) / dist) * d; //\todo USE CONTEXT LENGTH EPSILON
                IBoundingVolume::Touch();
            }
            return *this;
        }
    finline GSphere& Merge( const vec_type& point ) { return Merge(point,0); }
    finline GSphere& Merge( const GSphere& sphere ) { return Merge(sphere.m_Pos,sphere.m_Radius); }

public:
    vec_type m_Pos;
    Real m_Radius;
};

//! \name Usual dimensions
//@{
typedef GSphere<2> Sphere2;
typedef GSphere<3> Sphere3;

template <> struct bv_type_of<Sphere2> { enum { value = eBV_Sphere2 }; };
template <> struct bv_type_of<Sphere3> { enum { value = eBV_Sphere3 }; };
//@}

//Maximum distance between 2 points on the BV
template <unsigned D>
finline Real ComputeLargestExtent( const GSphere<D>& sphere )
{
    return 2*sphere.m_Radius;
}

//\todo Ugly way to support D==2 and D==3...
template <unsigned D>
finline Real ComputeVolume( const GSphere<D>& sphere )
{
    if( D == 2 )
        return mal::Pi<Real>()*mal::Sq(sphere.m_Radius);
    else /*D==3*/
        return Real(4.0f/3.0f)*mal::Pi<Real>()*mal::Sq(sphere.m_Radius)*sphere.m_Radius;
}

template <unsigned D>
finline Real ComputeVolume_Fast( const GSphere<D>& sphere ) { return ComputeVolume(sphere); }

template <unsigned D>
finline bool TestOverlap( const GSphere<D>& sphere1, const GSphere<D>& sphere2 )
{
    GEO_BV_STAT_INC( testoverlap.m_GSphere_GSphere );
    return ( mal::NormSq(sphere1.m_Pos-sphere2.m_Pos) <= mal::Sq(sphere1.m_Radius+sphere2.m_Radius) );
}

template <unsigned D>
finline bool TestRay( const GSphere<D>& sphere,
                      const mal::GVec<Real,D>& ray_pos, const mal::GVec<Real,D>& ray_dir, const mal::GInterval<Real>& interval )
{
#ifdef __GEO_ENABLE_STATS
    // if( D == 2 ) { GEO_BV_STAT_INC( testray.m_Sphere2 ); }
    // else { GEO_BV_STAT_INC( testray.m_Sphere2 ); }
#endif
    const mal::GVec<Real,D> ray_pos_rel( ray_pos - sphere.m_Pos );
    Real dot_pos_rel_dir( mal::Dot(ray_pos_rel,ray_dir) );
    Real dir_norm_sq( ray_dir.NormSq() );
    Real sqrt_arg( mal::Sq(dot_pos_rel_dir) - dir_norm_sq * ( ray_pos_rel.NormSq() - mal::Sq(sphere.m_Radius) ) );
    if( sqrt_arg >= Real(0) )
    {
        Real root( mal::Sqrt(sqrt_arg) );
        Real rcp_dir_norm_sq( mal::Rcp(dir_norm_sq) );
        Real lambda0( (-dot_pos_rel_dir-root)*rcp_dir_norm_sq );
        Real lambda1( (-dot_pos_rel_dir+root)*rcp_dir_norm_sq );
        if( lambda0 <= interval.Max() && lambda1 >= interval.Min() ) //implicitly lambda0 < lambda1
            return true;
    }
    return false;
}


}} // namespace geo::bv

#endif // GEO_BV_GSPHERE_H
