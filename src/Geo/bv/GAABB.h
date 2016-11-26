#ifndef GEO_BV_GAABB_H
#define GEO_BV_GAABB_H

#include <Geo/bv/BoundingVolume.h>
#include "stats.h"

namespace geo {
namespace bv {

template <unsigned D>
class GAABB: public GBoundingVolumeD<D>
{
public:
    enum EConstants { cDimension = D };
    typedef mal::GVec<Real,D> vec_type;

public:
    finline GAABB() : GBoundingVolumeD<D>((D==2)?eBV_AABB2:eBV_AABB3) {}
    finline GAABB( const vec_type &pos, const vec_type &half_sizes ) : GBoundingVolumeD<D>((D==2)?eBV_AABB2:eBV_AABB3),
                                                                       m_Pos(pos), m_HalfSizes(half_sizes) {}
    finline explicit GAABB( const vec_type &pos ) : GBoundingVolumeD<D>((D==2)?eBV_AABB2:eBV_AABB3),
                                                    m_Pos(pos), m_HalfSizes(0) {}
    finline ~GAABB() {}

    finline GAABB& Set( const vec_type &pos ) { m_Pos = pos; m_HalfSizes = vec_type::Zero(); IBoundingVolume::Touch(); return *this; }

    finline GAABB& SetMinMax( const vec_type &p0,
                              const vec_type &p1 ) { m_Pos=Real(0.5f)*(p0+p1); m_HalfSizes=Real(0.5f)*(p1-p0); IBoundingVolume::Touch(); return *this; }
    finline GAABB& SetPosHalfSizes( const vec_type &pos,
                                    const vec_type &half_sizes ) { m_Pos = pos; m_HalfSizes = half_sizes; IBoundingVolume::Touch(); return *this; }
    finline void GetMinMax( vec_type &pos0,
                            vec_type &pos1 ) const { pos0 = m_Pos - m_HalfSizes; pos1 = m_Pos + m_HalfSizes; }
    finline vec_type GetMin() const { return m_Pos-m_HalfSizes; }
    finline vec_type GetMax() const { return m_Pos+m_HalfSizes; }
    finline const vec_type &GetPos() const { return m_Pos; }
    finline const vec_type &GetHalfSizes() const { return m_HalfSizes; }

    finline GAABB& Extend( Real thickness ) { for( unsigned int i=0; i<cDimension; i++ ) m_HalfSizes[i] += thickness; IBoundingVolume::Touch(); return *this; }

    finline GAABB& Merge( const vec_type &point )
    {
        vec_type p0,p1;
        GetMinMax(p0,p1);
        for( unsigned int i=0; i<cDimension; i++ )
        {
            p0[i] = mal::Min( p0[i], point[i] );
            p1[i] = mal::Max( p1[i], point[i] );
        }
        SetMinMax(p0,p1); //!< Already touches
        return *this;
    }

    finline GAABB& Merge( const vec_type &point, Real radius )
    {
        vec_type p0,p1;
        GetMinMax(p0,p1);
        for( unsigned int i=0; i<cDimension; i++ )
        {
            p0[i] = mal::Min( p0[i], point[i] - radius );
            p1[i] = mal::Max( p1[i], point[i] + radius );
        }
        SetMinMax(p0,p1); //!< Already touches
        return *this;
    }

    finline GAABB& Merge( const GAABB &aabb )
    {
        vec_type p0,p1,q0,q1;
        GetMinMax(p0,p1);
        aabb.GetMinMax(q0,q1);
        for( unsigned int i=0; i<cDimension; i++ )
        {
            p0[i] = mal::Min( p0[i], q0[i] );
            p1[i] = mal::Max( p1[i], q1[i] );
        }
        SetMinMax(p0,p1); //!< Already touches
        return *this;
    }

public:
    vec_type m_Pos;
    vec_type m_HalfSizes;
};

//! \name Usual dimensions
//@{
typedef GAABB<2> AABB2;
typedef GAABB<3> AABB3;

template <> struct bv_type_of<AABB2> { enum { value = eBV_AABB2 }; };
template <> struct bv_type_of<AABB3> { enum { value = eBV_AABB3 }; };
//@}

//Maximum distance between 2 points on the BV
template <unsigned D>
finline Real ComputeLargestExtent( const GAABB<D>& aabb )
{
    Real half_diag_sq(0);
    for( unsigned int i=0; i<D; i++ ) half_diag_sq += mal::Sq(aabb.m_HalfSizes[i]);
    return 2*mal::Sqrt(half_diag_sq);
}

template <unsigned D>
finline Real ComputeVolume( const GAABB<D>& aabb )
{
    Real volume(1);
    for( unsigned int i=0; i<D; i++ ) volume *= 2*aabb.m_HalfSizes[i];
    return volume;
}
template <unsigned D>
finline Real ComputeVolume_Fast( const GAABB<D>& aabb ) { return ComputeVolume(aabb); }

template <unsigned D>
finline bool TestOverlap( const GAABB<D>& aabb1, const GAABB<D>& aabb2 )
{
    GEO_BV_STAT_INC( testoverlap.m_GAABB_GAABB );
    for( unsigned int i=0; i<D; i++ )
        if( mal::Abs(aabb1.m_Pos[i]-aabb2.m_Pos[i]) > aabb1.m_HalfSizes[i]+aabb2.m_HalfSizes[i] )
            return false;
    return true;
}

template <unsigned D>
finline bool TestRay( const GAABB<D>& aabb,
                      const mal::GVec<Real,D>& ray_pos, const mal::GVec<Real,D>& ray_dir, const mal::GInterval<Real>& interval )
{
#ifdef __GEO_ENABLE_STATS
    // if( D == 2 ) { GEO_BV_STAT_INC( testray.m_AABB2 ); }
    // else { GEO_BV_STAT_INC( testray.m_AABB3 ); }
#endif
    mal::GInterval<Real> overlap_interval( interval );
    mal::GVec<Real,D> ray_pos_rel( ray_pos - aabb.m_Pos );
    for( int it_axis=0; it_axis<int(D); it_axis++ )
    {
        // If parallel to slab, either overlaps for any lambda or for none.
        if( mal::Abs(ray_dir[it_axis]) < 1e-4 ) //\todo g_pDefaultContext->m_Epsilon_Dir )
        {
            if( mal::Abs(ray_pos_rel[it_axis]) > aabb.m_HalfSizes[it_axis] )
                return false;
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            // \todo THIS CODE IS LONG AND UGLY, try to simplify it...
            Real inv_divisor( mal::Rcp(ray_dir[it_axis]) );
            Real lambda0( ( -aabb.m_HalfSizes[it_axis] - ray_pos_rel[it_axis] ) * inv_divisor );
            Real lambda1( (  aabb.m_HalfSizes[it_axis] - ray_pos_rel[it_axis] ) * inv_divisor );
            if( lambda1 < lambda0 )
            {
                Real tmp( lambda0 );
                lambda0 = lambda1;
                lambda1 = tmp;
            }
            // Clip lambda-interval and update first-axis
            if( overlap_interval.Min() < lambda0 ) overlap_interval.Min() = lambda0;
            if( overlap_interval.Max() > lambda1 ) overlap_interval.Max() = lambda1;
            if( overlap_interval.IsEmpty() )
                return false;
        }
    }
    return true;
}

/*\todo MAYBE
template <unsigned D>
finline bool TestRay( const GAABB<D>& aabb,
        const mal::GVec<Real,D>& ray_pos, const mal::GVec<Real,D>& ray_dir, const mal::GInterval<Real>& interval, Real thickness )
*/

}} //namespace geo::bv

#endif // GEO_BV_GAABB_H
