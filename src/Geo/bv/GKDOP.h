#ifndef GEO_BV_GKDOP_H
#define GEO_BV_GKDOP_H

#include <Geo/bv/BoundingVolume.h>
#include <Geo/bv/GKDOP_Axis.h>
#include "stats.h"

namespace geo {
namespace bv {

template <unsigned D, unsigned K>
class GKDOP: public GBoundingVolumeD<D>
{
public:
    GEO_STATIC_ASSERT( K % 2 == 0 );
    enum EConstants { cDimension = D, cNumDirections = K/2 };
    typedef mal::GVec<Real,D> vec_type;
    Interval m_vecInterval[cNumDirections];

public:
    finline GKDOP()
    : GBoundingVolumeD<D>( (D==2)
                           ? ( (K==4) ? eBV_DOP2_K4 : eBV_DOP2_K8 )
                           : ( (K==6) ? eBV_DOP3_K6 : eBV_DOP3_K14 ) ) { /*Invervals are init to empty*/ }
    finline explicit GKDOP( const vec_type& pos ) : GBoundingVolumeD<D>( (D==2)
                                                                         ? ( (K==4) ? eBV_DOP2_K4 : eBV_DOP2_K8 )
                                                                         : ( (K==6) ? eBV_DOP3_K6 : eBV_DOP3_K14 ) ) { Set(pos); }
    finline ~GKDOP() {}

    finline void SetInterval( int index_d, const Interval& interval ) { m_vecInterval[index_d] = interval; IBoundingVolume::Touch(); }
    finline const Interval& GetInterval( int index_d ) const { return m_vecInterval[index_d]; }
    finline const vec_type GetPos() const { vec_type pos; for( int i=0; i<cDimension; i++ ) pos[i] = m_vecInterval[i].Mid(); return pos; }

    template <unsigned I> finline void SetInterval( const Interval& interval ) { m_vecInterval[I] = interval; IBoundingVolume::Touch(); }
    template <unsigned I> finline const Interval& GetInterval() const { return m_vecInterval[I]; }

    finline GKDOP& Set( const vec_type& point ) { for( unsigned int i=0; i<cNumDirections; i++ )
                                                      m_vecInterval[i].Set( GKDOP_Dot<D,K>(point,i) );
                                                  IBoundingVolume::Touch();
                                                  return *this; }

    finline GKDOP& Extend( Real thickness ) { for( unsigned int i=0; i<cNumDirections; i++ )
                                                  m_vecInterval[i].Extend( thickness );
                                              IBoundingVolume::Touch();
                                              return *this; }
    finline GKDOP& Merge( const GKDOP& other ) { for( unsigned int i=0; i<cNumDirections; i++ )
                                                     m_vecInterval[i].Merge( other.m_vecInterval[i] );
                                                 IBoundingVolume::Touch();
                                                 return *this; }
    finline GKDOP& Merge( const vec_type& point ) { for( unsigned int i=0; i<cNumDirections; i++ )
                                                        m_vecInterval[i].Merge( GKDOP_Dot<D,K>(point,i) );
                                                    IBoundingVolume::Touch();
                                                    return *this; }
    finline GKDOP& Merge( const vec_type& point, Real radius ) { GEO_ASSERT(false); return *this; }

};

//! \name Usual instantiations
//@{
//2D
typedef GKDOP<2,4> DOP2_K4;
typedef GKDOP<2,8> DOP2_K8;
//3D: Symmetric
typedef GKDOP<3,6> DOP3_K6;   //X,Y,Z (AABB)
typedef GKDOP<3,14> DOP3_K14; //X,Y,Z,X+Y+Z,X+Y-Z,X-Y+Z,X-Y-Z (AABB + cube vertex diagonals)
typedef GKDOP<3,18> DOP3_K18; //X,Y,Z,X+Y,X-Y,X+Z,X-Z,Y+Z,Y-Z (AABB + cube edge diagonals)
typedef GKDOP<3,24> DOP3_K24; //All: (AABB + cube vertex and edge diagonals)
//3D: Asymmetric
typedef GKDOP<3,2*(3+2)> DOP3_K10_XZ; //Axis: X,Y,Z,X+Z,X-Z (useful for 2.5D worlds with Y-up)

template <> struct bv_type_of<DOP2_K4> { enum { value = eBV_DOP2_K4 }; };
template <> struct bv_type_of<DOP2_K8> { enum { value = eBV_DOP2_K8 }; };
template <> struct bv_type_of<DOP3_K6> { enum { value = eBV_DOP3_K6 }; };
template <> struct bv_type_of<DOP3_K14> { enum { value = eBV_DOP3_K14 }; };
//@}

/*TODO Consider generic method(s) here... but that handle 2D and 3D differently
  template<unsigned K>
  ComputeSegments( const GKDOP<2,K>& dop2, std::pair<Vec2,Vec2>* vec_segments );
*/
unsigned int DOP2_K8_Segments( const DOP2_K8& dop2k8, std::pair<Vec2,Vec2>* vec_segments );

//Maximum distance between 2 points on the BV
//\todo By now we use GAABB<D> length, which is conservatively larger than the actual extent of an arbitrary KDOP
template <unsigned D, unsigned K>
finline Real ComputeLargestExtent( const GKDOP<D,K>& kdop )
{
    Real diag_sq(0);
    for( unsigned int i=0; i<D; i++ ) diag_sq += mal::Sq(kdop.m_vecInterval[i].Length());
    return mal::Sqrt(diag_sq);
}

//\todo ONLY implemented for DOP2_K4, DOP2_K8 and DOP3_K6 by now
template <unsigned D, unsigned K> Real GKDOP_ComputeVolume_Exact( const GKDOP<D,K>& kdop );

template <unsigned D, unsigned K>
finline Real GKDOP_ComputeVolume_Fast( const GKDOP<D,K>& kdop )
{
    /*The generic implementation cannot compute the actual volume,
      but a pseudo-volume by multiplying all K/2 interval lengths:
      - This is exact for AABB-like KDOP (K=4 int 2D, K=6 in 3D), but is
        clearly *not* a volume for other K. However, it provides a way
        to compare the relative "extension" of different KDOP with the
        same <D,K>.
      - Maybe we should try to prove some "triangular inequality" on
        this... or consider using sum of squared lengths instead of
        multiplication... we'll see.
    */
    Real pseudo_volume(1);
    for( unsigned int i=0; i<K/2; i++ ) pseudo_volume *= kdop.GetInterval(i).Length();
    return pseudo_volume;
}

finline Real ComputeVolume_Fast( const DOP2_K4& dop ) { return GKDOP_ComputeVolume_Fast(dop); }
finline Real ComputeVolume_Fast( const DOP2_K8& dop ) { return GKDOP_ComputeVolume_Fast(dop); }
finline Real ComputeVolume_Fast( const DOP3_K6& dop ) { return GKDOP_ComputeVolume_Fast(dop); }
finline Real ComputeVolume_Fast( const DOP3_K14& dop ) { return GKDOP_ComputeVolume_Fast(dop); }

Real ComputeVolume( const DOP2_K4& dop );
Real ComputeVolume( const DOP2_K8& dop );
Real ComputeVolume( const DOP3_K6& dop );
Real ComputeVolume( const DOP3_K14& dop );

template <unsigned D, unsigned K>
finline bool TestOverlap( const GKDOP<D,K>& kdop1, const GKDOP<D,K>& kdop2 )
{
#ifdef __GEO_ENABLE_STATS
    switch( kdop1.GetType() )
    {
    case eBV_DOP2_K4: GEO_BV_STAT_INC( testoverlap.m_DOP2K4_DOP2K4 ); break;
    case eBV_DOP2_K8: GEO_BV_STAT_INC( testoverlap.m_DOP2K8_DOP2K8 ); break;
    case eBV_DOP3_K6: GEO_BV_STAT_INC( testoverlap.m_DOP3K6_DOP3K6 ); break;
    case eBV_DOP3_K14: GEO_BV_STAT_INC( testoverlap.m_DOP3K14_DOP3K14 ); break;
    default: break;
    }
#endif
    for( unsigned int i=0; i<K/2; i++ )
        if( !kdop1.m_vecInterval[i].TestOverlap( kdop2.m_vecInterval[i] ) )
            return false;
    return true;
}

template <unsigned D, unsigned K>
finline bool TestRay( const GKDOP<D,K>& kdop,
                      const mal::GVec<Real,D>& ray_pos, const mal::GVec<Real,D>& ray_dir, const mal::GInterval<Real>& interval )
{
#ifdef __GEO_ENABLE_STATS
    // if( D == 2 ) { GEO_BV_STAT_INC( testray.m_DOP2 ); }
    // else { GEO_BV_STAT_INC( testray.m_AABB3 ); }
#endif
    mal::GInterval<Real> overlap_interval( interval );
    for( int it_axis=0; it_axis<K/2; it_axis++ )
    {
        // If parallel to slab, either overlaps for any lambda or for none.
        Real dot_kdop_raydir( GKDOP_Dot<D,K>(ray_dir,it_axis) );
        Real dot_kdop_raypos( GKDOP_Dot<D,K>(ray_pos,it_axis) );
        if( mal::Abs(dot_kdop_raydir) < 1e-4 ) //\todo g_pDefaultContext->m_Epsilon_Dir )
        {
            if( dot_kdop_raypos < kdop.m_vecInterval[it_axis].Min() || dot_kdop_raypos > kdop.m_vecInterval[it_axis].Max() )
                return false;
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            // \todo THIS CODE IS LONG AND UGLY, try to simplify it...
            Real inv_divisor( mal::Rcp(dot_kdop_raydir) );
            Real lambda0( ( kdop.m_vecInterval[it_axis].Min() - dot_kdop_raypos ) * inv_divisor );
            Real lambda1( ( kdop.m_vecInterval[it_axis].Max() - dot_kdop_raypos ) * inv_divisor );
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

}} //namespace geo::bv

#endif // GEO_BV_GKDOP_H
