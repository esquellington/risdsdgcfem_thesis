#include <Geo/Config.h>
#include "GKDOP.h"

namespace geo {
namespace bv {

/* p in Line = line_pos + \lambda in line_interval * line_dir;
   (x,y) in Slab- = Ax + By - slab_interval.min = 0
   (x,y) in Slab+ = Ax + By - slab_interval.max = 0

   Subst p into Slab-, solving for \lambda-
   \lambda- = (slab_interval.min - slab_dir*line_pos) / line_dir * slab_axis
*/
Interval ClipLineWithSlab( const Vec2& line_pos, const Vec2& line_dir, const Interval& line_interval,
                           const Vec2& slab_axis, const Interval& slab_interval )
{
    Real dot_axis_dir( mal::Dot( slab_axis, line_dir ) );
    Real dot_axis_pos( mal::Dot( slab_axis, line_pos ) );
    if( mal::Abs(dot_axis_dir) < 0.0001 )
    {
        if( slab_interval.TestOverlap( dot_axis_pos ) )
            return line_interval;
        else
            return Interval::Empty();
    }
    else
    {
        Real rcp_dot_axis_dir( mal::Rcp(dot_axis_dir) );
        Real lambda0 = rcp_dot_axis_dir*(slab_interval.Min() - dot_axis_pos);
        Real lambda1 = rcp_dot_axis_dir*(slab_interval.Max() - dot_axis_pos);
        if( lambda1 < lambda0 )
        {
            Real tmp( lambda0 );
            lambda0 = lambda1;
            lambda1 = tmp;
        }
        return Interval( mal::Max( line_interval.Min(), lambda0 ),
                         mal::Min( line_interval.Max(), lambda1 ) );
    }
}

/*
                 1
                 |      2
                 |     /
              /5---4\ /
             /       \
             6       3
             |       |--- 0
             7       2
             \       /
              \0---1/ \
                       \
                        3

Alg 2D O(n^2)
- Per cada parell de rectes d'un slab Si, clipejar-les
  contra TOTS els (altres) slabs.
Alg 3D O(n^3)
- Per cada parell d'slabs Si,Sj, trobar recta
  interseccio, i clipejar-la contra TOTS els (altres)
  slabs, fent early-out quan estigui completament fora
  d'algun d'ells.

  Com que totes les direccions son fixes, les uniques
  variables son els coeficients D dels plans/rectes, de
  manera que els denominadors de totes les divisions es
  poden guardar com a reciprocs constants.

  \todo General case, works for any K%2 == 0, generalize to any K, if
  DOP2_K>8 is required
*/
unsigned int DOP2_K8_Segments( const DOP2_K8& dop2k8, std::pair<Vec2,Vec2>* vec_segments )
{
    unsigned int num_segments(0);
    // For each slab axis d1
    for( int d1 = 0; d1<4; d1++ )
    {
        Vec2 n1( GKDOP_Axis<2,8>(d1) );
        Vec2 v1( mal::PerpendicularCW( n1 ) );
        // Consider both min and max segments perpendicular to slab axis d1
        Vec2 pA( dop2k8.GetInterval(d1).Min() * n1 );
        Vec2 pB( dop2k8.GetInterval(d1).Max() * n1 );
        Interval intervalA( Interval::Infinity() );
        Interval intervalB( Interval::Infinity() );
        // For each other slab axis d2, clip segments interval
        for( int d2 = 0; d2<4; d2++ )
        {
            if( d1 != d2 )
            {
                intervalA = ClipLineWithSlab( pA, v1, intervalA,
                                              GKDOP_Axis<2,8>(d2),
                                              dop2k8.GetInterval(d2) );
                intervalB = ClipLineWithSlab( pB, v1, intervalB,
                                              GKDOP_Axis<2,8>(d2),
                                              dop2k8.GetInterval(d2) );
            }
        }
        // Add clipped segments if not empty
        if( !intervalA.IsEmpty() )
        {
            vec_segments[num_segments].first = pA + intervalA.Min() * v1;
            vec_segments[num_segments].second = pA + intervalA.Max() * v1;
            num_segments++;
        }
        if( !intervalB.IsEmpty() )
        {
            vec_segments[num_segments].first = pB + intervalB.Min() * v1;
            vec_segments[num_segments].second = pB + intervalB.Max() * v1;
            num_segments++;
        }
    }
    return num_segments;
}

template <>
Real GKDOP_ComputeVolume_Exact( const GKDOP<2,4>& dop2k4 )
{
    return dop2k4.GetInterval(0).Length() * dop2k4.GetInterval(1).Length();
}

template <>
Real GKDOP_ComputeVolume_Exact( const GKDOP<2,8>& dop2k8 )
{
    // Get KDOP segments
    std::pair<Vec2,Vec2> vec_segments[8];
    unsigned int num_segments = DOP2_K8_Segments( dop2k8, vec_segments );
    Vec2 p0( vec_segments[0].first );
    Real twice_area(0);
    for( unsigned int it_segment=0; it_segment < num_segments; it_segment++ )
        twice_area += mal::Abs( mal::Dot(vec_segments[it_segment].first - p0,
                                         mal::PerpendicularCW(vec_segments[it_segment].second - p0) )); //fake cross-product norm
    return Real(0.5)*twice_area;
}

template <>
Real GKDOP_ComputeVolume_Exact( const GKDOP<3,6>& dop3k6 )
{
    return dop3k6.GetInterval(0).Length() * dop3k6.GetInterval(1).Length() * dop3k6.GetInterval(2).Length();
}

Real ComputeVolume( const DOP2_K4& dop ) { return GKDOP_ComputeVolume_Exact(dop); }
Real ComputeVolume( const DOP2_K8& dop ) { return GKDOP_ComputeVolume_Exact(dop); }
Real ComputeVolume( const DOP3_K6& dop ) { return GKDOP_ComputeVolume_Exact(dop); }
Real ComputeVolume( const DOP3_K14& dop ) { return GKDOP_ComputeVolume_Fast(dop); } //TEMP: Not exact

}} //namespace geo::bv
