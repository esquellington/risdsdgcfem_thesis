#ifndef GEO_BV_GBDOP_H
#define GEO_BV_GBDOP_H

#include <Mal/Config.h>
// #include <Geo/bv/BoundingVolume.h>
// #include "stats.h"
#include <Geo/bv/GAABB.h>

namespace geo {
namespace bv {

/* Barycentric Discrete Orientation Polytope:
   - Similar to an AABB but in barycentric coords, instead of cartesian ones.
   - MUST be as small as possible => Does NOT inherit from bv::BoundingVolume
   - ONLY offers a the minimal subset of BV ops required for DCR/BDOP-Trees
   IMPORTANT: BDOP is not a full BV
*/
template <unsigned B>
class GBDOP
{
public:
    inline GBDOP() {}
    const Interval& operator[](int i) const { return m_vecInterval[i]; }
    Interval& operator[](int i) { return m_vecInterval[i]; }
private:
    Interval m_vecInterval[B];
};

//! \name Common specializations
//@{
typedef GBDOP<2> BDOP2; //!< BDOP in a Segment
typedef GBDOP<3> BDOP3; //!< BDOP in a Triangle
typedef GBDOP<4> BDOP4; //!< BDOP in a Tetrahedron
//@}

/* Specialization of BV(BDOP) with tight AABB(BDOP) computation based
   on the papers:
   - [2007] Adaptive Deformations with Fast Tight Bounds
   - [2008] Tight and efficient surface bounds in meshless animation
   \todo Consider manual optimization
*/
inline AABB2 Compute_AABB2_From_BDOP3( const BDOP3& bdop, const Vec2* vec_pos )
{
    // sort node indices along X,Y
    int snie_x[3] = { 0, 1, 2 };
    if( vec_pos[snie_x[0]].x() > vec_pos[snie_x[1]].x() ) std::swap( snie_x[0], snie_x[1] );
    if( vec_pos[snie_x[1]].x() > vec_pos[snie_x[2]].x() ) std::swap( snie_x[1], snie_x[2] );
    if( vec_pos[snie_x[0]].x() > vec_pos[snie_x[1]].x() ) std::swap( snie_x[0], snie_x[1] );
    int snie_y[3] = { 0, 1, 2 };
    if( vec_pos[snie_y[0]].y() > vec_pos[snie_y[1]].y() ) std::swap( snie_y[0], snie_y[1] );
    if( vec_pos[snie_y[1]].y() > vec_pos[snie_y[2]].y() ) std::swap( snie_y[1], snie_y[2] );
    if( vec_pos[snie_y[0]].y() > vec_pos[snie_y[1]].y() ) std::swap( snie_y[0], snie_y[1] );
    // compute min/max
    return bv::AABB2().SetMinMax( Vec2( bdop[snie_x[0]].Max()*vec_pos[snie_x[0]].x() + bdop[snie_x[2]].Min()*vec_pos[snie_x[2]].x()
                                        + (1-bdop[snie_x[0]].Max()-bdop[snie_x[2]].Min()) * vec_pos[snie_x[1]].x(),
                                        bdop[snie_y[0]].Max()*vec_pos[snie_y[0]].y() + bdop[snie_y[2]].Min()*vec_pos[snie_y[2]].y()
                                        + (1-bdop[snie_y[0]].Max()-bdop[snie_y[2]].Min()) * vec_pos[snie_y[1]].y() ),
                                  Vec2( bdop[snie_x[2]].Max()*vec_pos[snie_x[2]].x() + bdop[snie_x[0]].Min()*vec_pos[snie_x[0]].x()
                                        + (1-bdop[snie_x[2]].Max()-bdop[snie_x[0]].Min()) * vec_pos[snie_x[1]].x(),
                                        bdop[snie_y[2]].Max()*vec_pos[snie_y[2]].y() + bdop[snie_y[0]].Min()*vec_pos[snie_y[0]].y()
                                        + (1-bdop[snie_y[2]].Max()-bdop[snie_y[0]].Min()) * vec_pos[snie_y[1]].y() ) );
}

/* Specialization of BV(BDOP) with tight AABB(BDOP) computation based
   on the papers:
   - [2007] Adaptive Deformations with Fast Tight Bounds
   - [2008] Tight and efficient surface bounds in meshless animation
   \todo Consider manual optimization
*/
inline AABB3 Compute_AABB3_From_BDOP4( const BDOP4& bdop, const Vec3* vec_pos )
{
    // sort node indices along X,Y,Z
    int snie_x[4] = { 0, 1, 2, 3 };
    int snie_y[4] = { 0, 1, 2, 3 };
    int snie_z[4] = { 0, 1, 2, 3 };
    //We use Sort4_Bubbble, 2x faster than generic std::sort() and simpler/shorter code than Sort4_HC while slightly (5%) slower
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( vec_pos[snie_x[j]].x() > vec_pos[snie_x[j+1]].x() ) std::swap(snie_x[j],snie_x[j+1]);
            if( vec_pos[snie_y[j]].y() > vec_pos[snie_y[j+1]].y() ) std::swap(snie_y[j],snie_y[j+1]);
            if( vec_pos[snie_z[j]].z() > vec_pos[snie_z[j+1]].z() ) std::swap(snie_z[j],snie_z[j+1]);
        }
    // compute min/max \sa Eq(8) in [2007] Adaptive Deformations with Fast Tight Bounds \todo Consider manual optimization
    return bv::AABB3().SetMinMax( Vec3( bdop[snie_x[0]].Max()*vec_pos[snie_x[0]].x() + bdop[snie_x[3]].Min()*vec_pos[snie_x[3]].x() + bdop[snie_x[2]].Min()*vec_pos[snie_x[2]].x()
                                        + (1-bdop[snie_x[0]].Max()-bdop[snie_x[3]].Min()-bdop[snie_x[2]].Min()) * vec_pos[snie_x[1]].x(),
                                        bdop[snie_y[0]].Max()*vec_pos[snie_y[0]].y() + bdop[snie_y[3]].Min()*vec_pos[snie_y[3]].y() + bdop[snie_y[2]].Min()*vec_pos[snie_y[2]].y()
                                        + (1-bdop[snie_y[0]].Max()-bdop[snie_y[3]].Min()-bdop[snie_y[2]].Min()) * vec_pos[snie_y[1]].y(),
                                        bdop[snie_z[0]].Max()*vec_pos[snie_z[0]].z() + bdop[snie_z[3]].Min()*vec_pos[snie_z[3]].z() + bdop[snie_z[2]].Min()*vec_pos[snie_z[2]].z()
                                        + (1-bdop[snie_z[0]].Max()-bdop[snie_z[3]].Min()-bdop[snie_z[2]].Min()) * vec_pos[snie_z[1]].z() ),
                                  Vec3( bdop[snie_x[3]].Max()*vec_pos[snie_x[3]].x() + bdop[snie_x[0]].Min()*vec_pos[snie_x[0]].x() + bdop[snie_x[1]].Min()*vec_pos[snie_x[1]].x()
                                        + (1-bdop[snie_x[3]].Max()-bdop[snie_x[0]].Min()-bdop[snie_x[1]].Min()) * vec_pos[snie_x[2]].x(),
                                        bdop[snie_y[3]].Max()*vec_pos[snie_y[3]].y() + bdop[snie_y[0]].Min()*vec_pos[snie_y[0]].y() + bdop[snie_y[1]].Min()*vec_pos[snie_y[1]].y()
                                        + (1-bdop[snie_y[3]].Max()-bdop[snie_y[0]].Min()-bdop[snie_y[1]].Min()) * vec_pos[snie_y[2]].y(),
                                        bdop[snie_z[3]].Max()*vec_pos[snie_z[3]].z() + bdop[snie_z[0]].Min()*vec_pos[snie_z[0]].z() + bdop[snie_z[1]].Min()*vec_pos[snie_z[1]].z()
                                        + (1-bdop[snie_z[3]].Max()-bdop[snie_z[0]].Min()-bdop[snie_z[1]].Min()) * vec_pos[snie_z[2]].z() ) );
}

/* Sort 3 node-in-element indices along X and Y
   \pre vec_nie_xy == { 0, 1, 2 ; 0, 1, 2 }
   \post vec_pos[ vec_nie_a[i] ].a() <= vec_pos[ vec_nie_a[i+1] ].a()

   \note For a generic GKDOP<D,K> the vec_nie_xy... array must
         allocate K/2 sub-arrays arrays of D node-in-element indices)
*/
inline void Compute_AABB2_From_BDOP3_PreSort( const Vec2* vec_pos, int* vec_nie_xy )
{
    // vec_nie_x[0] = 0; vec_nie_x[1] = 1; vec_nie_x[2] = 2;
    if( vec_pos[vec_nie_xy[0]].x() > vec_pos[vec_nie_xy[1]].x() ) std::swap( vec_nie_xy[0], vec_nie_xy[1] );
    if( vec_pos[vec_nie_xy[1]].x() > vec_pos[vec_nie_xy[2]].x() ) std::swap( vec_nie_xy[1], vec_nie_xy[2] );
    if( vec_pos[vec_nie_xy[0]].x() > vec_pos[vec_nie_xy[1]].x() ) std::swap( vec_nie_xy[0], vec_nie_xy[1] );
    // vec_nie_xy[0] = 0; vec_nie_xy[1] = 1; vec_nie_xy[2] = 2;
    const int y0(3); const int y1(4); const int y2(5);
    if( vec_pos[vec_nie_xy[y0]].y() > vec_pos[vec_nie_xy[y1]].y() ) std::swap( vec_nie_xy[y0], vec_nie_xy[y1] );
    if( vec_pos[vec_nie_xy[y1]].y() > vec_pos[vec_nie_xy[y2]].y() ) std::swap( vec_nie_xy[y1], vec_nie_xy[y2] );
    if( vec_pos[vec_nie_xy[y0]].y() > vec_pos[vec_nie_xy[y1]].y() ) std::swap( vec_nie_xy[y0], vec_nie_xy[y1] );
}

/* Compute AABB2 from BDOP and its pre-sorted nodes-in-element along X and Y */
inline AABB2 Compute_AABB2_From_BDOP3_PreSorted( const BDOP3& bdop, const Vec2* vec_pos, const int* vec_nie_xy )
{
    //const int y0(3); const int y1(4); const int y2(5);
    return bv::AABB2().SetMinMax( Vec2( bdop[vec_nie_xy[0]] .Max()*vec_pos[vec_nie_xy[0]] .x() + bdop[vec_nie_xy[2]] .Min()*vec_pos[vec_nie_xy[2]] .x()
                                        + (1-bdop[vec_nie_xy[0]] .Max()-bdop[vec_nie_xy[2]] .Min()) * vec_pos[vec_nie_xy[1]] .x(),
                                        bdop[vec_nie_xy[3]].Max()*vec_pos[vec_nie_xy[3]].y() + bdop[vec_nie_xy[5]].Min()*vec_pos[vec_nie_xy[5]].y()
                                        + (1-bdop[vec_nie_xy[3]].Max()-bdop[vec_nie_xy[5]].Min()) * vec_pos[vec_nie_xy[4]].y() ),
                                  Vec2( bdop[vec_nie_xy[2]] .Max()*vec_pos[vec_nie_xy[2]] .x() + bdop[vec_nie_xy[0]] .Min()*vec_pos[vec_nie_xy[0]] .x()
                                        + (1-bdop[vec_nie_xy[2]] .Max()-bdop[vec_nie_xy[0]] .Min()) * vec_pos[vec_nie_xy[1]] .x(),
                                        bdop[vec_nie_xy[5]].Max()*vec_pos[vec_nie_xy[5]].y() + bdop[vec_nie_xy[3]].Min()*vec_pos[vec_nie_xy[3]].y()
                                        + (1-bdop[vec_nie_xy[5]].Max()-bdop[vec_nie_xy[3]].Min()) * vec_pos[vec_nie_xy[4]].y() ) );
}

/* Sort 3 node-in-element indices along X,Y,Z
   \pre vec_nie_xyz == { 0, 1, 2, 3 ; 0, 1, 2, 3 ; 0, 1, 2, 3 }
   \post vec_pos[ vec_nie_a[i] ].a() <= vec_pos[ vec_nie_a[i+1] ].a()

   \note For a generic GKDOP<D,K> the vec_nie_xyz... array must
         allocate K/2 sub-arrays arrays of D node-in-element indices)
*/
inline void Compute_AABB3_From_BDOP4_PreSort( const Vec3* vec_pos, int* vec_nie_xyz )
{
    //We use Sort4_Bubbble, 2x faster than generic std::sort() and simpler/shorter code than Sort4_HC while slightly (5%) slower
    for( int i=0; i<3; i++ )
        for( int j=0; j<3; j++ )
        {
            if( vec_pos[vec_nie_xyz[j]].x() > vec_pos[vec_nie_xyz[j+1]].x() ) std::swap(vec_nie_xyz[j],vec_nie_xyz[j+1]);
            if( vec_pos[vec_nie_xyz[4+j]].y() > vec_pos[vec_nie_xyz[4+j+1]].y() ) std::swap(vec_nie_xyz[j+4],vec_nie_xyz[4+j+1]);
            if( vec_pos[vec_nie_xyz[8+j]].z() > vec_pos[vec_nie_xyz[8+j+1]].z() ) std::swap(vec_nie_xyz[j+8],vec_nie_xyz[8+j+1]);
        }
}

/* Compute AABB3 from BDOP and its pre-sorted nodes-in-element along X,Y,Z */
inline AABB3 Compute_AABB3_From_BDOP4_PreSorted( const BDOP4& bdop, const Vec3* vec_pos, const int* vec_nie_xyz )
{
    // const int y0(4); const int y1(5); const int y2(6); const int y3(7);
    // const int z0(8); const int z1(9); const int z2(10); const int z3(11);
    return bv::AABB3().SetMinMax( Vec3( bdop[vec_nie_xyz[0]].Max()*vec_pos[vec_nie_xyz[0]].x() + bdop[vec_nie_xyz[3]].Min()*vec_pos[vec_nie_xyz[3]].x() + bdop[vec_nie_xyz[2]].Min()*vec_pos[vec_nie_xyz[2]].x()
                                        + (1-bdop[vec_nie_xyz[0]].Max()-bdop[vec_nie_xyz[3]].Min()-bdop[vec_nie_xyz[2]].Min()) * vec_pos[vec_nie_xyz[1]].x(),
                                        bdop[vec_nie_xyz[4]].Max()*vec_pos[vec_nie_xyz[4]].y() + bdop[vec_nie_xyz[7]].Min()*vec_pos[vec_nie_xyz[7]].y() + bdop[vec_nie_xyz[6]].Min()*vec_pos[vec_nie_xyz[6]].y()
                                        + (1-bdop[vec_nie_xyz[4]].Max()-bdop[vec_nie_xyz[7]].Min()-bdop[vec_nie_xyz[6]].Min()) * vec_pos[vec_nie_xyz[5]].y(),
                                        bdop[vec_nie_xyz[8]].Max()*vec_pos[vec_nie_xyz[8]].z() + bdop[vec_nie_xyz[11]].Min()*vec_pos[vec_nie_xyz[11]].z() + bdop[vec_nie_xyz[10]].Min()*vec_pos[vec_nie_xyz[10]].z()
                                        + (1-bdop[vec_nie_xyz[8]].Max()-bdop[vec_nie_xyz[11]].Min()-bdop[vec_nie_xyz[10]].Min()) * vec_pos[vec_nie_xyz[9]].z() ),
                                  Vec3( bdop[vec_nie_xyz[3]].Max()*vec_pos[vec_nie_xyz[3]].x() + bdop[vec_nie_xyz[0]].Min()*vec_pos[vec_nie_xyz[0]].x() + bdop[vec_nie_xyz[1]].Min()*vec_pos[vec_nie_xyz[1]].x()
                                        + (1-bdop[vec_nie_xyz[3]].Max()-bdop[vec_nie_xyz[0]].Min()-bdop[vec_nie_xyz[1]].Min()) * vec_pos[vec_nie_xyz[2]].x(),
                                        bdop[vec_nie_xyz[7]].Max()*vec_pos[vec_nie_xyz[7]].y() + bdop[vec_nie_xyz[4]].Min()*vec_pos[vec_nie_xyz[4]].y() + bdop[vec_nie_xyz[5]].Min()*vec_pos[vec_nie_xyz[5]].y()
                                        + (1-bdop[vec_nie_xyz[7]].Max()-bdop[vec_nie_xyz[4]].Min()-bdop[vec_nie_xyz[5]].Min()) * vec_pos[vec_nie_xyz[6]].y(),
                                        bdop[vec_nie_xyz[11]].Max()*vec_pos[vec_nie_xyz[11]].z() + bdop[vec_nie_xyz[8]].Min()*vec_pos[vec_nie_xyz[8]].z() + bdop[vec_nie_xyz[9]].Min()*vec_pos[vec_nie_xyz[9]].z()
                                        + (1-bdop[vec_nie_xyz[11]].Max()-bdop[vec_nie_xyz[8]].Min()-bdop[vec_nie_xyz[9]].Min()) * vec_pos[vec_nie_xyz[10]].z() ) );
}

/* Test overlap between BDOP by casting them into an AABB3
   incrementally, which allows:
   - Per-axis early-out
   - Max2<Min1 and Max1<Min2 early-out

   The results are equivalent but around 5% faster than first casting into AABB and then testing AABB2AABB:
     const bv::AABB3 bv1_2( bv::Compute_AABB3_From_BDOP4_PreSorted( bdop1_1, vec_node_pos1_2, vec_sorted_nie1 ) );
     const bv::AABB3 bv2_2( bv::Compute_AABB3_From_BDOP4_PreSorted( bdop2_2, vec_node_pos2_2, vec_sorted_nie2 ) );
     return bv::TestOverlap( bv1_2, bv2_2 );
 */
inline bool TestOverlap_BDOP4_BDOP4_ComputeVolume( const BDOP4& bdop1, const Vec3* vec_pos1, const int* vec_nie_xyz1,
                                                   const BDOP4& bdop2, const Vec3* vec_pos2, const int* vec_nie_xyz2,
                                                   Real& volume1, Real& volume2 )
{
    GEO_BV_STAT_INC( testoverlap.m_GBDOP_GBDOP );
    Real acc_volume1(1), acc_volume2(1);
    for( int i=0; i<3; i++ )
    {
        const int i0(i*4); const int i1(i*4+1); const int i2(i*4+2); const int i3(i*4+3);
        Real interval1_min( bdop1[vec_nie_xyz1[i0]].Max()*vec_pos1[vec_nie_xyz1[i0]][i] + bdop1[vec_nie_xyz1[i3]].Min()*vec_pos1[vec_nie_xyz1[i3]][i] + bdop1[vec_nie_xyz1[i2]].Min()*vec_pos1[vec_nie_xyz1[i2]][i]
                            + (1-bdop1[vec_nie_xyz1[i0]].Max()-bdop1[vec_nie_xyz1[i3]].Min()-bdop1[vec_nie_xyz1[i2]].Min()) * vec_pos1[vec_nie_xyz1[i1]][i] );
        Real interval2_max( bdop2[vec_nie_xyz2[i3]].Max()*vec_pos2[vec_nie_xyz2[i3]][i] + bdop2[vec_nie_xyz2[i0]].Min()*vec_pos2[vec_nie_xyz2[i0]][i] + bdop2[vec_nie_xyz2[i1]].Min()*vec_pos2[vec_nie_xyz2[i1]][i]
                            + (1-bdop2[vec_nie_xyz2[i3]].Max()-bdop2[vec_nie_xyz2[i0]].Min()-bdop2[vec_nie_xyz2[i1]].Min()) * vec_pos2[vec_nie_xyz2[i2]][i] );
        if( interval2_max < interval1_min ) return false;
        Real interval1_max( bdop1[vec_nie_xyz1[i3]].Max()*vec_pos1[vec_nie_xyz1[i3]][i] + bdop1[vec_nie_xyz1[i0]].Min()*vec_pos1[vec_nie_xyz1[i0]][i] + bdop1[vec_nie_xyz1[i1]].Min()*vec_pos1[vec_nie_xyz1[i1]][i]
                            + (1-bdop1[vec_nie_xyz1[i3]].Max()-bdop1[vec_nie_xyz1[i0]].Min()-bdop1[vec_nie_xyz1[i1]].Min()) * vec_pos1[vec_nie_xyz1[i2]][i] );
        Real interval2_min( bdop2[vec_nie_xyz2[i0]].Max()*vec_pos2[vec_nie_xyz2[i0]][i] + bdop2[vec_nie_xyz2[i3]].Min()*vec_pos2[vec_nie_xyz2[i3]][i] + bdop2[vec_nie_xyz2[i2]].Min()*vec_pos2[vec_nie_xyz2[i2]][i]
                            + (1-bdop2[vec_nie_xyz2[i0]].Max()-bdop2[vec_nie_xyz2[i3]].Min()-bdop2[vec_nie_xyz2[i2]].Min()) * vec_pos2[vec_nie_xyz2[i1]][i] );
        if( interval1_max < interval2_min ) return false;
        acc_volume1 *= interval1_max - interval1_min;
        acc_volume2 *= interval2_max - interval2_min;
        /* Clearer code without per min/max early out, was SLIGHTLY slower
        Interval interval1( bdop1[vec_nie_xyz1[i0]].Max()*vec_pos1[vec_nie_xyz1[i0]][i] + bdop1[vec_nie_xyz1[i3]].Min()*vec_pos1[vec_nie_xyz1[i3]][i] + bdop1[vec_nie_xyz1[i2]].Min()*vec_pos1[vec_nie_xyz1[i2]][i]
                            + (1-bdop1[vec_nie_xyz1[i0]].Max()-bdop1[vec_nie_xyz1[i3]].Min()-bdop1[vec_nie_xyz1[i2]].Min()) * vec_pos1[vec_nie_xyz1[i1]][i],
                            bdop1[vec_nie_xyz1[i3]].Max()*vec_pos1[vec_nie_xyz1[i3]][i] + bdop1[vec_nie_xyz1[i0]].Min()*vec_pos1[vec_nie_xyz1[i0]][i] + bdop1[vec_nie_xyz1[i1]].Min()*vec_pos1[vec_nie_xyz1[i1]][i]
                            + (1-bdop1[vec_nie_xyz1[i3]].Max()-bdop1[vec_nie_xyz1[i0]].Min()-bdop1[vec_nie_xyz1[i1]].Min()) * vec_pos1[vec_nie_xyz1[i2]][i] );
        Interval interval2( bdop2[vec_nie_xyz2[i0]].Max()*vec_pos2[vec_nie_xyz2[i0]][i] + bdop2[vec_nie_xyz2[i3]].Min()*vec_pos2[vec_nie_xyz2[i3]][i] + bdop2[vec_nie_xyz2[i2]].Min()*vec_pos2[vec_nie_xyz2[i2]][i]
                            + (1-bdop2[vec_nie_xyz2[i0]].Max()-bdop2[vec_nie_xyz2[i3]].Min()-bdop2[vec_nie_xyz2[i2]].Min()) * vec_pos2[vec_nie_xyz2[i1]][i],
                            bdop2[vec_nie_xyz2[i3]].Max()*vec_pos2[vec_nie_xyz2[i3]][i] + bdop2[vec_nie_xyz2[i0]].Min()*vec_pos2[vec_nie_xyz2[i0]][i] + bdop2[vec_nie_xyz2[i1]].Min()*vec_pos2[vec_nie_xyz2[i1]][i]
                            + (1-bdop2[vec_nie_xyz2[i3]].Max()-bdop2[vec_nie_xyz2[i0]].Min()-bdop2[vec_nie_xyz2[i1]].Min()) * vec_pos2[vec_nie_xyz2[i2]][i] );
        if( !interval1.TestOverlap(interval2) ) return false;
        acc_volume1 *= interval1.Length();
        acc_volume2 *= interval2.Length();
        */
    }
    volume1 = acc_volume1;
    volume2 = acc_volume2;
    return true;
}

/* Test overlap along an arbitrary axis, given tetrahedron node
   positions, and their projections and order along such axis.
*/
inline bool TestOverlap_BDOP4_BDOP4_Axis( const BDOP4& bdop1, const Real* vec_proj1, const int* vec_nie1,
                                          const BDOP4& bdop2, const Real* vec_proj2, const int* vec_nie2 )
                                          // Real& volume1, Real& volume2 )
{
    GEO_BV_STAT_INC( testoverlap.m_GBDOP_GBDOP_Axis );
    Real interval1_min( bdop1[vec_nie1[0]].Max()*vec_proj1[vec_nie1[0]] + bdop1[vec_nie1[3]].Min()*vec_proj1[vec_nie1[3]] + bdop1[vec_nie1[2]].Min()*vec_proj1[vec_nie1[2]]
                        + (1-bdop1[vec_nie1[0]].Max()-bdop1[vec_nie1[3]].Min()-bdop1[vec_nie1[2]].Min()) * vec_proj1[vec_nie1[1]] );
    Real interval2_max( bdop2[vec_nie2[3]].Max()*vec_proj2[vec_nie2[3]] + bdop2[vec_nie2[0]].Min()*vec_proj2[vec_nie2[0]] + bdop2[vec_nie2[1]].Min()*vec_proj2[vec_nie2[1]]
                        + (1-bdop2[vec_nie2[3]].Max()-bdop2[vec_nie2[0]].Min()-bdop2[vec_nie2[1]].Min()) * vec_proj2[vec_nie2[2]] );
    if( interval2_max < interval1_min ) return false;
    Real interval1_max( bdop1[vec_nie1[3]].Max()*vec_proj1[vec_nie1[3]] + bdop1[vec_nie1[0]].Min()*vec_proj1[vec_nie1[0]] + bdop1[vec_nie1[1]].Min()*vec_proj1[vec_nie1[1]]
                        + (1-bdop1[vec_nie1[3]].Max()-bdop1[vec_nie1[0]].Min()-bdop1[vec_nie1[1]].Min()) * vec_proj1[vec_nie1[2]] );
    Real interval2_min( bdop2[vec_nie2[0]].Max()*vec_proj2[vec_nie2[0]] + bdop2[vec_nie2[3]].Min()*vec_proj2[vec_nie2[3]] + bdop2[vec_nie2[2]].Min()*vec_proj2[vec_nie2[2]]
                        + (1-bdop2[vec_nie2[0]].Max()-bdop2[vec_nie2[3]].Min()-bdop2[vec_nie2[2]].Min()) * vec_proj2[vec_nie2[1]] );
    if( interval1_max < interval2_min ) return false;
    // acc_volume1 *= interval1_max - interval1_min;
    // acc_volume2 *= interval2_max - interval2_min;
    return true;
}

}} //namespace geo::bv

#endif // GEO_BV_GBDOP_H
