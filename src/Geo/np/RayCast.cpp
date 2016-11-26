#include "RayCast.h"
#include <Geo/shape/TriSurfaceShape3.h>
#include <Geo/shape/MeshSolidShape2.h>
#include <Geo/shape/TetSolidShape3.h>
#include <Geo/np/Overlap.h>

namespace geo {
namespace np {

/* Coherent double-sided RC guarantees that the SAME rh will be
   generated for any permutation of triangle vertices. This is
   REQUIRED by Clip_Triangle3_Tetrahedron3() when used inside
   Clip_TriSurfaceShape3_TetSolidShape3() in order to produce the SAME
   intersection points for a given Tri against all TetSolid
   tetrahedrons, most of which share triangular faces with their
   neighbours.

   IMPORTANT: This is MUCH SLOWER (1.5x) than non-coherent version,
   use with care.
*/
bool RayCast_Triangle3_DoubleSided_COHERENT( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                             const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2,
                                             RayHit3& rh, RayCache3* p_rc,
                                             const Context* p_context )
{
    // GEO_NP_STAT_INC( raycast.m_Triangle3_DS_C );
#define __ENABLE_COHERENT
#ifdef __ENABLE_COHERENT
    std::vector<Vec3> vec_tri_pos = { triangle_p0, triangle_p1, triangle_p2 };
    std::sort( vec_tri_pos.begin(), vec_tri_pos.end(),
               []( const Vec3& a, const Vec3& b )
               {
                   bool less_x( a.x() < b.x() );
                   bool less_y( a.y() < b.y() );
                   bool less_z( a.z() < b.z() );
                   return less_x || (!less_x && less_y) || (!less_x && !less_y && less_z );
               } );
    return RayCast_Triangle3_DoubleSided( ray_pos, ray_dir, interval, vec_tri_pos[0], vec_tri_pos[1], vec_tri_pos[2], rh, p_rc );
#else
    return RayCast_Triangle3_DoubleSided( ray_pos, ray_dir, interval, triangle_p0, triangle_p1, triangle_p2, rh, p_rc );
#endif
}

/*! \todo Accept test flags to return first/closest-hit and single/double-sided */
bool RayCast_TriSurfaceShape3( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                               const TriSurfaceShape3* p_surface, const Transform3& surface_tr, const Real* p_surface_dof,
                               bool b_double_sided,
                               RayHit3& rh, RayCache3* p_rc,
                               const Context* p_context )
{
    GEO_NP_STAT_INC( raycast.m_TriSurface3 );
    const Vec3* vec_sdof( reinterpret_cast<const Vec3*>(p_surface_dof) );
    Vec3 ray_pos_local( surface_tr.Inverse()*ray_pos );
    Vec3 ray_dir_local( mal::Transposed(surface_tr.Rot())*ray_dir );
    rh.m_Interval = Interval::Empty(); //=>[inf,-inf]
    unsigned int num_hits(0);
#ifdef __GEO_TRISS_ENABLE_BVH__DISABLED_WHILE_DEVELOPING
    const BVH_TriSurfaceShape3* pBVH = p_surface->GetBVH(); //\todo ASSUME UP-TO-DATE!
    if( pBVH )
    {
        //\tod
    }
    else
#endif
    {
        RayHit3 tmp_rh;
        if( b_double_sided )
            for( unsigned int it_tri=0; it_tri < p_surface->GetNumT(); it_tri++ )
            {
                Vec3f p0( p_surface->V_Pos( p_surface->T_VID(it_tri,0), vec_sdof ) );
                Vec3f p1( p_surface->V_Pos( p_surface->T_VID(it_tri,1), vec_sdof ) );
                Vec3f p2( p_surface->V_Pos( p_surface->T_VID(it_tri,2), vec_sdof ) );
                if( RayCast_Triangle3_DoubleSided( ray_pos_local, ray_dir_local, interval, p0, p1, p2, tmp_rh ) )
                {
                    if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() ) rh = tmp_rh;
                    num_hits++;
                }
            }
        else
            for( unsigned int it_tri=0; it_tri < p_surface->GetNumT(); it_tri++ )
            {
                Vec3f p0( p_surface->V_Pos( p_surface->T_VID(it_tri,0), vec_sdof ) );
                Vec3f p1( p_surface->V_Pos( p_surface->T_VID(it_tri,1), vec_sdof ) );
                Vec3f p2( p_surface->V_Pos( p_surface->T_VID(it_tri,2), vec_sdof ) );
                if( RayCast_Triangle3_SingleSided( ray_pos_local, ray_dir_local, interval, p0, p1, p2, tmp_rh ) )
                {
                    if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() ) rh = tmp_rh;
                    num_hits++;
                }
            }
    }
    return num_hits > 0;
}


bool RayCast_DCR_E_DoubleSided( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                const DCR_MeshSolidShape2* p_dcr, uint32 eid,
                                const Mat3x3& Bs, const Mat3x3& invBs, const Mat3x3& Bs_invBm,
                                RayHit2& rh, RayCache2 *p_rc,
                                const Context* p_context )
{
    switch( p_context->m_RCDCR_Method )
    {
    case Context::eRCDCRM_BruteForce: return RayCast_DCR_E_DoubleSided_BruteForce( ray_pos, ray_dir, interval, p_dcr, eid, Bs, invBs, Bs_invBm, rh, p_rc, p_context ); break;
    case Context::eRCDCRM_BDT: return RayCast_DCR_E_DoubleSided_BDT( ray_pos, ray_dir, interval, p_dcr, eid, Bs, invBs, Bs_invBm, rh, p_rc, p_context ); break;
    default: return false; break;
    }
}

/* Performn RayCast on a DCR.E by testing all DCR.P.S inside
   - Only closest hit is reported
   - Computed in barycentric coords to avoid transforming DCR.E.V, but RayHit results in global coords as usual
*/
bool RayCast_DCR_E_DoubleSided_BruteForce( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                           const DCR_MeshSolidShape2* p_dcr, uint32 eid,
                                           const Mat3x3& Bs, const Mat3x3& invBs, const Mat3x3& Bs_invBm,
                                           RayHit2& rh, RayCache2 *p_rc,
                                           const Context* p_context )
{
    RayHit2 tmp_rh;
#define __USE_RC_DCR_E_DS_BF_Barycentric
#ifdef __USE_RC_DCR_E_DS_BF_Barycentric
    // Transform ray to barycentric coords (\todo Only b[1,2] are used, b[0] is implicit, could reduce to 2x2 MxV
    Vec3 b0( invBs * mal::Concat(1,ray_pos+interval.Min()*ray_dir) );
    Vec3 b1( invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
    Vec3 dir_b( b1-b0 );
    Interval interval_b(0,1);
    Mat3x3 invBm( invBs * Bs_invBm ); //TEMP: should be stored in DCR.E or passed as a param
    rh.m_Interval.Set(1);
#else
    rh.m_Interval.Set(interval.Max());
#endif

    const DCR_MeshSolidShape2::Element& dcr_e( p_dcr->m_vecE[eid] );
    for( unsigned int it_pie=0; it_pie<dcr_e.m_NumPatches; it_pie++ )
    {
        const DCR_MeshSolidShape2::Patch& dcr_p( p_dcr->m_vecP[dcr_e.m_FirstPID+it_pie] );
        for( unsigned int it_sip=0; it_sip<dcr_p.m_NumSegments; it_sip++ )
        {
            uint32 sid( dcr_p.m_FirstSID + it_sip );
#ifdef __USE_RC_DCR_E_DS_BF_Barycentric
            // Compute barycentric segment
            /* OPTIMIZATION \todo
               - This COULD be procomputed at the expense of A LOT OF extra memory in DCR.V... even if only b[1,2] are stored.
               - Also, a 2x2 MxV could be used because we NEVER use b[2]
            */
            Vec3 s0_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(0) ]) );
            Vec3 s1_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(1) ]) );
            // Perform raycast in barycentric coords using two independent coords b[1],b[2] (b[0]=1-b[1]-b[2])
            if( RayCast_Segment2_DoubleSided( mal::GRange<1,2>(b0), mal::GRange<1,2>(dir_b),
                                              Interval(0,rh.m_Interval.Min()),
                                              mal::GRange<1,2>(s0_b), mal::GRange<1,2>(s1_b),
                                              tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                rh = tmp_rh;
#else
            // IN GLOBAL COORDS
            Vec2 s0_0( mal::GRange<1,2>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(0) ]) ) );
            Vec2 s1_0( mal::GRange<1,2>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(1) ]) ) );
            //s0<->s1 swapped because p is assumed INSIDE p_dcr \todo This MAY NOT BE TRUE as p1 is an IB1.EP.m_AvgPos that MAY BE OUTSIDE DCR2
            if( RayCast_Segment2_DoubleSided( ray_pos, ray_dir,
                                              Interval(interval.Min(),rh.m_Interval.Min()), //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                                              s0_0, s1_0,
                                              tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                rh = tmp_rh;
#endif
        }
    }

#ifdef __USE_RC_DCR_E_DS_BF_Barycentric
    // Hit if any ray has actually hit a segment
    if( rh.m_Interval.Min() < interval_b.Max() )
    {
        // remap interval [0,1] to [min,max]
        rh.m_Interval.Min() = interval.Min() + rh.m_Interval.Min()*interval.Length();
        rh.m_Interval.Max() = interval.Min() + rh.m_Interval.Max()*interval.Length();
        // compute global and barycentric stuff
        rh.m_FeatureId = feature_id( eFT_Triangle, eid );
        Vec3 b( mal::Concat(1-mal::Sum(rh.m_Point),rh.m_Point) ); //b[0] = 1-b[1]-b[2]
        rh.m_Extra_BarycentricCoords = mal::Concat(b,0);
        rh.m_Point = mal::GRange<1,2>( Bs * b );
        //\todo rh.m_Normal = transform vector to global?!?!
        return true;
    }
    else
        return false;
#else
    // Hit if any ray has actually hit a segment
    if( rh.m_Interval.Min() < interval.Max() )
    {
        // compute global and barycentric stuff
        rh.m_FeatureId = feature_id( eFT_Triangle, eid );
        Vec3 b( invBs * mal::Concat(1,rh.m_Point) );
        rh.m_Extra_BarycentricCoords = mal::Concat(b,0);
        return true;
    }
    else
        return false;
#endif
}

bool Clip_BarycentricRay2_BDOP( const Vec3& ray_pos_b, const Vec3& ray_dir_b, const Interval& interval_b,
                                const bv::BDOP3& bdop,
                                Interval& clipped_interval_b,
                                const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( overlap.m_Segment2_Triangle2 ); //TEMP: To compare overlap/clip approaches
    clipped_interval_b = interval_b;
    // int first_hit_axis(0);
    for( int i=0; i<3; i++ )
    {
        if( mal::Abs(ray_dir_b[i]) < p_context->m_Epsilon_Dir )
        {
            if( ray_pos_b[i] < bdop[i].Min() || ray_pos_b[i] > bdop[i].Max() )
            {
                //Empty interval and return false
                clipped_interval_b = Interval::Empty();
                return false;
            }
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            Real inv_dir_b_i( mal::Rcp(ray_dir_b[i]) );
            Real lambda0( (bdop[i].Min()-ray_pos_b[i]) * inv_dir_b_i );
            Real lambda1( (bdop[i].Max()-ray_pos_b[i]) * inv_dir_b_i );
            if( lambda1 < lambda0 )
            {
                Real tmp( lambda0 );
                lambda0 = lambda1;
                lambda1 = tmp;
            }
            // Clip lambda-interval and update first-axis
            if( clipped_interval_b.Min() < lambda0 ) clipped_interval_b.Min() = lambda0;
            if( clipped_interval_b.Max() > lambda1 ) clipped_interval_b.Max() = lambda1;
            if( clipped_interval_b.IsEmpty() ) return false;
        }
    }
    return true;
}

/* Performn RayCast on a DCR.E using per-DCR.P BDOP-Tree
   - Only closest hit is reported
   - Computed in barycentric coords to avoid transforming DCR.E.V, but RayHit results in global coords as usual
   \todo Subst std::vector stack with std::array
*/
bool RayCast_DCR_E_DoubleSided_BDT( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                    const DCR_MeshSolidShape2* p_dcr, uint32 eid,
                                    const Mat3x3& Bs, const Mat3x3& invBs, const Mat3x3& Bs_invBm,
                                    RayHit2& rh, RayCache2 *p_rc,
                                    const Context* p_context )
{
    RayHit2 tmp_rh;
    // Transform ray to barycentric coords (\todo Only b[1,2] are used, b[0] is implicit, could reduce to 2x2 MxV
    Vec3 b0( invBs * mal::Concat(1,ray_pos+interval.Min()*ray_dir) );
    Vec3 b1( invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
    Vec3 dir_b( b1-b0 );
    Interval interval_b(0,1);
    Mat3x3 invBm( invBs * Bs_invBm ); //TEMP: should be stored in DCR.E or passed as a param
    rh.m_Interval.Set(1); //barycentric ray is a segment with lambda [0,1] that maps the euclidean segment interval [min,max]

    const DCR_MeshSolidShape2::Element& dcr_e( p_dcr->m_vecE[eid] );
    for( unsigned int it_pie=0; it_pie<dcr_e.m_NumPatches; it_pie++ )
    {
        const DCR_MeshSolidShape2::Patch& dcr_p( p_dcr->m_vecP[dcr_e.m_FirstPID+it_pie] );
        // Init stack
        typedef std::pair< DCR_MeshSolidShape2::Patch::bdt_node_sip_range,Interval> stack_entry_type;
        std::vector< stack_entry_type > stackBDTN;
        stackBDTN.push_back( stack_entry_type(DCR_MeshSolidShape2::Patch::bdt_node_sip_range(0,dcr_p.m_NumSegments),interval_b) );
        // GEO_LOG("DCR.P[%u].BDT",it_p);
        while( !stackBDTN.empty() )
        {
            // Pop PS subarray
            stack_entry_type se( stackBDTN.back() );
            stackBDTN.pop_back();
            const DCR_MeshSolidShape2::Patch::bdt_node_sip_range node_sr( se.first );
            const uint32 length_sr( node_sr.second - node_sr.first );
            if( length_sr <= p_context->m_RCDCR_BDT_MaxLeafSize )
            {
                // Test all segments in subtree
                for( unsigned int it_sin=0; it_sin<length_sr; it_sin++ )
                {
                    uint32 sid( dcr_p.m_FirstSID + node_sr.first + it_sin );
                    /* Compute barycentric segment
                       \todo OPTIMIZATION
                       - This COULD be procomputed at the expense of A LOT OF extra memory in DCR.V... even if only b[1,2] are stored.
                       - Also, a 2x2 MxV could be used because we NEVER use b[2]
                    */
                    Vec3 s0_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(0) ]) );
                    Vec3 s1_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecS[sid].GetVID(1) ]) );
                    // Perform raycast in barycentric coords using two independent coords b[1],b[2] (b[0]=1-b[1]-b[2])
                    if( RayCast_Segment2_DoubleSided( mal::GRange<1,2>(b0), mal::GRange<1,2>(dir_b),
                                                      Interval(0,rh.m_Interval.Min()),
                                                      mal::GRange<1,2>(s0_b), mal::GRange<1,2>(s1_b),
                                                      tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                        rh = tmp_rh;
                }
            }
            else
            {
                // Get node/segment data
                const DCR_MeshSolidShape2::Segment& bdtn( p_dcr->m_vecS[dcr_p.m_FirstSID+node_sr.first] );
                const bv::BDOP3 bdop( bdtn.BDTN_BDOPq() );
                // Cut recursive ray interval to current best hit
                if( se.second.Max() > rh.m_Interval.Min() ) se.second.Max() = rh.m_Interval.Min();
                // Test BDOP and clip interval
                Interval clipped_interval_b;
                if( Clip_BarycentricRay2_BDOP( b0, dir_b, se.second,
                                               bdop,
                                               clipped_interval_b ) ) //\todo IDEALLY, do this directly in BDOP3q space
                {
                    // Test node/segment
                    /* Compute barycentric segment
                       \todo OPTIMIZATION
                       - This COULD be procomputed at the expense of A LOT OF extra memory in DCR.V... even if only b[1,2] are stored.
                       - Also, a 2x2 MxV could be used because we NEVER use b[2]
                    */
                    Vec3 s0_b( invBm * mal::Concat(1,p_dcr->m_vecV[ bdtn.GetVID(0) ]) );
                    Vec3 s1_b( invBm * mal::Concat(1,p_dcr->m_vecV[ bdtn.GetVID(1) ]) );
                    // Perform raycast in barycentric coords using two independent coords b[1],b[2] (b[0]=1-b[1]-b[2])
                    if( RayCast_Segment2_DoubleSided( mal::GRange<1,2>(b0), mal::GRange<1,2>(dir_b),
                                                      Interval(0,rh.m_Interval.Min()),
                                                      mal::GRange<1,2>(s0_b), mal::GRange<1,2>(s1_b),
                                                      tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                        rh = tmp_rh;
                    // Recurse node using BDOP-clipped barycentric ray interval, so that all BCoords span shorten due to clipping
                    int remaining_sr( length_sr - 1 );
                    if( remaining_sr > 1 )
                    {
                        stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1, node_sr.first+1+remaining_sr/2 ), clipped_interval_b ) );
                        stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1+remaining_sr/2, node_sr.second ), clipped_interval_b ) );
                    }
                    else if( remaining_sr == 1 )
                        stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1, node_sr.first+2 ), clipped_interval_b ) );
                }
            }
        }
    }

    // Hit if any ray has actually hit a segment
    if( rh.m_Interval.Min() < interval_b.Max() )
    {
        // remap interval [0,1] to [min,max]
        rh.m_Interval.Min() = interval.Min() + rh.m_Interval.Min()*interval.Length();
        rh.m_Interval.Max() = interval.Min() + rh.m_Interval.Max()*interval.Length();
        // compute global and barycentric stuff
        rh.m_FeatureId = feature_id( eFT_Triangle, eid );
        Vec3 b( mal::Concat(1-mal::Sum(rh.m_Point),rh.m_Point) ); //b[0] = 1-b[1]-b[2]
        rh.m_Extra_BarycentricCoords = mal::Concat(b,0);
        rh.m_Point = mal::GRange<1,2>( Bs * b );
        //\todo rh.m_Normal = transform vector to global?!?!
        return true;
    }
    else
        return false;
}

bool RayCast_DCR_E_DoubleSided( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                const DCR_TetSolidShape3* p_dcr, uint32 eid,
                                const Mat4x4& Bs, const Mat4x4& invBs, const Mat4x4& Bs_invBm,
                                RayHit3& rh, RayCache3 *p_rc,
                                const Context* p_context )
{
    switch( p_context->m_RCDCR_Method )
    {
    case Context::eRCDCRM_BruteForce: return RayCast_DCR_E_DoubleSided_BruteForce( ray_pos, ray_dir, interval, p_dcr, eid, Bs, invBs, Bs_invBm, rh, p_rc, p_context ); break;
    case Context::eRCDCRM_BDT: return RayCast_DCR_E_DoubleSided_BDT( ray_pos, ray_dir, interval, p_dcr, eid, Bs, invBs, Bs_invBm, rh, p_rc, p_context ); break;
    default: return false; break;
    }
}

bool RayCast_DCR_E_DoubleSided_BruteForce( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                           const DCR_TetSolidShape3* p_dcr, uint32 eid,
                                           const Mat4x4& Bs, const Mat4x4& invBs, const Mat4x4& Bs_invBm,
                                           RayHit3& rh, RayCache3 *p_rc,
                                           const Context* p_context )
{
    RayHit3 tmp_rh;
#define __USE_RC_DCR_E_DS_BF_Barycentric //\todo This was used to test if barycentric RC worked properly, which does. KEEP it as reference. Non-Bruteforce BDT version uses barycentric exclusively
#ifdef __USE_RC_DCR_E_DS_BF_Barycentric
    // Transform ray to barycentric coords (\todo Only b[1,2] are used, b[0] is implicit, could reduce to 2x2 MxV
    Vec4 b0( invBs * mal::Concat(1,ray_pos+interval.Min()*ray_dir) );
    Vec4 b1( invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
    Vec4 dir_b( b1-b0 );
    Interval interval_b(0,1);
    Mat4x4 invBm( invBs * Bs_invBm ); //TEMP: should be stored in DCR.E or passed as a param
    rh.m_Interval.Set(1);
#else
    rh.m_Interval.Set(interval.Max());
#endif

    const DCR_TetSolidShape3::Element& dcr_e( p_dcr->m_vecE[eid] );
    for( unsigned int it_pie=0; it_pie<dcr_e.m_NumPatches; it_pie++ )
    {
        const DCR_TetSolidShape3::Patch& dcr_p( p_dcr->m_vecP[dcr_e.m_FirstPID+it_pie] );
        for( unsigned int it_tip=0; it_tip<dcr_p.m_NumTriangles; it_tip++ )
        {
            uint32 tid( dcr_p.m_FirstTID + it_tip );
#ifdef __USE_RC_DCR_E_DS_BF_Barycentric
            // Compute barycentric segment
            /* OPTIMIZATION \todo
               - This COULD be procomputed at the expense of A LOT OF extra memory in DCR.V... even if only b[1,2] are stored.
               - Also, a 2x2 MxV could be used because we NEVER use b[2]
            */
            Vec4 tri0_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(0) ]) );
            Vec4 tri1_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(1) ]) );
            Vec4 tri2_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(2) ]) );
            // Perform raycast in barycentric coords using two independent coords b[1],b[2] (b[0]=1-b[1]-b[2]-b[3])
            if( RayCast_Triangle3_DoubleSided( mal::GRange<1,3>(b0), mal::GRange<1,3>(dir_b),
                                               Interval(0,rh.m_Interval.Min()),
                                               mal::GRange<1,3>(tri0_b), mal::GRange<1,3>(tri1_b), mal::GRange<1,3>(tri2_b),
                                               tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                rh = tmp_rh;
#else
            // IN GLOBAL COORDS
            Vec3 tri0_0( mal::GRange<1,3>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(0) ]) ) );
            Vec3 tri1_0( mal::GRange<1,3>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(1) ]) ) );
            Vec3 tri2_0( mal::GRange<1,3>( Bs_invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(2) ]) ) );
            //s0<->s1 swapped because p is assumed INSIDE p_dcr \todo This MAY NOT BE TRUE as p1 is an IB1.EP.m_AvgPos that MAY BE OUTSIDE DCR2
            if( RayCast_Triangle3_DoubleSided( ray_pos, ray_dir,
                                               Interval(interval.Min(),rh.m_Interval.Min()), //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                                               tri0_0, tri1_0, tri2_0,
                                               tmp_rh, 0 ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                rh = tmp_rh;
#endif
        }
    }

#ifdef __USE_RC_DCR_E_DS_BF_Barycentric
    // Hit if any ray has actually hit a segment
    if( rh.m_Interval.Min() < interval_b.Max() )
    {
        // remap interval [0,1] to [min,max]
        rh.m_Interval.Min() = interval.Min() + rh.m_Interval.Min()*interval.Length();
        rh.m_Interval.Max() = interval.Min() + rh.m_Interval.Max()*interval.Length();
        // compute global and barycentric stuff
        rh.m_FeatureId = feature_id( eFT_Tetrahedron, eid );
        Vec4 b( mal::Concat(1-mal::Sum(rh.m_Point),rh.m_Point) ); //b[0] = 1-b[1]-b[2]-b[3]
        rh.m_Extra_BarycentricCoords = b;
        rh.m_Point = mal::GRange<1,3>( Bs * b );
        //\todo rh.m_Normal = transform vector to global?!?!
        return true;
    }
    else
        return false;
#else
    // Hit if any ray has actually hit a segment
    if( rh.m_Interval.Min() < interval.Max() )
    {
        // compute global and barycentric stuff
        rh.m_FeatureId = feature_id( eFT_Tetrahedron, eid );
        Vec4 b( invBs * mal::Concat(1,rh.m_Point) );
        rh.m_Extra_BarycentricCoords = b;
        return true;
    }
    else
        return false;
#endif
}

bool Clip_BarycentricRay3_BDOP( const Vec4& ray_pos_b, const Vec4& ray_dir_b, const Interval& interval_b,
                                const bv::BDOP4& bdop,
                                Interval& clipped_interval_b,
                                const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( overlap.m_Segment3_Tetrahedron3 ); //TEMP: To compare overlap/clip approaches
    clipped_interval_b = interval_b;
    // int first_hit_axis(0);
    for( int i=0; i<4; i++ )
    {
        if( mal::Abs(ray_dir_b[i]) < p_context->m_Epsilon_Dir )
        {
            if( ray_pos_b[i] < bdop[i].Min() || ray_pos_b[i] > bdop[i].Max() )
            {
                //Empty interval and return false
                clipped_interval_b = Interval::Empty();
                return false;
            }
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            Real inv_dir_b_i( mal::Rcp(ray_dir_b[i]) );
            Real lambda0( (bdop[i].Min()-ray_pos_b[i]) * inv_dir_b_i );
            Real lambda1( (bdop[i].Max()-ray_pos_b[i]) * inv_dir_b_i );
            if( lambda1 < lambda0 )
            {
                Real tmp( lambda0 );
                lambda0 = lambda1;
                lambda1 = tmp;
            }
            // Clip lambda-interval and update first-axis
            if( clipped_interval_b.Min() < lambda0 ) clipped_interval_b.Min() = lambda0;
            if( clipped_interval_b.Max() > lambda1 ) clipped_interval_b.Max() = lambda1;
            if( clipped_interval_b.IsEmpty() ) return false;
        }
    }
    return true;
}

/* Alternative ray-cast test using the TripleScalarProduct from RTCD 5.3.4 pg 186. */
inline bool RayCast_Triangle3_DoubleSided_TSP( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                               const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2,
                                               RayHit3& rh, RayCache3* p_rc = 0,
                                               const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Triangle3_DS );

    //\note Code from RTCD 5.3.6, pg 191
    Vec3 segment_p( ray_pos + ray_dir * interval.Min() );
    Vec3 segment_q( ray_pos + ray_dir * interval.Max() );
    Vec3 pq = segment_q - segment_p;
    Vec3 pa = triangle_p0 - segment_p;
    Vec3 pb = triangle_p1 - segment_p;
    Vec3 pc = triangle_p2 - segment_p;

    /* Single-sided
    Real u = mal::ScalarTriple(pq,pc,pb);
    if( u < 0 ) return false;
    Real v = mal::ScalarTriple(pq,pa,pc);
    if( v < 0 ) return false;
    Real w = mal::ScalarTriple(pq,pb,pa);
    if( w < 0 ) return false;
    */
    // Double-sided
    Real u = mal::ScalarTriple(pq,pc,pb);
    Real v = mal::ScalarTriple(pq,pa,pc);
    if( v*u < 0 ) return false;
    Real w = mal::ScalarTriple(pq,pb,pa);
    if( v*w < 0 || u*w < 0 ) return false;
    //END double-sided

    // If we arrive here, the ray hits the triangle, we compute a
    // double-sided GRayCast_Plane() to find hit lambda/point
    Vec3 triangle_n( mal::Cross( triangle_p1-triangle_p0, triangle_p2-triangle_p0 ) );
    //Test ray-plane
    Real n_dot_dir( mal::Dot(triangle_n,ray_dir) );
    Real n_dot_p_plus_d( mal::Dot(triangle_n,ray_pos - triangle_p0) );
    if( mal::Abs(n_dot_dir) > g_pDefaultContext->m_Epsilon_Dir ) // non-coplanar ray
    {
        Real lambda( -n_dot_p_plus_d / n_dot_dir );
        if( lambda >= interval.Min() && lambda <= interval.Max() )
        {
            rh.m_Interval = Interval(lambda,lambda);
            rh.m_Point = ray_pos + ray_dir * rh.m_Interval.Min();
            rh.m_Normal = mal::Normalized(triangle_n);
            rh.m_FeatureId = feature_id();
            rh.m_Extra_BarycentricCoords = Vec4::Zero(); //\todo Save barycentric coords according to FT
            return true;
        }
        else
            return false;
    }
    else if( mal::Abs(n_dot_p_plus_d) < g_pDefaultContext->m_Epsilon_Length ) // coplanar included ray, \lambda = any
    {
        rh.m_Interval = interval; //\todo CLIP IT TO TRIANGLE!!
        rh.m_Point = ray_pos + ray_dir*rh.m_Interval.Min();
        rh.m_Normal = mal::Normalized(triangle_n);
        rh.m_FeatureId = feature_id();
        return true;
    }
    else
        return false; // coplanar non-included ray
}

inline bool RayCast_Triangle3_DoubleSided_Fat( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                               const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2,
                                               RayHit3& rh, RayCache3* p_rc = 0,
                                               const Context* p_context = g_pDefaultContext )
{
    //\note Code from RTCD 5.3.6, pg 191
    Vec3 segment_p( ray_pos + ray_dir * interval.Min() );
    Vec3 segment_q( ray_pos + ray_dir * interval.Max() );
    Vec3 ab = triangle_p1 - triangle_p0;
    Vec3 ac = triangle_p2 - triangle_p0;
    Vec3 qp = segment_p - segment_q;
    Vec3 n = mal::Cross( ab, ac );
    Real d = mal::Dot( qp, n );
    // Coplanar?
    //if( mal::Abs(d) < g_pDefaultContext->m_Epsilon_Dir ) return false;
    if( d <= Real(0) )
    {
        // Invert triangle orientation
        n = -n;
        ab = triangle_p2 - triangle_p0;
        ac = triangle_p1 - triangle_p0;
        d = -d;
    }
    Vec3 ap = segment_p - triangle_p0;
    Real t = mal::Dot( ap, n );
    if( t < Real(-p_context->m_Epsilon_Length) ) return false;
    if( t > d+p_context->m_Epsilon_Length ) return false;
    Vec3 e = mal::Cross( qp, ap );
    Real v = mal::Dot( ac, e );
    if( v < Real(-p_context->m_Epsilon_Length) || v > d+p_context->m_Epsilon_Length ) return false;
    Real w = -mal::Dot( ab, e );
    if( w < Real(-p_context->m_Epsilon_Length) || v + w > d+p_context->m_Epsilon_Length ) return false;
    // Hit found, fill rh and return true
    t /= d;
    rh.m_Interval = Interval( interval.Min() + t * interval.Length() );
    rh.m_Point = ray_pos + rh.m_Interval.Min() * ray_dir;
    rh.m_Normal = mal::Normalized( n );
    //\todo THIS MAY NEED to consider inverted triangles (_DoubleSided)
    rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E/F sub-features
    rh.m_Extra_BarycentricCoords = Vec4::Zero(); //\todo Save barycentric coords according to FT
    return true;
}


/* Performn RayCast on a DCR.E using per-DCR.P BDOP-Tree
   - Only closest hit is reported
   - Computed in barycentric coords to avoid transforming DCR.E.V, but RayHit results in global coords as usual

   \todo OPTIMIZATION
   - Subst std::vector stack with std::array
   - Avoid Bs and related overhead per-triangle
   - Avoid calling Concat and GRange repeatedly per triangle
*/
bool RayCast_DCR_E_DoubleSided_BDT(const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                   const DCR_TetSolidShape3* p_dcr, uint32 eid,
                                   const Mat4x4& Bs, const Mat4x4& invBs, const Mat4x4& Bs_invBm,
                                   RayHit3& rh, RayCache3 *p_rc,
                                   const Context* p_context )
{
    RayHit3 tmp_rh;
// #define __ENABLE_OFFSET_RH //\todo THIS IS NOT WORKING and NOT ENOUGH
#ifdef __ENABLE_OFFSET_RH
    RayHit3 offset_rh;
    offset_rh.m_Interval.Set(1);
    // Real thickness_b( p_context->m_DCR2DCR_RC_Thickness );//\todo convert to bary, along ray_dir... should not require sqrt...
    Real thickness_b( mal::Norm( mal::GRange<1,1,3,3>(invBs)*(p_context->m_DCR2DCR_RC_Thickness*ray_dir) ) );
    // Real thickness_b_sq( mal::Det(mal::GRange<1,1,3,3>(invBs)) * mal::Sq(p_context->m_DCR2DCR_RC_Thickness) );
    // Real error( mal::Abs( thickness_b_sq - mal::Sq(thickness_b) ) );
    // GEO_LOG_WARNING_IF( error > 1e-6, "=>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ERROR= %f", error );
#endif
    // Transform ray to barycentric coords (\todo Only b[1,2] are used, b[0] is implicit, could reduce to 2x2 MxV
#ifdef __ENABLE_OFFSET_RH
    Vec4 b0( invBs * mal::Concat(1,ray_pos+(interval.Min()-p_context->m_DCR2DCR_RC_Thickness)*ray_dir) );
    Vec4 b1( invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
#else
    Vec4 b0( invBs * mal::Concat(1,ray_pos+interval.Min()*ray_dir) );
    Vec4 b1( invBs * mal::Concat(1,ray_pos+interval.Max()*ray_dir) );
#endif
    Vec4 dir_b( b1-b0 );
    Interval interval_b(0,1);
    Mat4x4 invBm( invBs * Bs_invBm ); //TEMP: should be stored in DCR.E or passed as a param
    rh.m_Interval.Set(1); //barycentric ray is a segment with lambda [0,1] that maps the euclidean segment interval [min,max]

    const DCR_TetSolidShape3::Element& dcr_e( p_dcr->m_vecE[eid] );
    for( unsigned int it_pie=0; it_pie<dcr_e.m_NumPatches; it_pie++ )
    {
        const DCR_TetSolidShape3::Patch& dcr_p( p_dcr->m_vecP[dcr_e.m_FirstPID+it_pie] );
        // Init stack
        typedef std::pair< DCR_TetSolidShape3::Patch::bdt_node_tip_range,Interval> stack_entry_type;
        std::vector< stack_entry_type > stackBDTN;
        stackBDTN.push_back( stack_entry_type(DCR_TetSolidShape3::Patch::bdt_node_tip_range(0,dcr_p.m_NumTriangles),interval_b) );
        // GEO_LOG("DCR.P[%u].BDT",it_p);
        while( !stackBDTN.empty() )
        {
            // Pop PS subarray
            stack_entry_type se( stackBDTN.back() );
            stackBDTN.pop_back();
            const DCR_TetSolidShape3::Patch::bdt_node_tip_range node_tr( se.first );
            const uint32 length_tr( node_tr.second - node_tr.first );
            if( length_tr <= p_context->m_RCDCR_BDT_MaxLeafSize )
            {
                // Test all triangles in subtree
                for( unsigned int it_sin=0; it_sin<length_tr; it_sin++ )
                {
                    uint32 tid( dcr_p.m_FirstTID + node_tr.first + it_sin );
                    /* Compute barycentric triangle
                       \todo OPTIMIZATION
                       - This COULD be procomputed at the expense of A LOT OF extra memory in DCR.V... even if only b[1,2] are stored.
                       - Also, a 2x2 MxV could be used because we NEVER use b[2]
                    */
                    Vec4 tri0_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(0) ]) );
                    Vec4 tri1_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(1) ]) );
                    Vec4 tri2_b( invBm * mal::Concat(1,p_dcr->m_vecV[ p_dcr->m_vecT[tid].GetVID(2) ]) );
                    // Perform raycast in barycentric coords using two independent coords b[1],b[2] (b[0]=1-b[1]-b[2]-b[3])
                    if( RayCast_Triangle3_DoubleSided( mal::GRange<1,3>(b0), mal::GRange<1,3>(dir_b),
                                                       Interval(0,rh.m_Interval.Min()),
                                                       mal::GRange<1,3>(tri0_b), mal::GRange<1,3>(tri1_b), mal::GRange<1,3>(tri2_b),
                                                       tmp_rh, 0, p_context ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                    {
#ifdef __ENABLE_OFFSET_RH
                        // rh is the closest hit FARTHER than the offset
                        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                            && tmp_rh.m_Interval.Min() >= thickness_b )
                            rh = tmp_rh;
                        // offset_rh is the closest hit to the offset on ANY side
                        if( mal::Abs( tmp_rh.m_Interval.Min() - thickness_b ) < mal::Abs( offset_rh.m_Interval.Min() - thickness_b ) )
                            offset_rh = tmp_rh;
#else
                        rh = tmp_rh;
#endif
                    }
                }
            }
            else
            {
                // Get node/triangle data
                const DCR_TetSolidShape3::Triangle& bdtn( p_dcr->m_vecT[dcr_p.m_FirstTID+node_tr.first] );
                const bv::BDOP4 bdop( bdtn.BDTN_BDOPq() );
                // Cut recursive ray interval to current best hit
                if( se.second.Max() > rh.m_Interval.Min() ) se.second.Max() = rh.m_Interval.Min();
                // Test BDOP and clip interval
                Interval clipped_interval_b;
                if( Clip_BarycentricRay3_BDOP( b0, dir_b, se.second,
                                               bdop,
                                               clipped_interval_b ) ) //\todo IDEALLY, do this directly in BDOP3q space
                {
                    // Test node/triangle
                    /* Compute barycentric triangle
                       \todo OPTIMIZATION
                       - This COULD be procomputed at the expense of A LOT OF extra memory in DCR.V... even if only b[1,2] are stored.
                       - Also, a 2x2 MxV could be used because we NEVER use b[2]
                    */
                    Vec4 tri0_b( invBm * mal::Concat(1,p_dcr->m_vecV[ bdtn.GetVID(0) ]) );
                    Vec4 tri1_b( invBm * mal::Concat(1,p_dcr->m_vecV[ bdtn.GetVID(1) ]) );
                    Vec4 tri2_b( invBm * mal::Concat(1,p_dcr->m_vecV[ bdtn.GetVID(2) ]) );
                    // Perform raycast in barycentric coords using two independent coords b[1],b[2] (b[0]=1-b[1]-b[2]-b[3])
                    if( RayCast_Triangle3_DoubleSided( mal::GRange<1,3>(b0), mal::GRange<1,3>(dir_b),
                                                       Interval(0,rh.m_Interval.Min()),
                                                       mal::GRange<1,3>(tri0_b), mal::GRange<1,3>(tri1_b), mal::GRange<1,3>(tri2_b),
                                                       tmp_rh, 0, p_context ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                    {
#ifdef __ENABLE_OFFSET_RH
                        // rh is the closest hit FARTHER than the offset
                        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min()
                            && tmp_rh.m_Interval.Min() >= thickness_b )
                            rh = tmp_rh;
                        // offset_rh is the closest hit to the offset on ANY side
                        if( mal::Abs( tmp_rh.m_Interval.Min() - thickness_b ) < mal::Abs( offset_rh.m_Interval.Min() - thickness_b ) )
                            offset_rh = tmp_rh;
#else
                        rh = tmp_rh;
#endif
                    }
                    /*
                    if( RayCast_Triangle3_DoubleSided( mal::GRange<1,3>(b0), mal::GRange<1,3>(dir_b),
                                                       Interval(0,rh.m_Interval.Min()),
                                                       mal::GRange<1,3>(tri0_b), mal::GRange<1,3>(tri1_b), mal::GRange<1,3>(tri2_b),
                                                       tmp_rh, 0, p_context ) ) //\note this guarantees: tmp_rh.m_Interval.Min() < rh.m_Interval.Min() on hit
                        rh = tmp_rh;
                    */
                    // Recurse node using BDOP-clipped barycentric ray interval, so that all BCoords span shorten due to clipping
                    int remaining_tr( length_tr - 1 );
                    if( remaining_tr > 1 )
                    {
                        stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr.first+1, node_tr.first+1+remaining_tr/2 ), clipped_interval_b ) );
                        stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr.first+1+remaining_tr/2, node_tr.second ), clipped_interval_b ) );
                    }
                    else if( remaining_tr == 1 )
                        stackBDTN.push_back( stack_entry_type( DCR_TetSolidShape3::Patch::bdt_node_tip_range( node_tr.first+1, node_tr.first+2 ), clipped_interval_b ) );
                }
            }
        }
    }

    // Hit if any ray has actually hit a triangle
    if( rh.m_Interval.Min() < interval_b.Max() )
    {
#ifdef __ENABLE_OFFSET_RH //BUG: Intervals are modified by the thickness, so we SHOULD undo thickness effect.
        // If offset_rh has a valid hit closer to ray_pos than rh, report that one instead, EVEN if behind ray_pos
        if( mal::Abs(offset_rh.m_Interval.Min()-thickness_b) < mal::Abs(rh.m_Interval.Min()-thickness_b) )
            rh = offset_rh;
        // remap interval [0,1] to [min,max]
        rh.m_Interval.Min() = interval.Min() + rh.m_Interval.Min()*interval.Length();
        rh.m_Interval.Max() = interval.Min() + rh.m_Interval.Max()*interval.Length();
#else
        // remap interval [0,1] to [min,max]
        rh.m_Interval.Min() = interval.Min() + rh.m_Interval.Min()*interval.Length();
        rh.m_Interval.Max() = interval.Min() + rh.m_Interval.Max()*interval.Length();
#endif
        // compute global and barycentric stuff
        rh.m_FeatureId = feature_id( eFT_Tetrahedron, eid );
        Vec4 b( mal::Concat(1-mal::Sum(rh.m_Point),rh.m_Point) ); //b[0] = 1-b[1]-b[2]-b[3]
        rh.m_Extra_BarycentricCoords = b;
        rh.m_Point = mal::GRange<1,3>( Bs * b );
        //\todo rh.m_Normal = transform vector to global?!?!
        return true;
    }
    else
        return false;
}

}} //namespace geo::np
