#include "Contact_X_Y_Stochastic.h"
#include <Geo/shape/IShape.h>
#include <Mal/GMatEigen.h>
#include <vector>
#include "stats.h"

#define __USE_CACHED_POB_DATA

#define __ENABLE_IM
#ifdef __ENABLE_IM
//#  define __ENABLE_IM_USE_BOUNDARIELS //NOT FOR A LOOONG TIME
#  define __ENABLE_IM_HIT_IP //check IP
//#  define __ENABLE_IM_ADVANCE_CLOSEST //otherwise __ENABLE_IM_ADVANCE_PERPENDICULAR, which is usually better
#  define __ENABLE_IM_BACKTRACK_AFTER_CROSSING //necessary
#  define __ENABLE_IM_CD
#  define __ENABLE_IM_CD_AREA_SCALE_DEPTH
//#  define __ENABLE_TRACE_IM
#endif

//#define __ENABLE_TRACE_PAIRS
//#define __ENABLE_TRACE_PAIRS_IP

#define __ENABLE_COMPARE_WITH_MESH2_MESH2_BRUTEFORCE
#ifdef __ENABLE_COMPARE_WITH_MESH2_MESH2_BRUTEFORCE
#  include <Geo/shape/MeshSolidShape2.h>
#  include "Contact_X_Y_BruteForce.h"
#endif

//#define __TRACE_PCA
#ifdef __TRACE_PCA
#  include <Mal/GSerialization.h>
#  include <iostream>
#endif

namespace geo {
namespace np {

class IntersectionMapping;

//\todo this should be in geo_params.h or something similar, and maybe be tunable through an Itemstream
struct PARAMS
{
    struct DEBUG
    {
        struct NP
        {
            struct CONTACT
            {
                struct STOCHASTIC
                {
                    Flags32 DDF;
                    enum EDDF { eDraw_Nothing = 0,
                                eDraw_IP      = 1<<0,
                                eDraw_CP      = 1<<1,
                                eDraw_NP      = 1<<2,
                                eDraw_RNP     = 1<<3,
                                eDraw_PCA     = 1<<4,
                                eDraw_IM      = 1<<5,
                                eDraw_Default = eDraw_IP | /*eDraw_NP  |*/ eDraw_IM
                    };
                    STOCHASTIC() : DDF(eDraw_Default) {}
                } stochastic;
            } contact;
        } np;
    } debug;
} params;

//--------------------------------------------------------------------------------------------------------------------------------
// Stochastic Cache
//--------------------------------------------------------------------------------------------------------------------------------

class StochasticSPC2: public ISpecificPairwiseCache2
{
public:
    finline StochasticSPC2( IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 )
    : m_pDS1(p_ds1), m_pDS2(p_ds2)
        {
            //\todo
        }
    finline ~StochasticSPC2()
        {
            //\todo dealloc from pool!!
            delete m_pDS1;
            delete m_pDS2;
        }

    finline void Clear()
        {
            /*\todo This deletes POB directlty regardless of their
              reference counts, which by now results in a LOG_ERROR in
              GIndexedRefCountedPoolDA... if it's considered safe, add
              a param to ClearPOB() that disables refcount checking in
              GIndexedRefCountedPoolDA::Clear() too
            */
            for( unsigned int it_ip=0; it_ip < m_vecIP.size(); it_ip++ )
                m_vecIP[it_ip].DiscardIM( m_pDS1, m_pDS2 );
            if(m_pDS1) m_pDS1->ClearPOB();
            if(m_pDS2) m_pDS2->ClearPOB();
            m_vecIP.clear(); m_vecCP.clear(); m_vecNP.clear();
#ifdef __USE_CACHED_POB_DATA
            m_vecP1_2.clear(); m_vecN1_2.clear();
            m_vecP2_2.clear(); m_vecN2_2.clear();
#endif
        }

public:
    IDomainSampler2 *m_pDS1;
    IDomainSampler2 *m_pDS2;

    struct Pair
    {
        Real m_DistSq;
        point_on_boundary_id_type m_POB1;
        point_on_boundary_id_type m_POB2;
        finline bool operator<( const Pair &other ) const { return m_DistSq < other.m_DistSq; }
        finline Pair() : m_DistSq(0), m_POB1(cInvalidPOB), m_POB2(cInvalidPOB) {}
        finline Pair( const Pair &pair ) : m_DistSq(pair.m_DistSq), m_POB1(pair.m_POB1), m_POB2(pair.m_POB2) {}
        finline Pair( point_on_boundary_id_type pob1, point_on_boundary_id_type pob2, Real dist_sq )
        : m_DistSq(dist_sq), m_POB1(pob1), m_POB2(pob2) {}
    };

    typedef Pair ClosestPair;
    typedef Pair NearPair;
    struct IntersectionPair: public Pair
    {
        bool m_bFlooded;
        IntersectionMapping *m_pIM;
        finline IntersectionPair( const Pair &pair ) : Pair(pair), m_bFlooded(false), m_pIM(0) {}

        finline bool IsFlooded() const { return m_bFlooded; }
        finline void SetFlooded( bool b_flooded ) { m_bFlooded = b_flooded; }

        finline bool HasIM() const { return 0 != m_pIM; }
        finline void SetIM( IntersectionMapping *p_im ) { GEO_ASSERT(!HasIM()); m_pIM = p_im; }
        finline IntersectionMapping *GetIM() { return m_pIM; }

        bool MergeIM( IntersectionPair &ip, IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
        void DiscardIM( IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
    };

    std::vector<IntersectionPair> m_vecIP; //Definitely intersecting pairs
    std::vector<ClosestPair> m_vecCP; //Closest pairs below NearDist
    std::vector<NearPair> m_vecNP; //Near pairs below NearDist

#ifdef __USE_CACHED_POB_DATA
    //\name Cached data for all POB
    //@{
    std::vector<Vec2> m_vecP1_2,m_vecN1_2;
    std::vector<Vec2> m_vecP2_2,m_vecN2_2;
    //@}
#endif

    //\todo declare global static POOL
};

/* \todo NO, do it lazily, check first one-way and, if not behind, AVOID retrieving the OTHER normal
finline bool AreBehind( const Vec2 &p1, const Vec2 &n1,
                        const Vec2 &p2, const Vec2 &n2 )
{
    return (p1-p2)*n2 < 0 && (p2-p1)*n1 < 0;
}
*/

#ifdef __ENABLE_IM
//! Segment2 vs Segment 2D => (a1,b1) vs (a2,b2)
inline bool Intersection_Segment2_Segment2( const Vec2 &a1, const Vec2 &b1,
                                            const Vec2 &a2, const Vec2 &b2,
                                            //\todo Real epsilon_length,
                                            Real &lambda1, Real &lambda2 )
{
    Vec2 d1( b1 - a1 );
    Vec2 d2( b2 - a2 );
    Real det( d1[1]*d2[0] - d1[0]*d2[1] );
    if( mal::Abs(det) > Real(0.00001f) )
    {
        Vec2 diff_ab( a2 - a1 );
        Real inv_det( mal::Rcp(det) );
        lambda1 = inv_det * (diff_ab[1]*d2[0] - diff_ab[0]*d2[1]);
        lambda2 = inv_det * (diff_ab[1]*d1[0] - diff_ab[0]*d1[1]);
        return lambda1 >= Real(0) && lambda1 <= Real(1)
            && lambda2 >= Real(0) && lambda2 <= Real(1);
    }
    else // \todo Parallel, may be partially coincident, we MUST handle it
    {
        /*\todo We'll need to discriminate between coincident cases
          CROSSING | INSIDE & NONCROSSING | OUTSIDE & NONCROSSING
          within a certain epsilon_length, because we're interested in
          shallow penetrations, where segments will tend to ALIGN
          while inside or outside (depending on contact response
          method)... this behaviour should be CONSISTENT, because
          incremental IM computation requires finding an outwards
          crossing for each inwards one.
        */
        return false;
    }
}

/* Fast test, accepts overlap if point is projected inside
   [-thickness, length+thickness] along segment direction and is
   perpendicularly closer than thickness.
   TEMP: Currently unsused
bool FastTest_Segment2_Point2( const Vec2 &p1, const Vec2 &p2, const Vec2 &point, Real thickness, Real &lambda )
{
    Vec2 u( p2 - p1 );
    Vec2 v( point - p1 );
    Real length( u.NormAndNormalize() );
    Real dot_u_v( mal::Dot( u, v ) );
    if( dot_u_v < -thickness //behind
        || dot_u_v > Real(length) + thickness //in-front
        ||  mal::NormSq( v - dot_u_v * u ) > mal::Sq(thickness) ) //perpendularly far
        return false;
    else
    {
        lambda = mal::Clamp01( dot_u_v / length );
        return true;
    }
}
*/

/* Forward segment fast test, accepts overlap if point is projected
   inside (+thickness, length+thickness] along segment direction and
   is perpendicularly closer than thickness.

   \note Forward segment ignores interval [0,thickness), it is tested
   by the previous forward segment.
*/
bool FastTest_ForwardSegment2_Point2( const Vec2 &p1, const Vec2 &p2, const Vec2 &point, Real thickness, Real &lambda )
{
    Vec2 u( p2 - p1 );
    Vec2 v( point - p1 );
    Real length( u.NormAndNormalize() );
    Real dot_u_v( mal::Dot( u, v ) );
    if( dot_u_v < thickness //after previous in-front
        || dot_u_v > Real(length) + thickness //in-front
        ||  mal::NormSq( v - dot_u_v * u ) > mal::Sq(thickness) ) //perpendularly far
        return false;
    else
    {
        lambda = mal::Clamp01( dot_u_v / length );
        return true;
    }
}

/* Test the segment against all IP, computing first hit lambda and
   flooding any hit IP

   \note We use Forward Segment tests to avoid hitting IP "backwards"
   and stopping unnecessarily. This is critical at the beginning of
   the flood, where redundant IP may exist.
*/
bool Test_Segment_IP( const Vec2 &a_2, const Vec2 &b_2,
                      std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int ignored_ip_index,
                      const Transform2 &tr12,
                      IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
                      Real thickness,
                      Real &min_lambda )
{
    min_lambda = Real(2); //>1
    for( unsigned int it_ip=0; it_ip<vec_ip.size(); it_ip++ )
        if( it_ip != ignored_ip_index )
        {
            StochasticSPC2::IntersectionPair &ip( vec_ip[it_ip] );
            Vec2 midpoint_ip_2( Real(0.5) * ( tr12 * p_ds1->POB_Position(ip.m_POB1) + p_ds2->POB_Position(ip.m_POB2) ) );
            Real lambda(-1);
            if( FastTest_ForwardSegment2_Point2( a_2, b_2, midpoint_ip_2, thickness /*+ mal::Sqrt(ip.m_DistSq)*/, lambda ) )
            {
                ip.SetFlooded(true);
                if( lambda < min_lambda) min_lambda = lambda;
#ifdef __ENABLE_TRACE_IM
                GEO_LOG( "IM::TSIP() Hit IP[%d] mp = (%f,%f) with l=%f", it_ip, midpoint_ip_2[0], midpoint_ip_2[1], lambda );
#endif
            }
            else if( FastTest_ForwardSegment2_Point2( a_2, b_2, midpoint_ip_2, 4*thickness, lambda ) )
            {
#ifdef __ENABLE_TRACE_IM
                GEO_LOG("IM::TSIP() Near-Miss IP[%d] mp = (%f,%f) is actually Hit with at l=%f with 4x thickness %f",
                        it_ip, midpoint_ip_2[0], midpoint_ip_2[1], lambda, thickness );
#endif
            }
        }
    return min_lambda <= Real(1);
}

//----- Heuristic "affinity" between boundariels
typedef Real (*HeuristicFunc)( const Vec2 &p1, const Vec2 &n1, Real r1,
                               const Vec2 &p2, const Vec2 &n2, Real r2 );

Real Heuristic_DistSq( const Vec2 &p1, const Vec2 &n1, Real r1,
                       const Vec2 &p2, const Vec2 &n2, Real r2 )
{
    return mal::NormSq( p2 - p1 );
}
Real Heuristic_DistSq_AngleSq( const Vec2 &p1, const Vec2 &n1, Real r1,
                               const Vec2 &p2, const Vec2 &n2, Real r2 )
{
    Vec2 v( p2 - p1 );
    return mal::Sq( mal::Dot( v, mal::PerpendicularCW(n1) ) )
        + mal::Sq( mal::Dot( v, mal::PerpendicularCW(n2) ) );
}
Real Heuristic_AngleSq( const Vec2 &p1, const Vec2 &n1, Real r1,
                        const Vec2 &p2, const Vec2 &n2, Real r2 )
{
    Vec2 v( p2 - p1 );
    Real norm_sq_v( mal::NormSq(v) );
    GEO_ASSERT( norm_sq_v > Real(0) );
    return mal::Rcp( norm_sq_v )
        * ( mal::Sq( mal::Dot( v, mal::PerpendicularCW(n1) ) )
            + mal::Sq( mal::Dot( v, mal::PerpendicularCW(n2) ) ) );
}


//----------------------------------------------------------------
//-- IM  Class
//----------------------------------------------------------------
#ifdef __ENABLE_IM_USE_BOUNDARIELS
/* IM data and operations
class IntersectionMapping
{
public:
    IntersectionMapping() {}

    void Clear() { m_IB1.clear(); m_IB2.clear(); m_IM.clear(); }
    bool Init( std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int seed_ip_index,
               const Vec2 &pos1_2, point_on_boundary_id_type neighbour_pob1, const Vec2 &neighbour_pos1_2,
               const Vec2 &pos2_2, point_on_boundary_id_type neighbour_pob2, const Vec2 &neighbour_pos2_2,
               const Transform2 &tr12, const Transform2 &tr2,
               IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
               Real crossing_neighbourhood_radius, Real step_length, Real step_length_epsilon_abs
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
               , ContactData2::VizData &vd
#endif
               );
    bool Update( std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int seed_ip_index,
                 IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
    void Destroy( IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
    bool ComputeContactData( ContactData2 &cd, const Transform2 &tr2, IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 ) const;

private:
    struct BoundaryElement: public GIDomainSamplerD<2>::boundariel_type
    {
        // TEMP: pos, normal, radius Inherited from boundariel_type
        point_on_boundary_id_type m_POB;
        Interval m_Interval; //For clipped boundariels, extents around pos reduced to [min*tangent, pos+max*tangent] (min,max may be negative)
        finline BoundaryElement( point_on_boundary_id_type pob, const Vec2 &pos_2, const Vec2 &normal, Real radius )
        : m_POB(pob), m_Pos_2(midpos_2), m_Normal(normal) m_Radius(radius) {}
        finline BoundaryElement() : m_POB(cInvalidPOB) {}
    };
    struct Pair
    {
        uint32 m_Id1;
        uint32 m_Id2;
        Real m_H;
        finline Pair( uint32 id1, uint32 id2, Real heuristic ) : m_Id1(id1), m_Id2(id2), m_H(heuristic) {}
    };
    TODOOOOO adapt all pob-segment based code to pob-boundariel based
};
*/
#else //__ENABLE_IM_USE_BOUNDARIELS

// IM data and operations
class IntersectionMapping
{
public:
    IntersectionMapping() {}

    void Clear() { m_IB1.clear(); m_IB2.clear(); m_IM.clear(); }
    bool Init( std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int seed_ip_index,
               const Vec2 &pos1_2, point_on_boundary_id_type neighbour_pob1, const Vec2 &neighbour_pos1_2,
               const Vec2 &pos2_2, point_on_boundary_id_type neighbour_pob2, const Vec2 &neighbour_pos2_2,
               const Transform2 &tr12, const Transform2 &tr2,
               IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
               Real crossing_neighbourhood_radius, Real step_length, Real step_length_epsilon_abs, Flags32 im_flags
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
               , ContactData2::VizData &vd
#endif
               );
    bool Update( std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int seed_ip_index,
                 IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
    void Destroy( IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
    bool ComputeContactData( ContactData2 &cd,
                             const Transform2 &tr2, IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
                             Real contact_depth_epsilon ) const;

private:
    struct BoundaryElement
    {
        point_on_boundary_id_type m_POB;
        Vec2 m_Pos_2;
        Vec2 m_Normal_2;
        Real m_Radius;
        uint8 m_Valence;
        finline BoundaryElement( point_on_boundary_id_type pob, const Vec2 &pos_2, Real radius = Real(0) )
        : m_POB(pob), m_Pos_2(pos_2), m_Normal_2(0,0), m_Radius(radius), m_Valence(0) {}
        finline BoundaryElement() : m_POB(cInvalidPOB), m_Valence(0) {}
    };
    struct Pair
    {
        uint32 m_Id1;
        uint32 m_Id2;
        Real m_H;
        finline Pair( uint32 id1, uint32 id2, Real heuristic ) : m_Id1(id1), m_Id2(id2), m_H(heuristic) {}
    };

private:

    //! Add normals to IB entries, needed by some remapping processes
    void ComputeNormals( const Transform2 &tr12, IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 );
    //! Count number of appearances of each BE in the mapping
    void ComputeValences();

    //\note Mapping Post-processes
    //@{
    void ReverseMapping();  //Reverse mapping
    void ReduceMapping();   //Remove redundant mappings
    void ProjectMapping( const Transform2 &tr12, const Transform2 &tr21,
                         IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
                         Real neighbourhood_radius ); //Project POB on opposite IB
    void FixMapping(); //\todo Remove degenerated mappings
    void RefineMapping() {};   //\todo Add intermediate POB where needed
    void OptimizeMapping() {}; //\todo Move existing POB to improve mapping

    // Total remappings
    void MinHeuristicPairRemapping( HeuristicFunc HF ); //Pure best heuristic global mapping, to compare
    //@}

    // Compute first intersection between Segment2 and IB2
    bool Intersection_Segment_IB( const Vec2 &a1, const Vec2 &b1,
                                  const std::vector<BoundaryElement> &vec_ibe2,
                                  //\todo Real epsilon_length,
                                  Real &lambda1, Real &lambda2, unsigned int &index_ibe2 );

private:
    std::vector< BoundaryElement > m_IB1;
    std::vector< BoundaryElement > m_IB2;
    std::vector< Pair > m_IM; //mapping between IB1 and IB2, first/second indices relative to ib1/ib2
};

bool IntersectionMapping::Init( std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int seed_ip_index,
                                const Vec2 &pos1_2, point_on_boundary_id_type neighbour_pob1, const Vec2 &neighbour_pos1_2,
                                const Vec2 &pos2_2, point_on_boundary_id_type neighbour_pob2, const Vec2 &neighbour_pos2_2,
                                const Transform2 &tr12, const Transform2 &tr2,
                                IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
                                Real crossing_neighbourhood_radius, Real step_length, Real step_length_epsilon_abs, Flags32 im_flags
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                                , ContactData2::VizData &vd
#endif
                                )
{
    GEO_ASSERT( !vec_ip[seed_ip_index].IsFlooded() );
    GEO_ASSERT( 0 == m_IB1.size() && 0 == m_IB2.size() &&0 == m_IM.size() );

#ifdef __ENABLE_TRACE_IM
    GEO_LOG( "IM::Init() IP[%d] = (%d,%d)", seed_ip_index, vec_ip[seed_ip_index].m_POB1, vec_ip[seed_ip_index].m_POB2 );
#endif

    // Init shortcuts
    const StochasticSPC2::IntersectionPair &ip( vec_ip[seed_ip_index] );
    std::vector< BoundaryElement > &ib1( m_IB1 );
    std::vector< BoundaryElement > &ib2( m_IB2 );
    std::vector< Pair > &im( m_IM );

    // Flood seed_ip
    vec_ip[seed_ip_index].SetFlooded(true);

    //-- Init Intersection Boundary stepping
    //IB1
    ib1.push_back( BoundaryElement( ip.m_POB1, pos1_2, 0 ) );
    p_ds1->POB_IncRef( ip.m_POB1 );
    ib1.push_back( BoundaryElement( neighbour_pob1, neighbour_pos1_2, Real(0.5) * mal::Norm(pos1_2 - neighbour_pos1_2) ) );
    p_ds1->POB_IncRef( neighbour_pob1 );
    //IB2
    ib2.push_back( BoundaryElement( ip.m_POB2, pos2_2, 0 ) );
    p_ds2->POB_IncRef( ip.m_POB2 );
    ib2.push_back( BoundaryElement( neighbour_pob2, neighbour_pos2_2, Real(0.5) * mal::Norm(pos2_2 - neighbour_pos2_2) ) );
    p_ds2->POB_IncRef( neighbour_pob2 );

    //-- Init Intersection Mapping
    im.push_back( Pair(0,0,0) );
    /* \todo This mapping is NOT correct if edges cross at a non-acute angle
       im.push_back( std::make_pair(1,1) );
       //\note If we add it, no lookahead POB exists at init, as the IM uses last POB in IB1/IB2
       */
    //-- Flood IM OLD: Incremental IB/IM with lookahead
    bool bCrossing1(false), bCrossing2(false);
    Real lambda1(-1), lambda2(-1);
    unsigned int index_ibe1(cInvalidFeatureIndex), index_ibe2(cInvalidFeatureIndex);
    const unsigned int cMaxStepsIM(1000);
    unsigned int num_steps_im(0);
#ifdef __ENABLE_IM_HIT_IP
    bool bHitIP1(false), bHitIP2(false);
    Real lambda_hit_ip1(-1), lambda_hit_ip2(-1);
#  ifdef __ENABLE_TRACE_IM
    GEO_LOG( "IM::Init() Pre-Flooding IP[%d]", seed_ip_index );
#  endif
    bHitIP1 = Test_Segment_IP( pos1_2, neighbour_pos1_2, vec_ip, seed_ip_index,
                               tr12, p_ds1, p_ds2,
                               crossing_neighbourhood_radius,
                               lambda_hit_ip1 );
    bHitIP2 = Test_Segment_IP( pos2_2, neighbour_pos2_2, vec_ip, seed_ip_index,
                               tr12, p_ds1, p_ds2,
                               crossing_neighbourhood_radius,
                               lambda_hit_ip2 );
#  ifdef __ENABLE_TRACE_IM
    GEO_LOG( "IM::Init() Flooding IP[%d], bHitIP = %d,%d", seed_ip_index, bHitIP1, bHitIP2 );
#  endif
#endif
    do
    {
        // Get lookahead boundariel geometry => segments (a1,b1) and (a2,b2)
        Vec2 a1_2( ib1[ib1.size()-2].m_Pos_2 );
        Vec2 b1_2( ib1.back().m_Pos_2 );
        Vec2 a2_2( ib2[ib2.size()-2].m_Pos_2 );
        Vec2 b2_2( ib2.back().m_Pos_2 );

        //-- Advance lookahead POB if necessary and check potential crossing incrementally
        //IB1
        bool bAdvance1( im.back().m_Id1 == ib1.size()-1 );
        if( bAdvance1
#ifdef __ENABLE_IM_HIT_IP
            && !bHitIP1
#endif
            )
        {
            Vec2 dir1_1( tr12.m_Rot.Transposed() * ( b1_2 - a1_2 ) );
            point_on_boundary_id_type spob = p_ds1->StepPOB( ib1.back().m_POB, dir1_1, step_length, step_length_epsilon_abs );
            a1_2 = b1_2;
            b1_2 = tr12 * p_ds1->POB_Position( spob );
            GEO_ASSERT( mal::NormSq(a1_2-b1_2) > 0 );
            ib1.push_back( BoundaryElement( spob, b1_2, Real(0.5) * mal::Norm(a1_2 - b1_2) ) ); //Ref == 1
            // Check crossing of the new segment against the whole boundary2
            bCrossing1 = Intersection_Segment_IB( a1_2, b1_2, ib2, lambda1, lambda2, index_ibe2 );
#ifdef __ENABLE_IM_HIT_IP
            // Test IP and flood any hit one
            bHitIP1 = Test_Segment_IP( a1_2, b1_2, vec_ip, seed_ip_index,
                                       tr12, p_ds1, p_ds2,
                                       crossing_neighbourhood_radius,
                                       lambda_hit_ip1 );
#  ifdef __ENABLE_TRACE_IM
            if( bCrossing1 || bHitIP1 )
            {
                GEO_LOG_WARNING( "IM::Init() bHitIP = %d,%d and bCrossing = %d,%d at #steps = %d",
                                 bHitIP1, bHitIP2, bCrossing1, bCrossing2, num_steps_im );
            }
#  endif
#endif
        }
        else
            bCrossing1 = false;
        //IB2
        bool bAdvance2( !bCrossing1 && im.back().m_Id2 == ib2.size()-1 );
        if( bAdvance2
#ifdef __ENABLE_IM_HIT_IP
            && !bHitIP2
#endif
            )
        {
            Vec2 dir2_2( b2_2 - a2_2 );
            point_on_boundary_id_type spob = p_ds2->StepPOB( ib2.back().m_POB, dir2_2, step_length, step_length_epsilon_abs );
            a2_2 = b2_2;
            b2_2 = p_ds2->POB_Position( spob );
            GEO_ASSERT( mal::NormSq(a2_2-b2_2) > 0 );
            ib2.push_back( BoundaryElement( spob, b2_2, Real(0.5) * mal::Norm(a2_2 - b2_2) ) ); //Ref == 1
            // Check crossing of the new segment against the whole boundary1, which MAY CONTAIN a new segment if bAdvance1 == true!
            bCrossing2 = Intersection_Segment_IB( a2_2, b2_2, ib1, lambda2, lambda1, index_ibe1 );
#ifdef __ENABLE_IM_HIT_IP
            // Test IP and flood any hit one
            bHitIP2 = Test_Segment_IP( a2_2, b2_2, vec_ip, seed_ip_index,
                                       tr12, p_ds1, p_ds2,
                                       crossing_neighbourhood_radius,
                                       lambda_hit_ip2 );
#  ifdef __ENABLE_TRACE_IM
            if( bCrossing2 || bHitIP2 )
            {
                GEO_LOG_WARNING( "IM::Init() bHitIP = %d,%d and bCrossing = %d,%d at #steps = %d",
                                 bHitIP1, bHitIP2, bCrossing1, bCrossing2, num_steps_im );
            }
#  endif
#endif
        }
        else
            bCrossing2 = false;

        //-- Advance Intersection Map, if no crossing found
        if( !bCrossing1 && !bCrossing2 )
        {
            //-- Compute heuristic for all potential mappings

            /*\todo Consider different heuristic/distance measures
              (POB to POB, POB to boundariel, boundariel midpoints,
              whole boundariels...)

              \todo Should we consider the h_aa case
              (if NOT already added to IM, of course)??
              MAY it happen that it's the best case but
              never added because not checked?... I don't
              think so... but we MUST be completely sure
            */
#ifdef __ENABLE_IM_ADVANCE_CLOSEST
            Real h_ab( mal::NormSq( a1_2 - b2_2 ) );
            Real h_ba( mal::NormSq( b1_2 - a2_2 ) );
            Real h_bb( mal::NormSq( b1_2 - b2_2 ) );
#else //__ENABLE_IM_ADVANCE_PERPENDICULAR
            Vec2 e1( b1_2 - a1_2 );
            Vec2 e2( b2_2 - a2_2 );
            Real inv_norm_sq_e1( mal::Rcp( mal::NormSq( e1 ) ) );
            Real inv_norm_sq_e2( mal::Rcp( mal::NormSq( e2 ) ) );
            /*\todo MAYBE no need to normalize e1,e2, their relative weight could
               decrease proportionally to length, actually, because it represents
               a smaller boundariel thus angle is less important... */
            Vec2 ab( a1_2 - b2_2 );
            Vec2 ba( b1_2 - a2_2 );
            Vec2 bb( b1_2 - b2_2 );
            //\todo if normalization is required in ab,ba,bb, divide h_xy by NormSq instead!!
            Real h_ab( mal::Sq( mal::Dot( ab, e1 ) ) * inv_norm_sq_e1 + mal::Sq( mal::Dot( ab, e2 ) ) * inv_norm_sq_e2 );
            Real h_ba( mal::Sq( mal::Dot( ba, e1 ) ) * inv_norm_sq_e1 + mal::Sq( mal::Dot( ba, e2 ) ) * inv_norm_sq_e2 );
            Real h_bb( mal::Sq( mal::Dot( bb, e1 ) ) * inv_norm_sq_e1 + mal::Sq( mal::Dot( bb, e2 ) ) * inv_norm_sq_e2 );
#endif

            //-- Select best mapping wrt heuristic
#ifdef __ENABLE_IM_HIT_IP
            if( !bHitIP1 && !bHitIP2 )
            {
                if( h_bb <= h_ab && h_bb <= h_ba ) im.push_back( Pair( ib1.size()-1, ib2.size()-1, h_bb ) );
                else if( h_ab < h_ba ) im.push_back( Pair( ib1.size()-2, ib2.size()-1, h_ab ) );
                else im.push_back( Pair( ib1.size()-1, ib2.size()-2, h_ba ) );
            }
            else if( bHitIP1 && !bHitIP2 )
                im.push_back( Pair( ib1.size()-2, ib2.size()-1, h_ab ) );
            else if( !bHitIP1 && bHitIP2 )
                im.push_back( Pair( ib1.size()-1, ib2.size()-2, h_ba ) );
            else
            {
                //GEO_LOG_ASSERT(false, "Cannot advance!" );

                /*\todo This CAN happen if one or both IB hit an IP in
                  such a way that the last IB segment does now
                  generate crossing with the other IB. Specifically,
                  it may happen if ANY of the IB reach the IP with a
                  lambda > 1, as may happen due to thickness in
                  FastTest_ForwardSegment2_Point2

                  \note The actual IP hit by each IB MAY be
                  different... unless we guarantee that there will
                  only be 1 IP within any thickness neighbourhood.

                  If BOTH bHitIP, there SHOULD be a crossing near the
                  hit IP... we could try to find it with a Line2Line
                  test instead of a Segment2Segment one...
                */
                GEO_LOG_WARNING("Cannot advance, forcing bCrossing1 && bCrossing2!");
                bCrossing1 = true;
                bCrossing2 = true;
                index_ibe1 = ib1.size()-2;
                lambda1 = 1.0;
                index_ibe2 = ib2.size()-2;
                lambda2 = 1.0;
            }
#else
            if( h_bb <= h_ab && h_bb <= h_ba ) //prioritize advancing both IB1 and IB2
                im.push_back( Pair( ib1.size()-1, ib2.size()-1, h_bb ) );
            else if( h_ab < h_ba ) //otherwise, either ab or ba is < bb
                im.push_back( Pair( ib1.size()-2, ib2.size()-1, h_ab ) );
            else // h_ba <= h_ab
                im.push_back( Pair( ib1.size()-1, ib2.size()-2, h_ba ) );
#endif //__ENABLE_IM_HIT_IP

            /*\todo detect backtracking, retreat IM as
              needed... backtracking would require access
              to ib.size()-3 POB, may not exist
            */
        }
        num_steps_im++;
    }
    while( !bCrossing1 && !bCrossing2 && num_steps_im < cMaxStepsIM );
    GEO_LOG_ASSERT( num_steps_im < cMaxStepsIM, "IM::Init() Too many steps > %d", cMaxStepsIM );

    /*
      After crossing, we may need to undo part of
      the mapping and unref unused POB, as the crossing
      may happen at boundariels in IB1/IB2 that are NOT
      the last added ones.

      \todo should compute crossing point, maybe
      add it to mapping, undo some mapping,
      discard potential IP near it, etc... by now
      we just stop
    */
    if( bCrossing1 )
    {
        //Compute crossing geometry
        Vec2 a1_2( ib1[ib1.size()-2].m_Pos_2 );
        Vec2 b1_2( ib1.back().m_Pos_2 );
        Vec2 cp_2( a1_2 + lambda1 * (b1_2 - a1_2) );

        //\todo Refine last entry in both IB to find the actual POB at the crossing point
        //\tood SUBSTITUTE last POB in IB with actual crossing POB
        //\todo Refine last mapping to include crossing POB
        //\todo Consume potentially unconsumed IP near crossing POB

#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
        if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
        {
            Vec2 a2_2( ib2[index_ibe2].m_Pos_2 );
            Vec2 b2_2( ib2[index_ibe2+1].m_Pos_2 );
            vd.m_vecIM_XSegment.push_back( std::make_pair( tr2 * a1_2, tr2 * b1_2 ) );
            vd.m_vecIM_XSegment.push_back( std::make_pair( tr2 * a2_2, tr2 * b2_2 ) );
            //vd.m_vecIM_XPoint.push_back( tr2 * cp_2 );
        }
#endif

#ifdef __ENABLE_IM_BACKTRACK_AFTER_CROSSING
#  ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
        if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
            for( unsigned int it_ib2=index_ibe2+1; it_ib2<ib2.size(); it_ib2++ )
                vd.m_vecIM_IB_Discarded.push_back( tr2 * ib2[it_ib2].m_Pos_2 );
#  endif
        //WIP: Backtrack ib2 and mapping up to crossing index_ibe2
        for( unsigned int it_ib2=index_ibe2+1; it_ib2<ib2.size(); it_ib2++ )
            p_ds2->POB_DecRef( ib2[it_ib2].m_POB );
        ib2.resize( index_ibe2+1 );
        unsigned int last_im( im.size()-1 );
        /*Remove mappings after crossing
          TEMP: This leaves orphan POB on IB1
          while( im[last_im].second > index_ibe2 ) last_im--;
          im.resize( last_im+1 );
        */
        // Remap all POB1 mapped to POB2 beyond crossing to crossing POB2
        while( im[last_im].m_Id2 > index_ibe2 ) im[last_im--].m_Id2 = index_ibe2;
#endif
    }
    else // bCrossing2
    {
        //Compute crossing geometry
        Vec2 a2_2( ib2[ib2.size()-2].m_Pos_2 );
        Vec2 b2_2( ib2.back().m_Pos_2 );
        Vec2 cp_2( a2_2 + lambda2 * (b2_2 - a2_2) );
        //\todo Backtrack ib1 and mapping up to crossing index_ibe1

#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
        if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
        {
            Vec2 a1_2( ib1[index_ibe1].m_Pos_2 );
            Vec2 b1_2( ib1[index_ibe1+1].m_Pos_2 );
            vd.m_vecIM_XSegment.push_back( std::make_pair( tr2 * a1_2, tr2 * b1_2 ) );
            vd.m_vecIM_XSegment.push_back( std::make_pair( tr2 * a2_2, tr2 * b2_2 ) );
            //vd.m_vecIM_XPoint.push_back( tr2 * cp_2 );
        }
#endif

#ifdef __ENABLE_IM_BACKTRACK_AFTER_CROSSING
#  ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
        if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
            for( unsigned int it_ib1=index_ibe1+1; it_ib1<ib1.size(); it_ib1++ )
                vd.m_vecIM_IB_Discarded.push_back( tr2 * ib1[it_ib1].m_Pos_2 );
#  endif
        //WIP: Backtrack ib1 and mapping up to crossing index_ibe1
        for( unsigned int it_ib1=index_ibe1+1; it_ib1<ib1.size(); it_ib1++ )
            p_ds1->POB_DecRef( ib1[it_ib1].m_POB );
        ib1.resize( index_ibe1+1 );
        unsigned int last_im( im.size()-1 );
        /*Remove mappings after crossing
          TEMP: This leaves orphan POB on IB2
          while( im[last_im].first > index_ibe1 ) last_im--;
          im.resize( last_im+1 );
        */
        // Remap all POB2 mapped to POB1 beyond crossing to crossing POB1
        while( im[last_im].m_Id1 > index_ibe1 ) im[last_im--].m_Id1 = index_ibe1;
#endif
    }

    // Viz Raw Mapping
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
        for( unsigned int it_im=0; it_im<im.size(); it_im++ )
            vd.m_vecIM_RMap.push_back( std::make_pair( tr2 * ib1[ im[it_im].m_Id1 ].m_Pos_2,
                                                       tr2 * ib2[ im[it_im].m_Id2 ].m_Pos_2 ) );
#endif

    //---- Post-processes IM (\note Application order IS relevant)
    if( im_flags.Test( Context::eIMF_Reverse ) )
    {
        ReverseMapping();
    }
    if( im_flags.Test( Context::eIMF_Reduce ) )
    {
        ReduceMapping();
    }
    if( im_flags.Test( Context::eIMF_MinDist ) )
    {
        ComputeNormals( tr12, p_ds1, p_ds2 );
        MinHeuristicPairRemapping( Heuristic_DistSq );
    }
    if( im_flags.Test( Context::eIMF_MinDistAngle ) )
    {
        ComputeNormals( tr12, p_ds1, p_ds2 );
        MinHeuristicPairRemapping( Heuristic_DistSq_AngleSq );
    }
    if( im_flags.Test( Context::eIMF_Project ) )
    {
        ComputeNormals( tr12, p_ds1, p_ds2 );
        ProjectMapping( tr12, tr12.Inverse(), p_ds1, p_ds2, step_length );
    }
    if( im_flags.Test( Context::eIMF_Fix ) )
    {
        ComputeNormals( tr12, p_ds1, p_ds2 );
        FixMapping();
    }

#ifdef __ENABLE_IM_CD_AREA_SCALE_DEPTH
    /* Modify BE radius to radius/valence, so that each the BE1
       appearance in the IM has a proper radius, distributed across
       all interacting BE2. This guarantees that each IB1 and IB2 area
       will only interact *once* with the other boundary, and that the
       total force integration area of IB1 and IB2 will be "correct".
    */
    ComputeValences();
    for( unsigned int it_im=0; it_im<m_IM.size(); it_im++ )
    {
        if( m_IB1[ m_IM[it_im].m_Id1 ].m_Valence > 0 ) m_IB1[ m_IM[it_im].m_Id1 ].m_Radius /= m_IB1[ m_IM[it_im].m_Id1 ].m_Valence;
        else m_IB1[ m_IM[it_im].m_Id1 ].m_Radius = 0;
        if( m_IB2[ m_IM[it_im].m_Id2 ].m_Valence > 0 ) m_IB2[ m_IM[it_im].m_Id2 ].m_Radius /= m_IB2[ m_IM[it_im].m_Id2 ].m_Valence;
        else m_IB2[ m_IM[it_im].m_Id2 ].m_Radius = 0;
    }
#endif

    //---- VIZ final IM
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
    ComputeNormals( tr12, p_ds1, p_ds2 ); //\todo Required to Viz
    //ComputeValences();  \todo Already done
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
    {
        Real length1(0);
        for( unsigned int it_ibe=0; it_ibe < ib1.size(); it_ibe++ ) length1 += ib1[it_ibe].m_Radius * ib1[it_ibe].m_Valence;
        Real length2(0);
        for( unsigned int it_ibe=0; it_ibe < ib2.size(); it_ibe++ ) length2 += ib2[it_ibe].m_Radius * ib2[it_ibe].m_Valence;
        vd.m_vecIM_XLength.push_back( std::make_pair( length1, length2 ) );
    }
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
        for( unsigned int it_im=0; it_im<im.size(); it_im++ )
            vd.m_vecIM_Map.push_back( std::make_pair( tr2 * ib1[ im[it_im].m_Id1 ].m_Pos_2,
                                                      tr2 * ib2[ im[it_im].m_Id2 ].m_Pos_2 ) );
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
        for( unsigned int it_ibe=0; it_ibe < ib1.size(); it_ibe++ )
        {
            vd.m_vecIM_IB1_Points.push_back( tr2 * ib1[it_ibe].m_Pos_2 );
            vd.m_vecIM_IB1_Normals.push_back( tr2.m_Rot * ib1[it_ibe].m_Normal_2 );
            vd.m_vecIM_IB1_Radii.push_back( ib1[it_ibe].m_Radius );
        }
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
        for( unsigned int it_ibe=0; it_ibe < ib2.size(); it_ibe++ )
        {
            vd.m_vecIM_IB2_Points.push_back( tr2 * ib2[it_ibe].m_Pos_2 );
            vd.m_vecIM_IB2_Normals.push_back( tr2.m_Rot * ib2[it_ibe].m_Normal_2 );
            vd.m_vecIM_IB2_Radii.push_back( ib2[it_ibe].m_Radius );
        }
#endif
    /*TEMP: if not stored, DecRef ALL POB in the IM!!
    for( unsigned int it_ibe=0; it_ibe < ib1.size(); it_ibe++ ) p_ds1->POB_DecRef( ib1[it_ibe].m_POB );
    for( unsigned int it_ibe=0; it_ibe < ib2.size(); it_ibe++ ) p_ds2->POB_DecRef( ib2[it_ibe].m_POB );
    */

    return true;
}

bool IntersectionMapping::Update( std::vector<StochasticSPC2::IntersectionPair> &vec_ip, unsigned int seed_ip_index,
                                  IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 )
{
    GEO_ASSERT( !vec_ip[seed_ip_index].IsFlooded() );

#ifdef __ENABLE_TRACE_IM
    GEO_LOG( "IM::Update IP[%d] = (%d,%d)", seed_ip_index, vec_ip[seed_ip_index].m_POB1, vec_ip[seed_ip_index].m_POB2 );
#endif

    //Update and flood IP as necessary
    vec_ip[seed_ip_index].SetFlooded(true);
    return false;
}

void IntersectionMapping::Destroy( IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 )
{
#ifdef __ENABLE_TRACE_IM
    GEO_LOG( "IM::Destroy IP = (%d,%d)", m_IB1[0].m_POB, m_IB2[0].m_POB );
#endif
    for( unsigned int it_ibe=0; it_ibe < m_IB1.size(); it_ibe++ ) p_ds1->POB_DecRef( m_IB1[it_ibe].m_POB );
    for( unsigned int it_ibe=0; it_ibe < m_IB2.size(); it_ibe++ ) p_ds2->POB_DecRef( m_IB2[it_ibe].m_POB );
    Clear();
}

bool IntersectionMapping::ComputeContactData( ContactData2 &cd, const Transform2 &tr2,
                                              IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
                                              Real contact_depth_epsilon ) const
{
#ifdef __ENABLE_IM_CD
    for( unsigned int it_im=0; it_im<m_IM.size(); it_im++ )
    {
        Vec2 p1_0( tr2 * m_IB1[ m_IM[it_im].m_Id1 ].m_Pos_2 );
        Vec2 p2_0( tr2 * m_IB2[ m_IM[it_im].m_Id2 ].m_Pos_2 );
        Vec2 n12_0( p2_0 - p1_0 );
        Real depth( mal::Norm(n12_0) );
        if( depth > contact_depth_epsilon )
        {
#ifdef __ENABLE_IM_CD_AREA_SCALE_DEPTH
            // CP area12 should combine r1 and r2 someway... by now we just average them...
            Real area12( Real(0.5) * ( m_IB1[ m_IM[it_im].m_Id1 ].m_Radius + m_IB2[ m_IM[it_im].m_Id2 ].m_Radius ) );
            cd.AddCP( p1_0, p2_0, n12_0 / depth, depth, area12 );
#else
            cd.AddCP( p1_0, p2_0, n12_0 / depth, depth );

#endif
            cd.AddPOF( p_ds1->POB_PointOnFeature( m_IB1[ m_IM[it_im].m_Id1 ].m_POB ),
                       p_ds2->POB_PointOnFeature( m_IB2[ m_IM[it_im].m_Id2 ].m_POB ) );
        }
        else
        {
            /*\todo Coincident contact should still yield a VALID
              normal, even with depth=0, to be used with
              constraint-based contact response (not with penalty.
              Normal could be POB1 or POB2 or a mix of both.
            */
            cd.AddCP( p1_0, p2_0, tr2.m_Rot * p_ds2->POB_Normal( m_IB2[ m_IM[it_im].m_Id2 ].m_POB ), 0 ); //\todo Using n2 by now
            cd.AddPOF( p_ds1->POB_PointOnFeature( m_IB1[ m_IM[it_im].m_Id1 ].m_POB ),
                       p_ds2->POB_PointOnFeature( m_IB2[ m_IM[it_im].m_Id2 ].m_POB ) );
        }
    }
#endif //__ENABLE_IM_CD
    return m_IM.size() > 0;
}

void IntersectionMapping::ComputeNormals( const Transform2 &tr12, IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 )
{
    for( unsigned int it_ib1=0; it_ib1<m_IB1.size(); it_ib1++ )
        m_IB1[it_ib1].m_Normal_2 = tr12.m_Rot * p_ds1->POB_Normal( m_IB1[it_ib1].m_POB );
    for( unsigned int it_ib2=0; it_ib2<m_IB2.size(); it_ib2++ )
        m_IB2[it_ib2].m_Normal_2 = p_ds2->POB_Normal( m_IB2[it_ib2].m_POB );
}

void IntersectionMapping::ComputeValences()
{
    for( unsigned int it_ib1=0; it_ib1<m_IB1.size(); it_ib1++ ) m_IB1[it_ib1].m_Valence = 0;
    for( unsigned int it_ib2=0; it_ib2<m_IB2.size(); it_ib2++ ) m_IB2[it_ib2].m_Valence = 0;
    for( unsigned int it_im=0; it_im<m_IM.size(); it_im++ )
    {
        m_IB1[ m_IM[it_im].m_Id1 ].m_Valence++;
        m_IB2[ m_IM[it_im].m_Id2 ].m_Valence++;
    }
}

void IntersectionMapping::ReverseMapping()
{
    //TEMP: GEO_ASSERT( ib1.size() > 1 && ib2.size() > 1 ); //BUG THIS
    //RAISES..... but makes no sense because on init we ALWAYS add
    //2POB to each IB... unless we REMOVE too many POB at
    //__ENABLE_IM_BACKTRACK_AFTER_CROSSING ???
    if( m_IB1.size() < 2 || m_IB2.size() < 2 ) return;
    m_IM.clear();

    /*TEMP: This doesn't work very well... bad mapping
      around crossing point, it's just to compare with
      direct mapping, but cannot replace it.

      \note Notice that "wrong" direct mappings tend to be
      corrected in reverse mappings, and viceversa.
    */

    int it_ib1( m_IB1.size()-2 );
    int it_ib2( m_IB2.size()-2 );
    // Init im with last POB
    m_IM.push_back( Pair( m_IB1.size()-1, m_IB2.size()-1, 0 ) );
    while( it_ib1 > 0 || it_ib2 > 0 )
    {
        // Get lookahead boundariel geometry => segments (a1,b1) and (a2,b2)
        Vec2 a1_2( m_IB1[ it_ib1+1 ].m_Pos_2 );
        Vec2 b1_2( m_IB1[ it_ib1   ].m_Pos_2 );
        Vec2 a2_2( m_IB2[ it_ib2+1 ].m_Pos_2 );
        Vec2 b2_2( m_IB2[ it_ib2   ].m_Pos_2 );

#ifdef __ENABLE_IM_ADVANCE_CLOSEST
        Real h_ab( mal::NormSq( a1_2 - b2_2 ) );
        Real h_ba( mal::NormSq( b1_2 - a2_2 ) );
        Real h_bb( mal::NormSq( b1_2 - b2_2 ) );
#else //__ENABLE_IM_ADVANCE_PERPENDICULAR
        Vec2 e1( b1_2 - a1_2 );
        Vec2 e2( b2_2 - a2_2 );
        Real inv_norm_sq_e1( mal::Rcp( mal::NormSq( e1 ) ) );
        Real inv_norm_sq_e2( mal::Rcp( mal::NormSq( e2 ) ) );
        /*\todo MAYBE no need to normalize e1,e2, their relative weight could
          decrease proportionally to length, actually, because it represents
          a smaller boundariel thus angle is less important... */
        Vec2 ab( a1_2 - b2_2 );
        Vec2 ba( b1_2 - a2_2 );
        Vec2 bb( b1_2 - b2_2 );
        //\todo if normalization is required in ab,ba,bb, divide h_xy by NormSq instead!!
        Real h_ab( mal::Sq( mal::Dot( ab, e1 ) ) * inv_norm_sq_e1 + mal::Sq( mal::Dot( ab, e2 ) ) * inv_norm_sq_e2 );
        Real h_ba( mal::Sq( mal::Dot( ba, e1 ) ) * inv_norm_sq_e1 + mal::Sq( mal::Dot( ba, e2 ) ) * inv_norm_sq_e2 );
        Real h_bb( mal::Sq( mal::Dot( bb, e1 ) ) * inv_norm_sq_e1 + mal::Sq( mal::Dot( bb, e2 ) ) * inv_norm_sq_e2 );
#endif

        if( h_bb <= h_ab && h_bb <= h_ba ) //prioritize advancing both IB1 and IB2
            m_IM.push_back( Pair( it_ib1--, it_ib2--, h_bb ) );
        else
        {
            if( it_ib2 > 0 && it_ib1 > 0 )
            {
                if( h_ab < h_ba ) m_IM.push_back( Pair( it_ib1+1, it_ib2--, h_ab ) );
                else m_IM.push_back( Pair( it_ib1--, it_ib2+1, h_ba ) );
            }
            else
            {
                if( it_ib2 > 0 ) m_IM.push_back( Pair( it_ib1+1, it_ib2--, h_ab ) );
                else m_IM.push_back( Pair( it_ib1--, it_ib2+1, h_ba ) );
            }
        }
        // Clamp to 0... simplifies previous conditional
        if( it_ib1 < 0 ) it_ib1 = 0;
        if( it_ib2 < 0 ) it_ib2 = 0;
    }
}

/*! Remove redundant mappings
  \note First/Last mappings are never considered redundant.
*/
void IntersectionMapping::ReduceMapping()
{
    //GEO_ASSERT( m_IM.size() > 1 );
    if( m_IM.size() < 2 ) return;
    // Copy non-redundant mappings to reduced IM
    std::vector< Pair > im;
    im.push_back( m_IM[0] );
    for( unsigned int it_im=1; it_im<m_IM.size()-1; it_im++ )
    {
        const Pair &prev( im.back() );
        const Pair &pair( m_IM[it_im] );
        const Pair &next( m_IM[it_im+1] );
        bool bRedundant( prev.m_Id1+1 == next.m_Id1 && prev.m_Id2+1 == next.m_Id2 );
        if( !bRedundant ) im.push_back( m_IM[it_im] );
    }
#ifdef __ENABLE_TRACE_IM
    if( im.size() < m_IM.size() ) { GEO_LOG("IM::ReduceMapping() reduced %d < %d", (int)im.size(), (int)m_IM.size() ); }
#endif
    // Replace IM with reduced IM
    std::swap( m_IM, im );
}

/*! Remove degenerated mappings
  Degenerated if:
  - v21 * n1 < 0
  - v21 * n2 > 0
  \todo Could accept some margin around <0 and >0 if necessary
*/
void IntersectionMapping::FixMapping()
{
    //GEO_ASSERT( m_IM.size() > 1 );
    if( m_IM.size() < 2 ) return;
    // Copy non-redundant mappings to reduced IM
    std::vector< Pair > im;
    im.push_back( m_IM[0] );
    for( unsigned int it_im=0; it_im<m_IM.size(); it_im++ )
    {
        const Pair &pair( m_IM[it_im] );
        Vec2 v21_2( m_IB1[pair.m_Id1].m_Pos_2 - m_IB2[pair.m_Id2].m_Pos_2 );
        Vec2 n1_2( m_IB1[pair.m_Id1].m_Normal_2 );
        Vec2 n2_2( m_IB2[pair.m_Id2].m_Normal_2 );
        bool bDegenerated( mal::Dot( v21_2, n1_2 ) < 0 || mal::Dot( v21_2, n2_2 ) > 0 );
        if( !bDegenerated ) im.push_back( m_IM[it_im] );
    }
#ifdef __ENABLE_TRACE_IM
    if( im.size() < m_IM.size() ) { GEO_LOG("IM:FixMapping() removed %d mappings", (int)m_IM.size() - (int)im.size() ); }
#endif
    // Replace IM with reduced IM
    std::swap( m_IM, im );
}

/* Project POB on the other IB
//\todo Use IDS::ClosestPOB() to project POB on the other IB
*/
void IntersectionMapping::ProjectMapping( const Transform2 &tr12, const Transform2 &tr21,
                                          IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2,
                                          Real neighbourhood_radius )
{
    if( m_IM.size() < 2 ) return;
    std::vector< Pair > im;
    std::vector< BoundaryElement > ib1;
    std::vector< BoundaryElement > ib2;

#define __ENABLE_PROJECT_MAPPING_ALWAYS_PROJECT_BOTH
#ifdef __ENABLE_PROJECT_MAPPING_ALWAYS_PROJECT_BOTH

#ifdef __ENABLE_IM_CD_AREA_SCALE_DEPTH
    ComputeValences();
    /* \todo Radius rescaling may not be completely correct. For new
       ClosestPOB() generated BE, it seems so, but for original BE
       that REMAIN in the IM, they preserve the original radius, while
       thei should also be split... not sure, I think I should track
       GLOBAL length and ensure it's preserved (approx) for all IM
       postprocesses.
    */
#endif

    // Copy IB
    ib1 = m_IB1;
    ib2 = m_IB2;
    // Init first mapping
    im.push_back( m_IM[0] );
    // Project and add all IB and IM unconditionally
    for( unsigned int it_im=1; it_im<m_IM.size(); it_im++ )
    {
        const Pair &pair( m_IM[it_im] );
        point_on_boundary_id_type closest_pob1 = p_ds1->ClosestPOB( m_IB1[ pair.m_Id1 ].m_POB,
                                                                    tr21 * m_IB2[ pair.m_Id2 ].m_Pos_2,
                                                                    neighbourhood_radius, 0 );
        point_on_boundary_id_type closest_pob2 = p_ds2->ClosestPOB( m_IB2[ pair.m_Id2 ].m_POB,
                                                                    m_IB1[ pair.m_Id1 ].m_Pos_2,
                                                                    neighbourhood_radius, 0 );
        if( closest_pob1 == m_IB1[ pair.m_Id1 ].m_POB ) p_ds1->POB_IncRef( closest_pob1 );
        if( closest_pob2 == m_IB2[ pair.m_Id2 ].m_POB ) p_ds2->POB_IncRef( closest_pob2 );
        // Add closest POB to IB
        ib1.push_back( BoundaryElement( closest_pob1, tr12 * p_ds1->POB_Position(closest_pob1),
#ifdef __ENABLE_IM_CD_AREA_SCALE_DEPTH //As we split each original mapping into 2, their area is also split in 2 parts
                                        Real(0.5) * m_IB1[pair.m_Id1].m_Radius / m_IB1[pair.m_Id1].m_Valence ) );
#else
                                        0 ) );
#endif
        ib2.push_back( BoundaryElement( closest_pob2, p_ds2->POB_Position(closest_pob2),
#ifdef __ENABLE_IM_CD_AREA_SCALE_DEPTH //As we split each original mapping into 2, their area is also split in 2 parts
                                        Real(0.5) * m_IB2[pair.m_Id2].m_Radius / m_IB2[pair.m_Id2].m_Valence ) );
#else
                                        0 ) );
#endif

#ifdef __ENABLE_IM_CD_AREA_SCALE_DEPTH //Remove from original BE the radius we've assigned to new projected BE
        ib1[ pair.m_Id1 ].m_Radius -= ib1.back().m_Radius;
        ib2[ pair.m_Id2 ].m_Radius -= ib2.back().m_Radius;
#endif
        // Add mappings with closest POB
        im.push_back( Pair( pair.m_Id1, ib2.size()-1, 0 ) );
        im.push_back( Pair( ib1.size()-1, pair.m_Id2, 0 ) );
        /* Copy original pair
        im.push_back( pair );
        */
    }

#ifdef __ENABLE_TRACE_IM
    GEO_LOG("IM:ProjectMapping() (ALWAYS) Mapping (#im,#ib1,#ib2) = (%d,%d,%d) projected to (%d,%d,%d)",
            (int)m_IM.size(), (int)m_IB1.size(), (int)m_IB2.size(),
            (int)im.size(), (int)ib1.size(), (int)ib2.size() );
#endif
    // Replace
    std::swap( m_IM, im );
    std::swap( m_IB1, ib1 );
    std::swap( m_IB2, ib2 );

#else

    // Init first mapping
    im.push_back( m_IM[0] );
    ib1.push_back( m_IB1[0] ); p_ds1->POB_IncRef( m_IB1[0].m_POB );
    ib2.push_back( m_IB2[0] ); p_ds2->POB_IncRef( m_IB2[0].m_POB );

    /* for each mapping in IM, check if "perpendicular enough".  If
       true, advance (copy IM to im, add IB to ib) Else, project one
       or both bdrl onto the other IB, insert new POB into IB, and
       adjust offset to keep remaining mappings valid.
    */
    unsigned int num_added1(0);
    unsigned int num_added2(0);
    const Real cMaxProjectionAngleThreshold( mal::Deg2Rad(10) );
    const Real cMinProjectionCosSqThreshold( mal::Sq( mal::Cos( cMaxProjectionAngleThreshold ) ) );
    for( unsigned int it_im=1; it_im<m_IM.size()-1; it_im++ )
    {
        const Pair &pair( m_IM[it_im] );
        Vec2 v( m_IB1[pair.m_Id1].m_Pos_2 - m_IB2[pair.m_Id2].m_Pos_2 );
        Real inv_norm_sq( mal::Rcp( mal::NormSq( v ) ) );
        Real cos_sq_wrt_ib1( inv_norm_sq * mal::Sq( mal::Dot( v, m_IB1[pair.m_Id1].m_Normal_2 ) ) );
        Real cos_sq_wrt_ib2( inv_norm_sq * mal::Sq( mal::Dot( v, m_IB2[pair.m_Id2].m_Normal_2 ) ) );
        bool bProjectOnIB1( cos_sq_wrt_ib1 < cMinProjectionCosSqThreshold );
        bool bProjectOnIB2( cos_sq_wrt_ib2 < cMinProjectionCosSqThreshold );

        //\todo CONSIDER projecting ONLY THE WORST-angle POB on
        //the best-angle boundariel, to ensure at least one
        //ortogonality, instead of projecting both... so that !(bProjectOnIB1 && bProjectOnIB2)

        //--  Project required POB and re-check if the projection returns different POB on the other IB
        point_on_boundary_id_type closest_pob1( cInvalidPOB );
        // Find POB1 on IB1 closest to POB2, using current POB1 as an approximation
        if( bProjectOnIB1 ) closest_pob1 = p_ds1->ClosestPOB( m_IB1[ pair.m_Id1 ].m_POB,
                                                              tr21 * m_IB2[ pair.m_Id2 ].m_Pos_2,
                                                              neighbourhood_radius, 0 );
        bProjectOnIB1 = bProjectOnIB1 && closest_pob1 != m_IB1[ pair.m_Id1 ].m_POB;
        // Find POB2 on IB2 closest to POB1, using current POB2 as an approximation
        point_on_boundary_id_type closest_pob2( cInvalidPOB );
        if( bProjectOnIB2 ) closest_pob2 =  p_ds2->ClosestPOB( m_IB2[ pair.m_Id2 ].m_POB,
                                                               m_IB1[ pair.m_Id1 ].m_Pos_2,
                                                               neighbourhood_radius, 0 );
        bProjectOnIB2 = bProjectOnIB2 && closest_pob2 != m_IB2[ pair.m_Id2 ].m_POB;

        //-- If not both IB need projection, copy current IM pair
        if( !bProjectOnIB1 || !bProjectOnIB2 )
        {
            // Copy IM
            im.push_back( pair );
            // Reindex IB
            im.back().m_Id1 += num_added1;
            im.back().m_Id2 += num_added2;
            // If new IM references IB not yet added, add them now
            if( im.back().m_Id1 == ib1.size() ) ib1.push_back( m_IB1[ pair.m_Id1 ] ); COMPUTE RADIUS
            if( im.back().m_Id2 == ib2.size() ) ib2.push_back( m_IB2[ pair.m_Id2 ] ); COMPUTE RADIUS
            im.back().m_H = 0; //\todo COMPUTE HEURISTIC
        }

        //\todo CONSIDER projecting ONLY THE WORST-angle POB on
        //the best-angle boundariel, to ensure at least one
        //ortogonality, instead of projecting both...
        if( bProjectOnIB1 )
        {
            // Add new closest POB2 to IB2 and add pair POB1,CPOB2
            ib1.push_back( BoundaryElement( closest_pob1, tr12 * p_ds1->POB_Position(closest_pob1), COMPUTE RADIUS ) );
            num_added1++;
            im.push_back( Pair( ib1.size()-1, pair.m_Id2 + num_added2, 0 ) ); //\todo Compute heuristic!
            // If new IM references IB not yet added, add them now
            if( im.back().m_Id2 == ib2.size() ) ib2.push_back( m_IB2[ pair.m_Id2 ] ); COMPUTE RADIUS
        }
        //else //Existing pair already added, do nothing

        if( bProjectOnIB2 )
        {
            // Add new closest POB2 to IB2 and add pair POB1,CPOB2
            ib2.push_back( BoundaryElement( closest_pob2, p_ds2->POB_Position(closest_pob2), COMPUTE RADIUS ) );
            num_added2++;
            im.push_back( Pair( pair.m_Id1 + num_added1, ib2.size()-1, 0 ) ); //\todo Compute heuristic!
            // If new IM references IB not yet added, add them now
            if( im.back().m_Id1 == ib1.size() ) ib1.push_back( m_IB1[ pair.m_Id1 ] ); COMPUTE RADIUS
        }
        //else //Existing pair already added, do nothing

        //\todo MUST decref all discarded (uncopied) POB on any IB!!!
        if( bProjectOnIB1 && bProjectOnIB2 )
        {
            p_ds1->POB_DecRef( m_IB1[pair.m_Id1].m_POB );
            p_ds2->POB_DecRef( m_IB2[pair.m_Id2].m_POB );
        }
    }
    //\todo ADD OR FIX LAST MAPPING, WHICH TENDS TO BE wrong because it does not use exact crossing points.

#ifdef __ENABLE_TRACE_IM
    if( num_added1 + num_added2 > 0 )
    {
        GEO_LOG("IM:ProjectMapping() Mapping (#im,#ib1,#ib2) = (%d,%d,%d) projected to (%d,%d,%d)",
                (int)m_IM.size(), (int)m_IB1.size(), (int)m_IB2.size(),
                (int)im.size(), (int)ib1.size(), (int)ib2.size() );
    }
#endif

    // Replace
    std::swap( m_IM, im );
    std::swap( m_IB1, ib1 );
    std::swap( m_IB2, ib2 );

#endif //__ENABLE_PROJECT_MAPPING_ALWAYS_PROJECT_BOTH
}

/* Remapping according to min heuristic exclusively.
   \note It's slow, O(n^2), just for comparison
   \note Non-sequential/crossing mappings generated
   \note No duplicated mappings are generated
*/
void IntersectionMapping::MinHeuristicPairRemapping( HeuristicFunc HF )
{
    // Clear mapping and valences
    m_IM.clear();
    for( unsigned int it_ib1=0; it_ib1<m_IB1.size(); it_ib1++ ) m_IB1[it_ib1].m_Valence = 0;
    for( unsigned int it_ib2=0; it_ib2<m_IB2.size(); it_ib2++ ) m_IB2[it_ib2].m_Valence = 0;
    // Map IB1 to IB2
    for( unsigned int it_ib1=0; it_ib1<m_IB1.size(); it_ib1++ )
    {
        unsigned int best_ib2(0);
        Real min_h( HF( m_IB1[it_ib1].m_Pos_2, m_IB1[it_ib1].m_Normal_2, m_IB1[it_ib1].m_Radius,
                        m_IB2[0].m_Pos_2, m_IB2[0].m_Normal_2, m_IB2[0].m_Radius ) );
        for( unsigned int it_ib2=0; it_ib2<m_IB2.size(); it_ib2++ )
        {
            Real h( HF( m_IB1[it_ib1].m_Pos_2, m_IB1[it_ib1].m_Normal_2, m_IB1[it_ib1].m_Radius,
                        m_IB2[it_ib2].m_Pos_2, m_IB2[it_ib2].m_Normal_2, m_IB2[it_ib2].m_Radius ) );
            if( h < min_h )
            {
                best_ib2 = it_ib2;
                min_h = h;
            }
        }
        m_IM.push_back( Pair( it_ib1, best_ib2, min_h ) );
        m_IB1[it_ib1].m_Valence = 1;
        GEO_ASSERT( m_IB2[best_ib2].m_Valence < 255 );
        m_IB2[best_ib2].m_Valence++;
    }
    // Map unmapped IB2 to IB1
    for( unsigned int it_ib2=0; it_ib2<m_IB2.size(); it_ib2++ )
    {
        if( 0 == m_IB2[it_ib2].m_Valence )
        {
            unsigned int best_ib1(0);
            Real min_h( HF( m_IB1[0].m_Pos_2, m_IB1[0].m_Normal_2, m_IB1[0].m_Radius,
                            m_IB2[it_ib2].m_Pos_2, m_IB2[it_ib2].m_Normal_2, m_IB2[it_ib2].m_Radius ) );
            for( unsigned int it_ib1=0; it_ib1<m_IB1.size(); it_ib1++ )
            {
                Real h( HF( m_IB1[it_ib1].m_Pos_2, m_IB1[it_ib1].m_Normal_2, m_IB1[it_ib1].m_Radius,
                            m_IB2[it_ib2].m_Pos_2, m_IB2[it_ib2].m_Normal_2, m_IB2[it_ib2].m_Radius ) );
                if( h < min_h )
                {
                    best_ib1 = it_ib1;
                    min_h = h;
                }
            }
            m_IM.push_back( Pair( best_ib1, it_ib2, min_h ) );
            GEO_ASSERT( m_IB1[best_ib1].m_Valence < 255 && m_IB2[it_ib2].m_Valence < 255 );
            m_IB1[best_ib1].m_Valence++;
            m_IB2[it_ib2].m_Valence++;
        }
    }
}

// Compute first intersection between Segment2 and IB2
bool IntersectionMapping::Intersection_Segment_IB( const Vec2 &a1, const Vec2 &b1,
                                                   const std::vector<BoundaryElement> &vec_ibe2,
                                                   //\todo Real epsilon_length,
                                                   Real &lambda1, Real &lambda2, unsigned int &index_ibe2 )
{
    //\todo We MUST handle coincidences, as incremental mapping termination REQUIRES finding a crossing that is known to EXIST
    for( unsigned int it_ibe=0; it_ibe<vec_ibe2.size()-1; it_ibe++ )
        if( Intersection_Segment2_Segment2( a1, b1,
                                            vec_ibe2[it_ibe].m_Pos_2, vec_ibe2[it_ibe+1].m_Pos_2,
                                            //\todo epsilon_length,
                                            lambda1, lambda2 ) )
        {
            index_ibe2 = it_ibe;
            return true;
        }
    return false;
}

#endif //__ENABLE_IM_USE_BOUNDARIELS
#endif //__ENABLE_IM


//----------------------------------------------------------------
//-- TEMP: IP implementation that depends on IntersectionMapping...
//----------------------------------------------------------------
bool StochasticSPC2::IntersectionPair::MergeIM( StochasticSPC2::IntersectionPair &ip,
                                                IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 )
{
    /*\todo Merge is possible when at most one pair in
      (unique_pair,pair) has an IM (just copy) and when
      both have IM and these are simple enough to merge.
    */
    if( HasIM() && ip.HasIM() ) { /*\todo return TryToMerge;*/ return false; }
    else if( HasIM() && !ip.HasIM() ) return true;
    else if( !HasIM() && ip.HasIM() ) { m_pIM = ip.m_pIM; ip.m_pIM = 0; return true; }
    else return true;
}

void StochasticSPC2::IntersectionPair::DiscardIM( IDomainSampler2 *p_ds1, IDomainSampler2 *p_ds2 )
{
    m_bFlooded = false;
#ifdef __ENABLE_IM
    if( 0 != m_pIM )
    {
        m_pIM->Destroy( p_ds1, p_ds2 );
        delete m_pIM;
        m_pIM = 0;
    }
#endif
}

//----------------------------------------------------------------
//-- Pair-processing functions
//----------------------------------------------------------------

#ifdef __USE_CACHED_POB_DATA
void UpdatePOB( const IDomainSampler2 *p_ds,
                std::vector<Vec2> &vec_pos, std::vector<Vec2> &vec_normal );
void TransformPOB( std::vector<Vec2> &vec_pos, std::vector<Vec2> &vec_normal, const Transform2 &tr );

template <typename PairT>
void GUpdatePairs( const Vec2 *vec_p1_2, const Vec2 *vec_p2_2,
                   std::vector< PairT > &vec_pair );
/*TEMP: Debug version
void UpdatePairs( const Vec2 *vec_p1_2, const Vec2 *vec_p2_2,
                  const IDomainSampler2 *p_ds1,
                  const IDomainSampler2 *p_ds2,
                  const Transform2 &tr12,
                  std::vector< StochasticSPC2::Pair > &vec_pair );
*/
#else
void UpdatePairs( const IDomainSampler2 *p_ds1,
                  const IDomainSampler2 *p_ds2,
                  const Transform2 &tr12,
                  std::vector< StochasticSPC2::Pair > &vec_pair );
#endif

void RefinePairs_MinDist( IDomainSampler2 *p_ds1,
                          IDomainSampler2 *p_ds2,
                          const Transform2 &tr12,
                          Real neighbourhood_radius,
                          unsigned int num_neighbours,
                          std::vector< StochasticSPC2::Pair > &vec_pair
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                          , ContactData2::VizData &vd
#endif
    );
template <typename PairT>
void GRefinePairs_MinDist_Adaptive( IDomainSampler2 *p_ds1,
                                    IDomainSampler2 *p_ds2,
                                    const Transform2 &tr12,
                                    unsigned int num_neighbours,
                                    std::vector< PairT > &vec_pair
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                                      , ContactData2::VizData &vd
#endif
    );

void ReclassifyPairs_NP( IDomainSampler2 *p_ds1,
                         IDomainSampler2 *p_ds2,
                         const Transform2 &tr12,
                         StochasticSPC2 *p_spc,
                         const Context *p_context );
void ReclassifyPairs_IP_to_IP_or_NP( IDomainSampler2 *p_ds1,
                                     IDomainSampler2 *p_ds2,
                                     const Transform2 &tr12,
                                     StochasticSPC2 *p_spc,
                                     const Context *p_context );
void ReclassifyPairs_CP_to_CP_or_NP( IDomainSampler2 *p_ds1,
                                     IDomainSampler2 *p_ds2,
                                     const Transform2 &tr12,
                                     StochasticSPC2 *p_spc,
                                     const Context *p_context );

/* Could be a methods of StochasticSPC2... but that would make the name
   "Cache" a bit incoherent, as it StochasticSPC2 would be an object
   that does actual work, not just cached data...
*/
void MatchPairs_BF( IDomainSampler2 *p_ds1,  unsigned int num_pob1,
                    const point_on_boundary_id_type *vec_pob1, const Vec2 *vec_p1_2, const Vec2 *vec_n1_2,
                    IDomainSampler2 *p_ds2,  unsigned int num_pob2,
                    const point_on_boundary_id_type *vec_pob2, const Vec2 *vec_p2_2, const Vec2 *vec_n2_2,
                    const Transform2 &tr12,
                    StochasticSPC2 *p_spc,
                    const Context *p_context );

void SimplifyPairs_Discard( IDomainSampler2 *p_ds1,
                            IDomainSampler2 *p_ds2,
                            const Transform2 &tr12,
                            unsigned int max_pairs,
                            std::vector< StochasticSPC2::Pair > &vec_pair );
void SimplifyPairs_KeepClosest( IDomainSampler2 *p_ds1,
                                IDomainSampler2 *p_ds2,
                                const Transform2 &tr12,
                                unsigned int max_pairs,
                                std::vector< StochasticSPC2::Pair > &vec_pair );

void SimplifyPairs_IP( IDomainSampler2 *p_ds1,
                       IDomainSampler2 *p_ds2,
                       const Transform2 &tr12,
                       Real epsilon,
                       unsigned int max_pairs,
                       std::vector< StochasticSPC2::IntersectionPair > &vec_pair );

template <typename PairT>
void GSimplifyPairs_RemoveDuplicatedAndKeepClosest( IDomainSampler2 *p_ds1,
                                                   IDomainSampler2 *p_ds2,
                                                   const Transform2 &tr12,
                                                   Real epsilon,
                                                   unsigned int max_pairs,
                                                   std::vector< PairT > &vec_pair );

/*! \todo Should work on POD vectors and on sample-pair std::vector, choosing one or the other point POB1/POB2 and avoiding reevaluation
 */
bool ComputePCA( unsigned int num_points, const Vec2 *vec_points,
                 Vec2 &avg_position,
                 Real &eigen_value0, Vec2 &eigen_vector0 );


//--------------------------------------------------------------------------------------------------------------------------------
// IShape Vs IShape
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_IShape2_IShape2_Stochastic( const IShape2* p_shape1, const Transform2 &tr1, const Real *vec_dof1,
                                         const IShape2* p_shape2, const Transform2 &tr2, const Real *vec_dof2,
                                         ContactData2 &cd, ContactCache2 *p_cc,
                                         const Context *p_context )
{

#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
    cd.m_VD.Clear();
#endif

    GEO_ASSERT( 0 != p_cc ); //NEEDS a valid cache

    Transform2 tr12( tr2.Inverse() * tr1 );

    const unsigned int cNumRP( p_context->m_Stochastic_Num_RP );
    const unsigned int cNumIRP( p_context->m_Stochastic_Num_IRP ); //Init RP
    unsigned int num_rp(0);

    //\todo Determine potential overlap region (POR) between shapes (AABB or kDOP are the simplest POR)

    // Retrieve StochasticSPC2 (Create it if not available in p_cc)
    StochasticSPC2 *pSPC( 0 );
    if( 0 == p_cc->GetSPC() )
    {
        //\todo I/N/F distance should be based in global max_vel criteria or relative to object size
        pSPC = new StochasticSPC2( p_shape1->CreateDomainSampler(), p_shape2->CreateDomainSampler() );
        p_cc->SetSPC( pSPC );
        // Initialize cache with enough sample pairs
        num_rp = cNumIRP;
    }
    else
    {
        pSPC = static_cast<StochasticSPC2*>( p_cc->GetSPC() );
        num_rp = cNumRP;
        if( !p_context->m_Stochastic_Persistent ) pSPC->Clear();
    }
    // Shortcuts to DS1/DS2
    IDomainSampler2 *pDS1( pSPC->m_pDS1 );
    IDomainSampler2 *pDS2( pSPC->m_pDS2 );
    GEO_ASSERT( 0 != pDS1 && 0 != pDS2 );

    // Update DS DOF, shape won't change during CD \todo CONSIDER doing this only once for rigid objects, as their DOF do not change...
    pDS1->UpdateDOF( vec_dof1 );
    pDS2->UpdateDOF( vec_dof2 );

    //TEMP: To check IM deterministically, without Stochastic IP
    if( p_context->m_Stochastic_Fallback
        && eShape_MeshSolid2 == p_shape1->GetType()
        && eShape_MeshSolid2 == p_shape2->GetType() )
    {
        pSPC->Clear(); //NO PERSISTENCE HERE... BY NOW
        ContactData2 bf_cd;
        if( Contact_MeshSolid2_MeshSolid2_BruteForce( static_cast<const MeshSolidShape2*>(p_shape1), tr1, vec_dof1,
                                                      static_cast<const MeshSolidShape2*>(p_shape2), tr2, vec_dof2,
                                                      bf_cd, p_cc, p_context ) )
        {
            GEO_ASSERT( bf_cd.HasPOF() );
            // Fill vecIP with deterministic crossings
            for( unsigned int i=0; i<bf_cd.Size(); i++ )
            {
                pSPC->m_vecIP.push_back( StochasticSPC2::Pair( pDS1->CreatePOB( bf_cd.GetPOF1( i ) ),
                                                               pDS2->CreatePOB( bf_cd.GetPOF2( i ) ),
                                                               0 ) );
                // No need to inc/decref anything
            }
        }
    }
    else
    {
        //---- BEGIN OF STOCHASTIC PAIR SEARCH
        /*\todo To improve memory access pattern, consider:
          - Keep old pairs unchanged (IP, NP, whatever)
          - Update old pairs and store them, if < far, into a SINGLE temporary Pair array vecUP (updated pairs)
          - Process vecUP, refine necessary pairs individually and save refined to vecRP (if 1 UP => K RP) or just overwrite UP (if 1 UP => 1 RP)
          - Process vecRP, simplify, clusterize, whatever, and finally save them, overwritting the old pairs
        */
        //\todo Update, Refine and Reclassify pairs, filter with POR
        // Update all pairs locally

        /*\todo If POB were NOT owned by IDomainSampler2 but by each pair
          ISpecificPairwiseCache2, they could be updated ONLY ONCE,
          regardless of pairs. To do so, either IDomainSampler2 does NOT
          allocate all its samples but manage them in an external
          ISampleSet that can be iterated and cleared from outside OR, A
          LOT SIMPLER, IDomainSampler are EXCLUSIVE to the
          ISpecificPairwiseCache2 that use them, so they are NOT saved in
          IShape, but created there and returned..... This second option
          is a LOT easier and cleaner... may have some overhead if
          specific IDomainSampler creation is expensive (ex: build
          spatial classification strutcures..., but a lot cheaper to
          update and clear without interfering constantly with POB from
          other IDomainSampler "listeners" in multi-contact scenarios.

          pDS1->UpdatePOB( tr12 ); -> Update points (and normals?!) in given tr12 ref sys
          pDS2->UpdatePOB( Identity ); -> Update points (and normals?!) and leave in local ref sys 2
          => POB_Position/Normal COULD BE INLINE if stored explicitly in GIDomainSamplerD<2>... this would save A LOT of virtual calls...
          also, could do pDS1->All_POB_Position( vec_pos, tr12 )
          Better -> pDS1->AllPositions( vec_pos1 ) => local to ref sys 1, mult by tr12 if necessary
          -> pDS2->AllPositions( vec_pos2 ) => local to ref sys 2, no need to transform
          UpdatePairs(xxxx) -> Just udpate pair-specific data, assuming POB are up-to-date
        */
#ifdef __USE_CACHED_POB_DATA
        UpdatePOB( pDS1, pSPC->m_vecP1_2, pSPC->m_vecN1_2 );
        TransformPOB( pSPC->m_vecP1_2, pSPC->m_vecN1_2, tr12 ); //\todo Transform POBD1 with tr12
        UpdatePOB( pDS2, pSPC->m_vecP2_2, pSPC->m_vecN2_2 );
        GUpdatePairs<StochasticSPC2::IntersectionPair>( &pSPC->m_vecP1_2[0], &pSPC->m_vecP2_2[0], pSPC->m_vecIP );
        GUpdatePairs<StochasticSPC2::ClosestPair>( &pSPC->m_vecP1_2[0], &pSPC->m_vecP2_2[0], pSPC->m_vecCP );
        GUpdatePairs<StochasticSPC2::NearPair>( &pSPC->m_vecP1_2[0], &pSPC->m_vecP2_2[0], pSPC->m_vecNP );
        /*TEMP: Debug version
          UpdatePairs( &pSPC->m_vecP1_2[0], &pSPC->m_vecP2_2[0], pDS1, pDS2, tr12, pSPC->m_vecIP );
          UpdatePairs( &pSPC->m_vecP1_2[0], &pSPC->m_vecP2_2[0], pDS1, pDS2, tr12, pSPC->m_vecCP );
          UpdatePairs( &pSPC->m_vecP1_2[0], &pSPC->m_vecP2_2[0], pDS1, pDS2, tr12, pSPC->m_vecNP );
        */
#else
        UpdatePairs( pDS1, pDS2, tr12, pSPC->m_vecIP );
        UpdatePairs( pDS1, pDS2, tr12, pSPC->m_vecCP );
        UpdatePairs( pDS1, pDS2, tr12, pSPC->m_vecNP );
#endif
        // Refine all pairs WITHOUT reclassifying them
        //IP: Refine aggressively \todo OPT: Could iterate only the ones with an associated IM (important crossing IP), non-IM IP are not that critical
        for( unsigned int it_rip=0; it_rip < p_context->m_Stochastic_Refine_Num_Neighbours_IP; it_rip++ )
            GRefinePairs_MinDist_Adaptive<StochasticSPC2::IntersectionPair>( pDS1, pDS2, tr12,
                                                                             p_context->m_Stochastic_Refine_Num_Neighbours_IP, pSPC->m_vecIP
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                                                                             , cd.m_VD
#endif
                );
        //CP: unused by now...
        GRefinePairs_MinDist_Adaptive<StochasticSPC2::ClosestPair>( pDS1, pDS2, tr12,
                                                                    p_context->m_Stochastic_Refine_Num_Neighbours_CP, pSPC->m_vecCP
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                                                                    , cd.m_VD
#endif
            );
        //NP: Single refinement pass by now
        GRefinePairs_MinDist_Adaptive<StochasticSPC2::NearPair>( pDS1, pDS2, tr12,
                                                                 p_context->m_Stochastic_Refine_Num_Neighbours_NP, pSPC->m_vecNP
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                                                                 , cd.m_VD
#endif
            );

        //\todo MAYBE: Filter pairs with POR
        // Reclassify pairs
        ReclassifyPairs_IP_to_IP_or_NP( pDS1, pDS2, tr12, pSPC, p_context );
        ReclassifyPairs_CP_to_CP_or_NP( pDS1, pDS2, tr12, pSPC, p_context );
        ReclassifyPairs_NP( pDS1, pDS2, tr12, pSPC, p_context );

        /* Generate new random pairs using IShape2 domain samplers and add
           to NP if in range

           \todo MAYBE: Filter POB with POR

           \todo Consider resolution-adaptive sampling on pDS1 and
           pDS2, instead of fixed sample count (larger domains require
           more samples to achieve the same resolution)
        */
        point_on_boundary_id_type vec_pob1[num_rp];
        point_on_boundary_id_type vec_pob2[num_rp];
        pDS1->CreateRandomPOB( vec_pob1, num_rp );
        pDS2->CreateRandomPOB( vec_pob2, num_rp );
        // Get all POB_Position and transform to common refsys, using tr12 to transform from tr1 to tr2, and avoid transforming BOTH
        //\todo Consider using APOB_Position to get all vec_pX_X[] so that 2*num_rp virtual calls to POB_Position() are avoided!)
        Vec2 vec_p1_2[num_rp];
        Vec2 vec_p2_2[num_rp];
        Vec2 vec_n1_2[num_rp];
        Vec2 vec_n2_2[num_rp];
        for( unsigned int it_rp=0; it_rp < num_rp; it_rp++ )
        {
            vec_p1_2[it_rp] = tr12 * pDS1->POB_Position( vec_pob1[it_rp] );
            vec_p2_2[it_rp] = pDS2->POB_Position( vec_pob2[it_rp] );
            vec_n1_2[it_rp] = tr12.m_Rot * pDS1->POB_Normal( vec_pob1[it_rp] );
            vec_n2_2[it_rp] = pDS2->POB_Normal( vec_pob2[it_rp] );
        }
        // Match pairs, increasing matched POB refcount
        unsigned int num_prev_np( pSPC->m_vecNP.size() );
        MatchPairs_BF( pDS1, num_rp, vec_pob1, vec_p1_2, vec_n1_2,
                       pDS2, num_rp, vec_pob2, vec_p2_2, vec_n2_2,
                       tr12,
                       pSPC,
                       p_context );
        GEO_NP_STAT_SET( stochastic.m_Ratio_NP_div_RP, float(pSPC->m_vecNP.size() - num_prev_np) / num_rp );

        // Unref all generated POB to remove unmatched ones
        pDS1->APOB_DecRef( vec_pob1, num_rp );
        pDS2->APOB_DecRef( vec_pob2, num_rp );

        /*\todo Limit total IP and NP by discarding or merging/clustering
          them, possibly with different thresholds

          \note First merge POB on the same shape, then reindex and merge
          pairs, updating their distance.

          A) O(n log n) If there are still too many IP or NP, sort them
          and discard farthest ones, but AFTER merging the ones that are
          very close, OTHERWISE we would just keep a narrow band of
          redundant close shapes instead of a representativa band
          [Near,Far] of potential local minima neighbourhoods.

          B) O(n) Alternatively, classify NP into distance-range buckets,
          NOT NECESSARILY LINEAR, in O(n) time, and keep a REPRESENTATIVE
          subset, not just the closest ones. A Pair is representative if
          there is NO OTHER pair in the same distance-range THAT is also
          VERY CLOSE spatially AND has similar normals (so that they may
          be in the same refinement/cluster neighbourhood). Non-linear
          interval distribution in buckets could ensure more resolution
          around Near and less towards Far.
        */
        SimplifyPairs_IP( pDS1, pDS2,
                          tr12,
                          p_context->m_Stochastic_SimplifyPairs_EpsilonRel_IP * p_context->m_Stochastic_ZeroDist,
                          p_context->m_Stochastic_Max_IP,
                          pSPC->m_vecIP );
        GSimplifyPairs_RemoveDuplicatedAndKeepClosest<StochasticSPC2::ClosestPair>( pDS1, pDS2,
                                                                                    tr12,
                                                                                    p_context->m_Stochastic_SimplifyPairs_EpsilonRel_CP
                                                                                    * p_context->m_Stochastic_NearDist,
                                                                                    p_context->m_Stochastic_Max_CP,
                                                                                    pSPC->m_vecCP );
        GSimplifyPairs_RemoveDuplicatedAndKeepClosest<StochasticSPC2::NearPair>( pDS1, pDS2,
                                                                                 tr12,
                                                                                 p_context->m_Stochastic_SimplifyPairs_EpsilonRel_NP
                                                                                 * p_context->m_Stochastic_NearDist,
                                                                                 p_context->m_Stochastic_Max_NP,
                                                                                 pSPC->m_vecNP );
        //---- END OF STOCHASTIC PAIR SEARCH
    }

    //-- BEGIN Contact Manifold computation
#ifdef __ENABLE_IM
    Real cCrossingNeighbourRadius( Real(2) * p_context->m_Stochastic_ZeroDist ); //\todo Maybe should also be relative to ACTUAL object sizes...

    // Unflood all IP
    for( unsigned int it_ip=0; it_ip<pSPC->m_vecIP.size(); it_ip++ ) pSPC->m_vecIP[it_ip].SetFlooded(false);
    // Update existing IM, flooding any reachable IP
    for( unsigned int it_ip=0; it_ip<pSPC->m_vecIP.size(); it_ip++ )
        if( pSPC->m_vecIP[it_ip].HasIM() )
        {
            pSPC->m_vecIP[it_ip].DiscardIM( pDS1, pDS2 ); //TEMP: Destroy to force re-init, by now
            //pSPC->m_vecIP[it_ip].GetIM()->Update( pSPC->m_vecIP, it_ip, pDS1, pDS2 );
        }
    // Create new IM for unflooded IP
    for( unsigned int it_ip=0; it_ip<pSPC->m_vecIP.size(); it_ip++ )
    {
        StochasticSPC2::IntersectionPair &ip( pSPC->m_vecIP[it_ip] );
        if( !ip.IsFlooded() )
        {
            // Find actual crossing
            //-- Determine EXACT crossing and its local direction in both objects
            Vec2 pos1_2( tr12 * pDS1->POB_Position( ip.m_POB1 ) );
            Vec2 pos2_2( pDS2->POB_Position( ip.m_POB2 ) );
            // Get 2 surrounding points
            //neighbour pob1
            point_on_boundary_id_type vec_npob1[2];
            pDS1->CreateNeighbourPOB( &vec_npob1[0], 2,
                                      ip.m_POB1,
                                      cCrossingNeighbourRadius, p_context->m_IM_StepLength_Epsilon_Abs );
            Vec2 vec_npos1_2[2];
            vec_npos1_2[0] = tr12 * pDS1->POB_Position( vec_npob1[0] );
            vec_npos1_2[1] = tr12 * pDS1->POB_Position( vec_npob1[1] );
            GEO_ASSERT( mal::NormSq( vec_npos1_2[0] - vec_npos1_2[1] ) > Real(0) );
            //neighbour pob2
            point_on_boundary_id_type vec_npob2[2];
            pDS2->CreateNeighbourPOB( &vec_npob2[0], 2,
                                      ip.m_POB2,
                                      cCrossingNeighbourRadius, p_context->m_IM_StepLength_Epsilon_Abs );
            Vec2 vec_npos2_2[2];
            vec_npos2_2[0] = pDS2->POB_Position( vec_npob2[0] );
            vec_npos2_2[1] = pDS2->POB_Position( vec_npob2[1] );
            GEO_ASSERT( mal::NormSq( vec_npos2_2[0] - vec_npos2_2[1] ) > Real(0) );

            // Intersect segments/parabolas
            /*\todo Crossing can also be tested checking if either
              npos1_2[0] o npos1_2[1] is behind pos2_2 local plane
              (using its local normal OR a fitted plane that accounts
              for npos1_2[0],npos1_2[1] ).

              More accurately, both npos1_2 could be checked
              against an implicit parabola fitted to 3 points on
              object2 (pos1_2, npos2_2[0], npos2_2[1]) to determine if
              they are inside or outside.

              \note The symmetric case would be checked too.
            */
            Real lambda1, lambda2;
            if( Intersection_Segment2_Segment2( vec_npos1_2[0], vec_npos1_2[1],
                                                vec_npos2_2[0], vec_npos2_2[1],
                                                lambda1, lambda2 ) )
            {
                /* \todo Consider computing approx intersection point,
                   but it wouldn't be efficient if we needed to
                   compute its POB... probably unnecessary, just use
                   IP as the starting pair for the IM assuming it MAY
                   be slightly exterior.
                */

                /* Compute normals \todo It'd be more robust and save
                   virtual calls if crossing normals were aprox from
                   neighbourhood points, maybe returned in
                   CreateNeighbourPOB() as CW/CCW orientation of
                   neighbours is not guaranteed...
                */
                Vec2 n1_2( tr12.m_Rot * pDS1->POB_Normal( ip.m_POB1 ) );
                Vec2 n2_2( pDS2->POB_Normal( ip.m_POB2 ) );

                // Select stepping directions
                int s1(-1);
                if( mal::Dot( vec_npos1_2[0] - pos2_2, n2_2 ) < 0 ) s1 = 0;
                else if( mal::Dot( vec_npos1_2[1] - pos2_2, n2_2 ) < 0 ) s1 = 1;
                else GEO_LOG_ASSERT(false,"Cannot select s1");
                int s2(-1);
                if( mal::Dot( vec_npos2_2[0] - pos1_2, n1_2 ) < 0 ) s2 = 0;
                else if( mal::Dot( vec_npos2_2[1] - pos1_2, n1_2 ) < 0 ) s2 = 1;
                else GEO_LOG_ASSERT(false,"Cannot select s2");

#  ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IM ) )
                    cd.m_VD.m_vecIM_XPoint.push_back( tr2 * ( vec_npos1_2[0] + lambda1 * (vec_npos1_2[1]-vec_npos1_2[0]) ) );
#  endif

                /* Flood any unfloodedd IP in a
                   cCrossingNeighbourRadius neighbourhood.

                   \note We do NOT rely on SimplifyPairs_IP to avoid
                   flooding nearly-coincident IP, instead, once we
                   find an unflooded IP that actually represents a
                   crossing, we flood any other unflooded IP in a
                   cCrossingNeighbourRadius neighbourhood around it.
                */
                Vec2 midpoint_ip_2( Real(0.5) * ( pos1_2 + pos2_2 ) );
                for( unsigned int it_unflooded_ip=0; it_unflooded_ip<pSPC->m_vecIP.size(); it_unflooded_ip++ )
                    if( it_unflooded_ip != it_ip && !pSPC->m_vecIP[it_unflooded_ip].IsFlooded() )
                    {
                        StochasticSPC2::IntersectionPair &unflooded_ip( pSPC->m_vecIP[it_unflooded_ip] );
                        if( mal::NormSq( Real(0.5) * ( tr12 * pDS1->POB_Position( unflooded_ip.m_POB1 )
                                                       + pDS2->POB_Position( unflooded_ip.m_POB2 ) ) - midpoint_ip_2 )
                            < mal::Sq(cCrossingNeighbourRadius) )
                            pSPC->m_vecIP[it_unflooded_ip].SetFlooded(true);
                    }
                // Create new IM
                IntersectionMapping *pIM = new IntersectionMapping();
                pIM->Init( pSPC->m_vecIP, it_ip,
                           pos1_2, vec_npob1[s1], vec_npos1_2[s1],
                           pos2_2, vec_npob2[s2], vec_npos2_2[s2],
                           tr12, tr2,
                           pDS1, pDS2,
                           cCrossingNeighbourRadius, p_context->m_IM_StepLength, p_context->m_IM_StepLength_Epsilon_Abs, p_context->m_IM_Flags
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                           , cd.m_VD
#endif
                           );
                ip.SetIM( pIM );
            }
            else
            {
                //Skip IP if no actual crossing found, it's a very-close pair, but not intersecting
                ip.SetFlooded(true);
                GEO_NP_STAT_INC( stochastic.m_Num_FalseIP );
#ifdef __ENABLE_TRACE_PAIRS_IP
                GEO_LOG("IP with no actual crossing ignored, dist_sq = %e", ip.m_DistSq );
#endif
            }

            // DecRef neighbour POB
            pDS1->APOB_DecRef( vec_npob1, 2 ); //only 1-s1 will be destroyed, s1 was just IncRef inside IM::Init()
            pDS2->APOB_DecRef( vec_npob2, 2 ); //only 1-s2 will be destroyed, s2 was just IncRef inside IM::Init()
        }
    }

    // Compute ContactData from IM
    unsigned int num_im(0);
    cd.Begin();
    for( unsigned int it_ip=0; it_ip<pSPC->m_vecIP.size(); it_ip++ )
        if( pSPC->m_vecIP[it_ip].HasIM() )
        {
            pSPC->m_vecIP[it_ip].GetIM()->ComputeContactData( cd, tr2, pDS1, pDS2, p_context->m_Epsilon_ContactDepth );
            num_im++;
        }
    cd.SetNumDisjointManifolds( num_im );
    cd.End();

#else //__ENABLE_IM

    //TEMP: Compute contact data from stochastic pairs only
    cd.Begin();
    //CP based contact points
    for( unsigned int it_cp=0; it_cp<pSPC->m_vecCP.size(); it_cp++ )
    {
        //\todo Consider saving POB_Positions for IP so that they do not need to be re-queried and re-transformed to create the CP
        Vec2 p1_0( tr1 * pDS1->POB_Position( pSPC->m_vecCP[it_cp].m_POB1 ) );
        Vec2 p2_0( tr2 * pDS2->POB_Position( pSPC->m_vecCP[it_cp].m_POB2 ) );
        cd.AddCP( p1_0, p2_0, mal::Normalized(p1_0-p2_0), mal::Norm(p1_0-p2_0) );
        cd.AddFeaturePair( pDS1->POB_FeatureId( pSPC->m_vecCP[it_cp].m_POB1 ) ,
                           pDS2->POB_FeatureId( pSPC->m_vecCP[it_cp].m_POB2 ) );
    }
    cd.End();
#endif

#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_IP ) )
        for( unsigned int it_ip=0; it_ip<pSPC->m_vecIP.size(); it_ip++ )
            cd.m_VD.m_vecIP.push_back( std::make_pair( tr1 * pDS1->POB_Position( pSPC->m_vecIP[it_ip].m_POB1 ),
                                                       tr2 * pDS2->POB_Position( pSPC->m_vecIP[it_ip].m_POB2 ) ) );
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_CP ) )
        for( unsigned int it_cp=0; it_cp<pSPC->m_vecCP.size(); it_cp++ )
            cd.m_VD.m_vecCP.push_back( std::make_pair( tr1 * pDS1->POB_Position( pSPC->m_vecCP[it_cp].m_POB1 ),
                                                       tr2 * pDS2->POB_Position( pSPC->m_vecCP[it_cp].m_POB2 ) ) );
    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_NP ) )
        for( unsigned int it_np=0; it_np<pSPC->m_vecNP.size(); it_np++ )
        {
            std::pair<Vec2,Vec2> pob_pair( std::make_pair( tr1 * pDS1->POB_Position( pSPC->m_vecNP[it_np].m_POB1 ),
                                                           tr2 * pDS2->POB_Position( pSPC->m_vecNP[it_np].m_POB2 ) ) );
            if( pSPC->m_vecNP[it_np].m_DistSq <= mal::Sq( p_context->m_Stochastic_NearDist ) )
                cd.m_VD.m_vecNP.push_back( pob_pair );
            else
                cd.m_VD.m_vecFP.push_back( pob_pair );
        }

    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_RNP ) )
        for( unsigned int it_rnp=0; it_rnp<cd.m_VD.m_vecRNP.size(); it_rnp++ )
        {
            // Globalize RNP
            cd.m_VD.m_vecRNP[it_rnp].first = tr1 * cd.m_VD.m_vecRNP[it_rnp].first;
            cd.m_VD.m_vecRNP[it_rnp].second = tr2 * cd.m_VD.m_vecRNP[it_rnp].second;
        }

    if( params.debug.np.contact.stochastic.DDF.Test( PARAMS::DEBUG::NP::CONTACT::STOCHASTIC::eDraw_PCA ) )
    {
        // Testing PCA
        std::vector< Vec2 > vec_points;
        for( unsigned int it_ip=0; it_ip<pSPC->m_vecIP.size(); it_ip++ )
        {
            vec_points.push_back( tr1 * pDS1->POB_Position( pSPC->m_vecIP[it_ip].m_POB1 ) );
            vec_points.push_back( tr2 * pDS2->POB_Position( pSPC->m_vecIP[it_ip].m_POB2 ) );
        }
        for( unsigned int it_cp=0; it_cp<pSPC->m_vecCP.size(); it_cp++ )
        {
            vec_points.push_back( tr1 * pDS1->POB_Position( pSPC->m_vecCP[it_cp].m_POB1 ) );
            vec_points.push_back( tr2 * pDS2->POB_Position( pSPC->m_vecCP[it_cp].m_POB2 ) );
        }
        for( unsigned int it_np=0; it_np<pSPC->m_vecNP.size(); it_np++ )
        {
            vec_points.push_back( tr1 * pDS1->POB_Position( pSPC->m_vecNP[it_np].m_POB1 ) );
            vec_points.push_back( tr2 * pDS2->POB_Position( pSPC->m_vecNP[it_np].m_POB2 ) );
        }
        Vec2 avg_pos;
        Real pca_value0;
        Vec2 pca_dir0;
        if( ComputePCA( vec_points.size(), &vec_points[0], avg_pos, pca_value0, pca_dir0 ) )
        {
            //GEO_LOG_WARNING("Adding PCA");
            cd.m_VD.m_vecPCA.push_back( std::make_pair( avg_pos, pca_dir0 * 10 ) );
        }
    }
#endif //__ENABLE_GEO_NP_CONTACTDATA_VIZDATA

#ifdef __ENABLE_TRACE_PAIRS
    GEO_LOG( "Generated %d RP, kept %d IP + %d CP + %d NP",
             num_rp, (int)pSPC->m_vecIP.size(), (int)pSPC->m_vecCP.size(), (int)pSPC->m_vecNP.size() );
    pDS1->Trace();
    pDS2->Trace();
#endif

#ifdef __ENABLE_COMPARE_WITH_MESH2_MESH2_BRUTEFORCE
    if( eShape_MeshSolid2 == p_shape1->GetType()
        && eShape_MeshSolid2 == p_shape2->GetType() )
    {
        bool bContactST( cd.Size() > 0 );
        ContactData2 bf_cd;
        bool bContactBF = Contact_MeshSolid2_MeshSolid2_BruteForce( static_cast<const MeshSolidShape2*>(p_shape1), tr1, vec_dof1,
                                                                    static_cast<const MeshSolidShape2*>(p_shape2), tr2, vec_dof2,
                                                                    bf_cd, p_cc, p_context );
        /*TEMP: verbose...
        if( bContactBF != bContactST )
        {
            GEO_LOG_ERROR("Reference contact MISSED");
        }
        else if( bf_cd.m_NumDisjointManifolds != cd.m_NumDisjointManifolds )
        {
            GEO_LOG_WARNING("Reference contact mismatch: #ref.dm = %d != #stoch.dm = %d", bf_cd.m_NumDisjointManifolds, cd.m_NumDisjointManifolds );
        }
        */
        GEO_NP_STAT_SET( stochastic.m_Num_UDCM, bf_cd.m_NumDisjointManifolds - cd.m_NumDisjointManifolds );
    }
#endif

    GEO_NP_STAT_SET( stochastic.m_Num_NP, pSPC->m_vecNP.size() );
    GEO_NP_STAT_SET( stochastic.m_Num_IP, pSPC->m_vecIP.size() );
    GEO_NP_STAT_SET( stochastic.m_Num_POB, pDS1->GetNumPOB() + pDS2->GetNumPOB() );
    GEO_NP_STAT_SET( stochastic.m_Num_Alloc_POB, pDS1->GetNumAllocatedPOB() + pDS2->GetNumAllocatedPOB() );

    return cd.Size() > 0;
}

//----------------------------------------------------------------
//-- Pair-processing functions
//----------------------------------------------------------------

#ifdef __USE_CACHED_POB_DATA
void UpdatePOB( const IDomainSampler2 *p_ds, std::vector<Vec2> &vec_pos, std::vector<Vec2> &vec_normal )
{
    vec_pos.resize( p_ds->GetNumAllocatedPOB() );
    vec_normal.resize( p_ds->GetNumAllocatedPOB() );
    p_ds->AllPositionsAndNormals( &vec_pos[0], &vec_normal[0] );
}
void TransformPOB( std::vector<Vec2> &vec_pos, std::vector<Vec2> &vec_normal, const Transform2 &tr )
{
    //\todo This transforms all POB unconditionally, even the free ones... but it's not a serious problem if #free << #valid
    GEO_ASSERT( vec_pos.size() == vec_normal.size() );
    for( unsigned int i=0; i<vec_pos.size(); i++ )
    {
        vec_pos[i] = tr * vec_pos[i];
        vec_normal[i] = tr.m_Rot * vec_normal[i];
    }
}

template <typename PairT>
void GUpdatePairs( const Vec2 *vec_p1_2, const Vec2 *vec_p2_2,
                   std::vector< PairT > &vec_pair )
{
    for( unsigned int it_pairs=0; it_pairs<vec_pair.size(); it_pairs++ )
    {
        PairT &pair( vec_pair[it_pairs] );
        pair.m_DistSq = mal::NormSq( vec_p1_2[pair.m_POB1] - vec_p2_2[pair.m_POB2] );
    }
}
/*TEMP: Debug version
void UpdatePairs( const Vec2 *vec_p1_2, const Vec2 *vec_p2_2,
                  const IDomainSampler2 *p_ds1,
                  const IDomainSampler2 *p_ds2,
                  const Transform2 &tr12,
                  std::vector< StochasticSPC2::Pair > &vec_pair )
{
    for( unsigned int it_pairs=0; it_pairs<vec_pair.size(); it_pairs++ )
    {
        StochasticSPC2::Pair &pair( vec_pair[it_pairs] );
        Vec2 p1_2( tr12 * p_ds1->POB_Position( pair.m_POB1 ) );
        Vec2 p2_2( p_ds2->POB_Position( pair.m_POB2 ) );
        Real dist_sq( mal::NormSq(p1_2 - p2_2) );
        pair.m_DistSq = dist_sq;
        //TEMP: Check
        Vec2 tmp_p1_2( vec_p1_2[pair.m_POB1] );
        Vec2 tmp_p2_2( vec_p2_2[pair.m_POB2] );
        Real tmp_dist_sq( mal::NormSq( tmp_p1_2 - tmp_p2_2 ) );
        //GEO_ASSERT( pair.m_DistSq == tmp_dist_sq ); //\todo TOO RESTRICTIVE, use epsilon
        pair.m_DistSq = tmp_dist_sq;
    }
}
*/
#else
void UpdatePairs( const IDomainSampler2 *p_ds1,
                  const IDomainSampler2 *p_ds2,
                  const Transform2 &tr12,
                  std::vector< StochasticSPC2::Pair > &vec_pair )
{
    for( unsigned int it_pairs=0; it_pairs<vec_pair.size(); it_pairs++ )
    {
        StochasticSPC2::Pair &pair( vec_pair[it_pairs] );
        Vec2 p1_2( tr12 * p_ds1->POB_Position( pair.m_POB1 ) );
        Vec2 p2_2( p_ds2->POB_Position( pair.m_POB2 ) );
        Real dist_sq( mal::NormSq(p1_2 - p2_2) );
        pair.m_DistSq = dist_sq;
        //\todo Update positiions if explicitly stored!
    }
}
#endif

/*! BruteForce closest pair, could use 1D sorting/bucketing or 2D/3D spatial hashing to reduce O(n²) */
Real ComputeClosestPair_BF( unsigned int num_pob1, const point_on_boundary_id_type *vec_pob1, const Vec2 *vec_p1_2,
                            unsigned int num_pob2, const point_on_boundary_id_type *vec_pob2, const Vec2 *vec_p2_2,
                            point_on_boundary_id_type &pob1, point_on_boundary_id_type &pob2 )
{
    pob1 = vec_pob1[0];
    pob2 = vec_pob2[0];
    Real min_dist_sq( mal::NormSq( vec_p1_2[0] - vec_p2_2[0] ) );
    for( unsigned int it_pob1 = 0; it_pob1 < num_pob1; it_pob1++ )
        for( unsigned int it_pob2 = 0; it_pob2 < num_pob2; it_pob2++ )
        {
            Real dist_sq( mal::NormSq( vec_p1_2[it_pob1] - vec_p2_2[it_pob2] ) );
            if( dist_sq < min_dist_sq )
            {
                min_dist_sq = dist_sq;
                pob1 = vec_pob1[it_pob1];
                pob2 = vec_pob2[it_pob2];
            }
        }
    return min_dist_sq;
}


/* Refine pairs using min-distance criterion

\todo Consider using pair.m_DistSq as the neighbourhood radius,
instead of a fixed value, specially for NP, which can get arbitrarily
closer and local refinement does nor make sense beyond its actual distance!!

\todo At least for FP, use actual m_DistSq, not always FarDist, which
can be QUITE BIG
*/
void RefinePairs_MinDist( IDomainSampler2 *p_ds1,
                          IDomainSampler2 *p_ds2,
                          const Transform2 &tr12,
                          Real neighbourhood_radius,
                          unsigned int num_neighbours,
                          std::vector< StochasticSPC2::Pair > &vec_pair
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                          , ContactData2::VizData &vd
#endif
    )
{
    //\todo By now, we generate the same number of neighbours on both DS
    const unsigned int num_pob( num_neighbours + 1 );
    // Iterate over original IP and refine (substitute) them with local optimal
    unsigned int num_pairs( vec_pair.size() );
    for( unsigned int it_pair=0; it_pair<num_pairs; it_pair++ )
    {
        StochasticSPC2::Pair &pair( vec_pair[it_pair] );
        // Generate K random samples around PAIR in a given neighbourhood
        point_on_boundary_id_type vec_pob1[num_pob];
        point_on_boundary_id_type vec_pob2[num_pob];
        p_ds1->CreateRandomNeighbourPOB( vec_pob1, num_neighbours,
                                         neighbourhood_radius, pair.m_POB1 );
        p_ds2->CreateRandomNeighbourPOB( vec_pob2, num_neighbours,
                                         neighbourhood_radius, pair.m_POB2 );
        // Add current POB
        vec_pob1[num_pob-1] = pair.m_POB1;
        vec_pob2[num_pob-1] = pair.m_POB2;
        // Compute POB positions wrt tr2 (\todo Optimize with APOB_Position)
        Vec2 vec_p1_2[num_pob];
        Vec2 vec_p2_2[num_pob];
        for( unsigned int it_pob=0; it_pob < num_pob; it_pob++ )
        {
            vec_p1_2[it_pob] = tr12 * p_ds1->POB_Position( vec_pob1[it_pob] );
            vec_p2_2[it_pob] = p_ds2->POB_Position( vec_pob2[it_pob] );
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
            /*TEMP: expensive viz, avoid it if possible
            vd.m_vecRNP.push_back( std::make_pair( p_ds1->POB_Position( vec_pob1[it_pob] ),
                                                   p_ds2->POB_Position( vec_pob2[it_pob] ) ) );
            */
#endif
        }
        /*TEMPORAL: by now, we JUST keep the best refined PAIR, not all
          < m_ZeroDistSq, this may change...
          \todo: Maybe match ALL refined sample-pairs < neighbourhood_radius,
          create new PAIR if found, add to end. This would be important
          in 3D, where PAIR would not be isolated but part of spatial
          curves.
        */
        // Compute closest pair and save it
        Real dist_sq = ComputeClosestPair_BF( num_pob, vec_pob1, vec_p1_2,
                                              num_pob, vec_pob2, vec_p2_2,
                                              pair.m_POB1, pair.m_POB2 );
        /*TEMP
        if( dist_sq < pair.m_DistSq ) { GEO_LOG_WARNING( "Refined pair %f => %f", pair.m_DistSq, dist_sq ); }
        else { GEO_ASSERT( pair.m_DistSq >= dist_sq ); }
        TEMP*/
        pair.m_DistSq = dist_sq;
        p_ds1->POB_IncRef( pair.m_POB1 );
        p_ds2->POB_IncRef( pair.m_POB2 );
        //\todo UPDATE pair POSITIONS if they are stored

        // Unref all generated POB, as well as old PAIR to remove unmatched ones
        p_ds1->APOB_DecRef( vec_pob1, num_pob );
        p_ds2->APOB_DecRef( vec_pob2, num_pob );
    }
}

/*! Refine pairs using adaptive min-distance criterion, with current
   distance as the neighbourhood radius.

   \note This has been observed to improve convergence, specially for
   NP, which can get arbitrarily closer and local refinement does nor
   make sense beyond its actual distance!!

   \todo As it requires an sqrt() per refined pair, MAYBE we should
   only use RefinePairs_MinDist_Adaptive() for specific/important
   pairs or pair-sets and use RefinePairs_MinDist() for the other
   ones...

   \todo CONSIDER using a variable number of neighbours on EACH
   IDomainSampler2 to ensure a minimum sampling density in their
   adaptive local neighbourhood.
*/
template <typename PairT>
void GRefinePairs_MinDist_Adaptive( IDomainSampler2 *p_ds1,
                                    IDomainSampler2 *p_ds2,
                                    const Transform2 &tr12,
                                    unsigned int num_neighbours,
                                    std::vector< PairT > &vec_pair
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                                    , ContactData2::VizData &vd
#endif
    )
{
    //\todo By now, we generate the same number of neighbours on both DS
    const unsigned int num_pob( num_neighbours + 1 );
    // Iterate over original pairs and refine (substitute) them with local optimal
    unsigned int num_pairs( vec_pair.size() );
    for( unsigned int it_pair=0; it_pair<num_pairs; it_pair++ )
    {
        PairT &pair( vec_pair[it_pair] );
        // Generate K random samples around PAIR in a given neighbourhood
        point_on_boundary_id_type vec_pob1[num_pob];
        point_on_boundary_id_type vec_pob2[num_pob];
        Real neighbourhood_radius( mal::Sqrt( pair.m_DistSq ) );
        p_ds1->CreateRandomNeighbourPOB( vec_pob1, num_neighbours,
                                         neighbourhood_radius, pair.m_POB1 );
        p_ds2->CreateRandomNeighbourPOB( vec_pob2, num_neighbours,
                                         neighbourhood_radius, pair.m_POB2 );
        // Add current POB
        vec_pob1[num_pob-1] = pair.m_POB1;
        vec_pob2[num_pob-1] = pair.m_POB2;
        // Compute POB positions wrt tr2 (\todo Optimize with APOB_Position)
        Vec2 vec_p1_2[num_pob];
        Vec2 vec_p2_2[num_pob];
        for( unsigned int it_pob=0; it_pob < num_pob; it_pob++ )
        {
            vec_p1_2[it_pob] = tr12 * p_ds1->POB_Position( vec_pob1[it_pob] );
            vec_p2_2[it_pob] = p_ds2->POB_Position( vec_pob2[it_pob] );
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
            /*TEMP: expensive viz, avoid it if possible
            vd.m_vecRNP.push_back( std::make_pair( p_ds1->POB_Position( vec_pob1[it_pob] ),
                                                   p_ds2->POB_Position( vec_pob2[it_pob] ) ) );
            */
#endif
        }
        /*TEMPORAL: by now, we JUST keep the best refined PAIR, not all
          < m_ZeroDistSq, this may change...
          \todo: Maybe match ALL refined sample-pairs < m_ZeroDistSq,
          create new PAIR if found, add to end. This would be important
          in 3D, where PAIR would not be isolated but part of spatial
          curves.
        */
        // Compute closest pair and save it
        Real dist_sq = ComputeClosestPair_BF( num_pob, vec_pob1, vec_p1_2,
                                              num_pob, vec_pob2, vec_p2_2,
                                              pair.m_POB1, pair.m_POB2 );
        /*TEMP
        if( dist_sq < pair.m_DistSq ) { GEO_LOG_WARNING( "Refined pair %f => %f", pair.m_DistSq, dist_sq ); }
        else { GEO_ASSERT( pair.m_DistSq >= dist_sq ); }
        TEMP*/
        pair.m_DistSq = dist_sq;
        p_ds1->POB_IncRef( pair.m_POB1 );
        p_ds2->POB_IncRef( pair.m_POB2 );
        //\todo UPDATE pair POSITIONS if they are stored

        // Unref all generated POB, as well as old PAIR to remove unmatched ones
        p_ds1->APOB_DecRef( vec_pob1, num_pob );
        p_ds2->APOB_DecRef( vec_pob2, num_pob );
    }
}

void ReclassifyPairs_NP( IDomainSampler2 *p_ds1,
                         IDomainSampler2 *p_ds2,
                         const Transform2 &tr12,
                         StochasticSPC2 *p_spc,
                         const Context *p_context )
{
    Real zero_dist_sq( mal::Sq( p_context->m_Stochastic_ZeroDist ) );
    Real near_dist_sq( mal::Sq( p_context->m_Stochastic_NearDist ) );
    Real far_dist_sq( mal::Sq( p_context->m_Stochastic_FarDist ) );
    unsigned int it_np(0);
    while( it_np < p_spc->m_vecNP.size() )
    {
        StochasticSPC2::Pair &np( p_spc->m_vecNP[it_np] );
        if( np.m_DistSq <= far_dist_sq )
        {
            bool bRemove(true);
            if( np.m_DistSq <= zero_dist_sq )
                p_spc->m_vecIP.push_back( np ); // Convert to IP
            else
            {
                // Check if opposed and behind, in tr2 ref sys
                Vec2 p1( tr12 * p_ds1->POB_Position( np.m_POB1 ) );
                Vec2 p2( p_ds2->POB_Position( np.m_POB2 ) );
                Vec2 n1( tr12.m_Rot * p_ds1->POB_Normal( np.m_POB1 ) );
                Vec2 n2( p_ds2->POB_Normal( np.m_POB2 ) );
                bool bOpposed( n1 * n2 < 0 );
                bool bInFront( (p1-p2) * n2 > 0 && (p2-p1) * n1 > 0 );
                if( bOpposed && bInFront && np.m_DistSq <= near_dist_sq ) // Convert to CP
                    p_spc->m_vecCP.push_back( np );
                else // Keep NP
                {
                    ++it_np;
                    bRemove = false;
                }
            }
            if( bRemove )
            {
                np = p_spc->m_vecNP.back();
                p_spc->m_vecNP.pop_back();
            }
        }
        else
        {
            // Unref and remove NP
            p_ds1->POB_DecRef( np.m_POB1 );
            p_ds2->POB_DecRef( np.m_POB2 );
            np = p_spc->m_vecNP.back();
            p_spc->m_vecNP.pop_back();
        }
    }
}

/* IP => IP | NP | Discard
*/
void ReclassifyPairs_IP_to_IP_or_NP( IDomainSampler2 *p_ds1,
                                     IDomainSampler2 *p_ds2,
                                     const Transform2 &tr12,
                                     StochasticSPC2 *p_spc,
                                     const Context *p_context )
{
    Real zero_dist_sq( mal::Sq( p_context->m_Stochastic_ZeroDist ) );
    Real near_dist_sq( mal::Sq( p_context->m_Stochastic_NearDist ) );
    unsigned int it_ip(0);
    while( it_ip < p_spc->m_vecIP.size() )
    {
        StochasticSPC2::IntersectionPair &ip( p_spc->m_vecIP[it_ip] );
        if( ip.m_DistSq <= zero_dist_sq ) //Keep IP
            ++it_ip;
        else
        {
            if( ip.m_DistSq <= near_dist_sq )
                p_spc->m_vecNP.push_back( ip );
            else
            {
                p_ds1->POB_DecRef( ip.m_POB1 );
                p_ds2->POB_DecRef( ip.m_POB2 );
            }
            // Remove IP
            ip.DiscardIM( p_ds1, p_ds2 ); //Discard potential IM and its internal IB
            ip = p_spc->m_vecIP.back();
            p_spc->m_vecIP.pop_back();
        }
    }
}

/* CP => CP | NP | Discard
*/
void ReclassifyPairs_CP_to_CP_or_NP( IDomainSampler2 *p_ds1,
                                     IDomainSampler2 *p_ds2,
                                     const Transform2 &tr12,
                                     StochasticSPC2 *p_spc,
                                     const Context *p_context )
{
    Real zero_dist_sq( mal::Sq( p_context->m_Stochastic_ZeroDist ) );
    Real near_dist_sq( mal::Sq( p_context->m_Stochastic_NearDist ) );
    unsigned int it_cp(0);
    while( it_cp < p_spc->m_vecCP.size() )
    {
        StochasticSPC2::Pair &cp( p_spc->m_vecCP[it_cp] );
        if( cp.m_DistSq <= near_dist_sq )
        {
            // Check (lazily) if opposed and in front, in tr2 ref sys
            Vec2 n1( tr12.m_Rot * p_ds1->POB_Normal( cp.m_POB1 ) );
            Vec2 n2( p_ds2->POB_Normal( cp.m_POB2 ) );
            bool bOpposed( n1 * n2 < 0 );
            bool bKeep(false);
            if( bOpposed )
            {
                Vec2 p1( tr12 * p_ds1->POB_Position( cp.m_POB1 ) );
                Vec2 p2( p_ds2->POB_Position( cp.m_POB2 ) );
                bool bInFront( (p1-p2) * n2 > 0 && (p2-p1) * n1 > 0 );
                bKeep = bInFront;
            }
            // Keep CP or convert to NP
            if( bKeep && cp.m_DistSq > zero_dist_sq )
                ++it_cp;
            else
            {
                p_spc->m_vecNP.push_back( cp );
                cp = p_spc->m_vecCP.back();
                p_spc->m_vecCP.pop_back();
            }
        }
        else
        {
            // Unref and remove CP
            p_ds1->POB_DecRef( cp.m_POB1 );
            p_ds2->POB_DecRef( cp.m_POB2 );
            cp = p_spc->m_vecCP.back();
            p_spc->m_vecCP.pop_back();
        }
    }
}

/*! Generate best NP in both directions.
  \note BruteForce implementation, could use 1D sorting/bucketing or 2D/3D spatial hashing to reduce O(n²)
  \note CAN use the same POB multiple times
  \note Current implementation may generate the same pair twice
  \note Written for clarity, not for speed

  \todo Generating IP is highly improbable so we don't check if < m_ZeroDist, we'll refine and reclassify it in the next frame if so
  \todo Should check bOpposed && bInFront for NP
*/
//#define __USE_EXISTING_IP_FOR_RP //\todo If disabled, matching is FASTER (less candidates and NO virtual calls too) but less effective to exploit cached POB
//#define __USE_EXISTING_NP_FOR_RP  //\todo If disabled, matching is FASTER (less candidates and NO virtual calls too) but less effective to exploit cached POB
//#define __REQUIRE_RP_OPPOSED_A_PRIORI //\todo This discards some RP which is actually a good thing, because we DO NOT WANT IP/NP/PP not opposed in ANY CASE, they CANNOT converge to a valid IP/NP
//#define __REQUIRE_RP_IN_FRONT_A_PRIORI //\todo This reduces RP effectiveness CONSIDERABLY because around half of the sample POB are inside and are not given the opportunity to converge...
//#define __ENABLE_REQUIRE_RP_OPPOSED_AND_IN_FRONT_A_POSTERIORI //\todo This completely kills RP effectiveness, just don't enable it
void MatchPairs_BF( IDomainSampler2 *p_ds1,  unsigned int num_pob1,
                    const point_on_boundary_id_type *vec_pob1, const Vec2 *vec_p1_2, const Vec2 *vec_n1_2,
                    IDomainSampler2 *p_ds2,  unsigned int num_pob2,
                    const point_on_boundary_id_type *vec_pob2, const Vec2 *vec_p2_2, const Vec2 *vec_n2_2,
                    const Transform2 &tr12,
                    StochasticSPC2 *p_spc,
                    const Context *p_context )
{
    Real near_dist_sq( mal::Sq( p_context->m_Stochastic_NearDist ) );
    // Match 1 to 2
    for( unsigned int it_pob1=0; it_pob1 < num_pob1; it_pob1++ )
    {
        //Find matching samples across Random, Near and Intersecting POB
        point_on_boundary_id_type min_dist_pob2(0);
        Real min_dist_sq( 2*near_dist_sq ); // Larger than largest distance we accept
        for( unsigned int it_pob2=0; it_pob2<num_pob2; it_pob2++ )
        {
            Real dist_sq( mal::NormSq( vec_p1_2[it_pob1] - vec_p2_2[it_pob2] ) );
            if( dist_sq < min_dist_sq )
            {
#ifdef __REQUIRE_RP_OPPOSED_A_PRIORI
                Vec2 n1( vec_n1_2[it_pob1] );
                Vec2 n2( vec_n2_2[it_pob2] );
                bool bOpposed( n1 * n2 < 0 );
                if( bOpposed )
                {
#ifdef __REQUIRE_RP_IN_FRONT_A_PRIORI  //THIS reduces effectiveness A LOT
                    Vec2 p1( vec_p1_2[it_pob1] );
                    Vec2 p2( vec_p2_2[it_pob2] );
                    bool bInFront( (p1 - p2) * n2 > 0 && (p2 - p1) * n1 > 0 );
                    if( bInFront )
#endif
                    {
                        min_dist_sq = dist_sq;
                        min_dist_pob2 = vec_pob2[it_pob2];
                    }
                }
#else
                min_dist_sq = dist_sq;
                min_dist_pob2 = vec_pob2[it_pob2];
#endif
            }
        }
#ifdef __USE_EXISTING_IP_FOR_RP
        for( unsigned int it_ip=0; it_ip<p_spc->m_vecIP.size(); it_ip++ )
        {
            point_on_boundary_id_type pob2( p_spc->m_vecIP[it_ip].m_POB2 );
            Real dist_sq( mal::NormSq( vec_p1_2[it_pob1] - p_ds2->POB_Position( pob2 ) ) );
            if( dist_sq < min_dist_sq )
            {
                min_dist_sq = dist_sq;
                min_dist_pob2 = pob2;
            }
        }
#endif
#ifdef __USE_EXISTING_NP_FOR_RP
        for( unsigned int it_np=0; it_np<p_spc->m_vecNP.size(); it_np++ )
        {
            point_on_boundary_id_type pob2( p_spc->m_vecNP[it_np].m_POB2 );
            Real dist_sq( mal::NormSq( vec_p1_2[it_pob1] - p_ds2->POB_Position( pob2 ) ) );
            if( dist_sq < min_dist_sq )
            {
                min_dist_sq = dist_sq;
                min_dist_pob2 = pob2;
            }
        }
#endif
        // Generate only best pair
        if( min_dist_sq <= near_dist_sq )
        {
#ifdef __ENABLE_REQUIRE_RP_OPPOSED_AND_IN_FRONT_A_POSTERIORI
            // Check if opposed and in front
            Vec2 n1( vec_n1_2[it_pob1] );
            Vec2 n2( p_ds2->POB_Normal( min_dist_pob2 ) );
            bool bOpposed( n1 * n2 < 0 );
            if( bOpposed )
            {
                Vec2 p1( vec_p1_2[it_pob1] );
                Vec2 p2( p_ds2->POB_Position( min_dist_pob2 ) );
                bool bInFront( (p1 - p2) * n2 > 0 && (p2 - p1) * n1 > 0 );
                if( bInFront )
                {
                    StochasticSPC2::Pair np;
                    np.m_DistSq = min_dist_sq;
                    np.m_POB1 = vec_pob1[it_pob1];
                    np.m_POB2 = min_dist_pob2;
                    p_ds1->POB_IncRef( np.m_POB1 );
                    p_ds2->POB_IncRef( np.m_POB2 );
                    p_spc->m_vecNP.push_back( np );
                }
            }
#else
            StochasticSPC2::Pair np;
            np.m_DistSq = min_dist_sq;
            np.m_POB1 = vec_pob1[it_pob1];
            np.m_POB2 = min_dist_pob2;
            p_ds1->POB_IncRef( np.m_POB1 );
            p_ds2->POB_IncRef( np.m_POB2 );
            p_spc->m_vecNP.push_back( np );
#endif
        }
    }
    // Match 2 to 1
    for( unsigned int it_pob2=0; it_pob2 < num_pob2; it_pob2++ )
    {
        //Find matching samples across Random, Near and Intersecting POB
        point_on_boundary_id_type min_dist_pob1(0);
        Real min_dist_sq( 2*near_dist_sq ); // Larger than largest distance we accept
        for( unsigned int it_pob1=0; it_pob1<num_pob1; it_pob1++ )
        {
            Real dist_sq( mal::NormSq( vec_p2_2[it_pob2] - vec_p1_2[it_pob1] ) );
            if( dist_sq < min_dist_sq )
            {
#ifdef __REQUIRE_RP_OPPOSED_A_PRIORI
                Vec2 n1( vec_n1_2[it_pob1] );
                Vec2 n2( vec_n2_2[it_pob2] );
                bool bOpposed( n1 * n2 < 0 );
                if( bOpposed )
                {
#ifdef __REQUIRE_RP_IN_FRONT_A_PRIORI  //THIS reduces effectiveness A LOT
                    Vec2 p1( vec_p1_2[it_pob1] );
                    Vec2 p2( vec_p2_2[it_pob2] );
                    bool bInFront( (p1 - p2) * n2 > 0 && (p2 - p1) * n1 > 0 );
                    if( bInFront )
#endif
                    {
                        min_dist_sq = dist_sq;
                        min_dist_pob1 = vec_pob1[it_pob1];
                    }
                }
#else
                min_dist_sq = dist_sq;
                min_dist_pob1 = vec_pob1[it_pob1];
#endif
            }
        }
#ifdef __USE_EXISTING_IP_FOR_RP
        for( unsigned int it_ip=0; it_ip<p_spc->m_vecIP.size(); it_ip++ )
        {
            point_on_boundary_id_type pob1( p_spc->m_vecIP[it_ip].m_POB1 );
            Real dist_sq( mal::NormSq( vec_p2_2[it_pob2] - tr12 * p_ds1->POB_Position( pob1 ) ) );
            if( dist_sq < min_dist_sq )
            {
                min_dist_sq = dist_sq;
                min_dist_pob1 = pob1;
            }
        }
#endif
#ifdef __USE_EXISTING_NP_FOR_RP
        for( unsigned int it_np=0; it_np<p_spc->m_vecNP.size(); it_np++ )
        {
            point_on_boundary_id_type pob1( p_spc->m_vecNP[it_np].m_POB1 );
            Real dist_sq( mal::NormSq( vec_p2_2[it_pob2] - tr12 * p_ds1->POB_Position( pob1 ) ) );
            if( dist_sq < min_dist_sq )
            {
                min_dist_sq = dist_sq;
                min_dist_pob1 = pob1;
            }
        }
#endif
        // Generate only best pair
        if( min_dist_sq <= near_dist_sq )
        {
#ifdef __ENABLE_REQUIRE_RP_OPPOSED_AND_IN_FRONT_A_POSTERIORI
            // Check if opposed and in front
            Vec2 n1( tr12 * p_ds1->POB_Normal( min_dist_pob1 ) );
            Vec2 n2( vec_n2_2[it_pob2] );
            bool bOpposed( n1 * n2 < 0 );
            if( bOpposed )
            {
                Vec2 p1( tr12 * p_ds1->POB_Position( min_dist_pob1 ) );
                Vec2 p2( vec_p2_2[it_pob2] );
                bool bInFront( (p1 - p2) * n2 > 0 && (p2 - p1) * n1 > 0 );
                if( bInFront )
                {
                    StochasticSPC2::Pair np;
                    np.m_DistSq = min_dist_sq;
                    np.m_POB1 = min_dist_pob1;
                    np.m_POB2 = vec_pob2[it_pob2];
                    p_ds1->POB_IncRef( np.m_POB1 );
                    p_ds2->POB_IncRef( np.m_POB2 );
                    p_spc->m_vecNP.push_back( np );
                }
            }
#else
            StochasticSPC2::Pair np;
            np.m_DistSq = min_dist_sq;
            np.m_POB1 = min_dist_pob1;
            np.m_POB2 = vec_pob2[it_pob2];
            p_ds1->POB_IncRef( np.m_POB1 );
            p_ds2->POB_IncRef( np.m_POB2 );
            p_spc->m_vecNP.push_back( np );
#endif
        }
    }
}

//! Discard exceess pairs, no specific criteria
void SimplifyPairs_Discard( IDomainSampler2 *p_ds1,
                            IDomainSampler2 *p_ds2,
                            const Transform2 &tr12,
                            unsigned int max_pairs,
                            std::vector< StochasticSPC2::Pair > &vec_pair )
{
    if( vec_pair.size() <= max_pairs ) return;
    for( unsigned int it_pair=max_pairs; it_pair<vec_pair.size(); it_pair++ )
    {
        p_ds1->POB_DecRef( vec_pair[it_pair].m_POB1 );
        p_ds2->POB_DecRef( vec_pair[it_pair].m_POB2 );
    }
    vec_pair.resize(max_pairs);
}

//! Sort and keep <= max_pairs closest ones
void SimplifyPairs_KeepClosest( IDomainSampler2 *p_ds1,
                                IDomainSampler2 *p_ds2,
                                const Transform2 &tr12,
                                unsigned int max_pairs,
                                std::vector< StochasticSPC2::Pair > &vec_pair )
{
    if( vec_pair.size() <= max_pairs ) return;
    std::sort( vec_pair.begin(), vec_pair.end() );
#ifdef __ENABLE_TRACE_PAIRS
    GEO_LOG( "Discarding pairs [%f..%f,%f]", vec_pair.front().m_DistSq, vec_pair[max_pairs].m_DistSq, vec_pair.back().m_DistSq );
#endif
    for( unsigned int it_pair=max_pairs; it_pair<vec_pair.size(); it_pair++ )
    {
        p_ds1->POB_DecRef( vec_pair[it_pair].m_POB1 );
        p_ds2->POB_DecRef( vec_pair[it_pair].m_POB2 );
    }
    vec_pair.resize(max_pairs);
}

//! Sort, remove duplicated and keep <= max_pairs closest ones
template <typename PairT>
void GSimplifyPairs_RemoveDuplicatedAndKeepClosest( IDomainSampler2 *p_ds1,
                                                    IDomainSampler2 *p_ds2,
                                                    const Transform2 &tr12,
                                                    Real epsilon,
                                                    unsigned int max_pairs,
                                                    std::vector< PairT > &vec_pair )
{
    if( 0 == vec_pair.size() ) return;
    Real epsilon_sq( mal::Sq(epsilon) );
    // Sort by distance
    std::sort( vec_pair.begin(), vec_pair.end() );
    std::vector< PairT > vec_unique_pair;
    // Remove duplicates using distance-sort, stop if unique max_pairs is reached. First pair is always unique
    vec_unique_pair.push_back( vec_pair.front() );
    for( unsigned int it_pair=1; it_pair < vec_pair.size(); it_pair++ )
    {
        PairT &pair( vec_pair[it_pair] );
        // If there's still room for unique pairs,
        if( vec_unique_pair.size() < max_pairs )
        {
            PairT &unique_pair( vec_unique_pair.back() );
            // Check if pair is potentially redundant comparing its dist_sq to last unique pair
            if( pair.m_DistSq - unique_pair.m_DistSq > epsilon_sq )
                vec_unique_pair.push_back( pair ); //refcounts do not change as we remove from vec_pair and add to vec_unique_pair
            else
            {
                // Merge POB and pairs, if within epsilon_sq
                Real dist_sq1( mal::NormSq( p_ds1->POB_Position( pair.m_POB1 )
                                            - p_ds1->POB_Position( unique_pair.m_POB1 ) ) );
                Real dist_sq2( mal::NormSq( p_ds2->POB_Position( pair.m_POB2 )
                                            - p_ds2->POB_Position( unique_pair.m_POB2 ) ) );
                // Merge POB if close enough
                if( dist_sq1 < epsilon_sq )
                {
                    // Merge POB1
                    p_ds1->POB_DecRef( pair.m_POB1 );
                    p_ds1->POB_IncRef( unique_pair.m_POB1 );
                    pair.m_POB1 = unique_pair.m_POB1;
                }
                if( dist_sq2 < epsilon_sq )
                {
                    // Merge POB2
                    p_ds2->POB_DecRef( pair.m_POB2 );
                    p_ds2->POB_IncRef( unique_pair.m_POB2 );
                    pair.m_POB2 = unique_pair.m_POB2;
                }
                // Save if unique, with AT MOST ONE potentially merged POB
                if( dist_sq1 >= epsilon_sq || dist_sq2 >= epsilon_sq )
                    vec_unique_pair.push_back( pair ); //refcounts do not change as we remove from vec_pair and add to vec_unique_pair
                else
                {
                    //DecRef, as the pair is not unique
                    p_ds1->POB_DecRef( pair.m_POB1 );
                    p_ds2->POB_DecRef( pair.m_POB2 );
                }
                /*\todo We do NOT update pair distances as they'll only be
                  off by 2*epsilon_sq at most, and we won't (probably) use
                  them any more during this iteration anyway */
            }
        }
        else //DecRef if there's already too many unique pairs
        {
            p_ds1->POB_DecRef( pair.m_POB1 );
            p_ds2->POB_DecRef( pair.m_POB2 );
        }
    }
#ifdef __ENABLE_TRACE_PAIRS
    if( vec_unique_pair.size() != vec_pair.size() )
        GEO_LOG( "Discarded pairs [%f..%f,%f] (%d<%d)",
                 vec_unique_pair.front().m_DistSq, vec_unique_pair.back().m_DistSq, vec_pair.back().m_DistSq,
                 (int)vec_unique_pair.size(), (int)vec_pair.size() );
#endif
    // Keep only unique pairs
    std::swap( vec_pair, vec_unique_pair );
}

/*! Sort, remove duplicated and keep <= max_pairs closest ones
  - Remove keeps IM
  - IP with IM beyond max_pairs are Kept
*/
void SimplifyPairs_IP( IDomainSampler2 *p_ds1,
                       IDomainSampler2 *p_ds2,
                       const Transform2 &tr12,
                       Real epsilon,
                       unsigned int max_pairs,
                       std::vector< StochasticSPC2::IntersectionPair > &vec_pair )
{
    if( 0 == vec_pair.size() ) return;
    Real epsilon_sq( mal::Sq(epsilon) );
    // Sort by distance
    std::sort( vec_pair.begin(), vec_pair.end() );
    std::vector< StochasticSPC2::IntersectionPair > vec_unique_pair;
    // Remove duplicates using distance-sort, stop if unique max_pairs is reached. First pair is always unique
    vec_unique_pair.push_back( vec_pair.front() );
    unsigned int num_discarded(0);
    for( unsigned int it_pair=1; it_pair < vec_pair.size(); it_pair++ )
    {
        StochasticSPC2::IntersectionPair &pair( vec_pair[it_pair] );
        StochasticSPC2::IntersectionPair &unique_pair( vec_unique_pair.back() );
        bool bDiscard(false);
        // Check if pair is potentially redundant comparing its increasing dist_sq to last unique pair
        if( pair.m_DistSq - unique_pair.m_DistSq > epsilon_sq )
        {
            // Keep it if there's still room for unique pairs or it has a valid IM
            if( vec_unique_pair.size() < max_pairs || pair.HasIM() )
                vec_unique_pair.push_back( pair );
            else
                bDiscard = true;
        }
        else
        {
            // Merge POB and pairs, if within epsilon_sq
            Real dist_sq1( mal::NormSq( p_ds1->POB_Position( pair.m_POB1 )
                                        - p_ds1->POB_Position( unique_pair.m_POB1 ) ) );
            Real dist_sq2( mal::NormSq( p_ds2->POB_Position( pair.m_POB2 )
                                        - p_ds2->POB_Position( unique_pair.m_POB2 ) ) );
            // Merge POB if close enough
            if( dist_sq1 < epsilon_sq )
            {
                // Merge POB1
                p_ds1->POB_DecRef( pair.m_POB1 );
                p_ds1->POB_IncRef( unique_pair.m_POB1 );
                pair.m_POB1 = unique_pair.m_POB1;
            }
            if( dist_sq2 < epsilon_sq )
            {
                // Merge POB2
                p_ds2->POB_DecRef( pair.m_POB2 );
                p_ds2->POB_IncRef( unique_pair.m_POB2 );
                pair.m_POB2 = unique_pair.m_POB2;
            }
            // Save if unique, with AT MOST ONE potentially merged POB
            if( dist_sq1 >= epsilon_sq || dist_sq2 >= epsilon_sq )
            {
                // Keep it if there's still room for unique pairs or it has a valid IM
                if( vec_unique_pair.size() < max_pairs || pair.HasIM() )
                    vec_unique_pair.push_back( pair );
                else
                    bDiscard = true;
            }
            else
            {
                // Try to merge, it it's not possible, Keep both IP
                if( !unique_pair.MergeIM( pair, p_ds1, p_ds2 ) )
                    vec_unique_pair.push_back( pair );
                else
                    bDiscard = true;
            }
            /*\todo We do NOT update pair distances as they'll only be
              off by 2*epsilon_sq at most, and we won't (probably) use
              them any more during this iteration anyway */
        }
        if( bDiscard )
        {
            p_ds1->POB_DecRef( pair.m_POB1 );
            p_ds2->POB_DecRef( pair.m_POB2 );
            num_discarded++;
        }
    }

#ifdef __ENABLE_TRACE_PAIRS_IP
    if( vec_unique_pair.size() < vec_pair.size() )
    {
        if( vec_unique_pair.size() <= max_pairs )
        {
            GEO_LOG( "Discarded %d pairs [%f..%f,%f], kept %d < %d", num_discarded,
                     vec_unique_pair.front().m_DistSq, vec_unique_pair.back().m_DistSq, vec_pair.back().m_DistSq,
                     (int)vec_unique_pair.size(), (int)vec_pair.size() );
        }
        else
        {
            GEO_LOG( "Discarded %d pairs [%f..%f,%f], kept %d > %d due to unmerged IM", num_discarded,
                     vec_unique_pair.front().m_DistSq, vec_unique_pair.back().m_DistSq, vec_pair.back().m_DistSq,
                     (int)vec_unique_pair.size(), max_pairs );
        }
    }
#endif
    // Keep only unique pairs
    std::swap( vec_pair, vec_unique_pair );
}

Vec2 ComputeAvgPos( unsigned int num_points, const Vec2 *vec_points )
{
    Vec2 acc_pos( Vec2::Zero() );
    for( unsigned int i=0; i<num_points; i++ ) acc_pos += vec_points[i];
    if( num_points > 0 ) return acc_pos / Real(num_points);
    else return Vec2::Zero();
}

Mat2x2 ComputeCovarianceMatrix2( unsigned int num_points, const Vec2 *vec_points, const Vec2 &avg_pos )
{
    Mat2x2 acc_covariance( Mat2x2::Zero() );
    for( unsigned int i=0; i<num_points; i++ )
    {
        Vec2 d( vec_points[i] - avg_pos );
        acc_covariance += Mat2x2( mal::Sq(d[0]), d[0]*d[1],
                                  d[1]*d[0], mal::Sq(d[1]) );
    }
    if( num_points > 1 ) return acc_covariance * mal::Rcp( Real(num_points-1) );
    else return Mat2x2::Zero();
}

/*! \todo Should work on POD vectors and on sample-pair std::vector, choosing one or the other point POB1/POB2 and avoiding reevaluation
 */
bool ComputePCA( unsigned int num_points, const Vec2 *vec_points,
                 Vec2 &avg_position,
                 Real &eigen_value0, Vec2 &eigen_vector0 )
{
    // Compute average
    avg_position = ComputeAvgPos( num_points, vec_points );
    // Compute covariance matrix
    Mat2x2 covariance = ComputeCovarianceMatrix2( num_points, vec_points, avg_position );
    // Compute largest eigenstuff
    Real error = mal::GComputeLargestEigenvalue_Symmetric( covariance,
                                                           eigen_value0, eigen_vector0,
                                                           10, Real(0.001f) );
#ifdef __TRACE_PCA
    std::cout << "avg_position = " << avg_position << std::endl
              << " covariance = " << covariance << std::endl
              << " err " << error << std::endl
              << " e0 " << eigen_value0 << std::endl
              << " v0 " << eigen_vector0 << std::endl;
#endif

    return error >= 0;
}

}} //namespace geo::np
