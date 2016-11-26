#ifndef GEO_NP_STATS_H
#define GEO_NP_STATS_H

#include <Geo/Config.h>

#ifdef __GEO_ENABLE_STATS
#include <util/Chronos.h>

namespace geo {
namespace np {

struct Stats
{
    Stats() { Reset(); }
    inline void Reset() { for( unsigned int i=0; i<sizeof(Stats) ; i++ ) reinterpret_cast<char*>(this)[i] = 0; }

    struct StochasticStats
    {
        // Cost
        uint32 m_Num_NP;
        uint32 m_Num_IP;
        uint32 m_Num_POB;
        uint32 m_Num_Alloc_POB;
        // Efficiency
        uint32 m_Num_FalseIP; //IP in close proximity but no crossingg
        float32 m_Ratio_NP_div_RP;
        uint32 m_AvgLifetime_NP;
        // Correctness
        int32 m_Num_UDCM; //Undetected Disjoint Contact Manifolds
        //-- IntersectionMapping
        //float32 m_IM_RelMaxError_Depth;
        //float32 m_IM_RelAvgError_Depth;
    } stochastic;

    /*
    struct BruteForceStats
    {
        ...
    } bruteforce;
    */

    struct RayCastStats
    {
        uint32 m_GPlane;
        uint32 m_GHalfSpace;
        uint32 m_GSphere;
        uint32 m_GCenteredAABB;
        uint32 m_GBox;
        uint32 m_Cylinder2;
        uint32 m_GCapsule;
        uint32 m_Segment2_SS;
        uint32 m_Segment2_DS;
        uint32 m_Triangle3_SS;
        uint32 m_Triangle3_DS;
        // uint32 m_Triangle3_DS_C; //\todo No need to detail it, use m_Triangle3_DS
        uint32 m_TriSurface3;
    } raycast;

    struct IntersectionStats
    {
        uint32 m_Segment2_Segment2;
        uint32 m_Triangle3_Triangle3;
    } intersection;

    struct OverlapStats
    {
        uint32 m_Point2_Triangle2;
        uint32 m_Point2_Polygonal2;
        uint32 m_Point3_Tetrahedron3;
        uint32 m_Point3_TriSurface3;
        uint32 m_Segment2_Triangle2;
        //uint32 m_Segment2_Triangle2_BDOP; \todo No need to detail it, use m_Segment2_Triangle2
        uint32 m_Segment3_Tetrahedron3;
        uint32 m_Triangle3_Triangle3;
        uint32 m_Triangle3_Tetrahedron3;
        uint32 m_Triangle3_TriSurface3;
    } overlap;

    struct ContactStats
    {
        uint32 m_Num_Vertices_Transformed; //Normal Tr * v
        uint32 m_Num_Vertices_Transformed_Barycentric; //Barycentric transform
        util::Chronos::time_type m_Profile_MSS2_MSS2_ms; //\todo Consider "milliseconds" type in PropertyTreeControl to display profiling in its relevant scale, seconds are WAY TOO COARSE!!
        util::Chronos::time_type m_Profile_TSS3_TSS3_ms; //\todo Consider "milliseconds" type in PropertyTreeControl to display profiling in its relevant scale, seconds are WAY TOO COARSE!!
    } contact;

    struct MIG2015
    {
        util::Chronos::time_type m_Profile_TSS3_TSS3_ms;
        util::Chronos::time_type m_Profile_RefitBVH_ms;
        util::Chronos::time_type m_Profile_OverlapBVH_ms;
        util::Chronos::time_type m_Profile_FindITP_ms;
        util::Chronos::time_type m_Profile_MergeIC_ms;
        util::Chronos::time_type m_Profile_FindCP_ms;
        util::Chronos::time_type m_Profile_FloodIB_CT_ms;
        util::Chronos::time_type m_Profile_FloodIB_IT_ms;
        util::Chronos::time_type m_Profile_RaycastBVH_ms;
        util::Chronos::time_type m_Profile_RaycastDCR_ms;
        util::Chronos::time_type m_Profile_TSS3_Plane_ms;

        util::Chronos::time_type m_Profile_TSS3_TSS3_ms_ACC;
        util::Chronos::time_type m_Profile_RefitBVH_ms_ACC;
        util::Chronos::time_type m_Profile_OverlapBVH_ms_ACC;
        util::Chronos::time_type m_Profile_FindITP_ms_ACC;
        util::Chronos::time_type m_Profile_MergeIC_ms_ACC;
        util::Chronos::time_type m_Profile_FindCP_ms_ACC;
        util::Chronos::time_type m_Profile_FloodIB_CT_ms_ACC;
        util::Chronos::time_type m_Profile_FloodIB_IT_ms_ACC;
        util::Chronos::time_type m_Profile_RaycastBVH_ms_ACC;
        util::Chronos::time_type m_Profile_RaycastDCR_ms_ACC;
        util::Chronos::time_type m_Profile_TSS3_Plane_ms_ACC;

        uint32 m_Num_ITP;
        uint32 m_Num_IC;
        uint32 m_Num_IB_CT;
        uint32 m_Num_IB_FT;
        uint32 m_Num_IB_FP;
        uint32 m_Num_IB;
        uint32 m_Num_CP;

    } mig2015;
};
extern Stats *g_pStats;

}} //namespace geo::np

#  define GEO_NP_STAT_SET( name, value ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) = (value); } }
#  define GEO_NP_STAT_INC( name ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name)++; } }
#  define GEO_NP_STAT_DEC( name ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name)--; } }
#  define GEO_NP_STAT_ADD( name, value ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) += (value); } }
#  define GEO_NP_STAT_ACC( name, value ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) += (value); } } //\todo Same as ADD, but assumed to ACC during a whole timestep

#  define GEO_NP_BEGIN_PROFILE_BLOCK_MS( name ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) = util::Chronos::GetSysTime_ms(); } }
#  define GEO_NP_END_PROFILE_BLOCK_MS( name ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) = util::Chronos::GetSysTime_ms() - (geo::np::g_pStats->name); } }

#  define GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( name ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) = util::Chronos::GetSysTime_ms(); } }
#  define GEO_NP_END_PROFILE_BLOCK_ACC_MS( name ) { if(geo::np::g_pStats) { (geo::np::g_pStats->name) = util::Chronos::GetSysTime_ms() - (geo::np::g_pStats->name); (geo::np::g_pStats->name##_ACC) += (geo::np::g_pStats->name); } }

#else //__GEO_ENABLE_STATS

#  define GEO_NP_STAT_SET( name, value )
#  define GEO_NP_STAT_INC( name )
#  define GEO_NP_STAT_DEC( name )
#  define GEO_NP_STAT_ADD( name, value )
#  define GEO_NP_STAT_ACC( name, value )

#endif //__GEO_ENABLE_STATS

#endif // GEO_NP_STATS_H
