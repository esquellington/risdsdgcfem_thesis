#include "geo_stats.h"

#include <Geo/Config.h>

#ifdef __GEO_ENABLE_STATS
#  include <util/Archetype.h>
//#  include "bp/stats.h"
//#  include "mp/stats.h"
#  include "np/stats.h"
#  include "bv/stats.h"
#  include <string>
#endif

namespace geo
{

#ifdef __GEO_ENABLE_STATS
struct Stats
{
    Stats() {}
    void Reset() { np.Reset(); bv.Reset(); }
    //bp::Stats bp;
    //mp::Stats mp;
    np::Stats np;
    bv::Stats bv;
    //...
};
/* \todo Consider splitting IArchetypeInstance stuff into a subclass,
 so that accross np:: only Stats needs to be included, and the
 IArchetypeInstance is just an optional extension, defined in some
 .cpp file not included from ANYWHERE ELSE
*/
struct Stats_Archetype: public Stats, public util::IArchetypeInstance
{
    bool Rebuild() { return true; }
    void SetName( const char *name ) {}
    const char *GetName() const { return "default geo::Stats"; }
    static void InitArchetype( util::ArchetypeLibrary &al );
};
static util::ArchetypeLibrary s_StatsAL;
static Stats_Archetype s_Stats;
np::Stats *np::g_pStats = 0;
bv::Stats *bv::g_pStats = 0;

#endif //__GEO_ENABLE_STATS

bool Init_Stats()
{
#ifdef __GEO_ENABLE_STATS
    //\todo s_StatsAL.Init();
    Stats_Archetype::InitArchetype( s_StatsAL );
    //bp::g_pStats = &s_Stats.BP;
    //mp::g_pStats = &s_Stats.MP;
    np::g_pStats = &s_Stats.np;
    bv::g_pStats = &s_Stats.bv;
    return true;
#else
    return false;
#endif
}

bool ShutDown_Stats()
{
#ifdef __GEO_ENABLE_STATS
    //\todo remove archetype or destroy local library...
    //\todo s_StatsAL.Clear();
    //bp::g_pStats = 0;
    //mp::g_pStats = 0,
    np::g_pStats = 0;
    bv::g_pStats = 0;
    return true;
#else
    return false;
#endif
}

// Begin/End stats collection
void BeginStats()
{
    s_Stats.Reset();
}

void EndStats()
{
}

util::ItemStream::ItemItRW QueryStats( util::ItemStream &stats_is )
{
#ifdef __GEO_ENABLE_STATS
    stats_is.BeginComplex( "geo::Stats", eType_Property_Group );
    {
        s_StatsAL.ExportInstance( "Archetype_geo_Stats", &s_Stats, stats_is );
        //s_Stats.Reset();
    }
    return stats_is.EndComplex();
#else
    return util::ItemStream::ItemItRW();
#endif
}

void SyncStats( util::ItemStream::ItemItRW stats_it )
{
#ifdef __GEO_ENABLE_STATS
    GEO_ASSERT( std::string("geo::Stats") == stats_it.GetName() );
    s_StatsAL.SyncInstance( "Archetype_geo_Stats", &s_Stats, stats_it.GetSubItem().Find("default geo::Stats") );
    //s_Stats.Reset();
#endif
}

//----------------------------------------------------------------
void Stats_Archetype::InitArchetype( util::ArchetypeLibrary &al )
{
    Stats_Archetype sa;
    al.BeginArchetype( "Archetype_geo_Stats" );
    {
        //\todo <BP>...
        al.BeginProperty_Group( "<NP>" );
        {
            al.BeginProperty_Group( "<MIG2015 Prof>" );
            {
                al.AddProperty( "TSS3vsTSS3 (ms)", sa.np.mig2015.m_Profile_TSS3_TSS3_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_TSS3_TSS3_ms_ACC ) );
                al.AddProperty( "+RefitBVH", sa.np.mig2015.m_Profile_RefitBVH_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_RefitBVH_ms_ACC ) );
                al.AddProperty( "+OverlapBVH", sa.np.mig2015.m_Profile_OverlapBVH_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_OverlapBVH_ms_ACC ) );
                al.AddProperty( "+FindITP", sa.np.mig2015.m_Profile_FindITP_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_FindITP_ms_ACC ) );
                al.AddProperty( "+MergeIC", sa.np.mig2015.m_Profile_MergeIC_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_MergeIC_ms_ACC ) );
                al.AddProperty( "+FindCP", sa.np.mig2015.m_Profile_FindCP_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_FindCP_ms_ACC ) );
                al.AddProperty( "++FloodIB.CT", sa.np.mig2015.m_Profile_FloodIB_CT_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_FloodIB_CT_ms_ACC ) );
                al.AddProperty( "++FloodIB.IT", sa.np.mig2015.m_Profile_FloodIB_IT_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_FloodIB_IT_ms_ACC ) );
                al.AddProperty( "+RaycastBVH", sa.np.mig2015.m_Profile_RaycastBVH_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_RaycastBVH_ms_ACC ) );
                al.AddProperty( "+RaycastDCR", sa.np.mig2015.m_Profile_RaycastDCR_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_RaycastDCR_ms_ACC ) );
                // al.AddProperty( "TSS3vsPlane (ms)", sa.np.mig2015.m_Profile_TSS3_Plane_ms_ACC, archetype_offset_of( sa, np.mig2015.m_Profile_TSS3_Plane_ms_ACC ) );//\todo DISABLED because it changes MIG counter indices in .out files
            }
            al.EndProperty_Group();
            al.BeginProperty_Group( "<MIG2015 Count>" );
            {
                al.AddProperty( "#ITP", sa.np.mig2015.m_Num_ITP, archetype_offset_of( sa, np.mig2015.m_Num_ITP ) );
                al.AddProperty( "#IC", sa.np.mig2015.m_Num_IC, archetype_offset_of( sa, np.mig2015.m_Num_IC ) );
                al.AddProperty( "#IB_CT", sa.np.mig2015.m_Num_IB_CT, archetype_offset_of( sa, np.mig2015.m_Num_IB_CT ) );
                al.AddProperty( "#IB_FT", sa.np.mig2015.m_Num_IB_FT, archetype_offset_of( sa, np.mig2015.m_Num_IB_FT ) );
                al.AddProperty( "#IB_FP", sa.np.mig2015.m_Num_IB_FP, archetype_offset_of( sa, np.mig2015.m_Num_IB_FP ) );
                al.AddProperty( "#CP", sa.np.mig2015.m_Num_CP, archetype_offset_of( sa, np.mig2015.m_Num_CP ) );
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<RayCast>" );
            {
                al.AddProperty( "GPlane", sa.np.raycast.m_GPlane, archetype_offset_of( sa, np.raycast.m_GPlane ) );
                al.AddProperty( "GHalfSpace", sa.np.raycast.m_GHalfSpace, archetype_offset_of( sa, np.raycast.m_GHalfSpace ) );
                al.AddProperty( "GSphere", sa.np.raycast.m_GSphere, archetype_offset_of( sa, np.raycast.m_GSphere ) );
                al.AddProperty( "GCenteredAABB", sa.np.raycast.m_GCenteredAABB, archetype_offset_of( sa, np.raycast.m_GCenteredAABB ) );
                al.AddProperty( "GBox", sa.np.raycast.m_GBox, archetype_offset_of( sa, np.raycast.m_GBox ) );
                al.AddProperty( "Cyl2", sa.np.raycast.m_Cylinder2, archetype_offset_of( sa, np.raycast.m_Cylinder2 ) );
                al.AddProperty( "GCapsule", sa.np.raycast.m_GCapsule, archetype_offset_of( sa, np.raycast.m_GCapsule ) );
                al.AddProperty( "Seg2SS", sa.np.raycast.m_Segment2_SS, archetype_offset_of( sa, np.raycast.m_Segment2_SS ) );
                al.AddProperty( "Seg2DS", sa.np.raycast.m_Segment2_DS, archetype_offset_of( sa, np.raycast.m_Segment2_DS ) );
                al.AddProperty( "Tri3SS", sa.np.raycast.m_Triangle3_SS, archetype_offset_of( sa, np.raycast.m_Triangle3_SS ) );
                al.AddProperty( "Tri3DS", sa.np.raycast.m_Triangle3_DS, archetype_offset_of( sa, np.raycast.m_Triangle3_DS ) );
                // al.AddProperty( "Tri3DSC", sa.np.raycast.m_Triangle3_DS_C, archetype_offset_of( sa, np.raycast.m_Triangle3_DS_C ) );
                al.AddProperty( "TriSurf3", sa.np.raycast.m_TriSurface3, archetype_offset_of( sa, np.raycast.m_TriSurface3 ) );
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<Stochastic>" );
            {
                /*TEMP: Disabled to avoid clutter
                al.AddProperty( "#NP", sa.np.stochastic.m_Num_NP, archetype_offset_of( sa, np.stochastic.m_Num_NP ) );
                al.AddProperty( "#IP", sa.np.stochastic.m_Num_IP, archetype_offset_of( sa, np.stochastic.m_Num_IP ) );
                al.AddProperty( "#POB", sa.np.stochastic.m_Num_POB, archetype_offset_of( sa, np.stochastic.m_Num_POB ) );
                al.AddProperty( "#AllocPOB", sa.np.stochastic.m_Num_Alloc_POB, archetype_offset_of( sa, np.stochastic.m_Num_Alloc_POB ) );
                al.AddProperty( "#FalseIP", sa.np.stochastic.m_Num_FalseIP, archetype_offset_of( sa, np.stochastic.m_Num_FalseIP ) );
                al.AddProperty( "NP/RP", sa.np.stochastic.m_Ratio_NP_div_RP, archetype_offset_of( sa, np.stochastic.m_Ratio_NP_div_RP ) );
                al.AddProperty( "#UDCM", sa.np.stochastic.m_Num_UDCM, archetype_offset_of( sa, np.stochastic.m_Num_UDCM ) );
                */
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<Intersection>" );
            {
                al.AddProperty( "Seg2-Seg2", sa.np.intersection.m_Segment2_Segment2, archetype_offset_of( sa, np.intersection.m_Segment2_Segment2 ) );
                al.AddProperty( "Tri3-Tri3", sa.np.intersection.m_Triangle3_Triangle3, archetype_offset_of( sa, np.intersection.m_Triangle3_Triangle3 ) );
            }
            al.EndProperty_Group();
            //\todo BruteForce, GJK, SAT, etc...

            al.BeginProperty_Group( "<Overlap>" );
            {
                al.AddProperty( "Point2-Tri2", sa.np.overlap.m_Point2_Triangle2, archetype_offset_of( sa, np.overlap.m_Point2_Triangle2 ) );
                al.AddProperty( "Point2-Poly2", sa.np.overlap.m_Point2_Polygonal2, archetype_offset_of( sa, np.overlap.m_Point2_Polygonal2 ) );
                al.AddProperty( "Point3_Tet3", sa.np.overlap.m_Point3_Tetrahedron3, archetype_offset_of( sa, np.overlap.m_Point3_Tetrahedron3 ) );
                al.AddProperty( "Point3_TriSurf3", sa.np.overlap.m_Point3_TriSurface3, archetype_offset_of( sa, np.overlap.m_Point3_TriSurface3 ) );
                al.AddProperty( "Seg2-Tri2", sa.np.overlap.m_Segment2_Triangle2, archetype_offset_of( sa, np.overlap.m_Segment2_Triangle2 ) );
                al.AddProperty( "Seg3-Tet3", sa.np.overlap.m_Segment3_Tetrahedron3, archetype_offset_of( sa, np.overlap.m_Segment3_Tetrahedron3 ) );
                al.AddProperty( "Tri3-Tri3", sa.np.overlap.m_Triangle3_Triangle3, archetype_offset_of( sa, np.overlap.m_Triangle3_Triangle3 ) );
                al.AddProperty( "Tri3-Tet3", sa.np.overlap.m_Triangle3_Tetrahedron3, archetype_offset_of( sa, np.overlap.m_Triangle3_Tetrahedron3 ) );
                al.AddProperty( "Tri3-TriSurf3", sa.np.overlap.m_Triangle3_TriSurface3, archetype_offset_of( sa, np.overlap.m_Triangle3_TriSurface3 ) );
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<Contact>" );
            {
                //\todo uninteresting, by now... al.AddProperty( "#V Tr", sa.np.contact.m_Num_Vertices_Transformed, archetype_offset_of( sa, np.contact.m_Num_Vertices_Transformed ) );
                al.AddProperty( "#V TrB", sa.np.contact.m_Num_Vertices_Transformed_Barycentric, archetype_offset_of( sa, np.contact.m_Num_Vertices_Transformed_Barycentric ) );
                al.AddProperty( "T MSS2-MSS2 (ms)", sa.np.contact.m_Profile_MSS2_MSS2_ms, archetype_offset_of( sa, np.contact.m_Profile_MSS2_MSS2_ms ) );
                al.AddProperty( "T TSS3-TSS3 (ms)", sa.np.contact.m_Profile_TSS3_TSS3_ms, archetype_offset_of( sa, np.contact.m_Profile_TSS3_TSS3_ms ) );
            }
            al.EndProperty_Group();
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<BV>" );
        {
            al.BeginProperty_Group( "<Overlap>" );
            {
                al.AddProperty( "#GAABB-GAABB", sa.bv.testoverlap.m_GAABB_GAABB, archetype_offset_of( sa, bv.testoverlap.m_GAABB_GAABB ) );
                al.AddProperty( "#GSphere-GSphere", sa.bv.testoverlap.m_GSphere_GSphere, archetype_offset_of( sa, bv.testoverlap.m_GSphere_GSphere ) );
                al.AddProperty( "#DOP2K4-DOP2K4", sa.bv.testoverlap.m_DOP2K4_DOP2K4, archetype_offset_of( sa, bv.testoverlap.m_DOP2K4_DOP2K4 ) );
                al.AddProperty( "#DOP2K8-DOP2K8", sa.bv.testoverlap.m_DOP2K8_DOP2K8, archetype_offset_of( sa, bv.testoverlap.m_DOP2K8_DOP2K8 ) );
                al.AddProperty( "#DOP3K6-DOP3K6", sa.bv.testoverlap.m_DOP3K6_DOP3K6, archetype_offset_of( sa, bv.testoverlap.m_DOP3K6_DOP3K6 ) );
                al.AddProperty( "#DOP3K14-DOP3K14", sa.bv.testoverlap.m_DOP3K14_DOP3K14, archetype_offset_of( sa, bv.testoverlap.m_DOP3K14_DOP3K14 ) );
                al.AddProperty( "#GBDOP-GBDOP", sa.bv.testoverlap.m_GBDOP_GBDOP, archetype_offset_of( sa, bv.testoverlap.m_GBDOP_GBDOP ) );
                al.AddProperty( "#GBDOP-GBDOP-A", sa.bv.testoverlap.m_GBDOP_GBDOP_Axis, archetype_offset_of( sa, bv.testoverlap.m_GBDOP_GBDOP_Axis ) );
            }
            al.EndProperty_Group();
        }
        al.EndProperty_Group();
    }
    al.EndArchetype();
}

} //namespace geo
