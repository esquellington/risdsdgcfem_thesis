#ifndef GEO_NP_CONTEXT_H
#define GEO_NP_CONTEXT_H

#include <Geo/Config.h>

#ifdef  __GEO_ENABLE_PARAMS
#  include <util/Archetype.h>
#endif

#define __GEO_ENABLE_NP_SCRATCHPAD
#ifdef __GEO_ENABLE_NP_SCRATCHPAD
#  include <util/LinearAllocator.h>
#endif

namespace geo {
namespace np {

#ifdef __GEO_ENABLE_PARAMS
struct Context: public util::IArchetypeInstance
#else
struct Context
#endif
{
    enum EIMFlags {
        eIMF_None         = 0,
        eIMF_Reverse      = 1,
        eIMF_Reduce       = 2,
        eIMF_MinDist      = 4,
        eIMF_MinDistAngle = 8,
        eIMF_Project      = 16,
        eIMF_Fix          = 32,
        eIMF_Default      = ( eIMF_Reduce | eIMF_Project | eIMF_Fix )
    };

    Context()
    // Epsilons
    : m_Epsilon_Length( 1e-4 )
    , m_Epsilon_LengthSq( 1e-6 )
    , m_Epsilon_Dir( 1e-4 )
    , m_Epsilon_ContactDepth( 1e-5 ) //below this depth, no normalization can be done
    , m_Epsilon_ProximityDistance( mal::Epsilon<Real>() )
      // Log10 Eps
    , m_Epsilon_Length_Log10( mal::Log10(m_Epsilon_Length) )
    , m_Epsilon_Dir_Log10( mal::Log10(m_Epsilon_Dir) )
    // Stochastic
    , m_Stochastic_Fallback( true )
    , m_Stochastic_Persistent( true )
    , m_Stochastic_ZeroDist( 1e-3 )
    , m_Stochastic_NearDist( 0.1f )
    , m_Stochastic_FarDist( 0.2f )
    , m_Stochastic_Num_IRP(100)
    , m_Stochastic_Num_RP(20)
    , m_Stochastic_Refine_Num_Iter_IP( 1 )
    , m_Stochastic_Refine_Num_Neighbours_IP( 5 )
    , m_Stochastic_Refine_Num_Neighbours_CP( 5 )
    , m_Stochastic_Refine_Num_Neighbours_NP( 5 )
    , m_Stochastic_Max_IP(10)
    , m_Stochastic_Max_CP(10)
    , m_Stochastic_Max_NP(100)
    , m_Stochastic_SimplifyPairs_EpsilonRel_IP( 0.5f )
    , m_Stochastic_SimplifyPairs_EpsilonRel_CP( 0.1f )
    , m_Stochastic_SimplifyPairs_EpsilonRel_NP( 0.5f )
    // IntersectionMapping
    , m_IM_StepLength( 0.05f )
    , m_IM_StepLength_Epsilon_Rel( 0.1f )
    , m_IM_StepLength_Epsilon_Abs( 0.005f )
    , m_IM_Flags( eIMF_Default )
    , m_SS2SS_Method( eSS2SS_Default )
    , m_BVH_Method( eBVHM_BV_Default )
    , m_TestBVH_Method( eTBVHM_Default )
    , m_DCR2DCR_Log_Enabled(true)
    , m_DCR2DCR_Viz_Enabled(true)
    , m_DCR2DCR_CFP_Method( eDCR2DCR_CFP_Default )
    , m_DCR2DCR_IC_Enabled(true)
    , m_DCR2DCR_IB_Method( eDCR2DCR_IB_Default )
    , m_DCR2DCR_IB_EnableCT(true)
    , m_DCR2DCR_CFP_EpsilonIC_Log10(-5)
    , m_DCR2DCR_E2E_BDT_Pretransform(true)
    , m_DCR2DCR_E2E_BDT_CacheTransformedV(true)
    , m_DCR2DCR_E2E_BDT_Presort(true)
    , m_DCR2DCR_E2E_BDT_Unpacked(false)
    , m_DCR2DCR_E2E_BDT_DescendLarger(true)
    , m_DCR2DCR_E2E_BDT_Split3(true)
    , m_DCR2DCR_E2E_BDT_MaxLeafTests(8)
    , m_DCR2DCR_H_Method( eDCR2DCR_H_Default )
    , m_DCR2DCR_IM_Method( eDCR2DCR_IM_Default )
    , m_DCR2DCR_RC_Method( eDCR2DCR_RC_Default )
    , m_DCR2DCR_RC_UseOnlyIB(false)
    , m_DCR2DCR_RC_UseClosestOrGlobal(false)
    , m_DCR2DCR_RC_UseBasicIfNoHit(false)
    , m_DCR2DCR_RC_Thickness(0.025f) //TEMPORAL: This IS NOT an actual thickness but an OFFSET, just testing...
    , m_DCR2DCR_RC_FarHitThreshold(0.25f)
    , m_DCR2DCR_RC_MinAngle_Deg(0)
    , m_RCDCR_Method( eRCDCRM_Default )
    , m_RCDCR_BDT_MaxLeafSize(1)
#ifdef __GEO_ENABLE_NP_SCRATCHPAD
    , m_ScratchPad( 1<<20 )
#endif
        {}

    //! \name Params
    //@{
    //-- Epsilons
    Real m_Epsilon_Length;
    Real m_Epsilon_LengthSq; //derived
    Real m_Epsilon_Dir;
    Real m_Epsilon_ContactDepth;
    Real m_Epsilon_ProximityDistance;
    // Derived Log10 eps
    Real m_Epsilon_Length_Log10;
    Real m_Epsilon_LengthSq_Log10;
    Real m_Epsilon_Dir_Log10;
    //-- Stochastic
    // Develop
    bool m_Stochastic_Fallback; //Use deterministic fallback algorithm
    bool m_Stochastic_Persistent;
    // Thresholds
    Real m_Stochastic_ZeroDist; //Distance at which a NP becomes an IP
    Real m_Stochastic_NearDist;
    Real m_Stochastic_FarDist; //In the "Monte-Carlo collision detection" paper, m_FarDist = k * d_min + epsilon (k = 1.2, epsilon = 0.01*object_size)
    // Random search
    uint32 m_Stochastic_Num_IRP;
    uint32 m_Stochastic_Num_RP;
    uint32 m_Stochastic_Refine_Num_Iter_IP;
    uint32 m_Stochastic_Refine_Num_Neighbours_IP;
    uint32 m_Stochastic_Refine_Num_Neighbours_CP;
    uint32 m_Stochastic_Refine_Num_Neighbours_NP;
    // Persistence/Limits
    uint32 m_Stochastic_Max_IP;
    uint32 m_Stochastic_Max_CP;
    uint32 m_Stochastic_Max_NP;
    Real m_Stochastic_SimplifyPairs_EpsilonRel_IP;
    Real m_Stochastic_SimplifyPairs_EpsilonRel_CP;
    Real m_Stochastic_SimplifyPairs_EpsilonRel_NP;
    //-- IntersectionMapping
    Real m_IM_StepLength;
    Real m_IM_StepLength_Epsilon_Rel;
    Real m_IM_StepLength_Epsilon_Abs; //derived
    Flags32 m_IM_Flags;
    //-- SS2SS (SolidShape2SolidShape)
    enum ESS2SSMethod {
        //SS-only
        eSS2SS_Bruteforce_SS = 0,
        eSS2SS_Stochastic_SS,
        //DCR-only
        eSS2SS_Bruteforce_DCR,
        eSS2SS_BVH_DCR,
        //SS+DCR
        eSS2SS_Bruteforce_SS_Lazy_DCR,
        eSS2SS_BVH_SS_Lazy_DCR,
        //Common
        cNumSS2SS,
        eSS2SS_Default = eSS2SS_BVH_SS_Lazy_DCR
    } m_SS2SS_Method;
    enum EBVHMethod {
        eBVHM_BV_E = 0,
        eBVHM_BV_BSlab,
        eBVHM_BV_BDOP,
        eBVHM_BV_NoRefit,
        //Common
        cNumBVHM,
        eBVHM_BV_Default = eBVHM_BV_BDOP
    } m_BVH_Method;
    enum ETestBVHMethod {
        eTBVHM_BV_Vs_BV = 0,
        eTBVHM_BDOP_Vs_GSlabVertices,
        //Common
        cNumTBVHM,
        eTBVHM_Default = eTBVHM_BDOP_Vs_GSlabVertices
    } m_TestBVH_Method;
    //-- DCR2DCR
    bool m_DCR2DCR_Log_Enabled;
    bool m_DCR2DCR_Viz_Enabled;
    enum EDCR2DCR_CFP_Method {
        eDCR2DCR_CFP_BruteForce = 0,
        eDCR2DCR_CFP_BDT,
        //Common
        cNumDCR2DCR_CFP,
        eDCR2DCR_CFP_Default = eDCR2DCR_CFP_BDT
    } m_DCR2DCR_CFP_Method;
    bool m_DCR2DCR_IC_Enabled;
    enum EDCR2DCR_IB_Method {
        eDCR2DCR_IB_None = 0,
        eDCR2DCR_IB_FloodT,
        eDCR2DCR_IB_FloodP,
        //Common
        cNumDCR2DCR_IB,
        eDCR2DCR_IB_Default = eDCR2DCR_IB_None
    } m_DCR2DCR_IB_Method;
    bool m_DCR2DCR_IB_EnableCT;
    Real m_DCR2DCR_CFP_EpsilonIC_Log10;
    bool m_DCR2DCR_E2E_BDT_Pretransform;
    bool m_DCR2DCR_E2E_BDT_CacheTransformedV;
    bool m_DCR2DCR_E2E_BDT_Presort; //IMPORTANT: Presort is NOTICEABLY FASTER
    bool m_DCR2DCR_E2E_BDT_Unpacked;
    bool m_DCR2DCR_E2E_BDT_DescendLarger;
    bool m_DCR2DCR_E2E_BDT_Split3; //IMPORTANT: Split3 is A LOT FASTER
    uint32 m_DCR2DCR_E2E_BDT_MaxLeafTests;
    enum EDCR2DCR_H_Method {
        eDCR2DCR_H_DistSq = 0,
        eDCR2DCR_H_DistSq_AngleSq,
        eDCR2DCR_H_FormFactor,
        eDCR2DCR_H_Affinity,
        //Common
        cNumDCR2DCR_H,
        eDCR2DCR_H_Default = eDCR2DCR_H_DistSq
    } m_DCR2DCR_H_Method;
    enum EDCR2DCR_IM_Method {
        eDCR2DCR_IM_Basic = 0,
        eDCR2DCR_IM_Projection,
        eDCR2DCR_IM_Raycast,
        //Common
        cNumDCR2DCR_IM,
        eDCR2DCR_IM_Default = eDCR2DCR_IM_Raycast
    } m_DCR2DCR_IM_Method;
    enum EDCR2DCR_RC_Method {
        eDCR2DCR_RC_None = 0,
        eDCR2DCR_RC_Global,
        eDCR2DCR_RC_Direct,
        eDCR2DCR_RC_Inverse,
        eDCR2DCR_RC_Average,
        //Common
        cNumDCR2DCR_RC,
        eDCR2DCR_RC_Default = eDCR2DCR_RC_Global
    } m_DCR2DCR_RC_Method;

    bool m_DCR2DCR_RC_UseOnlyIB;
    bool m_DCR2DCR_RC_UseClosestOrGlobal;
    bool m_DCR2DCR_RC_UseBasicIfNoHit;
    Real m_DCR2DCR_RC_Thickness;
    Real m_DCR2DCR_RC_FarHitThreshold;
    Real m_DCR2DCR_RC_MinAngle_Deg;

    enum ERayCastDCR_Method {
        eRCDCRM_BruteForce = 0,
        eRCDCRM_BDT,
        //Common
        cNumRCDCRM,
        eRCDCRM_Default = eRCDCRM_BDT
    } m_RCDCR_Method;
    uint32 m_RCDCR_BDT_MaxLeafSize;
    //@}

    /* \todo ANY other potentially global stuff for all Test_X_Y() calls should be here
       ex:
       MemBlock m_ScratchPad;
       GPU/OpenCL resource handles...

       \todo Expose for interactive tweaking!
    */
#ifdef __GEO_ENABLE_PARAMS
    bool Rebuild();
    void SetName( const char *name ) {}
    const char *GetName() const { return "default_geo_np_Context"; } //TEMP: Cannot LOAD names with spaces!!
    static void InitArchetype( util::ArchetypeLibrary &al );
#endif

#ifdef __GEO_ENABLE_NP_SCRATCHPAD
    mutable util::LinearAllocator m_ScratchPad; //\todo Mutable is UGLY, but required for const g_pDefaultContext... this may change
#endif
};

extern Context *g_pDefaultContext;

}} //namespace geo::np

#endif // GEO_NP_CONTEXT_H
