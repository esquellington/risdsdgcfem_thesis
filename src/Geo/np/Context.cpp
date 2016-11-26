#include "Context.h"

namespace geo {
namespace np {

/*
static Context gs_DefaultContext;
Context *g_pDefaultContext( &gs_DefaultContext );
*/

#ifdef __GEO_ENABLE_PARAMS
void Context::InitArchetype( util::ArchetypeLibrary &al )
{
    al.BeginArchetype( "Archetype_geo_np_Context" );
    {
        Context context;
        al.BeginProperty_Group( "<Tolerances>" );
        {
            al.BeginProperty_Group( "<Common>" );
            {
                al.AddProperty_NIR<Real>( "Log10(eps.Len)", context.m_Epsilon_Length_Log10, -10.0f, -1.0f,
                                          archetype_offset_of(context,m_Epsilon_Length_Log10), IArchetypeInstance::NTPF_Rebuild );
                al.AddProperty_NIR<Real>( "Log10(eps.Dir)", context.m_Epsilon_Dir_Log10, -10.0f, -1.0f,
                                          archetype_offset_of(context,m_Epsilon_Dir_Log10), IArchetypeInstance::NTPF_Rebuild );
                // al.AddProperty_NIR<Real>( "eps.Length", context.m_Epsilon_Length, 0.0f, 0.1f,
                //                           archetype_offset_of(context,m_Epsilon_Length), IArchetypeInstance::NTPF_Rebuild );
                // al.AddProperty_NIR<Real>( "eps.Dir", context.m_Epsilon_Dir, 0.0f, 0.1f,
                //                           archetype_offset_of(context,m_Epsilon_Dir), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty_NIR<Real>( "eps.Depth", context.m_Epsilon_ContactDepth, 0.0f, 1.0f,
                                          archetype_offset_of(context,m_Epsilon_ContactDepth), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty_NIR<Real>( "eps.Proximity", context.m_Epsilon_ProximityDistance, 0.0f, 1.0f,
                                          archetype_offset_of(context,m_Epsilon_ProximityDistance), IArchetypeInstance::NTPF_Ignore );
            }
            al.EndProperty_Group();
        }
        al.EndProperty_Group();

        /*TEMP: Disabled to avoid clutter...
        al.BeginProperty_Group( "<Stochastic>" );
        {
            al.BeginProperty_Group( "<Develop>" );
                {
                    al.AddProperty( "Fallback", context.m_Stochastic_Fallback,
                                    archetype_offset_of(context,m_Stochastic_Fallback), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty( "Persistent", context.m_Stochastic_Persistent,
                                    archetype_offset_of(context,m_Stochastic_Persistent), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();
                // Thresholds
                al.BeginProperty_Group( "<Thresholds>" );
                {
                    al.AddProperty_NIR<Real>( "ZeroDist", context.m_Stochastic_ZeroDist, 0.0f, 0.1f,
                                              archetype_offset_of(context,m_Stochastic_ZeroDist), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<Real>( "NearDist", context.m_Stochastic_NearDist, 0.0f, 1.0f,
                                              archetype_offset_of(context,m_Stochastic_NearDist), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<Real>( "FarDist", context.m_Stochastic_FarDist, 0.0f, 10.0f,
                                              archetype_offset_of(context,m_Stochastic_FarDist), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();
                // Stochastic generation counts
                al.BeginProperty_Group( "<RandomSearch>" );
                {
                    al.AddProperty_NIR<uint32>( "#IRP", context.m_Stochastic_Num_IRP, 1, 1000,
                                                archetype_offset_of(context,m_Stochastic_Num_IRP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "#RP", context.m_Stochastic_Num_RP, 1, 100,
                                                archetype_offset_of(context,m_Stochastic_Num_RP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "#IterIP", context.m_Stochastic_Refine_Num_Iter_IP, 1, 10,
                                                archetype_offset_of(context,m_Stochastic_Refine_Num_Iter_IP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "#RNIP", context.m_Stochastic_Refine_Num_Neighbours_IP, 1, 100,
                                                archetype_offset_of(context,m_Stochastic_Refine_Num_Neighbours_IP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "#RNCP", context.m_Stochastic_Refine_Num_Neighbours_CP, 1, 100,
                                                archetype_offset_of(context,m_Stochastic_Refine_Num_Neighbours_CP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "#RNNP", context.m_Stochastic_Refine_Num_Neighbours_NP, 1, 100,
                                                archetype_offset_of(context,m_Stochastic_Refine_Num_Neighbours_NP), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();
                // Persistence
                al.BeginProperty_Group( "<Persistence>" );
                {
                    // Limits
                    al.AddProperty_NIR<uint32>( "max.IP", context.m_Stochastic_Max_IP, 1, 100,
                                                archetype_offset_of(context,m_Stochastic_Max_IP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "max.CP", context.m_Stochastic_Max_CP, 1, 100,
                                                archetype_offset_of(context,m_Stochastic_Max_CP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "max.NP", context.m_Stochastic_Max_NP, 1, 1000,
                                                archetype_offset_of(context,m_Stochastic_Max_NP), IArchetypeInstance::NTPF_Ignore );
                    // Simplify
                    al.AddProperty_NIR<Real>( "reps.SimplifyIP", context.m_Stochastic_SimplifyPairs_EpsilonRel_IP, 0.0f, 1.0f,
                                              archetype_offset_of(context,m_Stochastic_SimplifyPairs_EpsilonRel_IP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<Real>( "reps.SimplifyCP", context.m_Stochastic_SimplifyPairs_EpsilonRel_CP, 0.0f, 1.0f,
                                              archetype_offset_of(context,m_Stochastic_SimplifyPairs_EpsilonRel_CP), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<Real>( "reps.SimplifyNP", context.m_Stochastic_SimplifyPairs_EpsilonRel_NP, 0.0f, 1.0f,
                                              archetype_offset_of(context,m_Stochastic_SimplifyPairs_EpsilonRel_NP), IArchetypeInstance::NTPF_Ignore );
                }
            al.EndProperty_Group();
        }
        al.EndProperty_Group();
        */

        /*TEMP: Disable to avoid clutter
        al.BeginProperty_Group( "<IM>" );
        {
            al.AddProperty_NIR<Real>( "StepLength", context.m_IM_StepLength, 0.01f, 1.0f,
                                      archetype_offset_of(context,m_IM_StepLength), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR<Real>( "reps.StepLength", context.m_IM_StepLength_Epsilon_Rel, 0.0f, 1.0f,
                                      archetype_offset_of(context,m_IM_StepLength_Epsilon_Rel), IArchetypeInstance::NTPF_Rebuild );

            //\todo WARNING: Cannot serialize Flags32 by now
            const char *vec_names[] = { "Rev", "Red", "MinD", "MinDA", "Proj", "Fix" };
            int32 vec_values[] = { (int32)eIMF_Reverse,
                                   (int32)eIMF_Reduce,
                                   (int32)eIMF_MinDist,
                                   (int32)eIMF_MinDistAngle,
                                   (int32)eIMF_Project,
                                   (int32)eIMF_Fix };
            al.AddProperty_Flags32( "IM_Flags", context.m_IM_Flags,
                                    6, vec_names, vec_values,
                                    archetype_offset_of(context,m_IM_Flags),
                                    IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
        */

        al.BeginProperty_Group( "<SS2SS>" );
        {
            const char *vec_names_m2m[] = { "BFSS", "StochSS", "BFDCR", "BVHDCR", "BFSS_LDCR", "BVHSS_LDCR" };
            uint32 vec_values_m2m[] = { eSS2SS_Bruteforce_SS,
                                        eSS2SS_Stochastic_SS,
                                        eSS2SS_Bruteforce_DCR,
                                        eSS2SS_BVH_DCR,
                                        eSS2SS_Bruteforce_SS_Lazy_DCR,
                                        eSS2SS_BVH_SS_Lazy_DCR };
            al.AddProperty_Enum32( "Method", (uint32)context.m_SS2SS_Method,
                                   cNumSS2SS, vec_names_m2m, vec_values_m2m,
                                   archetype_offset_of(context,m_SS2SS_Method),
                                   IArchetypeInstance::NTPF_Ignore );
            const char *vec_names_bvhm[] = { "BV(E)", "BV(BSlab)", "BV(BDOP)", "!Refit" };
            uint32 vec_values_bvhm[] = { eBVHM_BV_E,
                                         eBVHM_BV_BSlab,
                                         eBVHM_BV_BDOP,
                                         eBVHM_BV_NoRefit };
            al.AddProperty_Enum32( "BVHM", (uint32)context.m_BVH_Method,
                                   cNumBVHM, vec_names_bvhm, vec_values_bvhm,
                                   archetype_offset_of(context,m_BVH_Method),
                                   IArchetypeInstance::NTPF_Ignore );
            const char *vec_names_tbvhm[] = { "BV2BV", "BDOP2GSlabV" };
            uint32 vec_values_tbvhm[] = { eTBVHM_BV_Vs_BV,
                                          eTBVHM_BDOP_Vs_GSlabVertices };
            al.AddProperty_Enum32( "TBVHM", (uint32)context.m_TestBVH_Method,
                                   cNumTBVHM, vec_names_tbvhm, vec_values_tbvhm,
                                   archetype_offset_of(context,m_TestBVH_Method),
                                   IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<DCR2DCR>" );
        {
            al.AddProperty( "Log?", context.m_DCR2DCR_Log_Enabled,
                            archetype_offset_of(context,m_DCR2DCR_Log_Enabled), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "Viz?", context.m_DCR2DCR_Viz_Enabled,
                            archetype_offset_of(context,m_DCR2DCR_Viz_Enabled), IArchetypeInstance::NTPF_Ignore );
            const char *vec_names_cfpm[] = { "BF", "BDT" };
            uint32 vec_values_cfpm[] = { eDCR2DCR_CFP_BruteForce,
                                         eDCR2DCR_CFP_BDT };
            al.AddProperty_Enum32( "CFP.Method", (uint32)context.m_DCR2DCR_CFP_Method,
                                   cNumDCR2DCR_CFP, vec_names_cfpm, vec_values_cfpm,
                                   archetype_offset_of(context,m_DCR2DCR_CFP_Method),
                                   IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "IC?", context.m_DCR2DCR_IC_Enabled,
                            archetype_offset_of(context,m_DCR2DCR_IC_Enabled), IArchetypeInstance::NTPF_Ignore );

            const char *vec_names_ibm[] = { "None", "FloodT", "FloodP" };
            uint32 vec_values_ibm[] = { eDCR2DCR_IB_None,
                                        eDCR2DCR_IB_FloodT,
                                        eDCR2DCR_IB_FloodP };
            al.AddProperty_Enum32( "IB.Method", (uint32)context.m_DCR2DCR_IB_Method,
                                   cNumDCR2DCR_IB, vec_names_ibm, vec_values_ibm,
                                   archetype_offset_of(context,m_DCR2DCR_IB_Method),
                                   IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "IB.CT?", context.m_DCR2DCR_IB_EnableCT,
                            archetype_offset_of(context,m_DCR2DCR_IB_EnableCT), IArchetypeInstance::NTPF_Ignore );
            al.BeginProperty_Group( "<CFP>" );
            {
                al.AddProperty_NIR( "log1(CFP.EpsIC)", context.m_DCR2DCR_CFP_EpsilonIC_Log10,
                                    -10.0f, -1.0f,
                                    archetype_offset_of(context,m_DCR2DCR_CFP_EpsilonIC_Log10),
                                    IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "E2E.BDT.PreTr", context.m_DCR2DCR_E2E_BDT_Pretransform,
                                archetype_offset_of(context,m_DCR2DCR_E2E_BDT_Pretransform), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "E2E.BDT.CacheTr", context.m_DCR2DCR_E2E_BDT_CacheTransformedV,
                                archetype_offset_of(context,m_DCR2DCR_E2E_BDT_CacheTransformedV), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "E2E.BDT.Presort", context.m_DCR2DCR_E2E_BDT_Presort,
                                archetype_offset_of(context,m_DCR2DCR_E2E_BDT_Presort), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "E2E.BDT.Unpacked", context.m_DCR2DCR_E2E_BDT_Unpacked,
                                archetype_offset_of(context,m_DCR2DCR_E2E_BDT_Unpacked), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "E2E.BDT.DscndLrgr", context.m_DCR2DCR_E2E_BDT_DescendLarger,
                                archetype_offset_of(context,m_DCR2DCR_E2E_BDT_DescendLarger), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "E2E.BDT.Split3", context.m_DCR2DCR_E2E_BDT_Split3,
                                archetype_offset_of(context,m_DCR2DCR_E2E_BDT_Split3), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty_NIR<uint32>( "E2E.BDT.MaxLeafTests", context.m_DCR2DCR_E2E_BDT_MaxLeafTests,
                                            1, 32,
                                            archetype_offset_of(context,m_DCR2DCR_E2E_BDT_MaxLeafTests),
                                            IArchetypeInstance::NTPF_Ignore );
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<CP Gen>" );
            {
                const char *vec_names_h[] = { "d^2", "d^2a^2", "ff", "aff" };
                uint32 vec_values_h[] = { eDCR2DCR_H_DistSq,
                                          eDCR2DCR_H_DistSq_AngleSq,
                                          eDCR2DCR_H_FormFactor,
                                          eDCR2DCR_H_Affinity };
                al.AddProperty_Enum32( "H.Method", (uint32)context.m_DCR2DCR_H_Method,
                                       cNumDCR2DCR_H, vec_names_h, vec_values_h,
                                       archetype_offset_of(context,m_DCR2DCR_H_Method),
                                       IArchetypeInstance::NTPF_Ignore );
                const char *vec_names_im[] = { "Basic", "Projection", "Raycast" };
                uint32 vec_values_im[] = { eDCR2DCR_IM_Basic,
                                           eDCR2DCR_IM_Projection,
                                           eDCR2DCR_IM_Raycast };
                al.AddProperty_Enum32( "IM.Method", (uint32)context.m_DCR2DCR_IM_Method,
                                       cNumDCR2DCR_IM, vec_names_im, vec_values_im,
                                       archetype_offset_of(context,m_DCR2DCR_IM_Method),
                                       IArchetypeInstance::NTPF_Ignore );
                const char *vec_names_rc[] = { "None", "Global", "Direct", "Inverse", "Average" };
                uint32 vec_values_rc[] = { eDCR2DCR_RC_None,
                                           eDCR2DCR_RC_Global,
                                           eDCR2DCR_RC_Direct,
                                           eDCR2DCR_RC_Average,
                                           eDCR2DCR_RC_Inverse };
                al.AddProperty_Enum32( "RC.Method", (uint32)context.m_DCR2DCR_RC_Method,
                                       cNumDCR2DCR_RC, vec_names_rc, vec_values_rc,
                                       archetype_offset_of(context,m_DCR2DCR_RC_Method),
                                       IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "RC.OnlyIB", context.m_DCR2DCR_RC_UseOnlyIB,
                                archetype_offset_of(context,m_DCR2DCR_RC_UseOnlyIB), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "RC.ClosestOrGlobal", context.m_DCR2DCR_RC_UseClosestOrGlobal,
                                archetype_offset_of(context,m_DCR2DCR_RC_UseClosestOrGlobal), IArchetypeInstance::NTPF_Ignore );
                al.AddProperty( "RC.BasicIfNoHit", context.m_DCR2DCR_RC_UseBasicIfNoHit,
                                archetype_offset_of(context,m_DCR2DCR_RC_UseBasicIfNoHit), IArchetypeInstance::NTPF_Ignore );

                al.AddProperty_NIR( "RC.Thickness", context.m_DCR2DCR_RC_Thickness,
                                    0.0f, 0.1f,
                                    archetype_offset_of(context,m_DCR2DCR_RC_Thickness),
                                    IArchetypeInstance::NTPF_Ignore );
                al.AddProperty_NIR( "RC.FarHitDist", context.m_DCR2DCR_RC_FarHitThreshold,
                                    0.0f, 0.5f,
                                    archetype_offset_of(context,m_DCR2DCR_RC_FarHitThreshold),
                                    IArchetypeInstance::NTPF_Ignore );
                al.AddProperty_NIR( "RC.MinAngle_Deg", context.m_DCR2DCR_RC_MinAngle_Deg,
                                    0.0f, 45.0f,
                                    archetype_offset_of(context,m_DCR2DCR_RC_MinAngle_Deg),
                                    IArchetypeInstance::NTPF_Ignore );

                const char *vec_names_rcdcrm[] = { "BF", "BDT" };
                uint32 vec_values_rcdcrm[] = { eRCDCRM_BruteForce,
                                               eRCDCRM_BDT };
                al.AddProperty_Enum32( "RCDCR.Method", (uint32)context.m_RCDCR_Method,
                                       cNumRCDCRM, vec_names_rcdcrm, vec_values_rcdcrm,
                                       archetype_offset_of(context,m_RCDCR_Method),
                                       IArchetypeInstance::NTPF_Ignore );
                al.AddProperty_NIR<uint32>( "RCDCR.BDT.MaxLeafSize", context.m_RCDCR_BDT_MaxLeafSize,
                                            1, 16,
                                            archetype_offset_of(context,m_RCDCR_BDT_MaxLeafSize),
                                            IArchetypeInstance::NTPF_Ignore );
            }
            al.EndProperty_Group();
        }
        al.EndProperty_Group();
    }
    al.EndArchetype();

    /*TEMPORAL: Testing property groups:
    // IMPORTANT: This requires parseable property names (no spaces
    // inside identifiers...) and complains aboud <GROUP> item names
    // and Flags32 properties but seems to work anyway...

    util::ItemStream tmp_is(1<<15,1<<15);
    Context tmp_context;
    UTIL_LOG_WARNING("TEMPORAL: Saving...");
    al.SaveInstance( "Archetype_geo_np_Context", &tmp_context, tmp_is );
    UTIL_LOG_WARNING("TEMPORAL: ...Saved!");
    tmp_is.SaveTxt("tmp_Save.txt");


    UTIL_LOG_WARNING("TEMPORAL: Loading...");
    tmp_is.Clear();
    tmp_is.LoadTxt("tmp_Save.txt");
    tmp_context.m_IM_StepLength = 0;
    al.LoadInstance( "Archetype_geo_np_Context", &tmp_context, tmp_is.Begin() );
    tmp_context.m_IM_StepLength = 10;
    UTIL_LOG_WARNING("TEMPORAL: ...Loaded!");

    tmp_is.Clear();
    al.SaveInstance( "Archetype_geo_np_Context", &tmp_context, tmp_is );
    tmp_is.SaveTxt("tmp_Load.txt");
    */
}

bool Context::Rebuild()
{
    // <Common>
    m_Epsilon_Length = mal::Exp10( m_Epsilon_Length_Log10 );
    m_Epsilon_LengthSq = mal::Sq( m_Epsilon_Length );
    m_Epsilon_Dir = mal::Exp10( m_Epsilon_Dir_Log10 );
    // <IntersectionMapping>
    m_IM_StepLength_Epsilon_Abs = m_IM_StepLength_Epsilon_Rel * m_IM_StepLength;
    return true;
}

#endif //__GEO_ENABLE_PARAMS

}} //namespace geo::np
