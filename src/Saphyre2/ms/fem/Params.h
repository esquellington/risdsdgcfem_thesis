#ifndef S2_MS_FEM_PARAMS_H
#define S2_MS_FEM_PARAMS_H

#include "../Config.h"

namespace S2 {
namespace ms {
namespace fem {

struct Params
{
    // Material
    Real m_YoungModulus;
    Real m_PoissonRatio;
    Real m_DampingRatio;
    // Rayleigh
    Real m_RayleighFreq1;
    Real m_RayleighFreq2;
    Vec2 m_RayleighCoeff;
    // Plasticity
    Real m_PlasticYield;
    Real m_PlasticMax;
    Real m_PlasticCreepPerSecond;
    // Misc
    Real m_FactorDetF;
    Real m_InvertedCompressibilityFactorDetF;
    Real m_DegenerateThresholdDetF;
    Real m_ThresholdIpolDetF;
    float32 m_DAPD_L_Factor, m_DAPD_NL_Factor, m_DAPD_NL_Exponent;
    float32 m_ECIE_e_threshold, m_ECIE_k_factor;
    // Methods
    enum EMethod { // Linear
                   eMethod_L = 0,
                   // Corotational
                   eMethod_C_QR,
                   eMethod_C_SVD,
                   eMethod_C_PD_Project,
                   eMethod_C_PD_Reflect,
                   eMethod_C_DAPD,
                   // Hyperelastic-Corotational
                   eMethod_H_LCM_PD,
                   eMethod_H_CCM_PD,
                   // Hyperelastic
                   eMethod_H_LCM,
                   eMethod_H_CCM,
                   eMethod_H_NH_C0,
                   eMethod_H_NH_C1,
                   cNumMethods } m_Method;

    enum EMaterialMethod { // Linear
                           eMM_L = 0,
                           // Corotational
                           eMM_C_WRP,
                           eMM_C_LCM,
                           eMM_C_LCMH, //Hybrid LCM-WRP
                           eMM_C_CCM,
                           // Hyperelastic
                           eMM_H_LCM,
                           eMM_H_CCM,
                           eMM_H_NH_C0,
                           eMM_H_NH_C1,
                           cNumMM } m_MM;
    enum ERotationMethod { eRM_Id = 0,
                           eRM_QR,
                           eRM_MSVD,
                           eRM_DAPD,
                           cNumRM } m_RM;
    enum EDifferentialsMethod { eDM_Truncated,
                                eDM_Exact,
                                eDM_Exact_Inv,
                                eDM_Numerical,
                                cNumDM } m_DM;

    inline Params()
    // Material
    : m_YoungModulus(1000.0f), m_PoissonRatio(0.25f)
    , m_DampingRatio(1.0f)
      // Rayleigh
    , m_RayleighFreq1(0.1f), m_RayleighFreq2(10000.0f), m_RayleighCoeff(0,0)
      // Plastic
    , m_PlasticYield(0), m_PlasticMax(0), m_PlasticCreepPerSecond(0) //\note per-second as in Muller04 to achieve dt-independence.
      //, m_PlasticYield(0.01), m_PlasticMax(2.0), m_PlasticCreepPerSecond(1) //TEMP:
      // Misc
    , m_FactorDetF(0)
    , m_InvertedCompressibilityFactorDetF(0)
    , m_DegenerateThresholdDetF(0.1f)
    , m_ThresholdIpolDetF(0)
    , m_DAPD_L_Factor(1), m_DAPD_NL_Factor(0), m_DAPD_NL_Exponent(2) //exp >=2
    , m_ECIE_e_threshold(0.4f), m_ECIE_k_factor(1.0f)
    // Methods
    , m_Method( eMethod_H_LCM )//CCM ) //eMethod_C_SVD )//_PD )// eMethod_C_PD_DAPD )
    , m_MM( eMM_C_WRP )
    , m_RM( eRM_MSVD )
    , m_DM( eDM_Truncated )
        {
            m_RayleighCoeff[0] = 0;
            m_RayleighCoeff[1] = 0;
        }
};

}}} //namespace S2::ms::fem

#endif //S2_MS_FEM_PARAMS_H
