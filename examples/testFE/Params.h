#ifndef TEST_FE_PARAMS_H
#define TEST_FE_PARAMS_H

#include "Config.h"
#include <util/Archetype.h>

const double cDefaultTimeStep(0.01f);

class Params: public util::IArchetypeInstance
{
public:
    enum EFixFMethod { eFFM_None = 0,
                       eFFM_Project,
                       eFFM_Project_Nearest,
                       eFFM_Reflect,
                       cNumFFM };

    enum EMethod { eMethod_Id = 0,
                   eMethod_QR_XY,
                   eMethod_QR_YX,
                   eMethod_PD,
                   eMethod_PD_Project,
                   eMethod_PD_Reflect,
                   eMethod_PD_Unified,
                   eMethod_PD_Fix,
                   eMethod_PD_Project_Nearest,
                   eMethod_PD_SVD,
                   eMethod_IHFSDM,
                   eMethod_ITF_LRM,
                   eMethod_ECIE_CLRM,
                   eMethod_ECIE_NHC0,
                   eMethod_ECIE_NHC1,
                   eMethod_PD_CLRM,
                   eMethod_DAPD_NS, //NonSymmetricStrain
                   eMethod_DAPD_SS, //SymmetricStrain
                   eMethod_DAPD_EX, //Exact(numeric)
                   cNumMethods };
    enum ETrajectoryFlags { eTF_Nothing          = 0,
                            eTF_Free_x0          = (1<<0),
                            eTF_Free_x1          = (1<<1),
                            eTF_Free_x2          = (1<<2),
                            eTF_Free_x3          = (1<<3),
                            eTF_TotalEnergy      = (1<<4),
                            eTF_ElasticEnergy    = (1<<5),
                            eTF_KineticEnergy    = (1<<6),
                            eTF_AngularMomentum  = (1<<7),
                            eTF_DeformationRatio = (1<<8),
                            eTF_Default          = eTF_Free_x0 };
    enum EDrawFlags { eDraw_Nothing  = 0,
                      eDraw_Axis     = (1<<0),
                      eDraw_Grid     = (1<<1),
                      eDraw_Element  = (1<<2),
                      eDraw_f0_x     = (1<<3),
                      eDraw_f1_x     = (1<<4),
                      eDraw_f2_x     = (1<<5),
                      eDraw_f3_x     = (1<<6),
                      eDraw_P_e      = (1<<7),
                      eDraw_R_x      = (1<<8),
                      eDraw_Det_F    = (1<<9),
                      eDraw_Drag_f0  = (1<<10),
                      eDraw_Drag_EE0 = (1<<11),
                      eDraw_Drag_t0  = (1<<12),
                      eDraw_Drag_R   = (1<<13),
                      eDraw_Drag_dR  = (1<<14),
                      eDraw_Drag_df  = (1<<15),
                      eDraw_Default  = eDraw_Element | eDraw_f0_x /*| eDraw_Drag_t0*/ };

public:
    Params()
    : m_YoungModulus(1000.0f), m_PoissonRatio(0.35f)
    , m_FactorDetF(0), m_InvertedCompressibilityFactorDetF(0)
    , m_DegenerateThresholdDetF(0.1f), m_ThresholdIpolDetF(0)
    , m_DiffR_Numerical_H(0.01), m_DiffR_dX12_Bias_Factor(0.5)
    , m_Alpha(0), m_Beta(0), m_Angle0(0), m_Dist12(2), m_Scale12(1)
    , m_FFM(eFFM_None)
    , m_Method(eMethod_DAPD_EX)
    , m_Unified_L_Factor(2), m_Unified_NL_Factor(0), m_Unified_NL_Exponent(2)
    , m_ECIE_e_threshold(0.4f), m_ECIE_k_factor(10.0f) //\note For small values the force field is quite bad and inversion-increasing, ECIE paper used k=20 AFAICR
    // , m_HalfSizeX(8), m_HalfSizeY(8), m_HalfSizeZ(8)
    , m_HalfSizeX(3), m_HalfSizeY(3), m_HalfSizeZ(3)
    // , m_HalfCellsX(49), m_HalfCellsY(49), m_HalfCellsZ(49)
    , m_HalfCellsX(20), m_HalfCellsY(20), m_HalfCellsZ(20)
    , m_bSingleSliceX0(false), m_bSingleSliceY0(false), m_bSingleSliceZ0(true)
    , m_SliceX0(-20), m_SliceX1(20), m_SliceY0(-20), m_SliceY1(20), m_SliceZ0(0), m_SliceZ1(0)
    , m_NID(0) //varying node id
    , m_NoC(0) //gui-defined NoC
    , m_DragPointZ(0)
    , m_bEnableNoC(true), m_TimeStep(0.001), m_Duration(1), m_Mass(1), m_DampingRatio(0.1)
    , m_TrajectoryFlags(Params::eTF_Default)
    , m_DrawFlags(Params::eDraw_Default)
        {}
    ~Params() {}

    //!\name IArchetypeInstance implementation
    //@{
    bool Rebuild()
        {
            //APP_LOG("Rebuilding...");
            if( m_bSingleSliceX0 ) m_SliceX1 = m_SliceX0;
            else if( m_SliceX0 > m_SliceX1 ) m_SliceX0 = m_SliceX1; //\note We don't know which one changed, just force them to be equal
            if( m_bSingleSliceY0 ) m_SliceY1 = m_SliceY0;
            else if( m_SliceY0 > m_SliceY1 ) m_SliceY0 = m_SliceY1; //\note We don't know which one changed, just force them to be equal
            if( m_bSingleSliceZ0 ) m_SliceZ1 = m_SliceZ0;
            else if( m_SliceZ0 > m_SliceZ1 ) m_SliceZ0 = m_SliceZ1; //\note We don't know which one changed, just force them to be equal
            return true;
        }
    void SetName( const char *name ) {}
    const char *GetName() const { return "Unnamed_Params"; }
    static void InitArchetype( util::ArchetypeLibrary &al )
        {
            al.BeginArchetype( "Archetype_SP" );
            {
                Params sp;
                al.BeginProperty_Group( "<Material>" );
                {
                    al.AddProperty_NIR<float32>( "young modulus", 1000.0f, 100.0f, 10000.0f, archetype_offset_of(sp,m_YoungModulus), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "poisson ratio", 0.25f, 0.0f, 0.49f, archetype_offset_of(sp,m_PoissonRatio), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<Geometry>" );
                {
                    al.AddProperty_NIR<float32>( "alpha", 0, -0.99f*mal::Pi<float32>(), 0.99f*mal::Pi<float32>(), archetype_offset_of(sp,m_Alpha), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "beta", 0, -0.99f*mal::Pi<float32>(), 0.99f*mal::Pi<float32>(), archetype_offset_of(sp,m_Beta), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "angle0", 0, -0.99f*mal::HalfPi<float32>(), 0.99f*mal::HalfPi<float32>(), archetype_offset_of(sp,m_Angle0), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "dist12", 1.0f, 0.1f, 10.0f, archetype_offset_of(sp,m_Dist12), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "scale12", 1.0f, 0.1f, 5.0f, archetype_offset_of(sp,m_Scale12), IArchetypeInstance::NTPF_Rebuild );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<Thresholds>" );
                {
                    al.AddProperty_NIR<float32>( "factor det(F)", 0.0f, 0.0f, 1.0f, archetype_offset_of(sp,m_FactorDetF), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "ICF det(F)", 0.0f, 0.0f, 5.0f, archetype_offset_of(sp,m_InvertedCompressibilityFactorDetF), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "eps.det(F)", 0.01f, 0.001f, 1.0f, archetype_offset_of(sp,m_DegenerateThresholdDetF), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<dR>" );
                {
                    al.AddProperty_NIR<float32>( "dR h", 0.01f, 0.001f, 0.1f, archetype_offset_of(sp,m_DiffR_Numerical_H), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "dR bf", 0.5f, 0.0f, 2000.0f, archetype_offset_of(sp,m_DiffR_dX12_Bias_Factor), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<Methods>" );
                {
                    const char *vec_names_ffm[] = { "None", "P", "PN", "R" };
                    uint32 vec_values_ffm[] = { eFFM_None,
                                                eFFM_Project,
                                                eFFM_Project_Nearest,
                                                eFFM_Reflect };
                    al.AddProperty_Enum32( "Fix(F)",
                                           (uint32)eFFM_None,
                                           (uint32)cNumFFM, vec_names_ffm, vec_values_ffm,
                                           archetype_offset_of(sp,m_FFM),
                                           IArchetypeInstance::NTPF_Ignore );

                    const char *vec_names_mt[] = { "Id", "QR_XY", "QR_YX",
                                                   "PD", "PD_P", "PD_R", "PD_U", "PD_F", "PD_PN", "PD_SVD",
                                                   "IHFSDM", "ITF_LRM", "ECIE_CLRM", "ECIE_NHC0", "ECIE_NHC1",
                                                   "PD_CLRM",
                                                   "DAPDNS", "DAPDSS", "DAPDEX" };
                    uint32 vec_values_mt[] = { eMethod_Id,
                                               eMethod_QR_XY,
                                               eMethod_QR_YX,
                                               eMethod_PD,
                                               eMethod_PD_Project,
                                               eMethod_PD_Reflect,
                                               eMethod_PD_Unified,
                                               eMethod_PD_Fix,
                                               eMethod_PD_Project_Nearest,
                                               eMethod_PD_SVD,
                                               eMethod_IHFSDM,
                                               eMethod_ITF_LRM,
                                               eMethod_ECIE_CLRM,
                                               eMethod_ECIE_NHC0,
                                               eMethod_ECIE_NHC1,
                                               eMethod_PD_CLRM,
                                               eMethod_DAPD_NS,
                                               eMethod_DAPD_SS,
                                               eMethod_DAPD_EX };
                    al.AddProperty_Enum32( "Method",
                                           (uint32)eMethod_Id,
                                           (uint32)cNumMethods, vec_names_mt, vec_values_mt,
                                           archetype_offset_of(sp,m_Method),
                                           IArchetypeInstance::NTPF_Ignore );

                    al.AddProperty_NIR<float32>( "ipol.det(F)", 0.0f, -5.0f, 0.0f, archetype_offset_of(sp,m_ThresholdIpolDetF), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "PD_U L factor", 1.0f, 1.0f, 10.0f, archetype_offset_of(sp,m_Unified_L_Factor), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "PD_U NL factor", 0.0f, 0.0f, 10.0f, archetype_offset_of(sp,m_Unified_NL_Factor), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "PD_U NL exponent", 1.0f, 1.0f, 10.0f, archetype_offset_of(sp,m_Unified_NL_Exponent), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "ECIE e threshold", 0.4, 0.01, 0.99, archetype_offset_of(sp,m_ECIE_e_threshold), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "ECIE k factor", 10.0f, 0.0f, 20.0f, archetype_offset_of(sp,m_ECIE_k_factor), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<Domain>" );
                {
                    al.AddProperty_NIR<float32>( "hs X", 5.0f, 1.0f, 50.0f, archetype_offset_of(sp,m_HalfSizeX), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "hs Y", 5.0f, 1.0f, 50.0f, archetype_offset_of(sp,m_HalfSizeY), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "hs Z", 5.0f, 1.0f, 50.0f, archetype_offset_of(sp,m_HalfSizeZ), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "hc X", 20, 10, 100, archetype_offset_of(sp,m_HalfCellsX), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "hc Y", 20, 10, 100, archetype_offset_of(sp,m_HalfCellsY), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<uint32>( "hc Z", 20, 10, 100, archetype_offset_of(sp,m_HalfCellsZ), IArchetypeInstance::NTPF_Ignore );

                    al.AddProperty<bool>( "single X0", false, archetype_offset_of(sp,m_bSingleSliceX0), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<int32>( "slice X0",  0,  -20, 20, archetype_offset_of(sp,m_SliceX0), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<int32>( "slice X1",  0,  -20, 20, archetype_offset_of(sp,m_SliceX1), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty<bool>( "single Y0", false, archetype_offset_of(sp,m_bSingleSliceY0), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<int32>( "slice Y0",  0,  -20, 20, archetype_offset_of(sp,m_SliceY0), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<int32>( "slice Y1",  0,  -20, 20, archetype_offset_of(sp,m_SliceY1), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty<bool>( "single Z0", true, archetype_offset_of(sp,m_bSingleSliceZ0), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<int32>( "slice Z0",  0,  -20, 20, archetype_offset_of(sp,m_SliceZ0), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<int32>( "slice Z1",  0,  -20, 20, archetype_offset_of(sp,m_SliceZ1), IArchetypeInstance::NTPF_Rebuild );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<Misc>" );
                {
                    al.AddProperty_NIR<uint32>( "NID", 0, 0, 2, archetype_offset_of(sp,m_NID), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<int32>( "NoC", 0, 0, 2, archetype_offset_of(sp,m_NoC), IArchetypeInstance::NTPF_Rebuild );
                    al.AddProperty_NIR<float32>( "drag Z", 0.0f, -1.0f, 1.0f, archetype_offset_of(sp,m_DragPointZ), IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<Trajectory>" );
                {
                    al.AddProperty<bool>( "NoC?", true, archetype_offset_of(sp,m_bEnableNoC), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "dt", 0.001f, 0.0001f, 0.01f, archetype_offset_of(sp,m_TimeStep), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "Duration", 1.0f, 0.1f, 10.0f, archetype_offset_of(sp,m_Duration), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "Mass", 1.0f, 0.1f, 10.0f, archetype_offset_of(sp,m_Mass), IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_NIR<float32>( "DampingRatio", 0.1f, 0.0f, 1.0f, archetype_offset_of(sp,m_DampingRatio), IArchetypeInstance::NTPF_Ignore );
                    const char *vec_names_tf[] = { "x0", "x1", "x2", "x3", "TE", "EE", "KE", "L", "det(F)" };
                    int32 vec_values_tf[] = { (int32)eTF_Free_x0,
                                              (int32)eTF_Free_x1,
                                              (int32)eTF_Free_x2,
                                              (int32)eTF_Free_x3,
                                              (int32)eTF_TotalEnergy,
                                              (int32)eTF_ElasticEnergy,
                                              (int32)eTF_KineticEnergy,
                                              (int32)eTF_AngularMomentum,
                                              (int32)eTF_DeformationRatio };
                    al.AddProperty_Flags32( "TF",
                                            Params::eTF_Default,
                                            9, vec_names_tf, vec_values_tf,
                                            archetype_offset_of(sp,m_TrajectoryFlags),
                                            IArchetypeInstance::NTPF_Ignore );

                }
                al.EndProperty_Group();

                // DDF
                const char *vec_names_df[] = { "Axis", "Grid", "Elem",
                                               "f0", "f1", "f2", "f3",
                                               "P(e)",
                                               "R(x)", "Det(F)",
                                               "f0(x0)",  "EE(x0)", "X(t)", "R(t)", "dR(x)", "df(x)" };
                int32 vec_values_df[] = { (int32)eDraw_Axis,
                                          (int32)eDraw_Grid,
                                          (int32)eDraw_Element,
                                          (int32)eDraw_f0_x,
                                          (int32)eDraw_f1_x,
                                          (int32)eDraw_f2_x,
                                          (int32)eDraw_f3_x,
                                          (int32)eDraw_P_e,
                                          (int32)eDraw_R_x,
                                          (int32)eDraw_Det_F,
                                          (int32)eDraw_Drag_f0,
                                          (int32)eDraw_Drag_EE0,
                                          (int32)eDraw_Drag_t0,
                                          (int32)eDraw_Drag_R,
                                          (int32)eDraw_Drag_dR,
                                          (int32)eDraw_Drag_df };
                al.AddProperty_Flags32( "Draw",
                                        Params::eDraw_Default,
                                        16, vec_names_df, vec_values_df,
                                        archetype_offset_of(sp,m_DrawFlags),
                                        IArchetypeInstance::NTPF_Ignore );
            }
            al.EndArchetype();
        }
    //@}

    std::string GetMethodName( EMethod method ) const
        {
            switch( method )
            {
            case eMethod_Id: return std::string("Id"); break;
            case eMethod_QR_XY: return std::string("QR_XY"); break;
            case eMethod_QR_YX: return std::string("QR_YX"); break;
            case eMethod_PD: return std::string("PD"); break;
            case eMethod_PD_Project: return std::string("PD_P"); break;
            case eMethod_PD_Reflect: return std::string("PD_R"); break;
            case eMethod_PD_Unified: return std::string("PD_U"); break;
            case eMethod_PD_Fix: return std::string("PD_F"); break;
            case eMethod_PD_Project_Nearest: return std::string("PD_PN"); break;
            case eMethod_PD_SVD: return std::string("PD_SVD"); break;
            case eMethod_IHFSDM: return std::string("IHFSDM"); break;
            case eMethod_ITF_LRM: return std::string("ITF_LRM"); break;
            case eMethod_ECIE_CLRM: return std::string("ECIE_CLRM"); break;
            case eMethod_ECIE_NHC0: return std::string("ECIE_NH0"); break;
            case eMethod_ECIE_NHC1: return std::string("ECIE_NH1"); break;
            case eMethod_PD_CLRM: return std::string("PD_CLRM"); break;
            case eMethod_DAPD_NS: return std::string("DAPD_NS"); break;
            case eMethod_DAPD_SS: return std::string("DAPD_SS"); break;
            case eMethod_DAPD_EX: return std::string("DAPD_EX"); break;
            default: return std::string("NONE"); break;
            }
        }

public:

    float32 m_YoungModulus;
    float32 m_PoissonRatio;

    float32 m_FactorDetF;
    float32 m_InvertedCompressibilityFactorDetF;
    float32 m_DegenerateThresholdDetF;
    float32 m_ThresholdIpolDetF;
    float32 m_DiffR_Numerical_H;
    float32 m_DiffR_dX12_Bias_Factor;

    float32 m_Alpha;
    float32 m_Beta;
    float32 m_Angle0;
    float32 m_Dist12;
    float32 m_Scale12;

    uint32 m_FFM;
    uint32 m_Method;
    float32 m_Unified_L_Factor, m_Unified_NL_Factor, m_Unified_NL_Exponent;
    float32 m_ECIE_e_threshold, m_ECIE_k_factor;

    float32 m_HalfSizeX, m_HalfSizeY, m_HalfSizeZ;
    uint32 m_HalfCellsX, m_HalfCellsY, m_HalfCellsZ;

    bool m_bSingleSliceX0,m_bSingleSliceY0,m_bSingleSliceZ0;
    int32 m_SliceX0, m_SliceX1;
    int32 m_SliceY0, m_SliceY1;
    int32 m_SliceZ0, m_SliceZ1;

    uint32 m_NID;
    int32 m_NoC;

    float32 m_DragPointZ;

    bool m_bEnableNoC;
    float32 m_TimeStep, m_Duration, m_Mass, m_DampingRatio;
    Flags32 m_TrajectoryFlags;
    Flags32 m_DrawFlags;
};

inline void LameParameters_From_YoungAndPoisson( Real young_modulus, Real poisson_ratio,
                                                 Real &lame_mu, Real &lame_lambda )
{
    //\todo THIS MAY BE DIFFERENT IN 2D ?!?!?!?!
    lame_mu = Real(0.5) * young_modulus / (1+poisson_ratio);
    lame_lambda = young_modulus * poisson_ratio / ( (1+poisson_ratio) * (1-2*poisson_ratio) );
}

#endif //TEST_FE_PARAMS_H
