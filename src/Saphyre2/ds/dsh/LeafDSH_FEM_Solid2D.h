#ifndef S2_DS_DSH_LEAF_FEM_SOLID2D_H
#define S2_DS_DSH_LEAF_FEM_SOLID2D_H

#include <Saphyre2/ds/dsh/GLeafDSH_Model.h>
#include <Saphyre2/ms/ParticleSystem.h>

#include <Saphyre2/ms/fem/LinearTriangle2.h>
#include <vector>

#ifdef __S2_DS_ENABLE_PARAMS
#  include <util/Archetype.h>
#endif

#define __S2_DS_FEM_ENABLE_TEST

#define __ENABLE_CONTINUOUS
#define __ENABLE_SAVE_FORCES //TEMP: for Debug/Viz purposes

#define __ENABLE_CONTACT
#ifdef __ENABLE_CONTACT
#  define __ENABLE_CONTACT_D2K
#  ifdef __ENABLE_CONTACT_D2K
#    define __ENABLE_CONTACT_D2K_PRE //PRE Penalty works WAY better than POST!
#    include <Geo/mp/TestContact.h>
#    include <util/GPoolDA.h>
#  endif
#  define __ENABLE_CONTACT_CIRCLE //TEMPORAL: Fast and dirty test...
#endif

//Shape access... TEMPORAL
#include <Geo/ObjectFactory.h>
#include <Geo/shape/MeshSolidShape2.h>
#include <Geo/util/Viz.h>

#ifdef __S2_DS_ENABLE_STATS
#include <util/GStatician.h>
#endif

#define __USE_TRIANGLE_BASE //FASTER than using 3x3 barycentric coords 2x2 sub-matrix

//\todo THIS should be enabled by default
#define __ENABLE_RAYLEIGH_DAMPING
#ifdef __ENABLE_RAYLEIGH_DAMPING
//#  define __ENABLE_RIGID_MODE_NO_DAMPING //\todo Seems to be correct, but if enabled SE is too wiggly even at high damping_ratio...
#endif //__ENABLE_RAYLEIGH_DAMPING

#define __ENABLE_KNC
#define __ENABLE_KINEMATIC_B
#ifdef __ENABLE_KNC
#  include <util/GPoolDA.h>
#endif

//#  define __ENABLE_STATIC_ELEMENTS
#ifdef __ENABLE_STATIC_ELEMENTS
     const unsigned int cNumStaticElements = 4;
     const unsigned int cNumStaticNodes = 3*cNumStaticElements;
#else
     const unsigned int cNumStaticElements = 0;
     const unsigned int cNumStaticNodes = 0;
#endif

#define __ENABLE_VIZ_DEGENERATE_ELEMENTS

//#define __ENABLE_TRACE_IE
//#define __ENABLE_TRACE_IE_VERBOSE
//#define __ENABLE_TRACE_EIGENSTUFF

#include <Saphyre2/ns/MatrixD.h>
#define __ENABLE_TRACE_ASSEMBLE

#define __ENABLE_EIGENSTUFF

namespace S2 { namespace ds {

class LeafDSH_FEM_Solid2D: public GLeafDSH_Model< ms::ParticleSystem2D > //TEMPORAL should be ms::FEM_Solid2D
{
public:
    typedef geo::MeshSolidShape2 geo_shape_type;
    typedef geo::GObjectES<geo_shape_type> geo_object_type;
    enum EConstants {
        cNumEigenValues = 1
    };

public:
    LeafDSH_FEM_Solid2D( uint32 uid, IDynamicSystemHierarchy *p_parent );
    ~LeafDSH_FEM_Solid2D();
    EDSHType GetDSHType() const { return eDSH_Leaf_Solid2D; }
    void SyncGO();
    void RecomputeBV();
    bool Create( const ParamIt &pit );

    void Step( Real dt );
    void FixedStep( Real dt );

    bool Edit( const ParamIt &pit );
    bool Internal( const ParamIt &pit, ReturnStream &rets );

    void DoViz( util::VizStream &vs ) const;

    void QueryStats( Flags32 flags, ReturnStream &rets ) const;
    void QueryParams( Flags32 flags, ReturnStream &rets ) const;
    void SyncParams( const ParamIt &pit, ReturnStream &rets );

#ifdef __ENABLE_CONTACT_D2K
    void NotifyPairBP( const geo::bp::Proxy *p_other, const geo::bp::Pair *p_pair );
#endif

private:
    geo_object_type *m_pMeshO;
    geo_shape_type *m_pMeshS;
    unsigned int m_NumElements;
    ms::fem::LinearTriangle2 *m_vecElements;
    unsigned int m_NumNodes;
    Transform2 *m_vecTe2r;
    Mat2x2 *m_vecF;
    Real *m_vecDetF;
#ifdef __ENABLE_SAVE_FORCES
    Vec2 *m_vecForces;
#endif

#ifdef __S2_DS_ENABLE_PARAMS
    struct Params: public util::IArchetypeInstance
#else
    struct Params
#endif
    {
        LeafDSH_FEM_Solid2D *m_pDSH;

        // DEVELOP
        bool m_bRun;
        enum EMethod { eMethod_Id = 0,
                       eMethod_QR,
                       //-- PD-Methods
                       eMethod_PD,
                       eMethod_PD_Reflect,
                       eMethod_PD_Fix,
                       eMethod_PD_Project,
                       eMethod_PD_Project_Nearest,
                       //-- SVD-Methods
                       eMethod_PD_SVD,
                       eMethod_IHFSDM,
                       eMethod_ITF_LRM,
                       eMethod_ECIE_CLRM,
                       eMethod_ECIE_NHC0,
                       eMethod_ECIE_NHC1,
                       eMethod_PD_CLRM,
                       cNumMethods };
        EMethod m_Method;

        enum ENodeCollapseMethod { eNCM_NDN = 0,
                                   eNCM_ToC = 1,
                                   eNCM_IHFSDM = 2,
                                   cNumNCM };
        ENodeCollapseMethod m_NCM;

        Real m_FactorDetF;
        Real m_InvertedCompressibilityFactorDetF;
        Real m_DegenerateThresholdDetF;
        Real m_ThresholdIpolDetF;

        float32 m_ECIE_e_threshold, m_ECIE_k_factor;
        enum EKinematicMode { eKM_Disabled = 0, eKM_Static = 1, eKM_Periodic = 2, eKM_Mouse = 3, cNumKM };
        EKinematicMode m_KinematicMode;

#ifdef __S2_DS_FEM_ENABLE_TEST
        // Test
        enum ETestMode { eTM_Disabled = 0,
                         eTM_Random,
                         eTM_Plane,
                         eTM_Sphere,
                         cNumTM };
        ETestMode m_TestMode;
        enum ETestNodes { eTN_All = 0,
                          eTN_Boundary,
                          cNumTN };
        ETestNodes m_TestNodes;
        Real m_TestFraction;
#endif

        // Integrator
        enum EIntegratorType { eIT_SymplecticEuler = 0, eIT_ImplicitEuler = 1, cNumIT };
        EIntegratorType m_IntegratorType;
        bool m_bEnableVariableDT;
        Real m_FixedDT;
        Real m_TimeScale;
        // Material
        Real m_Thickness;
        Real m_Density;
        Real m_YoungModulus;
        Real m_PoissonRatio;
        Real m_DampingRatio;
        // Rayleigh
        Real m_RayleighFreq1;
        Real m_RayleighFreq2;
        Vec2 m_RayleighCoeff;
        // Environment
        Real m_AirDragCoeff;
        Vec2 m_Gravity;
        // LS Solver
        enum ELinearSystemSolverType { eLSST_CG = 0, eLSST_GS = 1, eLSST_Jacobi = 2, eLSST_LU = 3, cNumLSST };
        ELinearSystemSolverType m_SolverLS_Type;
        uint32 m_SolverLS_MaxIter;
#ifdef __S2_DS_ENABLE_PARAMS
        Real m_SolverLS_Log10_Epsilon;
#endif
        Real m_SolverLS_Epsilon; //derived
        // Contact Solver
        enum EContactSolverType { eCST_None = 0, eCST_Penalty = 1, eCST_Hack = 2, eCST_HackGroundPlane = 3, cNumCST };
        EContactSolverType m_ContactSolver_Type;
        Real m_ContactSolver_DepthMargin;
        Real m_ContactSolver_Penalty_Ks;
        Real m_ContactSolver_Penalty_Kd;
        Real m_ContactSolver_Relaxation_Coeff;
        Real m_ContactSolver_Restitution_Coeff;
        Real m_ContactSolver_MaxDepth;
        bool m_ContactSolver_BreakOnError;

        inline Params( LeafDSH_FEM_Solid2D *p_dsh = 0 )
        : m_pDSH(p_dsh)
          // Develop
        , m_bRun(true)

        , m_Method(eMethod_PD_Project)
        , m_NCM(eNCM_ToC)
        , m_FactorDetF(0)
        , m_InvertedCompressibilityFactorDetF(0)
        , m_DegenerateThresholdDetF(0.01f)
        , m_ThresholdIpolDetF(0)
        , m_ECIE_e_threshold(0.4f), m_ECIE_k_factor(1.0f)
        , m_KinematicMode(eKM_Static)
#ifdef __S2_DS_FEM_ENABLE_TEST
          // Test
        , m_TestMode(eTM_Disabled)
        , m_TestNodes(eTN_All)
        , m_TestFraction(0)
#endif
          // Integrator
        , m_IntegratorType(eIT_ImplicitEuler)
        , m_bEnableVariableDT(false)
        , m_FixedDT(0.01f)
        , m_TimeScale(1.0f)
          // Material
        , m_Thickness(1.0f)
        , m_Density(1000.0f) //kg/m^3
        , m_YoungModulus(1000.0f), m_PoissonRatio(0.25f)
        , m_DampingRatio(1.0f)
          // Rayleigh
        , m_RayleighFreq1(0.1f), m_RayleighFreq2(10000.0f), m_RayleighCoeff(0,0)
          // Environment
        , m_AirDragCoeff(0.0f)
          //, m_Gravity(0,-9.8f)
        , m_Gravity(0,0) //\todo should be a Creation and Editable param
          // Solver Internal
        , m_SolverLS_Type(eLSST_CG)
        , m_SolverLS_MaxIter(10)
#ifdef __S2_DS_ENABLE_PARAMS
        , m_SolverLS_Log10_Epsilon(-3.0f) //0.001
        , m_SolverLS_Epsilon( mal::Exp10( m_SolverLS_Log10_Epsilon ) )
#else
        , m_SolverLS_Epsilon( mal::Exp10( -3.0f ) )
#endif
          // Contact Solver
        , m_ContactSolver_Type(eCST_Penalty)
#define __USE_DEFAUT_SCALE_1M
#ifdef __USE_DEFAUT_SCALE_1M
        , m_ContactSolver_DepthMargin(0.001f)
        , m_ContactSolver_Penalty_Ks(10000)
        , m_ContactSolver_Penalty_Kd(100)
        , m_ContactSolver_Relaxation_Coeff(0.9)
        , m_ContactSolver_Restitution_Coeff(0)
        , m_ContactSolver_MaxDepth(1.0f)
#else
        , m_ContactSolver_DepthMargin(0.01f)
        , m_ContactSolver_Penalty_Ks(10000)
        , m_ContactSolver_Penalty_Kd(100)
        , m_ContactSolver_Relaxation_Coeff(0.9)
        , m_ContactSolver_Restitution_Coeff(0)
        , m_ContactSolver_MaxDepth(2.0f)
#endif
        , m_ContactSolver_BreakOnError(true)
            {
#ifdef __S2_DS_ENABLE_PARAMS
                RebuildRayleighCoeffsFromFreqs(this);
#else
                RebuildRayleighCoeffsFromFreqs();
#endif
            }
        inline bool Rebuild()
            {
                //\todo Anything missing???
                m_pDSH->NotifyChangedParams();
                return true;
            }

#ifdef __S2_DS_ENABLE_PARAMS
        void SetName( const char *name ) {}
        const char *GetName() const { return "default s2d::Params"; }
        static void InitArchetype( util::ArchetypeLibrary &al );
        static util::ArchetypeLibrary s_LeafDSH_FEM_Solid2D_Params_ArchetypeLibrary; //TEMP:
#endif

#ifdef __S2_DS_ENABLE_PARAMS
        static bool RebuildSolverLS_Epsilon( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            params.m_SolverLS_Epsilon = mal::Exp10( params.m_SolverLS_Log10_Epsilon );
            return true;
        }
        static bool RebuildAirDrag( IArchetypeInstance *p_instance )
        {
            static_cast<Params*>(p_instance)->m_pDSH->RebuildAirDrag();
            return true;
        }
        static bool RebuildRayleighCoeffsFromFreqs( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            Mat2x2 Omega;
            Omega(0,0) = mal::Rcp(params.m_RayleighFreq1); Omega(0,1) = params.m_RayleighFreq1;
            Omega(1,0) = mal::Rcp(params.m_RayleighFreq2); Omega(1,1) = params.m_RayleighFreq2;
            Vec2 zeta( params.m_DampingRatio, params.m_DampingRatio );
            params.m_RayleighCoeff = mal::Inverse(Omega) * (2*zeta);
            return true;
        }
#else
        void RebuildRayleighCoeffsFromFreqs()
        {
            Mat2x2 Omega;
            Omega(0,0) = mal::Rcp(m_RayleighFreq1); Omega(0,1) = m_RayleighFreq1;
            Omega(1,0) = mal::Rcp(m_RayleighFreq2); Omega(1,1) = m_RayleighFreq2;
            Vec2 zeta( m_DampingRatio, m_DampingRatio );
            m_RayleighCoeff = mal::Inverse(Omega) * (2*zeta);
        }
#endif
    };
    Params m_Params;

#ifdef __ENABLE_KNC
    struct KinematicNodeConstraint
    {
        enum { eKNCT_Position, eKNCT_Velocity } m_Type;
        uint16 m_nid;
        Vec2 m_Pos;
        Vec2 m_Vel;
        finline KinematicNodeConstraint() : m_Type(eKNCT_Position), m_nid(), m_Pos(0,0), m_Vel(0,0) {}
        finline KinematicNodeConstraint( uint16 nid, const Vec2 &pos ) : m_Type(eKNCT_Position), m_nid(nid), m_Pos(pos), m_Vel(0,0) {}
    };
    typedef util::GPoolDA<KinematicNodeConstraint> PoolKNC;
    PoolKNC m_poolKNC;
    finline bool IsKNC( uint16 nid ) const
        {
            for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
                if( nid == it_knc->m_nid ) return true;
            return false;
        }
#endif

    //\name Derived params
    //@{
    Real m_TotalMass;
    Real m_TotalArea;
    Real m_DensityPerArea;
    Real m_AirDragPerSecond;
    //@}

    //\name State
    //@{
    Real m_TotalTime;
    Real m_AccTime;
    Real m_ElasticEnergy;
    Real m_PotentialEnergy;
    Real m_KineticEnergy;
    //@}

    RadialPressureAtPoint2D m_LastRPAP;
    geo::bv::AABB2 m_AABB0;

private:

    void RecomputeEnergy();
    inline Real GetElasticEnergy() const { return m_ElasticEnergy; }
    inline Real GetPotentialEnergy() const { return m_PotentialEnergy; }
    inline Real GetKineticEnergy() const { return m_KineticEnergy; }

#ifdef __S2_DS_ENABLE_STATS
    struct Stats
    {
        uint32 m_Num_Degenerate;
        float32 m_Sum_Degenerate_DetF;
        uint32 m_SolverLS_NumIter;
        float32 m_SolverLS_Prec;
        util::GStatician<float32> m_RelArea;
        void Begin() { /* m_Num_Degenerate = 0; m_Sum_Degenerate_DetF = 0; */ m_SolverLS_NumIter = 0; m_SolverLS_Prec = 0; m_RelArea.Clear(); }
        void End() { m_RelArea.End(); }
    };
    mutable Stats m_Stats;
#endif

    // Eigenstuff
    ns::RealD m_vecEigenValue[cNumEigenValues];
    ns::VectorD m_vecEigenVector[cNumEigenValues];

#ifdef __ENABLE_CONTINUOUS
    struct DegeneratedElementData
    {
        uint32 m_ElementId;
        uint16 m_NoC; //Node-of-Collapse, local node index inside element [0..2]
        uint16 m_Age; //How long it's been degenerated
        Real m_ToC; //Time-of-Collapse [0..1]
        finline void Init( uint32 element_id, uint16 noc, Real toc ) { m_ElementId = element_id; m_NoC = noc; m_Age = 0; m_ToC = toc; }
    };
    typedef util::GPoolDA<DegeneratedElementData> PoolDED;
    PoolDED m_poolDED;
    Vec2 *m_vecPos0;
    void UpdateDED();
#endif

    // Assembled matrices and vectors
    ns::MatrixD m_M;
    ns::MatrixD m_K;
    //ns::MatrixD m_C;
    ns::MatrixD m_A;
    ns::VectorD m_Fe;

#ifdef __ENABLE_CONTACT_D2K
    struct ContactPointConstraintD2K
    {
        /*! Specific D2K contact data, may include baricentric coords,
         whatever... We ONLY store info for the Solid, the Geom is
         Kinematic and therefore we don't need its specific contact
         info, only the part required to compute response on the
         Solid. */
        ContactPointConstraintD2K( const geo::PointOnFeature &pof1, const Vec2 &normal, Real depth, Real radius )
        : m_POF1(pof1), m_Normal(normal), m_Depth(depth), m_Radius(radius) {}
        geo::PointOnFeature m_POF1;
        Vec2 m_Normal;
        Real m_Depth;
        Real m_Radius;
    };

    struct ContactConstraintD2K
    {
        /*\todo We'd like to have a FIXED cMaxCPC, but it's not
          realistic unless we come up with a good contact reduction
          scheme, and it's MUCH MORE complicated than the rigid case,
          because CP inside the "convex hull" are ALSO relevant
          here... so by now we store them all!
        enum EConstants { cMaxCPC = 8 };
        ContactPointConstraintD2K m_vecCPC[cMaxCPC];

        \todo Actually, if individual CPC persistence and vanishing is
        required, a better allocator with cheaper new/delete and
        iteration could be algorithmically better (eg: a dynamic
        iterable pool), BUT, depending on how we use the CPC in the
        solver, if >1 passes are required, an array will surely be
        faster due to contiguous memory access.
        */

        inline void Create() { Destroy(); m_TmpCD.Clear(); }
        inline void Destroy() { m_Age = 0; m_CC.Invalidate(); Clear(); }

        inline void Clear() { m_vecCPC.clear(); }

        inline void Reset( const geo::np::ContactData2 &cd, Real epsilon_length )
            {
                //\todo Use epsilon_length to merge close points, ONLY IF NOT already done by default with cd.Optimize(epsilon_length)
                DS_ASSERT( cd.HasPOF() );
                m_Age = 0;
                m_vecCPC.clear();
                for( unsigned int it_cp=0; it_cp<cd.Size(); it_cp++ )
                    m_vecCPC.push_back( ContactPointConstraintD2K( cd.GetPOF1(it_cp),
                                                                   cd.GetCP(it_cp).m_Normal,
                                                                   cd.GetCP(it_cp).m_Depth,
                                                                   cd.GetCP(it_cp).m_Radius ) );
            }
        inline void Persist( const geo::np::ContactData2 &cd, Real epsilon_length, Real epsilon_direction )
            {
                // \todo Match and persist...
                m_Age++;
                Reset(cd,epsilon_length);
            }

        finline bool IsEmpty() const { return m_vecCPC.empty(); }
        finline unsigned int GetNumCPC() const { return m_vecCPC.size(); }
        finline const ContactPointConstraintD2K& GetCPC( int i ) const { return m_vecCPC[i]; }

        uint32 m_Age;
        geo::np::ContactCache2 m_CC;
        std::vector<ContactPointConstraintD2K> m_vecCPC;

        geo::np::ContactData2 m_TmpCD; //TEMPORAL: Single CD stored to Viz later, otherwise local var
    };
    typedef util::GPoolDA<ContactConstraintD2K> PoolCCD2K;
    PoolCCD2K m_poolCCD2K;
    bool CheckError_Contact() const;
#endif //__ENABLE_CONTACT_D2K

private:
    //\name Dynamics
    //@{
    void AccumulateExplicitForces( Real dt );
    void Solve_SymplecticEuler( Real dt );

    void Solve_ImplicitEuler( Real dt );

    // Corotational
    int FindNDN( unsigned int index_e ) const;
    int FindNoC( unsigned int index_e ) const;
    void Compute_F( Mat2x2 *vec_F, Real *vec_det_F ) const;
    void Compute_T( Transform2 *vec_te2r ) const;

    // Assembly
    Real Compute_Effective_PoissonRatio( Real det_F ) const;
    Mat6x6 Compute_Effective_Ke( unsigned int index_e ) const;
    void Assemble_K( ns::MatrixD &K ) const;
    void Assemble_M( ns::MatrixD &M ) const;
    void Assemble_InvM( ns::MatrixD &InvM ) const;
    void Assemble_u( ns::VectorD &u ) const;
    void Assemble_v( ns::VectorD &v ) const;
    void Assemble_x( ns::VectorD &x ) const;
    void Assemble_r( ns::VectorD &r ) const;
    void Assemble_CorotationalElasticForce( ns::VectorD &cef ) const;
    void Assemble_Fe( ns::VectorD &Fe ) const;
    void Assemble_A( ns::MatrixD &A, Real dt, const ns::MatrixD &M, const ns::MatrixD &K );

    //Piola-Kirchhoff stuff
    void Compute_Effective_LameParams( Real det_F, Real &lame_mu, Real &lame_lambda ) const;
    Vec2 ComputeDiagonalP_ITF_LRM( const Vec2 &vec_diag_F ) const;
    Vec2 ComputeDiagonalP_ECIE_CLRM( const Vec2 &vec_diag_F ) const;
    Vec2 ComputeDiagonalP_ECIE_NHC0( const Vec2 &vec_diag_F ) const;
    Vec2 ComputeDiagonalP_ECIE_NHC1( const Vec2 &vec_diag_F ) const;
    //@}

private:
    void RebuildAirDrag();
    void RebuildMass();
    void NotifyChangedParams();
    void RebuildEigenstuff();

private:
    unsigned int FindClosestNode( const Vec2 &pos ) const;
    Vec2 ComputePosCoM() const;

private:
    //\name Edition
    //@{
    void SetPosCoM( const Vec2 &pos );
    //@}

private:
    // Tests
#ifdef __S2_DS_FEM_ENABLE_TEST
    void InitTest();
    void UpdateTest( Real dt );
    //Array of test-owned KNC
    std::vector<KinematicNodeConstraint*> m_vecTestKNC;
    std::vector<Vec2> m_vec_TM_Random_TargetPos;
#endif

private:
    //! Draw Debug Flags
    enum EDDF { eDraw_Nothing        = 0,
                eDraw_Mesh           = 1<<0,
                eDraw_KNC            = 1<<1,
                eDraw_RPAP           = 1<<2,
                eDraw_Degenerate     = 1<<3,
                eDraw_Corotational   = 1<<4,
                eDraw_CoM            = 1<<5,
                eDraw_Contacts       = 1<<6,
                eDraw_Eigenbasis     = 1<<7,
                eDraw_Forces         = 1<<8,
                eDraw_Trajectory     = 1<<9,
                eDraw_DED            = 1<<10,
                eDraw_Everything     = 0xFFFFFFFF,
                eDraw_Default        = eDraw_Mesh | eDraw_KNC | eDraw_RPAP | eDraw_Contacts /*| eDraw_Degenerate | eDraw_Corotational | eDraw_CoM | eDraw_Eigenbasis */  | eDraw_Forces | eDraw_Trajectory | eDraw_DED
    };
    Flags32 m_DDF;
    Flags32 m_GeoObjectDDF;
    Flags32 m_ContactDDF;
};

} } // namespace S2::ds

#endif // S2_DS_DSH_LEAF_FEM_SOLID2D_H
