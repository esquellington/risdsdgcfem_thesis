#ifndef S2_DS_DSH_LEAF_SOLID3D_FEM_H
#define S2_DS_DSH_LEAF_SOLID3D_FEM_H

#include <Saphyre2/ds/dsh/GLeafDSH_Model.h>

#include <Saphyre2/ms/ParticleSystem.h>
#include <Saphyre2/ms/fem/TetrahedronElementSet3.h>

#include <vector>

#ifdef __S2_DS_ENABLE_PARAMS
#  include <util/Archetype.h>
#endif

#define __S3D_ENABLE_SAVE_FORCES //TEMP: for Debug/Viz purposes

//Shape access... TEMPORAL
#include <Geo/ObjectFactory.h>
#include <Geo/shape/TetSolidShape3.h>
#include <Geo/util/Viz.h>

#ifdef __S2_DS_ENABLE_STATS
#  include <util/GStatician.h>
#  include <util/Chronos.h>
#endif

//#define __S3D_ENABLE_VIZ_DEGENERATE_ELEMENTS

#define __S3D_ENABLE_KNC
#ifdef __S3D_ENABLE_KNC
#  define __S3D_S2_DS_FEM_ENABLE_TEST
#  define __S3D_ENABLE_KNC_IN_LS
#  include <Saphyre2/ms/constraint/KinematicNodeConstraintSet3.h>
#endif

#define __S3D_ENABLE_CONTACT
#ifdef __S3D_ENABLE_CONTACT
#  define __S3D_ENABLE_CONTACT_D2K
#  ifdef __S3D_ENABLE_CONTACT_D2K
#    include <Geo/mp/TestContact.h>
#    include <util/GPoolDA.h>
#  endif
#endif

#define __S3D_ENABLE_BVH
#ifdef __S3D_ENABLE_BVH
#  include <Geo/bv/GBoundingVolumeHierarchy.h>
#endif

// ILSS
#include <Saphyre2/ns/ILSS_CG.h>
#include <Saphyre2/ns/ILSS_CR.h>
#include <Saphyre2/ns/ILSS_SYMMQMR.h>
#include <Saphyre2/ns/ILSS_MINRES.h>
#include <Saphyre2/ns/ILSS_GMRES.h>

namespace S2 { namespace ds {

class LeafDSH_Solid3D_FEM: public GLeafDSH_Model< ms::ParticleSystem3D >
{
public:
    LeafDSH_Solid3D_FEM( uint32 uid, IDynamicSystemHierarchy *p_parent );
    ~LeafDSH_Solid3D_FEM();
    EDSHType GetDSHType() const { return eDSH_Leaf_Solid3D; }

    //!\name IEntity implementation
    //@{
    bool Create( const ParamIt &pit );
    bool Edit( const ParamIt &pit );
    bool Internal( const ParamIt &pit, ReturnStream &rets );
    void DoViz( util::VizStream &vs ) const;
    void QueryStats( Flags32 flags, ReturnStream &rets ) const;
    void QueryParams( Flags32 flags, ReturnStream &rets ) const;
    void SyncParams( const ParamIt &pit, ReturnStream &rets );
    //@}

    //!\name IDSH Local update sub-phases
    //@{
    void Step( Real dt );
    void SyncGO();
    void RecomputeBV();
    //@}

#ifdef __S3D_ENABLE_CONTACT_D2K
    //!\name Broad-Phase related interface, only used by SS
    //@{
    void NotifyPairBP( const geo::bp::Proxy *p_other, const geo::bp::Pair *p_pair );
    //@}
#endif

private:
#ifdef __S2_DS_ENABLE_PARAMS
    struct Params: public util::IArchetypeInstance
#else
    struct Params
#endif
    {
        LeafDSH_Solid3D_FEM *m_pDSH;

        // FEM params
        ms::fem::Params m_FEM;

        // Develop
        bool m_bRun;
        enum EKinematicMode { eKM_Disabled = 0, eKM_Static = 1, eKM_Periodic = 2, eKM_Mouse = 3, cNumKM };
        EKinematicMode m_KinematicMode;

        // Integrator
        enum EIntegratorType { eIT_SymplecticEuler = 0, eIT_QuasiImplicitEuler_v = 1, eIT_QuasiImplicitEuler_dx = 2, eIT_FullyImplicitEuler_dx = 3, eIT_FullyImplicitEuler_dx_QuasiStatic = 4, cNumIT };
        EIntegratorType m_IntegratorType;
        bool m_bEnableVariableDT;
        Real m_FixedDT;
        Real m_TimeScale;

        // Material
        Real m_Density;

        // Environment
        Real m_AirDragCoeff;
        Vec3 m_Gravity;

        // LS Solver
        enum ELinearSystemSolverType { eLSST_CG = 0,
                                       eLSST_CR,
                                       eLSST_SQMR,
                                       eLSST_MINRES,
                                       eLSST_GMRES,
                                       cNumLSST, //\todo avoid usless LSST by now
                                       eLSST_CGN,
                                       eLSST_CGS,
                                       eLSST_RI,
                                       eLSST_GS,
                                       eLSST_Jacobi,
                                       eLSST_LU };
        ELinearSystemSolverType m_SolverLS_Type;
        uint32 m_SolverLS_MaxIter;
        uint32 m_SolverLS_MaxRestart; //GMRES
        bool m_SolverLS_bEnableWarmstarting;
#ifdef __S2_DS_ENABLE_PARAMS
        Real m_SolverLS_Log10_RelEpsilon;
        Real m_SolverLS_Log10_AbsEpsilon;
#endif
        Real m_SolverLS_RelEpsilon; //derived
        Real m_SolverLS_AbsEpsilon; //derived
        // NR Solver
        uint32 m_SolverNR_MaxIter;
#ifdef __S2_DS_ENABLE_PARAMS
        Real m_SolverNR_Log10_RelEpsilon;
        Real m_SolverNR_Log10_AbsEpsilon;
#endif
        Real m_SolverNR_RelEpsilon; //derived
        Real m_SolverNR_AbsEpsilon; //derived

        // Constraints
        enum EKNCSolverType { eKNCST_None = 0,
                              eKNCST_Reaction = 1,
                              cNumKNCST } m_ConstraintSolver_KNCST;
        enum ENDCSolverType { eNDCST_None = 0,
                              eNDCST_Reaction = 1,
                              cNumNDCST } m_ConstraintSolver_NDCST;
#ifdef __S3D_ENABLE_CONTACT
        // Contact Solver
        enum EPosContactSolverType { ePCST_None = 0,
                                     ePCST_Penalty = 1,
                                     ePCST_Alteration = 2,
                                     ePCST_Hack = 3,
                                     cNumPCST } m_ContactSolver_PCST;
        enum EVelContactSolverType { eVCST_None = 0,
                                     eVCST_Penalty = 1,
                                     eVCST_Reaction = 2,
                                     eVCST_Hack = 3,
                                     cNumVCST } m_ContactSolver_VCST;
        Real m_ContactSolver_DepthMargin;
        Real m_ContactSolver_Penalty_Ks;
        Real m_ContactSolver_Penalty_Kd;
        Real m_ContactSolver_Relaxation_Coeff;
        Real m_ContactSolver_Restitution_Coeff;
        Real m_ContactSolver_StatictFriction_Coeff;
        Real m_ContactSolver_DynamicFriction_Coeff;
        Real m_ContactSolver_MaxDepth;
        bool m_ContactSolver_BreakOnError;
#endif

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
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
        // Ground
        bool m_bEnableHackedGroundPlane;
#endif
        // Viz
        Flags32 m_GeoObjectDDF;

        inline Params()
        : m_pDSH(0)
          // Develop
        , m_bRun(true)
        , m_KinematicMode(eKM_Mouse)
          // Integrator
        , m_IntegratorType(eIT_QuasiImplicitEuler_v)
        , m_bEnableVariableDT(false)
        , m_FixedDT(0.001f)
        , m_TimeScale(1.0f)
          // Material
        , m_Density(1000.0f) //kg/m^3
          // Environment
        , m_AirDragCoeff(0.0f)
          //, m_Gravity(0,-9.8f,0)
        , m_Gravity(0,0,0)
          // Solver Internal
        , m_SolverLS_Type(eLSST_CG)
        , m_SolverLS_MaxIter(20)
        , m_SolverLS_MaxRestart(10)
        , m_SolverLS_bEnableWarmstarting(true)
#ifdef __S2_DS_ENABLE_PARAMS
        , m_SolverLS_Log10_RelEpsilon(-3.0f)
        , m_SolverLS_Log10_AbsEpsilon(-6.0f)
        , m_SolverLS_RelEpsilon( mal::Exp10( m_SolverLS_Log10_RelEpsilon ) )
        , m_SolverLS_AbsEpsilon( mal::Exp10( m_SolverLS_Log10_AbsEpsilon ) )
#else
        , m_SolverLS_RelEpsilon( mal::Exp10( -3.0f ) )
        , m_SolverLS_AbsEpsilon( mal::Exp10( -6.0f ) )
#endif
        , m_SolverNR_MaxIter(3)
#ifdef __S2_DS_ENABLE_PARAMS
        , m_SolverNR_Log10_RelEpsilon(-3.0f)
        , m_SolverNR_Log10_AbsEpsilon(-6.0f)
        , m_SolverNR_RelEpsilon( mal::Exp10( m_SolverNR_Log10_RelEpsilon ) )
        , m_SolverNR_AbsEpsilon( mal::Exp10( m_SolverNR_Log10_AbsEpsilon ) )
#else
        , m_SolverNR_RelEpsilon( mal::Exp10( -3.0f ) )
        , m_SolverNR_AbsEpsilon( mal::Exp10( -6.0f ) )
#endif
        , m_ConstraintSolver_KNCST(eKNCST_Reaction)
        , m_ConstraintSolver_NDCST(eNDCST_None)
#ifdef __S3D_ENABLE_CONTACT
          // Contact Solver (\note for a natural_scale of 1m)
        , m_ContactSolver_PCST(ePCST_Alteration)//ePCST_Penalty)
        , m_ContactSolver_VCST(eVCST_Reaction)//eVCST_Penalty)
        , m_ContactSolver_DepthMargin(0.01)
        , m_ContactSolver_Penalty_Ks(1000)
        , m_ContactSolver_Penalty_Kd(10)
        , m_ContactSolver_Relaxation_Coeff(0.1)
        , m_ContactSolver_Restitution_Coeff(0)
        , m_ContactSolver_StatictFriction_Coeff(0)//.15)
        , m_ContactSolver_DynamicFriction_Coeff(0)//.1)
        , m_ContactSolver_MaxDepth(1.0f)
        , m_ContactSolver_BreakOnError(false)
#endif
#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
          // Test
        , m_TestMode(eTM_Disabled)
        , m_TestNodes(eTN_All)
        , m_TestFraction(0)
        , m_bEnableHackedGroundPlane(false)
#endif
        // Viz
        , m_GeoObjectDDF( geo::eODDF_Shape | geo::eSDDF_DCR ) //geo::eODDF_Default | geo::eSDDF_Default )
            {
#ifdef __S2_DS_ENABLE_PARAMS
                RebuildRayleighCoeffsFromFreqs(this);
#else
                RebuildRayleighCoeffsFromFreqs();
#endif
            }

        finline void SetDSH( LeafDSH_Solid3D_FEM *p_dsh )
            {
                m_pDSH = p_dsh;
            }

#ifdef __S2_DS_ENABLE_PARAMS
        void SetName( const char *name ) {}
        const char *GetName() const { return "default s3d::Params"; }
        static void InitArchetype( util::ArchetypeLibrary &al );
        static util::ArchetypeLibrary s_LeafDSH_Solid3D_FEM_Params_ArchetypeLibrary; //TEMP:
#endif

#ifdef __S2_DS_ENABLE_PARAMS
        static bool RebuildFEM( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            if( params.m_pDSH && params.m_pDSH->m_TES.IsValid() ) params.m_pDSH->m_TES.SetParams( params.m_FEM );
            return true;
        }
        static bool RebuildMass( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            if( params.m_pDSH ) params.m_pDSH->RebuildMass();
            return true;
        }
        static bool RebuildAirDrag( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            if( params.m_pDSH ) params.m_pDSH->RebuildAirDrag();
            return true;
        }
        static bool RebuildSolverLS_Epsilon( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            params.m_SolverLS_RelEpsilon = mal::Exp10( params.m_SolverLS_Log10_RelEpsilon );
            params.m_SolverLS_AbsEpsilon = mal::Exp10( params.m_SolverLS_Log10_AbsEpsilon );
            return true;
        }
        static bool RebuildSolverLS_MaxRestart( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            if( params.m_pDSH ) params.m_pDSH->m_ILSS_GMRES.SetMaxRestarts( params.m_SolverLS_MaxRestart );
            return true;
        }
        static bool RebuildSolverNR_Epsilon( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            params.m_SolverNR_RelEpsilon = mal::Exp10( params.m_SolverNR_Log10_RelEpsilon );
            params.m_SolverNR_AbsEpsilon = mal::Exp10( params.m_SolverNR_Log10_AbsEpsilon );
            return true;
        }
        static bool RebuildRayleighCoeffsFromFreqs( IArchetypeInstance *p_instance )
        {
            Params &params( *static_cast<Params*>(p_instance) );
            Mat2x2 Omega;
            Omega(0,0) = mal::Rcp(params.m_FEM.m_RayleighFreq1); Omega(0,1) = params.m_FEM.m_RayleighFreq1;
            Omega(1,0) = mal::Rcp(params.m_FEM.m_RayleighFreq2); Omega(1,1) = params.m_FEM.m_RayleighFreq2;
            Vec2 zeta( params.m_FEM.m_DampingRatio, params.m_FEM.m_DampingRatio );
            params.m_FEM.m_RayleighCoeff = mal::Inverse(Omega) * (2*zeta);
            RebuildFEM( p_instance ); //\note propagate changes to FEM
            return true;
        }
#  ifdef __S3D_S2_DS_FEM_ENABLE_TEST
        static bool RebuildTest( IArchetypeInstance *p_instance )
            {
                Params &params( *static_cast<Params*>(p_instance) );
                if( params.m_pDSH ) params.m_pDSH->InitTest(); //\todo This resets the test completely... consider doing it only when TestMode actually changes
                return true;
            }
#  endif
#else
        void RebuildRayleighCoeffsFromFreqs()
        {
            Mat2x2 Omega;
            Omega(0,0) = mal::Rcp(m_FEM.m_RayleighFreq1); Omega(0,1) = m_FEM.m_RayleighFreq1;
            Omega(1,0) = mal::Rcp(m_FEM.m_RayleighFreq2); Omega(1,1) = m_FEM.m_RayleighFreq2;
            Vec2 zeta( m_FEM.m_DampingRatio, m_FEM.m_DampingRatio );
            m_FEM.m_RayleighCoeff = mal::Inverse(Omega) * (2*zeta);
        }
#endif
    };

private:
    Params m_Params;

    const geo::TetSolidShape3* m_pMeshS;

    unsigned int m_NumNodes;
    unsigned int m_NumElements;
    ms::fem::EditableTetrahedronElementSet3 m_TES;
    ms::IForce3::ICache* m_pTESC;

    //\name Continuous
    //@{
    Vec3 *m_vecPos0;
    //@}

    // SE tmp
    Vec3 *m_vecForces;

    // LS tmp
    Vec3 *m_LS_vec_b;
    Vec3 *m_LS_vec_y;
    Real *m_LS_real_array_tmp0;
    Real *m_LS_real_array_tmp1;
    Real *m_LS_real_array_tmp2;

    // NL-IE tmp
    Vec3 *m_IE_vec_x_k;
    Vec3 *m_IE_vec_v_k;

    //\name Derived params
    //@{
    Real m_TotalMass;
    Real m_DensityPerVolume;
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

    geo::bv::AABB3 m_AABB0;

#ifdef __S3D_ENABLE_KNC
    ms::KinematicNodeConstraintSet3 m_KNCS;
#endif

#ifdef __S3D_ENABLE_CONTACT_D2K
    //TEMP, this classes  should be in an ms::ContactConstraintSetD2K or something
    struct ContactPointConstraintD2K
    {
        /*! Specific D2K contact data, may include baricentric coords,
         whatever... We ONLY store info for the Solid, the Geom is
         Kinematic and therefore we don't need its specific contact
         info, only the part required to compute response on the
         Solid. */
        ContactPointConstraintD2K( const geo::PointOnFeature& pof1, const Vec3& normal, Real depth, Real radius )
        : m_POF1(pof1)
        , m_Normal(normal)
        , m_Depth(depth)
        , m_Radius(radius)
        , m_LambdaN(0)
        , m_IsActive(true)
            {}
        geo::PointOnFeature m_POF1;
        Vec3 m_Normal;
        Real m_Depth;
        Real m_Radius;
        Real m_LambdaN; //Normal impulse
        bool m_IsActive;
    };

    struct ContactConstraintD2K
    {
        /*\todo We'd like to have a FIXED cMaxCPC, but it's not
                realistic unless we come up with a good contact
                reduction scheme, and it's MUCH MORE complicated than
                the rigid case, because CP inside the "convex hull"
                are ALSO relevant here... so by now we store them all!
                enum EConstants { cMaxCPC = 8 };
                ContactPointConstraintD2K m_vecCPC[cMaxCPC];

        \todo Actually, if individual CPC persistence and vanishing is
              required, an allocator with cheaper new/delete and
              iteration could be algorithmically better (eg: a dynamic
              iterable pool), BUT, depending on how we use the CPC in
              the solver, if >1 passes are required, an array will
              surely be faster due to contiguous memory access.
        */
        inline void Create() { Destroy(); m_TmpCD.Clear(); }
        inline void Destroy() { m_Age = 0; m_CC.Invalidate(); Clear(); }

        inline void Clear() { m_vecCPC.clear(); }

        inline void Reset( const geo::np::ContactData3& cd, Real epsilon_length )
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
        inline void Persist( const geo::np::ContactData3& cd, Real epsilon_length, Real epsilon_direction )
            {
                // \todo Match and persist...
                m_Age++;
                Reset(cd,epsilon_length);
            }

        finline bool IsEmpty() const { return m_vecCPC.empty(); }
        finline unsigned int GetNumCPC() const { return m_vecCPC.size(); }
        finline ContactPointConstraintD2K& GetCPC( int i ) { return m_vecCPC[i]; }
        finline const ContactPointConstraintD2K& GetCPC( int i ) const { return m_vecCPC[i]; }

        uint32 m_Age;
        geo::np::ContactCache3 m_CC;
        std::vector<ContactPointConstraintD2K> m_vecCPC;

        geo::np::ContactData3 m_TmpCD; //TEMPORAL: Single CD stored to Viz later, otherwise local var
    };
    typedef util::GPoolDA<ContactConstraintD2K> PoolCCD2K;
    PoolCCD2K m_poolCCD2K;

    /*\todo These methods COULD be part of ContactConstraintSetD2K
      and, if generalized, of IConstraintSet (sparse constraints on
      DOF), which would include
      KNC,Contacts,Non-Collapse(DAPD),DistanceConstraintSet,etc...
    */
    //@{
    bool CheckError_Contact() const;
    void ApplyContacts_PositionAlteration( Vec3* vec_pos ) const;
    void ApplyContacts_PositionAlteration_Smoothed( Vec3* vec_pos ) const;
    void ApplyContacts_Penalty( Vec3* vec_forces ) const; //IN-integration
    void ApplyContacts_Reaction( Vec3* vec_v ) const; //IN-integration
    void ApplyContacts_Hack(); //POST-integration
    void ApplyContacts_Friction_dv( const Vec3* vec_vel0, Vec3* vec_vel1 ) const; //POST-integration, FAILS

    // Normal/Friction impulses, Normal uses "reaction constraints"
    void ApplyContacts_NormalImpulse( Vec3* vec_v, Real factor_v_to_impulse ); //Sets cpc.m_LambdaN and applies it to v, factor converts to impulses
    void ApplyContacts_FrictionImpulse( Vec3* vec_v, Real factor_impulse_to_v ) const;

    void UpdateContacts_ActiveSet( const Vec3* vec_vel, const Vec3* vec_acc );
    //@}
#endif //__S3D_ENABLE_CONTACT_D2K

private:
    //\name Dynamics
    //@{
    void FixedStep( Real dt );

    void Step_SymplecticEuler( Real dt );
    void Step_QuasiImplicitEuler_v( Real dt );
    //void Step_QuasiImplicitEuler_dx( Real dt );
    void Step_FullyImplicitEuler_dx( Real dt );
    //@}

    void UpdateKNC( Real dt );

    //\name LS internal methods
    //@{
    void LS_ZeroKNC( Vec3* v ) const; //Set to 0 all node-entries corresponding to KNC
    void LS_A_times_d( Real* q, Vec3* vq,
                       Real factor_M, Real factor_K,
                       const Real* d, const Vec3* vd ) const;
    void LS_Solve( Real factor_M, Real factor_K, Vec3* vy, Vec3* vb, Real& ls_rel_prec, Real& ls_abs_prec, int& ls_iter ) const;
    //@}

    //\name ILSS
    //@{
    void ILSS_D_Av( Real* Av, const Real* v ) const;
    void ILSS_D_ApplyConstraints( Real* v ) const;
    mutable ns::ILSS_CG m_ILSS_CG;
    mutable ns::ILSS_CR m_ILSS_CR;
    mutable ns::ILSS_SYMMQMR m_ILSS_SQMR;
    mutable ns::ILSS_MINRES m_ILSS_MINRES;
    mutable ns::ILSS_GMRES m_ILSS_GMRES;
    mutable Real m_ILSS_FactorM, m_ILSS_FactorK;
    //@}

    //\name Energy
    //@{
    void RecomputeEnergy();
    inline Real GetElasticEnergy() const { return m_ElasticEnergy; }
    inline Real GetPotentialEnergy() const { return m_PotentialEnergy; }
    inline Real GetKineticEnergy() const { return m_KineticEnergy; }
    //@}

    //\name Rebuild Params
    //@{
    void RebuildAirDrag();
    void RebuildMass();
    //@}

    unsigned int FindClosestNode( const Vec3 &pos ) const;
    Vec3 ComputePosCoM() const;
    Mat3x3 ComputeRotCoM() const;

    //\name Edition
    //@{
    void SetPosCoM( const Vec3 &pos );
    void SetRotCoM( const Mat3x3 &rot );
    //@}

private:
    //! Draw Debug Flags
    enum EDDF { eDraw_Nothing        = 0,
                eDraw_Mesh           = 1<<0,
                eDraw_KNC            = 1<<1,
                eDraw_Degenerate     = 1<<2,
                eDraw_Corotational   = 1<<3,
                eDraw_CoM            = 1<<4,
                eDraw_Contacts       = 1<<5,
                eDraw_Forces         = 1<<6,
                eDraw_Trajectory     = 1<<7,
                eDraw_DoC            = 1<<8,
                eDraw_Everything     = 0xFFFFFFFF,
                eDraw_Default        = eDraw_Mesh | eDraw_KNC | eDraw_Contacts | eDraw_Degenerate /* | eDraw_Corotational | eDraw_CoM | eDraw_Eigenbasis | eDraw_Forces | eDraw_Trajectory */ | eDraw_DoC
    };
    Flags32 m_DDF;

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
private:
    friend class Params;
    void InitTest();
    void UpdateTest( Real dt );
    void ApplyHackedGroundPlane();
    //Array of test-owned KNC
    std::vector<ms::KinematicNodeConstraintSet3::knc_id_type> m_vecTestKNC;
    std::vector<Vec3> m_vec_TM_Random_TargetPos;
#endif

#ifdef __S2_DS_ENABLE_STATS
private:
    struct Stats
    {
        uint32 m_SolverLS_NumIter;
        float32 m_SolverLS_RelPrec;
        float32 m_SolverLS_AbsPrec;
        uint32 m_SolverNR_NumIter;
        float32 m_SolverNR_RelPrec;
        float32 m_SolverNR_AbsPrec;
        // Last-step stuff
        uint32 m_VariableStep_LS_NumIter;
        uint32 m_VariableStep_NR_NumIter;
        float32 m_VariableStep_NR_AbsPrec;
        float32 m_VariableStep_Energy;
        float32 m_VariableStep_Duration;
        util::Chronos m_Clock;
        Stats() { BeginStep(); m_Clock.ResetTime(); } //Reset all stats
        void BeginStep() { m_SolverLS_NumIter = 0; m_SolverLS_RelPrec = 0; m_SolverLS_AbsPrec = 0;
                           m_SolverNR_NumIter = 0; m_SolverNR_RelPrec = 0; m_SolverNR_AbsPrec = 0;
                           m_VariableStep_LS_NumIter = 0; m_VariableStep_NR_NumIter = 0; m_VariableStep_NR_AbsPrec = 0; m_VariableStep_Energy = 0;
                           m_VariableStep_Duration = m_Clock.GetTime(); }
        void EndStep() { m_VariableStep_Duration = m_Clock.GetTime() - m_VariableStep_Duration; }
    };
    mutable Stats m_Stats;
#endif

    //TEMPORAL
    friend class SS_AggregateDSH_Basic;
    void NR_Init( Real dt );
    void NR_End( Real dt );
    void NR_PreStep( Real dt );
    void NR_MidStep( Real dt );
    void NR_PostStep( Real dt );
    uint32 m_NR_nr_k;
};

} } // namespace S2::ds

#endif // S2_DS_DSH_LEAF_SOLID3D_FEM_H
