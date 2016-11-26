#include "LeafDSH_Solid3D_FEM.h"

#include <Saphyre2/ds/dsh/IDynamicSystemHierarchy.h>
#include <Saphyre2/ds/dsh/IGeom.h>

#include <Mal/GConversion.h>
#include <Mal/GRandom.h>

#include <Saphyre2/ns/real_array.h>
#include <boost/bind.hpp> //for ILSS

#include <Geo/bp/Proxy.h>
#include <Geo/bp/Pair.h>
#include <Geo/geo_properties.h>

#include <util/ItemStream.h>

#define __ENABLE_PERSISTENT_CONTACTS
#define __ENABLE_CONTACT_EFFECTIVE_MASS //\todo Works BETTERER than non-consistent version, should be default
#define __ENABLE_CONTACT_ACTIVE_SET_IGNORE_ACCELERATION //\todo Mostly sure this is required, because gravity ALWAYS pulls downwards...
//#define __ENABLE_ILSS_CHECK_INDEFINITE //TEMP: shouldn't be req,
//#define __USE_CG_PRECONDITIONER \todo Not properly implemented, I fear...

//#define __S2_DS_ENABLE_TRACE
#ifdef __S2_DS_ENABLE_TRACE
#include <util/SimpleTracer.h>
#include <Mal/GRandom.h>
#  define S2_DS_TRACE_INIT() UTIL_TRACE_INIT()
#  define S2_DS_TRACE_SHUTDOWN() UTIL_TRACE_SHUTDOWN()
#  define S2_DS_TRACE_BEGIN_SLICE(x) UTIL_TRACE_BEGIN_SLICE(x)
#  define S2_DS_TRACE_END_SLICE(b_flush) UTIL_TRACE_END_SLICE(b_flush)
#  define S2_DS_TRACE_BEGIN_SCOPE(name) UTIL_TRACE_BEGIN_SCOPE(name)
#  define S2_DS_TRACE_END_SCOPE() UTIL_TRACE_END_SCOPE()
#  define S2_DS_TRACE_LOCAL(name,value) UTIL_TRACE_LOCAL(name,value)
#  define S2_DS_TRACE_LOCAL_ARRAY(name,vec_value,num_values) UTIL_TRACE_LOCAL_ARRAY(name,vec_value,num_values)
#  define S2_DS_TRACE_GLOBAL(name,value) UTIL_TRACE_GLOBAL(name,value)
#  define S2_DS_TRACE_GLOBAL_ARRAY(name,vec_value,num_values) UTIL_TRACE_GLOBAL_ARRAY(name,vec_value,num_values)
#else
#  define S2_DS_TRACE_INIT()
#  define S2_DS_TRACE_SHUTDOWN()
#  define S2_DS_TRACE_BEGIN_SLICE(x)
#  define S2_DS_TRACE_END_SLICE(b_flush)
#  define S2_DS_TRACE_BEGIN_SCOPE(name)
#  define S2_DS_TRACE_END_SCOPE()
#  define S2_DS_TRACE_LOCAL(name,value)
#  define S2_DS_TRACE_LOCAL_ARRAY(name,vec_value,num_values)
#  define S2_DS_TRACE_GLOBAL(name,value)
#  define S2_DS_TRACE_GLOBAL_ARRAY(name,vec_value,num_values)
#endif

namespace S2 { namespace ds {

//TEMP: This is useful, consider putting it somewhere else
static uint32 MapNameToValue_Enum32( const char* name, const char** vec_names, uint32* vec_values, unsigned int num_values )
{
    std::string str_name(name);
    for( unsigned int i=0; i<num_values; i++ )
        if( str_name == std::string(vec_names[i]) )
            return vec_values[i];
    return vec_values[0];
}

LeafDSH_Solid3D_FEM::LeafDSH_Solid3D_FEM( uint32 uid, IDynamicSystemHierarchy *p_parent )
: GLeafDSH_Model< ms::ParticleSystem3D >(uid,p_parent)
, m_pMeshS(0)
, m_NumNodes(0)
, m_NumElements(0)
, m_pTESC(0)
, m_vecPos0(0)
  // SE tmp
, m_vecForces(0)
  // LS tmp
, m_LS_vec_b(0)
, m_LS_vec_y(0)
, m_LS_real_array_tmp0(0)
, m_LS_real_array_tmp1(0)
, m_LS_real_array_tmp2(0)
  // NL-IE tmp
, m_IE_vec_x_k(0)
, m_IE_vec_v_k(0)
  // Derived params
, m_DensityPerVolume(0)
, m_AirDragPerSecond(0)
  // State
, m_TotalTime(0)
, m_AccTime(0)
  // Energy
, m_ElasticEnergy(0)
, m_PotentialEnergy(0)
, m_KineticEnergy(0)
#ifdef __S3D_ENABLE_CONTACT_D2K //TEMP, this should be in an ms::ContactConstraintSetD2K or something
, m_poolCCD2K(8)
#endif
  // ILSSS
, m_ILSS_FactorM(0)
, m_ILSS_FactorK(0)
  // DDF
, m_DDF( eDraw_Default )
{
    m_Params.SetDSH( this );
#ifdef __S2_DS_ENABLE_PARAMS
    if( Params::s_LeafDSH_Solid3D_FEM_Params_ArchetypeLibrary.IsEmpty() )
    {
        Params::InitArchetype( Params::s_LeafDSH_Solid3D_FEM_Params_ArchetypeLibrary );
    }
#endif
}

LeafDSH_Solid3D_FEM::~LeafDSH_Solid3D_FEM()
{
    if( m_vecPos0 ) ns::real_array::dealloc_GVec( m_vecPos0 );
    if( m_vecForces ) ns::real_array::dealloc_GVec( m_vecForces );

    if( m_LS_vec_b ) ns::real_array::dealloc_GVec( m_LS_vec_b );
    if( m_LS_vec_y ) ns::real_array::dealloc_GVec( m_LS_vec_y );
    if( m_LS_real_array_tmp0 ) ns::real_array::dealloc( m_LS_real_array_tmp0 );
    if( m_LS_real_array_tmp1 ) ns::real_array::dealloc( m_LS_real_array_tmp1 );
    if( m_LS_real_array_tmp2 ) ns::real_array::dealloc( m_LS_real_array_tmp2 );

    if( m_IE_vec_x_k ) ns::real_array::dealloc_GVec( m_IE_vec_x_k );
    if( m_IE_vec_v_k ) ns::real_array::dealloc_GVec( m_IE_vec_v_k );
    if( m_pTESC ) m_TES.DestroyCache( m_pTESC );
}

bool LeafDSH_Solid3D_FEM::Create( const ParamIt &pit )
{
    m_Params.m_FixedDT = pit.Find("fixed_dt").Get<Real>();

    // Get TetSolidShape3 from already-bound GO
    GEO_ASSERT( m_pGO != 0 );
    const geo::IShape* pShape( m_pGO->GetShapeInterface() );
    GEO_ASSERT( pShape->GetType() == geo::eShape_TetSolid3 );
    m_pMeshS = static_cast<const geo::TetSolidShape3*>(pShape);

    // Init particle system with fake uniform mass
    m_NumNodes = m_pMeshS->GetNumV();
    m_TotalMass = Real(666);
    m_Model.SetDscr( m_NumNodes, m_TotalMass );
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_Model.SetPos( it_node, m_pMeshS->GetVecDefaultSDOF()[it_node] );
        m_Model.SetVel( it_node, Vec3::Zero() );
    }
    // Init previous configuration C_i-1 = C_0
    m_vecPos0 = ns::real_array::alloc_GVec<Real,3>( m_NumNodes );
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_vecPos0[it_node] = m_Model.GetPos(it_node);

    // Get material params
    m_Params.m_Density = pit.Find("density").Get<Real>();
    m_Params.m_FEM.m_YoungModulus = pit.Find("young_modulus").Get<Real>();
    m_Params.m_FEM.m_PoissonRatio = pit.Find("poisson_ratio").Get<Real>();
    m_Params.m_FEM.m_DampingRatio = pit.Find("damping_ratio").Get<Real>();
    m_Params.m_FEM.m_PlasticYield = pit.Find("plastic_yield").Get<Real>();
    m_Params.m_FEM.m_PlasticMax = pit.Find("plastic_max").Get<Real>();
    m_Params.m_FEM.m_PlasticCreepPerSecond = pit.Find("plastic_creep_per_second").Get<Real>();

    // Create TetrahedronElementSet3
    m_NumElements = m_pMeshS->GetNumT();
    m_TES.BeginEdition( m_NumNodes, m_pMeshS->GetVecDefaultSDOF(), m_NumElements, m_Params.m_FEM );
    for( unsigned int it_tet=0; it_tet<m_NumElements; it_tet++ )
    {
        m_TES.AddElement( m_pMeshS->T_VID(it_tet,0),
                          m_pMeshS->T_VID(it_tet,1),
                          m_pMeshS->T_VID(it_tet,2),
                          m_pMeshS->T_VID(it_tet,3) );
    }
    m_TES.EndEdition();

    m_pTESC = m_TES.CreateCache( m_Model.GetVecInvMass() );

    // Rebuild everything except LT
#ifdef __S2_DS_ENABLE_PARAMS
    Params::RebuildRayleighCoeffsFromFreqs(&m_Params);
#else
    m_Params.RebuildRayleighCoeffsFromFreqs();
#endif
    RebuildMass();
    RebuildAirDrag();

    // Init CoM position
    if( pit.Find("pos_com").IsValid() ) SetPosCoM( pit.Find("pos_com").Get<Vec3>() );
    if( pit.Find("rot_com").IsValid() ) SetRotCoM( pit.Find("rot_com").Get<Mat3x3>() );

    // Init object-size based params
    m_pMeshS->ComputeBVD( m_AABB0, geo::Transform3::Identity(), 0 );

    // SE temp elastic force vector
    //\todo ALLOC with SIMD alignment, ideally in a SINGLE BLOCK
    m_vecForces = ns::real_array::alloc_GVec<Real,3>( m_NumNodes );
    memset( m_vecForces, 0, sizeof(Vec3) * m_NumNodes );

    // IE temporaries
    m_LS_vec_b = ns::real_array::alloc_GVec<Real,3>( m_NumNodes );
    m_LS_vec_y = ns::real_array::alloc_GVec<Real,3>( m_NumNodes );
    m_LS_real_array_tmp0 = ns::real_array::alloc( 3*m_NumNodes );
    m_LS_real_array_tmp1 = ns::real_array::alloc( 3*m_NumNodes );
    m_LS_real_array_tmp2 = ns::real_array::alloc( 3*m_NumNodes );
    m_IE_vec_x_k = ns::real_array::alloc_GVec<Real,3>( m_NumNodes );
    m_IE_vec_v_k = ns::real_array::alloc_GVec<Real,3>( m_NumNodes );

    // ILSS init
    m_ILSS_CG.Init( 3*m_NumNodes );
    m_ILSS_CR.Init( 3*m_NumNodes );
    m_ILSS_SQMR.Init( 3*m_NumNodes );
    m_ILSS_MINRES.Init( 3*m_NumNodes );
    m_ILSS_GMRES.Init( 3*m_NumNodes, m_Params.m_SolverLS_MaxRestart );

#ifdef __S3D_ENABLE_BVH
    if( m_pMeshS->GetBVH() == 0 )
    {
        geo::BVH_TetSolidShape3* pBVH = new geo::BVH_TetSolidShape3;
        uint32 num_entries( m_pMeshS->GetDCR() != 0 ? m_pMeshS->GetDCR()->m_NumElements : m_pMeshS->GetNumT() ); //TEMP: If there's a DCR, only use those P in the BVH
        pBVH->Rebuild_TopDown( num_entries,
                               boost::bind<void>( &geo::GEBV_TetSolidShape3_E<geo::BVH_TetSolidShape3::entry_index_type,geo::BVH_TetSolidShape3::bv_type>,
                                                  m_pMeshS, Transform3::Identity(), m_pMeshS->GetVecDefaultSDOF(),
                                                  _1, _2) );
        const_cast<geo::TetSolidShape3*>(m_pMeshS)->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH );
    }
#endif

    return true;
}

void LeafDSH_Solid3D_FEM::Step( Real dt )
{
    if( !m_Params.m_bRun ) return;

    Real scaled_dt( m_Params.m_TimeScale * dt );
    m_TotalTime += scaled_dt;

#ifdef __S2_DS_ENABLE_STATS
    m_Stats.BeginStep();
#endif

    //\todo FixedDT THIS SHOULD BE PROBABLY HANDLED IN THE DSH OR THE SIMULATION SCHEME, NOT HERE
    m_AccTime += scaled_dt;
    Real scaled_fixed_dt( m_Params.m_TimeScale * m_Params.m_FixedDT );
    while( m_AccTime > scaled_fixed_dt )
    {
        FixedStep( scaled_fixed_dt );
        m_AccTime -= scaled_fixed_dt;

        // Acc VarStep stats
        S2_DS_STAT_ADD( m_Stats.m_VariableStep_LS_NumIter, m_Stats.m_SolverLS_NumIter );
        S2_DS_STAT_ADD( m_Stats.m_VariableStep_NR_NumIter, m_Stats.m_SolverNR_NumIter );
        S2_DS_STAT_ADD( m_Stats.m_VariableStep_NR_AbsPrec, m_Stats.m_SolverNR_AbsPrec );
        S2_DS_STAT_SET( m_Stats.m_VariableStep_Energy, float32( m_Model.ComputeKineticEnergy()
                                                                + m_Model.ComputePotentialEnergy( m_Params.m_Gravity )
                                                                + GetElasticEnergy() ) / 1000.0f );
    }

#ifdef __S2_DS_ENABLE_STATS
    m_Stats.EndStep();
#endif
}

void LeafDSH_Solid3D_FEM::SyncGO()
{
    if( 0 != base_type::m_pGO )
    {
        //\note Assumes that transform == identity, as all
        //information is in the SDOF
        m_pGO->SetTransform( Transform3::Identity() );
        m_Model.GetDOF( static_cast< geo::GObjectSDOF<3,Vec3>* >(m_pGO)->GetVecDOF_WriteOnly() );
#ifdef __S3D_ENABLE_BVH
        //\todo Use more accurate GEBV: BDOP or BSlab
        m_pMeshS->GetBVH()->Refit( boost::bind<void>( &geo::GEBV_TetSolidShape3_E<geo::BVH_TetSolidShape3::entry_index_type,geo::BVH_TetSolidShape3::bv_type>,
                                                      m_pMeshS, geo::Transform3::Identity(), static_cast< const geo::GObjectSDOF<3,Vec3>* >(m_pGO)->GetVecSDOF(), _1, _2) );
#endif
    }
}

void LeafDSH_Solid3D_FEM::RecomputeBV()
{
    //\todo Use up-to-date BVH instead of recomputing in O(n)!!
    if( 0 != base_type::m_pGO && 0 != base_type::m_pBV )
    {
        base_type::m_pGO->ComputeBV( *base_type::m_pBV );
        //\todo Maybe use polimorphyc BV::Extend()...
        const Real cNearDist(0.25f); //\todo this should come from np::Context::m_Stochastic_Near/FarDist, probably a factor of it
        switch( base_type::m_pBV->GetType() )
        {
        case geo::bv::eBV_AABB3: static_cast<geo::bv::AABB3*>(base_type::m_pBV)->Extend( cNearDist ); break;
        case geo::bv::eBV_Sphere3: static_cast<geo::bv::Sphere3*>(base_type::m_pBV)->Extend( cNearDist ); break;
        case geo::bv::eBV_LSS3: static_cast<geo::bv::LSS3*>(base_type::m_pBV)->Extend( cNearDist ); break;
        default: DS_ASSERT(false); break;
        }
    }
}

bool LeafDSH_Solid3D_FEM::Edit( const ParamIt &pit )
{
    // CoM pos/rot
    ParamIt pit2 = pit.Find("pos_com");
    if( pit2.IsValid() ) SetPosCoM( pit2.Get<Vec3>() );
    pit2 = pit.Find("rot_com");
    if( pit2.IsValid() ) SetRotCoM( pit2.Get<Mat3x3>() );
    return true;
}

bool LeafDSH_Solid3D_FEM::Internal( const ParamIt &pit, ReturnStream &rets )
{
#ifdef __S3D_ENABLE_KNC
    ParamIt create_kpc_it = pit.Find("create_kpc");
    if( create_kpc_it.IsValid() )
    {
        // Find closest node to solid_pos and create a KNC, knc_eid is a pointer to the new KNC
        DS_LOG_WARNING("---------------------------------------------------------------- CREATING KPC!!");

        //\todo if we ever support non-node constraints, kpc_eid could
        //be a ptr to a generic ConstraintPoint class with
        //specializations Node and Generic... whatever...
        unsigned int closest_nid = FindClosestNode( create_kpc_it.GetSubItem().Find("solid_pos").Get<Vec3>() );

        // Return existing or new KNC, properly setting new params
        ms::KinematicNodeConstraintSet3::knc_id_type knc_id = m_KNCS.FindKNC( closest_nid );
        if( 0 == knc_id )
            knc_id = m_KNCS.AddKNC( closest_nid,
                                    create_kpc_it.GetSubItem().Find("world_pos").Get<Vec3>(),
                                    create_kpc_it.GetSubItem().Find("world_vel").Get<Vec3>() );
        else
            m_KNCS.ResetKNC( knc_id,
                             create_kpc_it.GetSubItem().Find("world_pos").Get<Vec3>(),
                             create_kpc_it.GetSubItem().Find("world_vel").Get<Vec3>() );

        machine_uint_type kpc_eid( reinterpret_cast<machine_uint_type>( knc_id ) );

        // Return new KPC data
        rets.BeginComplex( 666, eRet_Internal );
        {
            rets.WriteItem( "uid", pit.Find("uid") );
            rets.BeginComplex( "create_kpc", eRet_Internal );
            {
                //DS_LOG_WARNING( "Creating KNC %lld", kpc_eid );
                rets.WriteItem( "kpc_uid", create_kpc_it.GetSubItem().Find("kpc_uid") );
                rets.Write( "kpc_eid", kpc_eid );
            }
            rets.EndComplex();
        }
        rets.EndComplex();
    }

    ParamIt edit_kpc_it = pit.Find("edit_kpc");
    if( edit_kpc_it.IsValid() )
    {
        machine_uint_type kpc_eid = edit_kpc_it.GetSubItem().Find("kpc_eid").Get<machine_uint_type>();
        //DS_LOG_WARNING( "Editing KNC %lld", kpc_eid );
        ms::KinematicNodeConstraintSet3::knc_id_type knc_id = reinterpret_cast<ms::KinematicNodeConstraintSet3::knc_id_type>( kpc_eid );
        ParamIt pos_it = edit_kpc_it.GetSubItem().Find("world_pos");
        if( pos_it.IsValid() )
        {
            //DS_LOG_WARNING( "Editing KNC Pos %lld", kpc_eid );
            m_KNCS.SetPosKNC( knc_id, pos_it.Get<Vec3>() );
        }
        else
        {
            ParamIt vel_it = edit_kpc_it.GetSubItem().Find("world_vel");
            if( vel_it.IsValid() )
            {
                //DS_LOG_WARNING( "Editing KNC Vel %lld", kpc_eid );
                m_KNCS.SetVelKNC( knc_id, vel_it.Get<Vec3>() );
            }
        }
        // No data to be returned
    }

    ParamIt destroy_kpc_it = pit.Find("destroy_kpc");
    if( destroy_kpc_it.IsValid() )
    {
        machine_uint_type kpc_eid = destroy_kpc_it.GetSubItem().Find("kpc_eid").Get<machine_uint_type>();
        //DS_LOG_WARNING( "Destroying KNC %lld", kpc_eid );
        m_KNCS.RemoveKNC( reinterpret_cast<ms::KinematicNodeConstraintSet3::knc_id_type>( kpc_eid ) );
        // No data to be returned
    }
#endif

#define __ENABLE_PARAMS_EXTRA
#ifdef __ENABLE_PARAMS_EXTRA
    ParamIt spe_it = pit.Find("set_params_extra");
    if( spe_it.IsValid() )
    {
        m_Params.m_Gravity = spe_it.GetSubItem().Find("gravity").SafeGet<Vec3>( m_Params.m_Gravity ); //IArchetypeInstance::NTPF_Ignore );
        if( spe_it.GetSubItem().Find("mm").IsValid() )
        {
            const char *vec_names_mm[] = { "L",
                                           "C_WRP", "C_LCM", "C_LCMH", "C_CCM",
                                           "H_LCM", "H_CCM", "H_NH0", "H_NH1" };
            uint32 vec_values_mm[] = { ms::fem::Params::eMM_L,
                                       // Corotational
                                       ms::fem::Params::eMM_C_WRP,
                                       ms::fem::Params::eMM_C_LCM,
                                       ms::fem::Params::eMM_C_LCMH,
                                       ms::fem::Params::eMM_C_CCM,
                                       // Hyperelastic
                                       ms::fem::Params::eMM_H_LCM,
                                       ms::fem::Params::eMM_H_CCM,
                                       ms::fem::Params::eMM_H_NH_C0,
                                       ms::fem::Params::eMM_H_NH_C1 };
            m_Params.m_FEM.m_MM = (S2::ms::fem::Params::EMaterialMethod)MapNameToValue_Enum32( spe_it.GetSubItem().Find("mm").Get<String32>().GetStr(),
                                                                                               vec_names_mm, vec_values_mm, ms::fem::Params::cNumMM );
            Params::RebuildFEM( &m_Params );
        }
        if( spe_it.GetSubItem().Find("rm").IsValid() )
        {
            const char *vec_names_rm[] = { "Id", "QR", "MSVD", "DAPD" };
            uint32 vec_values_rm[] = { ms::fem::Params::eRM_Id,
                                       ms::fem::Params::eRM_QR,
                                       ms::fem::Params::eRM_MSVD,
                                       ms::fem::Params::eRM_DAPD };
            m_Params.m_FEM.m_RM = (S2::ms::fem::Params::ERotationMethod)MapNameToValue_Enum32( spe_it.GetSubItem().Find("rm").Get<String32>().GetStr(),
                                                                                               vec_names_rm, vec_values_rm, ms::fem::Params::cNumRM );
            Params::RebuildFEM( &m_Params );
        }
        if( spe_it.GetSubItem().Find("dm").IsValid() )
        {
            const char *vec_names_dm[] = { "T", "E", "I", "N" };
            uint32 vec_values_dm[] = { ms::fem::Params::eDM_Truncated,
                                       ms::fem::Params::eDM_Exact,
                                       ms::fem::Params::eDM_Exact_Inv,
                                       ms::fem::Params::eDM_Numerical };
            m_Params.m_FEM.m_DM = (S2::ms::fem::Params::EDifferentialsMethod)MapNameToValue_Enum32( spe_it.GetSubItem().Find("dm").Get<String32>().GetStr(),
                                                                                                    vec_names_dm, vec_values_dm, ms::fem::Params::cNumDM );
            Params::RebuildFEM( &m_Params );
        }
        if( spe_it.GetSubItem().Find("integrator").IsValid() )
        {
            const char *vec_names_it[] = { "SE", "QIE v", "QIE dx", "FIE dx", "FIQS dx" };
            uint32 vec_values_it[] = { (uint32)Params::eIT_SymplecticEuler,
                                       (uint32)Params::eIT_QuasiImplicitEuler_v,
                                       (uint32)Params::eIT_QuasiImplicitEuler_dx,
                                       (uint32)Params::eIT_FullyImplicitEuler_dx,
                                       (uint32)Params::eIT_FullyImplicitEuler_dx_QuasiStatic };
            m_Params.m_IntegratorType = (Params::EIntegratorType)MapNameToValue_Enum32( spe_it.GetSubItem().Find("integrator").Get<String32>().GetStr(),
                                                                                        vec_names_it, vec_values_it, Params::cNumIT );
        }
        m_Params.m_SolverNR_MaxIter = spe_it.GetSubItem().Find("nr_iter").SafeGet<uint32>( m_Params.m_SolverNR_MaxIter ); //IArchetypeInstance::NTPF_Ignore );
        m_Params.m_SolverNR_Log10_RelEpsilon = mal::Log10( spe_it.GetSubItem().Find("nr_rel_epsilon").SafeGet<float>( m_Params.m_SolverNR_RelEpsilon ) ); Params::RebuildSolverNR_Epsilon( &m_Params );
        if( spe_it.GetSubItem().Find("ls_solver").IsValid() )
        {
            const char *vec_names_ls[] = { "CG", "CR", "SQMR", "MINRES", "GMRES",
                                           "CGN", "CGS",
                                           "RI",
                                           "GS", "Jacobi", "LU" };
            uint32 vec_values_ls[] = { (uint32)Params::eLSST_CG,
                                       (uint32)Params::eLSST_CR,
                                       (uint32)Params::eLSST_SQMR,
                                       (uint32)Params::eLSST_MINRES,
                                       (uint32)Params::eLSST_GMRES,
                                       (uint32)Params::eLSST_CGN,
                                       (uint32)Params::eLSST_CGS,
                                       (uint32)Params::eLSST_RI,
                                       (uint32)Params::eLSST_GS,
                                       (uint32)Params::eLSST_Jacobi,
                                       (uint32)Params::eLSST_LU };
            m_Params.m_SolverLS_Type = (Params::ELinearSystemSolverType)MapNameToValue_Enum32( spe_it.GetSubItem().Find("ls_solver").Get<String32>().GetStr(),
                                                                                               vec_names_ls, vec_values_ls, Params::cNumLSST );
        }
        m_Params.m_SolverLS_MaxIter = spe_it.GetSubItem().Find("ls_iter").SafeGet<uint32>( m_Params.m_SolverLS_MaxIter ); //IArchetypeInstance::NTPF_Ignore );
        m_Params.m_SolverLS_MaxRestart = spe_it.GetSubItem().Find("ls_restarts").SafeGet<uint32>( m_Params.m_SolverLS_MaxRestart ); Params::RebuildSolverLS_MaxRestart( &m_Params );
        m_Params.m_SolverLS_Log10_RelEpsilon = mal::Log10( spe_it.GetSubItem().Find("ls_rel_epsilon").SafeGet<float>( m_Params.m_SolverLS_RelEpsilon ) ); Params::RebuildSolverNR_Epsilon( &m_Params );

        //Extra MIG2015
        {
            if( spe_it.GetSubItem().Find("cs_posm").IsValid() )
            {
                const char *vec_names_pcst[] = { "None", "Penalty", "Alteration", "Hack" };
                uint32 vec_values_pcst[] = { (uint32)Params::ePCST_None,
                                             (uint32)Params::ePCST_Penalty,
                                             (uint32)Params::ePCST_Alteration,
                                             (uint32)Params::ePCST_Hack };
                m_Params.m_ContactSolver_PCST = (Params::EPosContactSolverType)MapNameToValue_Enum32( spe_it.GetSubItem().Find("cs_posm").Get<String32>().GetStr(),
                                                                                                      vec_names_pcst, vec_values_pcst, Params::cNumPCST );
            }
            if( spe_it.GetSubItem().Find("cs_velm").IsValid() )
            {
                const char *vec_names_vcst[] = { "None", "Penalty", "Reaction", "Hack" };
                uint32 vec_values_vcst[] = { (uint32)Params::eVCST_None,
                                             (uint32)Params::eVCST_Penalty,
                                             (uint32)Params::eVCST_Reaction,
                                             (uint32)Params::eVCST_Hack };
                m_Params.m_ContactSolver_VCST = (Params::EVelContactSolverType)MapNameToValue_Enum32( spe_it.GetSubItem().Find("cs_velm").Get<String32>().GetStr(),
                                                                                                      vec_names_vcst, vec_values_vcst, Params::cNumVCST );
            }
            m_Params.m_ContactSolver_DepthMargin = spe_it.GetSubItem().Find( "cs_depth_margin" ).SafeGet<float32>( m_Params.m_ContactSolver_DepthMargin );
            m_Params.m_ContactSolver_Penalty_Ks = spe_it.GetSubItem().Find( "cs_penalty_ks" ).SafeGet<float32>( m_Params.m_ContactSolver_Penalty_Ks );
            m_Params.m_ContactSolver_Penalty_Kd = spe_it.GetSubItem().Find( "cs_penalty_kd" ).SafeGet<float32>( m_Params.m_ContactSolver_Penalty_Kd );
            m_Params.m_ContactSolver_Relaxation_Coeff = spe_it.GetSubItem().Find( "cs_relax_coeff" ).SafeGet<float32>( m_Params.m_ContactSolver_Relaxation_Coeff );
            m_Params.m_ContactSolver_StatictFriction_Coeff = spe_it.GetSubItem().Find( "cs_friction_s" ).SafeGet<float32>( m_Params.m_ContactSolver_StatictFriction_Coeff );
            m_Params.m_ContactSolver_DynamicFriction_Coeff = spe_it.GetSubItem().Find( "cs_friction_d" ).SafeGet<float32>( m_Params.m_ContactSolver_DynamicFriction_Coeff );
        }
    }
#endif

    return true;
}

void LeafDSH_Solid3D_FEM::DoViz( util::VizStream &vs ) const
{
    if( m_DDF.Test( eDraw_Mesh ) )
        geo::VizObject( m_pGO, vs, m_Params.m_GeoObjectDDF ); //\note we assume GO already Sync'ed

    // Highlight KNC
#ifdef __S3D_ENABLE_KNC
    if( m_DDF.Test( eDraw_KNC ) )
    {
        if( m_Params.m_KinematicMode != Params::eKM_Disabled )
        {
            for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
            {
                VIZ_POINT3( vs, it_knc->m_Pos, 10, Vec4f(1,1,1,1) );
                VIZ_POINT3( vs, m_Model.GetPos( it_knc->m_NID ), 10, Vec4f(1,0,0,1) );
                /*
                  char name[16];
                  sprintf( name, " %d ", it_knc->m_nid );
                  VIZ_POINT2_NAMED( vs, name, it_knc->m_Pos, 10, Vec4f(1,1,1,1) );
                */
            }
        }
    }
#endif

#ifdef __S3D_ENABLE_VIZ_DEGENERATE_ELEMENTS
    // Element viz: Degenerate or Extremely deformed elements
    if( m_DDF.Test( eDraw_Degenerate ) )
    {
        for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
        {
            // Gather VID
            unsigned int vec_nid[4];
            m_TES.GetNID( it_e, vec_nid[0], vec_nid[1], vec_nid[2], vec_nid[3] );
            Vec3 x0( m_Model.GetPos( vec_nid[0] ) );
            Vec3 x1( m_Model.GetPos( vec_nid[1] ) );
            Vec3 x2( m_Model.GetPos( vec_nid[2] ) );
            Vec3 x3( m_Model.GetPos( vec_nid[3] ) );
            float volume = ms::fem::TetrahedronElement3::Compute_Volume( x0, x1, x2, x3 );
            float rel_volume_deformation = (volume / m_TES.Volume( it_e )) - 1.0f;
            if( volume < 0 ) // Inverted
            {
                VIZ_TETRAHEDRON3( vs, x0, x1, x2, x3, Vec4f(0.5,0,0.5,0.25), util::eVizStyle_Solid );
            }
        }
    }
#endif

    // Show CoM
    if( m_DDF.Test( eDraw_CoM ) )
    {
        VIZ_POINT3( vs, ComputePosCoM(), 10, Vec4f(0,0,1,1) );
    }

    // Show trajectories
    if( m_DDF.Test( eDraw_Trajectory ) )
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            VIZ_SEGMENT3( vs, m_vecPos0[it_node], m_Model.GetPos(it_node), 2.0f, Vec4f(0.25,0.25,0.75,1) );

    // Show corot transforms
    if( m_DDF.Test( eDraw_Corotational ) )
    {
        for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
        {
            // Compute Transform from element 2 world
            Transform3 Te2w( m_TES.GetTransformE2W( m_pTESC, it_e, m_Model.GetVecPos() ) );
            // Element transform viz
            VIZ_TRANSFORM3( vs, Te2w, 0.33f );
        }
    }

    // Draw node forces
    if( m_DDF.Test( eDraw_Forces ) )
    {
#ifdef __S3D_ENABLE_SAVE_FORCES
        switch( m_Params.m_IntegratorType )
        {
        case Params::eIT_SymplecticEuler:
            for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
                VIZ_SEGMENT3( vs, m_Model.GetPos(it_node), m_Model.GetPos(it_node) + m_vecForces[it_node], 1.0f, Vec4f(0.25,0.25,1,1) );
            break;
        case Params::eIT_QuasiImplicitEuler_v:
        case Params::eIT_QuasiImplicitEuler_dx:
        case Params::eIT_FullyImplicitEuler_dx:
        case Params::eIT_FullyImplicitEuler_dx_QuasiStatic:
            for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
                VIZ_SEGMENT3( vs, m_Model.GetPos(it_node), m_Model.GetPos(it_node) + m_Model.GetVel(it_node), 1.0f, Vec4f(0,1,0,1) );
            break;
        default: break;
        }
#endif
    }

    if( m_DDF.Test( eDraw_DoC ) )
    {
        for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
        {
            ms::fem::TetrahedronElement3::DoC doc( m_TES.GetDoC( m_pTESC, it_e ) );
            if( doc.IsVF() )
            {
                unsigned int vec_nid[4];
                m_TES.GetNID( it_e, vec_nid[0], vec_nid[1], vec_nid[2], vec_nid[3] );
                Vec3 point_doc_V( m_Model.GetPos( vec_nid[doc.m_VIT0] ) );
                Vec3 point_doc_F( Real(1.0f/3.0f)*( m_Model.GetPos( vec_nid[ (doc.m_VIT0+1)%4 ] )
                                                    + m_Model.GetPos( vec_nid[ (doc.m_VIT0+2)%4 ] )
                                                    + m_Model.GetPos( vec_nid[ (doc.m_VIT0+3)%4 ] ) ) );
                VIZ_SEGMENT3( vs, point_doc_V, point_doc_F, 4.0f, Vec4f(1,0,1,1) );
                VIZ_POINT3( vs, point_doc_V, 20.0f, Vec4f(1,0,1,1) );
            }
            else if( doc.IsEE() )
            {
                unsigned int vec_nid[4];
                m_TES.GetNID( it_e, vec_nid[0], vec_nid[1], vec_nid[2], vec_nid[3] );
                Vec3 point_doc_V0( m_Model.GetPos( vec_nid[doc.m_VIT0] ) );
                Vec3 point_doc_V1( m_Model.GetPos( vec_nid[doc.m_VIT1] ) );
                VIZ_SEGMENT3( vs, point_doc_V0, point_doc_V1, 4.0f, Vec4f(1,0,1,1) );
            }
            //else doc.IsInvalid()
        }
    }

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
    switch( m_Params.m_TestMode )
    {
    case Params::eTM_Disabled: break;
    case Params::eTM_Random:
        {
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                for( unsigned int it_tp=0; it_tp < m_vec_TM_Random_TargetPos.size(); it_tp++ )
                {
                    VIZ_SEGMENT3( vs, m_pMeshS->V_Pos_0( it_tp ), m_vec_TM_Random_TargetPos[it_tp], 1, Vec4f(1,0,1,0.5) );
                    VIZ_POINT3( vs, m_vec_TM_Random_TargetPos[it_tp], 5, Vec4f(1,0,1,1) );
                }
                break;
#ifdef __TO_PORT_TO_3D
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            VIZ_SEGMENT3( vs, m_pMeshS->V_Pos_0( vid ), m_vec_TM_Random_TargetPos[vid], 1, Vec4f(1,0,1,0.5) );
                            VIZ_POINT3( vs, m_vec_TM_Random_TargetPos[vid], 5, Vec4f(1,0,1,1) );
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
#endif
            default: break;
            }
        }
        break;
    case Params::eTM_Plane:
        {
            Real half_height( 2*m_AABB0.GetHalfSizes()[1]*(Real(1)-m_Params.m_TestFraction) );
            Real half_size( 10 );
            Vec3 normal(0,1,0);
            Vec3 p0(0,half_height,0);
            Vec3 u( mal::Sq(normal[0]+1), mal::Sq(normal[1]+2), mal::Sq(normal[2]+3) ); //obtain linearly indep vector
            Vec3 v( mal::Normalized( mal::Cross(u, normal) ) );
            u = mal::Cross(normal, v);
            Vec3 p1 = p0 - half_size * u - half_size * v;
            Vec3 p2 = p0 + half_size * u - half_size * v;
            Vec3 p3 = p0 + half_size * u + half_size * v;
            Vec3 p4 = p0 - half_size * u + half_size * v;
            VIZ_QUAD3( vs, p1, p2, p3, p4, 1.0f, Vec4f(1,0,1,0.5), util::eVizStyle_Solid );
            p1 -= 2*half_height*normal;
            p2 -= 2*half_height*normal;
            p3 -= 2*half_height*normal;
            p4 -= 2*half_height*normal;
            VIZ_QUAD3( vs, p1, p2, p3, p4, 1.0f, Vec4f(1,0,1,0.5), util::eVizStyle_Solid );
            /*
            VIZ_PLANE3( vs, Vec3(0,1,0), half_height, Vec4f(1,0,1,0.5), util::eVizStyle_Wire );
            VIZ_PLANE3( vs, Vec3(0,1,0), -half_height, Vec4f(1,0,1,0.5), util::eVizStyle_Wire );
            */
        }
        break;
    case Params::eTM_Sphere:
        {
            Real radius( 1.1f*m_AABB0.GetHalfSizes().Norm()*(Real(1)-m_Params.m_TestFraction) );
            VIZ_SPHERE3_NO_ROT( vs, Vec3(0,0,0), radius, Vec4f(1,0,1,0.5), util::eVizStyle_Solid );
        }
        break;
    default: break;
    }
#endif

#ifdef __S3D_ENABLE_CONTACT_D2K
    if( m_DDF.Test( eDraw_Contacts ) )
        for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
            geo::VizContactData( it_ccd2k->m_TmpCD, vs, geo::eContactDataDDF_Points | geo::eContactDataDDF_Normals );
#endif
}

void LeafDSH_Solid3D_FEM::QueryStats( Flags32 flags, ReturnStream &rets ) const
{
#ifdef __S2_DS_ENABLE_STATS
    rets.BeginComplex("stats", eType_Property_Object );
    {
        rets.BeginComplex( "<Energy>", eType_Property_Group );
        {
            Real T( GetKineticEnergy() );
            Real V( GetPotentialEnergy() );
            Real E( GetElasticEnergy() );
            rets.Write( "kinetic", T );
            rets.Write( "potential", V );
            rets.Write( "elastic", E );
            rets.Write( "total", Pair_f32_f32( m_TotalTime, T+V+E ) );

            //TEMP: PS energy is QUITE similar to FEM energy and oscilates comparably...
            Real T1( m_Model.ComputeKineticEnergy() );
            Real V1( m_Model.ComputePotentialEnergy( m_Params.m_Gravity ) );
            rets.Write( "k1", T1 );
            rets.Write( "p1", V1 );
            rets.Write( "t1", Pair_f32_f32( m_TotalTime, T1+V1+E ) );
        }
        rets.EndComplex();

        rets.BeginComplex( "<Misc>", eType_Property_Group );
        {
            rets.Write( "#nodes", m_NumNodes );
            rets.Write( "#elements", m_NumElements );
            rets.Write( "total_mass", m_Model.GetMass() );
            rets.Write( "total_volume", m_TES.TotalVolume() );
            rets.Write( "airdrag/s", m_AirDragPerSecond );
            rets.Write( "total_time", m_TotalTime );
        }
        rets.EndComplex();

        rets.BeginComplex( "<SolverLS>", eType_Property_Group );
        {
            rets.Write( "#iter", Pair_f32_f32( m_TotalTime, (float32)m_Stats.m_SolverLS_NumIter ) );
            rets.Write( "rel_prec", Pair_f32_f32( m_TotalTime, m_Stats.m_SolverLS_RelPrec ) );
            rets.Write( "abs_prec", Pair_f32_f32( m_TotalTime, m_Stats.m_SolverLS_AbsPrec ) );
        }
        rets.EndComplex();

        rets.BeginComplex( "<SolverNR>", eType_Property_Group );
        {
            rets.Write( "#iter", Pair_f32_f32( m_TotalTime, (float32)m_Stats.m_SolverNR_NumIter ) );
            rets.Write( "rel_prec", Pair_f32_f32( m_TotalTime, m_Stats.m_SolverNR_RelPrec ) );
            rets.Write( "abs_prec", Pair_f32_f32( m_TotalTime, m_Stats.m_SolverNR_AbsPrec ) );
        }
        rets.EndComplex();

        rets.BeginComplex( "<VarStep>", eType_Property_Group );
        {
            rets.Write( "CPU time (ms)", Pair_f32_f32( m_TotalTime, 1000.0*m_Stats.m_VariableStep_Duration ) );
            rets.Write( "LS iter", m_Stats.m_VariableStep_LS_NumIter );
            rets.Write( "NR iter", m_Stats.m_VariableStep_NR_NumIter );
            rets.Write( "NR abs prec", m_Stats.m_VariableStep_NR_AbsPrec );
            rets.Write( "Energy", m_Stats.m_VariableStep_Energy );
        }
        rets.EndComplex();

        rets.BeginComplex( "<Degeneration>", eType_Property_Group );
        {
            ms::fem::TetrahedronElementSet3::Stats stats;
            m_TES.ComputeStats( m_pTESC, stats );
            rets.Write( "#degenerate", Pair_f32_f32( m_TotalTime, (float32)stats.m_Num_Degenerate ) );
            rets.Write( "Sum det(F) < 0", Pair_f32_f32( m_TotalTime, stats.m_Sum_Degenerate_DetF ) );
        }
        rets.EndComplex();
    }
    rets.EndComplex();
#endif //__S2_DS_ENABLE_STATS
}

void LeafDSH_Solid3D_FEM::QueryParams( Flags32 flags, ReturnStream &rets ) const
{
#ifdef __S2_DS_ENABLE_PARAMS
    rets.BeginComplex("params", eType_Property_Object );
       Params::s_LeafDSH_Solid3D_FEM_Params_ArchetypeLibrary.ExportInstance( "Archetype_s2_ds_LDSH_FEMS3D_Params", &m_Params, rets );
    rets.EndComplex();
#endif
}

void LeafDSH_Solid3D_FEM::SyncParams( const ParamIt &pit, ReturnStream &rets )
{
#ifdef __S2_DS_ENABLE_PARAMS
    util::ItemStream::ItemItRW out_pit = rets.WriteItem( "params", pit );
    Params::s_LeafDSH_Solid3D_FEM_Params_ArchetypeLibrary.SyncInstance( "Archetype_s2_ds_LDSH_FEMS3D_Params", &m_Params, out_pit.GetSubItem() );
#endif
}

#ifdef __S3D_ENABLE_CONTACT_D2K
bool LeafDSH_Solid3D_FEM::CheckError_Contact() const
{
    //\todo compute some sensible upper limit from object sizes
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
        for( unsigned int it_cpc=0; it_cpc<it_ccd2k->GetNumCPC(); it_cpc++ )
            if( it_ccd2k->GetCPC(it_cpc).m_Depth > m_Params.m_ContactSolver_MaxDepth )
                return false;
    return true;
}
void LeafDSH_Solid3D_FEM::NotifyPairBP( const geo::bp::Proxy *p_other, const geo::bp::Pair *p_pair )
{
    //DS_LOG_WARNING("LeafDSH_FEM_Solid3D_Simplified::NotifyPairBP");
    IEntity *pOtherEntity( p_other->GetObj<IEntity*>() );
    DS_ASSERT( eEntity_Geom == pOtherEntity->GetEntityType() );
    IGeom *pGeom( static_cast<IGeom*>( pOtherEntity ) );

#ifdef __ENABLE_PERSISTENT_CONTACTS
    switch( p_pair->m_State )
    {
    case geo::bp::Pair::eNew:
        {
            //DS_LOG_WARNING("New BP");
            ContactConstraintD2K *pCCD2K( m_poolCCD2K.New() );
            pCCD2K->Create();
            p_pair->SetUserData( pCCD2K );
            if( geo::mp::TestContact( static_cast<const geo::IObject3*>(base_type::GetGO()),
                                      static_cast<const geo::IObject3*>(pGeom->GetGO()),
                                      pCCD2K->m_TmpCD, &pCCD2K->m_CC ) ) //\todo Avoid casts specializing GetGO() with dimension
                pCCD2K->Reset( pCCD2K->m_TmpCD, 0.001f );
        }
        break;
    case geo::bp::Pair::ePersistent:
        {
            //DS_LOG_WARNING("Persistent BP");
            ContactConstraintD2K *pCCD2K( p_pair->GetUserData<ContactConstraintD2K*>() );
            if( geo::mp::TestContact( static_cast<const geo::IObject3*>(base_type::GetGO()),
                                      static_cast<const geo::IObject3*>(pGeom->GetGO()),
                                      pCCD2K->m_TmpCD, &pCCD2K->m_CC ) ) //\todo Avoid casts specializing GetGO() with dimension
                pCCD2K->Persist( pCCD2K->m_TmpCD, 0.001f, 0.001f );
            else
                pCCD2K->Clear();
        }
        break;
    case geo::bp::Pair::eVanished:
        {
            //DS_LOG_WARNING("Vanished BP");
            ContactConstraintD2K *pCCD2K( p_pair->GetUserData<ContactConstraintD2K*>() );
            pCCD2K->Destroy(); //This invalidates the cache too
            m_poolCCD2K.Delete( pCCD2K );
            p_pair->SetUserData(0);
        }
        break;
    default: break;
    }
#else
    m_poolCCD2K.Clear();
    switch( p_pair->m_State )
    {
    case geo::bp::Pair::eNew:
    case geo::bp::Pair::ePersistent:
        {
            geo::np::ContactCache3 *p_cc(0);
            // \todo an IGeom3 should return a geo::IObject3, but an IGeom should return an geo::IObject, consider 2 different GetGO() methods...
            //if( geo::mp::TestContact( base_type::GetGO(), pGeom->GetGO(), m_TmpCD, p_cc ) )
            if( geo::mp::TestContact( static_cast<const geo::IObject3*>(base_type::GetGO()),
                                      static_cast<const geo::IObject3*>(pGeom->GetGO()),
                                      m_TmpCD, p_cc ) ) //TEMP: casts are ugly
            {
                // detect contact
                //\todo Create contact constraint, store in DSH, solve...
                ContactConstraintD2K *pCCD2K = m_poolCCD2K.New();
                pCCD2K->Reset( m_TmpCD, 0.001f );
            }
        }
        break;
    case geo::bp::Pair::eVanished: break;
    default: break;
    }
#endif //__ENABLE_PERSISTENT_CONTACTS
    if( m_Params.m_ContactSolver_BreakOnError && !CheckError_Contact() )
    {
        DS_LOG_ERROR("Breaking simulation due to Contact Error!");
        m_Params.m_bRun = false;
    }
}
#endif //__S3D_ENABLE_CONTACT_D2K

void LeafDSH_Solid3D_FEM::UpdateKNC( Real dt )
{
#ifdef __S3D_ENABLE_KNC
    switch( m_Params.m_KinematicMode )
    {
    case Params::eKM_Static:
        for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
            m_KNCS.SetPosKNC( it_knc.GetPtr(), m_pMeshS->GetVecDefaultSDOF()[ it_knc->m_NID ] );
        break;
    case Params::eKM_Periodic:
        for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
            m_KNCS.SetPosKNC( it_knc.GetPtr(), mal::Sin( m_TotalTime )*Vec3(0,1,0) + m_pMeshS->GetVecDefaultSDOF()[ it_knc->m_NID ] );
        break;
    case Params::eKM_Disabled:
    case Params::eKM_Mouse:
    default:
        break;
    }
#endif
}

void LeafDSH_Solid3D_FEM::FixedStep( Real dt )
{
    S2_DS_TRACE_INIT(); //\todo should be global someday, now it's fine here

    S2_DS_TRACE_BEGIN_SLICE( m_TotalTime );

    UpdateKNC( dt );

#ifdef __S3D_ENABLE_CONTACT_D2K
    /*\note This actually CHANGES node positions in m_Model to fix
      penetration. This change is available to the integration step,
      and therefore accounted for in BOTH external and elastic forces
      in each Begin/EndEvaluation phase. This is slightly DIFFERENT
      from the reference B&W "position alteration" idea, as they only
      anticipate node displacement effect in the elastic forces on the
      right hand of Adv=b).
    */
    if( m_Params.m_ContactSolver_PCST == Params::ePCST_Alteration )
        ApplyContacts_PositionAlteration_Smoothed( m_Model.GetVecPos_RW() ); //\todo Smoothed SHOULD work betterer...
        // ApplyContacts_PositionAlteration( m_Model.GetVecPos_RW() );

#endif

    /* Here we have previous configuration C_i-1 (m_vecPos0) and the
       current C_i (m_Model::GetPos), which MAY not fulfill current
       KNC set after the last FixedStep() call. These will NOT be
       enforced here, but during integration, and only fulfilled at
       the end of the timestep.
    */
    m_TES.BeginEvaluation( m_pTESC, m_vecPos0, m_Model.GetVecPos(), m_Model.GetVecVel(), m_TotalTime );
    {
        // Afterwards, we compute F,T due to the trajectory C_i-1 => C_i^+
        /*TEMP: Old code BEFORE save pos, substituted by BeginEvaluation()
          Compute_F( m_vecF, m_vecDetF );
          Compute_T( m_vecTe2r );
        */

        // Save C_i into C_i-1 for the next frame BEFORE time-integration
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            m_vecPos0[it_node] = m_Model.GetPos(it_node);

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
        UpdateTest(dt);
#endif

        // Integrate C_i^+ up to C_i+1
        switch( m_Params.m_IntegratorType )
        {
        case Params::eIT_SymplecticEuler:
            Step_SymplecticEuler(dt);
            break;
        case Params::eIT_QuasiImplicitEuler_v:
            Step_QuasiImplicitEuler_v(dt);
            break;
        case Params::eIT_QuasiImplicitEuler_dx:
            //Step_QuasiImplicitEuler_dx(dt);
            break;
        case Params::eIT_FullyImplicitEuler_dx:
            Step_FullyImplicitEuler_dx(dt);
            break;
        case Params::eIT_FullyImplicitEuler_dx_QuasiStatic:
            Step_FullyImplicitEuler_dx(0);
            break;
        default:
            DS_ASSERT(false);
        }

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
        if( m_Params.m_bEnableHackedGroundPlane ) ApplyHackedGroundPlane();
#endif

#ifdef __S3D_ENABLE_CONTACT_D2K
        /*IMPORTANT: Disabled for Thesis, we use Hack for the version that worked best, see elsewhere
        if( m_Params.m_ContactSolver_PCST == Params::ePCST_Hack
            || m_Params.m_ContactSolver_VCST == Params::eVCST_Hack ) ApplyContacts_Hack();
        */
#endif

        // Compute elastic energy
        m_ElasticEnergy = m_TES.V( m_pTESC );
    }
    m_TES.EndEvaluation( m_pTESC );

    // Air drag
    if( m_AirDragPerSecond > 0 )
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            m_Model.SetVel( it_node, (Real(1)-dt*m_AirDragPerSecond)*m_Model.GetVel(it_node) );

    // Final energy
    RecomputeEnergy();

    S2_DS_TRACE_END_SLICE(true);
}

/*! Integrate dynamics from (x_0,v_0) to (x_1,v_1)
 */
void LeafDSH_Solid3D_FEM::Step_SymplecticEuler( Real dt )
{
    // Reset forces
    memset( m_vecForces, 0, sizeof(Vec3) * m_NumNodes );

    // Accumulate gravity forces
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_vecForces[it_node] += m_Model.GetMass(it_node)*m_Params.m_Gravity;

#ifdef __S3D_ENABLE_CONTACT_D2K
    // Accumulate penalty forces
    if( m_Params.m_ContactSolver_PCST == Params::ePCST_Penalty
        || m_Params.m_ContactSolver_VCST == Params::eVCST_Penalty )
        ApplyContacts_Penalty( m_vecForces );
#endif

    // Accumulate elastic and damping forces
    m_TES.f_e( m_pTESC, m_vecForces );
    m_TES.f_d( m_pTESC, m_vecForces );

    // Apply total forces
    m_Model.ApplyForce( m_vecForces );

#ifdef __S3D_ENABLE_KNC
    for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
    {
        unsigned int nid( it_knc->m_NID );
        m_Model.ApplyForce( nid, -m_Model.GetAccForce(nid) );
        m_Model.SetPos( nid, it_knc->m_Pos );
        m_Model.SetVel( nid, it_knc->m_Vel );
    }
#endif

    // Update PS (internally symplectic euler)
    m_Model.Update(dt);
}

//---- Linear System Internal methods

/* Set all kinematic node entries to 0, leave dynamic node entries
   unchanged. This should have the same effect as explicitly
   assembling the KNC into A*x = b with:

   For i Kinematic, j Dynamic
     A'_ii = Id
     A'_ij = 0
     A'_ji = 0
     b_j -= A_ji * b_i => A_ji = 0

   By correcting the products A*v instead.
*/
void LeafDSH_Solid3D_FEM::LS_ZeroKNC( Vec3* v ) const
{
#ifdef __S3D_ENABLE_KNC
    if( m_Params.m_KinematicMode != Params::eKM_Disabled
        && m_Params.m_ConstraintSolver_KNCST == Params::eKNCST_Reaction )
        for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
            v[it_knc->m_NID] = Vec3::Zero();
#endif
}

// Computes q = A*d
void LeafDSH_Solid3D_FEM::LS_A_times_d( Real* q, Vec3* vq,
                                        Real factor_M, Real factor_K,
                                        const Real* d, const Vec3* vd ) const
{
    /* q <== A*d
         <== (1 + h rayleigh_m) * (m_i * d_i) - (h^2 + h rayleigh_k) * df_x(d)
         <== factor_M + m_i * d_i - factor_K * df_x(d)
      q0 <== 0
      q1 <== q0 + df_x(d);
      q2 <== - factor_K * q1;
      q  <== q2 + factor_M * (m_i * d_i)
    */
    const unsigned int cN( 3 * m_NumNodes );
    ns::real_array::zero( q, cN ); //q0
    m_TES.df_x( m_pTESC, vd, vq ); //q1 == vq += df_x(vd)
    ns::real_array::muleq( q, -factor_K, cN ); //q2
    for( unsigned int i = 0; i < m_NumNodes; i++ )
    {
        Real factor( m_Model.GetMass(i) * factor_M );
        q[3*i]   += factor * d[3*i];
        q[3*i+1] += factor * d[3*i+1];
        q[3*i+2] += factor * d[3*i+2];
    }
}

void LeafDSH_Solid3D_FEM::LS_Solve( Real factor_M, Real factor_K, Vec3* vy, Vec3* vb, Real& ls_rel_prec, Real& ls_abs_prec, int& ls_iter ) const
{
    // Setup Av params
    m_ILSS_FactorM = factor_M;
    m_ILSS_FactorK = factor_K;

#ifdef __ENABLE_ILSS_CHECK_INDEFINITE
    Real lambda_ratio = ns::ComputeMinMaxEigenvalueRatio_Symmetric( boost::bind( &LeafDSH_Solid3D_FEM::ILSS_D_Av, this, _1, _2 ), 3*m_NumNodes,
                                                                    m_LS_real_array_tmp0, m_LS_real_array_tmp1, m_LS_real_array_tmp2,
                                                                    25, 0.01 );
    if( lambda_ratio < 0 ) DS_LOG_ERROR("Indefinite A with lambda_ratio = %f", lambda_ratio );
#endif

#ifdef __S2_DS_ENABLE_TRACE
    //\todo This is some OLD code that traced definiteness and symmetry, temporally kept as a reference
    {
        /* Sparse eigenvalue test
           Real lambda_max(0);
           //Real *eigen_vector( reinterpret_cast<Real*>( m_CG_vec_r ) ); //r unused by now, used
           // Find largest EV of A using power method
           Real error = LS_ComputeLargestEigenValue_A( factor_M, factor_K, 100, 0.01, 0, lambda_max );//, eigen_vector ); //d,q unused here, used inside
           //S2_DS_TRACE_LOCAL( "A.lambda_max", lambda_max );
           //S2_DS_TRACE_LOCAL( "A.lambda_max.err", error );
           //S2_DS_TRACE_LOCAL_ARRAY( "A.v_max", eigen_vector, cN );
           // Find smallest EV of A using shifted power method (http://math.stackexchange.com/questions/271864/power-iteration-smallest-eigenvalue)
           Real lambda_min(0);
           error = LS_ComputeLargestEigenValue_A( factor_M, factor_K, 100, 0.01, lambda_max, lambda_min );
           lambda_min += lambda_max;
           //S2_DS_TRACE_LOCAL( "A.lambda_min", lambda_min );
           S2_DS_TRACE_LOCAL( "A.lambda_min/max", lambda_min/lambda_max ); //If this is negative, A is indefinite!
        */

        /* Alternative EV using assembled matrix, to check assembly-less version and export A to Octave
        ns::RealD eigen_value_asm(0);
        ns::VectorD eigen_vector_asm( cN );
        ns::MatrixD A_asm( cN, cN );
        LS_Assemble_A( A_asm, factor_M, factor_K );
        {
            // Ouput A in Octave format for easy copy-paste
            std::cout << "A = " << A_asm << std::endl;
            const Real cAsymmetryThreshold(1);
            Real asymmetry( ns::ComputeLargestAssymetry(A_asm) );
            if( asymmetry <= cAsymmetryThreshold ) std::cout << "|A - At| = " << asymmetry << " <= " << cAsymmetryThreshold << std::endl;
            else std::cout << "|A - At| = " << asymmetry << " > " << cAsymmetryThreshold << " ASSYMETRIC!" << std::endl;
        }
        ns::ComputeLargestEigenvalue_Symmetric( A_asm, eigen_value_asm, eigen_vector_asm, 100, 0.01 );
        S2_DS_TRACE_LOCAL( "A_asm.lambda_max", Real(eigen_value_asm) );
        */
    }
#endif

    // Setup solver
    ns::IIterativeLinearSystemSolver *pILSS(0);
    switch( m_Params.m_SolverLS_Type )
    {
    case Params::eLSST_CG: pILSS = &m_ILSS_CG; break;
    case Params::eLSST_CR: pILSS = &m_ILSS_CR; break;
    case Params::eLSST_SQMR: pILSS = &m_ILSS_SQMR; break;
    case Params::eLSST_MINRES: pILSS = &m_ILSS_MINRES; break;
    case Params::eLSST_GMRES: pILSS = &m_ILSS_GMRES; break;
    default: break;
    }
    // Call solver
    if( 0 != pILSS )
    {
        pILSS->Solve( boost::bind( &LeafDSH_Solid3D_FEM::ILSS_D_Av, this, _1, _2 ),
                      boost::bind( &LeafDSH_Solid3D_FEM::ILSS_D_ApplyConstraints, this, _1 ),
                      reinterpret_cast<Real*>(vy),
                      reinterpret_cast<const Real*>(vb),
                      ls_rel_prec, ls_abs_prec, ls_iter );
        ls_rel_prec = pILSS->GetRelResidual();
        ls_abs_prec = pILSS->GetAbsResidual();
        ls_iter = pILSS->GetNumIterations();
    }
}

void LeafDSH_Solid3D_FEM::ILSS_D_Av( Real* Av, const Real* v ) const
{
    LS_A_times_d( Av, reinterpret_cast<Vec3*>(Av),
                  m_ILSS_FactorM, m_ILSS_FactorK,
                  v, reinterpret_cast<const Vec3*>(v) );
}
void LeafDSH_Solid3D_FEM::ILSS_D_ApplyConstraints( Real* v ) const
{
    /* Constraint application order IS important if constraints are
       NOT strictly dissipative (eg: depend on relative vel and may
       INCREASE specific node vel in previously constrained
       directions). The highes priority the latter it must be applied.
    */
    // if( m_Params.m_ConstraintSolver_NDCST == Params::eNDCST_Reaction )
    //     m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, reinterpret_cast<Vec3*>(v) ); \todo
#ifdef __S3D_ENABLE_CONTACT_D2K
    if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
        ApplyContacts_Reaction( reinterpret_cast<Vec3*>(v) );
#endif
    LS_ZeroKNC( reinterpret_cast<Vec3*>(v) );
}


/*! Quasi-Implicit Euler with \Delta v:

  Linear QIE that solves for v^{n+1} and assumes constant df_x = K^n

  With:
    x^{n+1} = x^n + h v^{n+1}
    v^{n+1} = v^n + h M^{-1} f( x^{n+1}, v^{n+1} )

  Solve for:
    A y = b

  Where:
    A = (1 + h rayleigh_m) M
        + (h^2 + h rayleigh_k) K;
    y = v^{n+1}
    b = M v^n + h f^n
    f^n = f_ext + f_s^n //WITOUT fd^n
    h = \Delta t

  KNC:
    Instead of building KNC-support into each per-node computation,
    all per-node methods assume all nodes to be free, and we apply KNC
    afterwards by overwritting KNC entries in assembled per-node
    arrays.
*/
void LeafDSH_Solid3D_FEM::Step_QuasiImplicitEuler_v( Real dt )
{
    Real inv_dt( mal::Rcp(dt) );

    //---- Assemble b (\todo Could be done in a SINGLE PASS except for dt*f_s^n)
    // Reset b
    memset( m_LS_vec_b, 0, sizeof(Vec3) * m_NumNodes );
    // Accumulate External forces (Gravity)
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_LS_vec_b[it_node] += m_Model.GetMass(it_node)*m_Params.m_Gravity;
#ifdef __S3D_ENABLE_CONTACT_D2K
    // Accumulate penalty forces (explicit!)
    if( m_Params.m_ContactSolver_PCST == Params::ePCST_Penalty
        || m_Params.m_ContactSolver_VCST == Params::eVCST_Penalty )
        ApplyContacts_Penalty( m_LS_vec_b );
#endif
    // Accumulate Elastic forces WITHOUT damping
    m_TES.f_e( m_pTESC, m_LS_vec_b );
    // Multiply f^n = (f_ext^n + f_s^n) by dt
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_LS_vec_b[it_node] *= dt;
    // Accumulate M v^n
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_LS_vec_b[it_node] += m_Model.GetMass(it_node)*m_Model.GetVel(it_node);

    //---- Solve A y = b
    if( m_Params.m_SolverLS_bEnableWarmstarting ) memcpy( m_LS_vec_y, m_Model.GetVecVel(), sizeof(Vec3) * m_NumNodes );
    else memset( m_LS_vec_y, 0, sizeof(Vec3) * m_NumNodes );

    //TEMP: Random initial guess avoids breakdown, according to: "LANCZOS METHODS FOR THE SOLUTION OF NONSYMMETRIC SYSTEMS OF LINEAR EQUATIONS*"
    if( m_Params.m_SolverLS_Type == Params::eLSST_CGS ) for( unsigned int i=0; i<m_NumNodes; i++ ) m_LS_vec_y[i] = Vec3( mal::Random01<Real>(), mal::Random01<Real>(), mal::Random01<Real>() );

    /*TEMP: Apply contacts to b and y
      \note THIS is equivalent to the:
         x0 = Sy + (I-S)z
      initialization in the paper
      2003_OnTheModifiedConjugateGradientMethodInClothSimulation,
      however, for D2K constraints we're NOT SETTING any z (target
      values) (we DO for KNC, which completely overwrite x0 entries.
      \todo WHICH is the correct value of z for D2K??  \todo DO WE
      NEED to apply it to b?... it seems so, but actually CG or any
      LSS ALREADY applies constraints to b at init. which SHOULD BE
      equivalent but, for some reason, causes different behaviour
     */
#ifdef __S3D_ENABLE_CONTACT_D2K
    // Contacts
    if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
    {
        UpdateContacts_ActiveSet( m_Model.GetVecVel(), m_LS_vec_b );
        ApplyContacts_NormalImpulse( m_LS_vec_b, 1 );
        ApplyContacts_Reaction( m_LS_vec_y );
    }
#endif
    /*\todo Non-degeneration
    if( m_Params.m_ConstraintSolver_NDCST == Params::eNDCST_Reaction )
    {
        m_TES.UpdateConstraints_ActiveSet( m_pTESC, m_Model.GetVecVel(), m_LS_vec_b );
        m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, m_LS_vec_b );
        m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, m_LS_vec_y );
    }
    */

    // Apply KNC to b and y AFTER other constraints (and therefore has highest priority)
#ifdef __S3D_ENABLE_KNC_IN_LS
    if( m_Params.m_KinematicMode != Params::eKM_Disabled
        && m_Params.m_ConstraintSolver_KNCST == Params::eKNCST_Reaction )
    {
        //m_KNCS.UpdateConstraints_ActiveSet( m_Model.GetVecVel(), m_LS_vec_b ); //\todo by now KNC are ALWAYS active, but Update() call required anyway
        //m_KNCS.ApplyConstraints_Reaction(m_LS_vec_b);
        //m_KNCS.ApplyConstraints_Reaction(m_LS_vec_y);
        for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int nid( it_knc->m_NID );
            // If it's a position constraint, Compute vel from position
            // error, otherwise, it's already a velocity constraint
            if( it_knc->m_Type == ms::KinematicNodeConstraintSet3::KNC::eKNCT_Position )
            {
                Vec3 vk( ( it_knc->m_Pos - m_Model.GetPos(nid) ) * inv_dt );
                Real norm_vk( vk.Norm() );
                if( norm_vk > Real(0.001f) )
                    it_knc->m_Vel = mal::Min( norm_vk, 50.0f ) * vk.Normalized();
                else
                    it_knc->m_Vel = Vec3::Zero();
            }
            // Set b_knc to KNC vel
            m_LS_vec_b[nid] = it_knc->m_Vel;
            // Set y_knc to KNC vel
            m_LS_vec_y[nid] = it_knc->m_Vel;
        }
    }
#endif //__S3D_ENABLE_KNC_IN_LS

    Real ls_rel_prec( m_Params.m_SolverLS_RelEpsilon );
    Real ls_abs_prec( m_Params.m_SolverLS_AbsEpsilon );
    int ls_iter( m_Params.m_SolverLS_MaxIter );
    const Real factor_M( Real(1) + dt*m_Params.m_FEM.m_RayleighCoeff[0] );
    const Real factor_K( mal::Sq(dt) + dt*m_Params.m_FEM.m_RayleighCoeff[1] );
    LS_Solve( factor_M, factor_K, m_LS_vec_y, m_LS_vec_b, ls_rel_prec, ls_abs_prec, ls_iter );
    S2_DS_STAT_SET( m_Stats.m_SolverLS_NumIter, ls_iter );
    S2_DS_STAT_SET( m_Stats.m_SolverLS_RelPrec, ls_rel_prec );
    S2_DS_STAT_SET( m_Stats.m_SolverLS_AbsPrec, ls_abs_prec );

#ifdef __S3D_ENABLE_CONTACT_D2K
            // Apply friction to candidate vel BEFORE position integration... this may be wrong...
    if( m_Params.m_ContactSolver_VCST == Params::eVCST_Hack )
    {
        UpdateContacts_ActiveSet( m_LS_vec_y,
                                  0 ); //IMPORTANT: crash if activeset tries to use acceleration, it's forbidden, live with it
        ApplyContacts_NormalImpulse( m_LS_vec_y, dt ); //y units = \Delta x => NormalImpulse = NormalDisplacement added to \Delta x
        ApplyContacts_FrictionImpulse( m_LS_vec_y, dt ); //y units = \Delta x => FrictionImpulse = TangentialDisplacement removed from \Delta x
    }
    else if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
    {
        ApplyContacts_FrictionImpulse( m_LS_vec_y, dt );
    }
#endif

    // Integrate x_1 with v_1
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_Model.SetVel( it_node, m_LS_vec_y[it_node] );
        m_Model.SetPos( it_node, m_Model.GetPos(it_node) + dt*m_LS_vec_y[it_node] );
    }
}

/*! Fully-Implicit Euler with \Delta x:

  Non-Linear (Q)-FIE that solves for x^{n+1}

  With:
    x^{n+1} = x^n + h v^{n+1}
    v^{n+1} = v^n + h M^{-1} f( x^{n+1}, v^{n+1} )
    \Delta x = x^{n+1} - x^n
    \Delta v = (1/h) * \Delta x

  Solve for:
    A y = b

  Where:
    A = (1/h^2 + rayleigh_m/h) M
        + (1 + rayleigh_k/h) K(x_k^{n+1})
    y = \Delta x
    b = (1/h) * M (v^n - v_k^{n+1}) + f_k^{n+1}
    f = f_ext + f_s + f_d + f_p
    h = \Delta t

  KNC:
    Instead of building KNC-support into each per-node computation,
    all per-node methods assume all nodes to be free, and we apply KNC
    afterwards by overwritting KNC entries in assembled per-node
    arrays.
*/
void LeafDSH_Solid3D_FEM::Step_FullyImplicitEuler_dx( Real dt )
{
    DS_ASSERT( m_Params.m_SolverLS_Type != Params::eLSST_CGN ); //TEMP: NOT YET IMPLEMENTED!!

    S2_DS_STAT_SET( m_Stats.m_SolverLS_NumIter, 0 );

    bool bQuasistatic( dt == 0 );
    Real inv_dt( bQuasistatic ? 0 : mal::Rcp(dt) );

    // Init x_k^{n+1} = x^n, v_k^{n+1} = v^n at nr-iteration k = 0
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_IE_vec_x_k[it_node] = m_Model.GetPos(it_node);
        m_IE_vec_v_k[it_node] = m_Model.GetVel(it_node);
    }

    // Perform NR-iterations while not converged, up to max iterations
    unsigned int nr_k(0);
    Real nr_delta_new(0);
    Real nr_delta_0(-1);
    bool bConverged1(false);
    bool bConverged2(false);
    do
    {
        // Begin NR iteration
        //\note BeginEvaluation() is open from FixedStep()

        //---- Assemble b_k (\todo All m_TES f could be done in a SINGLE PASS)
        // Reset b
        memset( m_LS_vec_b, 0, sizeof(Vec3) * m_NumNodes );
        // Accumulate External forces (Gravity)
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            m_LS_vec_b[it_node] += m_Model.GetMass(it_node)*m_Params.m_Gravity;
#ifdef __S3D_ENABLE_CONTACT_D2K
        // Accumulate penalty forces
        if( m_Params.m_ContactSolver_PCST == Params::ePCST_Penalty
            || m_Params.m_ContactSolver_VCST == Params::eVCST_Penalty )
            ApplyContacts_Penalty( m_LS_vec_b );
#endif
        // Accumulate internal forces (f_e, f_d, f_p) at ( x_k^{n+1}, v_k^{n+1} )
        m_TES.f_e( m_pTESC, m_LS_vec_b );
        //m_TES.f_d( m_pTESC, m_LS_vec_b ); //\todo This should be included, according to FEM.pdf, BUT adds noticeable damping...
        m_TES.f_p( m_pTESC, m_LS_vec_b );
        // Accumulate (1/dt) * M (v^n - v_k^{n+1}), subs v_k^{n+1} with (x_k^{n+1} - x^n)/dt
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            m_LS_vec_b[it_node] += inv_dt * m_Model.GetMass(it_node) * ( m_Model.GetVel(it_node) - inv_dt * (m_IE_vec_x_k[it_node] - m_Model.GetPos(it_node) ) );

        //---- Solve A_k y_k = b_k
        // Warmstarting only the FIRST NR-iter with \Delta x_0 = y = dt*vel_0 may be too optimistic and overshoot, but seems to work well. Further iterations use y = 0
        if( nr_k == 0 && m_Params.m_SolverLS_bEnableWarmstarting ) for( unsigned int i=0; i<m_NumNodes; i++ ) m_LS_vec_y[i] = dt*m_Model.GetVel(i);
        else memset( m_LS_vec_y, 0, sizeof(Vec3) * m_NumNodes );

        //TEMP: Random initial guess avoids breakdown, according to: "LANCZOS METHODS FOR THE SOLUTION OF NONSYMMETRIC SYSTEMS OF LINEAR EQUATIONS*"
        if( m_Params.m_SolverLS_Type == Params::eLSST_CGS ) for( unsigned int i=0; i<m_NumNodes; i++ ) m_LS_vec_y[i] = Vec3( mal::Random01<Real>(), mal::Random01<Real>(), mal::Random01<Real>() );

#ifdef __S3D_ENABLE_CONTACT_D2K
        /*\todo THIS changes the active set of each NR iteration,
         which is basically a QP-like approach (see Faure's
         QPContacts)... consider pre-computing the active set only
         once OUTSIDE the loop, may have better convergence properties
         ("2008 Globally coupled collision handling using volume
         preserving impulses" performs N outer iterations with
         changing active-set and M inner iterations for each one.
        */
        // Contacts
        if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
        {
            UpdateContacts_ActiveSet( m_IE_vec_v_k, //m_Model.GetVecVel(),
                                      m_LS_vec_b );
            ApplyContacts_NormalImpulse( m_LS_vec_b, dt ); //\todo This ONLY fixed b, but did not acc normal impulses: ApplyContacts_Reaction( m_LS_vec_b );
            ApplyContacts_Reaction( m_LS_vec_y );
        }
#endif
        /*\todo Non-degeneration
        if( m_Params.m_ConstraintSolver_NDCST == Params::eNDCST_Reaction )
        {
            m_TES.UpdateConstraints_ActiveSet( m_pTESC, m_Model.GetVecVel(), m_LS_vec_b );
            m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, m_LS_vec_b );
            m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, m_LS_vec_y );
        }
        */

        // Apply KNC to b_k and y_k AFTER other constraints (and therefore has highest priority)

#ifdef __S3D_ENABLE_KNC_IN_LS
        //LS_ApplyPosKNC( m_LS_vec_y, m_LS_vec_b, m_IE_vec_x_k );
        if( m_Params.m_KinematicMode != Params::eKM_Disabled
            && m_Params.m_ConstraintSolver_KNCST == Params::eKNCST_Reaction )
        {
            for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
            {
                unsigned int nid( it_knc->m_NID );
                // If it's a position constraint, Compute vel from position
                // error, otherwise, it's already a velocity constraint
                if( it_knc->m_Type == ms::KinematicNodeConstraintSet3::KNC::eKNCT_Position )
                {
                    Vec3 vk( ( it_knc->m_Pos - m_Model.GetPos(nid) ) * inv_dt );
                    Real norm_vk( vk.Norm() );
                    if( norm_vk > Real(0.001f) )
                        it_knc->m_Vel = mal::Min( norm_vk, 50.0f ) * vk.Normalized();
                    else
                        it_knc->m_Vel = Vec3::Zero();
                }
                // Set b_knc to KNC \Delta x
                m_LS_vec_b[nid] = it_knc->m_Pos - m_IE_vec_x_k[nid]; //it_knc->m_Vel; \todo NOT SURE
                // Set y_knc to KNC \Delta x
                m_LS_vec_y[nid] = it_knc->m_Pos - m_IE_vec_x_k[nid]; //it_knc->m_Vel; \todo NOT SURE
            }
        }
#endif //__S3D_ENABLE_KNC_IN_LS

        // Check convergence b_k < epsilon \IMPORTANT: AFTER applying KNC, otherwise b is 0 for stable configurations REGARDLESS of KNC
        {
            const unsigned int cN( 3 * m_NumNodes );
            const Real *b( reinterpret_cast<const Real*>( m_LS_vec_b ) );
            nr_delta_new = ns::real_array::norm_2( b, cN );
            if( nr_delta_0 < 0 ) nr_delta_0 = nr_delta_new;
        }
        bConverged1 = nr_delta_new < m_Params.m_SolverNR_AbsEpsilon; //\todo HACKED VERSION that checks only y_k = \Delta x correction instead of RHS b_k+1
        if( !bConverged1 )
        {
            Real ls_rel_prec( m_Params.m_SolverLS_RelEpsilon );
            Real ls_abs_prec( m_Params.m_SolverLS_AbsEpsilon );
            int ls_iter( m_Params.m_SolverLS_MaxIter );
            const Real factor_M( mal::Sq(inv_dt) + inv_dt * m_Params.m_FEM.m_RayleighCoeff[0] );
            const Real factor_K( Real(1) + inv_dt * m_Params.m_FEM.m_RayleighCoeff[1] );
            //CG_Solve( factor_M, factor_K, m_LS_vec_y, m_LS_vec_b, ls_rel_prec, ls_abs_prec, ls_iter );
            LS_Solve( factor_M, factor_K, m_LS_vec_y, m_LS_vec_b, ls_rel_prec, ls_abs_prec, ls_iter );
            S2_DS_STAT_ADD( m_Stats.m_SolverLS_NumIter, ls_iter );
            S2_DS_STAT_SET( m_Stats.m_SolverLS_RelPrec, ls_rel_prec );
            S2_DS_STAT_SET( m_Stats.m_SolverLS_AbsPrec, ls_abs_prec );

#ifdef __S3D_ENABLE_CONTACT_D2K
            // Apply friction to candidate vel BEFORE position integration... this may be wrong...
            if( m_Params.m_ContactSolver_VCST == Params::eVCST_Hack )
            {
                UpdateContacts_ActiveSet( m_LS_vec_y,
                                          0 ); //IMPORTANT: crash if activeset tries to use acceleration, it's forbidden, live with it
                ApplyContacts_NormalImpulse( m_LS_vec_y, 1 ); //y units = \Delta x => NormalImpulse = NormalDisplacement added to \Delta x
                ApplyContacts_FrictionImpulse( m_LS_vec_y, 1 ); //y units = \Delta x => FrictionImpulse = TangentialDisplacement removed from \Delta x
            }
            else if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
            {
                ApplyContacts_FrictionImpulse( m_LS_vec_y, dt );
            }
#endif

            // Update iterates
            {
                const unsigned int cN( 3 * m_NumNodes );
                const Real *y( reinterpret_cast<Real*>( m_LS_vec_y ) );
                Real *x_k( reinterpret_cast<Real*>( m_IE_vec_x_k ) );
                Real *v_k( reinterpret_cast<Real*>( m_IE_vec_v_k ) );
                ns::real_array::addeq( x_k, y, cN );
                //TEMP: From Vega, THIS is the correct way to do it!!
                for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
                {
                    m_IE_vec_v_k[it_node] = inv_dt * (m_IE_vec_x_k[it_node] - m_Model.GetPos( it_node ));
                }
                //This is what SiggCourse suggests, I think... ns::real_array::addeq_scaled( v_k, inv_dt, y, cN );
                nr_delta_new = ns::real_array::norm_2( y, cN );
                if( nr_delta_0 < 0 ) nr_delta_0 = nr_delta_new;
            }

            // NR iteration ends
            nr_k++;

            //bConverged2 = nr_delta_new < m_Params.m_SolverNR_AbsEpsilon; TEMP THIS SHOULD NOT BE COMMENTED... TEST IT, I don't know WHY I commented it in the first place
            if( !bConverged2 )
            {
                // If not finished, End current Eval and Begin a new one for the potential displacement x_0 -> x_k+1
                m_TES.EndEvaluation( m_pTESC );
                m_TES.BeginEvaluation( m_pTESC, m_Model.GetVecPos(), m_IE_vec_x_k, m_IE_vec_v_k, m_TotalTime ); //\note This does NOT update plasticity, as m_TotalTime has not advanced from last invokation
                //\note If no more iter, the currently open BeginEvaluation() will be closed at the end of FixedStep()
            }
        }
    } while( !bConverged1
             && !bConverged2
             && nr_k < m_Params.m_SolverNR_MaxIter );

    // results
    unsigned int nr_iter = nr_k;
    Real nr_rel_prec = (nr_delta_0>0) ? nr_delta_new / nr_delta_0 : 0;
    Real nr_abs_prec = nr_delta_new;
    // Stats
    S2_DS_STAT_SET( m_Stats.m_SolverNR_NumIter, nr_iter );
    S2_DS_STAT_SET( m_Stats.m_SolverNR_RelPrec, nr_rel_prec );
    S2_DS_STAT_SET( m_Stats.m_SolverNR_AbsPrec, nr_abs_prec );

    // Set final iterates x_k, v_k to Model state
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_Model.SetPos( it_node, m_IE_vec_x_k[it_node] );
        m_Model.SetVel( it_node, m_IE_vec_v_k[it_node] );
    }
}

void LeafDSH_Solid3D_FEM::RecomputeEnergy()
{
    m_PotentialEnergy = Real(0);
    m_KineticEnergy = Real(0);
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        // Gather element data
        unsigned int vec_nid[4];
        m_TES.GetNID( it_e, vec_nid[0], vec_nid[1], vec_nid[2], vec_nid[3] );
        Vec3 x0( m_Model.GetPos( vec_nid[0] ) );
        Vec3 x1( m_Model.GetPos( vec_nid[1] ) );
        Vec3 x2( m_Model.GetPos( vec_nid[2] ) );
        Vec3 x3( m_Model.GetPos( vec_nid[3] ) );
        Real element_mass( m_DensityPerVolume * m_TES.Volume( it_e ) );

        // Consistent kinetic energy? see FEM.pdf
        Vec3 v0( m_Model.GetVel( vec_nid[0] ) );
        Vec3 v1( m_Model.GetVel( vec_nid[1] ) );
        Vec3 v2( m_Model.GetVel( vec_nid[2] ) );
        Vec3 v3( m_Model.GetVel( vec_nid[3] ) );
#ifdef __TO_PORT_TO_3D
        m_KineticEnergy += ( element_mass / Real(12) )
                           * ( mal::NormSq(v1) + mal::NormSq(v2) + mal::NormSq(v3)
                               + v1.x()*v2.x() + v2.x()*v3.x() + v1.x()*v3.x()
                               + v1.y()*v2.y() + v2.y()*v3.y() + v1.y()*v3.y() );
#endif

        // Consistent potential energy?, see FEM.pdf
        m_PotentialEnergy += ( element_mass / Real(4) ) *  mal::Dot( x0+x1+x2+x3, -m_Params.m_Gravity );
    }
}

void LeafDSH_Solid3D_FEM::RebuildAirDrag()
{
    m_AirDragPerSecond = -mal::Log(10e-6) * m_Params.m_AirDragCoeff;
}

void LeafDSH_Solid3D_FEM::RebuildMass()
{
    // Compute consistent lumped masses
    m_DensityPerVolume = m_Params.m_Density;
    m_TotalMass = m_DensityPerVolume * m_TES.TotalVolume();
#ifdef __TO_PORT_TO_3D
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Real weighted_a(0);
        for( geo::MeshSolidShape2::iterator_around_vertex_ccw it_ccw( m_pMeshS->GetIteratorAroundVertexCCW(it_node) );
             it_ccw.IsValid();
             ++it_ccw )
            weighted_area += mal::Rcp(Real(3)) * m_TES.Area( it_ccw.PID() );
        m_Model.SetMass( it_node, weighted_area * m_DensityPerArea );
    }
#else
    Real mass_per_node = m_TotalMass / m_NumNodes;
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.SetMass( it_node, mass_per_node );
#endif
}

unsigned int LeafDSH_Solid3D_FEM::FindClosestNode( const Vec3 &pos ) const
{
    DS_ASSERT( m_NumNodes > 0 );
    unsigned int closest_nid(0);
    float min_dist_sq( mal::NormSq( pos - m_Model.GetPos(0) ) );
    unsigned int it_node(1);
    while( it_node<m_NumNodes )
    {
        float dist_sq( mal::NormSq( pos - m_Model.GetPos(it_node) ) );
        if( dist_sq < min_dist_sq )
        {
            min_dist_sq = dist_sq;
            closest_nid = it_node;
        }
        ++it_node;
    }
    return closest_nid;
}

Vec3 LeafDSH_Solid3D_FEM::ComputePosCoM() const
{
    Vec3 acc_com(0,0,0);
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        acc_com += m_Model.GetMass( it_node ) * m_Model.GetPos( it_node );
    return acc_com / m_TotalMass;
}

Mat3x3 LeafDSH_Solid3D_FEM::ComputeRotCoM() const
{
    return Mat3x3::Identity(); //HACK
}

//---- Edition
void LeafDSH_Solid3D_FEM::SetPosCoM( const Vec3 &pos )
{
    // Compte CoM
    Vec3 current_pos( ComputePosCoM() );
    // Compute displacement
    Vec3 displacement( pos - current_pos );
    // Apply displacement
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.SetPos( it_node, m_Model.GetPos( it_node ) + displacement );
}

void LeafDSH_Solid3D_FEM::SetRotCoM( const Mat3x3 &rot )
{
    // Compte CoM
    Mat3x3 current_rot( ComputeRotCoM() );
    // Compute displacement
    Mat3x3 diff_rot( mal::Inverse(current_rot)*rot );
    // Apply rotation around CoM to all nodes
    Vec3 pos_com( ComputePosCoM() );
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.SetPos( it_node, pos_com + diff_rot * (m_Model.GetPos( it_node ) - pos_com) );
}

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
void LeafDSH_Solid3D_FEM::InitTest()
{
    // Remove existing test KNC, but ignore preexisting non-test KNC
    for( unsigned int it_tknc=0; it_tknc<m_vecTestKNC.size(); it_tknc++ )
        if( 0 != m_vecTestKNC[it_tknc] )
            m_KNCS.RemoveKNC( m_vecTestKNC[it_tknc] );
    m_vecTestKNC.clear();
    m_vec_TM_Random_TargetPos.clear();

    // Add test KNC according to TM
    switch( m_Params.m_TestMode )
    {
    case Params::eTM_Disabled: break;
    case Params::eTM_Random:
        {
            mal::SetRandomSeed( -666 );
            Vec3 hs( 2*m_AABB0.GetHalfSizes() );
            //\note We ALWAYS allocate all test KNC, even if eTN_Boundary is selected, unused KNC will remain 0
            m_vecTestKNC.resize( m_NumNodes, 0 );
            for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ ) m_vec_TM_Random_TargetPos.push_back( mal::RandomV(-hs,hs) );
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                {
                    for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ )
                    {
                        if( m_KNCS.IsKNC(it_node) ) //Pre-Existing KNC are NOT readded or changed
                            m_vecTestKNC[it_node] = 0;
                        else
                            m_vecTestKNC[it_node] = m_KNCS.AddKNC( it_node, m_pMeshS->V_Pos_0(it_node), Vec3::Zero() );
                    }
                }
                break;
#ifdef __TO_PORT_TO_3D
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            if( m_KNCS.IsKNC(vid) ) //Pre-Existing KNC are NOT readded or changed
                                m_vecTestKNC[vid] = 0;
                            else
                                m_vecTestKNC[vid] = m_KNCS.AddKNC( vid, m_pMeshS->V_Pos_0(vid), Vec3::Zero() );
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
#endif
            default: break;
            }
        }
        break;
    case Params::eTM_Plane:
    case Params::eTM_Sphere:
        {
            // Create TestKNC slots, but the actual KNC are created during UpdateTest() when the constraint is violated
            m_vecTestKNC.resize( m_NumNodes, 0 );
        }
        break;
    default: break;
    };
}
void LeafDSH_Solid3D_FEM::UpdateTest( Real dt )
{
    switch( m_Params.m_TestMode )
    {
    case Params::eTM_Disabled: break;
    case Params::eTM_Random:
        {
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                {
                    for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ )
                        if( 0 != m_vecTestKNC[it_node] )
                        {
                            m_vecTestKNC[it_node]->m_Pos = (Real(1)-m_Params.m_TestFraction)*m_pMeshS->V_Pos_0(it_node) + m_Params.m_TestFraction*m_vec_TM_Random_TargetPos[it_node];
                            //TEMP: Avoid NR seeing Test::Random as large displacement induced vel!
                            m_vecPos0[it_node] = m_vecTestKNC[it_node]->m_Pos;
                            m_Model.SetPos( it_node, m_vecTestKNC[it_node]->m_Pos );
                            m_Model.SetVel( it_node, Vec3::Zero() );
                        }
                }
                break;
#ifdef __TO_PORT_TO_3D
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            if( 0 != m_vecTestKNC[vid] )
                            {
                                m_vecTestKNC[vid]->m_Pos = (Real(1)-m_Params.m_TestFraction)*m_pMeshS->V_Pos_0(vid) + m_Params.m_TestFraction*m_vec_TM_Random_TargetPos[vid];
                                //TEMP: Avoid NR seeing Test::Random as large displacement induced vel!
                                m_vecPos0[vid] = m_vecTestKNC[vid]->m_Pos;
                                m_Model.SetPos( vid, m_vecTestKNC[vid]->m_Pos );
                                m_Model.SetVel( vid, Vec3::Zero() );
                            }
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
#endif
            default: break;
            }
        }
        break;
    case Params::eTM_Plane:
        {
            Vec3 normal( 0, 1, 0 );
            Real half_height( 2*m_AABB0.GetHalfSizes()[1]*(Real(1)-m_Params.m_TestFraction) );
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                {
                    for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ )
                    {
                        unsigned int vid( it_node );
                        Vec3 pos( m_Model.GetPos(vid) );
                        Real dot_p_n( mal::Dot( pos, normal ) );
                        Vec3 projected_pos( (dot_p_n > 0)
                                            ? pos - (dot_p_n-half_height) * normal
                                            : pos + (-dot_p_n-half_height) * normal );
                        if( 0 != m_vecTestKNC[vid] ) //KNC exists
                        {
                            if( mal::Abs( dot_p_n ) >= half_height ) //Update if behind plane
                                m_KNCS.SetPosKNC( m_vecTestKNC[vid], projected_pos );
                            else //Remove if not
                            {
                                m_KNCS.RemoveKNC( m_vecTestKNC[vid] );
                                m_vecTestKNC[vid] = 0;
                            }
                        }
                        else if( mal::Abs( dot_p_n ) >= half_height ) //KNC does not exist, add it if behind plane
                        {
                            m_vecTestKNC[vid] = m_KNCS.AddKNC( vid, projected_pos, Vec3::Zero() );
                        }
                    }
                }
                break;
#ifdef __TO_PORT_TO_3D
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            Vec3 pos( m_Model.GetPos(vid) );
                            Real dot_p_n( mal::Dot( pos, normal ) );
                            Vec3 projected_pos( (dot_p_n > 0)
                                                ? pos - (dot_p_n-half_height) * normal
                                                : pos + (-dot_p_n-half_height) * normal );
                            if( 0 != m_vecTestKNC[vid] ) //KNC exists
                            {
                                if( mal::Abs( dot_p_n ) >= half_height ) //Update if behind plane
                                    m_KNCS.SetPosKNC( m_vecTestKNC[vid], projected_pos );
                                else //Remove if not
                                {
                                    m_KNCS.RemoveKNC( m_vecTestKNC[vid] );
                                    m_vecTestKNC[vid] = 0;
                                }
                            }
                            else if( mal::Abs( dot_p_n ) >= half_height ) //KNC does not exist, add it if behind plane
                            {
                                m_vecTestKNC[vid] = m_KNCS.AddKNC( vid, projected_pos, Vec3::Zero() );
                            }
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
#endif
            default: break;
            }
        }
        break;
    case Params::eTM_Sphere:
        {
            Real radius( 1.1f*m_AABB0.GetHalfSizes().Norm()*(Real(1)-m_Params.m_TestFraction) );
            Real radius_sq( mal::Sq(radius) );
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                {
                    for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ )
                    {
                        unsigned int vid( it_node );
                        Vec3 pos( m_Model.GetPos(vid) );
                        Real dist_sq( pos.NormSq() );
                        Vec3 projected_pos( (dist_sq > 0.000001f)
                                            ? radius * pos / mal::Sqrt(dist_sq)
                                            : pos );
                        if( 0 != m_vecTestKNC[vid] ) //KNC exists
                        {
                            if( dist_sq >= radius_sq ) //Update if behind plane
                                m_KNCS.SetPosKNC( m_vecTestKNC[vid], projected_pos );
                            else //Remove if not
                            {
                                m_KNCS.RemoveKNC( m_vecTestKNC[vid] );
                                m_vecTestKNC[vid] = 0;
                            }
                        }
                        else if( dist_sq >= radius_sq ) //KNC does not exist, add it if behind plane
                        {
                            m_vecTestKNC[vid] = m_KNCS.AddKNC( vid, projected_pos, Vec3::Zero() );
                        }
                    }
                }
                break;
#ifdef __TO_PORT_TO_3D
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            Vec3 pos( m_Model.GetPos(vid) );
                            Real dist_sq( pos.NormSq() );
                            Vec3 projected_pos( (dist_sq > 0.000001f)
                                                ? radius * pos / mal::Sqrt(dist_sq)
                                                : pos );
                            if( 0 != m_vecTestKNC[vid] ) //KNC exists
                            {
                                if( dist_sq >= radius_sq ) //Update if behind plane
                                    m_KNCS.SetPosKNC( m_vecTestKNC[vid], projected_pos );
                                else //Remove if not
                                {
                                    m_KNCS.RemoveKNC( m_vecTestKNC[vid] );
                                    m_vecTestKNC[vid] = 0;
                                }
                            }
                            else if( dist_sq >= radius_sq ) //KNC does not exist, add it if behind plane
                            {
                                m_vecTestKNC[vid] = m_KNCS.AddKNC( vid, projected_pos, Vec3::Zero() );
                            }
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
#endif
            default: break;
            }
        }
        break;
    default: break;
    }
}

void LeafDSH_Solid3D_FEM::ApplyHackedGroundPlane()
{
    const Real cContactSolver_DepthMargin(0.001f);
    const Real cContactSolver_Relaxation_Coeff(0.9);
    const Real cContactSolver_Restitution_Coeff(0);
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Vec3 normal( 0, 1, 0 );
        Vec3 pos( m_Model.GetPos(it_node) );
        Real depth( -pos.y()-1 ); //GroundPlane at Y=-1
        if( depth > cContactSolver_DepthMargin )
        {
            m_Model.SetPos( it_node, m_Model.GetPos(it_node) + (cContactSolver_Relaxation_Coeff * depth) * normal );
            // Could apply impulse instead of vel change...
            Vec3 vel( m_Model.GetVel(it_node) );
            Real vel_n( mal::Dot( vel, normal ) );
            if( vel_n < Real(0) ) m_Model.SetVel( it_node, vel - (Real(1)+cContactSolver_Restitution_Coeff) * vel_n * normal );
            // Constrait force
            Real force_n( mal::Dot( m_Model.GetAccForce(it_node), normal ) );
            if( force_n < Real(0) ) m_Model.ApplyForce( it_node, -force_n*normal );
        }
    }
}
#endif

#ifdef __S3D_ENABLE_CONTACT_D2K
/* Given the barycentric weights a,b of a
   point-on-feature P on the inter-node segment
   (p1,p2), find the displacements delta_p1,
   delta_p2

   As in Jakobsen's "Advanced Character Physics",
   we constrain node displacements to be:
   - Along contact normal N.
   - Proportional to their barycentric weight on P.

   Facts:
     p = a*p1 + b*p2
     p' = a*p1' + b*p2'
     delta_p = p'-p
     p1' = p1 + delta_p1
     p2' = p2 + delta_p2
     => delta_p = a*delta_p1 + b*delta_p2 (1)
   Constraint: All displacements along n
     delta_p = lambda*n
     delta_p1 = lambda1*n
     delta_p2 = lambda2*n
     => lambda = a*lambda1 + b*lambda2 (2)
   Constraint: Displacements proportional to barycentric weight
     => lambda1/lambda2 = a/b (3)
   Solution: Substitute (3) into (2) to find lambda1, then find lambda2 from (3)
     lambda1 = a/(a^2+b^2)*lambda
     lambda2 = b/(a^2+b^2)*lambda
   This generalizes to the triangle (a,b,c) and tetrahedron (a,b,c,d) cases.
   \todo NO, better use effective mass for consistency!

   \todo Actually, we SHOULD REEVALUATE depth, as we're effectively changing it with each displacement, GS-like
*/
void LeafDSH_Solid3D_FEM::ApplyContacts_PositionAlteration( Vec3* vec_pos ) const
{
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_POF1.m_FeatureId.IsVertex() )
                {
                    geo::tetsolid3_feature_index_type vid0 = cpc.m_POF1.m_FeatureId.AsVertex();
                    //vec_pos[vid0] += cpc.m_Depth * cpc.m_Normal;
                    vec_pos[vid0] += m_Params.m_ContactSolver_Relaxation_Coeff * (cpc.m_Depth-m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal;
                    // GEO_LOG_WARNING("V[%d] = (%f,%f), d = (%f,%f)", vid1, vec_pos[vid1].x(), vec_pos[vid1].y(), d.x(), d.y() )
                }
                else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
                {
                    geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                    uint32 vec_vid[4];
                    m_pMeshS->T_VecVID( eid, vec_vid );
                    // Compute change at P_cp \todo Actually, we SHOULD REEVALUATE depth, as we're effectively changing it with each displacement, GS-like
                    Vec3 delta_p( m_Params.m_ContactSolver_Relaxation_Coeff * (cpc.m_Depth-m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal );
#ifdef __ENABLE_CONTACT_EFFECTIVE_MASS
                    // Enforce \Delta p distributing it among p1..p4 using effective mass at p
                    Real effective_mass( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * m_Model.GetInvMass(vec_vid[0])
                                                   + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * m_Model.GetInvMass(vec_vid[1])
                                                   + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * m_Model.GetInvMass(vec_vid[2])
                                                   + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * m_Model.GetInvMass(vec_vid[3]) ) );
                    vec_pos[vec_vid[0]] += ( cpc.m_POF1.m_BarycentricCoords[0] * effective_mass * m_Model.GetInvMass(vec_vid[0]) ) * delta_p;
                    vec_pos[vec_vid[1]] += ( cpc.m_POF1.m_BarycentricCoords[1] * effective_mass * m_Model.GetInvMass(vec_vid[1]) ) * delta_p;
                    vec_pos[vec_vid[2]] += ( cpc.m_POF1.m_BarycentricCoords[2] * effective_mass * m_Model.GetInvMass(vec_vid[2]) ) * delta_p;
                    vec_pos[vec_vid[3]] += ( cpc.m_POF1.m_BarycentricCoords[3] * effective_mass * m_Model.GetInvMass(vec_vid[3]) ) * delta_p;
#else
                    // Apply \Delta p with barycentric weights
                    Real rcp_denom( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] )
                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] )
                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] )
                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) ) );
                    vec_pos[vec_vid[0]] += cpc.m_POF1.m_BarycentricCoords[0] * rcp_denom * delta_p;
                    vec_pos[vec_vid[1]] += cpc.m_POF1.m_BarycentricCoords[1] * rcp_denom * delta_p;
                    vec_pos[vec_vid[2]] += cpc.m_POF1.m_BarycentricCoords[2] * rcp_denom * delta_p;
                    vec_pos[vec_vid[3]] += cpc.m_POF1.m_BarycentricCoords[3] * rcp_denom * delta_p;
#endif
                }
                else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
            }
        }
    }
}
// \todo Jacobi-like smoothed position alteration
void LeafDSH_Solid3D_FEM::ApplyContacts_PositionAlteration_Smoothed( Vec3* vec_pos ) const
{
    // Use tmp arrays to avoid new/delete
    Vec3* vec_dp( m_LS_vec_y );
    Real* vec_num_dp( m_LS_real_array_tmp0 );
    memset( vec_dp, 0, sizeof(Vec3) * m_NumNodes );
    memset( vec_num_dp, 0, sizeof(Real) * m_NumNodes );

    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_POF1.m_FeatureId.IsVertex() )
                {
                    geo::tetsolid3_feature_index_type vid0 = cpc.m_POF1.m_FeatureId.AsVertex();
                    vec_dp[vid0] += m_Params.m_ContactSolver_Relaxation_Coeff * (cpc.m_Depth-m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal;
                    vec_num_dp[vid0]++;
                }
                else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
                {
                    geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                    uint32 vec_vid[4];
                    m_pMeshS->T_VecVID( eid, vec_vid );
                    // Compute change at P_cp
                    Vec3 delta_p( m_Params.m_ContactSolver_Relaxation_Coeff * (cpc.m_Depth-m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal );
#ifdef __ENABLE_CONTACT_EFFECTIVE_MASS
                    // Enforce \Delta p distributing it among p1..p4 using effective mass at p
                    Real effective_mass( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * m_Model.GetInvMass(vec_vid[0])
                                                   + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * m_Model.GetInvMass(vec_vid[1])
                                                   + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * m_Model.GetInvMass(vec_vid[2])
                                                   + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * m_Model.GetInvMass(vec_vid[3]) ) );
                    for( int i=0; i<4; i++ )
                    {
                        vec_dp[vec_vid[i]] += ( cpc.m_POF1.m_BarycentricCoords[i] * effective_mass * m_Model.GetInvMass(vec_vid[i]) ) * delta_p;
                        vec_num_dp[vec_vid[i]]++;
                    }
#else
                    // Apply \Delta p with barycentric weights
                    Real rcp_denom( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] )
                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] )
                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] )
                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) ) );
                    vec_dp[vec_vid[0]] += cpc.m_POF1.m_BarycentricCoords[0] * rcp_denom * delta_p;
                    vec_dp[vec_vid[1]] += cpc.m_POF1.m_BarycentricCoords[1] * rcp_denom * delta_p;
                    vec_dp[vec_vid[2]] += cpc.m_POF1.m_BarycentricCoords[2] * rcp_denom * delta_p;
                    vec_dp[vec_vid[3]] += cpc.m_POF1.m_BarycentricCoords[3] * rcp_denom * delta_p;
                    for( int i=0; i<4; i++ )
                    {
                        vec_dp[vec_vid[i]] += ( cpc.m_POF1.m_BarycentricCoords[0] * rcp_denom ) * delta_p;
                        vec_num_dp[vec_vid[i]]++;
                    }
#endif
                }
                else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
            }
        }
    }
    // Apply smoothed displacements
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
        if( vec_num_dp[it_node] > 0 )
            vec_pos[it_node] += vec_dp[it_node] * mal::Rcp<Real>(vec_num_dp[it_node]);
}

// Accumulate penalty forces per-node on vec_forces \todo THIS IS EXPLICIT and causes vibration/instability
void LeafDSH_Solid3D_FEM::ApplyContacts_Penalty( Vec3* vec_forces ) const
{
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        //if( !ccd2k.IsEmpty() ) DS_LOG_WARNING( "Contact %llx", (machine_uint_type)&ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_POF1.m_FeatureId.IsVertex() )
                {
                    geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
                    Vec3 vel( m_Model.GetVel(vid) );
                    Real pressure_n( m_Params.m_ContactSolver_Penalty_Ks * (cpc.m_Depth - m_Params.m_ContactSolver_DepthMargin)
                                     - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                    Real force_n( cpc.m_Radius * pressure_n );
                    vec_forces[vid] += force_n * cpc.m_Normal;
                }
                /*
                else if( cpc.m_POF1.m_FeatureId.IsSegment() )
                {
                    geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsSegment();
                    uint32 vec_vid[2] = { m_pMeshS->HE_OriginVID( eid ), m_pMeshS->HE_FinalVID( eid ) };
                    // Interpolate V_cp from node velocities
                    Vec2 vel( cpc.m_POF1.m_BarycentricCoords[0] * m_Model.GetVel(vec_vid[0])
                              + cpc.m_POF1.m_BarycentricCoords[1] * m_Model.GetVel(vec_vid[1]) );
                    Real pressure_n( m_Params.m_ContactSolver_Penalty_Ks * (cpc.m_Depth - m_Params.m_ContactSolver_DepthMargin)
                                     - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                    Real force_n( cpc.m_Radius * pressure_n );
                    // Distribute F_cp accross nodes
                    Real force_n1( cpc.m_POF1.m_BarycentricCoords[0] * force_n );
                    Real force_n2( cpc.m_POF1.m_BarycentricCoords[1] * force_n );
                    vec_forces[vec_vid[0]] += force_n1 * cpc.m_Normal;
                    vec_forces[vec_vid[1]] += force_n2 * cpc.m_Normal;
                }
                else if( cpc.m_POF1.m_FeatureId.IsTriangle() )
                {
                    geo::tetsolid3_feature_index_type tid = cpc.m_POF1.m_FeatureId.AsTriangle();
                    uint32 vec_vid[3];
                    m_pMeshS->P_VecVID( tid, vec_vid, 3 );
                    // Interpolate V_cp from node velocities
                    Vec2 vel( cpc.m_POF1.m_BarycentricCoords[0] * m_Model.GetVel(vec_vid[0])
                              + cpc.m_POF1.m_BarycentricCoords[1] * m_Model.GetVel(vec_vid[1])
                              + cpc.m_POF1.m_BarycentricCoords[2] * m_Model.GetVel(vec_vid[2]) );
                    Real pressure_n( m_Params.m_ContactSolver_Penalty_Ks * (cpc.m_Depth - m_Params.m_ContactSolver_DepthMargin)
                                     - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                    Real force_n( cpc.m_Radius * pressure_n );
                    // Distribute F_cp accross nodes
                    Real force_n0( cpc.m_POF1.m_BarycentricCoords[0] * force_n );
                    Real force_n1( cpc.m_POF1.m_BarycentricCoords[1] * force_n );
                    Real force_n2( cpc.m_POF1.m_BarycentricCoords[2] * force_n );
                    vec_forces[vec_vid[0]] += force_n0 * cpc.m_Normal;
                    vec_forces[vec_vid[1]] += force_n1 * cpc.m_Normal;
                    vec_forces[vec_vid[2]] += force_n2 * cpc.m_Normal;
                }
                */
                else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
                {
                    geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                    uint32 vec_vid[4];
                    m_pMeshS->T_VecVID( eid, vec_vid );
                    // Interpolate V_cp from node velocities
                    Vec3 vel( cpc.m_POF1.m_BarycentricCoords[0] * m_Model.GetVel(vec_vid[0])
                              + cpc.m_POF1.m_BarycentricCoords[1] * m_Model.GetVel(vec_vid[1])
                              + cpc.m_POF1.m_BarycentricCoords[2] * m_Model.GetVel(vec_vid[2])
                              + cpc.m_POF1.m_BarycentricCoords[3] * m_Model.GetVel(vec_vid[3]) );
                    Real pressure_n( m_Params.m_ContactSolver_Penalty_Ks * (cpc.m_Depth - m_Params.m_ContactSolver_DepthMargin)
                                     - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                    Real force_n( cpc.m_Radius * pressure_n );
                    // Distribute F_cp accross nodes
                    Real force_n0( cpc.m_POF1.m_BarycentricCoords[0] * force_n );
                    Real force_n1( cpc.m_POF1.m_BarycentricCoords[1] * force_n );
                    Real force_n2( cpc.m_POF1.m_BarycentricCoords[2] * force_n );
                    Real force_n3( cpc.m_POF1.m_BarycentricCoords[3] * force_n );
                    vec_forces[vec_vid[0]] += force_n0 * cpc.m_Normal;
                    vec_forces[vec_vid[1]] += force_n1 * cpc.m_Normal;
                    vec_forces[vec_vid[2]] += force_n2 * cpc.m_Normal;
                    vec_forces[vec_vid[3]] += force_n3 * cpc.m_Normal;
                    // GEO_LOG("Contact E[%u]: p_b = (%f,%f,%f,%f), n_0 = (%f,%f,%f), r = %f, d = %f",
                    //         eid,
                    //         cpc.m_POF1.m_BarycentricCoords[0], cpc.m_POF1.m_BarycentricCoords[1], cpc.m_POF1.m_BarycentricCoords[2], cpc.m_POF1.m_BarycentricCoords[3],
                    //         cpc.m_Normal[0], cpc.m_Normal[1], cpc.m_Normal[2],
                    //         cpc.m_Radius, cpc.m_Depth );
                }
                else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
            }
        }
    }
}

void LeafDSH_Solid3D_FEM::ApplyContacts_Reaction( Vec3* vec_v ) const
{
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K& ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K& cpc( ccd2k.GetCPC(it_cpc) );
            //if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_IsActive ) //active = collision || contact = v*n < 0 || \Delta v*n < 0
                {
                    if( cpc.m_POF1.m_FeatureId.IsVertex() )
                    {
                        uint32 vid = cpc.m_POF1.m_FeatureId.AsVertex();
                        vec_v[vid] -= mal::Dot(vec_v[vid],cpc.m_Normal) * cpc.m_Normal;
                    }
                    else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
                    {
                        geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                        uint32 vec_vid[4];
                        m_pMeshS->T_VecVID( eid, vec_vid );
                        // Interpolate V_cp from node velocities
                        Vec3 v( cpc.m_POF1.m_BarycentricCoords[0] * vec_v[vec_vid[0]]
                                + cpc.m_POF1.m_BarycentricCoords[1] * vec_v[vec_vid[1]]
                                + cpc.m_POF1.m_BarycentricCoords[2] * vec_v[vec_vid[2]]
                                + cpc.m_POF1.m_BarycentricCoords[3] * vec_v[vec_vid[3]] );
                        // Compute change at V_cp
                        Vec3 delta_v( -mal::Dot(v,cpc.m_Normal) * cpc.m_Normal );
                        // Barycentric distribution \note If seen as an impulse, this is the CORRECT way to distribute it among nodes
                        vec_v[vec_vid[0]] += cpc.m_POF1.m_BarycentricCoords[0] * delta_v;
                        vec_v[vec_vid[1]] += cpc.m_POF1.m_BarycentricCoords[1] * delta_v;
                        vec_v[vec_vid[2]] += cpc.m_POF1.m_BarycentricCoords[2] * delta_v;
                        vec_v[vec_vid[3]] += cpc.m_POF1.m_BarycentricCoords[3] * delta_v;
                    }
                    else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
                }
            }
        }
    }
}

void LeafDSH_Solid3D_FEM::ApplyContacts_Hack()
{
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        //if( !ccd2k.IsEmpty() ) DS_LOG_WARNING( "Contact %llx", (machine_uint_type)&ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_POF1.m_FeatureId.IsVertex() )
                {
                    geo::tetsolid3_feature_index_type vid1 = cpc.m_POF1.m_FeatureId.AsVertex();
                    // pos,vel,acc correction constraints
                    // Position correction \todo AIXO FA GUANYAR ENERGIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    m_Model.SetPos( vid1, m_Model.GetPos(vid1) + (m_Params.m_ContactSolver_Relaxation_Coeff * cpc.m_Depth) * cpc.m_Normal );
                    // Could apply impulse instead of vel change...
                    Vec3 vel( m_Model.GetVel(vid1) );
                    Real vel_n( mal::Dot( vel, cpc.m_Normal ) );
                    if( vel_n < Real(0) ) m_Model.SetVel( vid1, vel - (Real(1)+m_Params.m_ContactSolver_Restitution_Coeff) * vel_n * cpc.m_Normal );
                    // Constrait force
                    Real force_n( mal::Dot( m_Model.GetAccForce(vid1), cpc.m_Normal ) );
                    if( force_n < Real(0) ) m_Model.ApplyForce( vid1, -force_n*cpc.m_Normal );
                }
                else if( cpc.m_POF1.m_FeatureId.IsSegment() )
                {
                    DS_LOG_ERROR("Unsupported feature_id IsSegment()");
                }
                else
                {
                    DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type );
                }
            }
        }
    }
}

/* Using Bridson's formula in 2002_RobustTreatmentOfCollisionsContactAndFrictionForClothAnimation, sec 7.3
  \todo Should use RELATIVE normal vel in D2K and D2D contacts (not in D2S, as Static have vel=0)
  \note dv_n is the effect on velocity of the applied contact normal impulse (dv = j/m), j impulse, m mass

  IMPORTANT: This fails, because for resting contacts dv is always 0
  after the Ax=b solve, due to Reaction constraints.
*/
void LeafDSH_Solid3D_FEM::ApplyContacts_Friction_dv( const Vec3* vec_vel0, Vec3* vec_vel1 ) const
{
#ifdef __DISABLED_WHILE_DEVELOPING
                        \todo DEPRECATED METHOD, NOT USED IN 2D EITHER
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_IsActive )
            {
                if( cpc.m_POF1.m_FeatureId.IsVertex() )
                {
                    geo::tetsolid3_feature_index_type vid1 = cpc.m_POF1.m_FeatureId.AsVertex();
                    Vec3 v1t( vec_vel1[vid1] - mal::Dot(vec_vel1[vid1],cpc.m_Normal)*cpc.m_Normal );
                    Real norm_sq_v1t( mal::NormSq(v1t) );
                    if( norm_sq_v1t > 1e-12 )
                    {
                        Real dv_n( mal::Dot(vec_vel1[vid1]-vec_vel0[vid1],cpc.m_Normal) );
                        if( dv_n > 0 ) //\todo May be unnecessary if we already check IsActive...?
                        {
                            Real norm_v1t( mal::Sqrt(norm_sq_v1t) );
                            Real post_v1t_over_norm_v1t( mal::Max(0.0f,1.0f-m_Params.m_ContactSolver_DynamicFriction_Coeff*(dv_n/norm_v1t)) );
                            // Final vel1 includes tangential friction but keeps normal component unchanged
                            vec_vel1[vid1] = post_v1t_over_norm_v1t * v1t
                                             + mal::Dot(vec_vel1[vid1],cpc.m_Normal)*cpc.m_Normal;
                            // DS_LOG("V[%d] FF = %f",vid1,post_v1t_over_norm_v1t);
                        }
                        // else
                        //     DS_LOG("V[%d] dv_n <= 0",vid1);
                    }
                    // else
                    //     DS_LOG("V[%d] |v1t| == 0",vid1);
                }
                else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
            }
        }
    }
#endif //__DISABLED_WHILE_DEVELOPING
}

void LeafDSH_Solid3D_FEM::ApplyContacts_NormalImpulse( Vec3* vec_v, Real factor_v_to_impulse )
{
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        ContactConstraintD2K &ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            //if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_IsActive ) //active = collision || contact = v*n < 0 || \Delta v*n < 0
                {
                    if( cpc.m_POF1.m_FeatureId.IsVertex() )
                    {
                        geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
                        cpc.m_LambdaN = -mal::Dot(vec_v[vid],cpc.m_Normal);
                        vec_v[vid] += cpc.m_LambdaN * cpc.m_Normal;
                        cpc.m_LambdaN *= factor_v_to_impulse; //Apply factor to convert input magnitude v to impulse
                    }
                    else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
                    {
                        geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                        uint32 vec_vid[4];
                        m_pMeshS->T_VecVID( eid, vec_vid );
                        // Interpolate V_cp from node velocities
                        Vec3 v( cpc.m_POF1.m_BarycentricCoords[0] * vec_v[vec_vid[0]]
                                + cpc.m_POF1.m_BarycentricCoords[1] * vec_v[vec_vid[1]]
                                + cpc.m_POF1.m_BarycentricCoords[2] * vec_v[vec_vid[2]]
                                + cpc.m_POF1.m_BarycentricCoords[3] * vec_v[vec_vid[3]] );
                        // Compute change at V_cp
                        cpc.m_LambdaN = -mal::Dot(v,cpc.m_Normal);
                        Vec3 delta_v( cpc.m_LambdaN * cpc.m_Normal );
#ifdef __ENABLE_CONTACT_EFFECTIVE_MASS
                        // Enforce \Delta v distributing it among v1..v4 using effective mass at v
                        Real effective_mass( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * m_Model.GetInvMass(vec_vid[0])
                                                       + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * m_Model.GetInvMass(vec_vid[1])
                                                       + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * m_Model.GetInvMass(vec_vid[2])
                                                       + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * m_Model.GetInvMass(vec_vid[3]) ) );
                        vec_v[vec_vid[0]] += ( cpc.m_POF1.m_BarycentricCoords[0] * effective_mass * m_Model.GetInvMass(vec_vid[0]) ) * delta_v;
                        vec_v[vec_vid[1]] += ( cpc.m_POF1.m_BarycentricCoords[1] * effective_mass * m_Model.GetInvMass(vec_vid[1]) ) * delta_v;
                        vec_v[vec_vid[2]] += ( cpc.m_POF1.m_BarycentricCoords[2] * effective_mass * m_Model.GetInvMass(vec_vid[2]) ) * delta_v;
                        vec_v[vec_vid[3]] += ( cpc.m_POF1.m_BarycentricCoords[3] * effective_mass * m_Model.GetInvMass(vec_vid[3]) ) * delta_v;
                        cpc.m_LambdaN *= factor_v_to_impulse; //Apply factor to convert input magnitude v to impulse
#else
                        // Apply \Delta v with barycentric weights
                        vec_v[vec_vid[0]] += cpc.m_POF1.m_BarycentricCoords[0] * delta_v;
                        vec_v[vec_vid[1]] += cpc.m_POF1.m_BarycentricCoords[1] * delta_v;
                        vec_v[vec_vid[2]] += cpc.m_POF1.m_BarycentricCoords[2] * delta_v;
                        vec_v[vec_vid[3]] += cpc.m_POF1.m_BarycentricCoords[3] * delta_v;
                        cpc.m_LambdaN *= factor_v_to_impulse; //Apply factor to convert input magnitude v to impulse
#endif
                    }
                    else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
                }
            }
        }
    }
}

void LeafDSH_Solid3D_FEM::ApplyContacts_FrictionImpulse( Vec3* vec_v, Real factor_impulse_to_v ) const
{
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_POF1.m_FeatureId.IsVertex() )
            {
                geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
                Vec3 vt( vec_v[vid] - mal::Dot(vec_v[vid],cpc.m_Normal)*cpc.m_Normal );
                Real norm_sq_vt( mal::NormSq(vt) );
                if( norm_sq_vt > 1e-8 )
                {
                    // Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN * m_Model.GetInvMass(vid) ); //Apply factor to convert cached impulse to v magnitude
                    Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN ); //Apply factor to convert cached impulse to v magnitude
                    if( delta_vn > 0 )
                    {
                        Real norm_vt( mal::Sqrt(norm_sq_vt) );
                        Real post_vt_over_norm_vt( mal::Max(0.0f,1.0f-m_Params.m_ContactSolver_DynamicFriction_Coeff*(delta_vn/norm_vt)) );
                        // Final v includes tangential friction but keeps normal component unchanged
                        vec_v[vid] = post_vt_over_norm_vt * vt + mal::Dot(vec_v[vid],cpc.m_Normal)*cpc.m_Normal;
                    }
                }
            }
            else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
            {
                geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                uint32 vec_vid[4];
                m_pMeshS->T_VecVID( eid, vec_vid );
                // Interpolate V_cp from node velocities
                Vec3 v( cpc.m_POF1.m_BarycentricCoords[0] * vec_v[vec_vid[0]]
                        + cpc.m_POF1.m_BarycentricCoords[1] * vec_v[vec_vid[1]]
                        + cpc.m_POF1.m_BarycentricCoords[2] * vec_v[vec_vid[2]]
                        + cpc.m_POF1.m_BarycentricCoords[3] * vec_v[vec_vid[3]] );
                Real m( cpc.m_POF1.m_BarycentricCoords[0] * m_Model.GetMass(vec_vid[0])
                        + cpc.m_POF1.m_BarycentricCoords[1] * m_Model.GetMass(vec_vid[1])
                        + cpc.m_POF1.m_BarycentricCoords[2] * m_Model.GetMass(vec_vid[2])
                        + cpc.m_POF1.m_BarycentricCoords[3] * m_Model.GetMass(vec_vid[3]) );
                Vec3 vt( v - mal::Dot(v,cpc.m_Normal)*cpc.m_Normal );
                Real norm_sq_vt( mal::NormSq(vt) );
                if( norm_sq_vt > 1e-8 )
                {
                    // Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN / m ); //Apply factor to convert cached impulse to v magnitude \todo NOT SURE about barycentric mass...
                    Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN );
                    if( delta_vn > 0 )
                    {
                        Real norm_vt( mal::Sqrt(norm_sq_vt) );
                        Real post_vt_over_norm_vt( mal::Max(0.0f,1.0f-m_Params.m_ContactSolver_DynamicFriction_Coeff*(delta_vn/norm_vt)) );
                        // Final v includes tangential friction but keeps normal component unchanged
                        Vec3 v_post( post_vt_over_norm_vt * vt + mal::Dot(v,cpc.m_Normal)*cpc.m_Normal );
                        Vec3 delta_v( v_post - v );
#ifdef __ENABLE_CONTACT_EFFECTIVE_MASS
                        // Enforce \Delta v distributing it among v1..v4 using effective mass at v
                        Real effective_mass( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * m_Model.GetInvMass(vec_vid[0])
                                                       + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * m_Model.GetInvMass(vec_vid[1])
                                                       + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * m_Model.GetInvMass(vec_vid[2])
                                                       + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * m_Model.GetInvMass(vec_vid[3]) ) );
                        vec_v[vec_vid[0]] += ( cpc.m_POF1.m_BarycentricCoords[0] * effective_mass * m_Model.GetInvMass(vec_vid[0]) ) * delta_v;
                        vec_v[vec_vid[1]] += ( cpc.m_POF1.m_BarycentricCoords[1] * effective_mass * m_Model.GetInvMass(vec_vid[1]) ) * delta_v;
                        vec_v[vec_vid[2]] += ( cpc.m_POF1.m_BarycentricCoords[2] * effective_mass * m_Model.GetInvMass(vec_vid[2]) ) * delta_v;
                        vec_v[vec_vid[3]] += ( cpc.m_POF1.m_BarycentricCoords[3] * effective_mass * m_Model.GetInvMass(vec_vid[3]) ) * delta_v;
#else
                        // Barycentric distribution \note If seen as an impulse, this is the CORRECT way to distribute it among nodes
                        //\todo I THINK de distribution should also account for node masses... if it's an impulse, it should
                        vec_v[vec_vid[0]] += cpc.m_POF1.m_BarycentricCoords[0] * delta_v;
                        vec_v[vec_vid[1]] += cpc.m_POF1.m_BarycentricCoords[1] * delta_v;
                        vec_v[vec_vid[2]] += cpc.m_POF1.m_BarycentricCoords[2] * delta_v;
                        vec_v[vec_vid[3]] += cpc.m_POF1.m_BarycentricCoords[3] * delta_v;
#endif
                    }
                }
            }
            else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
        }
    }
}

void LeafDSH_Solid3D_FEM::UpdateContacts_ActiveSet( const Vec3* vec_vel, const Vec3* vec_acc )
{
    DS_ASSERT( vec_vel != 0 );
#ifdef __ENABLE_CONTACT_ACTIVE_SET_IGNORE_ACCELERATION
#else
    DS_ASSERT( vec_acc != 0 );
#endif
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        ContactConstraintD2K& ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            ContactPointConstraintD2K& cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_POF1.m_FeatureId.IsVertex() )
            {
                geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
#ifdef __ENABLE_CONTACT_ACTIVE_SET_IGNORE_ACCELERATION
                cpc.m_IsActive = mal::Dot(vec_vel[vid],cpc.m_Normal) < 0;     //Collision
#else
                cpc.m_IsActive = mal::Dot(vec_vel[vid],cpc.m_Normal) < 0     //Collision
                                 || mal::Dot(vec_acc[vid],cpc.m_Normal) < 0; //Contact
#endif
            }
            else if( cpc.m_POF1.m_FeatureId.IsTetrahedron() )
            {
                geo::tetsolid3_feature_index_type eid = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                uint32 vec_vid[4];
                m_pMeshS->T_VecVID( eid, vec_vid );
                // Interpolate V_cp, A_cp from node v,a
                Vec3 v( cpc.m_POF1.m_BarycentricCoords[0] * vec_vel[vec_vid[0]]
                        + cpc.m_POF1.m_BarycentricCoords[1] * vec_vel[vec_vid[1]]
                        + cpc.m_POF1.m_BarycentricCoords[2] * vec_vel[vec_vid[2]]
                        + cpc.m_POF1.m_BarycentricCoords[3] * vec_vel[vec_vid[3]] );
#ifdef __ENABLE_CONTACT_ACTIVE_SET_IGNORE_ACCELERATION
                cpc.m_IsActive = mal::Dot(v,cpc.m_Normal) < 0;     //Collision
#else
                Vec3 a( cpc.m_POF1.m_BarycentricCoords[0] * vec_acc[vec_vid[0]]
                        + cpc.m_POF1.m_BarycentricCoords[1] * vec_acc[vec_vid[1]]
                        + cpc.m_POF1.m_BarycentricCoords[2] * vec_acc[vec_vid[2]]
                        + cpc.m_POF1.m_BarycentricCoords[3] * vec_acc[vec_vid[3]] );
                cpc.m_IsActive = mal::Dot(v,cpc.m_Normal) < 0     //Collision
                                 || mal::Dot(a,cpc.m_Normal) < 0; //Contact
#endif
            }
            else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
        }
    }
}
#endif //__S3D_ENABLE_CONTACT_D2K

#ifdef __S2_DS_ENABLE_PARAMS
void LeafDSH_Solid3D_FEM::Params::InitArchetype( util::ArchetypeLibrary &al )
{
    Params params;
    al.BeginArchetype( "Archetype_s2_ds_LDSH_FEMS3D_Params" );
    {
        al.BeginProperty_Group( "<Develop>" );
        {
            al.AddProperty( "run", params.m_bRun,
                            archetype_offset_of(params,m_bRun), util::IArchetypeInstance::NTPF_Ignore );

            const char *vec_names_km[] = { "Disabled", "Static", "Periodic", "Mouse" };
            uint32 vec_values_km[] = { (uint32)Params::eKM_Disabled,
                                       (uint32)Params::eKM_Static,
                                       (uint32)Params::eKM_Periodic,
                                       (uint32)Params::eKM_Mouse };
            al.AddProperty_Enum32( "kinematic", (uint32)params.m_KinematicMode,
                                   Params::cNumKM, vec_names_km, vec_values_km,
                                   archetype_offset_of(params,m_KinematicMode), IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<Environment>" );
        {
            al.AddProperty_NIR( "gravity", params.m_Gravity[1], -100.0f, 0.0f,
                                archetype_offset_of(params,m_Gravity[1]), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "air_drag", params.m_AirDragCoeff, 0.0f, 1.0f,
                                archetype_offset_of(params,m_AirDragCoeff), Params::RebuildAirDrag );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<Material>" );
        {
            al.AddProperty_NIR( "density", params.m_Density, 1.0f, 1000.0f,
                                archetype_offset_of(params,m_Density), Params::RebuildMass );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<FEM>" );
        {
            al.BeginProperty_Group( "<Develop>" );
            {
                // Material
                const char *vec_names_mm[] = { "L",
                                               "C_WRP", "C_LCM", "C_LCMH", "C_CCM",
                                               "H_LCM", "H_CCM", "H_NH0", "H_NH1" };
                uint32 vec_values_mm[] = { ms::fem::Params::eMM_L,
                                           // Corotational
                                           ms::fem::Params::eMM_C_WRP,
                                           ms::fem::Params::eMM_C_LCM,
                                           ms::fem::Params::eMM_C_LCMH,
                                           ms::fem::Params::eMM_C_CCM,
                                           // Hyperelastic
                                           ms::fem::Params::eMM_H_LCM,
                                           ms::fem::Params::eMM_H_CCM,
                                           ms::fem::Params::eMM_H_NH_C0,
                                           ms::fem::Params::eMM_H_NH_C1 };
                al.AddProperty_Enum32( "Material",
                                       (uint32)ms::fem::Params::eMM_L,
                                       (uint32)ms::fem::Params::cNumMM, vec_names_mm, vec_values_mm,
                                       archetype_offset_of(params,m_FEM.m_MM),
                                       Params::RebuildFEM );
                // Rotation
                const char *vec_names_rm[] = { "Id", "QR", "MSVD", "DAPD" };
                uint32 vec_values_rm[] = { ms::fem::Params::eRM_Id,
                                           ms::fem::Params::eRM_QR,
                                           ms::fem::Params::eRM_MSVD,
                                           ms::fem::Params::eRM_DAPD };
                al.AddProperty_Enum32( "Rotation",
                                       (uint32)ms::fem::Params::eRM_Id,
                                       (uint32)ms::fem::Params::cNumRM, vec_names_rm, vec_values_rm,
                                       archetype_offset_of(params,m_FEM.m_RM),
                                       Params::RebuildFEM );
                // Differentials
                const char *vec_names_dm[] = { "T", "E", "I", "N" };
                uint32 vec_values_dm[] = { ms::fem::Params::eDM_Truncated,
                                           ms::fem::Params::eDM_Exact,
                                           ms::fem::Params::eDM_Exact_Inv,
                                           ms::fem::Params::eDM_Numerical };
                al.AddProperty_Enum32( "Differential",
                                       (uint32)ms::fem::Params::eDM_Truncated,
                                       (uint32)ms::fem::Params::cNumDM, vec_names_dm, vec_values_dm,
                                       archetype_offset_of(params,m_FEM.m_DM),
                                       Params::RebuildFEM );

                al.AddProperty_NIR<float32>( "factor det(F)", 0.0f, 0.0f, 1.0f, archetype_offset_of(params,m_FEM.m_FactorDetF), Params::RebuildFEM );
                al.AddProperty_NIR<float32>( "ICF det(F)", 0.0f, 0.0f, 5.0f, archetype_offset_of(params,m_FEM.m_InvertedCompressibilityFactorDetF), Params::RebuildFEM );
                al.AddProperty_NIR<float32>( "eps.det(F)", 0.01f, 0.001f, 1.0f, archetype_offset_of(params,m_FEM.m_DegenerateThresholdDetF), Params::RebuildFEM );
                al.AddProperty_NIR<float32>( "ipol.det(F)", 0.0f, -5.0f, 0.0f, archetype_offset_of(params,m_FEM.m_ThresholdIpolDetF), Params::RebuildFEM );

                al.AddProperty_NIR<float32>( "PD_U L factor", 1.0f, 1.0f, 10.0f, archetype_offset_of(params,m_FEM.m_DAPD_L_Factor), Params::RebuildFEM );
                al.AddProperty_NIR<float32>( "PD_U NL factor", 0.0f, 0.0f, 10.0f, archetype_offset_of(params,m_FEM.m_DAPD_NL_Factor), Params::RebuildFEM );
                al.AddProperty_NIR<float32>( "PD_U NL exponent", 1.0f, 1.0f, 10.0f, archetype_offset_of(params,m_FEM.m_DAPD_NL_Exponent), Params::RebuildFEM );

                al.AddProperty_NIR<float32>( "ECIE e threshold", 0.4, 0.01, 0.99, archetype_offset_of(params,m_FEM.m_ECIE_e_threshold), Params::RebuildFEM );
                al.AddProperty_NIR<float32>( "ECIE k factor", 1.0f, 0.0f, 20.0f, archetype_offset_of(params,m_FEM.m_ECIE_k_factor), Params::RebuildFEM );
            }
            al.EndProperty_Group();

#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
            al.BeginProperty_Group( "<Test>" );
            {
                // Test
                const char *vec_names_tm[] = { "Disabled", "Random", "Plane", "Sphere" };
                uint32 vec_values_tm[] = { (uint32)Params::eTM_Disabled,
                                           (uint32)Params::eTM_Random,
                                           (uint32)Params::eTM_Plane,
                                           (uint32)Params::eTM_Sphere };
                al.AddProperty_Enum32( "test mode", (uint32)params.m_TestMode,
                                       Params::cNumTM, vec_names_tm, vec_values_tm,
                                       archetype_offset_of(params,m_TestMode), Params::RebuildTest );
                const char *vec_names_tn[] = { "All", "Boundary" };
                uint32 vec_values_tn[] = { (uint32)Params::eTN_All,
                                           (uint32)Params::eTN_Boundary };
                al.AddProperty_Enum32( "test nodes", (uint32)params.m_TestNodes,
                                       Params::cNumTN, vec_names_tn, vec_values_tn,
                                       archetype_offset_of(params,m_TestNodes), Params::RebuildTest );
                al.AddProperty_NIR( "test fraction", params.m_TestFraction, 0.0f, 1.0f,
                                    archetype_offset_of(params,m_TestFraction), IArchetypeInstance::NTPF_Ignore );
                // Ground
                al.AddProperty( "hacked ground?", params.m_bEnableHackedGroundPlane,
                                archetype_offset_of(params,m_bEnableHackedGroundPlane), IArchetypeInstance::NTPF_Ignore );
            }
            al.EndProperty_Group();
#endif

            al.BeginProperty_Group( "<Material>" );
            {
                al.AddProperty_NIR( "young_modulus", params.m_FEM.m_YoungModulus, 1.0f, 10000.0f,
                                    archetype_offset_of(params,m_FEM.m_YoungModulus), Params::RebuildFEM );
                al.AddProperty_NIR( "poisson_ratio", params.m_FEM.m_PoissonRatio, 0.0f, 0.49f,
                                    archetype_offset_of(params,m_FEM.m_PoissonRatio), Params::RebuildFEM );
                al.AddProperty_NIR( "damping_ratio", params.m_FEM.m_DampingRatio, 0.0f, 20.0f,
                                    archetype_offset_of(params,m_FEM.m_DampingRatio), Params::RebuildRayleighCoeffsFromFreqs );
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<Plasticity>" );
            {
                al.AddProperty_NIR( "yield", params.m_FEM.m_PlasticYield, 0.0f, 1.0f,
                                    archetype_offset_of(params,m_FEM.m_PlasticYield), Params::RebuildFEM );
                al.AddProperty_NIR( "max", params.m_FEM.m_PlasticMax, 0.0f, 10.0f,
                                    archetype_offset_of(params,m_FEM.m_PlasticMax), Params::RebuildFEM );
                al.AddProperty_NIR( "creep/s", params.m_FEM.m_PlasticCreepPerSecond, 0.0f, 100.0f,
                                    archetype_offset_of(params,m_FEM.m_PlasticCreepPerSecond), Params::RebuildFEM );
            }
            al.EndProperty_Group();

            al.BeginProperty_Group( "<Rayleigh>" );
            {
                al.AddProperty_NIR( "freq1", params.m_FEM.m_RayleighFreq1, 0.1f, 10.0f,
                                    archetype_offset_of(params,m_FEM.m_RayleighFreq1), Params::RebuildRayleighCoeffsFromFreqs );
                al.AddProperty_NIR( "freq2", params.m_FEM.m_RayleighFreq2, 10.0f, 20000.0f,
                                    archetype_offset_of(params,m_FEM.m_RayleighFreq2), Params::RebuildRayleighCoeffsFromFreqs );
                al.AddProperty_NIR( "coeff_m", params.m_FEM.m_RayleighCoeff[0], 0.0f, 10.0f,
                                    archetype_offset_of(params,m_FEM.m_RayleighCoeff[0]), Params::RebuildFEM );
                al.AddProperty_NIR( "coeff_k", params.m_FEM.m_RayleighCoeff[1], 0.0f, 1.0f,
                                    archetype_offset_of(params,m_FEM.m_RayleighCoeff[1]), Params::RebuildFEM );
            }
            al.EndProperty_Group();
        }
        al.EndProperty_Group();

        //---- Solvers
        al.BeginProperty_Group( "<Integrator>" );
        {
            const char *vec_names[] = { "SE", "QIE v", "QIE dx", "FIE dx", "FIQS dx" };
            uint32 vec_values[] = { (uint32)Params::eIT_SymplecticEuler,
                                    (uint32)Params::eIT_QuasiImplicitEuler_v,
                                    (uint32)Params::eIT_QuasiImplicitEuler_dx,
                                    (uint32)Params::eIT_FullyImplicitEuler_dx,
                                    (uint32)Params::eIT_FullyImplicitEuler_dx_QuasiStatic };
            al.AddProperty_Enum32( "integrator", (uint32)params.m_IntegratorType,
                                   Params::cNumIT, vec_names, vec_values,
                                   archetype_offset_of(params,m_IntegratorType), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "variable_dt", params.m_bEnableVariableDT,
                            archetype_offset_of(params,m_bEnableVariableDT), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "dt", params.m_FixedDT, 0.001f, 0.03333f,
                                archetype_offset_of(params,m_FixedDT), IArchetypeInstance::NTPF_Rebuild );
            al.AddProperty_NIR( "time scale", params.m_TimeScale, 0.01f, 1.0f,
                                archetype_offset_of(params,m_TimeScale), IArchetypeInstance::NTPF_Rebuild );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<SolverLS>" );
        {
            const char *vec_names[] = { "CG", "CR", "SQMR", "MINRES","GMRES",
                                        "CGN", "CGS",
                                        "RI",
                                        "GS", "Jacobi", "LU" };
            uint32 vec_values[] = { (uint32)Params::eLSST_CG,
                                    (uint32)Params::eLSST_CR,
                                    (uint32)Params::eLSST_SQMR,
                                    (uint32)Params::eLSST_MINRES,
                                    (uint32)Params::eLSST_GMRES,
                                    (uint32)Params::eLSST_CGN,
                                    (uint32)Params::eLSST_CGS,
                                    (uint32)Params::eLSST_RI,
                                    (uint32)Params::eLSST_GS,
                                    (uint32)Params::eLSST_Jacobi,
                                    (uint32)Params::eLSST_LU };
            al.AddProperty_Enum32( "method", (uint32)params.m_SolverLS_Type,
                                   cNumLSST, vec_names, vec_values,
                                   archetype_offset_of(params,m_SolverLS_Type), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "max_iter", params.m_SolverLS_MaxIter, uint32(1), uint32(500),
                                archetype_offset_of(params,m_SolverLS_MaxIter), IArchetypeInstance::NTPF_Ignore );

            al.AddProperty_NIR( "max_restart", params.m_SolverLS_MaxRestart, uint32(1), uint32(500),
            archetype_offset_of(params,m_SolverLS_MaxRestart), Params::RebuildSolverLS_MaxRestart );

            al.AddProperty( "warmstarting?", params.m_SolverLS_bEnableWarmstarting,
                            archetype_offset_of(params,m_SolverLS_bEnableWarmstarting), IArchetypeInstance::NTPF_Ignore );

            al.AddProperty_NIR( "log10(rel_eps)", params.m_SolverLS_Log10_RelEpsilon, -15.0f, -1.0f,
                                archetype_offset_of(params,m_SolverLS_Log10_RelEpsilon), Params::RebuildSolverLS_Epsilon );
            al.AddProperty_NIR( "log10(abs_eps)", params.m_SolverLS_Log10_AbsEpsilon, -15.0f, -1.0f,
                                archetype_offset_of(params,m_SolverLS_Log10_AbsEpsilon), Params::RebuildSolverLS_Epsilon );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<SolverNR>" );
        {
            al.AddProperty_NIR( "max_iter", params.m_SolverNR_MaxIter, uint32(1), uint32(100),
                                archetype_offset_of(params,m_SolverNR_MaxIter), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "log10(rel_eps)", params.m_SolverNR_Log10_RelEpsilon, -15.0f, -1.0f,
                                archetype_offset_of(params,m_SolverNR_Log10_RelEpsilon), Params::RebuildSolverNR_Epsilon );
            al.AddProperty_NIR( "log10(abs_eps)", params.m_SolverNR_Log10_AbsEpsilon, -15.0f, -1.0f,
                                archetype_offset_of(params,m_SolverNR_Log10_AbsEpsilon), Params::RebuildSolverNR_Epsilon );
        }
        al.EndProperty_Group();

        al.BeginProperty_Group( "<ConstraintSolver>" );
        {
            const char *vec_names_kncst[] = { "None", "Reaction" };
            const char *vec_names_ndcst[] = { "None", "Reaction" };
            uint32 vec_values_kncst[] = { (uint32)Params::eKNCST_None,
                                          (uint32)Params::eKNCST_Reaction };
            uint32 vec_values_ndcst[] = { (uint32)Params::eNDCST_None,
                                          (uint32)Params::eNDCST_Reaction };
            al.AddProperty_Enum32( "KNCST", (uint32)params.m_ConstraintSolver_KNCST,
                                   Params::cNumKNCST, vec_names_kncst, vec_values_kncst,
                                   archetype_offset_of(params,m_ConstraintSolver_KNCST), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_Enum32( "NDCST", (uint32)params.m_ConstraintSolver_NDCST,
                                   Params::cNumNDCST, vec_names_ndcst, vec_values_ndcst,
                                   archetype_offset_of(params,m_ConstraintSolver_NDCST), IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();

#ifdef __S3D_ENABLE_CONTACT_D2K
        al.BeginProperty_Group( "<ContactSolver>" );
        {
            const char *vec_names_pcst[] = { "None", "Penalty", "Alteration", "Hack" };
            const char *vec_names_vcst[] = { "None", "Penalty", "Reaction", "Hack" };
            uint32 vec_values_pcst[] = { (uint32)Params::ePCST_None,
                                         (uint32)Params::ePCST_Penalty,
                                         (uint32)Params::ePCST_Alteration,
                                         (uint32)Params::ePCST_Hack };
            uint32 vec_values_vcst[] = { (uint32)Params::eVCST_None,
                                         (uint32)Params::eVCST_Penalty,
                                         (uint32)Params::eVCST_Reaction,
                                         (uint32)Params::eVCST_Hack };
            al.AddProperty_Enum32( "PosCST", (uint32)params.m_ContactSolver_PCST,
                                   Params::cNumPCST, vec_names_pcst, vec_values_pcst,
                                   archetype_offset_of(params,m_ContactSolver_PCST), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_Enum32( "VelCST", (uint32)params.m_ContactSolver_VCST,
                                   Params::cNumVCST, vec_names_vcst, vec_values_vcst,
                                   archetype_offset_of(params,m_ContactSolver_VCST), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Depth Margin", params.m_ContactSolver_DepthMargin, 0.0f, 0.1f,
                                archetype_offset_of(params,m_ContactSolver_DepthMargin), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Penalty Ks", params.m_ContactSolver_Penalty_Ks, 0.0f, 10000.0f,
                                archetype_offset_of(params,m_ContactSolver_Penalty_Ks), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Penalty Kd", params.m_ContactSolver_Penalty_Kd, 0.0f, 100.0f,
                                archetype_offset_of(params,m_ContactSolver_Penalty_Kd), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Relaxation Coeff", params.m_ContactSolver_Relaxation_Coeff, 0.0f, 1.0f,
                                archetype_offset_of(params,m_ContactSolver_Relaxation_Coeff), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Restitution Coeff", params.m_ContactSolver_Restitution_Coeff, 0.0f, 1.0f,
                                archetype_offset_of(params,m_ContactSolver_Restitution_Coeff), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "S.Friction", params.m_ContactSolver_StatictFriction_Coeff, 0.0f, 2.0f,
                                archetype_offset_of(params,m_ContactSolver_StatictFriction_Coeff), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "D.Friction", params.m_ContactSolver_DynamicFriction_Coeff, 0.0f, 2.0f,
                                archetype_offset_of(params,m_ContactSolver_DynamicFriction_Coeff), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Max Depth", params.m_ContactSolver_MaxDepth, 0.0f, 10.0f,
                                archetype_offset_of(params,m_ContactSolver_MaxDepth), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "BoE", params.m_ContactSolver_BreakOnError,
                            archetype_offset_of(params,m_ContactSolver_BreakOnError), util::IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
#endif

        al.BeginProperty_Group( "<Viz>" );
        {
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_GeoObjectDDF, archetype_offset_of(params,m_GeoObjectDDF), IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
    }
    al.EndArchetype();
}
util::ArchetypeLibrary LeafDSH_Solid3D_FEM::Params::s_LeafDSH_Solid3D_FEM_Params_ArchetypeLibrary;
#endif


//---------------------------------------------------------------- TEMPORAL
void LeafDSH_Solid3D_FEM::NR_Init( Real dt )
{
    UpdateKNC( dt );
#ifdef __S3D_ENABLE_CONTACT_D2K
    if( m_Params.m_ContactSolver_PCST == Params::ePCST_Alteration )
        ApplyContacts_PositionAlteration_Smoothed( m_Model.GetVecPos_RW() );
#endif
    m_TES.BeginEvaluation( m_pTESC, m_vecPos0, m_Model.GetVecPos(), m_Model.GetVecVel(), m_TotalTime );
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_vecPos0[it_node] = m_Model.GetPos(it_node);
#ifdef __S3D_S2_DS_FEM_ENABLE_TEST
    UpdateTest(dt);
#endif
}

void LeafDSH_Solid3D_FEM::NR_End( Real dt )
{
    // Set final iterates x_k, v_k to Model state
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_Model.SetPos( it_node, m_IE_vec_x_k[it_node] );
        m_Model.SetVel( it_node, m_IE_vec_v_k[it_node] );
    }
    m_TES.EndEvaluation( m_pTESC );
    if( m_AirDragPerSecond > 0 )
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            m_Model.SetVel( it_node, (Real(1)-dt*m_AirDragPerSecond)*m_Model.GetVel(it_node) );
    RecomputeEnergy();
    S2_DS_TRACE_END_SLICE(true);
}

void LeafDSH_Solid3D_FEM::NR_PreStep( Real dt )
{
    DS_ASSERT( m_Params.m_SolverLS_Type != Params::eLSST_CGN ); //TEMP: NOT YET IMPLEMENTED!!

    S2_DS_STAT_SET( m_Stats.m_SolverLS_NumIter, 0 );

    bool bQuasistatic( dt == 0 );
    Real inv_dt( bQuasistatic ? 0 : mal::Rcp(dt) );

    // Init x_k^{n+1} = x^n, v_k^{n+1} = v^n at nr-iteration k = 0
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_IE_vec_x_k[it_node] = m_Model.GetPos(it_node);
        m_IE_vec_v_k[it_node] = m_Model.GetVel(it_node);
    }
    m_NR_nr_k = 0;
}

void LeafDSH_Solid3D_FEM::NR_MidStep( Real dt )
{
    bool bQuasistatic( dt == 0 );
    Real inv_dt( bQuasistatic ? 0 : mal::Rcp(dt) );

    // Begin NR iteration
    //\note BeginEvaluation() is open from FixedStep()

    //---- Assemble b_k (\todo All m_TES f could be done in a SINGLE PASS)
    // Reset b
    memset( m_LS_vec_b, 0, sizeof(Vec3) * m_NumNodes );
    // Accumulate External forces (Gravity)
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_LS_vec_b[it_node] += m_Model.GetMass(it_node)*m_Params.m_Gravity;
#ifdef __S3D_ENABLE_CONTACT_D2K
    // Accumulate penalty forces
    if( m_Params.m_ContactSolver_PCST == Params::ePCST_Penalty
        || m_Params.m_ContactSolver_VCST == Params::eVCST_Penalty )
        ApplyContacts_Penalty( m_LS_vec_b );
#endif
    // Accumulate internal forces (f_e, f_d, f_p) at ( x_k^{n+1}, v_k^{n+1} )
    m_TES.f_e( m_pTESC, m_LS_vec_b );
    //m_TES.f_d( m_pTESC, m_LS_vec_b ); //\todo This should be included, according to FEM.pdf, BUT adds noticeable damping...
    m_TES.f_p( m_pTESC, m_LS_vec_b );
    // Accumulate (1/dt) * M (v^n - v_k^{n+1}), subs v_k^{n+1} with (x_k^{n+1} - x^n)/dt
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_LS_vec_b[it_node] += inv_dt * m_Model.GetMass(it_node) * ( m_Model.GetVel(it_node) - inv_dt * (m_IE_vec_x_k[it_node] - m_Model.GetPos(it_node) ) );

    //---- Solve A_k y_k = b_k
    // Warmstarting only the FIRST NR-iter with \Delta x_0 = y = dt*vel_0 may be too optimistic and overshoot, but seems to work well. Further iterations use y = 0
    if( m_NR_nr_k == 0 && m_Params.m_SolverLS_bEnableWarmstarting ) for( unsigned int i=0; i<m_NumNodes; i++ ) m_LS_vec_y[i] = dt*m_Model.GetVel(i);
    else memset( m_LS_vec_y, 0, sizeof(Vec3) * m_NumNodes );

    //TEMP: Random initial guess avoids breakdown, according to: "LANCZOS METHODS FOR THE SOLUTION OF NONSYMMETRIC SYSTEMS OF LINEAR EQUATIONS*"
    if( m_Params.m_SolverLS_Type == Params::eLSST_CGS ) for( unsigned int i=0; i<m_NumNodes; i++ ) m_LS_vec_y[i] = Vec3( mal::Random01<Real>(), mal::Random01<Real>(), mal::Random01<Real>() );

#ifdef __S3D_ENABLE_CONTACT_D2K
    /*\todo THIS changes the active set of each NR iteration,
      which is basically a QP-like approach (see Faure's
      QPContacts)... consider pre-computing the active set only
      once OUTSIDE the loop, may have better convergence properties
      ("2008 Globally coupled collision handling using volume
      preserving impulses" performs N outer iterations with
      changing active-set and M inner iterations for each one.
    */
    // Contacts
    if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
    {
        UpdateContacts_ActiveSet( m_IE_vec_v_k, //m_Model.GetVecVel(),
                                  m_LS_vec_b );
        ApplyContacts_NormalImpulse( m_LS_vec_b, dt ); //\todo This ONLY fixed b, but did not acc normal impulses: ApplyContacts_Reaction( m_LS_vec_b );
        ApplyContacts_Reaction( m_LS_vec_y );
    }
#endif
    /*\todo Non-degeneration
      if( m_Params.m_ConstraintSolver_NDCST == Params::eNDCST_Reaction )
      {
      m_TES.UpdateConstraints_ActiveSet( m_pTESC, m_Model.GetVecVel(), m_LS_vec_b );
      m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, m_LS_vec_b );
      m_TES.ApplyConstraints_NoC_Reaction( m_pTESC, m_LS_vec_y );
      }
    */

    // Apply KNC to b_k and y_k AFTER other constraints (and therefore has highest priority)

#ifdef __S3D_ENABLE_KNC_IN_LS
    //LS_ApplyPosKNC( m_LS_vec_y, m_LS_vec_b, m_IE_vec_x_k );
    if( m_Params.m_KinematicMode != Params::eKM_Disabled
        && m_Params.m_ConstraintSolver_KNCST == Params::eKNCST_Reaction )
    {
        for( ms::KinematicNodeConstraintSet3::PoolKNC::iterator it_knc=m_KNCS.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int nid( it_knc->m_NID );
            // If it's a position constraint, Compute vel from position
            // error, otherwise, it's already a velocity constraint
            if( it_knc->m_Type == ms::KinematicNodeConstraintSet3::KNC::eKNCT_Position )
            {
                Vec3 vk( ( it_knc->m_Pos - m_Model.GetPos(nid) ) * inv_dt );
                Real norm_vk( vk.Norm() );
                if( norm_vk > Real(0.001f) )
                    it_knc->m_Vel = mal::Min( norm_vk, 50.0f ) * vk.Normalized();
                else
                    it_knc->m_Vel = Vec3::Zero();
            }
            // Set b_knc to KNC \Delta x
            m_LS_vec_b[nid] = it_knc->m_Pos - m_IE_vec_x_k[nid]; //it_knc->m_Vel; \todo NOT SURE
            // Set y_knc to KNC \Delta x
            m_LS_vec_y[nid] = it_knc->m_Pos - m_IE_vec_x_k[nid]; //it_knc->m_Vel; \todo NOT SURE
        }
    }
#endif //__S3D_ENABLE_KNC_IN_LS

    // Check convergence b_k < epsilon \IMPORTANT: AFTER applying KNC, otherwise b is 0 for stable configurations REGARDLESS of KNC
    // {
    //     const unsigned int cN( 3 * m_NumNodes );
    //     const Real *b( reinterpret_cast<const Real*>( m_LS_vec_b ) );
    //     nr_delta_new = ns::real_array::norm_2( b, cN );
    //     if( nr_delta_0 < 0 ) nr_delta_0 = nr_delta_new;
    // }

    //\todo Convergence analysis disabled for MultiDSH, consider restoring it locally or globally...
    // m_NR_bConverged1 = nr_delta_new < m_Params.m_SolverNR_AbsEpsilon; //\todo HACKED VERSION that checks only y_k = \Delta x correction instead of RHS b_k+1
    // if( !m_NR_bConverged1 )
    {
        Real ls_rel_prec( m_Params.m_SolverLS_RelEpsilon );
        Real ls_abs_prec( m_Params.m_SolverLS_AbsEpsilon );
        int ls_iter( m_Params.m_SolverLS_MaxIter );
        const Real factor_M( mal::Sq(inv_dt) + inv_dt * m_Params.m_FEM.m_RayleighCoeff[0] );
        const Real factor_K( Real(1) + inv_dt * m_Params.m_FEM.m_RayleighCoeff[1] );
        //CG_Solve( factor_M, factor_K, m_LS_vec_y, m_LS_vec_b, ls_rel_prec, ls_abs_prec, ls_iter );
        LS_Solve( factor_M, factor_K, m_LS_vec_y, m_LS_vec_b, ls_rel_prec, ls_abs_prec, ls_iter );
        S2_DS_STAT_ADD( m_Stats.m_SolverLS_NumIter, ls_iter );
        S2_DS_STAT_SET( m_Stats.m_SolverLS_RelPrec, ls_rel_prec );
        S2_DS_STAT_SET( m_Stats.m_SolverLS_AbsPrec, ls_abs_prec );

#ifdef __S3D_ENABLE_CONTACT_D2K
        // Apply friction to candidate vel BEFORE position integration... this may be wrong...
        if( m_Params.m_ContactSolver_VCST == Params::eVCST_Hack )
        {
            UpdateContacts_ActiveSet( m_LS_vec_y,
                                      0 ); //IMPORTANT: crash if activeset tries to use acceleration, it's forbidden, live with it
            ApplyContacts_NormalImpulse( m_LS_vec_y, 1 ); //y units = \Delta x => NormalImpulse = NormalDisplacement added to \Delta x
            ApplyContacts_FrictionImpulse( m_LS_vec_y, 1 ); //y units = \Delta x => FrictionImpulse = TangentialDisplacement removed from \Delta x
        }
        else if( m_Params.m_ContactSolver_VCST == Params::eVCST_Reaction )
        {
            ApplyContacts_FrictionImpulse( m_LS_vec_y, dt );
        }
#endif

    }
}

void LeafDSH_Solid3D_FEM::NR_PostStep( Real dt )
{
    bool bQuasistatic( dt == 0 );
    Real inv_dt( bQuasistatic ? 0 : mal::Rcp(dt) );

    // Update iterates
    {
        const unsigned int cN( 3 * m_NumNodes );
        const Real *y( reinterpret_cast<Real*>( m_LS_vec_y ) );
        Real *x_k( reinterpret_cast<Real*>( m_IE_vec_x_k ) );
        Real *v_k( reinterpret_cast<Real*>( m_IE_vec_v_k ) );
        ns::real_array::addeq( x_k, y, cN );
        //TEMP: From Vega, THIS is the correct way to do it!!
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        {
            m_IE_vec_v_k[it_node] = inv_dt * (m_IE_vec_x_k[it_node] - m_Model.GetPos( it_node ));
        }
        //This is what SiggCourse suggests, I think... ns::real_array::addeq_scaled( v_k, inv_dt, y, cN );
        // nr_delta_new = ns::real_array::norm_2( y, cN );
        // if( nr_delta_0 < 0 ) nr_delta_0 = nr_delta_new;
    }

    // NR iteration ends
    m_NR_nr_k++;

    //bConverged2 = nr_delta_new < m_Params.m_SolverNR_AbsEpsilon; TEMP THIS SHOULD NOT BE COMMENTED... TEST IT, I don't know WHY I commented it in the first place
    // if( !bConverged2 )
    {
        // If not finished, End current Eval and Begin a new one for the potential displacement x_0 -> x_k+1
        m_TES.EndEvaluation( m_pTESC );
        m_TES.BeginEvaluation( m_pTESC, m_Model.GetVecPos(), m_IE_vec_x_k, m_IE_vec_v_k, m_TotalTime ); //\note This does NOT update plasticity, as m_TotalTime has not advanced from last invokation
        //\note If no more iter, the currently open BeginEvaluation() will be closed at the end of FixedStep()
    }
}

} } // namespace S2::ds
