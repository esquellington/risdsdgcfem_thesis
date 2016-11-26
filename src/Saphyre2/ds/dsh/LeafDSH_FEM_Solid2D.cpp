#include "LeafDSH_FEM_Solid2D.h"

#include <Saphyre2/ds/dsh/IDynamicSystemHierarchy.h>
#include <Saphyre2/ds/dsh/IGeom.h>

#include <Mal/GConversion.h>
#include <Mal/GMatDecomposition.h>

#include <Geo/bp/Proxy.h>
#include <Geo/bp/Pair.h>

#include <util/ItemStream.h>

#ifdef __ENABLE_EIGENSTUFF
#  include <Mal/GSerialization.h>
#  include <Mal/GMatUtils.h>
#endif //__ENABLE_EIGENSTUFF

#ifdef __ENABLE_CONTINUOUS
#  include <Mal/GSolvePolynomialEq.h>
#endif

namespace S2 { namespace ds {

//TEMPORAL: Move to GMatOps.h when it works, and add some tests in testMal for randomly generated matrices!!
inline Mat2x2 RLerp( const Mat2x2 &R0, const Mat2x2 &R1, Real lambda01 )
{
    Real w1( lambda01 );
    Real w0( Real(1)-w1 );
    Mat2x2 LR( w0*R0 + w1*R1 );
    return GRotation2x2_PolarDecomposition( LR, mal::Det(LR) ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix
}

inline Mat2x2 RCerp( const Mat2x2 &R0, const Mat2x2 &R1, Real lambda01 )
{
    Real w1( mal::Sq(lambda01) * (3-2*lambda01) ); //h01 polynomial from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
    Real w0( Real(1) - w1 );
    Mat2x2 CR( w0*R0 + w1*R1 );
    return GRotation2x2_PolarDecomposition( CR, mal::Det(CR) ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix
}
//TEMPORAL

LeafDSH_FEM_Solid2D::LeafDSH_FEM_Solid2D( uint32 uid, IDynamicSystemHierarchy *p_parent )
: GLeafDSH_Model< ms::ParticleSystem2D >(uid,p_parent) //TEMPORAL should be ms::FEM_Solid2D
, m_NumElements(0)
, m_vecElements(0)
, m_NumNodes(0)
, m_vecTe2r(0)
, m_vecF(0)
, m_vecDetF(0)
#ifdef __ENABLE_SAVE_FORCES
, m_vecForces(0)
#endif
, m_Params(this)
#ifdef __ENABLE_KNC
, m_poolKNC(16)
#endif
  // Derived params
, m_TotalArea(0)
, m_DensityPerArea(0)
, m_AirDragPerSecond(0)
  // State
, m_TotalTime(0)
, m_AccTime(0)
  // Energy
, m_ElasticEnergy(0)
, m_PotentialEnergy(0)
, m_KineticEnergy(0)
#ifdef __ENABLE_CONTINUOUS
, m_poolDED(8)
, m_vecPos0(0)
#endif
#ifdef __ENABLE_CONTACT_D2K
, m_poolCCD2K(8)
#endif
, m_DDF( eDraw_Default )
, m_GeoObjectDDF( geo::eODDF_Shape | geo::eSDDF_Boundary | geo::eSDDF_Vertices | geo::eSDDF_Interior )
, m_ContactDDF( geo::eContactDataDDF_Points | geo::eContactDataDDF_Normals )
{

#ifdef __S2_DS_ENABLE_PARAMS
    if( Params::s_LeafDSH_FEM_Solid2D_Params_ArchetypeLibrary.IsEmpty() )
    {
        Params::InitArchetype( Params::s_LeafDSH_FEM_Solid2D_Params_ArchetypeLibrary );
    }
#endif

}

LeafDSH_FEM_Solid2D::~LeafDSH_FEM_Solid2D()
{
    if( m_vecElements ) delete[] m_vecElements;
#ifdef __ENABLE_CONTINUOUS
    if( m_vecPos0 ) delete[] m_vecPos0;
#endif
    if( m_vecTe2r ) delete[] m_vecTe2r;
    if( m_vecF ) delete[] m_vecF;
    if( m_vecDetF ) delete[] m_vecDetF;
#ifdef __ENABLE_SAVE_FORCES
    if( m_vecForces ) delete[] m_vecForces;
#endif
}

void LeafDSH_FEM_Solid2D::SyncGO()
{
    if( 0 != base_type::m_pGO )
    {
        //\note Assumes that transform == identity, as all
        //information is in the SDOF
        m_pGO->SetTransform( Transform2::Identity() );
        m_Model.GetDOF( static_cast< geo::GObjectSDOF<2,Vec2>* >(m_pGO)->GetVecDOF_WriteOnly() );
    }
}

void LeafDSH_FEM_Solid2D::RecomputeBV()
{
    if( 0 != base_type::m_pGO && 0 != base_type::m_pBV )
    {
        base_type::m_pGO->ComputeBV( *base_type::m_pBV );
        //\todo Maybe use polimorphyc BV::Extend()...
        const Real cNearDist(1.0f); //\todo this should come from np::Context::m_Stochastic_Near/FarDist, probably a factor of it
        switch( base_type::m_pBV->GetType() )
        {
        case geo::bv::eBV_AABB2: static_cast<geo::bv::AABB2*>(base_type::m_pBV)->Extend( cNearDist ); break;
        case geo::bv::eBV_Sphere2: static_cast<geo::bv::Sphere2*>(base_type::m_pBV)->Extend( cNearDist ); break;
        case geo::bv::eBV_LSS2: static_cast<geo::bv::LSS2*>(base_type::m_pBV)->Extend( cNearDist ); break;
        default: DS_ASSERT(false); break;
        }
    }
}

bool LeafDSH_FEM_Solid2D::Create( const ParamIt &pit )
{
    m_Params.m_FixedDT = pit.Find("fixed_dt").Get<Real>();

    ParamIt shape_def = pit.Find("shape_def");
    m_pMeshO = static_cast<geo_object_type*>( geo::ObjectFactory::STATIC_CreateES( shape_def ) );
    m_pMeshO->ResetDOF();
    m_pMeshS = m_pMeshO->GetShape();

    // Init particle system with fake uniform mass
    m_NumNodes = m_pMeshS->GetNumV();
    m_TotalMass = Real(666);
    m_Model.SetDscr( m_NumNodes, m_TotalMass );
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_Model.SetPos( it_node, m_pMeshO->GetSDOF(it_node) );
        m_Model.SetVel( it_node, Vec2(0,0) );
    }

    // Get material params
    m_Params.m_Thickness = pit.Find("thickness").Get<Real>();
    m_Params.m_Density = pit.Find("density").Get<Real>();
    m_Params.m_YoungModulus = pit.Find("young_modulus").Get<Real>();
    m_Params.m_PoissonRatio = pit.Find("poisson_ratio").Get<Real>();
    m_Params.m_DampingRatio = pit.Find("damping_ratio").Get<Real>();

    // Create FEM LT
    m_TotalArea = 0;
    m_NumElements = m_pMeshS->GetNumP();
    m_vecElements = new ms::fem::LinearTriangle2[m_NumElements];
    for( unsigned int it_p=0; it_p<m_pMeshS->GetNumP(); it_p++ )
    {
#ifdef __ENABLE_TRACE_FEM
        std::cout << "------ LT " << it_p << std::endl;
#endif
        unsigned int it_he( m_pMeshS->P_FirstHEID(it_p) );
        unsigned int vid0 = m_pMeshS->HE_OriginVID(it_he); it_he = m_pMeshS->HE_Next(it_he);
        unsigned int vid1 = m_pMeshS->HE_OriginVID(it_he); it_he = m_pMeshS->HE_Next(it_he);
        unsigned int vid2 = m_pMeshS->HE_OriginVID(it_he); it_he = m_pMeshS->HE_Next(it_he);
        DS_PARANOID_ASSERT( it_he == m_pMeshS->P_FirstHEID(it_p) );
        Vec2 r1 = m_pMeshS->V_Pos_0( vid0 );
        Vec2 r2 = m_pMeshS->V_Pos_0( vid1 );
        Vec2 r3 = m_pMeshS->V_Pos_0( vid2 );
        ms::fem::LinearTriangle2 &element( m_vecElements[it_p] );
        element.Init( vid0, vid1, vid2, r1, r2, r3, m_Params.m_YoungModulus, m_Params.m_PoissonRatio );
        m_TotalArea += element.Area();
    }

    // Precompute static VID
#ifdef __ENABLE_STATIC_ELEMENTS
#  ifdef __ENABLE_KNC
    for( unsigned int it_e=0; it_e < m_NumElements && it_e<cNumStaticElements; it_e++ )
        for( int i=0; i<3; i++ )
        {
            uint16 nid( m_vecElements[it_e].nid(i) );
            if( !IsKNC( nid ) )
                *m_poolKNC.New() = KinematicNodeConstraint( nid, m_pMeshS->V_Pos_0(nid) );
        }
#  endif
#else
#  ifdef __ENABLE_KNC
    if( cNumStaticNodes > 0 )
        *m_poolKNC.New() = KinematicNodeConstraint( 0, m_pMeshS->V_Pos_0(0) );
    if( cNumStaticNodes > 1 )
        *m_poolKNC.New() = KinematicNodeConstraint( m_pMeshS->GetNumV()-1, m_pMeshS->V_Pos_0(m_pMeshS->GetNumV()-1) );
#  endif
#endif //__ENABLE_STATIC_ELEMENTS

#ifdef __ENABLE_CONTINUOUS
    m_vecPos0 = new Vec2[m_NumNodes];
    // Init previous configuration C_i-1 = C_0
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_vecPos0[it_node] = m_Model.GetPos(it_node);
#endif

    // Init F and Transforms
    m_vecTe2r = new Transform2[ m_NumElements ];
    m_vecF = new Mat2x2[ m_NumElements ];
    m_vecDetF = new Real[ m_NumElements ];
    Compute_F( m_vecF, m_vecDetF );
    Compute_T( m_vecTe2r );

    // Rebuild everything except LT
    RebuildMass();
    RebuildAirDrag();
    RebuildEigenstuff();

    // Init RPAP
    m_LastRPAP.m_Pos = Vec2(0,0);
    m_LastRPAP.m_Radius = 0;
    m_LastRPAP.m_Pressure = 0;

    // Init CoM position
    if( pit.Find("pos_com").IsValid() ) SetPosCoM( pit.Find("pos_com").Get<Vec2>() );

    // Init object-size based params
    m_pMeshS->ComputeBVD( m_AABB0, geo::Transform2::Identity(), 0 );
    m_Params.m_ContactSolver_MaxDepth = mal::Max( m_AABB0.GetHalfSizes() );

#ifdef __ENABLE_SAVE_FORCES
    m_vecForces = new Vec2[ m_NumNodes ];
    for( unsigned int i=0; i<m_NumNodes; i++ ) m_vecForces[i] = Vec2::Zero();
#endif

    return true;
}

void LeafDSH_FEM_Solid2D::Step( Real dt )
{
    if( !m_Params.m_bRun ) return;

    Real scaled_dt( m_Params.m_TimeScale * dt );
    m_TotalTime += scaled_dt;

#ifdef __S2_DS_ENABLE_STATS
    m_Stats.Begin();
#endif

    if( m_Params.m_bEnableVariableDT )
    {
        m_AccTime += dt;
        Real stable_dt( Real(2) / mal::Sqrt(m_vecEigenValue[0]) );
        stable_dt *= 0.95f;
        int num_substeps( mal::Ceil( m_AccTime / stable_dt ) );
        if( num_substeps > 0 )
        {
            Real variable_dt( m_AccTime / num_substeps );
            DS_ASSERT( variable_dt < stable_dt );
            for( int i=0; i<num_substeps; i++ )
                FixedStep(variable_dt);
            m_Params.m_FixedDT = variable_dt;
        }
        m_AccTime = 0;
    }
    else
    {
        //\todo FixedDT THIS SHOULD BE PROBABLY HANDLED IN THE DSH OR THE SIMULATION SCHEME, NOT HERE
        m_AccTime += scaled_dt;
        Real scaled_fixed_dt( m_Params.m_TimeScale * m_Params.m_FixedDT );
        while( m_AccTime > scaled_fixed_dt )
        {
            FixedStep( scaled_fixed_dt );
            m_AccTime -= scaled_fixed_dt;
        }
    }

#ifdef __S2_DS_ENABLE_STATS
    m_Stats.End();
#endif
}

void LeafDSH_FEM_Solid2D::FixedStep( Real dt )
{
    /* Here we have previous configuration C_i-1 (m_vecPos0) and the
       current C_i (m_Model::GetPos), which MAY not fulfill current
       KNC set after the last FixedStep() call. These will NOT be
       enforced here, but during integration, and only fulfilled at
       the end of the timestep.
    */

    // Afterwards, we compute F,T due to the trajectory C_i-1 => C_i^+
    Compute_F( m_vecF, m_vecDetF );

    UpdateDED();

    Compute_T( m_vecTe2r );

#ifdef __ENABLE_CONTINUOUS
    // Save C_i into C_i-1 for the next frame
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_vecPos0[it_node] = m_Model.GetPos(it_node);
#endif

#ifdef __S2_DS_FEM_ENABLE_TEST
    UpdateTest(dt);
#endif

    // Integrate C_i^+ up to C_i+1
    switch( m_Params.m_IntegratorType )
    {
    case Params::eIT_ImplicitEuler:
        Solve_ImplicitEuler(dt);
        break;
    case Params::eIT_SymplecticEuler:
        Solve_SymplecticEuler(dt);
        break;
    default:
        DS_ASSERT(false);
    }

    // Exponential air drag
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.SetVel( it_node, (Real(1)-dt*m_AirDragPerSecond)*m_Model.GetVel(it_node) );

    // Set mesh DOF
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_pMeshO->GetSDOF(it_node) = m_Model.GetPos(it_node);

    // Final energy
    RecomputeEnergy();
}

bool LeafDSH_FEM_Solid2D::Edit( const ParamIt &pit )
{
    // Radial pressure
    ParamIt pit2 = pit.Find("pressure");
    if( pit2.IsValid() )
    {
        m_LastRPAP = pit2.Get< RadialPressureAtPoint2D >();
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        {
            Vec2 r( m_Model.GetPos(it_node)-m_LastRPAP.m_Pos );
            Real dist = r.Norm();
            if( dist < m_LastRPAP.m_Radius )
                m_Model.ApplyForce( it_node,
                                    (r / dist) * (1-dist/m_LastRPAP.m_Radius) * m_LastRPAP.m_Pressure );
        }
    }
    // CoM position
    pit2 = pit.Find("pos_com");
    if( pit2.IsValid() ) SetPosCoM( pit2.Get<Vec2>() );
    return true;
}

bool LeafDSH_FEM_Solid2D::Internal( const ParamIt &pit, ReturnStream &rets )
{
    ParamIt create_kpc_it = pit.Find("create_kpc");
    if( create_kpc_it.IsValid() )
    {
        // Find closest node to solid_pos and create a KNC, knc_eid is a pointer to the new KNC

        //\todo if we ever support non-node constraints, kpc_eid could
        //be a ptr to a generic ConstraintPoint class with
        //specializations Node and Generic... whatever...
        unsigned int closest_nid = FindClosestNode( create_kpc_it.GetSubItem().Find("solid_pos").Get<Vec2>() );

        // Return existing or new KNC
        KinematicNodeConstraint *pKNC(0);
        for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid() && 0 == pKNC; ++it_knc )
            if( it_knc->m_nid == closest_nid )
                pKNC = it_knc.GetPtr();
        if( 0 == pKNC ) pKNC = m_poolKNC.New();

        // Init KNC
        *pKNC = KinematicNodeConstraint( closest_nid, create_kpc_it.GetSubItem().Find("world_pos").Get<Vec2>() );
        pKNC->m_Vel = create_kpc_it.GetSubItem().Find("world_vel").Get<Vec2>();
        machine_uint_type kpc_eid( reinterpret_cast<machine_uint_type>( pKNC ) );

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
        KinematicNodeConstraint *pKNC = reinterpret_cast<KinematicNodeConstraint*>( kpc_eid );
        ParamIt pos_it = edit_kpc_it.GetSubItem().Find("world_pos");
        if( pos_it.IsValid() )
        {
            //DS_LOG_WARNING( "Editing KNC Pos %lld", kpc_eid );
            pKNC->m_Pos = pos_it.Get<Vec2>();
            pKNC->m_Type = KinematicNodeConstraint::eKNCT_Position;
        }
        else
        {
            ParamIt vel_it = edit_kpc_it.GetSubItem().Find("world_vel");
            if( vel_it.IsValid() )
            {
                //DS_LOG_WARNING( "Editing KNC Vel %lld", kpc_eid );
                pKNC->m_Vel = vel_it.Get<Vec2>();
                pKNC->m_Type = KinematicNodeConstraint::eKNCT_Velocity;
            }
        }
        // No data to be returned
    }

    // destroy_kpc REQUIRES KNC pool instead of array
    ParamIt destroy_kpc_it = pit.Find("destroy_kpc");
    if( destroy_kpc_it.IsValid() )
    {
        machine_uint_type kpc_eid = destroy_kpc_it.GetSubItem().Find("kpc_eid").Get<machine_uint_type>();
        //DS_LOG_WARNING( "Destroying KNC %lld", kpc_eid );
        KinematicNodeConstraint *pKNC = reinterpret_cast<KinematicNodeConstraint*>( kpc_eid );
        m_poolKNC.Delete( pKNC );
        // No data to be returned
    }

    return true;
}

void LeafDSH_FEM_Solid2D::DoViz( util::VizStream &vs ) const
{
    if( m_DDF.Test( eDraw_Mesh ) )
        geo::VizObject( m_pMeshO, vs, m_GeoObjectDDF );
    // Highlight KNC
#ifdef __ENABLE_KNC
    if( m_DDF.Test( eDraw_KNC ) )
    {
        if( m_Params.m_KinematicMode != Params::eKM_Disabled )
        {
            for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
            {
                VIZ_POINT2( vs, it_knc->m_Pos, 10, Vec4f(1,1,1,1) );
                VIZ_POINT2( vs, m_Model.GetPos( it_knc->m_nid ), 10, Vec4f(1,0,0,1) );
                /*
                  char name[16];
                  sprintf( name, " %d ", it_knc->m_nid );
                  VIZ_POINT2_NAMED( vs, name, it_knc->m_Pos, 10, Vec4f(1,1,1,1) );
                */
            }
        }
    }
#endif

#ifdef __S2_DS_FEM_ENABLE_TEST
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
                    VIZ_SEGMENT2( vs, m_pMeshS->V_Pos_0( it_tp ), m_vec_TM_Random_TargetPos[it_tp], 1, Vec4f(1,0,1,0.5) );
                    VIZ_POINT2( vs, m_vec_TM_Random_TargetPos[it_tp], 5, Vec4f(1,0,1,1) );
                }
                break;
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            VIZ_SEGMENT2( vs, m_pMeshS->V_Pos_0( vid ), m_vec_TM_Random_TargetPos[vid], 1, Vec4f(1,0,1,0.5) );
                            VIZ_POINT2( vs, m_vec_TM_Random_TargetPos[vid], 5, Vec4f(1,0,1,1) );
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
            default: break;
            }
        }
        break;
    case Params::eTM_Plane:
        {
            Real half_height( 2*m_AABB0.GetHalfSizes()[1]*(Real(1)-m_Params.m_TestFraction) );
            VIZ_SEGMENT2( vs, Vec2(-5,half_height), Vec2(5,half_height), 1, Vec4f(1,0,1,0.5) );
            VIZ_SEGMENT2( vs, Vec2(-5,-half_height), Vec2(5,-half_height), 1, Vec4f(1,0,1,0.5) );
        }
        break;
    case Params::eTM_Sphere:
        {
            Real radius( 1.1f*m_AABB0.GetHalfSizes().Norm()*(Real(1)-m_Params.m_TestFraction) );
            VIZ_DISK2_NO_ROT( vs, Vec2(0,0), radius, Vec4f(1,0,1,0.5), util::eVizStyle_Wire );
        }
        break;
    default: break;
    }
#endif

    // RPAP
    if( m_DDF.Test( eDraw_RPAP ) )
        VIZ_DISK2_NO_ROT( vs, m_LastRPAP.m_Pos, m_LastRPAP.m_Radius, Vec4f(1,1,1,1), util::eVizStyle_Wire );

#ifdef __ENABLE_VIZ_DEGENERATE_ELEMENTS
    // Element viz: Degenerate or Extremely deformed elements
    if( m_DDF.Test( eDraw_Degenerate ) )
    {
        for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
        {
            const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
            // Gather VID
            Vec2 x0( m_Model.GetPos( element.nid(0) ) );
            Vec2 x1( m_Model.GetPos( element.nid(1) ) );
            Vec2 x2( m_Model.GetPos( element.nid(2) ) );
            float area = element.Compute_Area(x0,x1,x2);
            float rel_area_deformation = (area / element.Area()) - 1.0f;
            if( area < 0 ) // Inverted
            {
                VIZ_TRIANGLE2( vs, x0, x1, x2, 1, Vec4f(0.5,0,0.5,0.25), util::eVizStyle_Solid );
                /*TEMP: Te2r is also inverted
                if( m_Params.m_bEnableCorotational )
                {
                    // Compute Transform from element 2 world
                    Transform2 Te2w( m_vecTe2r[it_e] );
                    Te2w.m_Pos += element.Barycenter0();
                    // Element transform viz
                    VIZ_TRANSFORM2( vs, Te2w, 1.0f );
                }
                */
            }
            /*TEMP: non-degenerated not drawn to avoid clutter
            else// if( mal::Abs( rel_area_deformation ) > 0.1 ) //rel_area < 10% || rel_area > 190%
            {
                // relative area deformation: compressed -> resting -> dilated => red -> white -> blue
                VIZ_TRIANGLE2( vs, x0, x1, x2,
                               Vec4f(1-rel_area_deformation,1-mal::Abs(rel_area_deformation),1+rel_area_deformation,1),
                               util::eVizStyle_Solid );
            }
            */
        }
    }
#endif

    // Show CoM
    if( m_DDF.Test( eDraw_CoM ) )
    {
        VIZ_POINT2( vs, ComputePosCoM(), 10, Vec4f(0,0,1,1) );
    }

    // Show corot transforms
    if( m_DDF.Test( eDraw_Corotational ) )
    {
        for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
        {
            const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
            // Compute Transform from element 2 world
            Transform2 Te2w( m_vecTe2r[it_e] );
            Te2w.m_Pos += element.Barycenter0();
            // Element transform viz
            VIZ_TRANSFORM2( vs, Te2w, 0.33f );
            /* Local deformation viz
               VIZ_TRIANGLE2( vs,
               Te2w * ( element.r(0) - element.Barycenter0() ),
               Te2w * ( element.r(1) - element.Barycenter0() ),
               Te2w * ( element.r(2) - element.Barycenter0() ),
               Vec4f(0,1,0,1), util::eVizStyle_Solid );
            */
        }
    }

    // Draw material coords SVD
    if( m_DDF.Test( eDraw_Eigenbasis ) )
    {
        for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
        {
            const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );

            Mat2x2 U,Vt;
            Vec2 diag_S;
            Transform2 Te2w;

            // Draw eigenbasis of the undeformed elements => ortogonal deformation directions?!
            Mat2x2 P;
            P(0,0) = element.r(1).x() - element.r(0).x(); P(0,1) = element.r(2).x() - element.r(0).x();
            P(1,0) = element.r(1).y() - element.r(0).y(); P(1,1) = element.r(2).y() - element.r(0).y();
            mal::GSingularValueDecomposition_USVt( P, U, diag_S, Vt );
            Te2w.m_Rot = Vt.Transposed();
            Te2w.m_Pos = m_vecTe2r[it_e].m_Pos + element.Barycenter0();
            VIZ_TRANSFORM2( vs, Te2w, 0.25f );

            // Draw eigenbasis of the deformed elements (smaller)
            Vec2 x0( m_Model.GetPos( element.nid(0) ) );
            Vec2 x1( m_Model.GetPos( element.nid(1) ) );
            Vec2 x2( m_Model.GetPos( element.nid(2) ) );
            Mat2x2 F;
            element.ComputeF( x0, x1, x2, F );
            mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
            if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
            Te2w.m_Rot = Vt.Transposed();
            Te2w.m_Pos = m_vecTe2r[it_e].m_Pos + element.Barycenter0();
            VIZ_TRANSFORM2( vs, Te2w, 0.15f );
        }
    }

    // Draw node forces
    if( m_DDF.Test( eDraw_Forces ) )
    {
#ifdef __ENABLE_SAVE_FORCES
        switch( m_Params.m_IntegratorType )
        {
        case Params::eIT_ImplicitEuler:
            for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
                VIZ_SEGMENT2( vs, m_Model.GetPos(it_node), m_Model.GetPos(it_node) + m_Model.GetVel(it_node), 1.0f, Vec4f(0,1,0,1) );
            break;
        case Params::eIT_SymplecticEuler:
            for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
                VIZ_SEGMENT2( vs, m_Model.GetPos(it_node), m_Model.GetPos(it_node) + m_vecForces[it_node], 1.0f, Vec4f(0.25,0.25,1,1) );
            break;
        default: break;
        }
#endif
    }

#ifdef __ENABLE_CONTINUOUS
    if( m_DDF.Test( eDraw_Trajectory ) )
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
            VIZ_SEGMENT2( vs, m_vecPos0[it_node], m_Model.GetPos(it_node), 2.0f, Vec4f(0.25,0.25,0.75,1) );
    if( m_DDF.Test( eDraw_DED ) )
    {
        for( PoolDED::iterator it_ded = m_poolDED.Begin(); it_ded.IsValid(); ++it_ded )
        {
            Vec2 point_noc( m_Model.GetPos( m_vecElements[it_ded->m_ElementId].nid(it_ded->m_NoC) ) );
            Vec2 point_noc_opposite( 0.5f*( m_Model.GetPos( m_vecElements[it_ded->m_ElementId].nid((it_ded->m_NoC+1)%3) ) + m_Model.GetPos( m_vecElements[it_ded->m_ElementId].nid((it_ded->m_NoC+2)%3) ) ) );
            VIZ_SEGMENT2( vs, point_noc, point_noc_opposite, 2.0f, Vec4f(1.0,0.5,0.5,1) );
            //VIZ_POINT2( vs, point_noc, 10.0f, Vec4f(1.0,0.25,0.25,1) );
            int ndn( FindNDN(it_ded->m_ElementId) );
            if( ndn != geo::cInvalidFeatureIndex )
            {
                int noc( (ndn+2)%3 );
                Vec2 point_ndn( m_Model.GetPos( m_vecElements[it_ded->m_ElementId].nid( noc ) ) );
                VIZ_POINT2( vs, point_ndn, 10.0f, Vec4f(0.25,1.0,0.25,1) );
            }
        }
    }
#endif

#ifdef __ENABLE_CONTACT_D2K
    if( m_DDF.Test( eDraw_Contacts ) )
        for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
            geo::VizContactData( it_ccd2k->m_TmpCD, vs, m_ContactDDF );
#endif
}

void LeafDSH_FEM_Solid2D::QueryStats( Flags32 flags, ReturnStream &rets ) const
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
            rets.Write( "elastic", Pair_f32_f32( m_TotalTime, E ) );//E );
            rets.Write( "total", T+V+E );

            //TEMP: PS energy is QUITE similar to FEM energy and oscilates comparably...
            Real T1( m_Model.ComputeKineticEnergy() );
            Real V1( m_Model.ComputePotentialEnergy( m_Params.m_Gravity ) );
            rets.Write( "k1", T1 );
            rets.Write( "p1", V1 );
            rets.Write( "t1", T1+V1+E );
        }
        rets.EndComplex();

        rets.BeginComplex( "<Misc>", eType_Property_Group );
        {
            rets.Write( "total_mass", m_Model.GetMass() );
            rets.Write( "total_area", m_TotalArea );
            rets.Write( "airdrag/s", m_AirDragPerSecond );
#ifdef __ENABLE_EIGENSTUFF
            for( int i=0; i<cNumEigenValues; i++ )
            {
                rets.Write( "freq[0]", (float)mal::Sqrt(m_vecEigenValue[i]) );
                rets.Write( "dt[0]",  float(2)/mal::Sqrt(m_vecEigenValue[i]) );
            }
#endif //__ENABLE_EIGENSTUFF
        }
        rets.EndComplex();

        rets.BeginComplex( "<RelArea>", eType_Property_Group );
        {
            rets.Write( "min", m_Stats.m_RelArea.m_Min );
            rets.Write( "max", m_Stats.m_RelArea.m_Max );
            rets.Write( "median", m_Stats.m_RelArea.ComputeMedian() );
            rets.Write( "avg", m_Stats.m_RelArea.m_Avg );
            rets.Write( "stdev", m_Stats.m_RelArea.m_Stdev );
        }
        rets.EndComplex();

        rets.BeginComplex( "<SolverLS>", eType_Property_Group );
        {
            rets.Write( "#iter", m_Stats.m_SolverLS_NumIter );
            rets.Write( "prec", m_Stats.m_SolverLS_Prec );
        }
        rets.EndComplex();

        rets.BeginComplex( "<Degeneration>", eType_Property_Group );
        {
            //rets.Write( "#degenerate", m_Stats.m_Num_Degenerate );
            rets.Write( "#degenerate", Pair_f32_f32( m_TotalTime, (float32)m_Stats.m_Num_Degenerate ) );
            rets.Write( "Sum det(F) < 0", Pair_f32_f32( m_TotalTime, m_Stats.m_Sum_Degenerate_DetF ) );
        }
        rets.EndComplex();
    }
    rets.EndComplex();
#endif //__S2_DS_ENABLE_STATS
}

void LeafDSH_FEM_Solid2D::QueryParams( Flags32 flags, ReturnStream &rets ) const
{
#ifdef __S2_DS_ENABLE_PARAMS
    rets.BeginComplex("params", eType_Property_Object );
       Params::s_LeafDSH_FEM_Solid2D_Params_ArchetypeLibrary.ExportInstance( "Archetype_s2_ds_LDSH_FEMS2D_Params", &m_Params, rets );
    rets.EndComplex();
#endif
}

void LeafDSH_FEM_Solid2D::SyncParams( const ParamIt &pit, ReturnStream &rets )
{
#ifdef __S2_DS_ENABLE_PARAMS
    util::ItemStream::ItemItRW out_pit = rets.WriteItem( "params", pit );
    Params::s_LeafDSH_FEM_Solid2D_Params_ArchetypeLibrary.SyncInstance( "Archetype_s2_ds_LDSH_FEMS2D_Params", &m_Params, out_pit.GetSubItem() );
#endif
}


#ifdef __ENABLE_CONTACT_D2K

bool LeafDSH_FEM_Solid2D::CheckError_Contact() const
{
    //\todo compute some sensible upper limit from object sizes
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
        for( unsigned int it_cpc=0; it_cpc<it_ccd2k->GetNumCPC(); it_cpc++ )
            if( it_ccd2k->GetCPC(it_cpc).m_Depth > m_Params.m_ContactSolver_MaxDepth )
                return false;
    return true;
}

void LeafDSH_FEM_Solid2D::NotifyPairBP( const geo::bp::Proxy *p_other, const geo::bp::Pair *p_pair )
{
    //DS_LOG_WARNING("LeafDSH_FEM_Solid2D::NotifyPairBP");
    IEntity *pOtherEntity( p_other->GetObj<IEntity*>() );
    DS_ASSERT( eEntity_Geom == pOtherEntity->GetEntityType() );
    IGeom *pGeom( static_cast<IGeom*>( pOtherEntity ) );

#define __ENABLE_PERSISTENT_CONTACTS
#ifdef __ENABLE_PERSISTENT_CONTACTS
    switch( p_pair->m_State )
    {
    case geo::bp::Pair::eNew:
        {
            //DS_LOG_WARNING("New BP");
            ContactConstraintD2K *pCCD2K( m_poolCCD2K.New() );
            pCCD2K->Create();
            p_pair->SetUserData( pCCD2K );
            if( geo::mp::TestContact( static_cast<const geo::IObject2*>(base_type::GetGO()),
                                      static_cast<const geo::IObject2*>(pGeom->GetGO()),
                                      pCCD2K->m_TmpCD, &pCCD2K->m_CC ) ) //\todo Avoid casts specializing GetGO() with dimension
                pCCD2K->Reset( pCCD2K->m_TmpCD, 0.001f );
        }
        break;
    case geo::bp::Pair::ePersistent:
        {
            //DS_LOG_WARNING("Persistent BP");
            ContactConstraintD2K *pCCD2K( p_pair->GetUserData<ContactConstraintD2K*>() );
            if( geo::mp::TestContact( static_cast<const geo::IObject2*>(base_type::GetGO()),
                                      static_cast<const geo::IObject2*>(pGeom->GetGO()),
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
            //geo::np::ContactData2 cd;
            geo::np::ContactCache2 *p_cc(0);
            // \todo an IGeom2 should return a geo::IObject2, but an IGeom should return an geo::IObject, consider 2 different GetGO() methods...
            //if( geo::mp::TestContact( base_type::GetGO(), pGeom->GetGO(), m_TmpCD, p_cc ) )
            if( geo::mp::TestContact( static_cast<const geo::IObject2*>(base_type::GetGO()),
                                      static_cast<const geo::IObject2*>(pGeom->GetGO()),
                                      m_TmpCD, p_cc ) ) //TEMP: casts are ugly
            {
                //DS_LOG_WARNING("Contact_MeshSolid2_Plane2 => %d", m_TmpCD.Size() );
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
#endif //__ENABLE_CONTACT_D2K

Real LeafDSH_FEM_Solid2D::Compute_Effective_PoissonRatio( Real det_F ) const
{
    if( det_F >= 0 || m_Params.m_InvertedCompressibilityFactorDetF < 0.001 ) return m_Params.m_PoissonRatio;
    else //det_F < 0
    {
        //Real lambda01( mal::Clamp01(1.0 + m_InvertedCompressibilityFactorDetF*det_F) ); //\todo old, not exactly what I wanted...
        Real lambda01( mal::Clamp01( Real(1) + det_F/m_Params.m_InvertedCompressibilityFactorDetF ) ); //l01 = 1 if det_F >= 0, l01 = 0 if det_F < -ICFDF
        // Linear ipol
        //return m_Params.m_PoissonRatio * lambda01; // + (1-lambda01)*0
        // Cubic ipol
        return m_Params.m_PoissonRatio * mal::Sq(lambda01) * (3-2*lambda01); //h01 polynomial from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
    }
}

Mat6x6 LeafDSH_FEM_Solid2D::Compute_Effective_Ke( unsigned int index_e ) const
{
    const ms::fem::LinearTriangle2 &element( m_vecElements[index_e] );
    Real det_F( m_vecDetF[index_e] );
    if( det_F >= 0 ) return element.K();
    else //det_F < 0
    {
        if( m_Params.m_InvertedCompressibilityFactorDetF > 0.001 )
        {
#define __USE_PROPER_EFFECTIVE_Ke
#ifdef __USE_PROPER_EFFECTIVE_Ke
            Real poisson_ratio( Compute_Effective_PoissonRatio(det_F) );
            Real a( m_Params.m_YoungModulus / ((1+poisson_ratio)*(1-2*poisson_ratio)) );
            Real a_times_nu( a * poisson_ratio );
            Mat3x3 E;
            E(0,0) = a - a_times_nu; E(0,1) = a_times_nu;     E(0,2) = Real(0);
            E(1,0) = a_times_nu;     E(1,1) = a - a_times_nu; E(1,2) = Real(0);
            E(2,0) = 0;              E(2,1) = 0;              E(2,2) = Real(0.5)*a - a_times_nu;
            Mat3x6 B;
            element.Compute_B( element.r(0), element.r(1), element.r(2), B );
            return element.Area() * ( B.Transposed() * (E * B) );
#else
            // Compressible Ke => Increased young, Zero poisson
            Real compressible_young_factor( 1/*+2*m_Params.m_PoissonRatio*/ ); //increase multiplier for steeper force increase wrt poisson_ratio
            Real a( compressible_young_factor * m_Params.m_YoungModulus );
            // Comute B from node positions
            Mat3x6 B;
            element.Compute_B( element.r(0), element.r(1), element.r(2), B );
            /*Setting young_modulus to 0 in E yields
            Mat3x3 E;
            E(0,0) = a; E(0,1) = 0; E(0,2) = Real(0);
            E(1,0) = 0; E(1,1) = a; E(1,2) = Real(0);
            E(2,0) = 0; E(2,1) = 0; E(2,2) = 0.5*a;
            Mat6x6 KeC( element.Area() * ( B.Transposed() * (E * B) ) );
            */
            //Assuming E(2,2) == a instead of 0.5*a actually yields better force directions, less inverted-incompressible (=> E = a*Id3)
            Mat6x6 KeC( (element.Area() * a) * (B.Transposed() * B) );
            // Linear interpolation wrt det_F is C0 by now, could be made C1 or C2 with quadratic or cubic interpolation
            Real lambda01( mal::Clamp01( Real(1) + m_Params.m_InvertedCompressibilityFactorDetF*det_F ) );
            return lambda01*element.K() + (Real(1)-lambda01)*KeC; //This works, but may be SLOW to compute in RT if KeC is not precomputed, and memory-inefficient if it is...
#endif
        }
        else
        {
            //\todo Implicitly, det_F multiplies both young and poisson coeffs... so any k_det_F should be a dimensionless fraction 0..1, NOT a material param like lame lambda|mu
            return (Real(1)-m_Params.m_FactorDetF*det_F)*element.K();
        }
    }
}

/*! Compute deformations F and their determinants in a single pass
 */
void LeafDSH_FEM_Solid2D::Compute_F( Mat2x2 *vec_F, Real *vec_det_F ) const
{
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
        element.ComputeF( m_Model.GetPos( element.nid(0) ),
                          m_Model.GetPos( element.nid(1) ),
                          m_Model.GetPos( element.nid(2) ),
                          vec_F[it_e] );
        vec_det_F[it_e] = mal::Det( vec_F[it_e] );
    }
}

/*! Find Non-Degenerate-Neighbour of element index_e
  \todo CONSIDER if neighbours with 0 << det(F) < m_DegenerateThresholdDetF should be considered Non-Degenerate here...
*/
int LeafDSH_FEM_Solid2D::FindNDN( unsigned int index_e ) const
{
    const ms::fem::LinearTriangle2 &element( m_vecElements[index_e] );
    int vec_ndn[3];
    int num_ndn(0);
    geo::feature_index_type vec_neighbour[3];
    unsigned int it_he( m_pMeshS->P_FirstHEID(index_e) );
    vec_neighbour[0] = m_pMeshS->HE_RightPID(it_he); it_he = m_pMeshS->HE_Next(it_he);
    vec_neighbour[1] = m_pMeshS->HE_RightPID(it_he); it_he = m_pMeshS->HE_Next(it_he);
    vec_neighbour[2] = m_pMeshS->HE_RightPID(it_he); it_he = m_pMeshS->HE_Next(it_he);
    for( int it_neighbour=0; it_neighbour<3; it_neighbour++ )
        if( vec_neighbour[it_neighbour] != geo::cInvalidFeatureIndex )
            if( m_vecDetF[ vec_neighbour[it_neighbour] ] > m_Params.m_DegenerateThresholdDetF )
                vec_ndn[num_ndn++] = it_neighbour;
    if( num_ndn > 0 )
    {
        //\todo Consider selecting best...
        return vec_ndn[0];
    }
    else
        return geo::cInvalidFeatureIndex;
}

int LeafDSH_FEM_Solid2D::FindNoC( unsigned int index_e ) const
{
    switch( m_Params.m_NCM )
    {
    case Params::eNCM_NDN:
        {
            int ndn = FindNDN(index_e);
            if( ndn != geo::cInvalidFeatureIndex ) return (ndn + 2) % 3;
            else return geo::cInvalidFeatureIndex;
        }
        break;
    case Params::eNCM_ToC:
    case Params::eNCM_IHFSDM:
        {
            for( PoolDED::iterator it_ded = m_poolDED.Begin(); it_ded.IsValid(); ++it_ded )
                if( it_ded->m_ElementId == index_e )
                    return it_ded->m_NoC;
            return geo::cInvalidFeatureIndex;
        }
        break;
    default: return geo::cInvalidFeatureIndex; break;
    }
}


void LeafDSH_FEM_Solid2D::UpdateDED()
{
#ifdef __ENABLE_CONTINUOUS
    // Update previously degenerated elements
    PoolDED::iterator it_ded( m_poolDED.Begin() );
    while( it_ded.IsValid() )
    {
        if( m_vecDetF[ it_ded->m_ElementId ] > m_Params.m_DegenerateThresholdDetF )
        {
            //DS_LOG_WARNING( "Fixed Degenerated %d age %d", it_ded->m_ElementId, it_ded->m_Age );
            DegeneratedElementData *p_deleted( it_ded.GetPtr() );
            ++it_ded;
            m_poolDED.Delete( p_deleted );
        }
        else
        {
            //DS_LOG_WARNING( "Persistent Degenerated %d age %d", it_ded->m_ElementId, it_ded->m_Age );
            it_ded->m_Age++;
            ++it_ded;
        }
    }
    // Add new degenerated elements
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        if( m_vecDetF[it_e] <= m_Params.m_DegenerateThresholdDetF )
        {
            const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
            Mat2x2 F0;
            element.ComputeF( m_vecPos0[element.nid(0)], m_vecPos0[element.nid(1)], m_vecPos0[element.nid(2)], F0 );
            Real det_F0( mal::Det(F0) );
            if( det_F0 > m_Params.m_DegenerateThresholdDetF ) //if not degenerated before, add \note it is VERY IMPORTANT to use EXACTLY the SAME criteria here that when the DED was added
            {
                // Gather node positions and displacements
                Vec2 p0( m_vecPos0[ element.nid(0) ] );
                Vec2 p1( m_vecPos0[ element.nid(1) ] );
                Vec2 p2( m_vecPos0[ element.nid(2) ] );
                Vec2 d0( m_Model.GetPos( element.nid(0) ) - p0 );
                Vec2 d1( m_Model.GetPos( element.nid(1) ) - p1 );
                Vec2 d2( m_Model.GetPos( element.nid(2) ) - p2 );
                // Compute ToC, when det(F) = m_DegenerateThresholdDetF => det(Q) = m_DegenerateThresholdDetF * det(P) => det(Q) = m_DegenerateThresholdDetF * 2 * Area(P)
                Real D( m_Params.m_DegenerateThresholdDetF * 2 * element.Area() );
                Real a(   d1.x() * d2.y()
                        + d2.x() * d0.y()
                        + d0.x() * d1.y()
                        - d1.x() * d0.y()
                        - d0.x() * d2.y()
                        - d2.x() * d1.y() );
                Real b( (   d1.x() * p2.y() + p1.x() * d2.y() )
                        + ( d2.x() * p0.y() + p2.x() * d0.y() )
                        + ( d0.x() * p1.y() + p0.x() * d1.y() )
                        - ( d1.x() * p0.y() + p1.x() * d0.y() )
                        - ( d0.x() * p2.y() + p0.x() * d2.y() )
                        - ( d2.x() * p1.y() + p2.x() * d1.y() ) );
                Real c(   p1.x() * p2.y()
                        + p2.x() * p0.y()
                        + p0.x() * p1.y()
                        - p1.x() * p0.y()
                        - p0.x() * p2.y()
                        - p2.x() * p1.y()
                        - D );
                Real toc1(0), toc2(1);
                int num_roots( mal::GSolvePolynomialEq2<Real>( a, b, c, toc1, toc2 ) );
                if( num_roots > 0 )
                {
                    //DS_LOG_WARNING( "Potentially degenerated %d with toc = %f, %f", it_e, toc1, toc2 );
                    DS_ASSERT( toc1 <= toc2 );
                    //\todo We knot that det(F(0)) > DTDF, and that det(F(1)) <= DTDF, so there MUST BE STRICTLY 1 crossing in toc = [0,1]
                    Real toc( toc1 > 0 ? toc1 : toc2 );
                    if( toc > 1 )
                    {
                        DS_LOG_ERROR( "toc %f > 1", toc );
                    }
                    /* analyze inversion interval:
                       0 <= toc1 <= toc2 < 1 => Collapse->Invert->Uninvert->Uncollapse
                       toc1 < 0 <= toc2 < 1 => Collapse->Invert
                       toc1 <= toc2 < 0 => Extrapolated collapse outside [0,1], no actual collapse
                       1 <= toc1 <= toc2 < 1 => Extrapolated collapse outside [0,1], no actual collapse
                    */
                    /*TEMP
                    if( !(toc2 < 0 || toc1 > 1) //Inside [0,1]
                        &&
                        (toc1 < 0 || toc2 > 1) //2 roots inside [0,1] //\todo MAY UNINVERT ACROSS AXIS DIFERENT than collapsed one!!
                        )
                    */
                    {
                        // Compute NoC
                        uint32 noc( geo::cInvalidFeatureIndex );
                        switch( m_Params.m_NCM )
                        {
                        case Params::eNCM_IHFSDM:
                            {
                                //\todo ACTUALLY, COULD PERFORM ANALYSIS AT COLLAPSE CONFIGURATION INSTEAD OF CURRENT ONE!!!!!! WOULD BE MORE ACCURATE
                                // Compute current node positions and F
                                Vec2 vec_pos[3];
                                vec_pos[0] = m_Model.GetPos( element.nid(0) );
                                vec_pos[1] = m_Model.GetPos( element.nid(1) );
                                vec_pos[2] = m_Model.GetPos( element.nid(2) );
                                Mat2x2 F;
                                element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], F );
                                // Compute SVD
                                Mat2x2 U, Vt;
                                Vec2 diag_F;
                                mal::GSingularValueDecomposition_USVt( F, U, diag_F, Vt );
                                // Find *geometrically* shortest inversion direction in reference configuration
                                //\todo NOTICE that the selected inversion direction changes with dist12
                                Vec2 n_c[3];
                                n_c[0] = mal::PerpendicularCW( vec_pos[2] - vec_pos[1] );
                                n_c[1] = mal::PerpendicularCW( vec_pos[0] - vec_pos[2] );
                                n_c[2] = mal::PerpendicularCW( vec_pos[1] - vec_pos[0] );
                                Vec2 v[2];
                                v[0] = mal::GRow<0>(Vt);
                                v[1] = mal::GRow<1>(Vt);
                                Real min_lambda( mal::Infinity<Real>() );
                                int min_it_v(0);
                                //Real lambda_CxV[3][2];
                                for( int it_c=0; it_c<3; it_c++ )
                                {
                                    for( int it_v=0; it_v<2; it_v++ )
                                    {
                                        Real dot_c_v = mal::Dot( v[it_v], n_c[it_c] );
                                        Real lambda_c_v = mal::Abs(dot_c_v) > 0.0001f
                                                          ? mal::Dot( vec_pos[(it_c+1)%3] - vec_pos[it_c], n_c[it_c] ) / dot_c_v
                                                          : mal::Infinity<Real>();
                                        lambda_c_v = mal::Abs(lambda_c_v);
                                        //lambda_CxV[it_c][it_v] = lambda_c_v;
                                        if( lambda_c_v < min_lambda )
                                        {
                                            min_lambda = lambda_c_v;
                                            min_it_v = it_v;
                                            noc = it_c;
                                        }
                                    }
                                }
                            }
                            break;
                        case Params::eNCM_ToC:
                        default: //\note Default => ToC, used ony for Viz if m_NCM != eNCM_ToC
                            {
                                Vec2 q0( p0 + toc*d0 );
                                Vec2 q1( p1 + toc*d1 );
                                Vec2 q2( p2 + toc*d2 );
                                Real sql01( mal::NormSq( q1-q0 ) );
                                Real sql12( mal::NormSq( q2-q1 ) );
                                Real sql20( mal::NormSq( q0-q2 ) );
                                // Choose NoC, vertex opposite to longest edge at ToC
                                noc = (sql01 >= sql12)
                                      ? (sql01 >= sql20) ? 2 : 1
                                      : (sql12 >= sql20) ? 0 : 1;
                                /*TEMP: Debug
                                std::cout << "From "
                                          << p0 << ","
                                          << p1 << ","
                                          << p2 << " to "
                                          << m_Model.GetPos( element.nid(0) ) << ","
                                          << m_Model.GetPos( element.nid(1) ) << ","
                                          << m_Model.GetPos( element.nid(2) )
                                          << " ToC " << toc1 << "," << toc2
                                          << " NoC " << noc
                                          << " PoC " << ( (sql01 >= sql12) ? (sql01 >= sql20) ? q2 : q1 : (sql12 >= sql20) ? q0 : q1 )
                                          << " Det(Q(ToC)) " << mal::Det( mal::GMat2x2_From_Columns(q1-q0,q2-q0) ) << " ?= D = " << D << std::endl;
                                */
                            }
                            break;
                        }
                        DegeneratedElementData *pDED = m_poolDED.New();
                        pDED->Init( it_e, noc, toc1 );
                        //DS_LOG_WARNING( "Actually degenerated %d with noc = %d, toc = %f", it_e, noc, toc );
                    }
                    /*TEMP
                    else
                    {
                        DS_LOG_WARNING( "No 0 or 2 degenerations in range [0,1]" );
                    }
                    */
                }
                else
                {
                    DS_LOG_ERROR( "New Degenerated %d CANNOT SOLVE Eq2", it_e );
                }
            }
        }
    }
#endif
}

//---- Corotational
void LeafDSH_FEM_Solid2D::Compute_T( Transform2 *vec_te2r ) const
{
    unsigned int num_degenerate(0);
    float32 sum_degenerate_detF(0);
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
        // Gather node positions
        Vec2 vec_pos[3];
        vec_pos[0] = m_Model.GetPos( element.nid(0) );
        vec_pos[1] = m_Model.GetPos( element.nid(1) );
        vec_pos[2] = m_Model.GetPos( element.nid(2) );
        // Compute Te transform {e}->{0}
        Vec2 Be( mal::Rcp(Real(3))*(vec_pos[0]+vec_pos[1]+vec_pos[2]) );
        vec_te2r[it_e].m_Pos = Be - element.Barycenter0();
        DS_ASSERT( !mal::IsNaN( vec_te2r[it_e].m_Pos ) );

        Mat2x2 F( m_vecF[it_e] );
        Real det_F( m_vecDetF[it_e] );

        // Save degeneration state
        if( det_F <= m_Params.m_DegenerateThresholdDetF )
        {
            num_degenerate++;
            sum_degenerate_detF += det_F - m_Params.m_DegenerateThresholdDetF; //distance-to-undegenerate
        }

        // Compute T.R according to EMethod, fixing degeneration if necessary
        Mat2x2 R;//( Compute_R( x0, x1, x2, F, det_F ) )

        switch( m_Params.m_Method )
        {
        case Params::eMethod_Id:
            R = Mat2x2::Identity();
            break;
        case Params::eMethod_QR:
            {
                R = mal::GRotation2x2_GramSchmidtOrthonormalization_YX( F, Mat2x2::Identity() ); //\note YX to show the glitch with "aligned_romboid" scenes
                DS_ASSERT( !mal::IsNaN( R ) );
            }
            break;
        case Params::eMethod_PD:
            {
                if( mal::Abs(det_F) >= m_Params.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
                else R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, det_F ); //\note Handles collapse but NOT inversion
                DS_ASSERT( !mal::IsNaN( R ) );
            }
            break;
            // InvertedAxis cases
        case Params::eMethod_PD_Reflect:
        case Params::eMethod_PD_Fix:
        case Params::eMethod_PD_Project:
        case Params::eMethod_IHFSDM:
            {
                if( det_F > m_Params.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
                else
                {
                    //\todo Searching for NoC is SLOW, specially if eNCM_NoC is enabled, consider a specific pass for DED, ignoring them in an "undegenerated" first pass.
                    int noc = FindNoC( it_e );
                    if( noc != geo::cInvalidFeatureIndex )
                    {
                        //DS_LOG_WARNING("Using NoC %d", noc );
                        Vec2 &x0( vec_pos[ noc ] ); //Inverted vtx
                        const Vec2 &x1( vec_pos[ (noc + 1) % 3 ] );
                        const Vec2 &x2( vec_pos[ (noc + 2) % 3 ] );
                        Vec2 axis12( mal::PerpendicularCW( x2 - x1 ));
                        Real dist12( mal::Norm(axis12) );
                        if( dist12 > 0.001f )  //\todo use sensible epsilon
                        {
                            axis12 /= dist12;

                            // Use axis according to Method
                            switch( m_Params.m_Method )
                            {
                            case Params::eMethod_PD_Reflect:
                                {
                                    if( det_F <= -m_Params.m_DegenerateThresholdDetF ) //Inverted & !Collapsed
                                    {
                                        R = mal::GRotation2x2_PolarDecomposition( mal::GReflection2x2_From_Axis( axis12 ) * F, -det_F ); //Safe uncollapsed simple PD
                                        DS_ASSERT( !mal::IsNaN( R ) );
                                    }
                                    else if( det_F < 0 ) //-m_DegenerateThresholdDetF < det_F < 0 => Inverted & Collapsed
                                    {
                                        R = mal::GRotation2x2_PolarDecomposition_From_SVD( mal::GReflection2x2_From_Axis( axis12 ) * F, -det_F ); //\todo SHOULD Handle det_F = 0 properly
                                        DS_ASSERT( !mal::IsNaN( R ) );
                                    }
                                    else // 0 <= det_F < m_DegenerateThresholdDetF => !Inverted && Collapsed
                                    {
                                        R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, det_F ); //\todo SHOULD Handle det_F = 0 properly, but yields NaN sometimes
                                        DS_ASSERT( !mal::IsNaN( R ) );
                                    }
                                }
                                break;
                            case Params::eMethod_PD_Fix:
                                {
                                    Vec2 x12( x2 - x1 );
                                    Real dist12_sq( mal::NormSq(x12) );
                                    if( dist12_sq > 0.000001f )
                                    {
                                        // Compute rotation from r12 to x12, which is the same as global R when we rotate the resting state R=Id rigidly with x12
                                        Vec2 r1( element.r( (noc + 1) % 3 ) );
                                        Vec2 r2( element.r( (noc + 2) % 3 ) );
                                        Vec2 r12( r2-r1 );
                                        //\todo In 3D, rotR would be similarly computed from the NORMAL of the crossed face
                                        Mat2x2 rotR( mal::GRotation2x2_Vec2Vec( mal::Normalized(r12) , mal::Normalized(x12) ) ); //\todo Normalized(r12) could be precomputed
//#define __ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
#ifdef __ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
                                        if( det_F > 0 )
                                        {
                                            Vec2 n12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                                            /* Compute delta so that projected x0 yields an element with det(F) = DTDF
                                               det(F) = det(Q) / det(P) = DTDF
                                               => det(Q) = DTDF * det(P)
                                               ( det(Q) = 2*Area(Q) = delta * dist(x1,x2) )
                                               => delta * dist(x1,x2) = DTDF * det(P)
                                               => delta = DTDF * det(P) / dist(x1,x2)
                                               => delta = DTDF * 2*Area(P) / dist(x1,x2)
                                            */
                                            // Faster code that AVOIDS sqrt (using unnormalized axis12 and dist12_sq)
                                            Real delta_div_dist12( m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist12_sq ); //== delta / dist12
                                            x0 -= ( mal::Dot(x0-x1,n12)/dist12_sq ) * n12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
                                            x0 += delta_div_dist12*n12; //Move past collapse

                                            // Compute Proj(F)
                                            Mat2x2 projF;
                                            element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], projF );
                                            Real det_projF( mal::Det(projF) );
                                            DS_ASSERT( det_projF >= 0 );
                                            if( det_projF < 0.99f*m_Params.m_DegenerateThresholdDetF )
                                                DS_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_Params.m_DegenerateThresholdDetF );
                                            Mat2x2 projR = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                                            //R = RLerp( rotR, projR, det_F / m_Params.m_DegenerateThresholdDetF );
                                            R = RCerp( rotR, projR, det_F / m_Params.m_DegenerateThresholdDetF );
                                        }
                                        else //det_F <= 0
                                            R = rotR;
#else //__ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
                                        if( det_F > m_Params.m_ThresholdIpolDetF )
                                        {
                                            Vec2 n12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                                            /* Compute delta so that projected x0 yields an element with det(F) = DTDF
                                               det(F) = det(Q) / det(P) = DTDF
                                               => det(Q) = DTDF * det(P)
                                               ( det(Q) = 2*Area(Q) = delta * dist(x1,x2) )
                                               => delta * dist(x1,x2) = DTDF * det(P)
                                               => delta = DTDF * det(P) / dist(x1,x2)
                                               => delta = DTDF * 2*Area(P) / dist(x1,x2)
                                            */
                                            // Faster code that AVOIDS sqrt (using unnormalized axis12 and dist12_sq)
                                            Real delta_div_dist12( m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist12_sq ); //== delta / dist12
                                            x0 -= ( mal::Dot(x0-x1,n12)/dist12_sq ) * n12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
                                            x0 += delta_div_dist12*n12; //Move past collapse

                                            // Compute Proj(F)
                                            Mat2x2 projF;
                                            element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], projF );
                                            Real det_projF( mal::Det(projF) );
                                            DS_ASSERT( det_projF >= 0 );
                                            if( det_projF < 0.99f*m_Params.m_DegenerateThresholdDetF )
                                                DS_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_Params.m_DegenerateThresholdDetF );
                                            Mat2x2 projR = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                                            Real lambda01( (det_F - m_Params.m_ThresholdIpolDetF) / (m_Params.m_DegenerateThresholdDetF - m_Params.m_ThresholdIpolDetF) );
                                            R = RCerp( rotR, projR, lambda01 );
                                        }
                                        else //det_F <= m_Params.m_ThresholdIpolDetF
                                            R = rotR;
#endif //__ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
                                    }
                                    else
                                        R = Mat2x2::Identity();
                                    DS_ASSERT( !mal::IsNaN( R ) );
                                }
                                break;
                            case Params::eMethod_PD_Project:
                                {
                                    //\pre det_F <= m_DegenerateThresholdDetF , therefore Inverted | Collapsed
                                    /* Compute delta so that projected x0 yields an element with det(F) = DTDF
                                       det(F) = det(Q) / det(P) = DTDF
                                       => det(Q) = DTDF * det(P)
                                       ( det(Q) = 2*Area(Q) = delta * dist(x1,x2) )
                                       => delta * dist(x1,x2) = DTDF * det(P)
                                       => delta = DTDF * det(P) / dist(x1,x2)
                                       => delta = DTDF * 2*Area(P) / dist(x1,x2)
                                    */
                                    Real delta( m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist12 );
                                    x0 -= mal::Dot(x0-x1,axis12) * axis12; //Project on-to axis
                                    x0 += delta * axis12; //Move past collapse area, to ensure numerically safe
                                    // Recompute H and R from it
                                    Mat2x2 projF;
                                    element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], projF );
                                    Real det_projF = mal::Det(projF);
                                    DS_ASSERT( det_projF >= 0 );
                                    if( !mal::ApproxEq(det_projF,m_Params.m_DegenerateThresholdDetF,0.001f) )
                                    {
                                        DS_LOG_WARNING( "det(F) %f != %f, delta %f", det_projF, m_Params.m_DegenerateThresholdDetF, delta );
                                    }
                                    if( det_projF < 0.99f*m_Params.m_DegenerateThresholdDetF )
                                        DS_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_Params.m_DegenerateThresholdDetF );
                                    R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                                    DS_ASSERT( !mal::IsNaN( R ) );
                                }
                                break;
                                //\todo
                            case Params::eMethod_IHFSDM:
                                {
                                    // Compute SVD
                                    Mat2x2 U, Vt;
                                    Vec2 diag_S;
                                    mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
//#define __ENABLE_HACKISH_IHFSDM
#ifdef __ENABLE_HACKISH_IHFSDM
                                    //\note This yields a reasonable force field, but I'm not sure why
                                    if( mal::Abs( mal::Dot( mal::GColumn<0>(U), axis12 ) ) < mal::Abs( mal::Dot( mal::GColumn<1>(U), axis12 ) ) )
                                    {
                                        if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                                        if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                                    }
                                    else
                                    {
                                        if( mal::Det(U) < 0 ) mal::GSetColumn<0>( U, -mal::GColumn<0>(U) );
                                        if( mal::Det(Vt) < 0 ) mal::GSetRow<0>( Vt, -mal::GRow<0>(Vt) );
                                    }
#else
                                    // Find *geometrically* shortest inversion direction in reference configuration
                                    //\todo NOTICE that the selected inversion direction changes with dist12
                                    Vec2 n_c( mal::PerpendicularCW( element.r((noc+2)%3) - element.r((noc+1)%3) ) );
                                    Vec2 v[2];
                                    v[0] = mal::GRow<0>(Vt);
                                    v[1] = mal::GRow<1>(Vt);
                                    Real min_lambda( mal::Infinity<Real>() );
                                    int min_it_v(0);
                                    for( int it_v=0; it_v<2; it_v++ )
                                    {
                                        Real dot_c_v = mal::Dot( v[it_v], n_c );
                                        Real lambda_c_v = mal::Abs(dot_c_v) > 0.0001f
                                                          ? mal::Dot( element.r((noc+1)%3) - element.r(noc), n_c ) / dot_c_v
                                                          : mal::Infinity<Real>();
                                        lambda_c_v = mal::Abs(lambda_c_v);
                                        if( lambda_c_v < min_lambda )
                                        {
                                            min_lambda = lambda_c_v;
                                            min_it_v = it_v;
                                        }
                                    }
                                    // Invert selected direction
                                    if( 0 == min_it_v )
                                    {
                                        if( mal::Det(U) < 0 ) mal::GSetColumn<0>( U, -mal::GColumn<0>(U) );
                                        if( mal::Det(Vt) < 0 ) mal::GSetRow<0>( Vt, -mal::GRow<0>(Vt) );
                                    }
                                    else
                                    {
                                        if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                                        if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                                    }
#endif
                                    R = U*Vt;
                                    DS_ASSERT( !mal::IsNaN( R ) );
                                }
                            default: break;
                            }
                        }
                        else
                        {
                            DS_LOG_ERROR( "dist12 = %f is TOO SMALL, NDN is collapsed? this should NEVER HAPPEN unless element is completely collapsed, using R = Id", dist12 );
                            R = Mat2x2::Identity();
                        }
                    }
                    else
                    {
                        //\todo CANNOT USE Propagation method, consider using Local method (QR or SVD based)
                        //
                        R = Mat2x2::Identity();
                        //
                        DS_LOG_ERROR("eMethod_XXX: NoC not found, all neighbours collapsed? Using R = Id");
                        /*
                        //TEMP: Try Local SVD... but this may induce singularities/stable inversions
                        Mat2x2 U, Vt;
                        Vec2 diag_S;
                        mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
                        // Force pure rotation U,Vt by fixing potential inversion/reflection
                        if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                        if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                        R = U*Vt;
                        DS_ASSERT( !mal::IsNaN( R ) );
                        */
                    }
                }
                DS_ASSERT( !mal::IsNaN( R ) );
            }
            break;
        case Params::eMethod_PD_Project_Nearest:
            {
                if( det_F > m_Params.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
                else
                {
                    // Compute all possible axis
                    Vec2 axis12( mal::PerpendicularCW( vec_pos[2] - vec_pos[1] ) );
                    Real dist12( mal::Norm( axis12 ) );
                    bool bCollapsed12( dist12 < 0.0001f );

                    Vec2 axis20( mal::PerpendicularCW( vec_pos[0] - vec_pos[2] ) );
                    Real dist20( mal::Norm( axis20 ) );
                    bool bCollapsed20( dist20 < 0.0001f );

                    Vec2 axis01( mal::PerpendicularCW( vec_pos[1] - vec_pos[0] ) );
                    Real dist01( mal::Norm( axis01 ) );
                    bool bCollapsed01( dist01 < 0.0001f );
                    if( !bCollapsed12 || !bCollapsed20 || !bCollapsed01 )
                    {
                        // Compute all possible projections, discarding potentially collapsed axis
                        Real delta12(0);
                        if( !bCollapsed12 )
                        {
                            axis12 /= dist12;
                            delta12 = m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist12;
                        }
                        Vec2 projected_x0( vec_pos[0] - mal::Dot(vec_pos[0]-vec_pos[1],axis12)*axis12 ); //Project onto collapse line
                        projected_x0 += delta12*axis12; //Move past collapse

                        Real delta20(0);
                        if( !bCollapsed20 )
                        {
                            axis20 /= dist20;
                            delta20 = m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist20;
                        }
                        Vec2 projected_x1( vec_pos[1] - mal::Dot(vec_pos[1]-vec_pos[0],axis20)*axis20 ); //Project onto collapse line
                        projected_x1 += delta20*axis20; //Move past collapse

                        Real delta01(0);
                        if( !bCollapsed01 )
                        {
                            axis01 /= dist01;
                            delta01 = m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist01;
                        }
                        Vec2 projected_x2( vec_pos[2] - mal::Dot(vec_pos[2]-vec_pos[0],axis01)*axis01 ); //Project onto collapse line
                        projected_x2 += delta01*axis01; //Move past collapse

                        // Compute projection distances, to pick nearest valid (uncollapsed axis) one
                        Real d0( mal::NormSq( vec_pos[0] - projected_x0 ) );
                        Real d1( mal::NormSq( vec_pos[1] - projected_x1 ) );
                        Real d2( mal::NormSq( vec_pos[2] - projected_x2 ) );
                        Real max_d( d0 + d1 + d2 );
                        if( bCollapsed12 ) d0 = max_d;
                        if( bCollapsed20 ) d1 = max_d;
                        if( bCollapsed01 ) d2 = max_d;

                        // Choose nearest projection axis and compute projF and its R accordingly
                        Vec2 x0( vec_pos[0] );
                        Vec2 x1( vec_pos[1] );
                        Vec2 x2( vec_pos[2] );
                        if( d0 <= d1 )
                        {
                            if( d0 <= d2 ) x0 = projected_x0;
                            else x2 = projected_x2;
                        }
                        else
                        {
                            if( d1 <= d2 ) x1 = projected_x1;
                            else x2 = projected_x2;
                        }
                        Mat2x2 projF;
                        element.ComputeF( x0, x1, x2, projF );
                        Real det_projF( mal::Det(projF) );
                        DS_ASSERT( det_projF >= 0 );
                        if( det_projF < 0.99f*m_Params.m_DegenerateThresholdDetF )
                            DS_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_Params.m_DegenerateThresholdDetF );
                        R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                    }
                    else
                    {
                        R = Mat2x2::Identity();
                        DS_LOG_ERROR("eMethod_PD_Project_Nearest: all edges collapsed? Using R = Id");
                        //\todo R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, det_F ); \note det_F can be small, nonzero BUT NEGATIVE, so if we decide to try to extract R, we must account for -eps < det_F < eps
                    }
                }
                DS_ASSERT( !mal::IsNaN( R ) );
            }
            break;
        case Params::eMethod_PD_SVD:
            {
                //\note Undegenerated includes det_F == 1, which FAILS for GSingularValueDecomposition_USVt() due to repeated eigenvalues, resulting in NaN
                if( det_F > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
                {
                    R = mal::GRotation2x2_PolarDecomposition( F, det_F );
                    DS_ASSERT( !mal::IsNaN( R ) );
                }
                else //Degenerate or Inverted
                {
                    // Compute SVD
                    Mat2x2 U, Vt;
                    Vec2 diag_S;
                    mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
                    // Force pure rotation U,Vt by fixing potential inversion/reflection
                    if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                    if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                    R = U*Vt;
                    DS_ASSERT( !mal::IsNaN( R ) );
                }
            }
            break;
            // \todo
        case Params::eMethod_ITF_LRM:
        case Params::eMethod_ECIE_CLRM:
        case Params::eMethod_ECIE_NHC0:
        case Params::eMethod_ECIE_NHC1:
        case Params::eMethod_PD_CLRM:
            R = Mat2x2::Identity();
            break;
        default:
            R = Mat2x2::Identity();
            break;
        }
        DS_ASSERT( !mal::IsNaN( R ) );
        vec_te2r[it_e].m_Rot = R;
    }
    S2_DS_STAT_SET( m_Stats.m_Num_Degenerate, num_degenerate );
    S2_DS_STAT_SET( m_Stats.m_Sum_Degenerate_DetF, sum_degenerate_detF );
}

//---- Assembly
void LeafDSH_FEM_Solid2D::Assemble_K( ns::MatrixD &K ) const
{
    // Reset
    K.Resize( 2*m_NumNodes, 2*m_NumNodes );
    K.Zero();
    // For each element i
    Mat6x6 Ke;
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );

        // Get element K^e
        Mat6x6 Re( Mat6x6::Zero() );
        for( int i=0; i<3; i++ )
        {
            Re(2*i,  2*i  ) = m_vecTe2r[it_e].m_Rot(0,0);
            Re(2*i,  2*i+1) = m_vecTe2r[it_e].m_Rot(0,1);
            Re(2*i+1,2*i  ) = m_vecTe2r[it_e].m_Rot(1,0);
            Re(2*i+1,2*i+1) = m_vecTe2r[it_e].m_Rot(1,1);
        }
        Ke = Re * Compute_Effective_Ke(it_e) * Re.Transposed();

#ifdef __ENABLE_TRACE_IE_VERBOSE
        std::cout << "Ke[" << it_e << "] = ( " << element.nid(0) << ", " << element.nid(1) << ", " << element.nid(2) << " )" << std::endl;
        std::cout << element.K() << std::endl;
#endif
        // Acc K^e_ij on global K
        for( int i=0; i<3; i++ )
            for( int j=0; j<3; j++ )
            {
                int vi( element.nid(i) );
                int vj( element.nid(j) );
                // Accumulate 2x2 submatrix for (vi,vj)
                K( 2*vi   , 2*vj   ) += Ke( 2*i   , 2*j   );
                K( 2*vi   , 2*vj+1 ) += Ke( 2*i   , 2*j+1 );
                K( 2*vi+1 , 2*vj   ) += Ke( 2*i+1 , 2*j   );
                K( 2*vi+1 , 2*vj+1 ) += Ke( 2*i+1 , 2*j+1 );
            }
    }
    /* TEMPORAL: Current implementation works by accumulating all Ke
       on global K, but it may be faster to write just 1/2 matrix and
       explicitly symmetrize it later or store it without redundancy.
    */

    // Assemble KPC
#ifdef __ENABLE_KNC
    if( m_Params.m_KinematicMode != Params::eKM_Disabled )
    {
        for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int i( it_knc->m_nid );
            for( unsigned int j=0; j<m_NumNodes; j++ )
            {
                K( 2*i  , 2*j   ) = 0;
                K( 2*i  , 2*j+1 ) = 0;
                K( 2*i+1, 2*j   ) = 0;
                K( 2*i+1, 2*j+1 ) = 0;
                /* Symm
                   K( 2*j  , 2*i   ) = 0;
                   K( 2*j+1, 2*i   ) = 0;
                   K( 2*j  , 2*i+1 ) = 0;
                   K( 2*j+1, 2*i+1 ) = 0;
                */
            }
            /* Diag NO, no forces generated on static vid whatever the displacements or velocities
               K( 2*i  , 2*i   ) = 1;
               K( 2*i+1, 2*i+1 ) = 1;
            */
        }
    }
#endif
}

void LeafDSH_FEM_Solid2D::Assemble_M( ns::MatrixD &M ) const
{
    // Reset
    M.Resize( 2*m_NumNodes, 2*m_NumNodes );
    M.Zero();
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Real mass( m_Model.GetMass(it_node) );
        M(2*it_node,2*it_node) = mass;
        M(2*it_node+1,2*it_node+1) = mass;
    }
}

void LeafDSH_FEM_Solid2D::Assemble_InvM( ns::MatrixD &InvM ) const
{
    // Reset
    InvM.Resize( 2*m_NumNodes, 2*m_NumNodes );
    InvM.Zero();
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Real inv_m( mal::Rcp(m_Model.GetMass(it_node)) );
        InvM(2*it_node,2*it_node) = inv_m;
        InvM(2*it_node+1,2*it_node+1) = inv_m;
    }
}

void LeafDSH_FEM_Solid2D::Assemble_v( ns::VectorD &v ) const
{
    v.Resize( 2*m_NumNodes );
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        v[ 2*it_node   ] = m_Model.GetVel(it_node)[0];
        v[ 2*it_node+1 ] = m_Model.GetVel(it_node)[1];
    }
}

void LeafDSH_FEM_Solid2D::Assemble_u( ns::VectorD &u ) const
{
    u.Resize( 2*m_NumNodes );
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Vec2 r0 = m_pMeshS->V_Pos_0( it_node );
        u[ 2*it_node   ] = m_Model.GetPos(it_node)[0] - r0[0];
        u[ 2*it_node+1 ] = m_Model.GetPos(it_node)[1] - r0[1];
    }
}

void LeafDSH_FEM_Solid2D::Assemble_x( ns::VectorD &x ) const
{
    x.Resize( 2*m_NumNodes );
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        x[ 2*it_node   ] = m_Model.GetPos(it_node)[0];
        x[ 2*it_node+1 ] = m_Model.GetPos(it_node)[1];
    }
}

void LeafDSH_FEM_Solid2D::Assemble_r( ns::VectorD &r ) const
{
    r.Resize( 2*m_NumNodes );
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Vec2 r0( m_pMeshS->V_Pos_0( it_node ) );
        r[ 2*it_node   ] = r0[0];
        r[ 2*it_node+1 ] = r0[1];
    }
}

void LeafDSH_FEM_Solid2D::Assemble_CorotationalElasticForce( ns::VectorD &cef ) const
{
    cef.Resize( 2*m_NumNodes );
    cef.Zero();
    // For each element e
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
    {
        const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
        // Build R^e (\todo could be avoided by multiplying while accumulating final components of cef_e UNROTATED)
        Mat6x6 Re( Mat6x6::Zero() );
        for( int i=0; i<3; i++ )
        {
            Re(2*i,  2*i  ) = m_vecTe2r[it_e].m_Rot(0,0);
            Re(2*i,  2*i+1) = m_vecTe2r[it_e].m_Rot(0,1);
            Re(2*i+1,2*i  ) = m_vecTe2r[it_e].m_Rot(1,0);
            Re(2*i+1,2*i+1) = m_vecTe2r[it_e].m_Rot(1,1);
        }
        Mat6x6 Ke( Compute_Effective_Ke(it_e) );
        /* NO, DO NOT Remove dependency on static nodes
        for( int i=0; i<3; i++ )
            if( IsStaticVID( element.nid(i) ) )
            {
                for( unsigned int j=0; j<3; j++ )
                {
                    Ke( 2*i  , 2*j   ) = 0;
                    Ke( 2*i  , 2*j+1 ) = 0;
                    Ke( 2*i+1, 2*j   ) = 0;
                    Ke( 2*i+1, 2*j+1 ) = 0;
                    //Symm MUST NOT be removed, because non-static nodes DO depend on static ones
                    //Ke( 2*j  , 2*i   ) = 0;
                    //Ke( 2*j+1, 2*i   ) = 0;
                    //Ke( 2*j  , 2*i+1 ) = 0;
                    //Ke( 2*j+1, 2*i+1 ) = 0;
                }
            }
        */
        // Assemble r^e
        Vec6 r_e;
        for( int i=0; i<3; i++ )
        {
            r_e[2*i  ] = element.r(i)[0];
            r_e[2*i+1] = element.r(i)[1];
        }
        // Assemble x^e in reference coords
        Vec6 x_e;
        for( int i=0; i<3; i++ )
        {
            Vec2 p( m_vecTe2r[it_e].Rot().Transposed() * m_Model.GetPos( element.nid(i) ) );
            x_e[2*i  ] = p[0];
            x_e[2*i+1] = p[1];
        }
        // Compute cef^e
        Vec6 cef_e( - Re * Ke * ( x_e - r_e ) );
        // Accumulate cef^e_ij on global cef
        for( int i=0; i<3; i++ )
        {
            int vi( element.nid(i) );
            //\todo HACK: I think IsStaticVID() exception should happen automatically due to De2r, but maybe not... maybe De2r SHOULD displace x_e even if static?!?! If needed, better add it as a post-step and iterate only on static VID, FASTER
#ifdef __ENABLE_KNC
            if( !IsKNC(vi) )
#endif
            {
                cef[2*vi  ] += cef_e[2*i  ];
                cef[2*vi+1] += cef_e[2*i+1];
            }
        }
    }
}

void LeafDSH_FEM_Solid2D::Assemble_Fe( ns::VectorD &Fe ) const
{
    Fe.Resize( 2*m_NumNodes );
    Fe.Zero();
    // Apply gravity
    for( uint32 it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Fe[ 2*it_node   ] = m_Model.GetMass(it_node) * m_Params.m_Gravity[0];
        Fe[ 2*it_node+1 ] = m_Model.GetMass(it_node) * m_Params.m_Gravity[1];
    }
    // Assemble KNC
#ifdef __ENABLE_KNC
    if( m_Params.m_KinematicMode != Params::eKM_Disabled )
    {
        for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int i( it_knc->m_nid );
            Fe[ 2*i   ] = 0;
            Fe[ 2*i+1 ] = 0;
        }
    }
#endif

#ifdef __ENABLE_CONTACT_D2K_PRE
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        //if( !ccd2k.IsEmpty() ) DS_LOG_WARNING( "Contact %llx", (machine_uint_type)&ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                switch( m_Params.m_ContactSolver_Type )
                {
                case Params::eCST_None: break;
                case Params::eCST_Penalty:
                    {
                        if( cpc.m_POF1.m_FeatureId.IsVertex() )
                        {
                            geo::feature_index_type vid1 = cpc.m_POF1.m_FeatureId.AsVertex();
                            Vec2 vel( m_Model.GetVel(vid1) );
                            Real pressure_n( m_Params.m_ContactSolver_Penalty_Ks * (cpc.m_Depth - m_Params.m_ContactSolver_DepthMargin)
                                             - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                            Real force_n( cpc.m_Radius * pressure_n );
                            Fe[ 2*vid1   ] += force_n * cpc.m_Normal[0];
                            Fe[ 2*vid1+1 ] += force_n * cpc.m_Normal[1];
                        }
                        else if( cpc.m_POF1.m_FeatureId.IsSegment() )
                        {
                            geo::feature_index_type eid = cpc.m_POF1.m_FeatureId.AsSegment();
                            geo::feature_index_type vid1 = m_pMeshS->HE_OriginVID( eid );
                            geo::feature_index_type vid2 = m_pMeshS->HE_FinalVID( eid );
                            // Interpolate V_cp from node velocities
                            Vec2 vel( cpc.m_POF1.m_BarycentricCoords[0] * m_Model.GetVel(vid1)
                                      + cpc.m_POF1.m_BarycentricCoords[1] * m_Model.GetVel(vid2) );
                            Real pressure_n( m_Params.m_ContactSolver_Penalty_Ks * (cpc.m_Depth - m_Params.m_ContactSolver_DepthMargin)
                                             - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                            Real force_n( cpc.m_Radius * pressure_n );
                            // Distribute F_cp accross nodes
                            Real force_n1( cpc.m_POF1.m_BarycentricCoords[0] * force_n );
                            Real force_n2( cpc.m_POF1.m_BarycentricCoords[1] * force_n );
                            Fe[ 2*vid1   ] += force_n1 * cpc.m_Normal[0];
                            Fe[ 2*vid1+1 ] += force_n1 * cpc.m_Normal[1];
                            Fe[ 2*vid2   ] += force_n2 * cpc.m_Normal[0];
                            Fe[ 2*vid2+1 ] += force_n2 * cpc.m_Normal[1];
                        }
                        else { DS_LOG_ASSERT( false, "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
                    }
                    break;
                case Params::eCST_Hack: break; //\note It only works as a POST
                default: break;
                }
            }
        }
    }
#endif //__ENABLE_CONTACT_D2K_PRE

}

void LeafDSH_FEM_Solid2D::Assemble_A( ns::MatrixD &A, Real dt, const ns::MatrixD &M, const ns::MatrixD &K )
{
    //A = M + dt*C + mal::Sq(dt)*K; C = r[0]*M + r[1]*K
    if( m_Params.m_RayleighCoeff[0] != Real(0) && m_Params.m_RayleighCoeff[1] != Real(0) )
        A = (Real(1) + dt*m_Params.m_RayleighCoeff[0]) * M + (mal::Sq(dt) + dt*m_Params.m_RayleighCoeff[1]) * K;
    else if( m_Params.m_RayleighCoeff[0] != Real(0) )
        A = (Real(1) + dt*m_Params.m_RayleighCoeff[0]) * M + mal::Sq(dt) * K; //Optimization: M is diagonal
    else if( m_Params.m_RayleighCoeff[1] != Real(0) )
        A = M + (mal::Sq(dt) + dt*m_Params.m_RayleighCoeff[1]) * K;
    else
        A = M + mal::Sq(dt) * K;

    // Assemble KNC
#ifdef __ENABLE_KNC
    if( m_Params.m_KinematicMode != Params::eKM_Disabled )
    {
        for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int i( it_knc->m_nid );
            for( unsigned int j=0; j<m_NumNodes; j++ )
            {
                A( 2*i  , 2*j   ) = 0;
                A( 2*i  , 2*j+1 ) = 0;
                A( 2*i+1, 2*j   ) = 0;
                A( 2*i+1, 2*j+1 ) = 0;
                /* Symm NO, because j DOES DEPEND on i
                   A( 2*j  , 2*i   ) = 0;
                   A( 2*j+1, 2*i   ) = 0;
                   A( 2*j  , 2*i+1 ) = 0;
                   A( 2*j+1, 2*i+1 ) = 0;
                */
            }
            // Diag = Id
            A( 2*i  , 2*i   ) = 1;
            A( 2*i  , 2*i+1 ) = 0;
            A( 2*i+1, 2*i   ) = 0;
            A( 2*i+1, 2*i+1 ) = 1;
        }
    }
#endif
}

//TEMP
void inline LameParameters_From_YoungAndPoisson( Real young_modulus, Real poisson_ratio,
                                                 Real &lame_mu, Real &lame_lambda )
{
    //\todo THIS MAY BE DIFFERENT IN 2D ?!?!?!?!
    lame_mu = Real(0.5) * young_modulus / (1+poisson_ratio);
    lame_lambda = young_modulus * poisson_ratio / ( (1+poisson_ratio) * (1-2*poisson_ratio) );
}

void LeafDSH_FEM_Solid2D::Compute_Effective_LameParams( Real det_F, Real &lame_mu, Real &lame_lambda ) const
{
    if( m_Params.m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
        LameParameters_From_YoungAndPoisson( m_Params.m_YoungModulus, Compute_Effective_PoissonRatio(det_F), lame_mu, lame_lambda );
    else
        LameParameters_From_YoungAndPoisson( m_Params.m_YoungModulus, m_Params.m_PoissonRatio, lame_mu, lame_lambda );
}

Vec2 LeafDSH_FEM_Solid2D::ComputeDiagonalP_ITF_LRM( const Vec2 &vec_diag_F ) const
{
    /* TEMP: Literal implementation from ITF, slower
       Mat2x2r diag_F( vec_diag_F[0], 0,
       0, vec_diag_F[1] );
       Mat2x2r E( diag_F - Mat2x2r::Identity() );
       Mat2x2r diag_P( 2*m_LameMu * E
       + m_LameLambda * mal::Trace( E ) * Mat2x2r::Identity() );
    */
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real lame_mu,lame_lambda;
    Compute_Effective_LameParams( det_F, lame_mu, lame_lambda );
    Real trE( vec_diag_F[0] + vec_diag_F[1] - 2 );
    return Vec2( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * trE,
                 2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * trE );
}

Vec2 LeafDSH_FEM_Solid2D::ComputeDiagonalP_ECIE_CLRM( const Vec2 &vec_diag_F ) const
{
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real lame_mu,lame_lambda;
    Compute_Effective_LameParams( det_F, lame_mu, lame_lambda );
    Real J( det_F );
    return Vec2( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * (J-1) * vec_diag_F[1],
                 2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * (J-1) * vec_diag_F[0] );
}

Vec2 LeafDSH_FEM_Solid2D::ComputeDiagonalP_ECIE_NHC0( const Vec2 &vec_diag_F ) const
{
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real lame_mu,lame_lambda;
    Compute_Effective_LameParams( det_F, lame_mu, lame_lambda );
    Real ecie_e( m_Params.m_ECIE_e_threshold );
    Vec2 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                             mal::Max( ecie_e, vec_diag_F[1] ) );
    Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
    return Vec2( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                 lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
}

Vec2 LeafDSH_FEM_Solid2D::ComputeDiagonalP_ECIE_NHC1( const Vec2 &vec_diag_F ) const
{
    Real det_F( vec_diag_F[0] * vec_diag_F[1] );
    Real lame_mu,lame_lambda;
    Compute_Effective_LameParams( det_F, lame_mu, lame_lambda );
    Real ecie_e( m_Params.m_ECIE_e_threshold );
    Real ecie_k( m_Params.m_ECIE_k_factor * m_Params.m_YoungModulus );
    bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
    bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
    Vec2 vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                             mal::Max( ecie_e, vec_diag_F[1] ) );
    Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
    Vec2 vec_diag_P( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                     lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
    Real h01( lame_lambda / clamped_J );
    //\todo Cases can be generalized if delta0,delta1 are 0 for non-clamped axis
    if( bIsDegenerate0 && bIsDegenerate1 )
    {
        Real delta0( vec_diag_F[0] - ecie_e );
        Real delta1( vec_diag_F[1] - ecie_e );
        vec_diag_P[0] += h01*delta1 + 2*ecie_k*delta0;
        vec_diag_P[1] += h01*delta0 + 2*ecie_k*delta1;
    }
    else if( bIsDegenerate0 && !bIsDegenerate1 )
    {
        Real delta0( vec_diag_F[0] - ecie_e );
        vec_diag_P[0] += 2*ecie_k*delta0;
        vec_diag_P[1] += h01*delta0;
    }
    else if( !bIsDegenerate0 && bIsDegenerate1 )
    {
        Real delta1( vec_diag_F[1] - ecie_e );
        vec_diag_P[0] += h01*delta1;
        vec_diag_P[1] += 2*ecie_k*delta1;
    }
    return vec_diag_P;
}

Vec6 Compute_PiolaKirchhoff_Forces( const ms::fem::LinearTriangle2 &element, const Mat2x2 &diag_P, const Mat2x2 &U, const Mat2x2 &Vt )
{
    Mat2x2 Bm( element.InvDm() );
    Mat2x2 H( - element.Area() * U * diag_P * Vt * Bm.Transposed() ); // Negation is present in the Siggraph2012 course notes, but not in ITF...
    Vec2 f0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
    Vec6 f( f0.x(), f0.y(),
            H(0,0), H(1,0),
            H(0,1), H(1,1) );
    return f;
}

void LeafDSH_FEM_Solid2D::AccumulateExplicitForces( Real dt )
{
    // Apply gravity forces
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.ApplyForce( it_node, m_Model.GetMass(it_node)*m_Params.m_Gravity );
    // Apply FEM forces
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
        // Compute Fe according to method
        switch( m_Params.m_Method )
        {
            // LRM methods
        case Params::eMethod_Id:
        case Params::eMethod_QR:
        case Params::eMethod_PD:
        case Params::eMethod_PD_Reflect:
        case Params::eMethod_PD_Fix:
        case Params::eMethod_PD_Project:
        case Params::eMethod_PD_Project_Nearest:
        case Params::eMethod_PD_SVD:
        case Params::eMethod_IHFSDM:
            {
                // Gather nodes
                Vec2 x0( m_Model.GetPos( element.nid(0) ) );
                Vec2 x1( m_Model.GetPos( element.nid(1) ) );
                Vec2 x2( m_Model.GetPos( element.nid(2) ) );
                // Compute Te transform {e}->{0}
                Transform2 Te( m_vecTe2r[it_e] );
                // Compute Ue in {e} coords
                Transform2 InvTe( Te.Inverse() );
                Vec2 u0( InvTe * x0 - element.r(0) );
                Vec2 u1( InvTe * x1 - element.r(1) );
                Vec2 u2( InvTe * x2 - element.r(2) );
                Vec6 Ue( u0.x(), u0.y(),  u1.x(), u1.y(),  u2.x(), u2.y() );

#  ifdef __ENABLE_RAYLEIGH_DAMPING //\todo Convert into m_Params.m_Flags.Test( eFlag_RayleighDamping )
                // Mass matrix
                Real m0( m_Model.GetMass( element.nid(0) ) );
                Real m1( m_Model.GetMass( element.nid(1) ) );
                Real m2( m_Model.GetMass( element.nid(2) ) );
                Mat6x6 M( Mat6x6::Zero() );
                M(0,0) = m0; M(1,1) = m0;
                M(2,2) = m1; M(3,3) = m1;
                M(4,4) = m2; M(5,5) = m2;
                // Damping matrix
                Mat6x6 C( m_Params.m_RayleighCoeff[0] * M + m_Params.m_RayleighCoeff[1] * element.K() ); //\todo we ignore Compute_Effective_Ke() here?.... IF NOT, damping will be MORE BISED THAN ELASTICITY....
                // Local node velocities
                Vec2 v0( InvTe.m_Rot * m_Model.GetVel( element.nid(0) ) );
                Vec2 v1( InvTe.m_Rot * m_Model.GetVel( element.nid(1) ) );
                Vec2 v2( InvTe.m_Rot * m_Model.GetVel( element.nid(2) ) );
                Vec6 Ve( v0.x(), v0.y(),  v1.x(), v1.y(),  v2.x(), v2.y() );

#    ifdef __ENABLE_RIGID_MODE_NO_DAMPING //\todo Convert into a function RemoveRigidModeVel(m0,m1,m2,r0,r1,r2,u0,u1,u2,v0,v1,v2) and call if m_Params.m_Flags.Test( eFlag_NoRigidModeDamping )
                //---- Substract rigid mode velocity from Ve to avoid damping it unrealistically
                Real inv_m( mal::Rcp( m0 + m1 + m2 ) );
                Vec2 r_cm( inv_m * ( m0 * (element.r(0) + u0)
                                     + m1 * (element.r(1) + u1)
                                     + m2 * (element.r(2) + u2) ) );
                Vec2 v_cm( inv_m * ( m0*v0 + m1*v1 + m2*v2 ) );
                // Substract rigid CM vel
                v0 -= v_cm;
                v1 -= v_cm;
                v2 -= v_cm;
                // Compute node displacements from CM
                Vec2 d0p( mal::PerpendicularCW( element.r(0) + u0 - r_cm ) );
                Vec2 d1p( mal::PerpendicularCW( element.r(1) + u1 - r_cm ) );
                Vec2 d2p( mal::PerpendicularCW( element.r(2) + u2 - r_cm ) );
                // Compute node angular momentum L_i
                Real l0( m0 * (d0p * v0) );
                Real l1( m1 * (d1p * v1) );
                Real l2( m2 * (d2p * v2) );
                // Compute rigid rotation vel \omega_cm = L_cm / I_cm
                Real omega_cm( ( l0 + l1 + l2 )
                               /
                               ( m0*mal::NormSq(d0p) + m1*mal::NormSq(d1p) + m2*mal::NormSq(d2p) ) );
                // Substract rigid rotation vel
                v0 -= d0p*omega_cm;
                v1 -= d1p*omega_cm;
                v2 -= d2p*omega_cm;
#    endif //__ENABLE_RIGID_MODE_NO_DAMPING

                Vec6 Vd( v0.x(), v0.y(),  v1.x(), v1.y(),  v2.x(), v2.y() );
                // Compute Fe in {e} coords
                Vec6 Fd( -(C * Vd) );
                Vec6 Fs( -(Compute_Effective_Ke(it_e) * Ue) );
                Vec6 Fe( Fd + Fs );

#else  //__ENABLE_RAYLEIGH_DAMPING
                // Compute Fe in {e} coords
                Vec6 Fe( -(Compute_Effective_Ke() * Ue) );
#endif //__ENABLE_RAYLEIGH_DAMPING
                Vec2 f0( Fe[0], Fe[1] );
                Vec2 f1( Fe[2], Fe[3] );
                Vec2 f2( Fe[4], Fe[5] );
                // Apply Fe in {0} coords
                m_Model.ApplyForce( element.nid(0), Te.Rot()*f0 );
                m_Model.ApplyForce( element.nid(1), Te.Rot()*f1 );
                m_Model.ApplyForce( element.nid(2), Te.Rot()*f2 );
            }
            break;
            // P-K methods
        case Params::eMethod_ITF_LRM:
        case Params::eMethod_ECIE_CLRM:
        case Params::eMethod_ECIE_NHC0:
        case Params::eMethod_ECIE_NHC1:
            {
                Mat2x2 F( m_vecF[it_e] );
                Real det_F( m_vecDetF[it_e] );
                bool b_fix_f( false );//true ); //\todo param
#define __PROJECT_F
#ifdef __PROJECT_F //\todo Param Clamp/Reflect
                // Project F wrt NoC
                if( b_fix_f && det_F < m_Params.m_DegenerateThresholdDetF )
                {
                    int noc = FindNoC( it_e );
                    if( noc != geo::cInvalidFeatureIndex )
                    {
                        // Gather pos
                        Vec2 vec_pos[3];
                        vec_pos[0] = m_Model.GetPos( element.nid(0) );
                        vec_pos[1] = m_Model.GetPos( element.nid(1) );
                        vec_pos[2] = m_Model.GetPos( element.nid(2) );

                        //DS_LOG_WARNING("Using NoC %d", noc );
                        Vec2 &x0( vec_pos[ noc ] ); //Inverted vtx
                        const Vec2 &x1( vec_pos[ (noc + 1) % 3 ] );
                        const Vec2 &x2( vec_pos[ (noc + 2) % 3 ] );
                        Vec2 axis12( mal::PerpendicularCW( x2 - x1 ));
                        Real dist12( mal::Norm(axis12) );
                        if( dist12 > 0.0001f )
                        {
                            axis12 /= dist12;
                            // Project
                            Real delta( m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist12 );
                            x0 -= mal::Dot(x0-x1,axis12) * axis12; //Project on-to axis
                            x0 += delta * axis12; //Move past collapse area, to ensure numerically safe
                            // Recompute H and R from it
                            element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], F );
                            det_F = mal::Det(F);
                        }
                        else
                        {
                            F = Mat2x2::Identity();
                            det_F = 1;
                        }
                    }
                    else
                    {
                        F = Mat2x2::Identity();
                        det_F = 1;
                    }
                }
#else //__REFLECT_F //\todo REFLECT GENERATES REFLECTED FORCES, which is WRONG...
                if( b_fix_f && det_F < 0 )
                {
                    int noc = FindNoC( it_e );
                    if( noc != geo::cInvalidFeatureIndex )
                    {
                        // Gather pos
                        Vec2 vec_pos[3];
                        vec_pos[0] = m_Model.GetPos( element.nid(0) );
                        vec_pos[1] = m_Model.GetPos( element.nid(1) );
                        vec_pos[2] = m_Model.GetPos( element.nid(2) );

                        //DS_LOG_WARNING("Using NoC %d", noc );
                        Vec2 &x0( vec_pos[ noc ] ); //Inverted vtx
                        const Vec2 &x1( vec_pos[ (noc + 1) % 3 ] );
                        const Vec2 &x2( vec_pos[ (noc + 2) % 3 ] );
                        Vec2 axis12( mal::PerpendicularCW( x2 - x1 ));
                        Real dist12( mal::Norm(axis12) );
                        if( dist12 > 0.0001f )
                        {
                            axis12 /= dist12;
                            // Reflect
                            F = mal::GReflection2x2_From_Axis( axis12 ) * F;
                            det_F = mal::Det(F);
                        }
                        else
                        {
                            F = Mat2x2::Identity();
                            det_F = 1;
                        }
                    }
                    else
                    {
                        F = Mat2x2::Identity();
                        det_F = 1;
                    }
                }
#endif
                // Compute SVD
                Mat2x2 U, Vt;
                Vec2 vec_diag_F;
                mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
                // Force pure rotation U,Vt by fixing potential inversion/reflection
                if( mal::Det(F) < 0 ) //\todo Unnecessary if
                {
                    if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                    if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
                }
                // Compute P and its derived force
                Vec2 vec_diag_P;
                switch( m_Params.m_Method )
                {
                case Params::eMethod_ITF_LRM: vec_diag_P = ComputeDiagonalP_ITF_LRM(vec_diag_F); break;
                case Params::eMethod_ECIE_CLRM: vec_diag_P = ComputeDiagonalP_ECIE_CLRM(vec_diag_F); break;
                case Params::eMethod_ECIE_NHC0: vec_diag_P = ComputeDiagonalP_ECIE_NHC0(vec_diag_F); break;
                case Params::eMethod_ECIE_NHC1: vec_diag_P = ComputeDiagonalP_ECIE_NHC1(vec_diag_F); break;
                default: vec_diag_P = Vec2::Zero(); DS_ASSERT(false); break;
                }
                Mat2x2 diag_P( mal::GMatNxN_From_Diagonal( vec_diag_P ) );
                Vec6 fe( Compute_PiolaKirchhoff_Forces( element, diag_P, U, Vt ) );
                Vec2 f0( fe[0], fe[1] );
                Vec2 f1( fe[2], fe[3] );
                Vec2 f2( fe[4], fe[5] );
                // Apply Fe in {0} coords
                m_Model.ApplyForce( element.nid(0), f0 );
                m_Model.ApplyForce( element.nid(1), f1 );
                m_Model.ApplyForce( element.nid(2), f2 );
            }
            break;
        case Params::eMethod_PD_CLRM:
            {
                Mat2x2 F( m_vecF[it_e] );
                Real det_F( m_vecDetF[it_e] );
                Mat2x2 fixedF( F );
                Real det_fixedF( det_F );
                bool b_fix_f( true ); //\todo param
#define __PROJECT_F
#ifdef __PROJECT_F
                // Project F wrt NoC
                if( b_fix_f && det_F < m_Params.m_DegenerateThresholdDetF )
                {
                    int noc = FindNoC( it_e );
                    if( noc != geo::cInvalidFeatureIndex )
                    {
                        // Gather pos
                        Vec2 vec_pos[3];
                        vec_pos[0] = m_Model.GetPos( element.nid(0) );
                        vec_pos[1] = m_Model.GetPos( element.nid(1) );
                        vec_pos[2] = m_Model.GetPos( element.nid(2) );

                        //DS_LOG_WARNING("Using NoC %d", noc );
                        Vec2 &x0( vec_pos[ noc ] ); //Inverted vtx
                        const Vec2 &x1( vec_pos[ (noc + 1) % 3 ] );
                        const Vec2 &x2( vec_pos[ (noc + 2) % 3 ] );
                        Vec2 axis12( mal::PerpendicularCW( x2 - x1 ));
                        Real dist12( mal::Norm(axis12) );
                        if( dist12 > 0.0001f )
                        {
                            axis12 /= dist12;
                            // Project
                            Real delta( m_Params.m_DegenerateThresholdDetF * 2*element.Area() / dist12 );
                            x0 -= mal::Dot(x0-x1,axis12) * axis12; //Project on-to axis
                            x0 += delta * axis12; //Move past collapse area, to ensure numerically safe
                            // Recompute H and R from it
                            element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], fixedF );
                            det_fixedF = mal::Det(fixedF);
                        }
                        else
                        {
                            fixedF = Mat2x2::Identity();
                            det_fixedF = 1;
                        }
                    }
                    else
                    {
                        fixedF = Mat2x2::Identity();
                        det_fixedF = 1;
                    }
                }
#else //__REFLECT_F //\todo REFLECT GENERATES REFLECTED FORCES, which is WRONG...
                if( b_fix_f && det_F < 0 )
                {
                    int noc = FindNoC( it_e );
                    if( noc != geo::cInvalidFeatureIndex )
                    {
                        // Gather pos
                        Vec2 vec_pos[3];
                        vec_pos[0] = m_Model.GetPos( element.nid(0) );
                        vec_pos[1] = m_Model.GetPos( element.nid(1) );
                        vec_pos[2] = m_Model.GetPos( element.nid(2) );

                        //DS_LOG_WARNING("Using NoC %d", noc );
                        Vec2 &x0( vec_pos[ noc ] ); //Inverted vtx
                        const Vec2 &x1( vec_pos[ (noc + 1) % 3 ] );
                        const Vec2 &x2( vec_pos[ (noc + 2) % 3 ] );
                        Vec2 axis12( mal::PerpendicularCW( x2 - x1 ));
                        Real dist12( mal::Norm(axis12) );
                        if( dist12 > 0.0001f )
                        {
                            axis12 /= dist12;
                            // Reflect
                            fixed_F = mal::GReflection2x2_From_Axis( axis12 ) * F;
                            det_fixedF = -det_F;
                        }
                        else
                        {
                            fixedF = Mat2x2::Identity();
                            det_fixedF = 1;
                        }
                    }
                    else
                    {
                        fixedF = Mat2x2::Identity();
                        det_fixedF = 1;
                    }
                }
#endif
                Real lame_mu,lame_lambda;
                Compute_Effective_LameParams( det_F, lame_mu, lame_lambda );

                // PD of fixed F
                Mat2x2 R( mal::GRotation2x2_PolarDecomposition( fixedF, det_fixedF ) );
                Mat2x2 S( R.Transposed() * F ); //F = R*S
                Real det_S( mal::Det(S) ); //\note == det_F
                Mat2x2 P( R * ( 2*lame_mu*(S-Mat2x2::Identity())
                                + lame_lambda * (det_S-1) * Mat2x2( S(1,1), -S(1,0),
                                                                    -S(0,1), S(0,0) ) ) );
                // From Compute_PiolaKirchhoff_Forces()
                Mat2x2 Bm( element.InvDm() );
                Mat2x2 H( - element.Area() * P * Bm.Transposed() );
                Vec2 f1( H(0,0), H(1,0) );
                Vec2 f2( H(0,1), H(1,1) );
                Vec2 f0( - f1 - f2 ); //h0 = -h1-h2;
                // Apply Fe in {0} coords
                m_Model.ApplyForce( element.nid(0), f0 );
                m_Model.ApplyForce( element.nid(1), f1 );
                m_Model.ApplyForce( element.nid(2), f2 );
            }
            break;
        default: break;
        }

#ifdef __S2_DS_ENABLE_STATS
        {
            Vec2 x0( m_Model.GetPos( element.nid(0) ) );
            Vec2 x1( m_Model.GetPos( element.nid(1) ) );
            Vec2 x2( m_Model.GetPos( element.nid(2) ) );
            m_Stats.m_RelArea.Add( mal::Abs( element.Compute_Area(x0,x1,x2) / element.Area() ) );
        }
#endif

    }

#ifdef __ENABLE_SAVE_FORCES
    for( unsigned int i=0; i<m_NumNodes; i++ ) m_vecForces[i] = m_Model.GetAccForce(i);
#endif

    // Cancel some f^e
#ifdef __ENABLE_KNC
    if( m_Params.m_KinematicMode != Params::eKM_Disabled )
    {
        for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int nid( it_knc->m_nid );
            m_Model.ApplyForce( nid, -m_Model.GetAccForce(nid) );
            m_Model.SetPos( nid, it_knc->m_Pos );
            m_Model.SetVel( nid, it_knc->m_Vel );
        }
    }
#endif

#ifdef __ENABLE_CONTACT_CIRCLE
    if( m_LastRPAP.m_Radius > 0 )
    {
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        {
            Vec2 r( m_Model.GetPos(it_node)-m_LastRPAP.m_Pos );
            Real dist = r.Norm();
            if( dist < m_LastRPAP.m_Radius )
            {
                Vec2 u( r / dist );
                // Position correction
                m_Model.SetPos( it_node, m_Model.GetPos(it_node) + m_Params.m_ContactSolver_Relaxation_Coeff * (m_LastRPAP.m_Radius-dist) * u );
                // Could apply impulse instead of vel change...
                Vec2 vel( m_Model.GetVel(it_node) );
                Real vel_n( mal::Dot( vel, u ) );
                if( vel_n < Real(0) ) m_Model.SetVel( it_node, vel - (Real(1)+m_Params.m_ContactSolver_Restitution_Coeff) * vel_n * u );
                // Constrait force
                Real force_n( mal::Dot(m_Model.GetAccForce(it_node),u) );
                if( force_n < Real(0) ) m_Model.ApplyForce( it_node, -force_n*u );
            }
        }
    }
#endif

#ifdef __ENABLE_CONTACT_D2K
    /*
    DS_LOG_ASSERT( false, "__ENABLE_CONTACT_D2K for explicit integration usupported, must migrate to POF instead of FID" );
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            //TEMP
            //geo::feature_index_type vid1( cpc.m_FID1.AsVertex() );
            geo::feature_index_type vid1(geo::cInvalidFeatureIndex);
            if( cpc.m_FID1.IsVertex() ) vid1 = cpc.m_FID1.AsVertex();
            else if( cpc.m_FID1.IsEdge() ) vid1 = m_pMeshS->HE_OriginVID( cpc.m_FID1.AsEdge() );
            else { DS_LOG_ASSERT( false, "Unsupported feature_id" ); }
            //END TEMP
            // Position correction \todo AIXO FA GUANYAR ENERGIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            m_Model.SetPos( vid1, m_Model.GetPos(vid1) + (cCoeffRelaxation * cpc.m_Depth) * cpc.m_Normal );
            // Could apply impulse instead of vel change...
            Vec2 vel( m_Model.GetVel(vid1) );
            Real vel_n( mal::Dot( vel, cpc.m_Normal ) );
            if( vel_n < Real(0) ) m_Model.SetVel( vid1, vel - (Real(1)+cCoeffRestitution) * vel_n * cpc.m_Normal );
            // Constrait force
            Real force_n( mal::Dot( m_Model.GetAccForce(vid1), cpc.m_Normal ) );
            if( force_n < Real(0) ) m_Model.ApplyForce( vid1, -force_n*cpc.m_Normal );
        }
    }
    */
#endif
}

void LeafDSH_FEM_Solid2D::RecomputeEnergy()
{
    m_ElasticEnergy = Real(0);
    m_PotentialEnergy = Real(0);
    m_KineticEnergy = Real(0);
    for( unsigned int it_e=0; it_e<m_NumElements; it_e++ )
    {
        const ms::fem::LinearTriangle2 &element( m_vecElements[it_e] );
        // Gather element data
        Vec2 x0( m_Model.GetPos( element.nid(0) ) );
        Vec2 x1( m_Model.GetPos( element.nid(1) ) );
        Vec2 x2( m_Model.GetPos( element.nid(2) ) );
        Real element_mass( m_DensityPerArea * element.Area() );
        Transform2 InvTe( m_vecTe2r[it_e].Inverse() );

        // Elastic energy
        Vec2 u0( InvTe * x0 - element.r(0) );
        Vec2 u1( InvTe * x1 - element.r(1) );
        Vec2 u2( InvTe * x2 - element.r(2) );
        Vec6 Ue( u0.x(), u0.y(),  u1.x(), u1.y(),  u2.x(), u2.y() );
        m_ElasticEnergy += Real(0.5)*(Ue * (element.K() * Ue)); //== \frac{1}{2} (u^e)^T * K^e * u^e //\todo We IGORE Compute_Effective_Ke() here?...

        // Consistent kinetic energy? see FEM.pdf
        Vec2 v0( m_Model.GetVel( element.nid(0) ) );
        Vec2 v1( m_Model.GetVel( element.nid(1) ) );
        Vec2 v2( m_Model.GetVel( element.nid(2) ) );
        m_KineticEnergy += ( element_mass / Real(12) )
                           * ( mal::NormSq(v0) + mal::NormSq(v1) + mal::NormSq(v2)
                               + v0.x()*v1.x() + v1.x()*v2.x() + v0.x()*v2.x()
                               + v0.y()*v1.y() + v1.y()*v2.y() + v0.y()*v2.y() );

        // Consistent potential energy?, see FEM.pdf
        m_PotentialEnergy += ( element_mass / Real(3) ) *  mal::Dot( x0+x1+x2, -m_Params.m_Gravity );
    }
}

/*! Integrate dynamics from (x_0,v_0) to (x_1,v_1)
 */
void LeafDSH_FEM_Solid2D::Solve_SymplecticEuler( Real dt )
{
    // Accumulate forces
    AccumulateExplicitForces(dt);
    // Update PS (internally symplectic euler)
    m_Model.Update(dt);
}

/*! Integrate dynamics from (x_0,v_0) to (x_1,v_1)
 */
void LeafDSH_FEM_Solid2D::Solve_ImplicitEuler( Real dt )
{
#ifdef __ENABLE_CONTACT_CIRCLE
    if( m_LastRPAP.m_Radius > 0
        && m_Params.m_KinematicMode != Params::eKM_Mouse )
    {
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        {
            Vec2 r( m_Model.GetPos(it_node)-m_LastRPAP.m_Pos );
            Real dist = r.Norm();
            if( dist < m_LastRPAP.m_Radius )
            {
                Vec2 u( r / dist );
                // Position correction
                m_Model.SetPos( it_node, m_Model.GetPos(it_node) + m_Params.m_ContactSolver_Relaxation_Coeff * (m_LastRPAP.m_Radius-dist) * u );
                // Could apply impulse instead of vel change...
                Vec2 vel( m_Model.GetVel(it_node) );
                Real vel_n( mal::Dot( vel, u ) );
                if( vel_n < Real(0) ) m_Model.SetVel( it_node, vel - (Real(1)+m_Params.m_ContactSolver_Restitution_Coeff) * vel_n * u );
                /* Constrait force
                Real force_n( mal::Dot(m_Model.GetAccForce(it_node),u) );
                if( force_n < Real(0) ) m_Model.ApplyForce( it_node, -force_n*u );
                */
            }
        }
    }
#endif

    // Assemble common stuff \todo Some are invariant if !corotational
    Assemble_M( m_M );
    Assemble_K( m_K );
    Assemble_A( m_A, dt, m_M, /*m_C,*/ m_K );
    Assemble_Fe( m_Fe );
    ns::VectorD v_0;
    Assemble_v( v_0 );
    // Assemble b
    ns::VectorD b;
    if( m_Params.m_Method != Params::eMethod_Id ) //\todo By now, only LFEM does not use corotational, but this will become a switch(EMethod) someday
    {
        ns::VectorD cef_0;
        Assemble_CorotationalElasticForce( cef_0 );
        b = m_M*v_0 + dt * ( m_Fe + cef_0 );
    }
    else // Linear FEM
    {
        ns::VectorD u_0;
        Assemble_u( u_0 );
        b = m_M*v_0 + dt * ( m_Fe - m_K*u_0 );
    }

    // Assemble KNC
#ifdef __ENABLE_KNC
    if( m_Params.m_KinematicMode != Params::eKM_Disabled )
    {
        for( PoolKNC::iterator it_knc=m_poolKNC.Begin(); it_knc.IsValid(); ++it_knc )
        {
            unsigned int i( it_knc->m_nid );
            KinematicNodeConstraint *pKNC( it_knc.GetPtr() );
            /* No automatic control by now, all KNC positions set explicitly from API
            // TEMP: Compute target position
            switch( m_Params.m_KinematicMode )
            {
            case Params::eKM_Static:
                // KNC remains static
                break;
            case Params::eKM_Periodic:
                // Kinematic controller
                if( it_knc == 0 ) pKNC->m_Pos = Vec2( 10.0f * mal::Sin(m_TotalTime), 0.0f );
                break;
            case Params::eKM_Mouse:
                if( it_knc == 0 ) pKNC->m_Pos = m_LastRPAP.m_Pos;
                break;
            default: DS_ASSERT(false);
            }
            */

            // If it's a position constraint, Compute vel from position
            // error, otherwise, it's already a velocity constraint
            if( pKNC->m_Type == KinematicNodeConstraint::eKNCT_Position )
            {
                Vec2 vk( ( pKNC->m_Pos-m_Model.GetPos(i) ) / dt );
                Real norm_vk( vk.Norm() );
                if( norm_vk > Real(0.001f) )
                    pKNC->m_Vel = mal::Min( norm_vk, 50.0f ) * vk.Normalized();
                else
                    pKNC->m_Vel = Vec2(0,0);
            }

            // Set b to KNC vel
            b[2*i]   = pKNC->m_Vel[0];
            b[2*i+1] = pKNC->m_Vel[1];

#  ifdef __ENABLE_KINEMATIC_B
            // Add A_ji * Vi+1 to ALL non-kinematic particles j
            for( unsigned int j=0; j<m_NumNodes; j++ )
            {
                // Compute constant factor for i,j
                if( !IsKNC(j) )
                {
                    b[2*j  ] -= m_A( 2*j  , 2*i ) * b[2*i] + m_A( 2*j  , 2*i+1 ) * b[2*i+1];
                    b[2*j+1] -= m_A( 2*j+1, 2*i ) * b[2*i] + m_A( 2*j+1, 2*i+1 ) * b[2*i+1];
                }
                // Symmetrize A_ij for proper SolveCG
                if( j != i )
                {
                    m_A(2*j  ,2*i  ) = 0;
                    m_A(2*j  ,2*i+1) = 0;
                    m_A(2*j+1,2*i  ) = 0;
                    m_A(2*j+1,2*i+1) = 0;
                }
            }
#  endif

        }
    }
#endif //__ENABLE_KNC

#ifdef __ENABLE_TRACE_IE
    std::cout << "T = " << m_TotalTime << std::endl;
    static unsigned int s_CountIE(0);
    bool bOutputIE(false);
    bOutputIE = ( 0 == s_CountIE % 100 );
    s_CountIE++;
    if( bOutputIE )
    {
        std::cout << "M = " << std::endl << m_M << std::endl;
        std::cout << "K = " << std::endl << m_K << std::endl;
        //std::cout << "C = " << std::endl << C << std::endl;
        std::cout << "A = " << std::endl << m_A << std::endl;
        std::cout << "b = " << std::endl << b << std::endl;
        //std::cout << "u_0 = " << std::endl << u_0 << std::endl;
        //std::cout << "Fe_0 = " << std::endl << Fe_0 << std::endl;
        std::cout << "v_0 = " << std::endl << v_0 << std::endl;
    }
#endif

    // Solve for v_1
    ns::VectorD v_1;
    v_1.Resize( b.Size() );
    v_1.Zero();
    //
    //ns::VectorD v_1(v_0); //TEMPORAL: Warmstart v_1 as v_0, does NOT seem to reduce #iter...

    // SolverLS params
    ns::RealD ls_prec( m_Params.m_SolverLS_Epsilon );
    int ls_iter( m_Params.m_SolverLS_MaxIter );
    /*\todo Optimize using diagonal M and sparse K with combined timestep and damping factors
    SolveCG_M_K( factor_M, diagonal_M,
                 factor_K, K,
                 v_1,
                 b,
                 ls_prec, ls_iter );
    */
    switch( m_Params.m_SolverLS_Type )
    {
    case Params::eLSST_CG:
        ns::SolveCG( m_A, v_1, b, ls_prec, ls_iter );
        break;
    case Params::eLSST_GS:
        //\todo Works, but it's slower than CG
        ns::SolveGS( m_A, v_1, b, ls_prec, ls_iter );
        break;
    case Params::eLSST_Jacobi:
        //\todo Diverges miserably, bug or fails to converge?
        ns::SolveJacobi( m_A, v_1, b, ls_prec, ls_iter );
        break;
    case Params::eLSST_LU:
        //\todo Works, but is VERY VERY SLOW wrt CG
        SolveLU( m_A, v_1, b, ls_prec );
        break;
    default:
        DS_ASSERT(false);
    }

#ifdef __S2_DS_ENABLE_STATS
    m_Stats.m_SolverLS_NumIter = ls_iter;
    m_Stats.m_SolverLS_Prec = ls_prec;
#endif

#ifdef __ENABLE_TRACE_IE
    if( bOutputIE )
    {
        std::cout << "v_1 = " << std::endl << v_1 << std::endl;
        std::cout << "ls_iter = " << ls_iter << std::endl;
        std::cout << "ls_prec = " << ls_prec << std::endl;
    }
#endif

    /*\TODO HERE we could predict element inversion using p(t) and
      predicted p(t+dt), and avoid-it by means of a velocity filter
      (conditional incompressibility, probably iterative/PGS-like or
      even PBD relaxation-like...)

      Some kind of direct-global or iterative-local-propagation
      solution is required because velocity changes affect several
      elements in cascade... ideally only 1-ring or 2-ring
      neighbourhoods, so low iteration count should be enough.
    */

    // Integrate x_1 with v_1
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        m_Model.SetVel( it_node, Vec2( v_1[2*it_node], v_1[2*it_node+1] ) );
        m_Model.SetPos( it_node, m_Model.GetPos(it_node) + dt*m_Model.GetVel(it_node) );
    }

#ifdef __ENABLE_CONTACT_D2K_PRE //\todo Hack contact can only work as pre, not post...
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        //if( !ccd2k.IsEmpty() ) DS_LOG_WARNING( "Contact %llx", (machine_uint_type)&ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                switch( m_Params.m_ContactSolver_Type )
                {
                case Params::eCST_Hack:
                    {
                        if( cpc.m_POF1.m_FeatureId.IsVertex() )
                        {
                            geo::feature_index_type vid1 = cpc.m_POF1.m_FeatureId.AsVertex();
                            // pos,vel,acc correction constraints
                            // Position correction \todo AIXO FA GUANYAR ENERGIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            m_Model.SetPos( vid1, m_Model.GetPos(vid1) + (m_Params.m_ContactSolver_Relaxation_Coeff * cpc.m_Depth) * cpc.m_Normal );
                            // Could apply impulse instead of vel change...
                            Vec2 vel( m_Model.GetVel(vid1) );
                            Real vel_n( mal::Dot( vel, cpc.m_Normal ) );
                            if( vel_n < Real(0) ) m_Model.SetVel( vid1, vel - (Real(1)+m_Params.m_ContactSolver_Restitution_Coeff) * vel_n * cpc.m_Normal );
                            // Constrait force
                            Real force_n( mal::Dot( m_Model.GetAccForce(vid1), cpc.m_Normal ) );
                            if( force_n < Real(0) ) m_Model.ApplyForce( vid1, -force_n*cpc.m_Normal );
                        }
                        else if( cpc.m_POF1.m_FeatureId.IsSegment() )
                        {
                            DS_LOG_ERROR("Unsupported feature_id IsEdge()");
                        }
                        else
                        {
                            DS_LOG_ASSERT( false, "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type );
                        }
                    }
                    break;
                default: break;
                }
            }
        }
    }
    if( Params::eCST_HackGroundPlane == m_Params.m_ContactSolver_Type )
    {
        for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        {
            Vec2 normal( 0, 1 );
            Vec2 pos( m_Model.GetPos(it_node) );
            Real depth( -pos.y()-1 ); //GroundPlane at Y=-1
            if( depth > m_Params.m_ContactSolver_DepthMargin )
            {
                m_Model.SetPos( it_node, m_Model.GetPos(it_node) + (m_Params.m_ContactSolver_Relaxation_Coeff * depth) * normal );
                // Could apply impulse instead of vel change...
                Vec2 vel( m_Model.GetVel(it_node) );
                Real vel_n( mal::Dot( vel, normal ) );
                if( vel_n < Real(0) ) m_Model.SetVel( it_node, vel - (Real(1)+m_Params.m_ContactSolver_Restitution_Coeff) * vel_n * normal );
                // Constrait force
                Real force_n( mal::Dot( m_Model.GetAccForce(it_node), normal ) );
                if( force_n < Real(0) ) m_Model.ApplyForce( it_node, -force_n*normal );
            }
        }
    }
#endif

#ifdef __ENABLE_CONTACT_D2K_POST
    DS_LOG_ASSERT( false, "__ENABLE_CONTACT_D2K_POST usupported, must migrate to POF instead of FID" );
    /*
    for( PoolCCD2K::iterator it_ccd2k=m_poolCCD2K.Begin(); it_ccd2k.IsValid(); ++it_ccd2k )
    {
        const ContactConstraintD2K &ccd2k( *it_ccd2k );
        //if( !ccd2k.IsEmpty() ) DS_LOG_WARNING( "Contact %llx", (machine_uint_type)&ccd2k );
        for( unsigned int it_cpc=0; it_cpc<ccd2k.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2K &cpc( ccd2k.GetCPC(it_cpc) );
            // TEMP
            //geo::feature_index_type vid1( cpc.m_FID1.AsVertex() );
            geo::feature_index_type vid1;
            if( cpc.m_FID1.IsVertex() ) vid1 = cpc.m_FID1.AsVertex();
            else if( cpc.m_FID1.IsEdge() ) vid1 = m_pMeshS->HE_OriginVID( cpc.m_FID1.AsEdge() );
            else { DS_LOG_ASSERT( false, "Unsupported feature_id" ); }
            //END TEMP

            switch( m_Params.m_ContactSolver_Type )
            {
            case Params::eCST_None: break;
            case Params::eCST_Penalty:
                {
                    // penalty force
                    Vec2 vel( m_Model.GetVel(vid1) );
                    Real force_n( m_Params.m_ContactSolver_Penalty_Ks * cpc.m_Depth
                                  - m_Params.m_ContactSolver_Penalty_Kd * mal::Dot(vel,cpc.m_Normal) );
                    //m_Model.ApplyForce( vid1, force_n*cpc.m_Normal );
                    m_Model.SetVel( vid1, vel + dt*force_n*cpc.m_Normal );
                }
                break;
            case Params::eCST_Hack:
                {
                    // pos,vel,acc correction constraints
                    // Position correction \todo AIXO FA GUANYAR ENERGIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    m_Model.SetPos( vid1, m_Model.GetPos(vid1) + (cCoeffRelaxation * cpc.m_Depth) * cpc.m_Normal );
                    // Could apply impulse instead of vel change...
                    Vec2 vel( m_Model.GetVel(vid1) );
                    Real vel_n( mal::Dot( vel, cpc.m_Normal ) );
                    if( vel_n < Real(0) ) m_Model.SetVel( vid1, vel - (Real(1)+cCoeffRestitution) * vel_n * cpc.m_Normal );
                    // Constrait force
                    Real force_n( mal::Dot( m_Model.GetAccForce(vid1), cpc.m_Normal ) );
                    if( force_n < Real(0) ) m_Model.ApplyForce( vid1, -force_n*cpc.m_Normal );
                }
            default: break;
            }
        }
    }
    */
#endif //__ENABLE_CONTACT_D2K_POST
}


void LeafDSH_FEM_Solid2D::RebuildAirDrag()
{
    m_AirDragPerSecond = -mal::Log(10e-6) * m_Params.m_AirDragCoeff;
}

void LeafDSH_FEM_Solid2D::RebuildMass()
{
    // Compute consistent lumped masses
    //\todo For given mass: m_DensityPerArea = m_TotalMass / m_TotalArea;
    m_DensityPerArea = m_Params.m_Density * m_Params.m_Thickness;
    m_TotalMass = m_DensityPerArea * m_TotalArea;
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
    {
        Real weighted_area(0);
        for( geo::MeshSolidShape2::iterator_polygons_around_vertex_ccw it_ccw( m_pMeshS->GetIterator_PolygonsAroundVertexCCW(it_node) );
             it_ccw.IsValid();
             ++it_ccw )
            weighted_area += mal::Rcp(Real(3)) * m_vecElements[ it_ccw.PID() ].Area();
        m_Model.SetMass( it_node, weighted_area * m_DensityPerArea );
    }

    /*TEMPORAL HACK Uniform mass
    Real mass_per_node = m_TotalMass / m_pMeshS->GetNumV();
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.SetMass( it_node, mass_per_node );
    */

    RebuildEigenstuff();
}

void LeafDSH_FEM_Solid2D::NotifyChangedParams()
{
#ifdef __S2_DS_FEM_ENABLE_TEST
    InitTest(); //\todo This resets the test completely... consider doing it only when TestMode actually changes
#endif
    // Set material params
    for( unsigned int it_e=0; it_e < m_NumElements; it_e++ )
        m_vecElements[it_e].SetMaterialParams( m_Params.m_YoungModulus, m_Params.m_PoissonRatio );
    RebuildMass();
    RebuildAirDrag();
    RebuildEigenstuff();
}

void LeafDSH_FEM_Solid2D::RebuildEigenstuff()
{
#ifdef __ENABLE_EIGENSTUFF
    // Recompute eigenstuff
    ns::MatrixD K, InvM;
    Assemble_K( K );
    Assemble_InvM( InvM );
    ns::MatrixD A( InvM * K ); //\todo Include damping C in stability criterion!!
# ifdef __ENABLE_TRACE_EIGENSTUFF
    std::cout << "InvM = " << std::endl << InvM << std::endl;
    std::cout << "K = " << std::endl << K << std::endl;
    std::cout << "A = " << std::endl << A << std::endl;
# endif
    //Real eps_error = ns::ComputeLargestEigenvalue_Symmetric( A, m_vecEigenValue[0], m_vecEigenVector[0], 30, Real(0.001f) );
    Real eps_error = ns::ComputeEigenValues_Symmetric( A, m_vecEigenValue, m_vecEigenVector, cNumEigenValues, 100, Real(0.001f) );
    DS_ASSERT( eps_error >= Real(0) );
# ifdef __ENABLE_TRACE_EIGENSTUFF
    if( mal::Abs(eps_error/m_vecEigenValue[0]) < Real(0.1f) ) DS_LOG_WARNING("Large relative error i largest eigenvalue... maybe because A=0");
    std::cout << "eps_error = " << eps_error << std::endl;
# endif
#endif
}

unsigned int LeafDSH_FEM_Solid2D::FindClosestNode( const Vec2 &pos ) const
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

Vec2 LeafDSH_FEM_Solid2D::ComputePosCoM() const
{
    Vec2 acc_com(0,0);
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        acc_com += m_Model.GetMass( it_node ) * m_Model.GetPos( it_node );
    return acc_com / m_TotalMass;
}

//---- Edition
void LeafDSH_FEM_Solid2D::SetPosCoM( const Vec2 &pos )
{
    // Compte CoM
    Vec2 current_pos( ComputePosCoM() );
    // Compute displacement
    Vec2 displacement( pos - current_pos );
    // Apply displacement
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_Model.SetPos( it_node, m_Model.GetPos( it_node ) + displacement );
    // Set mesh DOF
    for( unsigned int it_node=0; it_node<m_NumNodes; it_node++ )
        m_pMeshO->GetSDOF(it_node) = m_Model.GetPos(it_node);
}


#ifdef __S2_DS_FEM_ENABLE_TEST
void LeafDSH_FEM_Solid2D::InitTest()
{
    // Remove existing test KNC, but ignore preexisting non-test KNC
    for( unsigned int it_tknc=0; it_tknc<m_vecTestKNC.size(); it_tknc++ )
        if( 0 != m_vecTestKNC[it_tknc] )
            m_poolKNC.Delete( m_vecTestKNC[it_tknc] );
    m_vecTestKNC.clear();
    m_vec_TM_Random_TargetPos.clear();

    // Add test KNC according to TM
    switch( m_Params.m_TestMode )
    {
    case Params::eTM_Disabled: break;
    case Params::eTM_Random:
        {
            mal::SetRandomSeed( -666 );
            Vec2 hs( 2*m_AABB0.GetHalfSizes() );
            //\note We ALWAYS allocate all test KNC, even if eTN_Boundary is selected, unused KNC will remain 0
            m_vecTestKNC.resize( m_NumNodes, static_cast<KinematicNodeConstraint*>(0) );
            for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ ) m_vec_TM_Random_TargetPos.push_back( mal::RandomV(-hs,hs) );
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                {
                    for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ )
                    {
                        if( IsKNC(it_node) ) //Pre-Existing KNC are NOT readded or changed
                            m_vecTestKNC[it_node] = 0;
                        else
                        {
                            m_vecTestKNC[it_node] = m_poolKNC.New();
                            *m_vecTestKNC[it_node] = KinematicNodeConstraint( it_node, m_pMeshS->V_Pos_0(it_node) );
                        }
                    }
                }
                break;
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            if( IsKNC(vid) ) //Pre-Existing KNC are NOT readded or changed
                                m_vecTestKNC[vid] = 0;
                            else
                            {
                                m_vecTestKNC[vid] = m_poolKNC.New();
                                *m_vecTestKNC[vid] = KinematicNodeConstraint( vid, m_pMeshS->V_Pos_0(vid) );
                            }
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
            default: break;
            }
        }
        break;
    case Params::eTM_Plane:
    case Params::eTM_Sphere:
        {
            // Create TestKNC slots, but the actual KNC are created during UpdateTest() when the constraint is violated
            m_vecTestKNC.resize( m_NumNodes, static_cast<KinematicNodeConstraint*>(0) );
        }
        break;
    default: break;
    };
}
void LeafDSH_FEM_Solid2D::UpdateTest( Real dt )
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
                            m_vecTestKNC[it_node]->m_Pos = (Real(1)-m_Params.m_TestFraction)*m_pMeshS->V_Pos_0(it_node) + m_Params.m_TestFraction*m_vec_TM_Random_TargetPos[it_node];
                }
                break;
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
                                m_vecTestKNC[vid]->m_Pos = (Real(1)-m_Params.m_TestFraction)*m_pMeshS->V_Pos_0(vid) + m_Params.m_TestFraction*m_vec_TM_Random_TargetPos[vid];
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
            default: break;
            }
        }
        break;
    case Params::eTM_Plane:
        {
            Vec2 normal( 0, 1 );
            Real half_height( 2*m_AABB0.GetHalfSizes()[1]*(Real(1)-m_Params.m_TestFraction) );
            switch( m_Params.m_TestNodes )
            {
            case Params::eTN_All:
                {
                    for( unsigned int it_node=0; it_node < m_NumNodes; it_node++ )
                    {
                        unsigned int vid( it_node );
                        Vec2 pos( m_Model.GetPos(vid) );
                        Real dot_p_n( mal::Dot( pos, normal ) );
                        Vec2 projected_pos( (dot_p_n > 0)
                                            ? pos - (dot_p_n-half_height) * normal
                                            : pos + (-dot_p_n-half_height) * normal );
                        if( 0 != m_vecTestKNC[vid] ) //KNC exists
                        {
                            if( mal::Abs( dot_p_n ) >= half_height ) //Update if behind plane
                                m_vecTestKNC[vid]->m_Pos = projected_pos;
                            else //Remove if not
                            {
                                m_poolKNC.Delete( m_vecTestKNC[vid] );
                                m_vecTestKNC[vid] = 0;
                            }
                        }
                        else if( mal::Abs( dot_p_n ) >= half_height ) //KNC does not exist, add it if behind plane
                        {
                            m_vecTestKNC[vid] = m_poolKNC.New();
                            *m_vecTestKNC[vid] = KinematicNodeConstraint( vid, projected_pos );
                        }
                    }
                }
                break;
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            Vec2 pos( m_Model.GetPos(vid) );
                            Real dot_p_n( mal::Dot( pos, normal ) );
                            Vec2 projected_pos( (dot_p_n > 0)
                                                ? pos - (dot_p_n-half_height) * normal
                                                : pos + (-dot_p_n-half_height) * normal );
                            if( 0 != m_vecTestKNC[vid] ) //KNC exists
                            {
                                if( mal::Abs( dot_p_n ) >= half_height ) //Update if behind plane
                                    m_vecTestKNC[vid]->m_Pos = projected_pos;
                                else //Remove if not
                                {
                                    m_poolKNC.Delete( m_vecTestKNC[vid] );
                                    m_vecTestKNC[vid] = 0;
                                }
                            }
                            else if( mal::Abs( dot_p_n ) >= half_height ) //KNC does not exist, add it if behind plane
                            {
                                m_vecTestKNC[vid] = m_poolKNC.New();
                                *m_vecTestKNC[vid] = KinematicNodeConstraint( vid, projected_pos );
                            }
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
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
                        Vec2 pos( m_Model.GetPos(vid) );
                        Real dist_sq( pos.NormSq() );
                        Vec2 projected_pos( (dist_sq > 0.000001f)
                                            ? radius * pos / mal::Sqrt(dist_sq)
                                            : pos );
                        if( 0 != m_vecTestKNC[vid] ) //KNC exists
                        {
                            if( dist_sq >= radius_sq ) //Update if behind plane
                                m_vecTestKNC[vid]->m_Pos = projected_pos;
                            else //Remove if not
                            {
                                m_poolKNC.Delete( m_vecTestKNC[vid] );
                                m_vecTestKNC[vid] = 0;
                            }
                        }
                        else if( dist_sq >= radius_sq ) //KNC does not exist, add it if behind plane
                        {
                            m_vecTestKNC[vid] = m_poolKNC.New();
                            *m_vecTestKNC[vid] = KinematicNodeConstraint( vid, projected_pos );
                        }
                    }
                }
                break;
            case Params::eTN_Boundary:
                {
                    // Iterate over boundary nodes only
                    for( unsigned int it_bp=0; it_bp < m_pMeshS->GetNumBoundaryP(); it_bp++ )
                    {
                        unsigned int it_he( m_pMeshS->BP_FirstHEID(it_bp) );
                        do
                        {
                            unsigned int vid( m_pMeshS->HE_OriginVID(it_he) );
                            Vec2 pos( m_Model.GetPos(vid) );
                            Real dist_sq( pos.NormSq() );
                            Vec2 projected_pos( (dist_sq > 0.000001f)
                                                ? radius * pos / mal::Sqrt(dist_sq)
                                                : pos );
                            if( 0 != m_vecTestKNC[vid] ) //KNC exists
                            {
                                if( dist_sq >= radius_sq ) //Update if behind plane
                                    m_vecTestKNC[vid]->m_Pos = projected_pos;
                                else //Remove if not
                                {
                                    m_poolKNC.Delete( m_vecTestKNC[vid] );
                                    m_vecTestKNC[vid] = 0;
                                }
                            }
                            else if( dist_sq >= radius_sq ) //KNC does not exist, add it if behind plane
                            {
                                m_vecTestKNC[vid] = m_poolKNC.New();
                                *m_vecTestKNC[vid] = KinematicNodeConstraint( vid, projected_pos );
                            }
                            it_he = m_pMeshS->HE_Next(it_he);
                        } while ( it_he != m_pMeshS->BP_FirstHEID(it_bp) );
                    }
                }
                break;
            default: break;
            }
        }
        break;
    default: break;
    }
}
#endif


#ifdef __S2_DS_ENABLE_PARAMS
void LeafDSH_FEM_Solid2D::Params::InitArchetype( util::ArchetypeLibrary &al )
{
    Params params;
    al.BeginArchetype( "Archetype_s2_ds_LDSH_FEMS2D_Params" );
    {
        al.BeginProperty_Group( "<Develop>" );
        {
            al.AddProperty( "run", params.m_bRun,
                            archetype_offset_of(params,m_bRun), util::IArchetypeInstance::NTPF_Ignore );

            const char *vec_names_mt[] = { "Id", "QR", "PD", "PD_R", "PD_F", "PD_P", "PD_PN", "PD_SVD", "IHFSDM", "ITF_LRM", "ECIE_CLRM", "ECIE_NHC0", "ECIE_NHC1", "PD_CLRM" };
            uint32 vec_values_mt[] = { eMethod_Id,
                                       eMethod_QR,
                                       eMethod_PD,
                                       eMethod_PD_Reflect,
                                       eMethod_PD_Fix,
                                       eMethod_PD_Project,
                                       eMethod_PD_Project_Nearest,
                                       eMethod_PD_SVD,
                                       eMethod_IHFSDM,
                                       eMethod_ITF_LRM,
                                       eMethod_ECIE_CLRM,
                                       eMethod_ECIE_NHC0,
                                       eMethod_ECIE_NHC1,
                                       eMethod_PD_CLRM };
            al.AddProperty_Enum32( "Method",
                                   (uint32)eMethod_Id,
                                   (uint32)cNumMethods, vec_names_mt, vec_values_mt,
                                   archetype_offset_of(params,m_Method),
                                   IArchetypeInstance::NTPF_Ignore );
            const char *vec_names_ncm[] = { "NDN", "ToC", "IHFSDM" };
            uint32 vec_values_ncm[] = { eNCM_NDN,
                                        eNCM_ToC,
                                        eNCM_IHFSDM };
            al.AddProperty_Enum32( "NCM",
                                   (uint32)eNCM_NDN,
                                   (uint32)cNumNCM, vec_names_ncm, vec_values_ncm,
                                   archetype_offset_of(params,m_NCM),
                                   IArchetypeInstance::NTPF_Ignore );

            al.AddProperty_NIR<float32>( "factor det(F)", 0.0f, 0.0f, 1.0f, archetype_offset_of(params,m_FactorDetF), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR<float32>( "ICF det(F)", 0.0f, 0.0f, 5.0f, archetype_offset_of(params,m_InvertedCompressibilityFactorDetF), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR<float32>( "eps.det(F)", 0.01f, 0.001f, 1.0f, archetype_offset_of(params,m_DegenerateThresholdDetF), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR<float32>( "ipol.det(F)", 0.0f, -5.0f, 0.0f, archetype_offset_of(params,m_ThresholdIpolDetF), IArchetypeInstance::NTPF_Ignore );

            al.AddProperty_NIR<float32>( "ECIE e threshold", 0.4, 0.01, 0.99, archetype_offset_of(params,m_ECIE_e_threshold), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR<float32>( "ECIE k factor", 1.0f, 0.0f, 20.0f, archetype_offset_of(params,m_ECIE_k_factor), IArchetypeInstance::NTPF_Ignore );

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
#ifdef __S2_DS_FEM_ENABLE_TEST
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
                                   archetype_offset_of(params,m_TestMode), IArchetypeInstance::NTPF_Rebuild );
            const char *vec_names_tn[] = { "All", "Boundary" };
            uint32 vec_values_tn[] = { (uint32)Params::eTN_All,
                                       (uint32)Params::eTN_Boundary };
            al.AddProperty_Enum32( "test nodes", (uint32)params.m_TestNodes,
                                   Params::cNumTN, vec_names_tn, vec_values_tn,
                                   archetype_offset_of(params,m_TestNodes), IArchetypeInstance::NTPF_Rebuild );
            al.AddProperty_NIR( "test fraction", params.m_TestFraction, 0.0f, 1.0f,
                                archetype_offset_of(params,m_TestFraction), IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
#endif
        al.BeginProperty_Group( "<Integrator>" );
        {
            const char *vec_names[] = { "SE", "IE" };
            uint32 vec_values[] = { (uint32)Params::eIT_SymplecticEuler, (uint32)Params::eIT_ImplicitEuler };
            al.AddProperty_Enum32( "integrator", (uint32)params.m_IntegratorType,
                                   2, vec_names, vec_values,
                                   archetype_offset_of(params,m_IntegratorType), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "variable_dt", params.m_bEnableVariableDT,
                            archetype_offset_of(params,m_bEnableVariableDT), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "dt", params.m_FixedDT, 0.001f, 0.03333f,
                                archetype_offset_of(params,m_FixedDT), IArchetypeInstance::NTPF_Rebuild );
            al.AddProperty_NIR( "time scale", params.m_TimeScale, 0.01f, 1.0f,
                                archetype_offset_of(params,m_TimeScale), IArchetypeInstance::NTPF_Rebuild );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<Material>" );
        {
            al.AddProperty_NIR( "young_modulus", params.m_YoungModulus, 1.0f, 10000.0f,
                                archetype_offset_of(params,m_YoungModulus), IArchetypeInstance::NTPF_Rebuild );
            al.AddProperty_NIR( "poisson_ratio", params.m_PoissonRatio, 0.0f, 0.49f,
                                archetype_offset_of(params,m_PoissonRatio), IArchetypeInstance::NTPF_Rebuild );
            al.AddProperty_NIR( "damping_ratio", params.m_DampingRatio, 0.0f, 20.0f,
                                archetype_offset_of(params,m_DampingRatio), Params::RebuildRayleighCoeffsFromFreqs );
            al.AddProperty_NIR( "air_drag", params.m_AirDragCoeff, 0.0f, 1.0f,
                                archetype_offset_of(params,m_AirDragCoeff), Params::RebuildAirDrag );
            al.AddProperty_NIR( "thickness", params.m_Thickness, 0.001f, 1.0f,
                                archetype_offset_of(params,m_Thickness), IArchetypeInstance::NTPF_Rebuild );
            al.AddProperty_NIR( "density", params.m_Density, 1.0f, 1000.0f,
                                archetype_offset_of(params,m_Density), IArchetypeInstance::NTPF_Rebuild );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<Rayleigh>" );
        {
            al.AddProperty_NIR( "freq1", params.m_RayleighFreq1, 0.1f, 10.0f,
                                archetype_offset_of(params,m_RayleighFreq1), Params::RebuildRayleighCoeffsFromFreqs );
            al.AddProperty_NIR( "freq2", params.m_RayleighFreq2, 10.0f, 20000.0f,
                                archetype_offset_of(params,m_RayleighFreq2), Params::RebuildRayleighCoeffsFromFreqs );
            al.AddProperty_NIR( "coeff_m", params.m_RayleighCoeff[0], 0.0f, 10.0f,
                                archetype_offset_of(params,m_RayleighCoeff[0]), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "coeff_k", params.m_RayleighCoeff[1], 0.0f, 1.0f,
                                archetype_offset_of(params,m_RayleighCoeff[1]), IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<Environment>" );
        {
            al.AddProperty_NIR( "gravity", params.m_Gravity[1], -100.0f, 0.0f,
                                archetype_offset_of(params,m_Gravity[1]), IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<SolverLS>" );
        {
            const char *vec_names[] = { "CG", "GS", "Jacobi", "LU" };
            uint32 vec_values[] = { (uint32)Params::eLSST_CG, (uint32)Params::eLSST_GS, (uint32)Params::eLSST_Jacobi, (uint32)Params::eLSST_LU };
            al.AddProperty_Enum32( "method", (uint32)params.m_SolverLS_Type,
                                   4, vec_names, vec_values,
                                   archetype_offset_of(params,m_SolverLS_Type), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "max_iter", params.m_SolverLS_MaxIter, uint32(1), uint32(100),
                                archetype_offset_of(params,m_SolverLS_MaxIter), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "log10(epsilon)", params.m_SolverLS_Log10_Epsilon, -15.0f, -1.0f,
                                archetype_offset_of(params,m_SolverLS_Log10_Epsilon), Params::RebuildSolverLS_Epsilon );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<ContactSolver>" );
        {
            const char *vec_names[] = { "None", "Penalty", "Hack", "HGP" };
            uint32 vec_values[] = { (uint32)Params::eCST_None, (uint32)Params::eCST_Penalty, (uint32)Params::eCST_Hack, (uint32)Params::eCST_HackGroundPlane };
            al.AddProperty_Enum32( "method", (uint32)params.m_ContactSolver_Type,
                                   Params::cNumCST, vec_names, vec_values,
                                   archetype_offset_of(params,m_ContactSolver_Type), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Depth Margin", params.m_ContactSolver_DepthMargin, 0.0f, 0.1f,
                                archetype_offset_of(params,m_ContactSolver_DepthMargin), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Penalty Ks", params.m_ContactSolver_Penalty_Ks, 0.0f, 100000.0f,
                                archetype_offset_of(params,m_ContactSolver_Penalty_Ks), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Penalty Kd", params.m_ContactSolver_Penalty_Kd, 0.0f, 1000.0f,
                                archetype_offset_of(params,m_ContactSolver_Penalty_Kd), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Relaxation Coeff", params.m_ContactSolver_Relaxation_Coeff, 0.0f, 1.0f,
                                archetype_offset_of(params,m_ContactSolver_Relaxation_Coeff), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Restitution Coeff", params.m_ContactSolver_Restitution_Coeff, 0.0f, 1.0f,
                                archetype_offset_of(params,m_ContactSolver_Restitution_Coeff), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty_NIR( "Max Depth", params.m_ContactSolver_MaxDepth, 0.0f, 10.0f,
                                archetype_offset_of(params,m_ContactSolver_MaxDepth), IArchetypeInstance::NTPF_Ignore );
            al.AddProperty( "BoE", params.m_ContactSolver_BreakOnError,
                            archetype_offset_of(params,m_ContactSolver_BreakOnError), util::IArchetypeInstance::NTPF_Ignore );
        }
        al.EndProperty_Group();
    }
    al.EndArchetype();
}
util::ArchetypeLibrary LeafDSH_FEM_Solid2D::Params::s_LeafDSH_FEM_Solid2D_Params_ArchetypeLibrary;
#endif

} } // namespace S2::ds
