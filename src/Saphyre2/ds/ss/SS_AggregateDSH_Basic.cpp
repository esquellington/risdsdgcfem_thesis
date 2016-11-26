#include "SS_AggregateDSH_Basic.h"
#include <Saphyre2/ds/dsh/IDynamicSystemHierarchy.h>
//#include <Saphyre2/ds/dsh/IConnector.h>
//#include <Saphyre2/ds/dsh/IConstraint.h>
#include <Saphyre2/ds/dsh/IGeom.h>

//#include <Geo/bp/IBroadPhase.h>
#include <Geo/bp/BasicBP.h> //BP type could be a creation param...
#include <Geo/bv/GSphere.h>
#include <Geo/bv/GAABB.h>
#include <Geo/util/Viz.h>

#include <Saphyre2/ds/DSG.h>

#ifdef __ENABLE_MULTI_SOLID3D
#  include <Saphyre2/ds/dsh/LeafDSH_Solid3D_FEM.h>
#endif

// #define __USE_OMP

namespace S2 {
namespace ds {

SS_AggregateDSH_Basic::SS_AggregateDSH_Basic()
: m_pDSH(0)
, m_pBP(0)
, m_PairsBP(100)
#ifdef __ENABLE_MULTI_SOLID3D
, m_AccTime(0)
, m_poolCCD2D(8)
#endif
{
}

SS_AggregateDSH_Basic::~SS_AggregateDSH_Basic()
{
}

bool SS_AggregateDSH_Basic::BindDSH( IDynamicSystemHierarchy* p_dsh )
{
    DS_ASSERT( p_dsh->IsAggregate() );
    m_pDSH = p_dsh;

    // Create Broad Phase
    m_pBP = new geo::bp::BasicBP( 5 ); //BUG: This leaks!

    // self-NotifyAdd all contained entities
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
        NotifyAddChild(it_dsh);
    for( ContainerGeoms::iterator it_geom=m_pDSH->GetGeomIterator(); it_geom.IsValid(); ++it_geom )
        NotifyAddGeom(it_geom);

    m_pBP->Update();

    return true;
}

/*! Predicted collisions and local constraint solving:
  - Integrate sub-dsh Vel
  - Predict sub-dsh state S'
  - Collision Detection at S'
    - Dispatch CD to proper handler
      - Sub-DSH if K2D
      - Merged-DSH SS if D2D
  - Undo prediction
  - Sub-DSH SolveConstraints at S->S' (SolveAcc/SolveVel)
  - IntegratePos
  - PostStab (SolvePos)
*/
void SS_AggregateDSH_Basic::Update( Real dt )
{
#ifdef __ENABLE_MULTI_SOLID3D
    bool bIsMultiDSH_Solid3D(true);
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
        bIsMultiDSH_Solid3D = bIsMultiDSH_Solid3D && it_dsh->GetDSHType() == eDSH_Leaf_Solid3D;
    if( bIsMultiDSH_Solid3D )
    {
        MultiDSH_Solid3D_Update(dt);
        return;
    }
#endif

    /*\todo from an SS::Update() we should be able to:
      - Integrate sub-dsh state with dt2
        - Pos/Vel independently or together
      - Predict sub-dsh state with dt2
        - Do/Undo prediction
      - Update sub-dsh GO and BV at any time
      - Perform full CD (bp,mp,np...)
        - Dispatch contacts to specific dsh (according to priorities,
          policy, etc...)
      - ...
    */

#ifdef __ENABLE_MERGE_DSH_PROTO
    /* Partitions: By now we'll fully rebuild partitions every
     * Update(), as favoured in Box2D, it's a LOT simpler than
     * persistent partitions in krm::phy and should not be noticeably
     * slower.
     */
    // Update each partition \todo Some DSH may have been destroyed
    // since last Update() that created the partitions!
    for( auto it_partition : m_vecPartitions )
        Partition_Update(*it_partition,dt);
#else
    // Update all children (and their GO and BV internally)
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
        it_dsh->Update(dt); //Step + SyncGO + RecomputeBV
    // Update all geoms
    // Update all BV NOT REQUIRED, as they are changed by their owners and Touched()
#endif
#ifdef __ENABLE_MERGE_DSH_PROTO
    //Unmerge everything
    m_vecPartitions.clear();
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
        it_dsh->SetPartition(0);
#endif

    // Perform collision detection and dispatching with new BV
    //---- CollisionDetection()
    m_pBP->Update();
    m_pBP->TestAll( m_PairsBP,
                    //( geo::bp::eTest_NonPersistent | geo::bp::eTest_Discrete
                    ( geo::bp::eTest_Persistent | geo::bp::eTest_Discrete
                      | geo::bp::eTest_Boundary | geo::bp::eTest_Internal ),
                    geo::bp::IBroadPhase::CBPF_True );

    for( geo::bp::PairContainer::iterator it_pair=m_PairsBP.Begin(); it_pair.IsValid(); ++it_pair )
    {
        DS_ASSERT( it_pair->IsValid1() && it_pair->IsValid2() );
        IEntity* pEntity1( it_pair->GetObj1<IEntity*>() );
        IEntity* pEntity2( it_pair->GetObj2<IEntity*>() );
        EEntityType et1( pEntity1->GetEntityType() );
        EEntityType et2( pEntity2->GetEntityType() );
        if( eEntity_DSH == et1 && eEntity_DSH == et2 )
        {
#ifdef __ENABLE_MERGE_DSH_PROTO
            //\todo dispatch to highest priority DSH or to merged solver...
            IDynamicSystemHierarchy* pDSH1 = static_cast<IDynamicSystemHierarchy*>(pEntity1);
            IDynamicSystemHierarchy* pDSH2 = static_cast<IDynamicSystemHierarchy*>(pEntity2);
            if( pDSH1->GetDSHType() == eDSH_Leaf_Solid3D && pDSH2->GetDSHType() == eDSH_Leaf_Solid3D )
            {
                DS_LOG_WARNING("NotifyPairBP S3D vs S3D");
                //\todo Non-persistent Merge DSH for next step
            }
#endif
        }
        else if( eEntity_DSH == et1 && eEntity_Geom == et2 )
            static_cast<IDynamicSystemHierarchy*>(pEntity1)->NotifyPairBP( it_pair->m_pProxy2, it_pair.GetPtr() );
        else if( eEntity_Geom == et1 && eEntity_DSH == et2 )
            static_cast<IDynamicSystemHierarchy*>(pEntity2)->NotifyPairBP( it_pair->m_pProxy1, it_pair.GetPtr() );
        else //both eEntity_Geom
        {
            DS_ASSERT( eEntity_Geom == et1 && eEntity_Geom == et2 );
            // Nothing to do...
        }
    }

    /* Finally, SyncGO and RecomputeBV of the whole AggregateDSH
    m_pDSH->SyncGO();
    m_pDSH->RecomputeBV();
    */
}

void SS_AggregateDSH_Basic::NotifyAddChild( IDynamicSystemHierarchy* p_child )
{
    // Create and Set DSH BV
    //\todo BV should be created in a pool!
    if( 2 == p_child->GetDimension() ) p_child->SetBV( new geo::bv::AABB2 );
    else if( 3 == p_child->GetDimension() ) p_child->SetBV( new geo::bv::AABB3 );
    else DS_ASSERT( false );

    // Add to BP
    const geo::bp::Proxy* pProxy = m_pBP->Add( reinterpret_cast<void*>( static_cast<IEntity*>(p_child) ),
                                               0,
                                               p_child->GetBV(),
                                               geo::bp::Proxy::eDefault,
                                               eSSPUF_Entity );
    //\todo IF added DSH has its own BP, we should add it as eNestedBP
    //with proper linking... this way it will be able to "go up" in
    //the BP hierarchy if required (eg: DSH Local Update)

    p_child->SetProxyBP( pProxy );
}

void SS_AggregateDSH_Basic::NotifyAddGeom( IGeom* p_geom )
{
    // Create and Set Geom BV
    //\todo BV should be created in a pool! BUG: These allocations leak!!
    if( 2 == p_geom->GetDimension() ) p_geom->SetBV( new geo::bv::AABB2 );
    else if( 3 == p_geom->GetDimension() ) p_geom->SetBV( new geo::bv::AABB3 );
    else DS_ASSERT( false );
    // Add to BP
    const geo::bp::Proxy* pProxy = m_pBP->Add( reinterpret_cast<void*>( static_cast<IEntity*>(p_geom) ),
                                               0,
                                               p_geom->GetBV(),
                                               geo::bp::Proxy::eBoundary,
                                               eSSPUF_Entity );
    p_geom->SetProxyBP( pProxy );
}


void SS_AggregateDSH_Basic::NotifyRemoveChild( IDynamicSystemHierarchy* p_child )
{
    // \todo Delete DSH anotations (BV, ProxyBP...)
}

void SS_AggregateDSH_Basic::NotifyRemoveGeom( IGeom* p_geom )
{
    // \todo Delete Geom anotations (BV, ProxyBP...)
}

void SS_AggregateDSH_Basic::DoViz( util::VizStream& vs ) const
{
    // BP structure
    //geo::VizBroadPhase( m_pBP, vs );

    // BP pairs
    for( geo::bp::PairContainer::iterator it_pair=m_PairsBP.Begin(); it_pair.IsValid(); ++it_pair )
    {
        geo::VizBoundingVolume( it_pair->m_pProxy1->m_pBV, vs );
        geo::VizBoundingVolume( it_pair->m_pProxy2->m_pBV, vs );
    }
}

#ifdef __ENABLE_MULTI_SOLID3D
inline LeafDSH_Solid3D_FEM* AsS3D( IDynamicSystemHierarchy* p_dsh ) { return static_cast<LeafDSH_Solid3D_FEM*>(p_dsh); }

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_Update( Real dt )
{
    // DS_LOG_WARNING("SS_AggregateDSH_Basic::Update_MultiDSH_Solid3D()");

    // Select smallest dt
    Real fixed_dt( mal::Infinity<Real>() );
    uint32 max_nr_iter( 0 );
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
    {
        fixed_dt = mal::Min( fixed_dt, AsS3D(it_dsh)->m_Params.m_FixedDT );
        max_nr_iter = mal::Max( max_nr_iter, AsS3D(it_dsh)->m_Params.m_SolverNR_MaxIter );
        AsS3D(it_dsh)->m_TotalTime += dt;
    }
    m_AccTime += dt;

    // Update all children (and their GO and BV internally)
    while( m_AccTime > fixed_dt )
    {
        MultiDSH_Solid3D_ApplyContacts_PositionAlteration_Smoothed();
        // NR init
        for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
            AsS3D(it_dsh)->NR_Init(dt); //Out of loop NR init + PositionAlteration
        // NR steps
        uint32 nr_iter(0);
        // bool bConverged(true);
        do
        {

#ifdef __USE_OMP
#  pragma omp parallel for
#endif
            for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
            {
                AsS3D(it_dsh)->NR_PreStep(dt);
                AsS3D(it_dsh)->NR_MidStep(dt); //includes K2D constraints
            }
            // DS_LOG_WARNING( "Should solve %d D2D contacts", m_poolCCD2D.Size() );
            MultiDSH_Solid3D_ApplyContacts_Impulse();
            // bConverged = true;

#ifdef __USE_OMP
#  pragma omp parallel for
#endif
            for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
                AsS3D(it_dsh)->NR_PostStep(dt);
        } while( /*bConverged || */ ++nr_iter < max_nr_iter );

#ifdef __USE_OMP
#  pragma omp parallel for
#endif
        for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
            AsS3D(it_dsh)->NR_End(dt); //Out of loop NR init

        m_AccTime -= fixed_dt;
    }
    // for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
    //     pS3D->EndStep(dt);


#ifdef __USE_OMP
#  pragma omp parallel for
#endif

    // Update all geoms and BV
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
    {
        // AsS3D(it_dsh)->FixedStep(fixed_dt);
        it_dsh->SyncGO();
        it_dsh->RecomputeBV();
    }

    // Perform collision detection and dispatching with new BV
    m_pBP->Update();
    m_pBP->TestAll( m_PairsBP,
                    //( geo::bp::eTest_NonPersistent | geo::bp::eTest_Discrete
                    ( geo::bp::eTest_Persistent | geo::bp::eTest_Discrete
                      | geo::bp::eTest_Boundary | geo::bp::eTest_Internal ),
                    geo::bp::IBroadPhase::CBPF_True );

    for( geo::bp::PairContainer::iterator it_pair=m_PairsBP.Begin(); it_pair.IsValid(); ++it_pair )
    {
        DS_ASSERT( it_pair->IsValid1() && it_pair->IsValid2() );
        IEntity* pEntity1( it_pair->GetObj1<IEntity*>() );
        IEntity* pEntity2( it_pair->GetObj2<IEntity*>() );
        EEntityType et1( pEntity1->GetEntityType() );
        EEntityType et2( pEntity2->GetEntityType() );
        if( eEntity_DSH == et1 && eEntity_DSH == et2 )
        {
            IDynamicSystemHierarchy* pDSH1 = static_cast<IDynamicSystemHierarchy*>(pEntity1);
            IDynamicSystemHierarchy* pDSH2 = static_cast<IDynamicSystemHierarchy*>(pEntity2);
            if( pDSH1->GetDSHType() == eDSH_Leaf_Solid3D && pDSH2->GetDSHType() == eDSH_Leaf_Solid3D )
            {
                // DS_LOG_WARNING("NotifyPairBP S3D vs S3D");
                MultiDSH_Solid3D_NotifyPairBP( it_pair.GetPtr() );
                //\todo Create D2D contact constraint so that it's handled in the next frame
            }
        }
        else if( eEntity_DSH == et1 && eEntity_Geom == et2 )
            static_cast<IDynamicSystemHierarchy*>(pEntity1)->NotifyPairBP( it_pair->m_pProxy2, it_pair.GetPtr() );
        else if( eEntity_Geom == et1 && eEntity_DSH == et2 )
            static_cast<IDynamicSystemHierarchy*>(pEntity2)->NotifyPairBP( it_pair->m_pProxy1, it_pair.GetPtr() );
        else //both eEntity_Geom
        {
            DS_ASSERT( eEntity_Geom == et1 && eEntity_Geom == et2 );
            // Nothing to do...
        }
    }

    /* Finally, SyncGO and RecomputeBV of the whole AggregateDSH
    m_pDSH->SyncGO();
    m_pDSH->RecomputeBV();
    */
}

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_NotifyPairBP( const geo::bp::Pair* p_pair )
{
    //DS_LOG_WARNING("LeafDSH_FEM_Solid3D_Simplified::NotifyPairBP");
    IEntity* pEntity1( p_pair->GetObj1<IEntity*>() );
    IEntity* pEntity2( p_pair->GetObj2<IEntity*>() );
    IDynamicSystemHierarchy* pDSH1 = static_cast<IDynamicSystemHierarchy*>(pEntity1);
    IDynamicSystemHierarchy* pDSH2 = static_cast<IDynamicSystemHierarchy*>(pEntity2);

    switch( p_pair->m_State )
    {
    case geo::bp::Pair::eNew:
        {
            //DS_LOG_WARNING("New BP");
            ContactConstraintD2D* pCCD2D( m_poolCCD2D.New() );
            pCCD2D->Create();
            p_pair->SetUserData( pCCD2D );
            if( geo::mp::TestContact( static_cast<const geo::IObject3*>(pDSH1->GetGO()),
                                      static_cast<const geo::IObject3*>(pDSH2->GetGO()),
                                      pCCD2D->m_TmpCD, &pCCD2D->m_CC ) ) //\todo Avoid casts specializing GetGO() with dimension
                pCCD2D->Reset( pDSH1, pDSH2, pCCD2D->m_TmpCD, 0.001f );
        }
        break;
    case geo::bp::Pair::ePersistent:
        {
            //DS_LOG_WARNING("Persistent BP");
            ContactConstraintD2D* pCCD2D( p_pair->GetUserData<ContactConstraintD2D*>() );
            if( geo::mp::TestContact( static_cast<const geo::IObject3*>(pDSH1->GetGO()),
                                      static_cast<const geo::IObject3*>(pDSH2->GetGO()),
                                      pCCD2D->m_TmpCD, &pCCD2D->m_CC ) ) //\todo Avoid casts specializing GetGO() with dimension
                pCCD2D->Reset( pDSH1, pDSH2, pCCD2D->m_TmpCD, 0.001f );
            else
                pCCD2D->Clear();
        }
        break;
    case geo::bp::Pair::eVanished:
        {
            //DS_LOG_WARNING("Vanished BP");
            ContactConstraintD2D* pCCD2D( p_pair->GetUserData<ContactConstraintD2D*>() );
            pCCD2D->Destroy(); //This invalidates the cache too
            m_poolCCD2D.Delete( pCCD2D );
            p_pair->SetUserData(0);
        }
        break;
    default: break;
    }
}

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_ApplyContacts_PositionAlteration_Smoothed()
{
    // Use tmp arrays to avoid new/delete
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
    {
        Vec3* vec_dp( AsS3D(it_dsh)->m_LS_vec_y );
        Real* vec_num_dp( AsS3D(it_dsh)->m_LS_real_array_tmp0 );
        memset( vec_dp, 0, sizeof(Vec3) * AsS3D(it_dsh)->m_NumNodes );
        memset( vec_num_dp, 0, sizeof(Real) * AsS3D(it_dsh)->m_NumNodes );
    }

    for( PoolCCD2D::iterator it_ccd2d=m_poolCCD2D.Begin(); it_ccd2d.IsValid(); ++it_ccd2d )
    {
        const ContactConstraintD2D& ccd2d( *it_ccd2d );
        if( ccd2d.IsEmpty() ) continue;

        LeafDSH_Solid3D_FEM* pS3D1 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH1);
        LeafDSH_Solid3D_FEM* pS3D2 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH2);

        Vec3* vec_dp1( pS3D1->m_LS_vec_y );
        Vec3* vec_dp2( pS3D2->m_LS_vec_y );
        Real* vec_num_dp1( pS3D1->m_LS_real_array_tmp0 );
        Real* vec_num_dp2( pS3D2->m_LS_real_array_tmp0 );

        for( unsigned int it_cpc=0; it_cpc<ccd2d.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2D& cpc( ccd2d.GetCPC(it_cpc) );
            if( cpc.m_Depth > pS3D1->m_Params.m_ContactSolver_DepthMargin )
            {
                // if( cpc.m_POF1.m_FeatureId.IsVertex() )
                // {
                //     geo::tetsolid3_feature_index_type vid0 = cpc.m_POF1.m_FeatureId.AsVertex();
                //     vec_dp[vid0] += m_Params.m_ContactSolver_Relaxation_Coeff * (cpc.m_Depth-m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal;
                //     vec_num_dp[vid0]++;
                // }
                // else
                if( cpc.m_POF1.m_FeatureId.IsTetrahedron() && cpc.m_POF2.m_FeatureId.IsTetrahedron() )
                {
                    geo::tetsolid3_feature_index_type eid1 = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                    geo::tetsolid3_feature_index_type eid2 = cpc.m_POF2.m_FeatureId.AsTetrahedron();
                    uint32 vec_vid1[4];
                    uint32 vec_vid2[4];
                    pS3D1->m_pMeshS->T_VecVID( eid1, vec_vid1 );
                    pS3D2->m_pMeshS->T_VecVID( eid2, vec_vid2 );
                    // Compute change at P_cp
                    //Vec3 delta_p( pS3D1->m_Params.m_ContactSolver_Relaxation_Coeff * (cpc.m_Depth-pS3D1->m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal );
                    Vec3 delta_p( (cpc.m_Depth-pS3D1->m_Params.m_ContactSolver_DepthMargin) * cpc.m_Normal ); //TEMP: apply full correction
                    // Enforce \Delta p distributing it among p1..p4 using effective mass at p
                    Real effective_mass1( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * pS3D1->m_Model.GetInvMass(vec_vid1[0])
                                                    + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * pS3D1->m_Model.GetInvMass(vec_vid1[1])
                                                    + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * pS3D1->m_Model.GetInvMass(vec_vid1[2])
                                                    + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * pS3D1->m_Model.GetInvMass(vec_vid1[3]) ) );
                    Real effective_mass2( mal::Rcp( mal::Sq( cpc.m_POF2.m_BarycentricCoords[0] ) * pS3D2->m_Model.GetInvMass(vec_vid2[0])
                                                    + mal::Sq( cpc.m_POF2.m_BarycentricCoords[1] ) * pS3D2->m_Model.GetInvMass(vec_vid2[1])
                                                    + mal::Sq( cpc.m_POF2.m_BarycentricCoords[2] ) * pS3D2->m_Model.GetInvMass(vec_vid2[2])
                                                    + mal::Sq( cpc.m_POF2.m_BarycentricCoords[3] ) * pS3D2->m_Model.GetInvMass(vec_vid2[3]) ) );
                    //\todo Apply displacements
                    Real w1( mal::Rcp(effective_mass1) / (mal::Rcp(effective_mass1) + mal::Rcp(effective_mass2)) );
                    Real w2( 1 - w1 );
                    for( int i=0; i<4; i++ )
                    {
                        vec_dp1[vec_vid1[i]] += ( w1 * cpc.m_POF1.m_BarycentricCoords[i] * effective_mass1 * pS3D1->m_Model.GetInvMass(vec_vid1[i]) ) * delta_p;
                        vec_dp2[vec_vid2[i]] -= ( w2 * cpc.m_POF2.m_BarycentricCoords[i] * effective_mass2 * pS3D2->m_Model.GetInvMass(vec_vid2[i]) ) * delta_p;
                        vec_num_dp1[vec_vid1[i]]++;
                        vec_num_dp2[vec_vid2[i]]++;
                    }
                }
                else { DS_LOG_ERROR( "Unsupported feature_id %d,%d", cpc.m_POF1.m_FeatureId.m_Type, cpc.m_POF2.m_FeatureId.m_Type ); }
            }
        }
    }
    // Apply smoothed displacements to ALL involved S3D
    // Use tmp arrays to avoid new/delete
    for( ContainerDSH::iterator it_dsh=m_pDSH->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
    {
        const Vec3* vec_dp( AsS3D(it_dsh)->m_LS_vec_y );
        const Real* vec_num_dp( AsS3D(it_dsh)->m_LS_real_array_tmp0 );
        Vec3* vec_pos( AsS3D(it_dsh)->m_Model.GetVecPos_RW() );
        for( uint32 it_node=0; it_node<AsS3D(it_dsh)->m_NumNodes; it_node++ )
            if( vec_num_dp[it_node] > 0 )
            {
                vec_pos[it_node] += vec_dp[it_node] * mal::Rcp<Real>(vec_num_dp[it_node]);
                // DS_LOG_WARNING("Displacement (%f,%f,%f)/%d",vec_dp[it_node].x(),vec_dp[it_node].y(),vec_dp[it_node].z(),(int)vec_num_dp[it_node]);
            }
    }
}

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_ApplyContacts_Impulse()
{
    MultiDSH_Solid3D_ApplyContacts_UpdateActiveSet();
    MultiDSH_Solid3D_ApplyContacts_NormalImpulse();
    MultiDSH_Solid3D_ApplyContacts_FrictionImpulse();
}

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_ApplyContacts_UpdateActiveSet()
{
    for( PoolCCD2D::iterator it_ccd2d=m_poolCCD2D.Begin(); it_ccd2d.IsValid(); ++it_ccd2d )
    {
        ContactConstraintD2D& ccd2d( *it_ccd2d );
        if( ccd2d.IsEmpty() ) continue;
        LeafDSH_Solid3D_FEM* pS3D1 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH1);
        LeafDSH_Solid3D_FEM* pS3D2 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH2);
        const Vec3* vec_vel1( pS3D1->m_LS_vec_y );
        const Vec3* vec_vel2( pS3D2->m_LS_vec_y );
        for( unsigned int it_cpc=0; it_cpc<ccd2d.GetNumCPC(); it_cpc++ )
        {
            ContactPointConstraintD2D& cpc( ccd2d.GetCPC(it_cpc) );
            // if( cpc.m_POF1.m_FeatureId.IsVertex() )
            // {
            //     geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
            //     cpc.m_IsActive = mal::Dot(vec_vel[vid],cpc.m_Normal) < 0;     //Collision
            // }
            // else
            if( cpc.m_POF1.m_FeatureId.IsTetrahedron() && cpc.m_POF2.m_FeatureId.IsTetrahedron() )
            {
                geo::tetsolid3_feature_index_type eid1 = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                geo::tetsolid3_feature_index_type eid2 = cpc.m_POF2.m_FeatureId.AsTetrahedron();
                uint32 vec_vid1[4];
                uint32 vec_vid2[4];
                pS3D1->m_pMeshS->T_VecVID( eid1, vec_vid1 );
                pS3D2->m_pMeshS->T_VecVID( eid2, vec_vid2 );
                // Interpolate V_cp, A_cp from node v,a
                Vec3 v1( cpc.m_POF1.m_BarycentricCoords[0] * vec_vel1[vec_vid1[0]]
                         + cpc.m_POF1.m_BarycentricCoords[1] * vec_vel1[vec_vid1[1]]
                         + cpc.m_POF1.m_BarycentricCoords[2] * vec_vel1[vec_vid1[2]]
                         + cpc.m_POF1.m_BarycentricCoords[3] * vec_vel1[vec_vid1[3]] );
                Vec3 v2( cpc.m_POF2.m_BarycentricCoords[0] * vec_vel2[vec_vid2[0]]
                         + cpc.m_POF2.m_BarycentricCoords[1] * vec_vel2[vec_vid2[1]]
                         + cpc.m_POF2.m_BarycentricCoords[2] * vec_vel2[vec_vid2[2]]
                         + cpc.m_POF2.m_BarycentricCoords[3] * vec_vel2[vec_vid2[3]] );
                cpc.m_IsActive = mal::Dot(v1-v2,cpc.m_Normal) < 0; //\todo
            }
            else { DS_LOG_ERROR( "Unsupported feature_id %d,%d", cpc.m_POF1.m_FeatureId.m_Type, cpc.m_POF2.m_FeatureId.m_Type ); }
        }
    }
}

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_ApplyContacts_NormalImpulse()
{
    for( PoolCCD2D::iterator it_ccd2d=m_poolCCD2D.Begin(); it_ccd2d.IsValid(); ++it_ccd2d )
    {
        ContactConstraintD2D &ccd2d( *it_ccd2d );
        if( ccd2d.IsEmpty() ) continue;
        LeafDSH_Solid3D_FEM* pS3D1 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH1);
        LeafDSH_Solid3D_FEM* pS3D2 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH2);
        Vec3* vec_v1( pS3D1->m_LS_vec_y );
        Vec3* vec_v2( pS3D2->m_LS_vec_y );

        for( unsigned int it_cpc=0; it_cpc<ccd2d.GetNumCPC(); it_cpc++ )
        {
            ContactPointConstraintD2D &cpc( ccd2d.GetCPC(it_cpc) );
            //if( cpc.m_Depth > m_Params.m_ContactSolver_DepthMargin )
            {
                if( cpc.m_IsActive ) //active = collision || contact = v*n < 0 || \Delta v*n < 0
                {
                    // if( cpc.m_POF1.m_FeatureId.IsVertex() )
                    // {
                    //     geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
                    //     cpc.m_LambdaN = -mal::Dot(vec_v[vid],cpc.m_Normal);
                    //     vec_v[vid] += cpc.m_LambdaN * cpc.m_Normal;
                    //     cpc.m_LambdaN *= factor_v_to_impulse; //Apply factor to convert input magnitude v to impulse
                    // }
                    // else
                    if( cpc.m_POF1.m_FeatureId.IsTetrahedron() && cpc.m_POF2.m_FeatureId.IsTetrahedron() )
                    {
                        geo::tetsolid3_feature_index_type eid1 = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                        geo::tetsolid3_feature_index_type eid2 = cpc.m_POF2.m_FeatureId.AsTetrahedron();
                        uint32 vec_vid1[4];
                        uint32 vec_vid2[4];
                        pS3D1->m_pMeshS->T_VecVID( eid1, vec_vid1 );
                        pS3D2->m_pMeshS->T_VecVID( eid2, vec_vid2 );
                        // Interpolate V_cp from node velocities
                        Vec3 v1( cpc.m_POF1.m_BarycentricCoords[0] * vec_v1[vec_vid1[0]]
                                 + cpc.m_POF1.m_BarycentricCoords[1] * vec_v1[vec_vid1[1]]
                                 + cpc.m_POF1.m_BarycentricCoords[2] * vec_v1[vec_vid1[2]]
                                 + cpc.m_POF1.m_BarycentricCoords[3] * vec_v1[vec_vid1[3]] );
                        Vec3 v2( cpc.m_POF2.m_BarycentricCoords[0] * vec_v2[vec_vid2[0]]
                                 + cpc.m_POF2.m_BarycentricCoords[1] * vec_v2[vec_vid2[1]]
                                 + cpc.m_POF2.m_BarycentricCoords[2] * vec_v2[vec_vid2[2]]
                                 + cpc.m_POF2.m_BarycentricCoords[3] * vec_v2[vec_vid2[3]] );
                        // Compute change at V_cp
                        cpc.m_LambdaN = -mal::Dot(v1-v2,cpc.m_Normal); //\todo
                        Vec3 delta_v( cpc.m_LambdaN * cpc.m_Normal );
                        // Enforce \Delta v distributing it among v1..v4 using effective mass at v
                        Real effective_mass1( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * pS3D1->m_Model.GetInvMass(vec_vid1[0])
                                                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * pS3D1->m_Model.GetInvMass(vec_vid1[1])
                                                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * pS3D1->m_Model.GetInvMass(vec_vid1[2])
                                                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * pS3D1->m_Model.GetInvMass(vec_vid1[3]) ) );
                        Real effective_mass2( mal::Rcp( mal::Sq( cpc.m_POF2.m_BarycentricCoords[0] ) * pS3D2->m_Model.GetInvMass(vec_vid2[0])
                                                        + mal::Sq( cpc.m_POF2.m_BarycentricCoords[1] ) * pS3D2->m_Model.GetInvMass(vec_vid2[1])
                                                        + mal::Sq( cpc.m_POF2.m_BarycentricCoords[2] ) * pS3D2->m_Model.GetInvMass(vec_vid2[2])
                                                        + mal::Sq( cpc.m_POF2.m_BarycentricCoords[3] ) * pS3D2->m_Model.GetInvMass(vec_vid2[3]) ) );
                        Real w1( mal::Rcp(effective_mass1) / (mal::Rcp(effective_mass1) + mal::Rcp(effective_mass2)) );
                        Real w2( 1 - w1 );
                        for( int i=0; i<4; i++ )
                        {
                            vec_v1[vec_vid1[i]] += ( w1 * cpc.m_POF1.m_BarycentricCoords[i] * effective_mass1 * pS3D1->m_Model.GetInvMass(vec_vid1[i]) ) * delta_v;
                            vec_v2[vec_vid2[i]] -= ( w2 * cpc.m_POF2.m_BarycentricCoords[i] * effective_mass2 * pS3D2->m_Model.GetInvMass(vec_vid2[i]) ) * delta_v;
                        }
                        cpc.m_LambdaN *= 1;//\todo proper factor from Step_FullyImplicitEuler_dx factor_v_to_impulse; //Apply factor to convert input magnitude v to impulse
                    }
                    else { DS_LOG_ERROR( "Unsupported feature_id %d,%d", cpc.m_POF1.m_FeatureId.m_Type, cpc.m_POF2.m_FeatureId.m_Type ); }
                }
            }
        }
    }
}

void SS_AggregateDSH_Basic::MultiDSH_Solid3D_ApplyContacts_FrictionImpulse()
{
    for( PoolCCD2D::iterator it_ccd2d=m_poolCCD2D.Begin(); it_ccd2d.IsValid(); ++it_ccd2d )
    {
        const ContactConstraintD2D &ccd2d( *it_ccd2d );
        if( ccd2d.IsEmpty() ) continue;
        LeafDSH_Solid3D_FEM* pS3D1 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH1);
        LeafDSH_Solid3D_FEM* pS3D2 = static_cast<LeafDSH_Solid3D_FEM*>(ccd2d.m_pDSH2);
        Vec3* vec_v1( pS3D1->m_LS_vec_y );
        Vec3* vec_v2( pS3D2->m_LS_vec_y );

        for( unsigned int it_cpc=0; it_cpc<ccd2d.GetNumCPC(); it_cpc++ )
        {
            const ContactPointConstraintD2D &cpc( ccd2d.GetCPC(it_cpc) );
            // if( cpc.m_POF1.m_FeatureId.IsVertex() )
            // {
            //     geo::tetsolid3_feature_index_type vid = cpc.m_POF1.m_FeatureId.AsVertex();
            //     Vec3 vt( vec_v[vid] - mal::Dot(vec_v[vid],cpc.m_Normal)*cpc.m_Normal );
            //     Real norm_sq_vt( mal::NormSq(vt) );
            //     if( norm_sq_vt > 1e-8 )
            //     {
            //         // Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN * m_Model.GetInvMass(vid) ); //Apply factor to convert cached impulse to v magnitude
            //         Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN ); //Apply factor to convert cached impulse to v magnitude
            //         if( delta_vn > 0 )
            //         {
            //             Real norm_vt( mal::Sqrt(norm_sq_vt) );
            //             Real post_vt_over_norm_vt( mal::Max(0.0f,1.0f-m_Params.m_ContactSolver_DynamicFriction_Coeff*(delta_vn/norm_vt)) );
            //             // Final v includes tangential friction but keeps normal component unchanged
            //             vec_v[vid] = post_vt_over_norm_vt * vt + mal::Dot(vec_v[vid],cpc.m_Normal)*cpc.m_Normal;
            //         }
            //     }
            // }
            // else
            if( cpc.m_POF1.m_FeatureId.IsTetrahedron() && cpc.m_POF2.m_FeatureId.IsTetrahedron() )
            {
                geo::tetsolid3_feature_index_type eid1 = cpc.m_POF1.m_FeatureId.AsTetrahedron();
                geo::tetsolid3_feature_index_type eid2 = cpc.m_POF2.m_FeatureId.AsTetrahedron();
                uint32 vec_vid1[4];
                uint32 vec_vid2[4];
                pS3D1->m_pMeshS->T_VecVID( eid1, vec_vid1 );
                pS3D2->m_pMeshS->T_VecVID( eid2, vec_vid2 );
                // Interpolate V_cp from node velocities
                Vec3 v1( cpc.m_POF1.m_BarycentricCoords[0] * vec_v1[vec_vid1[0]]
                         + cpc.m_POF1.m_BarycentricCoords[1] * vec_v1[vec_vid1[1]]
                         + cpc.m_POF1.m_BarycentricCoords[2] * vec_v1[vec_vid1[2]]
                         + cpc.m_POF1.m_BarycentricCoords[3] * vec_v1[vec_vid1[3]] );
                Vec3 v2( cpc.m_POF2.m_BarycentricCoords[0] * vec_v2[vec_vid2[0]]
                         + cpc.m_POF2.m_BarycentricCoords[1] * vec_v2[vec_vid2[1]]
                         + cpc.m_POF2.m_BarycentricCoords[2] * vec_v2[vec_vid2[2]]
                         + cpc.m_POF2.m_BarycentricCoords[3] * vec_v2[vec_vid2[3]] );
                Vec3 vt( v1 - v2 - mal::Dot(v1-v2,cpc.m_Normal)*cpc.m_Normal );
                Real norm_sq_vt( mal::NormSq(vt) );
                if( norm_sq_vt > 1e-8 )
                {
                    // Real delta_vn( factor_impulse_to_v * cpc.m_LambdaN / m ); //Apply factor to convert cached impulse to v magnitude \todo NOT SURE about barycentric mass...
                    Real delta_vn( /*factor_impulse_to_v * TEMP */ cpc.m_LambdaN );
                    if( delta_vn > 0 )
                    {
                        Real norm_vt( mal::Sqrt(norm_sq_vt) );
                        Real post_vt_over_norm_vt( mal::Max(0.0f,1.0f-pS3D1->m_Params.m_ContactSolver_DynamicFriction_Coeff*(delta_vn/norm_vt)) );
                        // Final v includes tangential friction but keeps normal component unchanged
                        Vec3 v_post( post_vt_over_norm_vt * vt + mal::Dot(v1-v2,cpc.m_Normal)*cpc.m_Normal );
                        Vec3 delta_v( v_post - v1 + v2 );
                        // Enforce \Delta v distributing it among v1..v4 using effective mass at v
                        Real effective_mass1( mal::Rcp( mal::Sq( cpc.m_POF1.m_BarycentricCoords[0] ) * pS3D1->m_Model.GetInvMass(vec_vid1[0])
                                                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[1] ) * pS3D1->m_Model.GetInvMass(vec_vid1[1])
                                                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[2] ) * pS3D1->m_Model.GetInvMass(vec_vid1[2])
                                                        + mal::Sq( cpc.m_POF1.m_BarycentricCoords[3] ) * pS3D1->m_Model.GetInvMass(vec_vid1[3]) ) );
                        Real effective_mass2( mal::Rcp( mal::Sq( cpc.m_POF2.m_BarycentricCoords[0] ) * pS3D2->m_Model.GetInvMass(vec_vid2[0])
                                                        + mal::Sq( cpc.m_POF2.m_BarycentricCoords[1] ) * pS3D2->m_Model.GetInvMass(vec_vid2[1])
                                                        + mal::Sq( cpc.m_POF2.m_BarycentricCoords[2] ) * pS3D2->m_Model.GetInvMass(vec_vid2[2])
                                                        + mal::Sq( cpc.m_POF2.m_BarycentricCoords[3] ) * pS3D2->m_Model.GetInvMass(vec_vid2[3]) ) );
                        Real w1( mal::Rcp(effective_mass1) / (mal::Rcp(effective_mass1) + mal::Rcp(effective_mass2)) );
                        Real w2( 1 - w1 );
                        for( int i=0; i<4; i++ )
                        {
                            vec_v1[vec_vid1[i]] += ( w1 * cpc.m_POF1.m_BarycentricCoords[i] * effective_mass1 * pS3D1->m_Model.GetInvMass(vec_vid1[i]) ) * delta_v;
                            vec_v2[vec_vid2[i]] -= ( w2 * cpc.m_POF2.m_BarycentricCoords[i] * effective_mass2 * pS3D2->m_Model.GetInvMass(vec_vid2[i]) ) * delta_v;
                        }
                    }
                }
            }
            else { DS_LOG_ERROR( "Unsupported feature_id %d", cpc.m_POF1.m_FeatureId.m_Type ); }
        }
    }
}

#endif

} } // namespace S2::ds
