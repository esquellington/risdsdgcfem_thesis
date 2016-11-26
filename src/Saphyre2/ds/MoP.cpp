#include <Saphyre2/ds/MoP.h>
#include <Saphyre2/ds/DSG.h>

#include <Saphyre2/ds/dsh/IDynamicSystemHierarchy.h>
#include <Saphyre2/ds/dsh/IGeom.h>
#include <Saphyre2/ds/dsh/IConnector.h>
#include <Saphyre2/ds/dsh/IConstraint.h>

#include <Saphyre2/ds/dsh/LeafDSH_Particle.h>
#include <Saphyre2/ds/dsh/LeafDSH_SPH_Fluid2D.h>
#include <Saphyre2/ds/dsh/LeafDSH_ParticleSystem.h>

//TEMP: by now all versions must work, but Solid2D_FEM will become the default one
//#include <Saphyre2/ds/dsh/LeafDSH_FEM_Solid2D.h> //Kept for DCFEM
//#include <Saphyre2/ds/dsh/LeafDSH_FEM_Solid2D_Simplified.h> //Simplfied version
#include <Saphyre2/ds/dsh/LeafDSH_Solid2D_FEM.h> //Generic version
#include <Saphyre2/ds/dsh/LeafDSH_Solid3D_FEM.h>

#include <Saphyre2/ds/dsh/GeomFactory.h>

#include <Saphyre2/ds/ss/SS_AggregateDSH_Basic.h>
#include <Saphyre2/ds/ss/SS_LeafDSH_ParticleSystem.h>

#include <Geo/shape/GSphereSetShape.h> //TEMPORAL
#include <Geo/mp/TestRay.h> //TEMPORAL
#include <Geo/bv/TestRay.h> //TEMPORAL

namespace S2 {
namespace ds {

//---- Some pointer hacks to cast EID to Entity and DSH
inline IEntity* Cast2Entity( machine_uint_type eid )
{
    return reinterpret_cast<IEntity*>(eid);
}
inline IDynamicSystemHierarchy* Cast2DSH( machine_uint_type eid )
{
    IEntity *pEntity( reinterpret_cast<IEntity*>(eid) );
    DS_ASSERT( pEntity->IsDSH() );
    return static_cast<IDynamicSystemHierarchy*>( pEntity );
}

//---- Specific MoP DSH subclasses
/*! MultiverseDSH is NOT, in essence, a DSH, just a top-level
    container of UniverseDSH3.
    \todo Consider a flat list, instead
*/
class MultiverseDSH: public AggregateDSH3
{
public:
    MultiverseDSH() : AggregateDSH3(0,0) {}
    EDSHType GetDSHType() const { return eDSH_Multiverse; }
    bool Create( const ParamIt &pit ) { return true; }
    bool Edit( const ParamIt &pit ) { return false; }
};

/*! UniverseDSH3 is an Aggregate with a hardcoded SS by now
 */
class UniverseDSH3: public AggregateDSH3
{
public:
    UniverseDSH3( machine_uint_type uid, IDynamicSystemHierarchy *p_parent )
    : AggregateDSH3(uid,p_parent)
    {
        BindSS(&m_SS);
    }
    EDSHType GetDSHType() const { return eDSH_Universe3D; }
    bool Create( const ParamIt &pit ) { return true; }
    bool Edit( const ParamIt &pit ) { return false; }

private:
    SS_AggregateDSH_Basic m_SS;
};

//---- Construction/Destruction
MoP::MoP()
: m_pRootDSH(0)
, m_pGeomFactory(0)
{
}

MoP::~MoP()
{
    ShutDown();
}

void MoP::Init()
{
    m_pRootDSH = new MultiverseDSH();
    m_pGeomFactory = new GeomFactory();
}

void MoP::ShutDown()
{
    if( m_pRootDSH ) delete m_pRootDSH;
    m_pRootDSH = 0;
    if( m_pGeomFactory ) delete m_pGeomFactory;
    m_pGeomFactory = 0;
}

//---- Public API for IDarkSide
bool MoP::Execute( const CommandStream &cmds, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Execute()...");
    bool bResult = true;
    for( CommandIt cmdit=cmds.Begin();
         bResult && cmdit.IsValid();
         ++cmdit )
        bResult = bResult && Execute( cmdit, rets );
    DS_ECMD(bResult);
    return bResult;
}


//---- Debug API for IDarkSide
void MoP::UpdateLog()
{
    //DSG::GetLogIS()->Clear();
    //m_pRootDSH->DoLog( *DSG::GetLogIS() );
}

void MoP::UpdateViz()
{
    DSG::GetVizIS()->Clear();
    m_pRootDSH->DoViz( *DSG::GetVizIS() );
}

void MoP::UpdateProf()
{
    //Just dump and clear existing Profiler, no need to call entities
    DSG::GetProfIS()->Clear();
    DSG::GetProfIS()->WriteArray( "ProfilerData",
                                  &DSG::GetProfiler().GetEntries()[0],
                                  DSG::GetProfiler().GetEntries().size() );
    DSG::GetProfiler().Clear();
}

//---- Internal methods
bool MoP::Execute( const CommandIt &cmdit, ReturnStream &rets )
{
    DS_ASSERT( cmdit.IsComplex() );
    switch( cmdit.GetType() )
    {
        /* Control commands
    case eCtrl_Log: return Ctrl_Log( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eCtrl_Viz: return Ctrl_Viz( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eCtrl_Prof: return Ctrl_Prof( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
        */

    // Entity commands
    case eCmd_Create: return Cmd_Create( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eCmd_Destroy: return Cmd_Destroy( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eCmd_Edit: return Cmd_Edit( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eCmd_Update: return Cmd_Update( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eCmd_Internal: return Cmd_Internal( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;

    // Entity queries
    case eQuery_RayCast: return Query_RayCast( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    //\todo case eQuery_Intersection: return Query_Intersection( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;

    // Entity debug commands
    case eDbg_Conf: return Dbg_Conf( cmdit.GetId(), cmdit.GetSubItem(), rets ); break;
    case eDbg_QueryStats: return Dbg_QueryStats( cmdit.GetId(), cmdit.GetSubItem() ); break;
    case eDbg_QueryParams: return Dbg_QueryParams( cmdit.GetId(), cmdit.GetSubItem() ); break;
    case eDbg_SyncParams: return Dbg_SyncParams( cmdit.GetId(), cmdit.GetSubItem() ); break;

    default: DS_INFO("Unknown Cmd!"); return false; break;
    }
    return false;
}

//---- Control commands
/*
bool MoP::Ctrl_Log( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Ctrl_Log()...");
    DS_ECMD(false);
    return false;
}

bool MoP::Ctrl_Viz( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Ctrl_Viz()...");
    DS_ECMD(false);
    return false;
}
bool MoP::Ctrl_Prof( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Ctrl_Prof()...");
    DS_ECMD(false);
    return false;
}
*/

//---- Entity commands

/*! Create a DS Entity:
  - DSH
    - Universe
    - Aggregate
    - Leaf
  - Connector
  - Constraint
  - Geom
*/
bool MoP::Cmd_Create( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Cmd_Create()...");
    bool bResult(false);

    //---- Try to create dsh entity
    IEntity *pEntity = 0;
    IDynamicSystemHierarchy *pParentDSH = 0;
    EEntityType entitytype = (EEntityType)pit.Find("entity_type").Get<uint32>();
    switch( entitytype )
    {
    case eEntity_DSH:
        {
            //---- Instantiate dsh entity
            IDynamicSystemHierarchy *pDSH = 0;
            EDSHType dsh_type = (EDSHType)pit.Find("dsh_type").Get<uint32>();
            switch( dsh_type )
            {
            //-- Aggregates
            case eDSH_Universe3D:
                {
                    DS_INFO("Creating eDSH_Universe3D...");
                    pParentDSH = m_pRootDSH; //!< Parent of Universe is Root Multiverse
                    pDSH = new UniverseDSH3(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = 0 != pDSH;
                }
                break;
            //-- 3D Models
            case eDSH_Leaf_Particle3D:
                {
                    DS_INFO("Creating eDSH_Leaf_Particle3D...");
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );
                    pDSH = new LeafDSH_Particle3D(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = 0 != pDSH;
                }
                break;
            //-- 2D Models
            case eDSH_Leaf_Particle2D:
                {
                    DS_INFO("Creating eDSH_Leaf_Particle2D...");
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );
                    pDSH = new LeafDSH_Particle2D(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = 0 != pDSH;
                }
                break;
            case eDSH_Leaf_ParticleSystem2D:
                {
                    DS_INFO("Creating eDSH_Leaf_ParticleSystem2D...");
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );
                    pDSH = new LeafDSH_ParticleSystem2D(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = pDSH->Create(pit);
                    pDSH->BindSS( new SS_LeafDSH_ParticleSystem() );
                    //TEMPORAL HACK: Should bind user-specified GO like with other Models...
                    {
                        geo::GObjectSS<geo::SphereSetShape2> *pGO( new geo::GObjectSS<geo::SphereSetShape2>() );
                        geo::SphereSetShape2 *pSS2( new geo::SphereSetShape2( pit.Find("count").Get<uint32>() ) );
                        pSS2->SetRadius( pit.Find("radius").Get<float32>() );
                        pGO->SetShape( pSS2 );
                    }
                    bResult = 0 != pDSH;
                }
                break;
            case eDSH_Leaf_Fluid2D:
                {
                    DS_INFO("Creating eDSH_Leaf_Fluid2D...");
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );
                    pDSH = new LeafDSH_SPH_Fluid2D(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = 0 != pDSH;
                }
                break;
            case eDSH_Leaf_Solid2D:
                {
                    DS_INFO("Creating eDSH_Leaf_Solid2D...");
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );
                    //TEMP: by now all versions must work, but Solid2D_FEM will become the default one
                    //pDSH = new LeafDSH_FEM_Solid2D(pit.Find("uid").Get<machine_uint_type>(),pParentDSH); //Kept for DCFEM
                    //pDSH = new LeafDSH_FEM_Solid2D_Simplified(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    pDSH = new LeafDSH_Solid2D_FEM(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = 0 != pDSH;
                }
                break;
            case eDSH_Leaf_Solid3D:
                {
                    DS_INFO("Creating eDSH_Leaf_Solid3D...");
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );
                    pDSH = new LeafDSH_Solid3D_FEM(pit.Find("uid").Get<machine_uint_type>(),pParentDSH);
                    bResult = 0 != pDSH;
                }
                break;
            default: break;
            }

            //---- BindGO if specified
            if( bResult )
            {
                // Create and bind geometry, if any
                ParamIt shape_id( pit.Find("shape_id") );
                if( shape_id.IsValid() )
                {
                    //\todo This comment is out of date?? DS_LOG_ERROR( "We DO NOT support shape_id BS->DS by now..." );
                    pDSH->BindGO( m_pGeomFactory->GetGeoObjFactory().CreateSS( shape_id.Get<uint32>() ) );
                }
                else
                {
                    ParamIt shape_def = pit.Find("shape_def");
                    if( shape_def.IsValid() )
                        pDSH->BindGO( m_pGeomFactory->GetGeoObjFactory().CreateES( shape_def ) );
                }

                // Link to parent
                pEntity = pDSH;
                if( pParentDSH ) pParentDSH->AddChild( pDSH );
            }
            else if( pDSH )
                delete pDSH;

            //---- Create the dsh entity (\note CAN use the already bound GO)
            if( bResult )
                bResult = pDSH->Create(pit);
        }
        break;
    case eEntity_Connector: break;
    case eEntity_Constraint: break;
    case eEntity_Interaction: break;
    case eEntity_Geom:
        {
            //---- Try to create geom entity
            IGeom *pGeom = 0;
            EGeomType geom_type = (EGeomType)pit.Find("geom_type").Get<uint32>();
            switch( geom_type )
            {
            case eGeom_Simple:
                {
                    pParentDSH = Cast2DSH( pit.Find("peid").Get<machine_uint_type>() );

                    // Get shape params
                    ParamIt shape_id( pit.Find("shape_id") );
                    if( shape_id.IsValid() )
                    {
                        DS_LOG_ERROR( "We DO NOT support shape_id BS->DS by now..." );
                        DS_INFO( "shape_id = " << shape_id.Get<uint32>() );
                        pGeom = m_pGeomFactory->Create( pit.Find("uid").Get<machine_uint_type>(), pParentDSH, shape_id.Get<uint32>() );
                    }
                    else
                    {
                        ParamIt shape_def = pit.Find("shape_def");
                        DS_INFO( "shape_def = " << shape_def );
                        pGeom = m_pGeomFactory->Create( pit.Find("uid").Get<machine_uint_type>(), pParentDSH, shape_def );
                    }

                    //bind transform and DOF functions to DS or Static according to binding type:
                    // DOF,
                    bResult = pGeom->Create(pit);
                }
                break;
            default: bResult = false; break;
            }

            //---- Analyze result and link/delete and return data accordingly
            if( bResult )
            {
                // Link to parent
                pEntity = pGeom;
                if( pParentDSH ) pParentDSH->AddGeom( pGeom );
            }
            else if( pGeom )
                delete pGeom;
        }
        break;
    default: break;
    }

    //---- Analyze result and link/delete and return data accordingly
    if( bResult )
    {
        // Return new entity (\todo Return state?)
        rets.BeginComplex( cid, (uint32)eRet_NewEntity );
        {
            rets.Write("uid",pEntity->GetUID());
            rets.Write("eid",(machine_uint_type)pEntity);
            rets.Write("entity_type",(uint32)entitytype); //TEMP unnecessary?
        }
        rets.EndComplex();
    }
    else
    {
        // Return and delete
        rets.Write( cid, (uint32)eRet_Error );
    }

    DS_ECMD(bResult);
    return bResult;
}

bool MoP::Cmd_Destroy( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Cmd_Destroy()...");
      IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());
      bool bResult = pEntity->Destroy(rets);
    DS_ECMD(bResult);
    return bResult;
}

bool MoP::Cmd_Edit( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Cmd_Edit()...");
      IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());
      bool bResult = pEntity->Edit(pit);
      if( !bResult ) rets.Write( cid, (uint32)eRet_Error );
    DS_ECMD(bResult);
    return bResult;
}

bool MoP::Cmd_Update( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Cmd_Update()...");
    IDynamicSystemHierarchy *pDSH = Cast2DSH(pit.Find("eid").Get<machine_uint_type>());
    // Update
    pDSH->Update( pit.Find("dt").Get<Real>() );
    // Return new state
    pDSH->QueryState(rets);
    DS_ECMD(true);
    return true;
}

bool MoP::Cmd_Internal( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Cmd_Internal()...");
      IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());
      bool bResult = pEntity->Internal(pit,rets);
      if( !bResult ) rets.Write( cid, (uint32)eRet_Error );
    DS_ECMD(bResult);
    return bResult;
}

void Write_RayHitReport2( machine_uint_type uid, const geo::np::RayHit2 &rh, ReturnStream &rets )
{
    rets.BeginComplex( "ray_hit_report", eType_Unknown ); //TEMP: Unknown should be "Generic", as in "valid-but-not-predefined"
    {
        rets.Write( "uid", uid );
        rets.Write( "interval", rh.m_Interval );
        rets.Write( "point", rh.m_Point );
        rets.Write( "normal", rh.m_Normal );
        rets.Write( "feature_type", (uint32)rh.m_FeatureId.m_Type );
        rets.Write( "feature_index", (geo::feature_index_type)rh.m_FeatureId.m_Index );
        rets.Write( "extra_bc", rh.m_Extra_BarycentricCoords );
    }
    rets.EndComplex();
}

void Write_RayHitReport3( machine_uint_type uid, const geo::np::RayHit3 &rh, ReturnStream &rets )
{
    rets.BeginComplex( "ray_hit_report", eType_Unknown ); //TEMP: Unknown should be "Generic", as in "valid-but-not-predefined"
    {
        rets.Write( "uid", uid );
        rets.Write( "interval", rh.m_Interval );
        rets.Write( "point", rh.m_Point );
        rets.Write( "normal", rh.m_Normal );
        rets.Write( "feature_type", (uint32)rh.m_FeatureId.m_Type );
        rets.Write( "feature_index", (geo::feature_index_type)rh.m_FeatureId.m_Index );
        rets.Write( "extra_bc", rh.m_Extra_BarycentricCoords );
    }
    rets.EndComplex();
}

/*\todo Raycast requires a BP? if so, it should be done by the
  SimulationScheme, not by the DSH. Or, Maybe, it could be done
  outside, using DSH/SS but not as a method, it's just a QUERY,
  does NOT change their state... maybe each LeafDSH/Geom will
  require an object-centric QueryRayCast(), but hirearchical
  query with optional BP is pretty generic.

  When querying nested BP, possible EXTERNAL proxies, as they are
  already tested in a higher level or must be ignored, if raycast
  is local!.

  If DSH or SS has BP
  - QueryRayCast on BP => Gather proxies => Get LeafDSH & Geom from proxies and test explicitly
  else
  - Recursively QueryRayCast on each sub-DSH => Gather LeafDSH & Geom and test explicitly
*/
unsigned int Query_RayCast2_DSH( const IDynamicSystemHierarchy *p_dsh, const geo::np::Ray2 &ray, ReturnStream &rets )
{
    geo::np::RayHit2 rh;
    unsigned int num_hits(0);
    if( 0 != p_dsh->GetSS() && 0 != p_dsh->GetSS()->GetBP() )
    {
        // Raycast BP to gather hit proxies
        std::vector<geo::bp::Proxy*> vec_hit_proxies;
        if( p_dsh->GetSS()->GetBP()->TestRay2( ray, vec_hit_proxies, geo::bp::eTest_Internal | geo::bp::eTest_Boundary ) > 0 )
        {
            //DS_LOG_WARNING( "Query_RayCast2_DSH: Using DSH->SS->BP! %d hits", vec_hit_proxies.size() );
            for( unsigned int it_proxy=0; it_proxy < vec_hit_proxies.size(); it_proxy++ )
            {
                if( vec_hit_proxies[it_proxy]->IsNestedBP() )
                {
                    //\todo Recursive call?
                    DS_LOG_WARNING( "Query_RayCast2_DSH: Nested BP not yet implemented" );
                }
                else if( vec_hit_proxies[it_proxy]->m_UserFlags.Test(ISimulationScheme::eSSPUF_Entity) )
                {
                    // Retrieve IEntity from proxy (safe as eSSPUF_Entity)
                    const IEntity *pEntity( vec_hit_proxies[it_proxy]->GetObj<const IEntity*>() );
                    // Retrieve GO and raycast it
                    const geo::IObject *pGO(0);
                    if( pEntity->IsGeom() ) pGO = static_cast<const IGeom*>(pEntity)->GetGO();
                    else if( pEntity->IsDSH() ) pGO = static_cast<const IDynamicSystemHierarchy*>(pEntity)->GetGO();
                    if( 0 != pGO
                        && pGO->GetDimension() == 2
                        && geo::mp::TestRay( static_cast<const geo::IObject2*>(pGO), ray, rh, 0 ) )
                    {
                        num_hits++;
                        Write_RayHitReport2( pEntity->GetUID(), rh, rets );
                    }
                    /*TEMPORAL
                    if( 0 == pGO ) DS_LOG_WARNING( "Query_RayCast2_DSH: DSH without a GO found" );
                    */
                    // Recursive call if DSH
                    if( pEntity->IsDSH() )
                        num_hits += Query_RayCast2_DSH( static_cast<const IDynamicSystemHierarchy*>(pEntity), ray, rets );
                }
                else
                {
                    DS_LOG_WARNING( "Query_RayCast2_DSH: Non eSSPUF_Entity BP proxy ignored" );
                }
            }
        }
        /*TEMPORAL
        else
            DS_LOG_WARNING( "Query_RayCast2_DSH: Using DSH->SS->BP! 0 hits" );
        */
    }
    else
    {
        // Iterate over Geoms, filter with optional BV
        for( ContainerGeoms::iterator it_geom=p_dsh->GetGeomIterator(); it_geom.IsValid(); ++it_geom )
        {
            const geo::bv::IBoundingVolume* pBV( it_geom->GetBV() );
            if( 0 == pBV
                || (pBV->GetDimension() == 2
                    && geo::bv::TestRay( static_cast<const geo::bv::BoundingVolume2*>(pBV), ray )) ) //\todo ClipRay could reduce interval instead
            {
                const geo::IObject *pGO( it_geom->GetGO() );
                if( pGO->GetDimension() == 2
                    && geo::mp::TestRay( static_cast<const geo::IObject2*>(pGO), ray, rh, 0 ) )
                {
                    num_hits++;
                    Write_RayHitReport2( it_geom->GetUID(), rh, rets );
                }
            }
            /*TEMPORAL
            else if( 0 != it_geom->GetBV() )
                DS_LOG_WARNING( "Query_RayCast2_DSH: Discarded hit by Geom BV!" );
            */
        }
        // Iterate over children DSH, filter with optional BV
        for( ContainerDSH::iterator it_dsh=p_dsh->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
        {
            const geo::bv::IBoundingVolume* pBV( it_dsh->GetBV() );
            if( 0 == pBV
                || (pBV->GetDimension() == 2
                    && geo::bv::TestRay( static_cast<const geo::bv::BoundingVolume2*>(pBV), ray )) ) //\todo ClipRay could reduce interval instead
            {
                const geo::IObject *pGO( it_dsh->GetGO() );
                if( pGO->GetDimension() == 2
                    && geo::mp::TestRay( static_cast<const geo::IObject2*>(pGO), ray, rh, 0 ) )
                {
                    num_hits++;
                    Write_RayHitReport2( it_dsh->GetUID(), rh, rets );
                }
                // Recursive call
                num_hits += Query_RayCast2_DSH( it_dsh, ray, rets );
            }
            /*TEMPORAL
            else if( 0 != it_dsh->GetBV() )
                DS_LOG_WARNING( "Query_RayCast2_DSH: Discarded hit by DSH BV!" );
            */
        }
    }
    return num_hits;
}

/*\todo Raycast requires a BP? if so, it should be done by the
  SimulationScheme, not by the DSH. Or, Maybe, it could be done
  outside, using DSH/SS but not as a method, it's just a QUERY,
  does NOT change their state... maybe each LeafDSH/Geom will
  require an object-centric QueryRayCast(), but hirearchical
  query with optional BP is pretty generic.

  When querying nested BP, possible EXTERNAL proxies, as they are
  already tested in a higher level or must be ignored, if raycast
  is local!.

  If DSH or SS has BP
  - QueryRayCast on BP => Gather proxies => Get LeafDSH & Geom from proxies and test explicitly
  else
  - Recursively QueryRayCast on each sub-DSH => Gather LeafDSH & Geom and test explicitly
*/
unsigned int Query_RayCast3_DSH( const IDynamicSystemHierarchy *p_dsh, const geo::np::Ray3 &ray, ReturnStream &rets )
{
    geo::np::RayHit3 rh;
    unsigned int num_hits(0);
    if( 0 != p_dsh->GetSS() && 0 != p_dsh->GetSS()->GetBP() )
    {
        // Raycast BP to gather hit proxies
        std::vector<geo::bp::Proxy*> vec_hit_proxies;
        if( p_dsh->GetSS()->GetBP()->TestRay3( ray, vec_hit_proxies, geo::bp::eTest_Internal | geo::bp::eTest_Boundary ) > 0 )
        {
            //DS_LOG_WARNING( "Query_RayCast3_DSH: Using DSH->SS->BP! %d hits", vec_hit_proxies.size() );
            for( unsigned int it_proxy=0; it_proxy < vec_hit_proxies.size(); it_proxy++ )
            {
                if( vec_hit_proxies[it_proxy]->IsNestedBP() )
                {
                    //\todo Recursive call?
                    DS_LOG_WARNING( "Query_RayCast3_DSH: Nested BP not yet implemented" );
                }
                else if( vec_hit_proxies[it_proxy]->m_UserFlags.Test(ISimulationScheme::eSSPUF_Entity) )
                {
                    // Retrieve IEntity from proxy (safe as eSSPUF_Entity)
                    const IEntity *pEntity( vec_hit_proxies[it_proxy]->GetObj<const IEntity*>() );
                    // Retrieve GO and raycast it
                    const geo::IObject *pGO(0);
                    if( pEntity->IsGeom() ) pGO = static_cast<const IGeom*>(pEntity)->GetGO();
                    else if( pEntity->IsDSH() ) pGO = static_cast<const IDynamicSystemHierarchy*>(pEntity)->GetGO();
                    if( 0 != pGO
                        && pGO->GetDimension() == 3
                        && geo::mp::TestRay( static_cast<const geo::IObject3*>(pGO), ray, rh, 0 ) )
                    {
                        num_hits++;
                        Write_RayHitReport3( pEntity->GetUID(), rh, rets );
                    }
                    /*TEMPORAL
                    if( 0 == pGO ) DS_LOG_WARNING( "Query_RayCast3_DSH: DSH without a GO found" );
                    */
                    // Recursive call if DSH
                    if( pEntity->IsDSH() )
                        num_hits += Query_RayCast3_DSH( static_cast<const IDynamicSystemHierarchy*>(pEntity), ray, rets );
                }
                else
                {
                    DS_LOG_WARNING( "Query_RayCast3_DSH: Non eSSPUF_Entity BP proxy ignored" );
                }
            }
        }
        /*TEMPORAL
        else
            DS_LOG_WARNING( "Query_RayCast3_DSH: Using DSH->SS->BP! 0 hits" );
        */
    }
    else
    {
        // Iterate over Geoms, filter with optional BV
        for( ContainerGeoms::iterator it_geom=p_dsh->GetGeomIterator(); it_geom.IsValid(); ++it_geom )
        {
            const geo::bv::IBoundingVolume* pBV( it_geom->GetBV() );
            if( 0 == pBV
                || (pBV->GetDimension() == 3
                    && geo::bv::TestRay( static_cast<const geo::bv::BoundingVolume3*>(pBV), ray )) ) //\todo ClipRay could reduce interval instead
            {
                const geo::IObject *pGO( it_geom->GetGO() );
                if( pGO->GetDimension() == 3
                    && geo::mp::TestRay( static_cast<const geo::IObject3*>(pGO), ray, rh, 0 ) )
                {
                    num_hits++;
                    Write_RayHitReport3( it_geom->GetUID(), rh, rets );
                }
            }
            /*TEMPORAL
            else if( 0 != it_geom->GetBV() )
                DS_LOG_WARNING( "Query_RayCast3_DSH: Discarded hit by Geom BV!" );
            */
        }
        // Iterate over children DSH, filter with optional BV
        for( ContainerDSH::iterator it_dsh=p_dsh->GetChildrenIterator(); it_dsh.IsValid(); ++it_dsh )
        {
            const geo::bv::IBoundingVolume* pBV( it_dsh->GetBV() );
            if( 0 == pBV
                || (pBV->GetDimension() == 3
                    && geo::bv::TestRay( static_cast<const geo::bv::BoundingVolume3*>(pBV), ray )) ) //\todo ClipRay could reduce interval instead
            {
                const geo::IObject *pGO( it_dsh->GetGO() );
                if( pGO->GetDimension() == 3
                    && geo::mp::TestRay( static_cast<const geo::IObject3*>(pGO), ray, rh, 0 ) )
                {
                    num_hits++;
                    Write_RayHitReport3( it_dsh->GetUID(), rh, rets );
                }
                // Recursive call
                num_hits += Query_RayCast3_DSH( it_dsh, ray, rets );
            }
            /*TEMPORAL
            else if( 0 != it_dsh->GetBV() )
                DS_LOG_WARNING( "Query_RayCast3_DSH: Discarded hit by DSH BV!" );
            */
        }
    }
    return num_hits;
}

bool MoP::Query_RayCast( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Query_RayCast()...");
    IEntity *pEntity = Cast2Entity( pit.Find("eid").Get<machine_uint_type>() );
    DS_ASSERT( pEntity->IsDSH() ); //\todo COULD BE GEOM if any entity can receive a RayCast query!!
    IDynamicSystemHierarchy *pDSH = static_cast<IDynamicSystemHierarchy*>( pEntity );
    Flags32 ray_flags = pit.Find("flags").Get<Flags32>();


    // if( pDSH->GetDimension() == 2 ) //\todo NO! because the global Universe is 3D, and the BP is dimension-agnostic... THIS IS SHIT!
    if( pit.Find("pos").GetType() == eType_Vec2f )
    {
        geo::np::Ray2 ray( pit.Find("pos").Get<Vec2f>(),
                           pit.Find("dir").Get<Vec2f>(),
                           pit.Find("interval").Get<Intervalf>(),
                           pit.Find("thickness").Get<float32>() );

        /* \todo Consider var-sized compact array, instead of nested complex per ray_hit...
           rets.BeginArray( "array_hit", eType_RayHit );
           rets.EndArray();
        */
        uint32 num_hits(0);
        rets.BeginComplex( cid, (uint32)eRet_RayCast );
        {
            rets.Write( "qid", pit.Find("qid").Get<machine_uint_type>() );
            rets.BeginComplex( "array_hit", eType_Unknown ); //TEMP: Unknown should be "Generic", as in "valid-but-not-predefined"
            {
                num_hits = Query_RayCast2_DSH( pDSH, ray, rets );
            }
            rets.EndComplex();
            rets.Write( "num_hits", num_hits );
        }
        rets.EndComplex();
    }
    else
    {
        // DS_LOG("RC 3D!!!");
        geo::np::Ray3 ray( pit.Find("pos").Get<Vec3f>(),
                           pit.Find("dir").Get<Vec3f>(),
                           pit.Find("interval").Get<Intervalf>(),
                           pit.Find("thickness").Get<float32>() );

        /* \todo Consider var-sized compact array, instead of nested complex per ray_hit...
           rets.BeginArray( "array_hit", eType_RayHit );
           rets.EndArray();
        */
        uint32 num_hits(0);
        rets.BeginComplex( cid, (uint32)eRet_RayCast );
        {
            rets.Write( "qid", pit.Find("qid").Get<machine_uint_type>() );
            rets.BeginComplex( "array_hit", eType_Unknown ); //TEMP: Unknown should be "Generic", as in "valid-but-not-predefined"
            {
                num_hits = Query_RayCast3_DSH( pDSH, ray, rets );
            }
            rets.EndComplex();
            rets.Write( "num_hits", num_hits );
        }
        rets.EndComplex();

        // DS_LOG("RC 3D %d hits", num_hits );
    }

    DS_ECMD(true);
    return true;
}

//---- Entity debug commands
//! Log/Viz/Profile
bool MoP::Dbg_Conf( CommandID cid, const ParamIt &pit, ReturnStream &rets )
{
    DS_BCMD("ds::MoP::Dbg_Conf()...");
    // Find entity (todo: check!!)
    IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());

    EDbgModule dbgm = (EDbgModule)pit.Find("dbg_module").Get<uint32>();
    Flags32 flags( pit.Find("flags").Get<Flags32>() );

    // Configure Log/Viz/Prof
    switch( dbgm )
    {
    case eDbgModule_Log: pEntity->SetLogFlags(flags); break;
    case eDbgModule_Viz: pEntity->SetVizFlags(flags); break;
    case eDbgModule_Prof: pEntity->SetProfFlags(flags); break;
    case eDbgModule_Stats: break;
    case eDbgModule_Params: break;
    default: break;
    }

    DS_ECMD(true);
    return true;
}

//! Stats
bool MoP::Dbg_QueryStats( CommandID cid, const ParamIt &pit )
{
    DS_BCMD("ds::MoP::Dbg_QueryStats()...");
    IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());
    DSG::GetStatIS()->BeginComplex( cid, (uint32)eRet_Stats );
    {
        DSG::GetStatIS()->Write( "uid", pEntity->GetUID() );
        pEntity->QueryStats( pit.Find("flags").Get<Flags32>(), *DSG::GetStatIS() );
    }
    DSG::GetStatIS()->EndComplex();
    DS_ECMD(true);
    return true;
}

//! Params
bool MoP::Dbg_QueryParams( CommandID cid, const ParamIt &pit )
{
    DS_BCMD("ds::MoP::Dbg_QueryParams()...");
    IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());
    DSG::GetParamIS()->BeginComplex( cid, (uint32)eRet_Params );
    {
        DSG::GetParamIS()->Write( "uid", pEntity->GetUID() );
        pEntity->QueryParams( pit.Find("flags").Get<Flags32>(), *DSG::GetParamIS() );
    }
    DSG::GetParamIS()->EndComplex();
    DS_ECMD(true);
    return true;
}

bool MoP::Dbg_SyncParams( CommandID cid, const ParamIt &pit )
{
    DS_BCMD("ds::MoP::Dbg_SyncParams()...");
    IEntity *pEntity = Cast2Entity(pit.Find("eid").Get<machine_uint_type>());
    DSG::GetParamIS()->BeginComplex( cid, (uint32)eRet_Params );
    {
        DSG::GetParamIS()->Write( "uid", pEntity->GetUID() );
        pEntity->SyncParams( pit.Find("params"), *DSG::GetParamIS() );
    }
    DSG::GetParamIS()->EndComplex();
    DS_ECMD(true);
    return true;
}

} } // namespace S2::ds
