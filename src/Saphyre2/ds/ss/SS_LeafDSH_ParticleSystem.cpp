#include "SS_LeafDSH_ParticleSystem.h"
#include <Saphyre2/ds/dsh/LeafDSH_ParticleSystem.h>
#include <Geo/bp/BasicBP.h>

namespace S2 {
namespace ds {

SS_LeafDSH_ParticleSystem::SS_LeafDSH_ParticleSystem()
: m_pDSH(0)
, m_pBP(0)
, m_PairsBP(100)
{
}

SS_LeafDSH_ParticleSystem::~SS_LeafDSH_ParticleSystem()
{
    //delete BP, etc...
    if( m_pBP ) delete m_pBP;

    if( m_vecProxy ) delete m_vecProxy;
    if( 2 == m_Dimension && m_vecBV2 ) delete m_vecBV2;
    else if( 3 == m_Dimension && m_vecBV3 ) delete m_vecBV3;
}

bool SS_LeafDSH_ParticleSystem::BindDSH( IDynamicSystemHierarchy *p_dsh )
{
    if( eDSH_Leaf_ParticleSystem2D != p_dsh->GetDSHType() && eDSH_Leaf_ParticleSystem3D != p_dsh->GetDSHType() )
        return false;

    m_Dimension = p_dsh->GetDimension();
    m_pDSH = p_dsh;
    if( 2 == m_Dimension )
    {
        m_pPS = &static_cast<LeafDSH_ParticleSystem2D*>(p_dsh)->GetModel();
        m_pPS2 = static_cast<ms::ParticleSystem2D *>(m_pPS);
    }
    else
    {
        m_pPS = &static_cast<LeafDSH_ParticleSystem3D*>(p_dsh)->GetModel();
        m_pPS3 = static_cast<ms::ParticleSystem3D *>(m_pPS);
    }

    //\todo Notify possibly existing Geoms & Constraints using NotifyXXX()

    // Init per-particle annotations
    m_vecProxy = new const geo::bp::Proxy*[ m_pPS->GetNumParticles() ];
    if( 2 == m_Dimension ) m_vecBV2 = new geo::bv::Sphere2[ m_pPS->GetNumParticles() ];
    else m_vecBV3 = new geo::bv::Sphere3[ m_pPS->GetNumParticles() ];

    // Create Broad Phase
    m_pBP = new geo::bp::BasicBP( m_pPS->GetNumParticles() );
    for( unsigned int i=0; i<m_pPS->GetNumParticles(); i++ )
        m_vecProxy[i] = m_pBP->Add( reinterpret_cast<void*>(i),
                                    0, //\todo part_id could be particle index, and the PS entitity could be the p_obj!!
                                    (2==m_Dimension)
                                    ? static_cast<const geo::bv::IBoundingVolume*>(&m_vecBV2[i])
                                    : static_cast<const geo::bv::IBoundingVolume*>(&m_vecBV3[i]),
                                    (i<3) ? geo::bp::Proxy::eDefault : geo::bp::Proxy::eBoundary,
                                    eSSPUF_None ); //\note Not added with eSSPF_Entity => NO RayCast!!
    RecomputePerParticleBV();
    m_pBP->Update();
    return true;
}

void SS_LeafDSH_ParticleSystem::Update( Real dt )
{
    /*TEMP: Local update (Integrate + SyncGO + RecomputeBV)
      \todo Instead of local update, here we should accept different
      integration methods using generic IForce (FEM, mass and spring,
      gravitation, ...)
    */
    static_cast<IBaseDSH*>(m_pDSH)->Update_Local(dt);

    // Collision Detection:
    //   Update BV (=>also BP proxies)
    RecomputePerParticleBV();
    //   Perform BP
    m_pBP->Update();
    m_pBP->TestAll( m_PairsBP,
                    geo::bp::eTest_NonPersistent
                    | geo::bp::eTest_Discrete
                    | geo::bp::eTest_Boundary
                    | geo::bp::eTest_Internal,
                    geo::bp::IBroadPhase::CBPF_True );

    //     Perform MP/NP => Create/Update Contacts
    // Constraints:
    //   Solve Contacts & Constraints
}

//---- Internal methods
void SS_LeafDSH_ParticleSystem::RecomputePerParticleBV()
{
    /*!\todo HARDCODED RADIUS, should come from LeafDSH bound
      shape/geom, NEVER from an hypothetical "GetRadius(i)" in
      ms::ParticleSystem because that'll break the Model/Shape
      decoupling
    */
    if( 2 == m_Dimension )
    {
        Real radius(1.0f);
        if( m_pDSH->GetGO()
            && geo::eShape_SphereSet2 == m_pDSH->GetGO()->GetShapeInterface()->GetType() )
            radius = static_cast<const geo::SphereSetShape2 *>(m_pDSH->GetGO()->GetShapeInterface())->GetRadius();
        for( unsigned int i=0; i<m_pPS2->GetNumParticles(); i++ )
            m_vecBV2[i].SetPosRadius( m_pPS2->GetPos(i), radius );
    }
    else
    {
        Real radius(1.0f);
        if( m_pDSH->GetGO()
            && geo::eShape_SphereSet3 == m_pDSH->GetGO()->GetShapeInterface()->GetType() )
            radius = static_cast<const geo::SphereSetShape3 *>(m_pDSH->GetGO()->GetShapeInterface())->GetRadius();
        for( unsigned int i=0; i<m_pPS3->GetNumParticles(); i++ )
            m_vecBV3[i].SetPosRadius( m_pPS3->GetPos(i), radius );
    }
}

void SS_LeafDSH_ParticleSystem::DoViz( util::VizStream &vs ) const
{
    // BP pairs
    if( 2 == m_Dimension )
        for( geo::bp::PairContainer::iterator it_pair=m_PairsBP.Begin(); it_pair.IsValid(); ++it_pair )
        {
            vs.BeginComplex(1,util::eType_VizSegment);
            {
                vs.Write("pos1", mal::CastDimension<3,float,2>( m_vecBV2[it_pair->GetObj1<machine_uint_type>()].GetPos() ) );
                vs.Write("pos2", mal::CastDimension<3,float,2>( m_vecBV2[it_pair->GetObj2<machine_uint_type>()].GetPos() ) );
                vs.BeginComplex("style",util::eType_VizStyle);
                {
                    vs.Write("color",Vec4f(1,0,0,1));
                    vs.Write("pen_size",1.0f);
                }
                vs.EndComplex();
            }
            vs.EndComplex();
        }
    else
        for( geo::bp::PairContainer::iterator it_pair=m_PairsBP.Begin(); it_pair.IsValid(); ++it_pair )
        {
            vs.BeginComplex(1,util::eType_VizSegment);
            {
                vs.Write("pos1", Vec3f( m_vecBV3[it_pair->GetObj1<machine_uint_type>()].GetPos() ) );
                vs.Write("pos2", Vec3f( m_vecBV3[it_pair->GetObj2<machine_uint_type>()].GetPos() ) );
                vs.BeginComplex("style",util::eType_VizStyle);
                {
                    vs.Write("color",Vec4f(1,0,0,1));
                    vs.Write("pen_size",1.0f);
                }
                vs.EndComplex();
            }
            vs.EndComplex();
        }
}

}} //namespace S2::ds
