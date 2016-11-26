#include "Scene.h"
#include "Params.h"

// S2 Includes
#include <Saphyre2/bs/BSG.h>
#include <Saphyre2/bs/Universe.h>

#include <string.h> //for strlen()
#include <sstream>

Scene::Scene( const Params& params )
: m_rParams(params)
, m_Time(0)
, m_pUniverse(0)
{}

Scene::~Scene()
{
    delete m_pUniverse; //\todo ugly... should be a S2::BSG method? or pimpl? whatever except raw delete!
}

void Scene::Init()
{
    // Create Universe
    m_pUniverse = S2::BSG::CreateUniverse();
    m_pUniverse->BeginDef();
    {
        ;//SetParams()
    }
    m_pUniverse->EndDef(true);

    // Init debug flags
    m_pUniverse->DbgSetFlags_Viz( 0xFFFFFFFF );
}

bool Scene::Update( double dt )
{
    m_Time += dt;
    if( dt != 0.0 )
    {
        // Apply forces to Particles, by now
        for( MapSyncEntity::iterator it_s2e=m_mapSyncEntity.begin();
             it_s2e != m_mapSyncEntity.end();
             ++it_s2e )
        {
            S2::ISyncEntity *pS2E = it_s2e->second;
            pS2E->Lock();
            switch(pS2E->GetType())
            {
            case S2::eEntity_Particle2D: static_cast<S2::Particle2D*>(pS2E)->ApplyForce( S2::Vec2(0,-1) ); break;
            case S2::eEntity_Particle3D: static_cast<S2::Particle3D*>(pS2E)->ApplyForce( S2::Vec3(0,-1,0) ); break;
            default: break;
            }
            pS2E->Unlock();
        }

        // Update whole universe
        m_pUniverse->Update(dt);
        /* [ALTERNATIVE] Update each Object separatedly... useful for
           individual object advance, but NOT for different
           time-speeds (use time_scale for that)
        */

        if( m_rParams.scene.m_bLog ) S2::BSG::UpdateLog();
        if( m_rParams.scene.m_bProf ) S2::BSG::UpdateProf();

        if( m_rParams.scene.m_bViz )
            S2::BSG::UpdateViz();
        else //Clear explicitly to avoid re-drawing last update Viz stream contents
            S2::BSG::ClearViz();

        //S2::BSG::UpdateStats();
        //S2::BSG::SyncParams(); //??

        return true;
    }
    else
        return false;
}

const char *Scene::Add( S2::ISyncEntity *p_sync_entity, const char *name )
{
    // Init debug flags
    p_sync_entity->DbgSetFlags_Viz( 0xFFFFFFFF );

    // Add entity to map
    MapSyncEntity::const_iterator inserted_s2e;
    if( 0 == name || 0 == strlen(name) )
    {
        // Generate unique name from ptr
        std::stringstream default_name;
        default_name << "s2e[0x" << std::hex << reinterpret_cast<machine_uint_type>(p_sync_entity) << "]";
        inserted_s2e = m_mapSyncEntity.insert( std::make_pair( default_name.str() , p_sync_entity ) ).first;
        p_sync_entity->SetName( default_name.str().c_str() );
    }
    else
    {
        inserted_s2e = m_mapSyncEntity.insert( std::make_pair( std::string(name), p_sync_entity ) ).first;
        p_sync_entity->SetName( name );
    }
    return inserted_s2e->first.c_str();
}

bool Scene::Clear()
{
    for( MapSyncEntity::iterator it_s2e=m_mapSyncEntity.begin();
         it_s2e != m_mapSyncEntity.end();
         ++it_s2e )
        it_s2e->second->Destroy();
    m_mapSyncEntity.clear();
    return true;
}

bool Scene::Kill( const char *name )
{
    MapSyncEntity::iterator it_s2e;
    it_s2e = m_mapSyncEntity.find(std::string(name));
    if ( it_s2e != m_mapSyncEntity.end() )
    {
        S2::ISyncEntity *pSE( it_s2e->second );
        m_mapSyncEntity.erase(it_s2e);
        return pSE->Destroy();
    }
    else
        return false;
}

//---- Simple object creation
/*
S2::Rigid * Scene::CreateRigid( const S2::Point3 &pos, const S2::Quat &rot,
                                const S2::Vec3 &vel, const S2::Vec3 &vel_rot,
                                S2::Real mass, const S2::Vec3 &diagonal_inertia_tensor )
{
    S2::Rigid *pRigid = m_UniverseOLD.CreateRigid( m_pSystem );
    pRigid->BeginDefinition();
      pRigid->SetMass( mass );
      pRigid->SetBodyInertiaTensor( diagonal_inertia_tensor );
      pRigid->SetPos( pos );
      pRigid->SetOrientation( rot );
      pRigid->SetVel( vel );
      pRigid->SetVelRot( vel_rot );
    pRigid->EndDefinition();

    m_vecRigids.push_back( pRigid );
    return pRigid;
}
*/

//---- Custom object creation
/*
S2::Particle *Scene::CreateParticleChain( unsigned int num_objects,
                                          float handle_mass,
                                          float mass_factor,
                                          const S2::Point3 &handle_pos,
                                          const S2::Vec3 &direction, float separation )
{
    if( num_objects > 0 )
    {
        // Add n particles
        float mass = handle_mass;
        unsigned int num_existing_particles = m_vecParticles.size();
        for( unsigned int i=num_existing_particles; i<num_existing_particles + num_objects; i++ )
        {
            m_vecParticles.push_back( m_UniverseOLD.CreateParticle( m_pSystem ) );
            m_vecParticles[i]->BeginDefinition();
              m_vecParticles[i]->SetMass( mass );
              m_vecParticles[i]->SetPos( handle_pos + float(i-num_existing_particles)*separation*direction );
              m_vecParticles[i]->SetVel( S2::Vec3::Zero() );
            m_vecParticles[i]->EndDefinition();
            mass *= mass_factor;
        }
        // Add n-1 distance constraints
        for( unsigned int i=0; i<num_objects-1; i++ )
        {
            S2::ConstraintFDP2P *pConstraint = m_UniverseOLD.CreateConstraintFDP2P( m_pSystem );
            pConstraint->BeginDefinition();
              pConstraint->SetStabilizationMode( DEFAULT_STAB_MODE );
              pConstraint->SetConnectorPoint1( m_vecParticles[num_existing_particles+i]->GetConnectorPoint() );
              pConstraint->SetConnectorPoint2( m_vecParticles[num_existing_particles+i+1]->GetConnectorPoint() );
              pConstraint->SetDist( (m_vecParticles[num_existing_particles+i]->GetPos()
                                     - m_vecParticles[num_existing_particles+i+1]->GetPos()).Norm() );
            pConstraint->EndDefinition();
        }
        // Return chain handle
        return m_vecParticles[num_existing_particles];
    }
    return NULL;
}
*/

/*
S2::Rigid *Scene::CreateRigidChain( unsigned int num_objects,
                                    float handle_mass, const S2::Vec3 &handle_diagonal_inertia_tensor,
                                    float mass_factor,
                                    const S2::Point3 &handle_pos, const S2::Quat &handle_orientation,
                                    const S2::Vec3 &direction, float separation )
{
    if( num_objects > 0 )
    {
        // Add n particles
        float mass = handle_mass;
        unsigned int num_existing_rigids = m_vecRigids.size();
        for( unsigned int i=num_existing_rigids; i<num_existing_rigids + num_objects; i++ )
        {
            m_vecRigids.push_back( m_UniverseOLD.CreateRigid( m_pSystem ) );
            m_vecRigids[i]->BeginDefinition();
              m_vecRigids[i]->SetMass( mass );
              m_vecRigids[i]->SetBodyInertiaTensor( handle_diagonal_inertia_tensor );
              m_vecRigids[i]->SetPos( handle_pos + float(i-num_existing_rigids)*separation*direction );
              m_vecRigids[i]->SetOrientation( handle_orientation );
              m_vecRigids[i]->SetVel( S2::Vec3::Zero() );
              m_vecRigids[i]->SetVelRot( S2::Vec3::Zero() );
            m_vecRigids[i]->EndDefinition();
            mass *= mass_factor;
        }
        // Add n-1 distance constraints
        for( unsigned int i=0; i<num_objects-1; i++ )
        {
            S2::ConstraintFDP2P *pConstraint = m_UniverseOLD.CreateConstraintFDP2P( m_pSystem );
            pConstraint->BeginDefinition();
              pConstraint->SetStabilizationMode( DEFAULT_STAB_MODE );
              pConstraint->SetConnectorPoint1( m_vecRigids[num_existing_rigids+i]->GetConnectorPointLocal( S2::Point3::Zero()) );
              pConstraint->SetConnectorPoint2( m_vecRigids[num_existing_rigids+i+1]->GetConnectorPointLocal( S2::Point3(0,3,0) ) );
              pConstraint->SetDist( (m_vecRigids[num_existing_rigids+i]->GetPos()
                                     - m_vecRigids[num_existing_rigids+i+1]->GetPos()).Norm() );
            pConstraint->EndDefinition();
        }
        // Return chain handle
        return m_vecRigids[num_existing_rigids];
    }
    return NULL;
}
*/

/*
S2::Particle *Scene::CreateParticleChainRigidUFO( unsigned int num_objects,
                                                  float handle_mass,
                                                  float mass_factor,
                                                  const S2::Point3 &handle_pos,
                                                  const S2::Vec3 &direction, float separation,
                                                  S2::Rigid *p_rigid )
{
    S2::Particle *pParticleHandle = CreateParticleChain( num_objects,
                                                         handle_mass, mass_factor,
                                                         handle_pos,
                                                         direction, separation );
    if( p_rigid )
    {
        S2::ConstraintFDP2P *pConstraint = m_UniverseOLD.CreateConstraintFDP2P( m_pSystem );
        pConstraint->BeginDefinition();
          pConstraint->SetStabilizationMode( DEFAULT_STAB_MODE );
          pConstraint->SetConnectorPoint1( m_vecParticles.back()->GetConnectorPoint() );
          pConstraint->SetConnectorPoint2( p_rigid->GetConnectorPointLocal( S2::Point3::Zero() ) );
          pConstraint->SetDist( separation );
        pConstraint->EndDefinition();
    }
    return pParticleHandle;
}
*/

#if defined(ENABLE_TREE)
void Scene::CreateParticleTree( unsigned int num_levels )
{
    m_vecTree = new S2::Particle*[ (1<<num_levels) - 1 ];

    // Root
    m_vecTree[0] = m_UniverseOLD.CreateParticle( m_pSystem );
    m_vecTree[0]->BeginDefinition();
    m_vecTree[0]->SetMass( PARTICLE_M1 );
    m_vecTree[0]->SetPos( PARTICLE0_POS0 );
    m_vecTree[0]->SetVel( PARTICLE0_VEL0 );
    m_vecTree[0]->EndDefinition();

    // Tree
    for( unsigned int i=2; i<num_levels+1; i++ )
    {
        unsigned int num_objects_level_i = 1<<(i-1);
        unsigned int num_objects_before_level_i = num_objects_level_i-1;

        // Fill level i
        for( unsigned int j=0; j<num_objects_level_i; j++ )
        {
            unsigned int k = num_objects_before_level_i + j;
            unsigned int parent = (k - 1) / 2;
            int child_displ = ( (num_objects_before_level_i + j - 1) % 2 ) ? 1 : -1;

            m_vecTree[k] = m_UniverseOLD.CreateParticle( m_pSystem );
            m_vecTree[k]->BeginDefinition();
            m_vecTree[k]->SetMass( PARTICLE_M1 );
            m_vecTree[k]->SetPos( m_vecTree[parent]->GetPos()
                                  + V3D( 0, -10, 0)
                                  + V3D( child_displ * 5, 0, 0 ) );
            m_vecTree[k]->SetVel( PARTICLE0_VEL0 );
            m_vecTree[k]->EndDefinition();

            // Constraint to parent
            S2::ConstraintFDP2P *pConstraint = m_UniverseOLD.CreateConstraintFDP2P( m_pSystem );
            pConstraint->BeginDefinition();
            pConstraint->SetSolveMode( S2::Constraint::SOLVE_COUPLED );
            pConstraint->SetStabilizationMode( DEFAULT_STAB_MODE );
            pConstraint->SetObjects( m_vecTree[parent], m_vecTree[k] );
            pConstraint->SetDist( (m_vecTree[parent]->GetPos()-m_vecTree[k]->GetPos()).Modul() );
            pConstraint->EndDefinition();
        }
    }
}
#endif
