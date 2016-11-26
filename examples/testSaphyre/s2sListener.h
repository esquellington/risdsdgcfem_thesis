#ifndef TEST_SAPHYRE_S2S_LISTENER_H
#define TEST_SAPHYRE_S2S_LISTENER_H

#include "Config.h"
#include "Scene.h"
#include "SceneRenderer.h"
#include "AppTask.h"
#include "AppRenderer.h"

#include <Safra/Safra.h>
#include <Safra/task/Console.h>
#include <Safra/task/ProfilerMonitor.h>

#include <Geo/shape/shape.h> //All shape types

#include <boost/bind.hpp>

#ifdef __ENABLE_TEST_ID32
#  include <util/Id32.h>
#endif

#define __ENABLE_SOLID_PARAMS_EXTRA

#include <string>
bool string_has_suffix( const std::string& str, const std::string& suffix )
{
    return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

struct AppModules
{
    AppModules( AppTask *p_task,
                Scene *p_scene,
                SceneRenderer *p_scene_renderer,
                sfr::IView *p_viz_view,
                sfr::ProfilerMonitor *p_profiler_monitor,
                sfr::Console *p_console )
    : m_pTask(p_task)
    , m_pScene(p_scene)
    , m_pSceneRenderer(p_scene_renderer)
    , m_pVizView(p_viz_view)
    , m_pProfilerMonitor(p_profiler_monitor)
    , m_pConsole(p_console)
    {}

    AppTask *m_pTask;
    Scene *m_pScene;

    //! \name Optional Visualization "Modules"
    //@{
    SceneRenderer *m_pSceneRenderer;
    sfr::IView *m_pVizView;
    //@}

    //! \name Optional Modules that MAY have an associated View
    //@{
    sfr::ProfilerMonitor *m_pProfilerMonitor;
    sfr::Console *m_pConsole;
    //VariableMonitor
    //Logger
    //@}
};

/*! Script Listener that dispatches script commands to suitable app
    modules and implements non-local script commands that involve >1
    module.

    \note All modules are optional.
*/
class s2sListener
{
public:
    s2sListener( AppModules &app_modules )
    : m_rAppModules(app_modules)
    , m_Interpreter(m_LanguageDef)
    {
        // Extend language
        ExtendLanguageDef_scn(m_LanguageDef,"scn_");
        ExtendLanguageDef_viz(m_LanguageDef,"viz_");
        ExtendLanguageDef_stats(m_LanguageDef,"stats_");
        ExtendLanguageDef_params(m_LanguageDef,"params_");
        ExtendLanguageDef_mig2015(m_LanguageDef,"mig2015_");
        m_Interpreter.GetUserContext().Add("p_s2s_listener",reinterpret_cast<void*>(this));

        // Extend language with module-specific stuff
        //if( m_pProfilerMonitor ) m_pProfilerMonitor->ExtendLanguageDef(m_LanguageDef,"prof_");

        // Setup console language and interpreter
        if( m_rAppModules.m_pConsole )
        {
            m_rAppModules.m_pConsole->ExtendLanguageDef(m_LanguageDef,""); //Console cmds have no prefix
            m_rAppModules.m_pConsole->SetInterpreter(&m_Interpreter);
        }
    }

    bool Execute( const char *script, unsigned int length = 0 )
        {
            return m_Interpreter.Execute(script,length);
        }

private:
    AppModules &m_rAppModules;
    skr::LanguageDef m_LanguageDef;
    skr::Interpreter m_Interpreter;

private:
    //! \name s2s Time Manipulation Methods
    //@{
    static bool Method_Scene_Run( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        double duration = params.Get<double>("t");
        AppTask *pTask = pListener->m_rAppModules.m_pTask;
        if( !pTask ) return false;
        if( duration < 0 )
            pTask->SetActive(true);
        else
        {
            pTask->SetActive(true);
            pTask->Update(duration);
            pTask->SetActive(false);
        }
        return true;
    }
    static bool Method_Scene_Step( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        double t = mal::Max<double>(0,params.Get<double>("t"));
        double dt = mal::Max<double>(0,params.Get<double>("dt"));
        if( dt <= 0 ) dt = t;
        //std::cout << "Execute t = " << t << " dt = " << dt << std::endl;
        AppTask *pTask = pListener->m_rAppModules.m_pTask;
        if( !pTask ) return false;
        pTask->SetActive(true);
        unsigned int num_steps(0);
        for( double time = 0; time < t; time +=dt, num_steps++ )
            pTask->Update(dt);
        pTask->SetActive(false);
        //std::cout << "Executed #num_steps " << num_steps << std::endl;
        return true;
    }
    static bool Method_Scene_Pause( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        AppTask *pTask = pListener->m_rAppModules.m_pTask;
        if( !pTask ) return false;
        pTask->SetActive(false);
        return true;
    }
    static bool Method_Scene_Reset( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;
        return false;
    }
    //@}

    //! \name s2s Scene Query methods
    //@{
    static bool Method_Scene_List( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        //\todo iterate over scene, write entity names
        std::cout << "Listing Entities..." << std::endl;
        for( Scene::EntityIterator it_entity=pScene->GetEntityIterator(); it_entity.IsValid(); ++it_entity )
            std::cout << "Entity " << (*it_entity)->GetName() << " Type " << (*it_entity)->GetType() << std::endl;
        //\todo MAYBE add to results ?? results.AddString("name",obj_name);
        return true;
    }
    //@{

    //! \name s2s Scene Manipulation Methods
    //@{
    static bool Method_Scene_Particle3D( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Particle3D using S2 API
        S2::Particle3D *pParticle3D = pScene->GetUniverse()->CreateParticle3D();
        pParticle3D->BeginDef();
        {
            S2::Particle3D::Params p3d_params;
            //p3d_params.m_Flags;
            p3d_params.m_Mass = params.Get<double>("m");
            p3d_params.m_Radius = params.Get<double>("r");
            //p3d_params.m_CoeffRestitution;
            pParticle3D->SetParams(p3d_params);
            //pParticle3D->AttachShape( geo::BoxShape3( Vec3f(3,2,1) ) );
            pParticle3D->AttachShape( geo::SphereShape3( 1.5f ) );

            pParticle3D->SetPos( S2::Point3(params.GetArrayPtr<double>("p")) );
            pParticle3D->SetVel( S2::Vec3(params.GetArrayPtr<double>("v")) );
        }
        pParticle3D->EndDef(true);

        const char *obj_name = pScene->Add( pParticle3D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pParticle3D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Particle2D( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Particle2D using S2 API
        S2::Particle2D *pParticle2D = pScene->GetUniverse()->CreateParticle2D();
        pParticle2D->BeginDef();
        {
            S2::Particle2D::Params p2d_params;
            //p2d_params.m_Flags;
            p2d_params.m_Mass = params.Get<double>("m");
            p2d_params.m_Radius = params.Get<double>("r");
            //p2d_params.m_CoeffRestitution;
            pParticle2D->SetParams(p2d_params);
            pParticle2D->AttachShape( geo::SphereShape2( 1 ) );
            //pParticle2D->AttachShape( geo::BoxShape2( Vec2f(2,1) ) );

            pParticle2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pParticle2D->SetVel( S2::Vec2(params.GetArrayPtr<double>("v")) );
        }
        pParticle2D->EndDef(true);

        const char *obj_name = pScene->Add( pParticle2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pParticle2D) );
        results.AddString("name",obj_name);
        return true;
    }


    //! \name s2s Scene Manipulation Methods
    //@{
    static bool Method_Scene_ParticleSystem2D( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a ParticleSys2D using S2 API
        S2::ParticleSys2D *pParticleSys2D = pScene->GetUniverse()->CreateParticleSys2D();
        pParticleSys2D->BeginDef();
        {
            S2::ParticleSys2D::Params ps2d_params;
            //p3d_params.m_Flags;
            ps2d_params.m_NumParticles = (uint32)params.Get<double>("np");
            ps2d_params.m_TotalMass = params.Get<double>("m");
            //ps2d_params.m_CoeffRestitution;
            pParticleSys2D->SetParams(ps2d_params);

            for( unsigned int i=0; i<ps2d_params.m_NumParticles; i++ )
                pParticleSys2D->DefineParticle( i, S2::Point2(i,0) );
        }
        pParticleSys2D->EndDef(true);

        pParticleSys2D->Lock();
           pParticleSys2D->SetPos(0,S2::Vec2(0,1));
           pParticleSys2D->SetVel(1,S2::Vec2(0,1));
           pParticleSys2D->ApplyForce(2,S2::Vec2(0,1));
           pParticleSys2D->ApplyImpulse(3,S2::Vec2(0,1));
        pParticleSys2D->Unlock(true);

        const char *obj_name = pScene->Add( pParticleSys2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pParticleSys2D) );
        results.AddString("name",obj_name);
        return true;
    }

    /*! Creates a sheet of water with dimensions: 1m x 0.5m
      mass = volume*density_vol = thickness*surface*density_vol = surface * density_surf
      density_surf = thickness * density_vol
      Assuming water density_vol, density_surf = 1000kg/m^3 * 0.01m = 10kg/m^2
    */
    static bool Method_Scene_Fluid2D( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Fluid2D using S2 API
        S2::Fluid2D *pFluid2D = pScene->GetUniverse()->CreateFluid2D();
        pFluid2D->BeginDef();
        {
            S2::Fluid2D::Params f2d_params;
            f2d_params.m_Flags = 1;//flags=eFlags_KeiserSPH;
            f2d_params.m_NumParticles = (int)params.Get<double>("np");
            f2d_params.m_Density = params.Get<double>("density");
            f2d_params.m_Thickness = params.Get<double>("thickness");
            const double *p_shape_aabb = params.GetArrayPtr<double>("shape_aabb");
            f2d_params.m_InitShape_AABB_PosMin = S2::Point2(p_shape_aabb[0],p_shape_aabb[1]);
            f2d_params.m_InitShape_AABB_PosMax = S2::Point2(p_shape_aabb[2],p_shape_aabb[3]);
            const double *p_bounds_aabb = params.GetArrayPtr<double>("bounds_aabb");
            f2d_params.m_Bounds_AABB_PosMin = S2::Point2(p_bounds_aabb[0],p_bounds_aabb[1]);
            f2d_params.m_Bounds_AABB_PosMax = S2::Point2(p_bounds_aabb[2],p_bounds_aabb[3]);
            pFluid2D->SetParams(f2d_params);
        }
        pFluid2D->EndDef(true);

        // add to scene
        const char *obj_name = pScene->Add( pFluid2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pFluid2D) );
        results.AddString("name",obj_name);
        return true;
    }

    //! Creates a Solid2D
    static bool Method_Scene_Solid2D( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        int num_subd( mal::Clamp<int>( (int)params.Get<double>("subd"), 0, 4 ) );

        // Create a Solid2D using S2 API
        S2::Solid2D *pSolid2D = pScene->GetUniverse()->CreateSolid2D();
        if( 0 == g_pSolid2D ) g_pSolid2D = pSolid2D; //TEMPORAL to apply pressure
        pSolid2D->BeginDef();
        {
            S2::Solid2D::Params s2d_params;
            s2d_params.m_Flags = 0;
            s2d_params.m_FixedDT = params.Get<double>("dt");
            s2d_params.m_Density = params.Get<double>("density");
            s2d_params.m_Thickness = params.Get<double>("thickness");
            s2d_params.m_YoungModulus = params.Get<double>("young_modulus");
            s2d_params.m_PoissonRatio = params.Get<double>("poisson_ratio");
            s2d_params.m_DampingRatio = params.Get<double>("damping_ratio");
            s2d_params.m_PlasticYield = params.Get<double>("plastic_yield");
            s2d_params.m_PlasticMax = params.Get<double>("plastic_max");
            s2d_params.m_PlasticCreepPerSecond = params.Get<double>("plastic_creep_per_second");
            pSolid2D->SetParams(s2d_params);

            // Create mesh
            const char *mesh_file = params.GetString("file");
            if( std::string(mesh_file) != "" ) // Load mesh file
            {
                geo::ShapeLibrary sl;
                sl.Reserve( 1<<20 );
                sl.Load( mesh_file, string_has_suffix( std::string(mesh_file), ".bin" ) );
                geo::ShapeID sid2 = sl.GetShapeIdByName("SIM");
                APP_LOG( "Creating SIM shape with sid %d", sid2 );
                geo::ShapeDef shape_def = sl.Lookup( sid2 );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                pSolid2D->CreateMeshGO( shape_id );
                //Load embedded shape EMB, if exists
                {
                    geo::ShapeID sid3 = sl.GetShapeIdByName("EMB");
                    if( sid3 != geo::cInvalidShapeId )
                    {
                        APP_LOG( "Creating EMB shape with sid %d", sid3 );
                        geo::ShapeDef shape_def = sl.Lookup( sid3 );
                        geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                        pSolid2D->CreateEmbeddedGO( shape_id )->Embed( pSolid2D->GetMeshGO(), geo::eEM_Barycentric );
                    }
                }
            }
            else if( params.GetArrayCount("v") > 0 ) // Create vertex-list mesh
            {
                const double *array_double = params.GetArrayPtr<double>("v");
                uint32 array_count = params.GetArrayCount("v");
                SFR_ASSERT( array_count%6 == 0 );
                uint32 polygon_count = array_count/6;
                geo::EditableMeshSolidShape2 emss2;
                emss2.Clear();
                emss2.BeginEdition();
                {
                    for( uint32 it_p = 0; it_p < polygon_count; it_p++ )
                    {
                        geo::feature_index_type vid0 = emss2.AddVertex( geo::Vec2( array_double[6*it_p+0], array_double[6*it_p+1]) );
                        geo::feature_index_type vid1 = emss2.AddVertex( geo::Vec2( array_double[6*it_p+2], array_double[6*it_p+3]) );
                        geo::feature_index_type vid2 = emss2.AddVertex( geo::Vec2( array_double[6*it_p+4], array_double[6*it_p+5]) );
                        emss2.AddPolygon3( vid0, vid1, vid2 );
                        //APP_LOG("Creating elem (%d,%d,%d)",vid0,vid1,vid2);
                    }
                }
                emss2.EndEdition();
                for( int it_subd=0; it_subd < num_subd; it_subd++ ) emss2.Subdivide();
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( emss2 );
                pSolid2D->CreateMeshGO( shape_id );
            }
            else // Create regular mesh
            {
                int nx = (int) params.Get<double>("nx");
                int ny = (int) params.Get<double>("ny");
                Vec2f half_sizes = S2::Vec2(params.GetArrayPtr<double>("h")); //half_sizes
                geo::EditableMeshSolidShape2 emss2;
                if( nx > 1 && ny > 1 )
                    geo::Make_MeshSolidShape2_Box( emss2, half_sizes, nx, ny );
                else //Single-triangle mesh test
                {
                    emss2.Clear();
                    emss2.BeginEdition();
                    {
                        /*
                        geo::feature_index_type vid0 = emss2.AddVertex( geo::Vec2(-half_sizes.x(),-half_sizes.y()) );
                        geo::feature_index_type vid1 = emss2.AddVertex( geo::Vec2( half_sizes.x(),-half_sizes.y()) );
                        geo::feature_index_type vid2 = emss2.AddVertex( geo::Vec2( half_sizes.x(), half_sizes.y()) );
                        emss2.AddPolygon3( vid0, vid1, vid2 );
                        */

                        //TEMP: 1-ring face-neighbours around polygon0
                        const float sin60deg = mal::Sin( mal::Deg2Rad( 60.0f ) );
                        geo::feature_index_type vid0 = emss2.AddVertex( geo::Vec2( 0, 0 ) );
                        geo::feature_index_type vid1 = emss2.AddVertex( geo::Vec2( 1, 0 ) );
                        geo::feature_index_type vid2 = emss2.AddVertex( geo::Vec2( 2, 0 ) );
                        geo::feature_index_type vid3 = emss2.AddVertex( geo::Vec2( 1.5, sin60deg ) );
                        geo::feature_index_type vid4 = emss2.AddVertex( geo::Vec2( 1, 2*sin60deg ) );
                        geo::feature_index_type vid5 = emss2.AddVertex( geo::Vec2( 0.5, sin60deg ) );
                        emss2.AddPolygon3( vid1, vid3, vid5 );
                        emss2.AddPolygon3( vid0, vid1, vid5 );
                        emss2.AddPolygon3( vid1, vid2, vid3 );
                        emss2.AddPolygon3( vid5, vid3, vid4 );
                    }
                    emss2.EndEdition();
                }
                for( int it_subd=0; it_subd < num_subd; it_subd++ ) emss2.Subdivide();
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( emss2 );
                pSolid2D->CreateMeshGO( shape_id );
            }
        }
        pSolid2D->EndDef(true);

        // Add optional static KPC
        util::ItemStream::ItemIt kpc_it = params.GetIS().Find("%vec_kpc");
        if( kpc_it.IsValid() )
        {
            if( kpc_it.IsArray()
                && kpc_it.GetArrayCount() % 2 == 0
                && kpc_it.GetType() == pla_type_id<double>::value )
            {
                // Gather KPC
                unsigned int num_kpc( kpc_it.GetArrayCount() / 2 );
                const double *array_kpc( kpc_it.GetArrayPtr<double>() );
                S2::Vec2 vec_kpc_points[num_kpc];
                for( unsigned int i=0; i<num_kpc; i++ )
                {
                    vec_kpc_points[i].x() = array_kpc[2*i];
                    vec_kpc_points[i].y() = array_kpc[2*i+1];
                }

                // Add KPC
                pSolid2D->Lock();
                {
                    for( unsigned int i=0; i<num_kpc; i++ )
                        pSolid2D->CreateKinematicPointConstraint( vec_kpc_points[i], vec_kpc_points[i], S2::Vec2::Zero() );
                }
                pSolid2D->Unlock(true);

                //Add subd KPC
                if( num_subd > 0 )
                {
                    // Get Mesh
                    const S2::Solid2D::geo_object_type *pMSSGO( pSolid2D->GetMeshGO() );
                    const geo::MeshSolidShape2* pMSS( pMSSGO->GetShape() );

                    //for each KPC edge, add as KPC all contained nodes
                    pSolid2D->Lock();
                    {
                        for( unsigned int i=0; i<num_kpc; i++ )
                        {
                            S2::Vec2 p0( vec_kpc_points[i] );
                            S2::Vec2 p1( vec_kpc_points[i+1] );
                            for( unsigned int it_node=0; it_node<pMSS->GetNumV(); it_node++ )
                            {
                                S2::Vec2 p( pMSS->V_Pos( it_node, pMSS->GetVecDefaultSDOF() ) );
                                S2::Vec2 dir( mal::Normalized(p1-p0) );
                                S2::Vec2 normal( mal::PerpendicularCW(dir) );
                                S2::Vec2 diff( p - p0 );
                                S2::Real lambda( mal::Dot(diff, dir) / mal::Norm(p1-p0) );
                                if( mal::Abs(mal::Dot(diff,normal)) < 0.001 //on edge line
                                    && lambda > 0.01 && lambda < 0.99 )  //strictly inside edge, NO a verted
                                    pSolid2D->CreateKinematicPointConstraint( p, p, S2::Vec2::Zero() );
                            }
                        }
                    }
                    pSolid2D->Unlock(true);
                }
            }
            else
                SFR_LOG_ERROR( "scn_s2d: parameter vec_kpc wrong type" );
        }

        pSolid2D->Lock();
          pSolid2D->SetPosCoM( S2::Point2(params.GetArrayPtr<double>("p")) );
        pSolid2D->Unlock(true);

#ifdef __ENABLE_SOLID_PARAMS_EXTRA
        // Set extra params
        {
            S2::Solid2D::Params_Extra s2d_params_extra;
            double default_double_array_zero_2[] = { 0,0 };
            s2d_params_extra.m_Gravity = S2::Vec2( params.SafeGetArrayPtr<double>("%gravity", default_double_array_zero_2 ) );
            s2d_params_extra.m_MM = params.SafeGetString("%mm", "C_LCMH");
            s2d_params_extra.m_RM = params.SafeGetString("%rm", "DAPD" );
            s2d_params_extra.m_DM = params.SafeGetString("%dm", "T" );
            s2d_params_extra.m_Integrator = params.SafeGetString("%integrator", "QIE v" ); //"FIE dx" );
            s2d_params_extra.m_NR_Iter = (uint32)params.SafeGet<double>("%nr_iter", 3.0 );
            s2d_params_extra.m_NR_RelEpsilon = params.SafeGet<double>("%nr_rel_epsilon", 0.001 );
            s2d_params_extra.m_LS_Solver = params.SafeGetString("%ls_solver", "CG" ); //"GMRES" );
            s2d_params_extra.m_LS_Iter = (uint32)params.SafeGet<double>("%ls_iter", 10.0 );
            s2d_params_extra.m_LS_Restarts = (uint32)params.SafeGet<double>("%ls_restarts", s2d_params_extra.m_LS_Iter );
            s2d_params_extra.m_LS_RelEpsilon = params.SafeGet<double>("%ls_rel_epsilon", 0.001 );
            pSolid2D->Lock();
            pSolid2D->SetParamsExtra(s2d_params_extra);
            pSolid2D->Unlock(true);
        }
#endif

        // add to scene
        const char *obj_name = pScene->Add( pSolid2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pSolid2D) );
        results.AddString("name",obj_name);
        return true;
    }

    //! Creates a Solid3D
    static bool Method_Scene_Solid3D( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        int num_subd( mal::Clamp<int>( (int)params.Get<double>("subd"), 0, 4 ) );

        // Create a Solid3D using S2 API
        S2::Solid3D *pSolid3D = pScene->GetUniverse()->CreateSolid3D();
        if( 0 == g_pSolid3D ) g_pSolid3D = pSolid3D; //TEMPORAL
        pSolid3D->BeginDef();
        {
            S2::Solid3D::Params s3d_params;
            s3d_params.m_Flags = 0;
            s3d_params.m_FixedDT = params.Get<double>("dt");
            s3d_params.m_Density = params.Get<double>("density");
            s3d_params.m_YoungModulus = params.Get<double>("young_modulus");
            s3d_params.m_PoissonRatio = params.Get<double>("poisson_ratio");
            s3d_params.m_DampingRatio = params.Get<double>("damping_ratio");
            s3d_params.m_PlasticYield = params.Get<double>("plastic_yield");
            s3d_params.m_PlasticMax = params.Get<double>("plastic_max");
            s3d_params.m_PlasticCreepPerSecond = params.Get<double>("plastic_creep_per_second");
            pSolid3D->SetParams(s3d_params);

            // Create mesh
            const char *mesh_file = params.GetString("file");
            if( std::string(mesh_file) != "" ) // Load mesh file
            {
                geo::ShapeLibrary sl;
                sl.Reserve( 1<<20 );
                sl.Load( mesh_file, string_has_suffix( std::string(mesh_file), ".bin" ) );
                geo::ShapeID sid1 = sl.GetShapeIdByName("SIM"); //\todo Assumes "bundle" format
                if( sid1 == geo::cInvalidShapeId ) sid1 = sl.GetShapeIdByName("#1"); //\todo OLD ids
                APP_LOG( "Creating shape 'SIM' or '#1' with sid %d", sid1 );
                geo::ShapeDef shape_def = sl.Lookup( sid1 );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                pSolid3D->CreateMeshGO( shape_id );
                //Load embedded shape EMB, if exists
                {
                    geo::ShapeID sid3 = sl.GetShapeIdByName("EMB");
                    if( sid3 != geo::cInvalidShapeId )
                    {
                        APP_LOG( "Creating EMB shape with sid %d", sid3 );
                        geo::ShapeDef shape_def = sl.Lookup( sid3 );
                        geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                        static int s_num_solid3d = 0;
                        if( true ) //s_num_solid3d == 0 ) //TEMP: This was for B vs MLS comparison in Thesis, disabled for MIG2015 MultiDSH tests
                        {
                            pSolid3D->CreateEmbeddedGO( shape_id )->Embed( pSolid3D->GetMeshGO(), geo::eEM_Barycentric );
                            s_num_solid3d++;
                        }
                        else
                            pSolid3D->CreateEmbeddedGO( shape_id )->Embed( pSolid3D->GetMeshGO(), geo::eEM_MLS );
                    }
                }
            }
            else if( false )
            {
                geo::EditableTetSolidShape3 etss3;
                geo::Make_TetSolidShape3_Box( etss3, geo::Vec3(2,2,2) );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( etss3 );
                pSolid3D->CreateMeshGO( shape_id );
            }
            else //TEMP: Single tetrahedron
            {
                geo::EditableTetSolidShape3 etss3;
                etss3.BeginEdition();
                {
                    etss3.AddVertex( S2::Vec3(0,0,0) );
                    etss3.AddVertex( S2::Vec3(1,0,0) );
                    etss3.AddVertex( S2::Vec3(0.45,1,0.3) );
                    etss3.AddVertex( S2::Vec3(0,0,1) );
                    etss3.AddTetrahedron( 0, 1, 2, 3 );
                }
                etss3.EndEdition();
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( etss3 );
                pSolid3D->CreateMeshGO( shape_id );
            }
        }
        pSolid3D->EndDef(true);

        pSolid3D->Lock();
          pSolid3D->SetPosCoM( S2::Point3(params.GetArrayPtr<double>("p")) );
          pSolid3D->SetRotCoM( mal::GRotation3x3_From(Vec3f(1,0,0),(float)params.Get<double>("ax"))
                               * mal::GRotation3x3_From(Vec3f(0,1,0),(float)params.Get<double>("ay"))
                               * mal::GRotation3x3_From(Vec3f(0,0,1),(float)params.Get<double>("az")) );
        pSolid3D->Unlock(true);

        // Add optional static KPC
        util::ItemStream::ItemIt kpc_it = params.GetIS().Find("%vec_kpc");
        if( kpc_it.IsValid() )
        {
            if( kpc_it.IsArray()
                && kpc_it.GetArrayCount() % 3 == 0
                && kpc_it.GetType() == pla_type_id<double>::value )
            {
                // Gather KPC
                unsigned int num_kpc( kpc_it.GetArrayCount() / 3 );
                const double *array_kpc( kpc_it.GetArrayPtr<double>() );
                S2::Vec3 vec_kpc_points[num_kpc];
                for( unsigned int i=0; i<num_kpc; i++ )
                {
                    vec_kpc_points[i].x() = array_kpc[3*i];
                    vec_kpc_points[i].y() = array_kpc[3*i+1];
                    vec_kpc_points[i].z() = array_kpc[3*i+2];
                }

                // Add KPC
                pSolid3D->Lock();
                {
                    for( unsigned int i=0; i<num_kpc; i++ )
                        pSolid3D->CreateKinematicPointConstraint( vec_kpc_points[i], vec_kpc_points[i], S2::Vec3::Zero() );
                }
                pSolid3D->Unlock(true);
            }
        }

#ifdef __ENABLE_SOLID_PARAMS_EXTRA
        // Set extra params
        {
            S2::Solid3D::Params_Extra s3d_params_extra;
            double default_double_array_zero_3[] = { 0,0,0 };
            s3d_params_extra.m_Gravity = S2::Vec3( params.SafeGetArrayPtr<double>("%gravity", default_double_array_zero_3 ) );
            s3d_params_extra.m_MM = params.SafeGetString("%mm", "C_LCMH");
            s3d_params_extra.m_RM = params.SafeGetString("%rm", "DAPD" );
            s3d_params_extra.m_DM = params.SafeGetString("%dm", "T" );
            s3d_params_extra.m_Integrator = params.SafeGetString("%integrator", "QIE v" ); //"FIE dx" );
            s3d_params_extra.m_NR_Iter = (uint32)params.SafeGet<double>("%nr_iter", 3.0 );
            s3d_params_extra.m_NR_RelEpsilon = params.SafeGet<double>("%nr_rel_epsilon", 0.001 );
            s3d_params_extra.m_LS_Solver = params.SafeGetString("%ls_solver", "MINRES" );//"CG" ); //"GMRES" );
            s3d_params_extra.m_LS_Iter = (uint32)params.SafeGet<double>("%ls_iter", 10.0 );
            s3d_params_extra.m_LS_Restarts = (uint32)params.SafeGet<double>("%ls_restarts", s3d_params_extra.m_LS_Iter );
            s3d_params_extra.m_LS_RelEpsilon = params.SafeGet<double>("%ls_rel_epsilon", 0.001 );
            //Extra MIG2015
            s3d_params_extra.m_CS_PosCST = params.SafeGetString("%cs_posm", "Alteration" );
            s3d_params_extra.m_CS_VelCST = params.SafeGetString("%cs_velm", "Reaction" );
            s3d_params_extra.m_CS_DepthMargin = params.SafeGet<double>("%cs_depth_margin", 0.01 );
            s3d_params_extra.m_CS_PenaltyKs = params.SafeGet<double>("%cs_penalty_ks", 1000 );
            s3d_params_extra.m_CS_PenaltyKd = params.SafeGet<double>("%cs_penalty_kd", 10 );
            s3d_params_extra.m_CS_RelaxCoeff = params.SafeGet<double>("%cs_relax_coeff", 0.1 );
            s3d_params_extra.m_CS_FrictionS = params.SafeGet<double>("%cs_friction_s", 0 );
            s3d_params_extra.m_CS_FrictionD = params.SafeGet<double>("%cs_friction_d", 0 );
            //Extra Thesis
            //\todo s3d_params_extra.m_TestType = params.SafeGetString("%test_type", "Random" );
            // s3d_params_extra.m_TestFraction = params.SafeGet<double>("%test_fraction", 0 );
            pSolid3D->Lock();
            pSolid3D->SetParamsExtra(s3d_params_extra);
            pSolid3D->Unlock(true);
        }
#endif

        // add to scene
        const char *obj_name = pScene->Add( pSolid3D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pSolid3D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Sphere( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pKine2D->SetAngle( params.Get<double>("a") );
            //pKine2D->AttachShape( 122323 );
            pKine2D->AttachShape( geo::SphereShape2( params.Get<double>("r") ) );
        }
        pKine2D->EndDef(true);

        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine3D_Sphere( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine3D using S2 API
        S2::Kine3D *pKine3D = pScene->GetUniverse()->CreateKine3D();
        pKine3D->BeginDef();
        {
            S2::Kine3D::Params k3d_params;
            //k3d_params.m_Flags;
            pKine3D->SetParams(k3d_params);
            pKine3D->SetPos( S2::Point3(params.GetArrayPtr<double>("p")) );
            //pKine3D->SetAngle( params.Get<double>("a") );
            //pKine3D->AttachShape( 122323 );
            pKine3D->AttachShape( geo::SphereShape3( params.Get<double>("r") ) );
        }
        pKine3D->EndDef(true);

        const char *obj_name = pScene->Add( pKine3D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine3D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Capsule( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pKine2D->SetAngle( params.Get<double>("a") );
            //pKine2D->AttachShape( 122323 );
            pKine2D->AttachShape( geo::CapsuleShape2( params.Get<double>("r"), params.Get<double>("hh") ) );
        }
        pKine2D->EndDef(true);

        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Box( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pKine2D->SetAngle( params.Get<double>("a") );
            pKine2D->AttachShape( geo::BoxShape2( S2::Vec2(params.GetArrayPtr<double>("h")) ) );
        }
        pKine2D->EndDef(true);

        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Plane( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(0,0) );
            pKine2D->SetAngle( 0 );
            pKine2D->AttachShape( geo::PlaneShape2( S2::Vec2(params.GetArrayPtr<double>("n")).Normalize(),
                                                    params.Get<double>("d"),
                                                    0 != params.Get<double>("hs") ) );
        }
        pKine2D->EndDef(true);

        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine3D_Plane( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine3D using S2 API
        S2::Kine3D *pKine3D = pScene->GetUniverse()->CreateKine3D();
        pKine3D->BeginDef();
        {
            S2::Kine3D::Params k3d_params;
            //k3d_params.m_Flags;
            pKine3D->SetParams(k3d_params);
            pKine3D->SetPos( S2::Point3(0,0,0) );
            // pKine3D->SetAngle( 0 );
            pKine3D->AttachShape( geo::PlaneShape3( S2::Vec3(params.GetArrayPtr<double>("n")).Normalize(),
                                                    params.Get<double>("d"),
                                                    0 != params.Get<double>("hs") ) );
        }
        pKine3D->EndDef(true);

        const char *obj_name = pScene->Add( pKine3D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine3D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Polygonal( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pKine2D->SetAngle( params.Get<double>("a") );

            // Load mesh or use explicit vertex array
            const char *poly_file = params.GetString("file");
            if( std::string(poly_file) != "" ) // Load poly file
            {
                geo::ShapeLibrary sl;
                sl.Reserve( 1<<20 );
                sl.Load( poly_file, string_has_suffix( std::string(poly_file), ".bin" ) );
                geo::ShapeID sid3 = sl.GetShapeIdByName("#3"); //\todo Assumes shape library that contains (??,??,polygonal2...)
                APP_LOG( "Creating shape '#3' with sid %d", sid3 );
                geo::ShapeDef shape_def = sl.Lookup( sid3 );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                pKine2D->AttachShape( shape_id );
            }
            else if( params.GetArrayCount("v") > 0 ) // Create vertex-list mesh
            {
                // Edit the polygonal...
                const double *array_double = params.GetArrayPtr<double>("v");
                uint32 array_count = params.GetArrayCount("v");
                SFR_ASSERT( array_count%2 == 0 );
                uint32 vtx_count = array_count/2;
                geo::EditablePolygonalShape2 eps2;
                eps2.Alloc( vtx_count );
                eps2.BeginEdition();
                {
                    for( unsigned int i=0; i<vtx_count; i++ )
                        eps2.AddPoint( geo::Vec2(array_double[2*i], array_double[2*i+1]) );
                    eps2.SetRadius( params.Get<double>("r") );
                    eps2.SetClosed( (bool)params.Get<double>("c") );
                }
                eps2.EndEdition();
                eps2.Subdivide( (int)params.Get<double>("s") );
                eps2.Distort( params.Get<double>("d") );
                pKine2D->AttachShape( eps2 );
            }
        }
        pKine2D->EndDef(true);

        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Path( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pKine2D->SetAngle( params.Get<double>("a") );
            // Edit the path...
            const double *array_double = params.GetArrayPtr<double>("v");
            uint32 array_count = params.GetArrayCount("v");
            SFR_ASSERT( array_count%2 == 0 );
            uint32 vtx_count = array_count/2;

            // Create mesh
            const char *path_file = params.GetString("file");
            if( std::string(path_file) != "" )
            {
                geo::ShapeLibrary sl;
                sl.Reserve( 1<<20 );
                sl.Load( path_file, string_has_suffix( std::string(path_file), ".bin" ) );
                geo::ShapeID sid1 = sl.GetShapeIdByName("#1");
                APP_LOG( "Creating shape '#1' with sid %d", sid1 );
                geo::ShapeDef shape_def = sl.Lookup( sid1 );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                pKine2D->AttachShape( shape_id );
            }
            else
            {
                //TEMPORAL
                geo::EditablePathShape2 eps2;
                geo::Make_PathShape2_Mushroom( eps2, geo::Vec2(2,4) );

                /* TEMP: These do not seem to work anymore... maybe too large?
                //geo::Make_PathShape2_SVG( eps2, "M -2.0203053 630.11842 C 65.536226 596.53831 58.64823 538.92729 112.12693 523.04225 c 86.22553 -25.61198 -22.223353 93.94419 42.42641 100.0051 64.64976 6.06092 83.84266 47.47717 52.52793 56.56854 -31.31473 9.09138 -153.543185 35.35534 -209.1015753 -49.49747 z" );

                //geo::Make_PathShape2_SVG( eps2, "M 80.812203 324.04219 34.450739 321.1832 22.843316 276.20735 62.030998 251.26974 97.85774 280.83329 z" );
                geo::Make_PathShape2_SVG( eps2, "m 412.14225 176.55992 -33.07005 -44.27136 -55.21761 2.14419 31.88537 -45.132085 -19.10243 -51.852476 52.77628 16.378206 43.41166 -34.190789 0.73217 55.254374 45.93232 30.72141 -52.32378 17.77087 z" );
                */

                //
                /* TEMP: reenable vertex definition
                   eps2.BeginEdition();
                   {
                   eps2.SetFirstPoint( geo::Vec2(array_double[0], array_double[1]) );
                   for( unsigned int i=1; i<vtx_count; i++ )
                   eps2.AddLineTo( geo::Vec2(array_double[2*i], array_double[2*i+1]) );
                   if( (bool)params.Get<double>("c") ) eps2.Close();
                   }
                   eps2.EndEdition();
                */
                /*
                  eps2.Subdivide( (int)params.Get<double>("s") );
                  eps2.Distort( params.Get<double>("d") );
                */
                pKine2D->AttachShape( eps2 );
            }
        }
        pKine2D->EndDef(true);

        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_MeshSolid( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(params.GetArrayPtr<double>("p")) );
            pKine2D->SetAngle( params.Get<double>("a") );
            // Load mesh file or create parametric rectangular mesh
            const char *mesh_file = params.GetString("file");
            if( std::string(mesh_file) != "" )
            {
                geo::ShapeLibrary sl;
                sl.Reserve( 1<<20 );
                sl.Load( mesh_file, string_has_suffix( std::string(mesh_file), ".bin" ) );
                geo::ShapeID sid2 = sl.GetShapeIdByName("SIM");
                APP_LOG( "Creating SIM shape with sid %d", sid2 );
                geo::ShapeDef shape_def = sl.Lookup( sid2 );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                pKine2D->AttachShape( shape_id );

                //Load embedded shape EMB, if exists
                // {
                //     geo::ShapeID sid3 = sl.GetShapeIdByName("EMB");
                //     if( sid3 != geo::cInvalidShapeId )
                //     {
                //         APP_LOG( "Creating EMB shape with sid %d", sid3 );
                //         geo::ShapeDef shape_def = sl.Lookup( sid3 );
                //         geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                //         pKine2D->CreateEmbeddedGO( shape_id )->Embed( pKine2D->GetMeshGO(), geo::eEM_Barycentric );
                //     }
                // }
            }
            else
            {
                // Edit the mesh
                geo::EditableMeshSolidShape2 emss2;
                /*
                  geo::Make_MeshSolidShape2_BoxWithHole( emss2,
                  S2::Vec2(params.GetArrayPtr<double>("h")), //half_sizes
                  0.5f*S2::Vec2(params.GetArrayPtr<double>("h")) ); //half_sizes
                */
                //
                geo::Make_MeshSolidShape2_Box( emss2,
                                               S2::Vec2(params.GetArrayPtr<double>("h")), //half_sizes
                                               (int) params.Get<double>("nx"),
                                               (int) params.Get<double>("ny") );
                //

                /* Test edition-modification
                   emss2.BeginEdition();
                   {
                   geo::feature_id vid = emss2.AddVertex( Vec2(0,-1) );
                   emss2.AddPolygon3( 1, 0, vid );
                   }
                   emss2.EndEdition();
                */
                pKine2D->AttachShape( emss2 );
            }
        }
        pKine2D->EndDef(true);

#if __cplusplus > 199711L //C++11 FTW!
        //TEMP: Add BVH because they're not yet serialized
        geo::BVH_MeshSolidShape2* pBVH = new geo::BVH_MeshSolidShape2;
        geo::MeshSolidShape2* pMSS = static_cast<geo::MeshSolidShape2*>( pKine2D->GetShape() );
        uint32 num_entries_bvh( pMSS->GetDCR() ? pMSS->GetDCR()->m_NumElements : pMSS->GetNumP() );
        pBVH->Rebuild_TopDown( num_entries_bvh,
                               boost::bind<void>( &geo::GEBV_MeshSolidShape2_E<geo::BVH_MeshSolidShape2::entry_index_type,geo::BVH_MeshSolidShape2::bv_type>,
                                                  pMSS, geo::Transform2::Identity(), pMSS->GetVecDefaultSDOF(),
                                                  _1, _2 ) );
        pMSS->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH );
#endif
        const char *obj_name = pScene->Add( pKine2D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine2D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine3D_TetSolid( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        // Create a Kine2D using S2 API
        S2::Kine3D *pKine3D = pScene->GetUniverse()->CreateKine3D();
        pKine3D->BeginDef();
        {
            S2::Kine3D::Params k3d_params;
            //k2d_params.m_Flags;
            pKine3D->SetParams(k3d_params);
            pKine3D->SetPos( S2::Point3(params.GetArrayPtr<double>("p")) );
            pKine3D->SetRot( mal::GRotation3x3_From(Vec3f(1,0,0),(float)params.Get<double>("ax"))
                             * mal::GRotation3x3_From(Vec3f(0,1,0),(float)params.Get<double>("ay"))
                             * mal::GRotation3x3_From(Vec3f(0,0,1),(float)params.Get<double>("az")) );

            // Load mesh file or create parametric rectangular mesh
            const char *mesh_file = params.GetString("file");
            if( std::string(mesh_file) != "" )
            {
                geo::ShapeLibrary sl;
                sl.Reserve( 1<<20 );
                sl.Load( mesh_file, string_has_suffix( std::string(mesh_file), ".bin" ) );
                geo::ShapeID sid1 = sl.GetShapeIdByName("SIM"); //\todo Assumes "bundle" format
                if( sid1 == geo::cInvalidShapeId ) sid1 = sl.GetShapeIdByName("#1"); //\todo OLD ids
                APP_LOG( "Creating shape 'SIM' or '#1' with sid %d", sid1 );
                geo::ShapeDef shape_def = sl.Lookup( sid1 );
                geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                /*\todo APP: Load embedded shape EMB, if exists
                pSolid3D->CreateMeshGO( shape_id );
                {
                    geo::ShapeID sid3 = sl.GetShapeIdByName("EMB");
                    if( sid3 != geo::cInvalidShapeId )
                    {
                        APP_LOG( "Creating EMB shape with sid %d", sid3 );
                        geo::ShapeDef shape_def = sl.Lookup( sid3 );
                        geo::ShapeID shape_id = S2::BSG::GetShapeLibrary().Register( shape_def );
                        pSolid3D->CreateEmbeddedGO( shape_id )->Embed( pSolid3D->GetMeshGO(), geo::eEM_Barycentric );
                    }
                }
                */
                pKine3D->AttachShape( shape_id );
            }
            else
            {
                // Edit the mesh
                geo::EditableTetSolidShape3 etss3;
                geo::Make_TetSolidShape3_Box( etss3,
                                              S2::Vec3(params.GetArrayPtr<double>("h")), //half_sizes
                                              (int) params.Get<double>("nx"),
                                              (int) params.Get<double>("ny"),
                                              (int) params.Get<double>("nz") );
                pKine3D->AttachShape( etss3 );
            }
        }
        pKine3D->EndDef(true);

        //TEMP: Check DCR
        {
            const geo::DCR_TetSolidShape3* pDCR( static_cast<const geo::TetSolidShape3*>(pKine3D->GetShape())->GetDCR() );
            if( pDCR )
            {
                APP_LOG( "Checking Obj[%s].DCR...", params.GetString("name") );
                if( geo::Check_DCR_TetSolidShape3(*pDCR) )
                {
                    APP_LOG( "...OK" );
                }
                else
                {
                    APP_LOG( "...ERRORS!" );
                }
            }
        }

#if __cplusplus > 199711L //C++11 FTW!
        //TEMP: Add BVH because they're not yet serialized
        geo::BVH_TetSolidShape3* pBVH = new geo::BVH_TetSolidShape3;
        geo::TetSolidShape3* pTSS = static_cast<geo::TetSolidShape3*>( pKine3D->GetShape() );
        uint32 num_entries_bvh( pTSS->GetDCR() ? pTSS->GetDCR()->m_NumElements : pTSS->GetNumT() );
        pBVH->Rebuild_TopDown( num_entries_bvh,
                               boost::bind<void>( &geo::GEBV_TetSolidShape3_E<geo::BVH_TetSolidShape3::entry_index_type,geo::BVH_TetSolidShape3::bv_type>,
                                                  pTSS, geo::Transform3::Identity(), pTSS->GetVecDefaultSDOF(),
                                                  _1, _2 ) );
        pTSS->SetBakedBVH_StrictlyNonshared_UglyHack( pBVH );
#endif
        const char *obj_name = pScene->Add( pKine3D, params.GetString("name") );
        results.Add("ptr", reinterpret_cast<void*>(pKine3D) );
        results.AddString("name",obj_name);
        return true;
    }

    static bool Method_Scene_Kine2D_Cradle( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;

        //---- Create 3 Kine2D planes using S2 API
        double cradle_width( params.Get<double>("w") );
        double cradle_depth( params.Get<double>("d") );
        double cradle_angle( params.Get<double>("a") );
        // Base
        S2::Kine2D *pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(0,0) );
            pKine2D->SetAngle( 0 );
            pKine2D->AttachShape( geo::PlaneShape2( S2::Vec2(0,1), cradle_depth, true ) );
        }
        pKine2D->EndDef(true);
        pScene->Add( pKine2D, "" );

        // Compute wall plane coeff_d
        double tan_cradle_angle( mal::Tan(cradle_angle) );
        double cradle_wall_d;
        if( mal::Abs(tan_cradle_angle) > 0.0001 )
            cradle_wall_d = (0.5*cradle_width + cradle_depth / mal::Tan(cradle_angle)) * mal::Sin(cradle_angle);
        else
            cradle_wall_d = 0.5 * cradle_width * mal::Sin(cradle_angle);

        // Left wall
        pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(0,0) );
            pKine2D->SetAngle( 0 );
            pKine2D->AttachShape( geo::PlaneShape2( S2::Vec2( mal::Sin(cradle_angle),
                                                              mal::Cos(cradle_angle) ),
                                                    cradle_wall_d, true ) );
        }
        pKine2D->EndDef(true);
        pScene->Add( pKine2D, "" );

        // Right wall
        pKine2D = pScene->GetUniverse()->CreateKine2D();
        pKine2D->BeginDef();
        {
            S2::Kine2D::Params k2d_params;
            //k2d_params.m_Flags;
            pKine2D->SetParams(k2d_params);
            pKine2D->SetPos( S2::Point2(0,0) );
            pKine2D->SetAngle( 0 );
            pKine2D->AttachShape( geo::PlaneShape2( S2::Vec2( -mal::Sin(cradle_angle),
                                                               mal::Cos(cradle_angle) ),
                                                    cradle_wall_d, true ) );
        }
        pKine2D->EndDef(true);
        pScene->Add( pKine2D, "" );

        return true;
    }

    static bool Method_Scene_Kill( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;
        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
            return pScene->Clear();
        else
            return pScene->Kill( obj_name );
    }

    static bool Method_Scene_SaveState( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;
        const char *obj_name = params.GetString("name");
        const char *file_name = params.GetString("file");
        util::ItemStream its( 1<<12, 1<<10 );
        bool bResult(false);
        if( 0 == strcmp(obj_name,"*") )
        {
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
            {
                its.BeginComplex( (*it)->GetName(), eType_Unknown );
                (*it)->Dbg_SaveState( its );
                its.EndComplex();
            }
            bResult = true;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE )
            {
                its.BeginComplex( pSE->GetName(), eType_Unknown );
                pSE->Dbg_SaveState( its );
                its.EndComplex();
                bResult = true;
            }
            else
                bResult = false;
        }
        // Save to file
        if( bResult )
        {
            its.SaveTxt( file_name );
        }
        return bResult;
    }

    //---- Viz methods
    static bool Method_Viz_On( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;
        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
        {
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
                (*it)->DbgSetFlags_Viz( (*it)->DbgGetFlags_Viz().Enable(0xFFFFFFFF) );
            return true;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE ) pSE->DbgSetFlags_Viz( pSE->DbgGetFlags_Viz().Enable(0xFFFFFFFF) );
            else return false;
        }
        return true;
    }
    static bool Method_Viz_Off( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;
        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
        {
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
                (*it)->DbgSetFlags_Viz( (*it)->DbgGetFlags_Viz().Disable(0xFFFFFFFF) );
            return true;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE ) pSE->DbgSetFlags_Viz( pSE->DbgGetFlags_Viz().Disable(0xFFFFFFFF) );
            else return false;
        }
        return true;
    }
    static bool Method_Viz_Toggle( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        if( !pScene ) return false;
        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
        {
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
                (*it)->DbgSetFlags_Viz( (*it)->DbgGetFlags_Viz().Toggle(0xFFFFFFFF) );
            return true;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE ) pSE->DbgSetFlags_Viz( pSE->DbgGetFlags_Viz().Toggle(0xFFFFFFFF) );
            else return false;
        }
        return true;
    }

    static void OnSyncStats( S2::ISyncEntity *p_entity, util::ItemStream::ItemItRW sit )
        {
            util::ItemStream::ItemItRW it = p_entity->Dbg_SyncStats();
            SFR_ASSERT( it == sit );
        }

    static bool Method_Stats_Show( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        sfr::IView *pVizView = pListener->m_rAppModules.m_pVizView;
        if( 0 == pScene || 0 == pVizView || 0 == pVizView->GetDesktop() )
        {
            SFR_LOG_WARNING( "s2sListener::Method_Stats_Show() Requires Scene, VizView and Desktop... ignoring" );
            return true;
        }

        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
        {
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
            {
                util::ItemStream::ItemItRW sit = (*it)->Dbg_QueryStats();
                if( sit.IsValid() )
                {
                    sfr::gui::IWidget *pPTW
                        = pVizView->GetDesktop()->CreateWidget_PropertyTree( (*it)->GetName(),
                                                                             sfr::gui::eDWF_Default,
                                                                             sfr::Vec2(50,50),
                                                                             sit,
                                                                             boost::bind( &s2sListener::OnSyncStats,
                                                                                          static_cast<S2::ISyncEntity*>(*it),
                                                                                          _1 ) );
                    pPTW->GetBasicAPI()->SetColor( 0.8,0.0,0.0,0.25 );
                    pPTW->GetBasicAPI()->Minimize();
                }
                else
                    SFR_LOG_WARNING("s2sListener::Method_Stats_Show() entity '%s' returned no valid stats, ignoring", (*it)->GetName() );
            }
            return true;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE )
            {
                util::ItemStream::ItemItRW sit = pSE->Dbg_QueryStats();
                if( sit.IsValid() )
                {
                    sfr::gui::IWidget *pPTW
                        = pVizView->GetDesktop()->CreateWidget_PropertyTree( pSE->GetName(),
                                                                             sfr::gui::eDWF_Default,
                                                                             sfr::Vec2(50,50),
                                                                             sit,
                                                                             boost::bind( &s2sListener::OnSyncStats,
                                                                                          pSE,
                                                                                          _1 ) );
                    pPTW->GetBasicAPI()->SetColor( 0.8,0.0,0.0,0.25 );
                    pPTW->GetBasicAPI()->Minimize();
                }
                else
                    SFR_LOG_WARNING("s2sListener::Method_Stats_Show() entity '%s' returned no valid stats, ignoring", obj_name );
            }
            else return false;
        }
        return true;
    }

    static void OnSyncParams( S2::ISyncEntity *p_entity, util::ItemStream::ItemItRW pit )
        {
            util::ItemStream::ItemItRW it = p_entity->Dbg_SyncParams();
            SFR_ASSERT( it == pit );
        }

    static bool Method_Params_Show( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        sfr::IView *pVizView = pListener->m_rAppModules.m_pVizView;

        if( 0 == pScene || 0 == pVizView || 0 == pVizView->GetDesktop() )
        {
            SFR_LOG_WARNING( "s2sListener::Method_Params_Show() Requires Scene, VizView and Desktop... ignoring" );
            return true;
        }

        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
        {
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
            {
                util::ItemStream::ItemItRW sit = (*it)->Dbg_QueryParams();
                if( sit.IsValid() )
                {
                    sfr::gui::IWidget *pPTW
                        = pVizView->GetDesktop()->CreateWidget_PropertyTree( (*it)->GetName(),
                                                                             sfr::gui::eDWF_Default,
                                                                             sfr::Vec2(50,50),
                                                                             sit,
                                                                             boost::bind( &s2sListener::OnSyncParams,
                                                                                          static_cast<S2::ISyncEntity*>(*it),
                                                                                          _1 ) );
                    pPTW->GetBasicAPI()->SetColor( 0.0,0.8,0.0,0.25 );
                    pPTW->GetBasicAPI()->Minimize();
                }
                else
                    SFR_LOG_WARNING("s2sListener::Method_Params_Show() entity '%s' returned no valid params, ignoring", (*it)->GetName() );
            }
            return true;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE )
            {
                util::ItemStream::ItemItRW sit = pSE->Dbg_QueryParams();
                if( sit.IsValid() )
                {
                    sfr::gui::IWidget *pPTW
                        = pVizView->GetDesktop()->CreateWidget_PropertyTree( pSE->GetName(),
                                                                             sfr::gui::eDWF_Default,
                                                                             sfr::Vec2(50,50),
                                                                             sit,
                                                                             boost::bind( &s2sListener::OnSyncParams,
                                                                                          pSE,
                                                                                          _1 ) );
                    pPTW->GetBasicAPI()->SetColor( 0.0,0.8,0.0,0.25 );
                    pPTW->GetBasicAPI()->Minimize();
                }
                else
                    SFR_LOG_WARNING("s2sListener::Method_Params_Show() entity '%s' returned no valid params, ignoring", obj_name );
            }
            else return false;
        }
        return true;
    }

    static void OnGetStats( S2::ISyncEntity *p_entity, util::ItemStream::ItemIt sit )
        {
            p_entity->Dbg_SyncStats(); //\todo This MAY be redundant if stats_show has been called, but avoiding it would require centralized Sync...
            util::ItemStream::ItemIt it = p_entity->Dbg_GetStats();
            SFR_ASSERT( it == sit );
        }

    static bool Method_Stats_Monitor( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
#ifdef __ENABLE_MONITOR
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        AppTask *pTask = pListener->m_rAppModules.m_pTask;
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        sfr::IView *pVizView = pListener->m_rAppModules.m_pVizView;
        if( 0 == pTask || 0 == pScene )
        {
            SFR_LOG_WARNING( "s2sListener::Method_Stats_Monitor() Requires Task and Scene... ignoring" );
            return true;
        }

        const char *obj_name = params.GetString("name");
        if( 0 == strcmp(obj_name,"*") )
        {
            /*
            for( Scene::EntityIterator it=pScene->GetEntityIterator(); it.IsValid(); ++it )
            {
                util::ItemStream::ItemItRW sit = (*it)->Dbg_QueryStats();
                if( sit.IsValid() )
                {
                    sfr::gui::IWidget *pPTW
                        = pVizView->GetDesktop()->CreateWidget_PropertyTree( (*it)->GetName(),
                                                                             sfr::gui::eDWF_Default,
                                                                             sfr::Vec2(50,50),
                                                                             sit,
                                                                             boost::bind( &s2sListener::OnGetStats,
                                                                                          static_cast<S2::ISyncEntity*>(*it),
                                                                                          _1 ) );
                    pPTW->GetBasicAPI()->SetColor( 0.8,0.0,0.0,0.25 );
                    pPTW->GetBasicAPI()->Minimize();
                }
                else
                    SFR_LOG_WARNING("s2sListener::Method_Stats_Monitor() entity '%s' returned no valid stats, ignoring", (*it)->GetName() );
            }
            return true;
            */
            SFR_LOG_ERROR("s2sListener::Method_Stats_Monitor() wildcard \"*\" unsupported by now..." );
            return false;
        }
        else
        {
            S2::ISyncEntity *pSE = pScene->Find(obj_name);
            if( pSE )
            {
                util::ItemStream::ItemItRW sit = pSE->Dbg_QueryStats();
                if( sit.IsValid() )
                {
                    pTask->GetMonitor().AddStatsMonitor( pSE->GetName(), sit,
                                                         boost::bind( &s2sListener::OnGetStats, pSE, _1 ),
                                                         (std::string(pSE->GetName()) + params.GetString("file_sufix")).c_str() );
                }
                else
                    SFR_LOG_WARNING("s2sListener::Method_Stats_Monitor() entity '%s' returned no valid stats, ignoring", obj_name );
            }
            else return false;
        }
        return true;
#else
        return false;
#endif
    }

    static bool Method_mig2015_AbortNonInteractive( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        if( 0 == pListener->m_rAppModules.m_pTask
            || 0 == pListener->m_rAppModules.m_pTask->GetAppView() )
        {
            SFR_LOG("ABORTING non-interactive run");
            exit(0);
        }
        return true;
    }

    static bool Method_Screenshot( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        if( 0 == pListener->m_rAppModules.m_pTask
            || 0 == pListener->m_rAppModules.m_pTask->GetAppView() )
        {
            SFR_LOG_WARNING( "s2sListener::Method_Screenshot() Requires interactive kernel (test_Saphyre -i true)" );
            return false;
        }
        const char *file_name = params.GetString("file");
        // TODO: Force draw to take proper screnshot: This is
        // IMPOSSIBLE, because all .skr commands are executed
        // sequentially from a single _run(), there's no way to
        // "pause" after a command, redisplay, and continue running
        // the .skr script by now
        // pKernel->ForceDisplay();
        system( (std::string("gnome-screenshot -w -B -f ") + file_name).c_str() );
        return true;
    }

    static bool Method_mig2015_SetContactEntities( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        AppTask *pTask = pListener->m_rAppModules.m_pTask;
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        const char *obj_name1 = params.GetString("name1");
        const char *obj_name2 = params.GetString("name2");
        S2::ISyncEntity *pSE1 = pScene->Find(obj_name1);
        S2::ISyncEntity *pSE2 = pScene->Find(obj_name2);
        if( pSE1 && pSE2 )
        {
            pTask->SetContactEntities3D_MIG2015( pSE1, pSE2 );
            return true;
        }
        else
            return false;
    }

    static bool Method_mig2015_Stats_Monitor( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
#ifdef __ENABLE_MONITOR
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        AppTask *pTask = pListener->m_rAppModules.m_pTask;
        Scene *pScene = pListener->m_rAppModules.m_pScene;
        sfr::IView *pVizView = pListener->m_rAppModules.m_pVizView;
        if( 0 == pTask || 0 == pScene )
        {
            SFR_LOG_WARNING( "s2sListener::Method_mig2015_Stats_Monitor() Requires Task and Scene... ignoring" );
            return true;
        }

        util::ItemStream* pIS_mig2015( new util::ItemStream( 1<<20, 1<<10, util::ItemStream::eRealloc_Identifiers ) ); //IMPORTANT: Data MUST NOT RESIZE or statsmonitor may crash
        util::ItemStream::ItemItRW sit = geo::QueryStats( *pIS_mig2015 );
        if( sit.IsValid() )
        {
            pTask->GetMonitor().AddStatsMonitor( "mig2015 Stats", sit,
                                                 []
                                                 ( util::ItemStream::ItemItRW it )
                                                 { geo::SyncStats( it ); },
                                                 params.GetString("file") );
            geo::np::g_pDefaultContext->m_DCR2DCR_Log_Enabled = false;
            geo::np::g_pDefaultContext->m_DCR2DCR_Viz_Enabled = false;
        }
        else
            SFR_LOG_WARNING("s2sListener::Method_mig2015_Stats_Monitor() returned no valid stats, ignoring" );
        return true;
#else
        return false;
#endif
    }

#ifdef __ENABLE_TEST_ID32
    static bool Method_Id32( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
    {
        s2sListener *pListener = reinterpret_cast<s2sListener*>(user_context.Get<void*>("p_s2s_listener"));
        const char *str_id32 = params.GetString("id");
        util::Id32 id32( str_id32 );
        char id32_to_str[10];
        id32.ToStr(id32_to_str);
        std::cout << "Id32( " << str_id32 << " ) => " << id32_to_str << std::endl; //"error" if !IsValid()
        return true;
    }
#endif
    //@}

    // Other module methods...
    /*
    static bool Method_Monitor( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params );
    static bool Method_Tweak( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params );
    */

    //! Scn commands
    void ExtendLanguageDef_scn( skr::LanguageDef &language_def, const char *prefix )
    {
        //---- time control
        language_def.DefineMethod( prefix, "run","Run simulation (t=duration, non-stop by default)",&Method_Scene_Run)
            .Add("t",double(-1.0f));
        language_def.DefineMethod( prefix, "step","Single Step simulation (t=duration, 1/60 by default, dt = timestep, by default t=dt )",&Method_Scene_Step)
            .Add("t",double(cDefaultTimeStep))
            .Add("dt",double(0));
        //TEMP scn_step(t=5) short cmd
        language_def.DefineMethod( "", "sq","Single Step of 0.25s... changes regularly to extract results for specific scenes",&Method_Scene_Step)
            .Add("t",double(0.25))
            .Add("dt",double(0.0333333));
        language_def.DefineMethod( "", "s1q","Single Step of 1s... changes regularly to extract results for specific scenes",&Method_Scene_Step)
            .Add("t",double(1))
            .Add("dt",double(0.01666666));
        language_def.DefineMethod( "", "s2q","Single Step of 2s... changes regularly to extract results for specific scenes",&Method_Scene_Step)
            .Add("t",double(2))
            .Add("dt",double(0.01666666));
        language_def.DefineMethod( "", "s3q","Single Step of 3s... changes regularly to extract results for specific scenes",&Method_Scene_Step)
            .Add("t",double(3))
            .Add("dt",double(0.01666666));
        //TEMP
        language_def.DefineMethod( prefix, "pause","Pause simulation",&Method_Scene_Pause);
        language_def.DefineMethod( prefix, "reset","Reset simulation (todo)",&Method_Scene_Reset);

        //---- Scene query
        language_def.DefineMethod( prefix, "ls","List scene entities",&Method_Scene_List);

        //---- Scene edition
        double zero_array[4] = {0,0,0,0};
        double default_half_size_array[4] = {1,2,3,4};
        double default_normal_array[3] = {0,1,0};
        double default_3_vertices_2d[3*2] = { 1, 0, 0, 0.5, 0, -0.5 };
        double aabb_0_0_1_1[4] = {0,0,1,1};
        double aabb_0_0_1_05[4] = {0,0,1,0.5};
        double default_solid2d_half_size_array[2] = {2.0,0.5};

        language_def.DefineMethod( prefix, "p3d","Create a Particle3D",&Method_Scene_Particle3D)
            .Add("m",double(1))         //mass
            .Add("r",double(0.1))       //radius
            .AddArray("p",zero_array,3) //pos
            .AddArray("v",zero_array,3) //vel
            .AddString("name","");

        language_def.DefineMethod( prefix, "p2d","Create a Particle2D",&Method_Scene_Particle2D)
            .Add("m",double(1))
            .Add("r",double(0.1))
            .AddArray("p",zero_array,2)
            .AddArray("v",zero_array,2)
            .AddString("name","");

        language_def.DefineMethod( prefix, "k2d_sphere","Create a Kine2D Sphere",&Method_Scene_Kine2D_Sphere)
            .Add("r",double(1.0))        //radius
            .AddArray("p",zero_array,2)  //pos
            .Add("a",double(0))          //angle
            .AddString("name","");

        language_def.DefineMethod( prefix, "k3d_sphere","Create a Kine3D Sphere",&Method_Scene_Kine3D_Sphere)
            .Add("r",double(1.0))        //radius
            .AddArray("p",zero_array,3)  //pos
            // .Add("a",double(0))          //angle
            .AddString("name","");

        language_def.DefineMethod( prefix, "k2d_capsule","Create a Kine2D Capsule",&Method_Scene_Kine2D_Capsule)
            .Add("r",double(1.0))        //radius
            .Add("hh",double(1.0))       //half_height
            .AddArray("p",zero_array,2)  //pos
            .Add("a",double(0))          //angle
            .AddString("name","");

        language_def.DefineMethod( prefix, "k2d_box","Create a Kine2D Box",&Method_Scene_Kine2D_Box)
            .AddArray("h",default_half_size_array,2)  //half_sizes
            .AddArray("p",zero_array,2)  //pos
            .Add("a",double(0))          //angle
            .AddString("name","");

        language_def.DefineMethod( prefix, "k2d_plane","Create a Kine2D Plane",&Method_Scene_Kine2D_Plane)
            .AddArray("n",default_normal_array,2)  //normal
            .Add("d",double(0))                    //coeff_d
            .Add("hs",double(1))                   //half_space
            .AddString("name","");

        language_def.DefineMethod( prefix, "k3d_plane","Create a Kine3D Plane",&Method_Scene_Kine3D_Plane)
            .AddArray("n",default_normal_array,3)  //normal
            .Add("d",double(0))                    //coeff_d
            .Add("hs",double(1))                   //half_space
            .AddString("name","");

        language_def.DefineMethod( prefix, "k2d_poly","Create a Kine2D Polygonal",&Method_Scene_Kine2D_Polygonal)
            .Add("c",double(true))       //closed
            .AddArray("v",default_3_vertices_2d,3*2) //vtx array
            .Add("r",double(0))          //radius
            .AddArray("p",zero_array,2)  //pos
            .Add("a",double(0))          //angle
            .Add("s",double(3))          //subdivision count
            .Add("d",double(0.1))        //distortion
            .AddString("name","")
            .AddString("file","");

        language_def.DefineMethod( prefix, "k2d_path","Create a Kine2D Path",&Method_Scene_Kine2D_Path)
            .Add("c",double(true))       //closed
            .AddArray("v",default_3_vertices_2d,3*2) //vtx array
            .AddArray("p",zero_array,2)  //pos
            .Add("a",double(0))          //angle
            .Add("s",double(3))          //subdivision count
            .Add("d",double(0.1))        //distortion
            .AddString("name","")
            .AddString("file","");

        language_def.DefineMethod( prefix, "k2d_cradle","Create a 2d cradle (3 planes)",&Method_Scene_Kine2D_Cradle)
            /*
            .Add("w",double(10))  // width
            .Add("d",double(5))   // depth
            .Add("a",mal::SixthPi<double>()) //wallangle
            */
            .Add("w",double(2))  // width
            .Add("d",double(1))   // depth
            .Add("a",mal::HalfPi<double>()) //wallangle
            .AddString("name","");

        language_def.DefineMethod( prefix, "k2d_mesh","Create a Kine2D MeshSolid",&Method_Scene_Kine2D_MeshSolid)
            .AddArray("p",zero_array,2)  //pos
            .Add("a",double(0))          //angle
            .AddArray("h",default_half_size_array,2)  //half_sizes
            .Add("nx",double(5))
            .Add("ny",double(5))
            .AddString("name","")
            .AddString("file","");

        language_def.DefineMethod( prefix, "k3d_tetsolid","Create a Kine3D TetSolid",&Method_Scene_Kine3D_TetSolid)
            .AddArray("p",zero_array,3)  //pos
            .Add("ax",double(0))          //angle X
            .Add("ay",double(0))          //angle Y
            .Add("az",double(0))          //angle Z
            .AddArray("h",default_half_size_array,3)  //half_sizes
            .Add("nx",double(5))
            .Add("ny",double(5))
            .Add("nz",double(5))
            .AddString("name","")
            .AddString("file","");

        language_def.DefineMethod( prefix, "ps2d","Create a ParticleSystem2D",&Method_Scene_ParticleSystem2D)
            .Add("np",double(5))
            .Add("m",double(1.0))
            .AddString("name","");

        language_def.DefineMethod( prefix, "f2d","Create a Fluid2D",&Method_Scene_Fluid2D)
            .Add("np",double(1000.0))
            .Add("density",double(1000.0)) // density(water) = 1000kg/m^3
            .Add("thickness", double(0.01)) // fluid sheet thickness in m
            .AddArray("shape_aabb",aabb_0_0_1_05,4)
            .AddArray("bounds_aabb",aabb_0_0_1_1,4)
            .AddString("name","");

        language_def.DefineMethod( prefix, "s2d","Create a Solid2D",&Method_Scene_Solid2D)
            .AddArray("p",zero_array,2)  //pos
            .AddArray("v",default_3_vertices_2d,0) //empty vtx array
            .AddArray("h",default_solid2d_half_size_array,2)  //half_sizes
            .Add("nx",double(9))
            .Add("ny",double(3))
            .Add("subd",double(0))
            .Add("dt", double(cDefaultTimeStep))
            .Add("density",double(1000.0)) // density = kg/m^3
            .Add("thickness", double(0.01)) // z-thickness in m
            .Add("young_modulus",double(1000.0)) // \todo units?
            .Add("poisson_ratio", double(0.25)) // \todo units?
            .Add("damping_ratio", double(1.0))
            .Add("plastic_yield", double(0))
            .Add("plastic_max", double(0))
            .Add("plastic_creep_per_second", double(0))
            .AddString("name","")
            .AddString("file","");

        language_def.DefineMethod( prefix, "s3d","Create a Solid3D",&Method_Scene_Solid3D)
            .AddArray("p",zero_array,3)  //pos
            .Add("ax",double(0))          //angle X
            .Add("ay",double(0))          //angle Y
            .Add("az",double(0))          //angle Z
            //.AddArray("v",default_3_vertices_2d,0) //empty vtx array
            //.AddArray("h",default_solid2d_half_size_array,2)  //half_sizes
            //.Add("nx",double(9))
            //.Add("ny",double(3))
            .Add("subd",double(0))
            .Add("dt", double(cDefaultTimeStep))
            .Add("density",double(1000.0)) // density = kg/m^3
            //.Add("thickness", double(0.01)) // z-thickness in m
            .Add("young_modulus",double(1000.0)) // \todo units?
            .Add("poisson_ratio", double(0.25)) // \todo units?
            .Add("damping_ratio", double(1.0))
            .Add("plastic_yield", double(0))
            .Add("plastic_max", double(0))
            .Add("plastic_creep_per_second", double(0))
            .AddString("name","")
            .AddString("file","");

        language_def.DefineMethod( prefix, "kill","Kill'em All",&Method_Scene_Kill )
            .AddString("name","*");

        language_def.DefineMethod( prefix, "save_state","Save raw entity state",&Method_Scene_SaveState )
            .AddString("name","*")
            .AddString("file","*");

#ifdef __ENABLE_TEST_ID32
        language_def.DefineMethod( "", "id","Parses an id32",&Method_Id32 )
            .AddString("id","abcd#2");
#endif
        // \todo Consider promoting to console command with '_' prefix, as it's always useful
        language_def.DefineMethod( "", "screenshot","Shoot the screen",&Method_Screenshot ).AddString("file","screenshot.png");
    }

    //! Viz commands
    void ExtendLanguageDef_viz( skr::LanguageDef &language_def, const char *prefix )
    {
        language_def.DefineMethod( prefix, "on","Enable Viz",&Method_Viz_On ).AddString("name","*");
        language_def.DefineMethod( prefix, "off","Disable Viz",&Method_Viz_Off ).AddString("name","*");
        language_def.DefineMethod( prefix, "tgl","Toggle Viz",&Method_Viz_Toggle ).AddString("name","*");
    }

    //! Stats commands
    void ExtendLanguageDef_stats( skr::LanguageDef &language_def, const char *prefix )
    {
        language_def.DefineMethod( prefix, "show","Show Stats",&Method_Stats_Show ).AddString("name","*");
        language_def.DefineMethod( prefix, "monitor","Monitor Stats",&Method_Stats_Monitor )
            .AddString("name","*")
            .AddString("magnitude","*")
            .AddString("file_sufix","_stats.out");
    }

    //! Stats commands
    void ExtendLanguageDef_params( skr::LanguageDef &language_def, const char *prefix )
    {
        language_def.DefineMethod( prefix, "show","Show Params",&Method_Params_Show ).AddString("name","*");
    }

    void ExtendLanguageDef_mig2015( skr::LanguageDef &language_def, const char *prefix )
    {
        // MIG2015 commands
        language_def.DefineMethod( prefix, "set_contact_entities","MIG2015",&Method_mig2015_SetContactEntities )
            .AddString("name1","UNNAMED")
            .AddString("name2","UNNAMED");
        language_def.DefineMethod( prefix, "monitor","Monitor mig2015 Stats",&Method_mig2015_Stats_Monitor )
            .AddString("name","*")
            .AddString("magnitude","*")
            .AddString("file","mig2015_stats.out");
        language_def.DefineMethod( prefix, "abort_non_interactive","Abort if not interactive",&Method_mig2015_AbortNonInteractive );
    }

};

#endif //TEST_SAPHYRE_S2S_LISTENER_H
