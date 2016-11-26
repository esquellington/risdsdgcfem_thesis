#ifndef TEST_SAPHYRE_APP_TASK_H
#define TEST_SAPHYRE_APP_TASK_H

#include "Config.h"
#include "Scene.h"
#include <Safra/Safra.h>
#include <Saphyre2/bs/RayCastQuery.h>
#include <Geo/shape/shape.h> //\todo req to use SphereShape2
//\todo req to test contacts
#include <Geo/np/Contact_X_Y_Stochastic.h>
#include <Geo/mp/DoubleDispatcher.h>

#include <Geo/geo_params.h>
#include <Geo/geo_stats.h>
#include <Mal/GConversion.h>

#define __ENABLE_WIND
#ifdef __ENABLE_WIND
#  include "Wind.h"
#endif

#include "Params.h"

#define __ENABLE_MONITOR
#ifdef __ENABLE_MONITOR
class Monitor
{
public:
    // typedef boost::function<void (util::ItemStream::ItemIt)> D_MonitorSyncCall;
    typedef boost::function<void (util::ItemStream::ItemItRW)> D_MonitorSyncCall;
private:
    class IMonitor
    {
    public:
        inline IMonitor() {}
        virtual ~IMonitor() {}
        virtual void Update( double dt ) = 0;
    };
    class StatsMonitor: public IMonitor
    {
    public:
        // inline StatsMonitor( const char* monitor_name, util::ItemStream::ItemIt sit, D_MonitorSyncCall sync_call, const char* file_name )
        inline StatsMonitor( const char* monitor_name, util::ItemStream::ItemItRW sit, D_MonitorSyncCall sync_call, const char* file_name )
        : m_Name(monitor_name)
        , m_SIT(sit)
        , m_SyncCall(sync_call)
        , m_NumSamples(0)
        , m_Time(0)
            {
                m_FileStream.open( file_name );
                m_FileStream << "# ---- Legend ----" << std::endl;
                int count(1); //First gnuplot column is 1
                //m_FileStream << "#" << count++ << ":sample_id" << std::endl; // Unnecessary
                m_FileStream << "#" << count++ << ":time" << std::endl;
                OutputLegend( m_Name, m_SIT, count );
                m_FileStream << "# Gnuplot Usage: \"gnuplot> plot 'file_name' using 1:C with lines\", where 1 = time column and C = magnitude column index to be plot" << std::endl;
                m_FileStream << "# ---- Data ----" << std::endl;

                Update(0); //First update to guarantee init state gets written
            }
        ~StatsMonitor()
            {
                m_FileStream.close();
            }
        void Update( double dt )
            {
                m_Time += dt; //IMonitor::Update gets called AFTER scene changes
                if( !m_SyncCall.empty() ) m_SyncCall( m_SIT );
                //m_FileStream << m_NumSamples << " "; // Unnecessary
                m_FileStream << m_Time << " ";
                OutputStats( m_SIT );
                m_FileStream << std::endl;
                m_NumSamples++;
            }
    private:
        void OutputLegend( std::string prefix, util::ItemStream::ItemIt it, int& count )
            {
                if( it.IsComplex() )
                {
                    std::string new_prefix = prefix + "." + std::string(it.GetName());
                    for( util::ItemStream::ItemIt it_sub=it.GetSubItem(); it_sub.IsValid(); ++it_sub )
                        OutputLegend( new_prefix, it_sub, count );
                }
                else //\todo Could filter by type, or something here...
                    m_FileStream << "#" << count++ << ":" << prefix << "." << it.GetName() << std::endl;
            }
        void OutputStats( util::ItemStream::ItemIt it )
            {
                if( it.IsComplex() )
                    for( util::ItemStream::ItemIt it_sub=it.GetSubItem(); it_sub.IsValid(); ++it_sub )
                        OutputStats( it_sub );
                else
                {
                    //TEMP: m_FileStream << it.GetName() << " = ";
                    switch( it.GetType() )
                    {
                    case eType_UInt32: m_FileStream << it.Get<uint32>() << " "; break;
                    case eType_Float32: m_FileStream << it.Get<float32>() << " "; break;
                    case eType_Float64: m_FileStream << it.Get<float64>() << " "; break;
                    case eType_Pair_Float32_Float32: m_FileStream << it.Get< GPair<float32,float32> >().m_Second << " "; break;
                    default: m_FileStream << "?!?!"; break;
                    }
                }
            }

    private:
        std::string m_Name;
        // util::ItemStream::ItemIt m_SIT;
        util::ItemStream::ItemItRW m_SIT;
        D_MonitorSyncCall m_SyncCall;
        std::ofstream m_FileStream;
        uint32 m_NumSamples;
        double m_Time;
    };

public:
    inline Monitor() {}
    ~Monitor()
        {
            //\todo std::foreach( m_Monitors.begin(), m_Monitors.end(), []( delete *element ) );
            for( std::vector<IMonitor*>::iterator it_m=m_Monitors.begin(); it_m != m_Monitors.end(); ++it_m )
                delete *it_m;
        }
    void Update( double dt )
        {
            //\todo std::foreach( m_Monitors.begin(), m_Monitors.end(), []( (*element)->Update(dt) ) );
            for( std::vector<IMonitor*>::iterator it_m=m_Monitors.begin(); it_m != m_Monitors.end(); ++it_m )
                (*it_m)->Update(dt);
        }
    void AddStatsMonitor( const char* monitor_name, util::ItemStream::ItemItRW sit, D_MonitorSyncCall d_sync, const char* file_name )
        {
            m_Monitors.push_back( new StatsMonitor(monitor_name,sit,d_sync,file_name) );
        }

private:
    std::vector<IMonitor*> m_Monitors;
};
#endif //__ENABLE_MONITOR

S2::Solid2D *g_pSolid2D(0); //TEMPORAL to apply pressure!!
S2::Solid3D *g_pSolid3D(0); //TEMPORAL to apply pressure!!
const double cDefaultTimeStep(1.0f/60.0f);

class AppTask: public sfr::IUpdateTarget
{
public:
    AppTask()
    : m_Scene( m_AppParams )
    , m_bIsActive(true)
    , m_pAppView(0)
    , m_Dimension(3)//2
    , m_bIsEnabledRayCast(false), m_pRCQ2D(0), m_RayCastPoint(5,5), m_RayCastDir(0,-1)
    , m_bIsEnabledPressure(false), m_PressurePoint(0,0), m_PressureRadius(2.0f)
    , m_bIsEnabledDrag(false), m_DragKPCID(-1), m_DragPoint2D(0,0), m_DragPoint3D(0,0,0), m_DragNodePoint3D(0,0,0), m_DragNodeDir3D(1,0,0)
    , m_bIsEnabledMove(false), m_bIsEnabledRotate(false)
    , m_MoveOrRotatePoint2D(0,0), m_pMoveOrRotateEntity2D(0)
    , m_MoveOrRotatePoint3D(0,0,0), m_MoveOrRotateDir3D(0,0,0), m_MoveOrRotateHitPoint3D(0,0,0), m_pMoveOrRotateEntity3D(0), m_pRCQ3D(0)
    , m_bIsEnabledContact2D(false), m_pContactKine2D_1(0), m_pContactKine2D_2(0)
    , m_bIsEnabledContact3D(false), m_pContactKine3D_1(0), m_pContactKine3D_2(0)
#ifdef __ENABLE_WIND
    , m_bIsEnabledWind(false), m_pWind(0)
#endif
    , m_ParamsIS( 4096, 4096 )
    , m_GeoIS( 1<<15, 1<<15 )
        {
            m_Scene.Init();
            m_pRCQ2D = m_Scene.GetUniverse()->CreateRayCastQuery2D();
            m_pRCQ3D = m_Scene.GetUniverse()->CreateRayCastQuery3D();
#ifdef __ENABLE_WIND
            m_pWind = new Wind2D( m_Scene );
#endif
        }

    void OnSync_AppParams( sfr::gui::PropertyIt pit )
        {
            m_ArchetypeLibrary.SyncInstance( "Archetype_AP", &m_AppParams, pit );//.GetSubItem() );
        }

    void OnSync_GeoParams( sfr::gui::PropertyIt pit )
        {
            //APP_LOG("OnSync_GeoParams BEGIN");
            //std::cout << pit << std::endl;
            geo::SyncParams( pit );
            //APP_LOG("OnSync_GeoParams END");
        }

    void OnSync_GeoStats( sfr::gui::PropertyIt pit )
        {
            //APP_LOG("OnSync_GeoParams BEGIN");
            //std::cout << pit << std::endl;
            geo::SyncStats( pit );
            //APP_LOG("OnSync_GeoParams END");
        }

    void SetView( sfr::IView *p_app_view )
        {
            m_pAppView = p_app_view;
            m_pAppView->GetDesktop()->Message( "Starting....", 2 );

            // Init archetype library
            {
                // Test widget from ExportInstance
                Params::InitArchetype( m_ArchetypeLibrary );
                m_ParamsIS.BeginComplex( "AppParams", eType_Property_Group );
                {
                    m_ArchetypeLibrary.ExportInstance( "Archetype_AP", &m_AppParams, m_ParamsIS );
                }
                m_ParamsIS.EndComplex();
                util::ItemStream::ItemItRW sit = m_ParamsIS.BeginRW().GetSubItem();
                sfr::gui::IWidget *pAppParamsPTW
                    = m_pAppView->GetDesktop()->CreateWidget_PropertyTree( "app",
                                                                           sfr::gui::eDWF_Default,
                                                                           sfr::Vec2(50,50),
                                                                           sit,
                                                                           sfr::gui::D_PropertyTreeSyncCall (
                                                                               std::bind1st( std::mem_fun( &AppTask::OnSync_AppParams ),
                                                                                             this )
                                                                               ) );
                pAppParamsPTW->GetBasicAPI()->Minimize();
            }

            // Geo params
            {
                util::ItemStream::ItemItRW geo_params_it = geo::QueryParams( m_GeoIS );
                sfr::gui::IWidget *pGeoParamsPTW
                    = m_pAppView->GetDesktop()->CreateWidget_PropertyTree( "geo",
                                                                           sfr::gui::eDWF_Default,
                                                                           sfr::Vec2(50,50),
                                                                           geo_params_it,
                                                                           sfr::gui::D_PropertyTreeSyncCall (
                                                                               std::bind1st( std::mem_fun( &AppTask::OnSync_GeoParams ),
                                                                                             this )
                                                                               ) );
                pGeoParamsPTW->GetBasicAPI()->Minimize();
            }

            // Geo stats
            {
                util::ItemStream::ItemItRW geo_stats_it = geo::QueryStats( m_GeoIS );
                sfr::gui::IWidget *pGeoStatsPTW
                    = m_pAppView->GetDesktop()->CreateWidget_PropertyTree( "geo",
                                                                           sfr::gui::eDWF_Default,
                                                                           sfr::Vec2(50,50),
                                                                           geo_stats_it,
                                                                           sfr::gui::D_PropertyTreeSyncCall (
                                                                               std::bind1st( std::mem_fun( &AppTask::OnSync_GeoStats ),
                                                                                             this )
                                                                               ) );
                pGeoStatsPTW->GetBasicAPI()->Minimize();
                //m_GeoIS.SaveTxt("geo_Export.txt"); //TEMPORAL
            }
        }

    inline sfr::IView *GetAppView() const { return m_pAppView; }

    bool Update( double dt )
        {
            if( m_bIsActive ) geo::BeginStats();

            if( GetDimension() == 2 )
            {
                // Raycast
                m_vecRH2D.clear();
                if( m_bIsEnabledRayCast ) RayCast();
                // Either run or force sync to ensure commands get flushed
                if( m_bIsActive )
                {
                    if( m_bIsEnabledPressure ) Pressure();
                    m_Scene.Update(dt);
                }
                /* \todo Consider ForceSync if needed...
                   else
                   m_Scene.GetUniverse()->ForceSync();
                */

                if( m_bIsActive && m_bIsEnabledContact2D )
                {
                    //\todo Kine2D do NOT have a geo::IObject, therefore we use default shape DOF
                    /*
                      bool bContact = geo::np::Contact_IShape2_IShape2_Stochastic( m_pContactKine1->GetShape(),
                      m_pContactKine1->GetTransform(), m_pContactKine1->GetShape()->GetDefaultDOF(),
                      m_pContactKine2->GetShape(),
                      m_pContactKine2->GetTransform(), m_pContactKine2->GetShape()->GetDefaultDOF(),
                      m_ContactData, &m_ContactCache );
                    */
                    m_ContactData2.Clear();
                    geo::mp::g_pDefaultDoubleDispatcher->TestContact( m_pContactKine2D_1->GetShape(),
                                                                      m_pContactKine2D_1->GetTransform(),
                                                                      m_pContactKine2D_1->GetShape()->GetVecDefaultDOF(),
                                                                      m_pContactKine2D_2->GetShape(),
                                                                      m_pContactKine2D_2->GetTransform(),
                                                                      m_pContactKine2D_2->GetShape()->GetVecDefaultDOF(),
                                                                      m_ContactData2, &m_ContactCache2 );
                }
#ifdef __ENABLE_WIND
                // Wind
                if( m_bIsEnabledWind ) m_pWind->Update( dt );
#endif
            }
            else //dim == 3
            {
                if( m_bIsActive )
                {
                    if( m_bIsEnabledPressure ) Pressure();
                    m_Scene.Update(dt);
                }

                if( m_bIsActive && m_bIsEnabledContact3D )
                {
                    m_ContactData3.Clear();
                    geo::mp::g_pDefaultDoubleDispatcher->TestContact( m_pContactKine3D_1->GetShape(),
                                                                      m_pContactKine3D_1->GetTransform(),
                                                                      m_pContactKine3D_1->GetShape()->GetVecDefaultDOF(),
                                                                      m_pContactKine3D_2->GetShape(),
                                                                      m_pContactKine3D_2->GetTransform(),
                                                                      m_pContactKine3D_2->GetShape()->GetVecDefaultDOF(),
                                                                      m_ContactData3, &m_ContactCache3 );
                }
            }

            if( m_bIsActive ) geo::EndStats();

            if( 0 != m_pAppView && 0 != m_pAppView->GetDesktop() ) m_pAppView->GetDesktop()->Update(dt);

#ifdef __ENABLE_MONITOR
            if( m_bIsActive ) m_Monitor.Update(dt);
#endif

            return true;
        }

    inline void SetActive( bool b_active ) { m_bIsActive = b_active; }
    inline bool IsActive() const { return m_bIsActive; }

    inline Scene &GetScene() { return m_Scene; }
    inline const Scene &GetScene() const { return m_Scene; }

    inline Params &GetParams() { return m_AppParams; }
    inline const Params &GetParams() const { return m_AppParams; }


#ifdef __ENABLE_MONITOR
    Monitor &GetMonitor() { return m_Monitor; }
    const Monitor &GetMonitor() const { return m_Monitor; }
#endif

    inline void SetDimension( int d ) { m_Dimension = d; }
    inline int GetDimension() const { return m_Dimension; }

    //\name Drag
    //@{
    //---- 2D
    void BeginDrag2D( const sfr::Vec2 &pos, bool b_nail )
        {
            if( 0 != g_pSolid2D )
            {
                m_bIsEnabledDrag = true;
                m_DragPoint2D = pos;
                SFR_ASSERT( m_DragKPCID == S2::Solid2D::kpc_id_type(-1) );
                g_pSolid2D->Lock();
                m_DragKPCID = g_pSolid2D->CreateKinematicPointConstraint( m_DragPoint2D, m_DragPoint2D, S2::Vec2::Zero() );
                g_pSolid2D->Unlock(true); //we force sync to ensure that KPC EID will be valid before Drag or EndDrag
            }
        }
    void Drag2D( const sfr::Vec2 &pos, bool b_nail )
        {
            if( 0 != g_pSolid2D && m_bIsEnabledDrag )
            {
                m_DragPoint2D = pos;
                g_pSolid2D->Lock();
                g_pSolid2D->SetKinematicPointConstraint_Pos( m_DragKPCID, m_DragPoint2D );
                //g_pSolid2D->SetKinematicPointConstraint_Vel( m_DragKPCID, S2::Vec2(1,1) ); //TEMP testing vel constraints
                g_pSolid2D->Unlock();
            }
        }
    void EndDrag2D( bool b_nail )
        {
            if( 0 != g_pSolid2D && m_bIsEnabledDrag )
            {
                m_bIsEnabledDrag = false;
                if( !b_nail )
                {
                    g_pSolid2D->Lock();
                    g_pSolid2D->DestroyKinematicPointConstraint( m_DragKPCID );
                    g_pSolid2D->Unlock();
                }
                m_DragKPCID = S2::Solid2D::kpc_id_type(-1);
            }
        }
    //---- 3D
    void BeginDrag3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir, bool b_nail )
        {
            if( 0 != g_pSolid3D )
            {
                SFR_ASSERT( m_DragKPCID == S2::Solid3D::kpc_id_type(-1) );
                m_DragPoint3D = pos;
                if( g_pSolid3D->RayCastNode( pos, dir, 0.25f, m_DragNodePoint3D ) )
                {
                    m_bIsEnabledDrag = true;
                    m_DragNodeDir3D = dir;
                    g_pSolid3D->Lock();
                    m_DragKPCID = g_pSolid3D->CreateKinematicPointConstraint( m_DragNodePoint3D, m_DragNodePoint3D, S2::Vec3::Zero() );
                    g_pSolid3D->Unlock(true); //we force sync to ensure that KPC EID will be valid before Drag or EndDrag
                    //APP_LOG("BeginDrag3D kpc_id %d", m_DragKPCID );
                }
            }
        }
    void Drag3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir, bool b_nail )
        {
            if( 0 != g_pSolid3D && m_bIsEnabledDrag )
            {
                //APP_LOG("Drag3D kpc_id %d", m_DragKPCID );
                geo::np::RayHit3 rh;
                if( geo::np::GRayCast_Plane<3>( pos, dir, Intervalf(0,1000),
                                                -m_DragNodeDir3D, -mal::Dot( m_DragNodePoint3D, -m_DragNodeDir3D ),
                                                rh, 0 ) )
                {
                    m_DragNodePoint3D = rh.m_Point;
                    m_DragNodeDir3D = dir;
                    m_DragPoint3D = pos;
                    g_pSolid3D->Lock();
                    g_pSolid3D->SetKinematicPointConstraint_Pos( m_DragKPCID, m_DragNodePoint3D );
                    //g_pSolid3D->SetKinematicPointConstraint_Vel( m_DragKPCID, S2::Vec2(1,1) ); //TEMP testing vel constraints
                    g_pSolid3D->Unlock();
                }
            }
        }
    void EndDrag3D( bool b_nail )
        {
            if( 0 != g_pSolid3D && m_bIsEnabledDrag )
            {
                m_bIsEnabledDrag = false;
                if( !b_nail )
                {
                    g_pSolid3D->Lock();
                    g_pSolid3D->DestroyKinematicPointConstraint( m_DragKPCID );
                    g_pSolid3D->Unlock();
                }
                //APP_LOG("EndDrag3D kpc_id %d", m_DragKPCID );
                m_DragKPCID = S2::Solid3D::kpc_id_type(-1);
            }
        }
    void FixBoundaryPlane3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir )
        {
            if( 0 != g_pSolid3D )
                g_pSolid3D->FixAllNodesInBoundaryPlane( pos, dir, 0.001f );
        }
    //@}

    //\name Move
    //@{
    //---- 2D
    void BeginMove2D( const sfr::Vec2 &pos )
        {
            if( !m_bIsEnabledMove )
            {
                m_pMoveOrRotateEntity2D = TryToSelectEntity2D(pos);
                if( 0 != m_pMoveOrRotateEntity2D )
                {
                    m_bIsEnabledMove = true;
                    m_MoveOrRotatePoint2D = pos;
                    {
                        char str[128];
                        snprintf( str, 128, "BeginMove2D() Selected entity %s", m_pMoveOrRotateEntity2D->GetName() );
                        GetAppView()->GetDesktop()->Message( str, 2 );
                    }
                }
            }
        }
    void BeginRotate2D( const sfr::Vec2 &pos )
        {
            if( !m_bIsEnabledRotate )
            {
                m_pMoveOrRotateEntity2D = TryToSelectEntity2D(pos);
                if( 0 != m_pMoveOrRotateEntity2D )
                {
                    m_bIsEnabledRotate = true;
                    m_MoveOrRotatePoint2D = pos;
                    {
                        char str[128];
                        snprintf( str, 128, "BeginRotate2D() Selected entity %s", m_pMoveOrRotateEntity2D->GetName() );
                        GetAppView()->GetDesktop()->Message( str, 2 );
                    }
                }
            }
        }
    void MoveOrRotate2D( const sfr::Vec2 &pos )
        {
            if( 0 != m_pMoveOrRotateEntity2D
                && (m_bIsEnabledMove || m_bIsEnabledRotate) )
            {
                S2::Vec2 diff( pos - m_MoveOrRotatePoint2D );
                switch( m_pMoveOrRotateEntity2D->GetType() )
                {
                case S2::eEntity_Kine2D:
                    {
                        S2::Kine2D *pK2D( static_cast<S2::Kine2D*>(m_pMoveOrRotateEntity2D) );
                        m_pMoveOrRotateEntity2D->Lock();
                        if( m_bIsEnabledMove )
                            pK2D->SetPos( pK2D->GetPos() + diff );
                        else if( m_bIsEnabledRotate )
                        {
                            S2::Real angle( -diff.x() * 0.25f );
                            pK2D->SetRot( mal::GRotation2x2_From(angle) * pK2D->GetRot() );
                        }
                        m_pMoveOrRotateEntity2D->Unlock();
                    }
                    break;
                case S2::eEntity_Particle2D:
                    {
                        S2::Particle2D *pP2D( static_cast<S2::Particle2D*>(m_pMoveOrRotateEntity2D) );
                        m_pMoveOrRotateEntity2D->Lock();
                        if( m_bIsEnabledMove )
                            pP2D->SetPos( pP2D->GetPos() + diff );
                        m_pMoveOrRotateEntity2D->Unlock();
                    }
                    break;
                    /*
                case S2::eEntity_Solid2D:
                    {
                        S2::Solid2D *pS2D( static_cast<S2::Solid2D*>(m_pMoveOrRotateEntity) );
                        m_pMoveOrRotateEntity->Lock();
                        if( m_bIsEnabledMove )
                            pS2D->SetPosCoM( pS2D->GetPosCoM() + diff ); //\todo GetPosCoM()
                        else if( m_bIsEnabledRotate )
                            pS2D->SetRotCoM( pS2D->GetRotCoM() + diff ); //\todo GetRotCoM()
                        m_pMoveOrRotateEntity->Unlock(true);
                    }
                    break;
                    */
                default: break;
                }
                m_MoveOrRotatePoint2D = pos;
            }
        }
    void EndMove2D( const sfr::Vec2 &pos )
        {
            if( m_bIsEnabledMove )
            {
                m_bIsEnabledMove = false;
                m_pMoveOrRotateEntity2D = 0;
            }
        }
    void EndRotate2D( const sfr::Vec2 &pos )
        {
            if( m_bIsEnabledRotate )
            {
                m_bIsEnabledRotate = false;
                m_pMoveOrRotateEntity2D = 0;
            }
        }
    //---- 3D
    void BeginMove3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir )
        {
            if( !m_bIsEnabledMove )
            {
                m_pMoveOrRotateEntity3D = TryToSelectEntity3D(pos,dir,m_MoveOrRotateHitPoint3D);
                if( 0 != m_pMoveOrRotateEntity3D )
                {
                    m_bIsEnabledMove = true;
                    m_MoveOrRotatePoint3D = pos;
                    m_MoveOrRotateDir3D = dir;
                    {
                        char str[128];
                        snprintf( str, 128, "BeginMove3D() Selected entity %s", m_pMoveOrRotateEntity3D->GetName() );
                        GetAppView()->GetDesktop()->Message( str, 2 );
                    }
                }
            }
        }
    void BeginRotate3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir )
        {
            if( !m_bIsEnabledRotate )
            {
                m_pMoveOrRotateEntity3D = TryToSelectEntity3D(pos,dir,m_MoveOrRotateHitPoint3D);
                if( 0 != m_pMoveOrRotateEntity3D )
                {
                    m_bIsEnabledRotate = true;
                    m_MoveOrRotatePoint3D = pos;
                    m_MoveOrRotateDir3D = dir;
                    {
                        char str[128];
                        snprintf( str, 128, "BeginRotate3D() Selected entity %s", m_pMoveOrRotateEntity3D->GetName() );
                        GetAppView()->GetDesktop()->Message( str, 2 );
                    }
                }
            }
        }
    void MoveOrRotate3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir )
        {
            if( 0 != m_pMoveOrRotateEntity3D
                && (m_bIsEnabledMove || m_bIsEnabledRotate) )
            {
                switch( m_pMoveOrRotateEntity3D->GetType() )
                {
                case S2::eEntity_Kine3D:
                    {
                        geo::np::RayHit3 rh;
                        if( geo::np::GRayCast_Plane<3>( pos, dir, Intervalf(0,1000),
                                                        -m_MoveOrRotateDir3D, -mal::Dot( m_MoveOrRotateHitPoint3D, -m_MoveOrRotateDir3D ),
                                                        rh, 0 ) )
                        {
                            S2::Vec3 diff( rh.m_Point - m_MoveOrRotateHitPoint3D );
                            m_MoveOrRotateHitPoint3D = rh.m_Point;
                            m_MoveOrRotatePoint3D = pos;
                            m_MoveOrRotateDir3D = dir;
                            S2::Kine3D *pK3D( static_cast<S2::Kine3D*>(m_pMoveOrRotateEntity3D) );
                            m_pMoveOrRotateEntity3D->Lock();
                            if( m_bIsEnabledMove )
                                pK3D->SetPos( pK3D->GetPos() + diff );
                            else if( m_bIsEnabledRotate )
                            {
                                S2::Real angle( -diff.x() * 0.25f );
                                //pK3D->SetRot( mal::GRotation3x3_From(Vec3f(0,0,1),angle) * pK3D->GetRot() );
                                pK3D->SetRot( mal::GRotation3x3_From(Vec3f(0,1,0),angle) * pK3D->GetRot() );
                            }
                            m_pMoveOrRotateEntity3D->Unlock();
                        }
                    }
                    break;
                default: break;
                }
            }
        }
    void EndMove3D( const sfr::Vec3 &pos )
        {
            if( m_bIsEnabledMove )
            {
                m_bIsEnabledMove = false;
                m_pMoveOrRotateEntity3D = 0;
            }
        }
    void EndRotate3D( const sfr::Vec3 &pos )
        {
            if( m_bIsEnabledRotate )
            {
                m_bIsEnabledRotate = false;
                m_pMoveOrRotateEntity3D = 0;
            }
        }
    //@}

    //\name Rotate
    //@{
    //@}

#ifdef __ENABLE_WIND
    //\name Wind
    //@{
    void BeginWind( const sfr::Vec2& pos )
        {
            if( !m_bIsEnabledWind )
            {
                m_bIsEnabledWind = true;
                m_pWind->Init( pos, mal::Normalized(Vec2f(1,0.5)), 5.0f, 10.0f,
                               50.0f,
                               50 );
            }
        }
    void Wind( const sfr::Vec2& pos )
        {
            if( m_bIsEnabledWind )
            {
                m_pWind->Init( pos, mal::Normalized(Vec2f(1,0.5)), 5.0f, 10.0f,
                               50.0f,
                               50 );
            }
        }
    void EndWind( const sfr::Vec2& pos )
        {
            if( m_bIsEnabledWind )
            {
                m_bIsEnabledWind = false;
            }
        }
    //@}
#endif

    //\name Contact
    //@{
    void SelectContactEntities2D( const sfr::Vec2 &pos )
        {
            S2::ISyncEntity *pEntity = TryToSelectEntity2D(pos);
            if( 0 == pEntity || S2::eEntity_Kine2D != pEntity->GetType() ) return;

            if( 0 == m_pContactKine2D_1 ) m_pContactKine2D_1 = static_cast<S2::Kine2D*>(pEntity);
            else if( 0 == m_pContactKine2D_2 && static_cast<S2::Kine2D*>(pEntity) != m_pContactKine2D_2 ) m_pContactKine2D_2 = static_cast<S2::Kine2D*>(pEntity);
            else //Unselect 2, select 1, to allow seqüential selection
            {
                m_pContactKine2D_1 = static_cast<S2::Kine2D*>(pEntity);
                m_pContactKine2D_2 = 0;
            }
            m_bIsEnabledContact2D = ( 0 != m_pContactKine2D_1 && 0 != m_pContactKine2D_2 && m_pContactKine2D_1 != m_pContactKine2D_2 );
        }
    void SelectContactEntities3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir )
        {
            sfr::Vec3 tmp;
            S2::ISyncEntity *pEntity = TryToSelectEntity3D(pos,dir,tmp);
            if( 0 == pEntity || S2::eEntity_Kine3D != pEntity->GetType() ) return;

            if( 0 == m_pContactKine3D_1 ) m_pContactKine3D_1 = static_cast<S2::Kine3D*>(pEntity);
            else if( 0 == m_pContactKine3D_2 && static_cast<S2::Kine3D*>(pEntity) != m_pContactKine3D_2 ) m_pContactKine3D_2 = static_cast<S2::Kine3D*>(pEntity);
            else //Unselect 2, select 1, to allow seqüential selection
            {
                m_pContactKine3D_1 = static_cast<S2::Kine3D*>(pEntity);
                m_pContactKine3D_2 = 0;
            }
            m_bIsEnabledContact3D = ( 0 != m_pContactKine3D_1 && 0 != m_pContactKine3D_2 && m_pContactKine3D_1 != m_pContactKine3D_2 );
        }
    void SetContactEntities3D_MIG2015( S2::ISyncEntity* pEntity1, S2::ISyncEntity* pEntity2 )
        {
            if( 0 == pEntity1
                || 0 == pEntity2
                || S2::eEntity_Kine3D != pEntity1->GetType()
                || S2::eEntity_Kine3D != pEntity2->GetType() ) return;
            m_pContactKine3D_1 = static_cast<S2::Kine3D*>(pEntity1);
            m_pContactKine3D_2 = static_cast<S2::Kine3D*>(pEntity2);
            m_bIsEnabledContact3D = ( 0 != m_pContactKine3D_1 && 0 != m_pContactKine3D_2 && m_pContactKine3D_1 != m_pContactKine3D_2 );
            SetDimension(3);
        }
    //@}

    void CreateParticle2D( const sfr::Vec2 &pos )
        {
            // Create a Particle2D using S2 API
            S2::Particle2D *pParticle2D = m_Scene.GetUniverse()->CreateParticle2D();
            pParticle2D->BeginDef();
            {
                S2::Particle2D::Params p2d_params;
                //p2d_params.m_Flags;
                p2d_params.m_Mass = 1.0f;
                p2d_params.m_Radius = 0.5f;
                //p2d_params.m_CoeffRestitution;
                pParticle2D->SetParams(p2d_params);
                pParticle2D->AttachShape( geo::SphereShape2( 0.5f ) );
                //pParticle2D->AttachShape( geo::BoxShape2( S2::Vec2(0.5f,1.0f) ) );

                pParticle2D->SetPos( pos );
                pParticle2D->SetVel( S2::Vec2(0,0) );
            }
            pParticle2D->EndDef(true);
            m_Scene.Add( pParticle2D, "" );
        }

private:
    S2::ISyncEntity *TryToSelectEntity2D( const sfr::Vec2 &pos )
        {
            m_pRCQ2D->Run( pos, Vec2f(0,-1), Intervalf(0,1), 0.0f, S2::eRCQF_Default );
            S2::RayCastQuery2D::Hit2D h2d;
            if( m_pRCQ2D->Closest(h2d) )
                return h2d.m_pEntity;
            else
                return 0;
        }
        S2::ISyncEntity *TryToSelectEntity3D( const sfr::Vec3 &pos, const sfr::Vec3 &dir, sfr::Vec3 &hit_point )
        {
            m_pRCQ3D->Run( pos, dir, Intervalf(0,1000), 0.0f, S2::eRCQF_Default );
            S2::RayCastQuery3D::Hit3D h3d;
            if( m_pRCQ3D->Closest(h3d) )
            {
                hit_point = h3d.m_RayHit.m_Point;
                return h3d.m_pEntity;
            }
            else
                return 0;
        }
    void RayCast()
        {
            m_pRCQ2D->Run( m_RayCastPoint, m_RayCastDir, Intervalf(0,1000), 0.0f, S2::eRCQF_Default );
            //m_pRCQ2D->Run( m_RayCastPoint, m_RayCastDir, Intervalf(-0.001,0.001), 0.0f, S2::RayCastQuery2D::eRCQF_Default );
            S2::RayCastQuery2D::Hit2D h2d;
            if( m_pRCQ2D->Closest(h2d) )
            {
                //std::cout << "RayHit at lambda = " << h2d.m_RayHit2.m_Interval.Min() << std::endl;
                m_vecRH2D.push_back( std::make_pair( h2d.m_pEntity, h2d.m_RayHit ) );
            }
            /*
            else
            std::cout << "No RayHit" << std::endl;
            */
        }
    void Pressure()
        {
            if( 0 != g_pSolid2D )
            {
                g_pSolid2D->Lock();
                g_pSolid2D->ApplyPressure( m_PressurePoint, m_PressureRadius, 10.0f );
                g_pSolid2D->Unlock();
            }
        }

public:
    Params m_AppParams;
    Scene m_Scene;
    bool m_bIsActive;

    sfr::IView *m_pAppView;

private:
    friend class AppInputListener;
    friend class AppRenderer;

    int m_Dimension;

    bool m_bIsEnabledRayCast;
    S2::RayCastQuery2D *m_pRCQ2D;
    Vec2f m_RayCastPoint;
    Vec2f m_RayCastDir;
    std::vector< std::pair< S2::ISyncEntity*, geo::np::RayHit2 > > m_vecRH2D;

    bool m_bIsEnabledPressure;
    Vec2f m_PressurePoint;
    float m_PressureRadius;

    bool m_bIsEnabledDrag, m_bIsEnabledDragNail;
    S2::Solid2D::kpc_id_type m_DragKPCID;
    Vec2f m_DragPoint2D;
    Vec3f m_DragPoint3D;
    Vec3f m_DragNodePoint3D;
    Vec3f m_DragNodeDir3D;

    bool m_bIsEnabledMove;
    bool m_bIsEnabledRotate;
    Vec2f m_MoveOrRotatePoint2D;
    S2::ISyncEntity* m_pMoveOrRotateEntity2D;
    Vec3f m_MoveOrRotatePoint3D;
    Vec3f m_MoveOrRotateDir3D;
    Vec3f m_MoveOrRotateHitPoint3D;
    S2::ISyncEntity* m_pMoveOrRotateEntity3D;
    S2::RayCastQuery3D *m_pRCQ3D;

    bool m_bIsEnabledContact2D;
    S2::Kine2D* m_pContactKine2D_1;
    S2::Kine2D* m_pContactKine2D_2;
    geo::np::ContactData2 m_ContactData2;
    geo::np::ContactCache2 m_ContactCache2;
    bool m_bIsEnabledContact3D;
    S2::Kine3D* m_pContactKine3D_1;
    S2::Kine3D* m_pContactKine3D_2;
    geo::np::ContactData3 m_ContactData3;
    geo::np::ContactCache3 m_ContactCache3;

#ifdef __ENABLE_WIND
    bool m_bIsEnabledWind;
    Wind2D* m_pWind;
#endif

    util::ArchetypeLibrary m_ArchetypeLibrary;
    util::ItemStream m_ParamsIS;
    util::ItemStream m_GeoIS;

#ifdef __ENABLE_MONITOR
    Monitor m_Monitor;
#endif
};

#endif //TEST_SAPHYRE_APP_TASK_H
