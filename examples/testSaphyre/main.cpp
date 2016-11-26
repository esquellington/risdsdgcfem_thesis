#include <Safra/Safra.h>
#include <Safra/gfx/VizRenderer.h>

#include <Safra/task/Console.h>
#include <Safra/task/ProfilerMonitor.h>

#include <iostream>
#include <string>

#include "AppTask.h"
#include "AppRenderer.h"
#include "AppInputListener.h"
#include "s2sListener.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//TEMPORAL: Required
#ifdef BOOST_NO_EXCEPTIONS
#  include <exception>
namespace boost { void throw_exception( std::exception const & e ) {} }
#endif

#define __USE_3D_CAMERA

int main( int argc, const char *argv[] )
{
    util::LogStream g_LogStream;
    std::string AppName;

#ifdef _DEBUG
    AppName = "testSaphyre_d";
#else
    AppName = "testSaphyre";
#endif

    g_LogStream.Open( AppName + ".log.xml", AppName + " Log" );
    S2::BSG::Init( AppName + "_S2BS.log.xml" );

    // Default arguments
    std::string skr_file_name("");
    bool bInteractive(true);
    bool bEnableViz(false);
    bool bEnableConsole(false);
    bool bEnableProfiler(false);

    // Parse arguments
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("interactive,i", po::value<bool>(&bInteractive)->default_value(bInteractive), "enable interactive visualization")
        ("viz,v", po::value<bool>(&bEnableViz)->default_value(bEnableViz), "enable VizRender")
        ("console,c", po::value<bool>(&bEnableConsole)->default_value(bEnableConsole), "enable Console")
        ("profiler,p", po::value<bool>(&bEnableProfiler)->default_value(bEnableProfiler), "enable Profiler")
        ("skr_file,f", po::value<std::string>(&skr_file_name)->default_value(skr_file_name), "Load .skr scene file" )
        ;

    // Parse params and init the AppTask accordingly
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if( vm.count("help") )
    {
        std::cout << desc << "\n";
        return 1;
    }

    AppTask app_task;
    app_task.GetParams().scene.m_bViz = bEnableViz;

    // Parse params
    if( bInteractive )
    {
        Safra::Init(Safra::eKernelGLUT);
        Safra::SetLogStream( &g_LogStream );

        // Scene Task
        Safra::AddTask( new sfr::RealtimeUT( &app_task, cDefaultTimeStep ) );

        // Scene View
        AppRenderer *pAppRenderer
            = new AppRenderer( app_task, Safra::GetDefaultRenderer() );
        sfr::IView *pView = Safra::CreateView( "Scene View",
                                               sfr::IView::eFlags_EnableDesktop | sfr::IView::eFlags_EnableFPS,
                                               600, 600 );
        /*
        pView->GetCamera()->InitOrtho( sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                       10.0, 1.0,  -100.0, 100.0,  600, 600 );
        */
        pView->GetCamera()->InitPerspective( sfr::Vec3(3,4,5.0), sfr::Vec3(0,2,0),
                                             // sfr::Vec3(5,7,10.0), sfr::Vec3(0,0,0),
                                             sfr::Vec3(0,1,0),
                                             45.0, 1.0, 0.1, 1000.0, 600, 600 );
        pView->SetRenderer( pAppRenderer );
        app_task.SetView( pView );

        // Set custom input listener
#ifdef __USE_3D_CAMERA
        AppInputListener *pAppInputListener = new AppInputListener( app_task, pView->GetCamera(), sfr::gui::eCameraController_Spherical );
#else
        AppInputListener *pAppInputListener = new AppInputListener( app_task, pView->GetCamera(), sfr::gui::eCameraController_Planar );
#endif
        pView->SetKeyboardListener( pAppInputListener );
        pView->SetMouseListener( pAppInputListener );
        Safra::AddTask( new sfr::RealtimeUT( pView, -cDefaultTimeStep ) ); //-0.01f ) );

        // VizStream View
        sfr::IView *pVizView(0);
        if(bEnableViz)
        {
            sfr::gfx::VizRenderer *pVizRenderer
                = new sfr::gfx::VizRenderer( S2::BSG::GetVizIS(), Safra::GetDefaultRenderer(), sfr::gfx::VizRenderer::eDraw_All );
            pVizView = Safra::CreateView( "Viz View", sfr::IView::eFlags_EnableDesktop, 600, 600 );
#ifdef __USE_3D_CAMERA
            // Perspective 3D
            // pVizView->GetCamera()->InitPerspective( sfr::Vec3(5,7,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
            //                                         45.0, 1.0, 0.1, 1000.0, 600, 600 );
            //
            // Ortho 3D
               pVizView->GetCamera()->InitOrtho( sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
               10.0, 1.0, -10.0, 500.0, 600, 600 );
            //
            pVizView->SetRenderer( pVizRenderer );
            pVizView->SetMouseListener( new sfr::gui::CameraControllerML( pVizView->GetCamera(), sfr::gui::eCameraController_Spherical ) );
#else
            // Ortho 2D
            pVizView->GetCamera()->InitOrtho( sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                              10.0, 1.0, 0.1, 1000.0,  600, 600 );
            pVizView->SetRenderer( pVizRenderer );
            pVizView->SetMouseListener( new sfr::gui::CameraControllerML( pVizView->GetCamera(), sfr::gui::eCameraController_Planar ) );
#endif

            //IMPORTANT: Viz freq seems to limit App draw freq...
            Safra::AddTask( new sfr::RealtimeUT( pVizView, -0.1f ) );
            //Safra::AddTask( new sfr::RealtimeUT( pVizView, -0.05f ) );
            // Safra::AddTask( new sfr::RealtimeUT( pVizView ) ); //var freq, as fast as it gets
        }

        // Console TextView and Task (sizes related to font width x height)
        //IMPORTANT: Console is EXTREMELY SLOW (TextView is just horrible) and has input problems... consider total rewrite inside a normal view, using Desktop() instead
        sfr::Console* pConsole(0);
        if( bEnableConsole )
        {
            sfr::ITextView *pViewConsole = Safra::CreateTextView( "Console", 0, 8*60, 13*30, 30, 60 );
            pViewConsole->SetFontSize(13); //8x13
            pViewConsole->Write("Consolaaaa't!!");
            pConsole = new sfr::Console( pViewConsole );
            Safra::AddTask( pConsole );
        }

        // Profiler IMPORTANT: The TextView reduces fps to 30, AVOID THIS THIS AT ALL COSTS
        if( bEnableProfiler )
        {
            sfr::ProfilerMonitor *pProfilerMonitor
                = new sfr::ProfilerMonitor( S2::BSG::GetProfIS(),
                                            Safra::CreateTextView( "Profiler", 0,
                                                                   400, 600, 50, 80 ) );
            Safra::AddTask( pProfilerMonitor );
        }

        // Modules and Script Listener
        AppModules app_modules( &app_task,
                                &app_task.GetScene(),
                                pAppRenderer,
                                pVizView,
                                0,//pProfilerMonitor,
                                pConsole );
        s2sListener s2s_listener( app_modules );

        if( !skr_file_name.empty() )
        {
            char str[1024];
            snprintf( str, 1024, "_run( file = \"%s\" );", skr_file_name.c_str() );
            std::cout << "Running script " << skr_file_name << std::endl;
            s2s_listener.Execute( str );
        }

        // run lola run
        Safra::Run();
    }
    else
    {
        Safra::Init(Safra::eKernelConsole);
        Safra::SetLogStream( &g_LogStream );

        Safra::AddTask( new sfr::IdleUT( &app_task, 0.02f ) );

        // Modules and Script Listener
        AppModules app_modules( &app_task, &app_task.GetScene(), 0,0,0,/*0,*/0 );
        s2sListener s2s_listener( app_modules );

        if( !skr_file_name.empty() )
        {
            char str[1024];
            snprintf( str, 1024, "_run( file = \"%s\" );", skr_file_name.c_str() );
            std::cout << "Running script " << skr_file_name << std::endl;
            s2s_listener.Execute( str );
        }

        // run lola run
        Safra::Run();
    }

    // THIS may never be executed... add some kind of "Ending Callback to Safra"
    g_LogStream.Close();
    S2::BSG::ShutDown();
    Safra::ShutDown();

    return 0;
}
