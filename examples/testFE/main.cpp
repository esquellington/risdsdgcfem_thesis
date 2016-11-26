#include "Config.h"
#include "AppTask.h"
#include "AppInputListener.h"
#include "AppRenderer.h"
#include <Safra/Safra.h>

int main( int argc, const char *argv[] )
{

    //util::LogStream g_LogStream;
    std::string AppName;

#ifdef _DEBUG
    AppName = "testFE_d";
#else
    AppName = "testFE";
#endif

    AppTask app_task;

    bool bUseVizKernel( argc <= 1 || std::string(argv[1]) != "-c" );

    // Parse params
    if( bUseVizKernel )
    {
        Safra::Init(Safra::eKernelGLUT);
        //Safra::SetLogStream( &g_LogStream );

        // Element2D Task
        Safra::AddTask( new sfr::RealtimeUT( &app_task, cDefaultTimeStep ) );

        // Element2D View
        AppRenderer *pAppRenderer
            = new AppRenderer( app_task, Safra::GetDefaultRenderer(), AppRenderer::eDraw_Default );
        sfr::IView *pView = Safra::CreateView( "Element2D View", sfr::IView::eFlags_EnableDesktop, 600, 600 );
        pView->GetCamera()->InitOrtho( sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                       10.0, 1.0,  1.0, 500.0,  600, 600 );
        pView->SetRenderer( pAppRenderer );
        app_task.SetView( pView );

        // Set custom input listener
        AppInputListener *pAppInputListener = new AppInputListener( app_task, pView->GetCamera(),
                                                                    sfr::gui::eCameraController_Spherical );
        pView->SetKeyboardListener( pAppInputListener );
        pView->SetMouseListener( pAppInputListener );
        Safra::AddTask( new sfr::RealtimeUT( pView, -0.05f ) );

        // run lola run
        Safra::Run();

    }
    else
    {
        Safra::Init(Safra::eKernelConsole);
        //Safra::SetLogStream( &g_LogStream );
        Safra::AddTask( new sfr::IdleUT( &app_task, 0.02f ) ); //\todo Actually, should run it just *once* and write output
        Safra::Run();
    }

    // THIS may never be executed... add some kind of "Ending Callback to Safra"
    //g_LogStream.Close();
    Safra::ShutDown();

    return 0;
}
