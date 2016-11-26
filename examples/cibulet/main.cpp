#include "Config.h"
#include "AppTask.h"
#include "AppInputListener.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main( int argc, const char *argv[] )
{
    std::string app_name = "cibulet";
    Params default_params;
    std::string svg_file_name("s2s/svg/shapes2d.svg");
    std::string svg_path_element_name("Patito");
    std::string output_name("cibulet_output");
    // Parse arguments
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("viz,v", "enable interactive visualization")
        ("svg_file,f", po::value<std::string>(&svg_file_name)->default_value(svg_file_name), "SVG input file name" )
        ("svg_path,p", po::value<std::string>(&svg_path_element_name)->default_value(svg_path_element_name), "SVG path element name" )
        ("output,o", po::value<std::string>(&output_name)->default_value(output_name), "output name (=> output_name.txt/.bin/.params)" )
        ("src_size", po::value<float32>(&default_params.m_SRC_NaturalSize)->default_value(default_params.m_SRC_NaturalSize), "Object natural size")
        ("viz_detail", po::value<float32>(&default_params.m_VIZ_Detail)->default_value(default_params.m_VIZ_Detail), "Visualization representation detail")
        ("viz_distortion_s", po::value<float32>(&default_params.m_VIZ_DistortionScale)->default_value(default_params.m_VIZ_DistortionScale), "Visualization distortion scale")
        ("viz_distortion_f", po::value<float32>(&default_params.m_VIZ_DistortionFreq)->default_value(default_params.m_VIZ_DistortionFreq), "Visualization distortion frequency")
        ("ext_detail", po::value<float32>(&default_params.m_EXT_Detail)->default_value(default_params.m_EXT_Detail), "External representation detail")
        ("ext_offset", po::value<float32>(&default_params.m_EXT_Offset)->default_value(default_params.m_EXT_Offset), "External offset (fraction of src_size)")
        ("sim_cdt_size", po::value<float32>(&default_params.m_SIM_CDT_Size)->default_value(default_params.m_SIM_CDT_Size), "Simulation CDT size criteria")
        ("sim_cdt_ratio", po::value<float32>(&default_params.m_SIM_CDT_Ratio)->default_value(default_params.m_SIM_CDT_Ratio), "SimulationCDT ratio criteria")
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
    app_task.Init( svg_file_name.c_str(), svg_path_element_name.c_str(), output_name.c_str(), default_params );

    if( vm.count("viz") == 0 )
    {
        //console tool
        app_task.RebuildAllPhasesAndSave();
    }
    else
    {
        Safra::Init(Safra::eKernelGLUT);
        Safra::AddTask( new sfr::RealtimeUT( &app_task, 0.05f ) );
        // VizStream View
        sfr::gfx::VizRenderer *pVizRenderer = new sfr::gfx::VizRenderer( &app_task.GetVizIS(),
                                                                         Safra::GetDefaultRenderer(),
                                                                         sfr::gfx::VizRenderer::eDraw_All );
        sfr::IView *pView = Safra::CreateView( app_name.c_str(), sfr::IView::eFlags_EnableDesktop, 600, 600 );
        pView->GetCamera()->InitOrtho( sfr::Vec3(0.0,0.0,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                       10.0, 1.0,  1.0, 500.0,  300, 300 );
        pView->SetRenderer( pVizRenderer );
        app_task.SetView(pView);

        // Set custom input listener
        AppInputListener *pAppInputListener = new AppInputListener( app_task, pView->GetCamera(),
                                                                    sfr::gui::eCameraController_Planar );
        pView->SetKeyboardListener( pAppInputListener );
        pView->SetMouseListener( pAppInputListener );
        Safra::AddTask( new sfr::RealtimeUT( pView, -0.05f ) );

        // run lola run
        Safra::Run();
        // THIS may never be executed... add some kind of "Ending Callback to Safra"
        Safra::ShutDown();
    }
    return 0;
}
