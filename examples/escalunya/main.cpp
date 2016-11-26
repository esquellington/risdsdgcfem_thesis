#include "Config.h"
#include "AppTask.h"
#include "AppInputListener.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main( int argc, const char *argv[] )
{
    std::string app_name = "escalunya";
    Params default_params;
    std::string src_file_name("");
    std::string output_name("escalunya_output");
    // Parse arguments
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("viz,v", "enable interactive visualization")
        ("src_file,f", po::value<std::string>(&src_file_name)->default_value(src_file_name), "SRC input file name" )
        ("output,o", po::value<std::string>(&output_name)->default_value(output_name), "output name (=> output_name.txt/.bin/.params)" )
        ("src_size", po::value<float32>(&default_params.m_SRC_NaturalSize)->default_value(default_params.m_SRC_NaturalSize), "Object natural size")
        ("viz_detail", po::value<float32>(&default_params.m_VIZ_Detail)->default_value(default_params.m_VIZ_Detail), "VIZ detail")
        /* EXT/SIM_CDT params are deprecated... kept for future comparison between CDT and FIT
        ("ext_detail", po::value<float32>(&default_params.m_EXT_Detail)->default_value(default_params.m_EXT_Detail), "EXT detail")
        ("ext_offset", po::value<float32>(&default_params.m_EXT_Offset)->default_value(default_params.m_EXT_Offset), "EXT offset (fraction of src_size)")
        ("sim_cdt_facet_angle", po::value<float32>(&default_params.m_SIM_CDT_Facet_Angle)->default_value(default_params.m_SIM_CDT_Facet_Angle), "SIM CDT facet angle criteria")
        ("sim_cdt_facet_size", po::value<float32>(&default_params.m_SIM_CDT_Facet_Size)->default_value(default_params.m_SIM_CDT_Facet_Size), "SIM CDT facet size criteria")
        ("sim_cdt_facet_dist", po::value<float32>(&default_params.m_SIM_CDT_Facet_Distance)->default_value(default_params.m_SIM_CDT_Facet_Distance), "SIM CDT facet distance criteria")
        ("sim_cdt_cell_size", po::value<float32>(&default_params.m_SIM_CDT_Cell_Size)->default_value(default_params.m_SIM_CDT_Cell_Size), "SIM CDT size criteria")
        ("sim_cdt_cell_ratio", po::value<float32>(&default_params.m_SIM_CDT_Cell_Ratio)->default_value(default_params.m_SIM_CDT_Cell_Ratio), "SIM CDT ratio criteria")
        ("sim_cdt_lloyd", po::value<bool>(&default_params.m_SIM_CDT_Lloyd)->default_value(default_params.m_SIM_CDT_Lloyd), "SIM CDT enable Lloyd")
        ("sim_cdt_odt", po::value<bool>(&default_params.m_SIM_CDT_Odt)->default_value(default_params.m_SIM_CDT_Odt), "SIM CDT enable odt")
        ("sim_cdt_perturb", po::value<bool>(&default_params.m_SIM_CDT_Perturb)->default_value(default_params.m_SIM_CDT_Perturb), "SIM CDT enable perturb")
        ("sim_cdt_exude", po::value<bool>(&default_params.m_SIM_CDT_Exude)->default_value(default_params.m_SIM_CDT_Exude), "SIM CDT enable exude")
        */
        ("sim_cell_size", po::value<float32>(&default_params.m_SIM_CDT_Cell_Size)->default_value(default_params.m_SIM_CDT_Cell_Size), "SIM cell size")
        ("sim_fit_odt", po::value<float32>(&default_params.m_SIM_Fit_ODT_RelaxationCoeff)->default_value(default_params.m_SIM_Fit_ODT_RelaxationCoeff), "SIM fit ODT relaxation coeff")
        ("sim_fit_lpc", po::value<float32>(&default_params.m_SIM_Fit_Lpc_RelaxationCoeff)->default_value(default_params.m_SIM_Fit_Lpc_RelaxationCoeff), "SIM fit Laplacian smoothing \\lambda coeff")
        ("sim_fit_iter", po::value<uint32>(&default_params.m_SIM_Fit_MaxIter)->default_value(default_params.m_SIM_Fit_MaxIter), "SIM fit iterations")
        ("clp_eps_length", po::value<float32>(&default_params.m_CLP_EpsilonLength)->default_value(default_params.m_CLP_EpsilonLength), "CLP close edge epsilon length")
        ("clp_close_holes_iter", po::value<uint32>(&default_params.m_CLP_MaxIter_CloseHoles)->default_value(default_params.m_CLP_MaxIter_CloseHoles), "CLP close holes iteration limit")
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
    app_task.Init( src_file_name, output_name, default_params );

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
        pView->GetCamera()->InitPerspective( sfr::Vec3(5,7,10.0), sfr::Vec3(0,0,0), sfr::Vec3(0,1,0),
                                             45.0, 1.0,  0.1, 1000.0,  600, 600 );
        pView->SetRenderer( pVizRenderer );
        app_task.SetView(pView);

        // Set custom input listener
        AppInputListener *pAppInputListener = new AppInputListener( app_task, pView->GetCamera(),
                                                                    sfr::gui::eCameraController_Spherical );
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
