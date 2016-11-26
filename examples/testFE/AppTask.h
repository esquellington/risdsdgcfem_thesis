#ifndef TEST_FE_APPTASK_H
#define TEST_FE_APPTASK_H

#include "Config.h"
#include "Params.h"
#include "Element2D.h"
#include "Element3D.h"

#include <Safra/Safra.h>

#include <util/Archetype.h>
#include <util/GStatician.h>
#include "SimpleProfiler.h"

#include <iostream>
#include <fstream> //to output plots
#include <string>

#include "DifferentialR.h"

template<typename OStreamT, typename T >
inline OStreamT &operator<<( OStreamT &o_stream, const util::GStatician<T> &statician )
{
    for( unsigned int it_sample=0; it_sample<statician.GetNumSamples(); it_sample++ )
        o_stream << it_sample << " " << statician.GetSample(it_sample) << std::endl;
    return o_stream;
}

class AppTask: public sfr::IUpdateTarget
{
public:
    AppTask()
    : m_Element2D( m_Params )
    , m_Element3D( m_Params )
    , m_bIsActive(true)
    , m_pAppView(0)
    , m_Dimension(2)
    , m_bIsEnabledDrag(false), m_DragPoint(0,0)
    , m_TmpIS( 1<<15, 1<<15 )
        {
            m_Element2D.RebuildFromParams();
            m_Element3D.RebuildFromParams();
        }

    void OnSync_Params( sfr::gui::PropertyIt pit )
        {
            //APP_LOG("OnSync_Params BEGIN");
            //std::cout << pit << std::endl;
            m_ArchetypeLibrary.SyncInstance( "Archetype_SP", &m_Params, pit );//.GetSubItem() );
            /*\todo If RebuildFromParams() was too expensive, we could only do it in this cases...
            if( m_Params.m_YoungModulus != m_Element2D.m_YoungModulus
                || m_Params.m_PoissonRatio != m_Element2D.m_PoissonRatio
                || m_Params.m_Alpha != m_Element2D.m_Alpha
                || m_Params.m_Beta != m_Element2D.m_Beta
                || m_Params.m_Angle0 != m_Element2D.m_Angle0
                || m_Params.m_Dist12 != m_Element2D.m_Dist12 )
                m_Element2D.RebuildFromParams();
            */
            m_Element2D.RebuildFromParams();
            m_Element3D.RebuildFromParams();

            /*TEMP: OLD copied params that do not need RebuildFromParams()
            m_Element2D.m_FactorDetF = m_Params.m_FactorDetF;
            m_Element2D.m_InvertedCompressibilityFactorDetF = m_Params.m_InvertedCompressibilityFactorDetF;
            m_Element2D.m_DegenerateThresholdDetF = m_Params.m_DegenerateThresholdDetF;
            m_Element2D.m_ThresholdIpolDetF = m_Params.m_ThresholdIpolDetF;
            m_Element2D.m_FFM = (Params::EFixFMethod)m_Params.m_FFM;
            m_Element2D.m_ECIE_e_threshold = m_Params.m_ECIE_e_threshold;
            m_Element2D.m_ECIE_k_factor = m_Params.m_ECIE_k_factor;
            */
            //APP_LOG("OnSync_Params END");
            //APP_LOG("Params = %f, %d", m_Params.m_Float32, m_Params.m_Enum32 );
        }

    void SetView( sfr::IView *p_app_view )
        {
            m_pAppView = p_app_view;
            m_pAppView->GetDesktop()->SetBackgroundColor(0.25,0.25,0.25,0.1);
            m_pAppView->GetDesktop()->Message( "Starting....", 2 );
            Params::InitArchetype( m_ArchetypeLibrary );
            m_TmpIS.BeginComplex( "Params", eType_Property_Group );
            {
                m_ArchetypeLibrary.ExportInstance( "Archetype_SP", &m_Params, m_TmpIS );
            }
            m_TmpIS.EndComplex();
            util::ItemStream::ItemItRW sit = m_TmpIS.BeginRW().GetSubItem();
            sfr::gui::IWidget *pParamsPTW
                = m_pAppView->GetDesktop()->CreateWidget_PropertyTree( "app",
                                                                       sfr::gui::eDWF_Default,
                                                                       sfr::Vec2(50,50),
                                                                       sit,
                                                                       sfr::gui::D_PropertyTreeSyncCall (
                                                                           std::bind1st( std::mem_fun( &AppTask::OnSync_Params ),
                                                                                         this )
                                                                           ) );
            pParamsPTW->GetBasicAPI()->Minimize();
            pParamsPTW->GetBasicAPI()->Unroll();
            pParamsPTW->GetBasicAPI()->AddFootButton( "Plot Stats", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotStats), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Plot e(x)", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotEnergyField), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Plot ne(x)", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotEnergyField_Numerical), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Plot f(x)", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotForceField), this), 0 );
            //pParamsPTW->GetBasicAPI()->AddFootButton( "PlotSlices", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotSlices), this), 0 );
            //pParamsPTW->GetBasicAPI()->AddFootButton( "Profile", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_Profile), this), 0 );
            //pParamsPTW->GetBasicAPI()->AddFootButton( "Plot Ellipsoid", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotEllipsoid), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Plot_L&W", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_Plot_Lambda_And_W), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Plot3D", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_Plot3D), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Profile2D", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_Profile2D), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Profile3D", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_Profile3D), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Error DAPD f", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotErrorForceDAPD), this), 0 );
            pParamsPTW->GetBasicAPI()->AddFootButton( "Error df_dx", std::bind1st( std::mem_fun(&AppTask::OnButtonCall_PlotErrorForceJacobian), this), 0 );
        }

    inline sfr::IView *GetAppView() const { return m_pAppView; }

    bool Update( double dt )
        {
            if( m_bIsActive ) ;//m_Element2D.Update(dt);
            if( 0 != m_pAppView && 0 != m_pAppView->GetDesktop() ) m_pAppView->GetDesktop()->Update(dt);
            return true;
        }

    void SetActive( bool b_active ) { m_bIsActive = b_active; }
    bool IsActive() const { return m_bIsActive; }

    int GetElementDimension() const { return m_Dimension; }
    Element2D &GetElement2D() { return m_Element2D; }
    const Element2D &GetElement2D() const { return m_Element2D; }
    Element3D &GetElement3D() { return m_Element3D; }
    const Element3D &GetElement3D() const { return m_Element3D; }

    const Params &GetParams() const { return m_Params; }
    util::SimpleProfiler &GetProfiler() { return m_Profiler; }

private:

    void SetElementDimension( int dimension ) { m_Dimension = dimension; }

    //\name Drag
    //@{
    void BeginDrag( const sfr::Vec2 &pos )
        {
            if( true /*0 != g_pSolid2D*/ )
            {
                m_bIsEnabledDrag = true;
                m_DragPoint = pos;
                /*
                SFR_ASSERT( m_DragKPCID == S2::Solid2D::kpc_id_type(-1) );
                g_pSolid2D->Lock();
                m_DragKPCID = g_pSolid2D->CreateKinematicPointConstraint( m_DragPoint, m_DragPoint, S2::Vec2::Zero() );
                g_pSolid2D->Unlock(true); //we force sync to ensure that KPC EID will be valid before Drag or EndDrag
                */
            }
        }
    void Drag( const sfr::Vec2 &pos )
        {
            if( /*0 != g_pSolid2D && */ m_bIsEnabledDrag )
            {
                m_DragPoint = pos;
                /*
                g_pSolid2D->Lock();
                g_pSolid2D->SetKinematicPointConstraint_Pos( m_DragKPCID, m_DragPoint );
                //g_pSolid2D->SetKinematicPointConstraint_Vel( m_DragKPCID, S2::Vec2(1,1) ); //TEMP testing vel constraints
                g_pSolid2D->Unlock();
                */
            }
        }
    void EndDrag()
        {
            if( /*0 != g_pSolid2D && */ m_bIsEnabledDrag )
            {
                m_bIsEnabledDrag = false;
                /*
                if( !b_nail )
                {
                    g_pSolid2D->Lock();
                    g_pSolid2D->DestroyKinematicPointConstraint( m_DragKPCID );
                    g_pSolid2D->Unlock();
                }
                m_DragKPCID = S2::Solid2D::kpc_id_type(-1);
                */
            }
        }
    //@}

    void OnButtonCall_PlotStats( void *p_user_data ) const
        {
            /*\todo Fer-ho a través d'una pipe enlloc de cridar system()
              FILE *pipe = popen("gnuplot -persist 2>/dev/null", "w");
              fprintf(pipe, "set term x11 enhanced \n");
              fprintf(pipe, "plot x^2 ti 'x^2' with lines\n");
              pclose(pipe);
            */
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            method_name = method_name + std::string( (m_Dimension == 2) ? "_2D" : "_3D" );
            if( m_Params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) )
            {
                std::string file_name( std::string("plot_testFE_TotalEnergy_") + method_name + ".txt" );
                std::ofstream fs;
                fs.open( file_name.c_str() );
                fs << m_TotalEnergyST;
                fs.close();
                std::cout << "TE in " << "[ " << m_TotalEnergyST.m_Min << " < " << m_TotalEnergyST.m_Max << " ]" << " = " << m_TotalEnergyST.m_Avg << " +- " << m_TotalEnergyST.m_Stdev << std::endl;
                system( (std::string("gnuplot -e \"plot '") + file_name + "' with lines; \" -p").c_str() );
            }
            if( m_Params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) )
            {
                std::string file_name( std::string("plot_testFE_ElasticEnergy_") + method_name + ".txt" );
                std::ofstream fs;
                fs.open( file_name.c_str() );
                fs << m_ElasticEnergyST;
                fs.close();
                std::cout << "EE in " << "[ " << m_ElasticEnergyST.m_Min << " < " << m_ElasticEnergyST.m_Max << " ]" << " = " << m_ElasticEnergyST.m_Avg << " +- " << m_ElasticEnergyST.m_Stdev << std::endl;
                system( (std::string("gnuplot -e \"plot '") + file_name + "' with lines; \" -p").c_str() );
            }
            if( m_Params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) )
            {
                std::string file_name( std::string("plot_testFE_KineticEnergy_") + method_name + ".txt" );
                std::ofstream fs;
                fs.open( file_name.c_str() );
                fs << m_KineticEnergyST;
                fs.close();
                std::cout << "KE in " << "[ " << m_KineticEnergyST.m_Min << " < " << m_KineticEnergyST.m_Max << " ]" << " = " << m_KineticEnergyST.m_Avg << " +- " << m_KineticEnergyST.m_Stdev << std::endl;
                system( (std::string("gnuplot -e \"plot '") + file_name + "' with lines; \" -p").c_str() );
            }
            if( m_Params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) )
            {
                std::string file_name( std::string("plot_testFE_AngularMomentum_") + method_name + ".txt" );
                std::ofstream fs;
                fs.open( file_name.c_str() );
                fs << m_AngularMomentumST;
                fs.close();
                std::cout << "L in " << "[ " << m_AngularMomentumST.m_Min << " < " << m_AngularMomentumST.m_Max << " ]" << " = " << m_AngularMomentumST.m_Avg << " +- " << m_AngularMomentumST.m_Stdev << std::endl;
                system( (std::string("gnuplot -e \"plot '") + file_name + "' with lines; \" -p").c_str() );
            }
            if( m_Params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) )
            {
                std::string file_name( std::string("plot_testFE_detF_") + method_name + ".txt" );
                std::ofstream fs;
                fs.open( file_name.c_str() );
                fs << m_DeformationRatioST;
                fs.close();
                std::cout << "det(F) in " << "[ " << m_DeformationRatioST.m_Min << " < " << m_DeformationRatioST.m_Max << " ]" << " = " << m_DeformationRatioST.m_Avg << " +- " << m_DeformationRatioST.m_Stdev << std::endl;
                system( (std::string("gnuplot -e \"plot '") + file_name + "' with lines; \" -p").c_str() );
            }
        }


    void OnButtonCall_PlotEnergyField( void *p_user_data ) const
        {
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            method_name = method_name + std::string( (m_Dimension == 2) ? "_2D" : "_3D" );
            switch( m_Dimension )
            {
            case 2:
                {
                    std::string file_name_xy( std::string("ef_testFE_xy_") + method_name + ".txt" );
                    std::ofstream fs_xy;
                    fs_xy.open( file_name_xy.c_str() );
                    int nid( 0 );
                    Vec2r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY );
                    int hcx( m_Params.m_HalfCellsX );
                    int hcy( m_Params.m_HalfCellsY );
                    Vec2r X[3] = { m_Element2D.m_X[0], m_Element2D.m_X[1], m_Element2D.m_X[2] };
                    Mat2x2r R;
                    Real energy;
                    for( int it_y = -hcy; it_y < hcy; it_y++ )
                    {
                        for( int it_x = -hcx; it_x < hcx; it_x++ )
                        {
                            X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                            m_Element2D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method, &R, &energy );
                            fs_xy << X[nid].x() << " " << X[nid].y() << " " << energy << " " << std::endl;
                        }
                        fs_xy << std::endl;
                    }
                    fs_xy.close();
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy + "' with image; \" -p").c_str() );
                }
                break;
            case 3:
                {
                    //\todo
                }
                break;
            default: break;
            }
        }

    void OnButtonCall_PlotEnergyField_Numerical( void *p_user_data ) const
        {
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            method_name = method_name + std::string( (m_Dimension == 2) ? "_2D" : "_3D" );
            switch( m_Dimension )
            {
            case 2:
                {
                    int nid( m_Params.m_NID );
                    Vec2r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY );
                    int hcx( m_Params.m_HalfCellsX );
                    int hcy( m_Params.m_HalfCellsY );
                    Vec2r X[3] = { m_Element2D.m_X[0], m_Element2D.m_X[1], m_Element2D.m_X[2] };
                    Mat2x2r R;
                    Real energy;
                    // Fill force/energy grids
                    int32 nx( 2*hcx );
                    int32 ny( 2*hcy );
                    int32 grid_size(nx*ny);
                    std::vector<Vec2r> grid_f( grid_size, Vec2r::Zero() );
                    std::vector<Real> grid_e( grid_size, 0 );
                    for( int it_y = -hcy; it_y < hcy; it_y++ )
                    {
                        int i( it_y + ny/2 );
                        for( int it_x = -hcx; it_x < hcx; it_x++ )
                        {
                            int j( it_x + nx/2 );
                            X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                            Vec6r f = m_Element2D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method, &R, &energy );
                            grid_f[i*nx+j] = Vec2r( f[2*nid], f[2*nid+1] );
                            grid_e[i*nx+j] = energy;
                        }
                    }
                    Real dx( half_sizes.x() / hcx );
                    Real dy( half_sizes.y() / hcy );

                    // We will sync energies at reference nid, where energy is expected to be 0
                    int32 ri( (ny/2+(m_Element2D.m_r[nid].y()/dy)) );
                    int32 rj( (nx/2+(m_Element2D.m_r[nid].x()/dx)) );
                    int32 reference_index( ri*nx + rj );

                    //---- Integrate energy from force (along X)
                    std::vector<Real> grid_ne_x( grid_size, 0 );
                    // // Starting from (-hcx,-hcy) dx,dy>0
                    // grid_ne_x[0] = 0; //base value
                    // //Fill i=0
                    // for( int32 j=1; j<nx; j++ )
                    //     grid_ne_x[j] = grid_ne_x[j-1] - 0.5*dx*( grid_f[j-1].x() + grid_f[j].x() );
                    // // Fill i=1..ny-1
                    // for( int32 j=0; j<nx; j++ )
                    //     for( int32 i=1; i<ny; i++ )
                    //         grid_ne_x[i*nx+j] = grid_ne_x[(i-1)*nx+j] - 0.5*dy*( grid_f[(i-1)*nx+j].y() + grid_f[i*nx+j].y() );

                    // // Starting from (hcx,hcy), dx,dy <0
                    // grid_ne_x[grid_size-1] = 0; //base value
                    // // Fill i=last
                    // for( int32 j=nx-2; j>=0; j-- )
                    //     grid_ne_x[(ny-1)*nx+j] = grid_ne_x[(ny-1)*nx+j+1] + 0.5*dx*( grid_f[(ny-1)*nx+j+1].x() + grid_f[(ny-1)*nx+j].x() );
                    // // Fill i=last..first
                    // for( int32 j=nx-1; j>=0; j-- )
                    //     for( int32 i=ny-2; i>=0; i-- )
                    //         grid_ne_x[i*nx+j] = grid_ne_x[(i+1)*nx+j] + 0.5*dy*( grid_f[(i+1)*nx+j].y() + grid_f[i*nx+j].y() );

                    //-- ACCURATE VERSION: Starts from base value at r1=(1,0)
                    // base value
                    grid_ne_x[reference_index] = 0;
                    // Fill i=ri, dx<0
                    for( int32 j=rj-1; j>=0; j-- )
                        grid_ne_x[ri*nx+j] = grid_ne_x[ri*nx+j+1] + 0.5*dx*( grid_f[ri*nx+j+1].x() + grid_f[ri*nx+j].x() );
                    // Fill i=ri, dx>0
                    for( int32 j=rj+1; j<nx; j++ )
                        grid_ne_x[ri*nx+j] = grid_ne_x[ri*nx+j-1] - 0.5*dx*( grid_f[ri*nx+j-1].x() + grid_f[ri*nx+j].x() );
                    // Starting from (ri,rj), dx,dy <0
                    for( int32 j=rj; j>=0; j-- )
                        for( int32 i=ri-1; i>=0; i-- )
                            grid_ne_x[i*nx+j] = grid_ne_x[(i+1)*nx+j] + 0.5*dy*( grid_f[(i+1)*nx+j].y() + grid_f[i*nx+j].y() );
                    // Starting from (ri,rj), dx<0, dy>0
                    for( int32 j=rj; j>=0; j-- )
                        for( int32 i=ri+1; i<ny; i++ )
                            grid_ne_x[i*nx+j] = grid_ne_x[(i-1)*nx+j] - 0.5*dy*( grid_f[(i-1)*nx+j].y() + grid_f[i*nx+j].y() );
                    // Starting from (ri,rj), dx,dy >0
                    for( int32 j=rj; j<nx; j++ )
                        for( int32 i=ri+1; i<ny; i++ )
                            grid_ne_x[i*nx+j] = grid_ne_x[(i-1)*nx+j] - 0.5*dy*( grid_f[(i-1)*nx+j].y() + grid_f[i*nx+j].y() );
                    // Starting from (ri,rj), dx>0, dy<0
                    for( int32 j=rj; j<nx; j++ )
                        for( int32 i=ri-1; i>=0; i-- )
                            grid_ne_x[i*nx+j] = grid_ne_x[(i+1)*nx+j] + 0.5*dy*( grid_f[(i+1)*nx+j].y() + grid_f[i*nx+j].y() );

                    //---- Integrate energy from force (along Y)
                    std::vector<Real> grid_ne_y( grid_size, 0 );
                    // // Starting from (-hcx,-hcy) dx,dy>0
                    // grid_ne_y[0] = 0; //base value
                    // // Fill j=0
                    // for( int32 i=1; i<ny; i++ )
                    //     grid_ne_y[i*nx] = grid_ne_y[(i-1)*nx] - 0.5*dy*( grid_f[(i-1)*nx].y() + grid_f[i*nx].y() );
                    // // Fill j=1..nx-1
                    // for( int32 i=0; i<ny; i++ )
                    //     for( int32 j=1; j<nx; j++ )
                    //         grid_ne_y[i*nx+j] = grid_ne_y[i*nx+j-1] - 0.5*dx*( grid_f[i*nx+j-1].x() + grid_f[i*nx+j].x() );

                    // Starting from (hcx,hcy) dx,dy<0
                    // grid_ne_y[grid_size-1] = 0; //base value
                    // // Fill j=0
                    // for( int32 i=ny-2; i>=0; i-- )
                    //     grid_ne_y[i*nx+(nx-1)] = grid_ne_y[(i+1)*nx+(nx-1)] + 0.5*dy*( grid_f[(i+1)*nx+(nx-1)].y() + grid_f[i*nx+(nx-1)].y() );
                    // // Fill j=1..nx-1
                    // for( int32 i=ny-1; i>=0; i-- )
                    //     for( int32 j=nx-2; j>=0; j-- )
                    //         grid_ne_y[i*nx+j] = grid_ne_y[i*nx+j+1] + 0.5*dx*( grid_f[i*nx+j+1].x() + grid_f[i*nx+j].x() );

                    //-- ACCURATE VERSION: Starts from base value at r1=(1,0)
                    grid_ne_y[reference_index] = 0; //base value
                    // Fill j=rj, dy<0
                    for( int32 i=ri-1; i>=0; i-- )
                        grid_ne_y[i*nx+rj] = grid_ne_y[(i+1)*nx+rj] + 0.5*dy*( grid_f[(i+1)*nx+rj].y() + grid_f[i*nx+rj].y() );
                    // Fill j=rj, dy>0
                    for( int32 i=ri; i<ny; i++ )
                        grid_ne_y[i*nx+rj] = grid_ne_y[(i-1)*nx+rj] - 0.5*dy*( grid_f[(i-1)*nx+rj].y() + grid_f[i*nx+rj].y() );
                    // Starting from (ri,rj), dx,dy <0
                    for( int32 i=ri; i>=0; i-- )
                        for( int32 j=rj-1; j>=0; j-- )
                            grid_ne_y[i*nx+j] = grid_ne_y[i*nx+j+1] + 0.5*dx*( grid_f[i*nx+j+1].x() + grid_f[i*nx+j].x() );
                    // Starting from (ri,rj), dx<0, dy>0
                    for( int32 i=ri; i<ny; i++ )
                        for( int32 j=rj-1; j>=0; j-- )
                            grid_ne_y[i*nx+j] = grid_ne_y[i*nx+j+1] + 0.5*dx*( grid_f[i*nx+j+1].x() + grid_f[i*nx+j].x() );
                    // Starting from (ri,rj), dx,dy >0
                    for( int32 i=ri; i<ny; i++ )
                        for( int32 j=rj+1; j<nx; j++ )
                            grid_ne_y[i*nx+j] = grid_ne_y[i*nx+j-1] - 0.5*dx*( grid_f[i*nx+j-1].x() + grid_f[i*nx+j].x() );
                    // Starting from (ri,rj), dx>0, dy<0
                    for( int32 i=ri; i>=0; i-- )
                        for( int32 j=rj+1; j<nx; j++ )
                            grid_ne_y[i*nx+j] = grid_ne_y[i*nx+j-1] - 0.5*dx*( grid_f[i*nx+j-1].x() + grid_f[i*nx+j].x() );

                    // Equalize reference position energy
                    Real e0( grid_e[reference_index] );
                    Real ne0_x_offset( e0 - grid_ne_x[reference_index] );
                    Real ne0_y_offset( e0 - grid_ne_y[reference_index] );
                    for( int32 i=0; i<ny; i++ )
                    {
                        for( int32 j=0; j<nx; j++ )
                        {
                            grid_ne_x[i*nx+j] += ne0_x_offset;
                            grid_ne_y[i*nx+j] += ne0_y_offset;
                        }
                    }
                    // Plot
                    std::string file_name_xy_nef_x( std::string("nef_x_testFE_xy_") + method_name + ".txt" );
                    std::string file_name_xy_nef_y( std::string("nef_y_testFE_xy_") + method_name + ".txt" );
                    std::string file_name_xy_err_x( std::string("nef_x_testFE_xy_") + method_name + "_ERROR.txt" );
                    std::string file_name_xy_err_y( std::string("nef_y_testFE_xy_") + method_name + "_ERROR.txt" );
                    std::string file_name_xy_diff_xy( std::string("nef_xy_testFE_xy_") + method_name + "_DIFF.txt" );
                    std::ofstream fs_xy_nef_x;
                    std::ofstream fs_xy_nef_y;
                    std::ofstream fs_xy_err_x;
                    std::ofstream fs_xy_err_y;
                    std::ofstream fs_xy_diff_xy;
                    fs_xy_nef_x.open( file_name_xy_nef_x.c_str() );
                    fs_xy_nef_y.open( file_name_xy_nef_y.c_str() );
                    fs_xy_err_x.open( file_name_xy_err_x.c_str() );
                    fs_xy_err_y.open( file_name_xy_err_y.c_str() );
                    fs_xy_diff_xy.open( file_name_xy_diff_xy.c_str() );
                    for( int it_y = -hcy; it_y < hcy; it_y++ )
                    {
                        int i( it_y + ny/2 );
                        for( int it_x = -hcx; it_x < hcx; it_x++ )
                        {
                            int j( it_x + nx/2 );
                            X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );

                            fs_xy_nef_x << X[nid].x() << " " << X[nid].y() << " " << grid_ne_x[i*nx+j] << " " << std::endl;
                            fs_xy_nef_y << X[nid].x() << " " << X[nid].y() << " " << grid_ne_y[i*nx+j] << " " << std::endl;

                            Real err_x( mal::Abs(grid_e[i*nx+j]) > 1e-6 ? mal::Abs(grid_e[i*nx+j]-grid_ne_x[i*nx+j])/mal::Abs(grid_e[i*nx+j]) : 0 );
                            Real err_y( mal::Abs(grid_e[i*nx+j]) > 1e-6 ? mal::Abs(grid_e[i*nx+j]-grid_ne_y[i*nx+j])/mal::Abs(grid_e[i*nx+j]) : 0 );
                            err_x = mal::Min(10.0f,err_x);
                            err_y = mal::Min(10.0f,err_y);
                            fs_xy_err_x << X[nid].x() << " " << X[nid].y() << " " << err_x << " " << std::endl;
                            fs_xy_err_y << X[nid].x() << " " << X[nid].y() << " " << err_y << " " << std::endl;

                            Real diff_xy( (mal::Abs(grid_ne_x[i*nx+j]) > 1e-6 || mal::Abs(grid_ne_y[i*nx+j]) > 1e-6)
                                          ? mal::Abs(grid_ne_x[i*nx+j]-grid_ne_y[i*nx+j])/mal::Max(mal::Abs(grid_ne_x[i*nx+j]),mal::Abs(grid_ne_y[i*nx+j]))
                                          : 0 );
                            fs_xy_diff_xy << X[nid].x() << " " << X[nid].y() << " " << diff_xy << " " << std::endl;
                        }
                        fs_xy_nef_x << std::endl;
                        fs_xy_nef_y << std::endl;
                        fs_xy_err_x << std::endl;
                        fs_xy_err_y << std::endl;
                        fs_xy_diff_xy << std::endl;
                    }
                    fs_xy_nef_x.close();
                    fs_xy_nef_y.close();
                    fs_xy_err_x.close();
                    fs_xy_err_y.close();
                    fs_xy_diff_xy.close();
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy_nef_x + "' with image; \" -p").c_str() );
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy_nef_y + "' with image; \" -p").c_str() );
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy_err_x + "' with image; \" -p").c_str() );
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy_err_y + "' with image; \" -p").c_str() );
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy_diff_xy + "' with image; \" -p").c_str() );
                }
                break;
            case 3:
                {
                    //\todo
                }
                break;
            default: break;
            }
        }

    void OnButtonCall_PlotForceField( void *p_user_data ) const
        {
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            method_name = method_name + std::string( (m_Dimension == 2) ? "_2D" : "_3D" );
            switch( m_Dimension )
            {
            case 2:
                {
                    std::string file_name_xy( std::string("ff_testFE_xy_") + method_name + ".txt" );
                    std::string file_name_xyz( std::string("ff_testFE_xyz_") + method_name + ".txt" );
                    std::ofstream fs_xy, fs_xyz;
                    fs_xy.open( file_name_xy.c_str() );
                    fs_xyz.open( file_name_xyz.c_str() );
                    int nid( 0 );
                    int idx( 0 );
                    Vec2r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY );
                    float max_norm_f( half_sizes.Norm() * m_Params.m_YoungModulus );
                    int hcx( m_Params.m_HalfCellsX );
                    int hcy( m_Params.m_HalfCellsY );
                    Vec2r X[3] = { m_Element2D.m_X[0], m_Element2D.m_X[1], m_Element2D.m_X[2] };
                    float length( 0.8f * half_sizes[0] / hcx );
                    for( int it_y = -hcy; it_y < hcy; it_y++ )
                    {
                        for( int it_x = -hcx; it_x < hcx; it_x++ )
                        {
                            X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                            Vec6r f( m_Element2D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method ) );
                            Vec2r f_i( f[2*idx], f[2*idx+1] );
                            fs_xy << X[nid].x() << " " << X[nid].y() << " " << f_i.x() << " " << f_i.y() << " " << std::endl;
                            float c_i( f_i.Norm() / max_norm_f );
                            Vec2r v_i( length * mal::SafeNormalized( f_i ) );
                            fs_xyz << X[nid].x() << " " << X[nid].y() << " " << c_i << " " << v_i.x() << " " << v_i.y() << " " << 0 << std::endl;
                        }
                        fs_xy << std::endl;
                        fs_xyz << std::endl;
                    }
                    fs_xy.close();
                    fs_xyz.close();
                    system( (std::string("gnuplot -e \"plot '") + file_name_xy + "' with vectors; \" -p").c_str() );
                    system( (std::string("gnuplot -e \"splot '") + file_name_xyz + "' with vectors; \" -p").c_str() );
                }
                break;
            case 3:
                {
                    int nid( 0 );
                    int idx( 0 );
                    Vec3r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY, m_Params.m_HalfSizeZ );
                    int hcx( m_Params.m_HalfCellsX );
                    int hcy( m_Params.m_HalfCellsY );
                    int hcz( m_Params.m_HalfCellsZ );
                    float length( 0.8f * half_sizes[0] / hcx );
                    float max_norm_f( half_sizes.Norm() * m_Params.m_YoungModulus );
                    Vec3r X[4] = { m_Element3D.m_X[0], m_Element3D.m_X[1], m_Element3D.m_X[2], m_Element3D.m_X[3] };

                    std::string file_name_zyx( std::string("ff_testFE_zyx_") + method_name + ".txt" );
                    std::ofstream fs_zyx;
                    fs_zyx.open( file_name_zyx.c_str() );
                    // ZXY
                    for( int it_z = m_Params.m_SliceZ0; it_z <= m_Params.m_SliceZ1; it_z++ )
                    {
                        for( int it_y = m_Params.m_SliceY0; it_y <= m_Params.m_SliceY1; it_y++ )
                        {
                            for( int it_x = m_Params.m_SliceX0; it_x <= m_Params.m_SliceX1; it_x++ )
                            {
                                X[nid] = Vec3r( half_sizes[0] * Real(it_x) / hcx,
                                                half_sizes[1] * Real(it_y) / hcy,
                                                half_sizes[2] * Real(it_z) / hcz );
                                Vec12r f( m_Element3D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method ) );
                                Vec3r f_i( f[3*idx], f[3*idx+1], f[3*idx+2] );
                                fs_zyx << X[nid].x() << " " << X[nid].y() << " " << X[nid].z() << " " << f_i.x() << " " << f_i.y() << " " << f_i.z() << std::endl;
                            }
                        }
                    }
                    fs_zyx.close();
                    system( (std::string("gnuplot -e \"splot '") + file_name_zyx + "' with vectors; \" -p").c_str() );
                }
                break;
            default: break;
            }
        }
    /*
    void OnButtonCall_PlotSlices( void *p_user_data ) const
        {
            if( 3 != m_Dimension ) return;
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            method_name = method_name + std::string( "_3D" );
            int nid( 0 );
            int idx( 0 );
            Vec3r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY, m_Params.m_HalfSizeZ );
            int hcx( m_Params.m_HalfCellsX );
            int hcy( m_Params.m_HalfCellsY );
            int hcz( m_Params.m_HalfCellsZ );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * m_Params.m_YoungModulus );
            Vec3r X[4] = { m_Element3D.m_X[0], m_Element3D.m_X[1], m_Element3D.m_X[2], m_Element3D.m_X[3] };
            // Slice X0
            std::string file_name_x0( std::string("ff_testFE_xyz_X0_") + method_name + ".txt" );
            std::ofstream fs_x0;
            fs_x0.open( file_name_x0.c_str() );
            for( int it_z = m_Params.m_SliceZ0; it_z <= m_Params.m_SliceZ1; it_z++ )
            {
                for( int it_y = m_Params.m_SliceY0; it_y <= m_Params.m_SliceY1; it_y++ )
                {
                    X[nid] = Vec3r( 0, //x=0
                                    half_sizes[1] * Real(it_y) / hcy,
                                    half_sizes[2] * Real(it_z) / hcz );
                    Vec12r f( m_Element3D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method ) );
                    Vec3r f_i( f[3*idx], f[3*idx+1], f[3*idx+2] );
                    fs_x0 << X[nid].y() << " " << X[nid].z() << " " << f_i.x() << " " << f_i.y() << " " << f_i.z() << std::endl;
                }
            }
            fs_x0.close();



            std::string file_name_x0( std::string("ff_testFE_xyz_X0_") + method_name + ".txt" );
            std::ofstream fs_x0;
            fs_x0.open( file_name_x0.c_str() );
            for( int it_z = m_Params.m_SliceZ0; it_z <= m_Params.m_SliceZ1; it_z++ )
            {
                for( int it_y = m_Params.m_SliceY0; it_y <= m_Params.m_SliceY1; it_y++ )
                {
                    X[nid] = Vec3r( 0, //x=0
                                    half_sizes[1] * Real(it_y) / hcy,
                                    half_sizes[2] * Real(it_z) / hcz );
                    Vec12r f( m_Element3D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method ) );
                    Vec3r f_i( f[3*idx], f[3*idx+1], f[3*idx+2] );
                    fs_x0 << X[nid].y() << " " << X[nid].z() << " " << f_i.x() << " " << f_i.y() << " " << f_i.z() << std::endl;
                }
            }
            fs_x0.close();

            system( (std::string("gnuplot -e \"splot '") + file_name_x0 + "' with vectors; \" -p").c_str() );
            system( (std::string("gnuplot -e \"plot '") + file_name_y0 + "' with vectors; \" -p").c_str() );
            system( (std::string("gnuplot -e \"plot '") + file_name_z0 + "' with vectors; \" -p").c_str() );
        }
    */

    void OnButtonCall_Plot3D( void *p_user_data ) const
        {
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            method_name = method_name + "_3D";
            int nid( 0 );
            int idx( 0 );
            Vec3r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY, m_Params.m_HalfSizeZ );
            int hcx( m_Params.m_HalfCellsX );
            int hcy( m_Params.m_HalfCellsY );
            int hcz( m_Params.m_HalfCellsZ );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( 0 );
            Vec3r X[4] = { m_Element3D.m_X[0], m_Element3D.m_X[1], m_Element3D.m_X[2], m_Element3D.m_X[3] };

            // Output X Y Z fx fy fz
            std::string file_name_zyx( std::string("ff_testFE_zyx_") + method_name + ".txt" );
            std::ofstream fs_zyx;
            fs_zyx.open( file_name_zyx.c_str() );
            for( int it_z = m_Params.m_SliceZ0; it_z <= m_Params.m_SliceZ1; it_z++ )
            {
                for( int it_y = m_Params.m_SliceY0; it_y <= m_Params.m_SliceY1; it_y++ )
                {
                    for( int it_x = m_Params.m_SliceX0; it_x <= m_Params.m_SliceX1; it_x++ )
                    {
                        X[nid] = Vec3r( half_sizes[0] * Real(it_x) / hcx,
                                        half_sizes[1] * Real(it_y) / hcy,
                                        half_sizes[2] * Real(it_z) / hcz );
                        Vec12r f( m_Element3D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method ) );
                        Vec3r f_i( f[3*idx], f[3*idx+1], f[3*idx+2] );
                        if( f_i.Norm() > max_norm_f ) max_norm_f = f_i.Norm();
                        fs_zyx << X[nid].x() << " " << X[nid].y() << " " << X[nid].z() << " " << f_i.x() << " " << f_i.y() << " " << f_i.z() << std::endl;
                    }
                }
            }
            fs_zyx.close();
            // Output gnuplot commands
            std::ofstream fs_gnp;
            fs_gnp.open( "plot3d.gnp" );
            {
                fs_gnp << "reset" << std::endl;
                fs_gnp << "set terminal wxt size 350,262 enhanced font 'Verdana,10' persist" << std::endl;
                fs_gnp << "set border 0" << std::endl;
                fs_gnp << "vector_length(x,y,z) = sqrt( x**2 + y**2 + z**2 )" << std::endl;
                fs_gnp << "scaling = .5" << std::endl;
                fs_gnp << "dx(x,y,z) = scaling * x/"<< max_norm_f << std::endl; //#/vector_length(x,y,z)" << std::endl;
                fs_gnp << "dy(x,y,z) = scaling * y/"<< max_norm_f << std::endl; //#/vector_length(x,y,z)" << std::endl;
                fs_gnp << "dz(x,y,z) = scaling * z/"<< max_norm_f << std::endl; //#/vector_length(x,y,z)" << std::endl;
                fs_gnp << "compute_color(x,y,z,vx,vy,vz) = vector_length(vx,vy,vz)" << std::endl;
                //fs_gnp << "compute_color(x,y,z,vx,vy,vz) = 1-abs(z)/" << half_sizes[2]*m_Params.m_SliceZ0/hcz << std::endl; //Gradient slices in Z, this works GREAT for an orthogonal view along Z
                //fs_gnp << "compute_color(x,y,z,vx,vy,vz) = vector_length(x+1,y,z)" << std::endl; //Gradient centered at the critical point

                fs_gnp << "set object 1 polygon from "
                       << m_Element3D.m_r[0][0] << "," << m_Element3D.m_r[0][1] << "," << m_Element3D.m_r[0][2] << " to "
                       << m_Element3D.m_r[1][0] << "," << m_Element3D.m_r[1][1] << "," << m_Element3D.m_r[1][2] << " to "
                       << m_Element3D.m_r[2][0] << "," << m_Element3D.m_r[2][1] << "," << m_Element3D.m_r[2][2] << " to "
                       << m_Element3D.m_r[3][0] << "," << m_Element3D.m_r[3][1] << "," << m_Element3D.m_r[3][2] << " to "
                       << m_Element3D.m_r[0][0] << "," << m_Element3D.m_r[0][1] << "," << m_Element3D.m_r[0][2] << " to "
                       << m_Element3D.m_r[2][0] << "," << m_Element3D.m_r[2][1] << "," << m_Element3D.m_r[2][2] << " to "
                       << m_Element3D.m_r[3][0] << "," << m_Element3D.m_r[3][1] << "," << m_Element3D.m_r[3][2] << " to "
                       << m_Element3D.m_r[1][0] << "," << m_Element3D.m_r[1][1] << "," << m_Element3D.m_r[1][2] << std::endl;
                fs_gnp << "set object 1 fc rgb \"black\" fillstyle solid 0.5 lw 2 front" << std::endl;

                //fs_gnp << "set palette defined ( 0 '#0000ff', 1 '#00ee22', 2 '#777700', 3 '#ffee00', 4 '#ff000000')" << std::endl;
                fs_gnp << "set palette defined ( 0 '#ffffff', 1 '#ffee00', 2 '#ff7000', 3 '#ee0000', 4 '#7f0000')" << std::endl;
                fs_gnp << "set xrange [" << half_sizes[0]*m_Params.m_SliceX0/hcx-0.1 << ":" << half_sizes[0]*m_Params.m_SliceX1/hcx+0.1 << "]" << std::endl;
                fs_gnp << "set yrange [" << half_sizes[1]*m_Params.m_SliceY0/hcy-0.1 << ":" << half_sizes[1]*m_Params.m_SliceY1/hcy+0.1 << "]" << std::endl;
                fs_gnp << "set zrange [" << half_sizes[2]*m_Params.m_SliceZ0/hcz-0.1 << ":" << half_sizes[2]*m_Params.m_SliceZ1/hcz+0.1 << "]" << std::endl;
                fs_gnp << "splot '" << file_name_zyx << "' using 1:2:3:(dx($4,$5,$6)):(dy($4,$5,$6)):(dz($4,$5,$6)):(compute_color($1,$2,$3,$4,$5,$6)) with vectors head size 0.1,20,60 filled lc palette lw 3" << std::endl;
                /*
                fs_gnp << "set palette defined ( 0 '#ff000000', 1 '#ffee00', 2 '#777700', 3 '#00ee22', 4 '#0000ff')" << std::endl;
                fs_gnp << "set zrange [ -0.1 : 0.1 ]" << std::endl;
                fs_gnp << "replot" << std::endl;
                */
            }
            fs_gnp.close();
            // Run gnuplot commands
            system( "gnuplot plot3d.gnp -p" );
        }

    void OnButtonCall_PlotEllipsoid( void *p_user_data ) const
        {
            std::string file_name_xyz( "ff_testFE_xyz_Ellipsoid3D.txt" );
            std::ofstream fs_xyz;
            fs_xyz.open( file_name_xyz.c_str() );
            Vec3r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY, m_Params.m_HalfSizeZ );
            int hcx( m_Params.m_HalfCellsX );
            int hcy( m_Params.m_HalfCellsY );
            int hcz( m_Params.m_HalfCellsZ );
            for( int it_z = m_Params.m_SliceZ0; it_z <= m_Params.m_SliceZ1; it_z++ )
            {
                for( int it_y = m_Params.m_SliceY0; it_y <= m_Params.m_SliceY1; it_y++ )
                {
                    for( int it_x = m_Params.m_SliceX0; it_x <= m_Params.m_SliceX1; it_x++ )
                    {
                        Vec3r pos( half_sizes[0] * Real(it_x) / hcx,
                                   half_sizes[1] * Real(it_y) / hcy,
                                   half_sizes[2] * Real(it_z) / hcz );
                        Vec3r f( pos.x(), pos.y(), 0 );
                        fs_xyz << pos.x() << " " << pos.y() << " " << pos.z() << " " << f.x() << " " << f.y() << " " << f.z() << std::endl;
                    }
                }
            }
            fs_xyz.close();
            system( (std::string("gnuplot -e \"splot '") + file_name_xyz + "' with vectors; \" -p").c_str() );
        }

    void OnButtonCall_Profile2D( void *p_user_data )
        {
            m_Profiler.Clear();
            m_Profiler.Begin("Trajectory");
            const unsigned int cNumProfileSamples(1000);
            for( unsigned int it_sample=0; it_sample<cNumProfileSamples; it_sample++ )
            {
                Vec2r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY );
                Mat2x2r R;
                Vec2r X[3] = { m_Element2D.m_X[0], m_Element2D.m_X[1], m_Element2D.m_X[2] };
                X[0] = Vec2r( m_DragPoint );
                Vec2r prevX[3] = { X[0], X[1], X[2] };
                Vec2r velX[3] = { Vec2r::Zero(), Vec2r::Zero(), Vec2r::Zero() };
                Real t(0);
                Real inv_node_mass( Real(3)/m_Params.m_Mass );
                Real damping_coeff( 0 );//m_Params.m_DampingRatio * 2 * mal::Sqrt( m_Params.m_YoungModulus * m_Params.m_Mass ) );
                // Init degeneration tracking and NoC
                Real det_F( mal::Det( m_Element2D.Compute_F( X, Params::eFFM_None ) ) );
                bool bDegenerate( det_F <= m_Params.m_DegenerateThresholdDetF );
                int trajectory_noc( bDegenerate ? m_Params.m_NoC : -1 ); //Assume NoC = 0 if starts degenerate
                do
                {
                    Vec6r f( m_Element2D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method, &R, 0 ) );
                    // Draw elastic force
                    Vec2r f0( f[0], f[1] );
                    Vec2r f1( f[2], f[3] );
                    Vec2r f2( f[4], f[5] );
                    if( m_Params.m_TrajectoryFlags.Test(Params::eTF_Free_x0) )
                    {
                        velX[0] += m_Params.m_TimeStep * inv_node_mass * (f0 - damping_coeff * velX[0]);
                        X[0] += m_Params.m_TimeStep * velX[0];
                    }
                    if( m_Params.m_TrajectoryFlags.Test(Params::eTF_Free_x1) )
                    {
                        velX[1] += m_Params.m_TimeStep * inv_node_mass * (f1 - damping_coeff * velX[1]);
                        X[1] += m_Params.m_TimeStep * velX[1];
                    }
                    if( m_Params.m_TrajectoryFlags.Test(Params::eTF_Free_x2) )
                    {
                        velX[2] += m_Params.m_TimeStep * inv_node_mass * (f2 - damping_coeff * velX[2]);
                        X[2] += m_Params.m_TimeStep * velX[2];
                    }
                    t += m_Params.m_TimeStep;
                    if( Params::eMethod_PD_Project == m_Params.m_Method
                        || Params::eMethod_PD_Reflect == m_Params.m_Method
                        || Params::eMethod_PD_Fix == m_Params.m_Method
                        || Params::eMethod_PD_Unified == m_Params.m_Method
                        || Params::eMethod_DAPD_NS == m_Params.m_Method
                        || Params::eMethod_DAPD_SS == m_Params.m_Method
                        || Params::eMethod_DAPD_EX == m_Params.m_Method )
                    {
                        // Compute ToC->NoC
                        det_F = mal::Det( m_Element2D.Compute_F( X, Params::eFFM_None ) );
                        if( !bDegenerate && det_F <= m_Params.m_DegenerateThresholdDetF )
                        {
                            bDegenerate = true;
                            if( m_Params.m_bEnableNoC ) trajectory_noc = ComputeNoC( prevX, X, m_Params.m_DegenerateThresholdDetF, m_Element2D.m_Area );
                            else trajectory_noc = m_Params.m_NoC;
                        }
                        else if( det_F > m_Params.m_DegenerateThresholdDetF )
                        {
                            bDegenerate = false;
                            trajectory_noc = -1;
                        }
                        m_Element2D.m_NoC = trajectory_noc;
                    }
                    // Advance
                    prevX[0] = X[0];
                    prevX[1] = X[1];
                    prevX[2] = X[2];
                } while( t < m_Params.m_Duration );
            }
            m_Profiler.End(); //"Trajectory"
            m_Profiler.Print();

            m_Element2D.m_NoC = m_Params.m_NoC;
        }

    void OnButtonCall_Profile3D( void *p_user_data )
        {
            m_Profiler.Clear();
            m_Profiler.Begin("Trajectory");
            const unsigned int cNumProfileSamples(1000);
            for( unsigned int it_sample=0; it_sample<cNumProfileSamples; it_sample++ )
            {
                Vec3r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY, m_Params.m_HalfSizeZ );
                Mat3x3r R;
                Vec3r X[4] = { m_Element3D.m_X[0], m_Element3D.m_X[1], m_Element3D.m_X[2], m_Element3D.m_X[3] };
                X[0] = Vec3r( m_DragPoint.x(), m_DragPoint.y(), m_Params.m_DragPointZ * m_Params.m_HalfSizeZ );
                Vec3r prevX[4] = { X[0], X[1], X[2], X[3] };
                Vec3r velX[4] = { Vec3r::Zero(), Vec3r::Zero(), Vec3r::Zero(), Vec3r::Zero() };
                Real t(0);
                Real inv_node_mass( Real(4)/m_Params.m_Mass );
                Real damping_coeff( 0 );//m_Params.m_DampingRatio * 2 * mal::Sqrt( m_Params.m_YoungModulus * m_Params.m_Mass ) );
                // Init degeneration tracking and NoC
                Real det_F( mal::Det( m_Element3D.Compute_F( X, Params::eFFM_None ) ) );
                bool bDegenerate( det_F <= m_Params.m_DegenerateThresholdDetF );
                int trajectory_noc( bDegenerate ? m_Params.m_NoC : -1 ); //Assume NoC = 0 if starts degenerate
                do
                {
                    Vec12r f( m_Element3D.ComputeElasticForce( X, (Params::EMethod)m_Params.m_Method, &R, 0 ) );
                    // Draw elastic force
                    Vec3r f0( mal::GRange<0,2>(f) );
                    Vec3r f1( mal::GRange<3,5>(f) );
                    Vec3r f2( mal::GRange<6,8>(f) );
                    Vec3r f3( mal::GRange<9,11>(f) );
                    if( m_Params.m_TrajectoryFlags.Test(Params::eTF_Free_x0) )
                    {
                        velX[0] += m_Params.m_TimeStep * inv_node_mass * (f0 - damping_coeff * velX[0]);
                        X[0] += m_Params.m_TimeStep * velX[0];
                    }
                    if( m_Params.m_TrajectoryFlags.Test(Params::eTF_Free_x1) )
                    {
                        velX[1] += m_Params.m_TimeStep * inv_node_mass * (f1 - damping_coeff * velX[1]);
                        X[1] += m_Params.m_TimeStep * velX[1];
                    }
                    if( m_Params.m_TrajectoryFlags.Test(Params::eTF_Free_x2) )
                    {
                        velX[2] += m_Params.m_TimeStep * inv_node_mass * (f2 - damping_coeff * velX[2]);
                        X[2] += m_Params.m_TimeStep * velX[2];
                    }
                    t += m_Params.m_TimeStep;
                    if( Params::eMethod_PD_Project == m_Params.m_Method
                        || Params::eMethod_PD_Reflect == m_Params.m_Method
                        || Params::eMethod_PD_Fix == m_Params.m_Method )
                    {
                        // Compute ToC->NoC
                        det_F = mal::Det( m_Element3D.Compute_F( X, Params::eFFM_None ) );
                        if( !bDegenerate && det_F <= m_Params.m_DegenerateThresholdDetF )
                        {
                            bDegenerate = true;
                            if( m_Params.m_bEnableNoC )
                            {
                                DoC doc = ComputeDoC( prevX[0], prevX[1], prevX[2], prevX[3],
                                                      X[0], X[1], X[2], X[3],
                                                      m_Params.m_DegenerateThresholdDetF,
                                                      m_Element3D.m_Volume );
                                if( doc.IsVF() ) trajectory_noc = doc.m_VIT0;
                                else trajectory_noc = -1;
                            }
                            else trajectory_noc = m_Params.m_NoC;
                        }
                        else if( det_F > m_Params.m_DegenerateThresholdDetF )
                        {
                            bDegenerate = false;
                            trajectory_noc = -1;
                        }
                        m_Element3D.m_NoC = trajectory_noc;
                    }
                    // Advance
                    prevX[0] = X[0];
                    prevX[1] = X[1];
                    prevX[2] = X[2];
                    prevX[3] = X[3];
                } while( t < m_Params.m_Duration );
            }
            m_Profiler.End(); //"Trajectory"
            m_Profiler.Print();

            m_Element3D.m_NoC = m_Params.m_NoC;
        }

    void OnButtonCall_Plot_Lambda_And_W( void *p_user_data ) const
        {
            std::string method_name( m_Params.GetMethodName( (Params::EMethod)m_Params.m_Method ) );
            std::string file_name_lambda( "ff_testFE_lambda_" + method_name + ".txt" );
            std::string file_name_w( "ff_testFE_w_" + method_name + ".txt" );
            std::ofstream fs_lambda;
            std::ofstream fs_w;
            fs_lambda.open( file_name_lambda.c_str() );
            fs_w.open( file_name_w.c_str() );
            Vec3r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY, m_Params.m_HalfSizeZ );
            for( Real w = -1; w <= 1.0f; w += 0.01f )
            {
                Real delta_h_alpha( m_Params.m_DegenerateThresholdDetF * 2 * m_Element2D.m_Area / m_Params.m_Dist12 );
                Real s( w / delta_h_alpha );
                Real t( Real(1) - s );
                Real lambda_p = (w < delta_h_alpha)
                                ? delta_h_alpha - w
                                : 0;
                Real lambda_r = (w < 0)
                                ? -2*w
                                : 0;
                Real lambda_u = ( w < delta_h_alpha )
                                ? ( ( w >= 0 ) ? m_Params.m_Unified_L_Factor*delta_h_alpha*( (2-t)*t*t ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                                    : m_Params.m_Unified_L_Factor*delta_h_alpha*( t + m_Params.m_Unified_NL_Factor*mal::Pow( -s, m_Params.m_Unified_NL_Exponent ) ) ) //lambda_u^L + lambda_u^NL
                                : 0;
                Real lambda(0);
                switch( m_Params.m_Method )
                {
                case Params::eMethod_PD_Project: lambda = lambda_p; break;
                case Params::eMethod_PD_Reflect: lambda = lambda_r; break;
                case Params::eMethod_PD_Unified: lambda = lambda_u; break;
                default: break;
                }
                fs_lambda << w << " " << lambda << std::endl;
                fs_w << w << " " << w + lambda << std::endl;
            }
            fs_lambda.close();
            fs_w.close();
            system( (std::string("gnuplot -e \"plot '") + file_name_lambda + "' with lines; \" -p").c_str() );
            system( (std::string("gnuplot -e \"plot '") + file_name_w + "' with lines; \" -p").c_str() );
        }

    void OnButtonCall_PlotErrorForceJacobian( void *p_user_data ) const
        {
            std::string file_name_WRP_SVD_E( "ff_testFE_error_WRP_SVD_E_2D.txt" );

            std::string file_name_LRM_SVD_T( "ff_testFE_error_LRM_SVD_T_2D.txt" );
            std::string file_name_WRP_SVD_T( "ff_testFE_error_WRP_SVD_T_2D.txt" );
            std::string file_name_LRM_SVD_H( "ff_testFE_error_LRM_SVD_H_2D.txt" );

            std::string file_name_LRM_PDU_T( "ff_testFE_error_LRM_PDU_T_2D.txt" );
            std::string file_name_WRP_PDU_T( "ff_testFE_error_WRP_PDU_T_2D.txt" );

            std::string file_name_CLRM_SVD_T( "ff_testFE_error_CLRM_SVD_T_2D.txt" );
            std::string file_name_CLRM_PDU_T( "ff_testFE_error_CLRM_PDU_T_2D.txt" );

            std::ofstream fs_WRP_SVD_E; fs_WRP_SVD_E.open( file_name_WRP_SVD_E.c_str() );
            std::ofstream fs_LRM_SVD_T; fs_LRM_SVD_T.open( file_name_LRM_SVD_T.c_str() );
            std::ofstream fs_WRP_SVD_T; fs_WRP_SVD_T.open( file_name_WRP_SVD_T.c_str() );
            std::ofstream fs_LRM_SVD_H; fs_LRM_SVD_H.open( file_name_LRM_SVD_H.c_str() );

            std::ofstream fs_LRM_PDU_T; fs_LRM_PDU_T.open( file_name_LRM_PDU_T.c_str() );
            std::ofstream fs_WRP_PDU_T; fs_WRP_PDU_T.open( file_name_WRP_PDU_T.c_str() );

            std::ofstream fs_CLRM_SVD_T; fs_CLRM_SVD_T.open( file_name_CLRM_SVD_T.c_str() );
            std::ofstream fs_CLRM_PDU_T; fs_CLRM_PDU_T.open( file_name_CLRM_PDU_T.c_str() );

            int nid( m_Params.m_NID );
            Vec2r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY );
            float max_norm_f( half_sizes.Norm() * m_Params.m_YoungModulus );
            int hcx( m_Params.m_HalfCellsX );
            int hcy( m_Params.m_HalfCellsY );
            Vec2r X[3] = { m_Element2D.m_X[0], m_Element2D.m_X[1], m_Element2D.m_X[2] };
            float length( 0.8f * half_sizes[0] / hcx );
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );

                    Mat6x6r Pf_Px_LRM_SVD_E( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_WRP_SVD_E( Mat6x6r::Identity() );

                    Mat6x6r Pf_Px_LRM_SVD_T( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_WRP_SVD_T( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_LRM_SVD_H( Mat6x6r::Identity() );

                    Mat6x6r Pf_Px_LRM_PDU_E( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_LRM_PDU_T( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_WRP_PDU_T( Mat6x6r::Identity() );

                    Mat6x6r Pf_Px_CLRM_SVD_E( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_CLRM_SVD_T( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_CLRM_PDU_E( Mat6x6r::Identity() );
                    Mat6x6r Pf_Px_CLRM_PDU_T( Mat6x6r::Identity() );

                    for( int i=0; i<6; i++ )
                    {
                        Vec6r dx( Vec6r::Zero() );
                        dx[i] = 1;
                        Vec2r dX[3] = { mal::GRange<0,1>(dx),
                                        mal::GRange<2,3>(dx),
                                        mal::GRange<4,5>(dx) };
                        // Compute common F
                        Mat2x2r F( m_Element2D.Compute_F( X, Params::eFFM_None ) );
                        // Computes R according to SVD/PDU
                        Mat2x2r R_SVD, R_PDU;
                        m_Element2D.ComputeElasticForce( X, Params::eMethod_PD_SVD, &R_SVD );
                        m_Element2D.ComputeElasticForce( X, Params::eMethod_PD_Unified, &R_PDU );
                        // Compute dR according to SVD/PDU
                        Mat2x2r dR_SVD( Mat2x2r::Zero() ), dR_PDU( Mat2x2r::Zero() );
                        bool bCorrect_dR_SVD = Compute_dR( m_Element2D,
                                                           dX[0], dX[1], dX[2],
                                                           m_Element2D.m_invDm, F, R_SVD, mal::Transposed(R_SVD),
                                                           dR_SVD );
                        bool bCorrect_dR_PDU = Compute_dR_PD_U( m_Element2D, m_Params,
                                                                X[0], X[1], X[2],
                                                                dX[0], dX[1], dX[2],
                                                                m_Element2D.m_invDm, F, R_PDU, mal::Transposed(R_PDU),
                                                                dR_PDU );
                        // Compute df_dx according to SVD/PDU LRM/WRP Approx/Exact
                        Vec6r df;
                        //--- LRM
                        // Exact df
                        df = m_Element2D.Compute_H_LRM_df_x( X, dX, R_SVD, dR_SVD );
                        mal::SetColumn( Pf_Px_LRM_SVD_E, i, df );
                        df = m_Element2D.Compute_C_df_x( X, dX, R_SVD, dR_SVD );
                        mal::SetColumn( Pf_Px_WRP_SVD_E, i, df );
                        df = m_Element2D.Compute_H_LRM_df_x( X, dX, R_PDU, dR_PDU );
                        mal::SetColumn( Pf_Px_LRM_PDU_E, i, df );
                        // Truncated df
                        Real tmp( m_Element2D.m_rParams.m_ECIE_k_factor ); //SAVE
                        // dP1
                        const_cast<Params&>(m_Element2D.m_rParams).m_ECIE_k_factor = 0;
                        df = m_Element2D.Compute_H_LRM_df_x( X, dX, R_SVD, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_LRM_SVD_T, i, df );
                        // dP2
                        df = m_Element2D.Compute_C_df_x( X, dX, R_SVD, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_WRP_SVD_T, i, df );
                        // dP Hybrid
                        const_cast<Params&>(m_Element2D.m_rParams).m_ECIE_k_factor = 1;
                        df = m_Element2D.Compute_H_LRM_df_x( X, dX, R_SVD, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_LRM_SVD_H, i, df );
                        const_cast<Params&>(m_Element2D.m_rParams).m_ECIE_k_factor = tmp; //RESTORE
                        // Truncated df PDU
                        df = m_Element2D.Compute_H_LRM_df_x( X, dX, R_PDU, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_LRM_PDU_T, i, df );
                        df = m_Element2D.Compute_C_df_x( X, dX, R_PDU, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_WRP_PDU_T, i, df );
                        //--- CLRM
                        // Exact df
                        df = m_Element2D.Compute_H_CLRM_df_x( X, dX, R_SVD, dR_SVD );
                        mal::SetColumn( Pf_Px_CLRM_SVD_E, i, df );
                        df = m_Element2D.Compute_H_CLRM_df_x( X, dX, R_PDU, dR_PDU );
                        mal::SetColumn( Pf_Px_CLRM_PDU_E, i, df );
                        // Truncated df
                        df = m_Element2D.Compute_H_CLRM_df_x( X, dX, R_SVD, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_CLRM_SVD_T, i, df );
                        df = m_Element2D.Compute_H_CLRM_df_x( X, dX, R_PDU, Mat2x2r::Zero() );
                        mal::SetColumn( Pf_Px_CLRM_PDU_T, i, df );
                    }
                    Real error_WRP_SVD_E = mal::NormF( Pf_Px_LRM_SVD_E - Pf_Px_WRP_SVD_E ) / mal::Max( mal::NormF(Pf_Px_LRM_SVD_E), mal::NormF(Pf_Px_WRP_SVD_E) );
                    Real error_LRM_SVD_T = mal::NormF( Pf_Px_LRM_SVD_E - Pf_Px_LRM_SVD_T ) / mal::NormF(Pf_Px_LRM_SVD_E);
                    Real error_WRP_SVD_T = mal::NormF( Pf_Px_LRM_SVD_E - Pf_Px_WRP_SVD_T ) / mal::NormF(Pf_Px_LRM_SVD_E);
                    Real error_LRM_SVD_H = mal::NormF( Pf_Px_LRM_SVD_E - Pf_Px_LRM_SVD_H ) / mal::NormF(Pf_Px_LRM_SVD_E);
                    Real error_LRM_PDU_T = mal::NormF( Pf_Px_LRM_PDU_E - Pf_Px_LRM_PDU_T ) / mal::NormF(Pf_Px_LRM_PDU_E);
                    Real error_WRP_PDU_T = mal::NormF( Pf_Px_LRM_PDU_E - Pf_Px_WRP_PDU_T ) / mal::NormF(Pf_Px_LRM_PDU_E);
                    Real error_CLRM_SVD_T = mal::NormF( Pf_Px_CLRM_SVD_E - Pf_Px_CLRM_SVD_T ) / mal::NormF(Pf_Px_CLRM_SVD_E);
                    Real error_CLRM_PDU_T = mal::NormF( Pf_Px_CLRM_PDU_E - Pf_Px_CLRM_PDU_T ) / mal::NormF(Pf_Px_CLRM_PDU_E);
                    // Real error_CLRM_PDU_T = mal::SignedNormInf( Pf_Px_LRM_SVD_E - Pf_Px_LRM_SVD_T ) / mal::Max( mal::NormInf( Pf_Px_LRM_SVD_E ),  mal::NormInf( Pf_Px_LRM_SVD_T ) ); //\todo TEMP: overestimation measure
                    fs_WRP_SVD_E << X[nid].x() << " " << X[nid].y() << " " << error_WRP_SVD_E << std::endl;

                    fs_LRM_SVD_T << X[nid].x() << " " << X[nid].y() << " " << error_LRM_SVD_T << std::endl;
                    fs_WRP_SVD_T << X[nid].x() << " " << X[nid].y() << " " << error_WRP_SVD_T << std::endl;
                    fs_LRM_SVD_H << X[nid].x() << " " << X[nid].y() << " " << error_LRM_SVD_H << std::endl;

                    fs_LRM_PDU_T << X[nid].x() << " " << X[nid].y() << " " << error_LRM_PDU_T << std::endl;
                    fs_WRP_PDU_T << X[nid].x() << " " << X[nid].y() << " " << error_WRP_PDU_T << std::endl;

                    fs_CLRM_SVD_T << X[nid].x() << " " << X[nid].y() << " " << error_CLRM_SVD_T << std::endl;
                    fs_CLRM_PDU_T << X[nid].x() << " " << X[nid].y() << " " << error_CLRM_PDU_T << std::endl;
                }
                fs_WRP_SVD_E << std::endl;
                fs_LRM_SVD_T << std::endl;
                fs_WRP_SVD_T << std::endl;
                fs_LRM_SVD_H << std::endl;
                fs_LRM_PDU_T << std::endl;
                fs_WRP_PDU_T << std::endl;
                fs_CLRM_SVD_T << std::endl;
                fs_CLRM_PDU_T << std::endl;
            }
            fs_WRP_SVD_E.close();
            fs_LRM_SVD_T.close();
            fs_WRP_SVD_T.close();
            fs_LRM_SVD_H.close();
            fs_LRM_PDU_T.close();
            fs_WRP_PDU_T.close();
            fs_CLRM_SVD_T.close();
            fs_CLRM_PDU_T.close();

            // FUCKING BLOODY SHIT nans stop gnuplot form working...
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_WRP_SVD_E).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_LRM_SVD_T).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_WRP_SVD_T).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_LRM_SVD_H).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_LRM_PDU_T).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_WRP_PDU_T).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_CLRM_SVD_T).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_CLRM_PDU_T).c_str() );

            //system( (std::string("gnuplot -e \"splot '") + file_name_WRP_SVD_E + "' with pm3d  \" -p").c_str() );
            system( (std::string("gnuplot -e \"splot '") + file_name_LRM_SVD_T + "' with pm3d title 'dP(1)' \" -p").c_str() );
            system( (std::string("gnuplot -e \"splot '") + file_name_WRP_SVD_T + "' with pm3d title 'dP(2)' \" -p").c_str() );
            system( (std::string("gnuplot -e \"splot '") + file_name_LRM_SVD_H + "' with pm3d title 'dP(H)' \" -p").c_str() );
            //system( (std::string("gnuplot -e \"splot '") + file_name_LRM_PDU_T + "' with pm3d  \" -p").c_str() );
            //system( (std::string("gnuplot -e \"splot '") + file_name_WRP_PDU_T + "' with pm3d  \" -p").c_str() );
            //system( (std::string("gnuplot -e \"splot '") + file_name_CLRM_SVD_T + "' with pm3d \" -p").c_str() );
            //system( (std::string("gnuplot -e \"splot '") + file_name_CLRM_PDU_T + "' with pm3d \" -p").c_str() );
        }

    void OnButtonCall_PlotErrorForceDAPD( void *p_user_data ) const
        {
            std::string file_name_abs( "ff_testFE_error_DAPD_SS_vs_EX_2D_Abs.txt" );
            std::string file_name_rel( "ff_testFE_error_DAPD_SS_vs_EX_2D_Rel.txt" );
            std::string file_name_mix( "ff_testFE_error_DAPD_SS_vs_EX_2D_Mix.txt" );
            std::string file_name_angle( "ff_testFE_error_DAPD_SS_vs_EX_2D_Angle.txt" );
            std::string file_name_norm( "ff_testFE_error_DAPD_SS_vs_EX_2D_Norm.txt" );
            std::ofstream fs_abs; fs_abs.open( file_name_abs.c_str() );
            std::ofstream fs_rel; fs_rel.open( file_name_rel.c_str() );
            std::ofstream fs_mix; fs_mix.open( file_name_mix.c_str() );
            std::ofstream fs_angle; fs_angle.open( file_name_angle.c_str() );
            std::ofstream fs_norm; fs_norm.open( file_name_norm.c_str() );
            int nid( m_Params.m_NID );
            Vec2r half_sizes( m_Params.m_HalfSizeX, m_Params.m_HalfSizeY );
            int hcx( m_Params.m_HalfCellsX );
            int hcy( m_Params.m_HalfCellsY );
            Vec2r X[3] = { m_Element2D.m_X[0], m_Element2D.m_X[1], m_Element2D.m_X[2] };
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Vec6r f_ss = m_Element2D.ComputeElasticForce( X, Params::eMethod_DAPD_SS );
                    Vec6r f_ex = m_Element2D.ComputeElasticForce( X, Params::eMethod_DAPD_EX );
                    // Real error = mal::Norm( f_ss - f_ex ) / mal::Max( mal::Norm(f_ss), mal::Norm(f_ex) );

                    // f_ss.Normalize();
                    // f_ex.Normalize();
                    // Real error( mal::Max<Real>(0, 1 - mal::Dot(f_ss,f_ex) ) ); //\todo Value -0.05 appears, meaning >90º, but I doubt it, can be a numerical derivation error
                    // Real error = mal::Rad2Deg(mal::ACos(mal::Max<Real>(0,mal::Dot(f_ss,f_ex))));

                    Vec2r f_ss1( mal::GRange<0,1>(f_ss) );
                    Vec2r f_ss2( mal::GRange<2,3>(f_ss) );
                    Vec2r f_ss3( mal::GRange<4,5>(f_ss) );
                    Vec2r f_ex1( mal::GRange<0,1>(f_ex) );
                    Vec2r f_ex2( mal::GRange<2,3>(f_ex) );
                    Vec2r f_ex3( mal::GRange<4,5>(f_ex) );

                    // Real error_std_max = mal::Max( mal::Max( (f_ss1-f_ex1).Norm()/f_ex1.Norm(),
                    //                                          (f_ss2-f_ex2).Norm()/f_ex2.Norm() ),
                    //                                (f_ss3-f_ex3).Norm()/f_ex3.Norm() );
                    // Real error_std_avg = Real(0.3333333f)*( (f_ss1-f_ex1).Norm()/f_ex1.Norm()
                    //                                         + (f_ss2-f_ex2).Norm()/f_ex2.Norm()
                    //                                         + (f_ss3-f_ex3).Norm()/f_ex3.Norm() );
                    // Real error_std_abs = mal::Max( mal::Max( (f_ss1-f_ex1).Norm(), (f_ss2-f_ex2).Norm() ),
                    //                                (f_ss3-f_ex3).Norm() );

                    // Real error_std_mix = mal::Max( mal::Max( (f_ss1-f_ex1).Norm() * (f_ss1-f_ex1).Norm()/f_ex1.Norm(),
                    //                                          (f_ss2-f_ex2).Norm() * (f_ss2-f_ex2).Norm()/f_ex2.Norm() ),
                    //                                (f_ss3-f_ex3).Norm() * (f_ss3-f_ex3).Norm()/f_ex3.Norm() );


                    Real error_std_rel = (f_ss-f_ex).Norm()/f_ex.Norm();
                    Real error_std_abs = (f_ss-f_ex).Norm();
                    Real error_std_mix = error_std_abs*error_std_rel;
                    // Real error_std_rel = (f_ss2-f_ex2).Norm()/f_ex2.Norm();
                    // Real error_std_abs = (f_ss2-f_ex2).Norm();
                    // Real error_std_mix = error_std_abs*error_std_rel;
                    Real error_norm_max = mal::Max( mal::Max( mal::Abs(f_ss1.Norm()-f_ex1.Norm())/f_ex1.Norm(),
                                                              mal::Abs(f_ss2.Norm()-f_ex2.Norm())/f_ex2.Norm() ),
                                                    mal::Abs( f_ss3.Norm()-f_ex3.Norm())/f_ex3.Norm() );
                    Real error_norm_avg = Real(0.3333333f)*( mal::Abs(f_ss1.Norm()-f_ex1.Norm())/f_ex1.Norm()
                                                             + mal::Abs(f_ss2.Norm()-f_ex2.Norm())/f_ex2.Norm()
                                                             + mal::Abs( f_ss3.Norm()-f_ex3.Norm())/f_ex3.Norm() );

                    f_ss1.Normalize();
                    f_ss2.Normalize();
                    f_ss3.Normalize();
                    f_ex1.Normalize();
                    f_ex2.Normalize();
                    f_ex3.Normalize();
                    Real error_angle_max = mal::Rad2Deg( mal::Max( mal::Max( mal::ACos(mal::Max<Real>(-1,mal::Dot(f_ss1,f_ex1))),
                                                                             mal::ACos(mal::Max<Real>(-1,mal::Dot(f_ss2,f_ex2)))),
                                                                   mal::ACos(mal::Max<Real>(-1,mal::Dot(f_ss3,f_ex3)))));
                    Real error_angle_avg = Real(0.3333333f) * (mal::Rad2Deg( mal::ACos(mal::Max<Real>(-1,mal::Dot(f_ss1,f_ex1))))
                                                               + mal::Rad2Deg( mal::ACos(mal::Max<Real>(-1,mal::Dot(f_ss1,f_ex1))))
                                                               + mal::Rad2Deg( mal::ACos(mal::Max<Real>(-1,mal::Dot(f_ss1,f_ex1)))) );

                    fs_abs << X[nid].x() << " " << X[nid].y() << " " << error_std_abs << std::endl;
                    fs_rel << X[nid].x() << " " << X[nid].y() << " " << error_std_rel << std::endl;
                    fs_mix << X[nid].x() << " " << X[nid].y() << " " << error_std_mix << std::endl;
                    fs_angle << X[nid].x() << " " << X[nid].y() << " " << error_angle_max << std::endl;
                    fs_norm << X[nid].x() << " " << X[nid].y() << " " << error_norm_max << std::endl;
                }
                fs_abs << std::endl;
                fs_rel << std::endl;
                fs_mix << std::endl;
                fs_angle << std::endl;
                fs_norm << std::endl;
            }
            fs_abs.close();
            fs_rel.close();
            fs_mix.close();
            fs_angle.close();
            fs_norm.close();
            // FUCKING BLOODY SHIT nans stop gnuplot form working...
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_abs).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_rel).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_mix).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_angle).c_str() );
            system( (std::string("sed -i.bak s/nan/0/g ") + file_name_norm).c_str() );
            system( (std::string("gnuplot -e \"splot '") + file_name_abs + "' with pm3d title 'abs'\" -p").c_str() );
            system( (std::string("gnuplot -e \"splot '") + file_name_rel + "' with pm3d title 'rel'\" -p").c_str() );
            system( (std::string("gnuplot -e \"splot '") + file_name_mix + "' with pm3d title 'mix'\" -p").c_str() );
            // system( (std::string("gnuplot -e \"splot '") + file_name_angle + "' with pm3d title 'angle'\" -p").c_str() );
            // system( (std::string("gnuplot -e \"splot '") + file_name_norm + "' with pm3d title 'norm'\" -p").c_str() );
        }

public:
    Params m_Params;
    Element2D m_Element2D;
    Element3D m_Element3D;

    bool m_bIsActive;
    sfr::IView *m_pAppView;

private:
    friend class AppInputListener;
    friend class AppRenderer;

    int m_Dimension;
    bool m_bIsEnabledDrag;
    Vec2f m_DragPoint;

    util::ArchetypeLibrary m_ArchetypeLibrary;
    util::ItemStream m_TmpIS;

    util::GStatician<float> m_TotalEnergyST;
    util::GStatician<float> m_ElasticEnergyST;
    util::GStatician<float> m_KineticEnergyST;
    util::GStatician<float> m_AngularMomentumST;
    util::GStatician<float> m_DeformationRatioST;

    util::SimpleProfiler m_Profiler;
};

#endif //TEST_FE_APPTASK_H
