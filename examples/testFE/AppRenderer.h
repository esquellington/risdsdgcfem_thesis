#ifndef TEST_FE_APPRENDERER_H
#define TEST_FE_APPRENDERER_H

#include "Config.h"
#include "AppTask.h"
#include <Safra/core/gfx/IRenderer.h>
#include "DifferentialR.h"

#define __ENABLE_COMPARE_SVD1_AND_LRM

//#define __ENABLE_COMPUTE_APPROX_JACOBIAN_ERROR
//#define __ENABLE_COMPUTE_APPROX_JACOBIAN_ERROR_EIGEN
#ifdef __ENABLE_COMPUTE_APPROX_JACOBIAN_ERROR_EIGEN
#include <Eigen/Dense>
#include <Eigen/LU>
Mat6x6r ComputeInverse( const Mat6x6r& M )
{
    typedef Eigen::Matrix<double,6,6> M66f;

    M66f M_eigen;
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            M_eigen( i, j ) = M(i,j);
    std::cout << "M_eigen = " << M_eigen << std::endl;

    Eigen::FullPivLU<M66f> lu(M_eigen);
    //M66f invM_eigen = M_eigen.inverse();
    M66f invM_eigen = lu.inverse();

    std::cout << "isInvertible = " << lu.isInvertible() << std::endl;
    std::cout << "invM_eigen = " << invM_eigen << std::endl;

    std::cout << "Id_eigen = " << M_eigen*invM_eigen << std::endl;

    Mat6x6r invM;
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            invM( i, j ) = invM_eigen(i,j);

    return invM;
}

Real ComputeSpectralRadius_Nonsymmetric( const Mat6x6r& M )
{
    return 0;
}
#elif defined(__ENABLE_COMPUTE_APPROX_JACOBIAN_ERROR_GSL)
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_permutation.h>
#  include <gsl/gsl_linalg.h>
#  include <gsl/gsl_cblas.h>
#  include <gsl/gsl_blas.h>
#  include <gsl/gsl_eigen.h>
Mat6x6r ComputeInverse( const Mat6x6r& M )
{
    // Alloc
    gsl_matrix* M_gsl = gsl_matrix_alloc(6,6);
    gsl_matrix* invM_gsl = gsl_matrix_alloc(6,6);
    gsl_permutation* perm_gsl = gsl_permutation_alloc(6);

    /* Fill the matrix M
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            gsl_matrix_set( M_gsl, i, j, M(i,j) );
    */

    /*
    double values[] = { -599.541, -13.5421, 313.31, -192.77, 286.231, 206.312,
                        -13.5421, -200.459, -192.77, 86.6883, 206.312, 113.771,
                        313.31, -192.77, -363.883, 209.694, 50.5735, -16.9241,
                        -192.77, 86.6883, 209.694, -636.115, -16.9241, 549.427,
                        286.231, 206.312, 50.5735, -16.9241, -336.804, -189.388,
                        206.312, 113.771, -16.9241, 549.427, -189.388, -663.197 };
    Mat6x6r M_tmp;
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            M_tmp(i,j) = values[i*6+j];
    std::cout << "M_tmp = " << M_tmp << std::endl;
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            gsl_matrix_set( M_gsl, i, j, values[i*6+j] );//M_tmp(i,j) );
    */

    /*WORKS
    double values[] = { -599.541, -13.5421, 313.31, -192.77, 286.231, 206.312,
                        -13.5421, -200.459, -192.77, 86.6883, 206.312, 113.771,
                        313.31, -192.77, -363.883, 209.694, 50.5735, -16.9241,
                        -192.77, 86.6883, 209.694, -636.115, -16.9241, 549.427,
                        286.231, 206.312, 50.5735, -16.9241, -336.804, -189.388,
                        206.312, 113.771, -16.9241, 549.427, -189.388, -663.197 };
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            gsl_matrix_set( M_gsl, i, j, values[i*6+j] );
    */

    //
    double values[6*6];
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            values[i*6+j] = M(i,j);
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            gsl_matrix_set( M_gsl, i, j, values[i*6+j] );

    // Make LU decomposition of matrix m
    int s;
    gsl_linalg_LU_decomp( M_gsl, perm_gsl, &s );
    // Invert the matrix m
    gsl_linalg_LU_invert( M_gsl, perm_gsl, invM_gsl );
    // Get result
    Mat6x6r invM;
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            invM(i,j) = gsl_matrix_get( invM_gsl, i, j );

    // Cleanup
    gsl_permutation_free( perm_gsl );
    gsl_matrix_free( invM_gsl );
    gsl_matrix_free( M_gsl );

    return invM;
}

Real ComputeSpectralRadius_Nonsymmetric( const Mat6x6r& M )
{
    // Fill the matrix M
    gsl_matrix* M_gsl = gsl_matrix_alloc(6,6);
    for( int i=0; i<6; i++ )
        for( int j=0; j<6; j++ )
            gsl_matrix_set( M_gsl, i, j, M(i,j) );
    // Compute eigen
    gsl_vector_complex* eval = gsl_vector_complex_alloc(6);
    gsl_matrix_complex* evec = gsl_matrix_complex_alloc(6,6);
    gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(6);
    gsl_eigen_nonsymmv( M_gsl, eval, evec, w );
    gsl_eigen_nonsymmv_free(w);
    gsl_eigen_nonsymmv_sort( eval, evec, GSL_EIGEN_SORT_ABS_DESC );
    // Get largest EV
    gsl_complex eval_max = gsl_vector_complex_get(eval,0);
    Real rho( mal::Sqrt( mal::Sq(GSL_REAL(eval_max)) + mal::Sq(GSL_IMAG(eval_max)) ) );
    /*
    for( int i = 0; i < 6; i++ )
    {
        gsl_complex eval_i = gsl_vector_complex_get(eval,i);
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,i);
        printf ("eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for( int j = 0; j < 6; ++j )
        {
            gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
        }
    }
    */
    // Cleanup
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_matrix_free(M_gsl);

    return rho;
}

void test_fucking_gsl_invert()
{
    // Define the dimension n of the matrix
    // and the signum s (for LU decomposition)
    int n = 6;
    int s;

    // Define all the used matrices
    gsl_matrix * m = gsl_matrix_alloc (n, n);
    gsl_matrix * inverse = gsl_matrix_alloc (n, n);
    gsl_matrix * identity = gsl_matrix_alloc (n, n);
    gsl_permutation * perm = gsl_permutation_alloc (n);

    /* Fill the matrix m
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            gsl_matrix_set( m, i, j, 1 + i + j + std::max(i,j)*10 + i*j ); //(i==j)?1:0 ); //float(4*i+1)/(j*j+1) );
    */
    double values[] = { -583.21, 80.2121, 211.393, -223.316, 371.817, 143.104,
                        80.2121, -216.79, -223.316, 188.607, 143.104, 28.1828,
                        211.393, -223.316, -282.38, 123.051, 70.9871, 100.265,
                        -223.316, 188.607, 123.051, -717.62, 100.265, 529.013,
                        371.817, 143.104, 70.9871, 100.265, -442.804, -243.369,
                        143.104, 28.1828, 100.265, 529.013, -243.369, -557.196 };
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            gsl_matrix_set( m, i, j, values[i*n+j] );

    std::cout << "m = " << std::endl;
    std::cout << "[ ";
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            std::cout << gsl_matrix_get( m, i, j ) << ((j<n-1) ? ", " : ";\n");
    std::cout << " ]" << std::endl;

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp (m, perm, &s);

    // Invert the matrix m
    gsl_linalg_LU_invert (m, perm, inverse);

    // identity = m*inverse
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, m, inverse,
                    0.0, identity );

    std::cout << "Inv = " << std::endl;
    std::cout << "[ ";
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            std::cout << gsl_matrix_get( inverse, i, j ) << ((j<n-1) ? ", " : ";\n");
    std::cout << " ]" << std::endl;

    std::cout << "Id = " << std::endl;
    std::cout << "[ ";
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            std::cout << gsl_matrix_get( identity, i, j ) << ((j<n-1) ? ", " : ";\n");
    std::cout << " ]" << std::endl;
}
#endif

#define __ENABLE_DRAW_DCLFEM

/*! App-specific Element2D Renderer with volatile data (raycast hits, intersections, etc...) */
class AppRenderer: public sfr::gfx::RendererPropagated
{
public:
    enum EDrawFlags {
        eDraw_Nothing = 0,
        eDraw_All     = 0xFFFFFFFF,
        eDraw_Default = eDraw_All
    };

public:
    AppRenderer( AppTask &app_task, sfr::gfx::IRenderer *p_renderer, int32 draw_flags )
    : RendererPropagated( p_renderer )
    , m_rAppTask( app_task )
    , m_DrawFlags( draw_flags )
        {}

    bool Render()
    {
        switch( m_rAppTask.GetElementDimension() )
        {
        case 2: return Render2D(); break;
        case 3: return Render3D(); break;
        default: return false; break;
        }
    }

    bool Render2D()
    {
        const Params &params( m_rAppTask.GetParams() );
        int nid( params.m_NID );

        if( params.m_DrawFlags.Test(Params::eDraw_Axis) )
            DrawRefSys( Vec3f::Zero(), Mat3x3f::Identity(), 25.0f, sfr::gfx::Style() );

        float grid_size = 50.0f;
        float half_grid_size = 0.5f * grid_size;
        if( params.m_DrawFlags.Test(Params::eDraw_Grid) )
            DrawGrid( Vec3f(-half_grid_size,-half_grid_size,0),
                      Vec3f(1,0,0), Vec3f(0,1,0),
                      grid_size, grid_size,
                      10, 10,
                      sfr::gfx::Style(sfr::gfx::Color(0.1,0.1,0.5,1),1) );

        // Element
        Element2D &element( m_rAppTask.GetElement2D() );
        if( params.m_DrawFlags.Test(Params::eDraw_Element) )
        {
#ifdef __ENABLE_DRAW_DCLFEM
            // Rest
            DrawSegment2( element.m_r[0], element.m_r[1], sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),2) );
            DrawSegment2( element.m_r[1], element.m_r[2], sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),2) );
            DrawSegment2( element.m_r[2], element.m_r[0], sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),2) );
#else
            // Rest
            DrawPoint2( element.m_r[0], sfr::gfx::Style(sfr::gfx::Color(1,0.25,0.25,1),10) );
            DrawPoint2( element.m_r[1], sfr::gfx::Style(sfr::gfx::Color(0.25,1,0.25,1),10) );
            DrawPoint2( element.m_r[2], sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,1,1),10) );
            DrawSegment2( element.m_r[0], element.m_r[1], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment2( element.m_r[1], element.m_r[2], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment2( element.m_r[2], element.m_r[0], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            // Rotated
            DrawPoint2( element.m_X[0], sfr::gfx::Style(sfr::gfx::Color(0.5,0.25,0.25,1),6) );
            DrawPoint2( element.m_X[1], sfr::gfx::Style(sfr::gfx::Color(0.25,0.5,0.25,1),6) );
            DrawPoint2( element.m_X[2], sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,0.5,1),6) );
            DrawSegment2( element.m_X[0], element.m_X[1], sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,1),2) );
            DrawSegment2( element.m_X[1], element.m_X[2], sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,1),2) );
            DrawSegment2( element.m_X[2], element.m_X[0], sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,1),2) );
#endif
            // Show eps det(F) threshold
            Real delta_h_alpha( params.m_DegenerateThresholdDetF * 2 * element.m_Area / params.m_Dist12 );
            DrawSegment2( Vec2r( delta_h_alpha, -params.m_HalfSizeY ),
                          Vec2r( delta_h_alpha,  params.m_HalfSizeY ),
                          sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0.5,1),1) );
        }

        // Draw Drag stuff
        if( true )//m_rAppTask.m_bIsEnabledDrag )
        {
            DrawPoint2( m_rAppTask.m_DragPoint, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),10) );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            X[nid] = Vec2r( m_rAppTask.m_DragPoint );

            /*TEMP: Draw SVD details
            // Compute F
            Mat2x2r Ds( mal::GMat2x2_From_Columns( X[1] - X[0], X[2] - X[0] ) );
            Mat2x2r F( Q * mal::Inverse( mal::GMat2x2_From_Columns( element.m_r[1] - element.m_r[0], element.m_r[2] - element.m_r[0] ) ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            // Draw eigenvalues
            char str[256];
            sprintf( str, "%f + %f = %f", vec_diag_F[0], vec_diag_F[1], vec_diag_F[0] + vec_diag_F[1] );
            DrawLabel2( X[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,1,1),1) );
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
            }
            Mat2x2r R( U*Vt );
            DrawRefSys2( X[nid], R, 1, sfr::gfx::Style() );
            DrawRefSys2( X[nid], U, 0.66, sfr::gfx::Style() );
            DrawRefSys2( X[nid], Vt.Transposed(), 0.33, sfr::gfx::Style() );
            */

            // Drag force/rotation
            if( params.m_DrawFlags.Test( Params::eDraw_Drag_f0
                                         | Params::eDraw_Drag_R
                                         | Params::eDraw_Drag_EE0 ) )
            {
                Mat2x2r R;
                Real elastic_energy(0);
                Vec6r f( element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R, &elastic_energy ) );
                if( params.m_DrawFlags.Test(Params::eDraw_Drag_f0) )
                {
                    // Draw elastic force
                    Vec2r f0( f[0], f[1] );
                    Vec2r f1( f[2], f[3] );
                    Vec2r f2( f[4], f[5] );
                    float max_norm_f( mal::Max( mal::Max( f0.Norm(), f1.Norm() ), f2.Norm() ) );
                    if( max_norm_f > mal::Epsilon<float>() );
                    {
                        Vec2r v0( f0 / max_norm_f );
                        Vec2r v1( f1 / max_norm_f );
                        Vec2r v2( f2 / max_norm_f );
                        DrawVector2( sfr::Vec2(X[0]), sfr::Vec2(v0), sfr::gfx::Style( sfr::gfx::Color( 1.0f, 0.25f, 0.25f,1 ), 2 ) );
                        DrawVector2( sfr::Vec2(X[1]), sfr::Vec2(v1), sfr::gfx::Style( sfr::gfx::Color( 0.25f, 1.0f, 0.25f,1 ), 2 ) );
                        DrawVector2( sfr::Vec2(X[2]), sfr::Vec2(v2), sfr::gfx::Style( sfr::gfx::Color( 0.25f, 0.25f, 1.0f,1 ), 2 ) );
                    }
                }
                if( params.m_DrawFlags.Test(Params::eDraw_Drag_EE0) )
                {
                    char str[128];
                    sprintf( str, "EE = %f", elastic_energy );
                    DrawLabel2( X[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),1) );
                }
                if( params.m_DrawFlags.Test(Params::eDraw_Drag_R) )
                {
                    // Draw inferred rotation
                    DrawRefSys2( X[nid], R, 1, sfr::gfx::Style() );
                }
            }

            // Drag trajectory
            if( params.m_DrawFlags.Test(Params::eDraw_Drag_t0) )
            {
                // Init statistics
                Real elastic_energy(0);
                Real *p_elastic_energy(0);
                unsigned int num_samples( mal::Ceil( params.m_Duration / params.m_TimeStep ) );
                if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) ) { m_rAppTask.m_TotalEnergyST.Begin( num_samples ); p_elastic_energy = &elastic_energy; }
                if( params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) ) { m_rAppTask.m_ElasticEnergyST.Begin( num_samples ); p_elastic_energy = &elastic_energy; }
                if( params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) ) m_rAppTask.m_KineticEnergyST.Begin( num_samples );
                if( params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) ) m_rAppTask.m_AngularMomentumST.Begin( num_samples );
                if( params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) ) m_rAppTask.m_DeformationRatioST.Begin( num_samples );

                // Init trajectory
                Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
                Mat2x2r R;
                Vec2r prevX[3] = { X[0], X[1], X[2] };
                Vec2r velX[3] = { Vec2r::Zero(), Vec2r::Zero(), Vec2r::Zero() };
                Real t(0);
                Real inv_node_mass( Real(3)/params.m_Mass );
                Real damping_coeff( params.m_DampingRatio * 2 * mal::Sqrt( params.m_YoungModulus * params.m_Mass ) );
                // Init degeneration tracking and NoC
                Real det_F( mal::Det( element.Compute_F( X, Params::eFFM_None ) ) );
                bool bDegenerate( det_F <= params.m_DegenerateThresholdDetF );
                int trajectory_noc( bDegenerate ? params.m_NoC : -1 ); //Assume NoC = 0 if starts degenerate
#ifdef __ENABLE_DRAW_DCLFEM
                unsigned int cNumGhostFrames(25);
                float acc_screenshot_time(0), screenshot_timestep( params.m_Duration / cNumGhostFrames );
#endif
                do
                {
                    Vec6r f( element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R, p_elastic_energy ) );
                    // Draw elastic force
                    Vec2r f0( f[0], f[1] );
                    Vec2r f1( f[2], f[3] );
                    Vec2r f2( f[4], f[5] );
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x0) )
                    {
                        velX[0] += params.m_TimeStep * inv_node_mass * (f0 - damping_coeff * velX[0]);
                        X[0] += params.m_TimeStep * velX[0];
#ifdef __ENABLE_DRAW_DCLFEM
                        DrawSegment2( sfr::Vec2(prevX[0]), sfr::Vec2(X[0]), sfr::gfx::Style( sfr::gfx::Color( 1.0f, 0.25f, 0.25f,1 ), 4 ) );
#else
                        DrawSegment2( sfr::Vec2(prevX[0]), sfr::Vec2(X[0]), sfr::gfx::Style( sfr::gfx::Color( 1.0f, 0.25f, 0.25f,1 ), 1 ) );
#endif
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x1) )
                    {
                        velX[1] += params.m_TimeStep * inv_node_mass * (f1 - damping_coeff * velX[1]);
                        X[1] += params.m_TimeStep * velX[1];
#ifdef __ENABLE_DRAW_DCLFEM
                        DrawSegment2( sfr::Vec2(prevX[1]), sfr::Vec2(X[1]), sfr::gfx::Style( sfr::gfx::Color( 0.25f, 1.0f, 0.25f,1 ), 4 ) );
#else
                        DrawSegment2( sfr::Vec2(prevX[1]), sfr::Vec2(X[1]), sfr::gfx::Style( sfr::gfx::Color( 0.25f, 1.0f, 0.25f,1 ), 1 ) );
#endif
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x2) )
                    {
                        velX[2] += params.m_TimeStep * inv_node_mass * (f2 - damping_coeff * velX[2]);
                        X[2] += params.m_TimeStep * velX[2];
#ifdef __ENABLE_DRAW_DCLFEM
                        DrawSegment2( sfr::Vec2(prevX[2]), sfr::Vec2(X[2]), sfr::gfx::Style( sfr::gfx::Color( 0.25f, 0.25f, 1.0f,1 ), 4 ) );
#else
                        DrawSegment2( sfr::Vec2(prevX[2]), sfr::Vec2(X[2]), sfr::gfx::Style( sfr::gfx::Color( 0.25f, 0.25f, 1.0f,1 ), 1 ) );
#endif
                    }
                    t += params.m_TimeStep;


#ifdef __ENABLE_DRAW_DCLFEM
                    acc_screenshot_time += params.m_TimeStep;
                    if( acc_screenshot_time > screenshot_timestep )
                    {
                        Real lambda01( t / params.m_Duration );
                        sfr::gfx::Color col( sfr::gfx::Color(lambda01,1-lambda01,0,1) );
                        /*
                        DrawSegment2( X[0], X[1], sfr::gfx::Style(col,4) );
                        DrawSegment2( X[1], X[2], sfr::gfx::Style(col,4) );
                        DrawSegment2( X[2], X[0], sfr::gfx::Style(col,4) );
                        */
                        DrawRefSys2( X[nid], R, 0.25, sfr::gfx::Style(col,2) );
                        acc_screenshot_time -= screenshot_timestep;
                    }
#endif

                    // Compute ToC->NoC
                    det_F = mal::Det( element.Compute_F( X, Params::eFFM_None ) );
                    if( !bDegenerate && det_F <= params.m_DegenerateThresholdDetF )
                    {
                        bDegenerate = true;
                        if( params.m_bEnableNoC ) trajectory_noc = ComputeNoC( prevX, X, params.m_DegenerateThresholdDetF, element.m_Area );
                        else trajectory_noc = params.m_NoC;
                    }
                    else if( det_F > params.m_DegenerateThresholdDetF )
                    {
                        bDegenerate = false;
                        trajectory_noc = -1;
                    }
                    element.m_NoC = trajectory_noc;

                    // Advance
                    prevX[0] = X[0];
                    prevX[1] = X[1];
                    prevX[2] = X[2];

                    // Gather statistics
                    if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy | Params::eTF_ElasticEnergy | Params::eTF_KineticEnergy ) )
                    {
                        //Real kinetic_energy( Real(0.5) * params.m_Mass * ( mal::NormSq(velX[0]) + mal::NormSq(velX[1]) + mal::NormSq(velX[2]) ) ); //Lumped KE
                        //Consistent KE seems more accurate than lumped M-KE... both behave similarly but the scale of consistent KE is closer to EE, which is expected...
                        Real kinetic_energy( ( params.m_Mass / Real(12) )
                                             * ( mal::NormSq(velX[0]) + mal::NormSq(velX[1]) + mal::NormSq(velX[2])
                                                 + velX[0].x()*velX[1].x() + velX[1].x()*velX[2].x() + velX[0].x()*velX[2].x()
                                                 + velX[0].y()*velX[1].y() + velX[1].y()*velX[2].y() + velX[0].y()*velX[2].y() ) );
                        //
                        if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) )
                            m_rAppTask.m_TotalEnergyST.Add( kinetic_energy + elastic_energy );
                        if( params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) )
                            m_rAppTask.m_ElasticEnergyST.Add( elastic_energy );
                        if( params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) )
                            m_rAppTask.m_KineticEnergyST.Add( kinetic_energy );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) )
                    {
                        Vec2r CoM( mal::Rcp<float>(3) * ( X[0] + X[2] + X[2] ) );
                        m_rAppTask.m_AngularMomentumST.Add( mal::Dot( velX[0], mal::PerpendicularCW( X[0] - CoM ) )
                                                            + mal::Dot( velX[1], mal::PerpendicularCW( X[1] - CoM ) )
                                                            + mal::Dot( velX[2], mal::PerpendicularCW( X[2] - CoM ) ) );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) )
                        m_rAppTask.m_DeformationRatioST.Add( det_F );
                } while( t < params.m_Duration );

                //Restore Element2D.NoC to Params.NoC, changed through GUI
                element.m_NoC = params.m_NoC;

                // End statistics and show results
                if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) ) m_rAppTask.m_TotalEnergyST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) ) m_rAppTask.m_ElasticEnergyST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) ) m_rAppTask.m_KineticEnergyST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) ) m_rAppTask.m_AngularMomentumST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) ) m_rAppTask.m_DeformationRatioST.End();
                // Final configuration
                {
#ifdef __ENABLE_DRAW_DCLFEM
                    sfr::gfx::Color col( sfr::gfx::Color(1,0,0,1) );
                    DrawSegment2( X[0], X[1], sfr::gfx::Style(col,4) );
                    DrawSegment2( X[1], X[2], sfr::gfx::Style(col,4) );
                    DrawSegment2( X[2], X[0], sfr::gfx::Style(col,4) );
#else
                    sfr::gfx::Color col( det_F > 0
                                         ? sfr::gfx::Color(0.25,0.75,0.25,1)
                                         : sfr::gfx::Color(0.75,0.25,0.25,1) );
                    DrawSegment2( X[0], X[1], sfr::gfx::Style(col,2) );
                    DrawSegment2( X[1], X[2], sfr::gfx::Style(col,2) );
                    DrawSegment2( X[2], X[0], sfr::gfx::Style(col,2) );
#endif
                }
            }
        }

        // Sample force field
        if( params.m_DrawFlags.Test(Params::eDraw_f0_x)
            || params.m_DrawFlags.Test(Params::eDraw_f1_x)
            || params.m_DrawFlags.Test(Params::eDraw_f2_x) )
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * params.m_YoungModulus );
            int idx( params.m_DrawFlags.Test(Params::eDraw_f0_x)
                     ? 0
                     : params.m_DrawFlags.Test(Params::eDraw_f1_x) ? 1 : 2 );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Vec6r f( element.ComputeElasticForce( X, (Params::EMethod)params.m_Method ) );
                    Vec2r f_i( f[2*idx], f[2*idx+1] );
                    float c_i( f_i.Norm() / max_norm_f );
                    Vec2r v_i( length * mal::SafeNormalized( f_i ) );
                    DrawVector( sfr::Vec3(X[nid].x(),X[nid].y(),c_i), sfr::Vec3(v_i.x(),v_i.y(),0),
                                sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c_i,
                                                                  0.25f,
                                                                  0.25f + 0.75f*(1.0-c_i),
                                                                  1 ),
                                                 1 ) );
                }
            }
        }

        // Sample energy field
        if( true )//params.m_DrawFlags.Test(Params::eDraw_e_x0) ) //\todo plot to file?
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * params.m_YoungModulus );
            int idx( params.m_DrawFlags.Test(Params::eDraw_f0_x)
                     ? 0
                     : params.m_DrawFlags.Test(Params::eDraw_f1_x) ? 1 : 2 );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            Mat2x2r R;
            Real energy(0);
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R, &energy );
                    DrawSegment( sfr::Vec3(X[nid].x(),X[nid].y(),0),
                                 sfr::Vec3(X[nid].x(),X[nid].y(),energy/5000),
                                 sfr::gfx::Style( sfr::gfx::Color( energy/5000, 0, 0, 1 ),1 ) );
                }
            }
        }

        if( params.m_DrawFlags.Test(Params::eDraw_P_e) )
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * params.m_YoungModulus );
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    Vec2r vec_diag_F( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Vec2r vec_diag_P( - element.ComputeDiagonalP( vec_diag_F, (Params::EMethod)params.m_Method ) ); //\todo MUST NEGATE to draw... just convention
                    float c( vec_diag_P.Norm() / max_norm_f );
                    Vec2r v( length * mal::SafeNormalized( vec_diag_P ) );
                    DrawVector( sfr::Vec3(vec_diag_F.x(),vec_diag_F.y(),c), sfr::Vec3(v.x(),v.y(),0),
                                sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c,
                                                                  0.25f,
                                                                  0.25f + 0.75f*(1.0-c),
                                                     1 ),
                                                 1 ) );
                }
            }
            DrawSegment2( sfr::Vec2(-half_sizes.x(), half_sizes.y()),
                          sfr::Vec2( half_sizes.x(),-half_sizes.y()),
                          sfr::gfx::Style( sfr::gfx::Color(0.5,0.5,0.5,1), 2 ) );
        }
        // Sample Rotation field
        if( params.m_DrawFlags.Test(Params::eDraw_R_x) )
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * params.m_YoungModulus );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Mat2x2r R;
                    element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R );
                    // Draw inferred rotation
                    DrawRefSys2( X[nid], R, length, sfr::gfx::Style() );
                    //DrawVector2( sfr::Vec2(X[nid]), sfr::Vec2(length*mal::GColumn<0>(R)), sfr::gfx::Style( sfr::gfx::Color( 0.66f, 0.25f, 0.25f ), 1 ) );
                }
            }
        }
        // Sample Det(F)
        if( params.m_DrawFlags.Test(Params::eDraw_Det_F) )
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Mat2x2r F = element.Compute_F( X, params.m_FFM );
                    Real det_F( mal::Det(F) );
                    float c( det_F );/// max_det_F );
                    DrawVector( sfr::Vec3(X[nid].x(),X[nid].y(),0), sfr::Vec3(0,0,c),
                                sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c,
                                                                  0.25f,
                                                                  0.25f + 0.75f*(1.0-c),
                                                                  1 ),
                                                 1 ) );
                }
            }
        }

        // Draw QR-GSO critical point
        {
            Real den( element.m_invDm(0,0) + element.m_invDm(1,0) );
            if( mal::Abs(den) > mal::Epsilon<Real>() )
            {
                Vec2r gso_critical_X0( (element.m_r[1].x() * element.m_invDm(0,0) + element.m_r[2].x() * element.m_invDm(1,0) ) / den,
                                       (element.m_r[1].y() * element.m_invDm(0,0) + element.m_r[2].y() * element.m_invDm(1,0) ) / den );
                Mat2x2r Ds( mal::GMat2x2_From_Columns( element.m_r[1] - gso_critical_X0, element.m_r[2] - gso_critical_X0 ) );
                Mat2x2r F( Ds * element.m_invDm );
                //TEMP: std::cout << mal::GColumn<0>(F) << std::endl;
                DrawPoint2( gso_critical_X0, sfr::gfx::Style(sfr::gfx::Color(1,0,1,1),20) );
            }
        }

        // Draw |dR_exact - dR_numerical| along radial dx[nid] directions
        if( params.m_DrawFlags.Test(Params::eDraw_Drag_dR) )
        {
            // Draw dR errors around drag point
            DrawPoint2( m_rAppTask.m_DragPoint, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),10) );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            X[nid] = Vec2r( m_rAppTask.m_DragPoint );
            const unsigned int cNumDirectionalDerivatives( 8 );
            for( unsigned int it_dd = 0; it_dd < cNumDirectionalDerivatives; it_dd++ )
            {
                Real angle_rad( mal::TwoPi<Real>() * float(it_dd) / cNumDirectionalDerivatives );
                Vec2r radial( Vec2r( mal::Cos( angle_rad ), mal::Sin( angle_rad ) ) );
                Vec2r perpendicular( mal::PerpendicularCW(radial) );
                //Vec2r dX[3] = { Vec2r(0,0), Vec2r(0,0), Vec2r(0,0) };
                Vec2r dX[3] = { Vec2r(0,0),
                                params.m_DiffR_dX12_Bias_Factor * perpendicular, //params.m_DiffR_dX12_Bias_Factor * mal::Normalized(radial+perpendicular),
                                -params.m_DiffR_dX12_Bias_Factor * perpendicular }; //params.m_DiffR_dX12_Bias_Factor * mal::Normalized(radial-perpendicular) };
                dX[nid] = radial;

                Mat2x2r F( element.Compute_F( X, Params::eFFM_None ) );
                Mat2x2r R;
                element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R ); //\note Computes R according to selected method
                Mat2x2r dR,dR_numerical;

                bool bCorrect_dR(false);
                switch( params.m_Method )
                {
                case Params::eMethod_Id: bCorrect_dR = true; dR = Mat2x2r::Identity(); break;
                case Params::eMethod_QR_XY: break;
                case Params::eMethod_QR_YX: break;
                case Params::eMethod_PD: break;
                case Params::eMethod_PD_Project: break;
                case Params::eMethod_PD_Reflect: break;
                case Params::eMethod_PD_Unified:
                    bCorrect_dR = Compute_dR_PD_U( element, params,
                                                   X[0], X[1], X[2],
                                                   dX[0], dX[1], dX[2],
                                                   element.m_invDm, F, R, mal::Transposed(R),
                                                   dR );
                    break;
                case Params::eMethod_PD_Fix: break;
                case Params::eMethod_PD_Project_Nearest: break;
                    // All SVD-based methods use the same dR from McAdams
                case Params::eMethod_PD_SVD:
                case Params::eMethod_IHFSDM:
                case Params::eMethod_ITF_LRM:
                case Params::eMethod_ECIE_CLRM:
                case Params::eMethod_ECIE_NHC0:
                case Params::eMethod_ECIE_NHC1:
                    bCorrect_dR = Compute_dR( element,
                                              dX[0], dX[1], dX[2],
                                              element.m_invDm, F, R, mal::Transposed(R),
                                              dR );
                    break;
                case Params::eMethod_PD_CLRM:
                    bCorrect_dR = Compute_dR_PD_U( element, params,
                                                   X[0], X[1], X[2],
                                                   dX[0], dX[1], dX[2],
                                                   element.m_invDm, F, R, mal::Transposed(R),
                                                   dR );
                    break;
                default: break;
                }

                bool bCorrect_dR_numerical = Compute_dR_Numerical( element, (Params::EMethod)params.m_Method,
                                                                   X[0], X[1], X[2],
                                                                   dX[0], dX[1], dX[2],
                                                                   element.m_invDm, params.m_DegenerateThresholdDetF,
                                                                   params.m_DiffR_Numerical_H,
                                                                   dR_numerical );
                if( bCorrect_dR && bCorrect_dR_numerical )
                {
                    Real norm_dR( mal::NormF( dR ) );
                    Real norm_dR_n( mal::NormF( dR_numerical ) );
                    Real error_abs( mal::NormF( dR - dR_numerical ) );
                    Real error_rel( error_abs / mal::Max( norm_dR, norm_dR_n ) );
                    Real norm_dx( mal::Sqrt( mal::NormSq(dX[0]) + mal::NormSq(dX[1]) + mal::NormSq(dX[2]) ) );
                    char str[128];
                    sprintf( str, "err %3.3f,%3.3f |dR|/|dx|=%1.2f", error_abs, error_rel, norm_dR / norm_dx );

                    /*
                    Real eigen_value_Ke;
                    Vec6r eigen_vector_Ke;
                    mal::GComputeLargestEigenvalue_Symmetric( element.m_Ke, eigen_value_Ke, eigen_vector_Ke, 20, 0.01f );
                    Mat6x6r Re(Real(0));
                    mal::GSetRange<0,0,1,1>( Re, R );
                    mal::GSetRange<2,2,3,3>( Re, R );
                    mal::GSetRange<4,4,5,5>( Re, R );
                    Real eigen_value_KeR;
                    Vec6r eigen_vector_KeR;
                    mal::GComputeLargestEigenvalue_Symmetric( Re*element.m_Ke*mal::Transposed(Re), eigen_value_KeR, eigen_vector_KeR, 20, 0.01f );
                    Mat6x6r dRe(Real(0));
                    mal::GSetRange<0,0,1,1>( dRe, dR );
                    mal::GSetRange<2,2,3,3>( dRe, dR );
                    mal::GSetRange<4,4,5,5>( dRe, dR );
                    Real eigen_value_dR;
                    Vec6r eigen_vector_dR;
                    mal::GComputeLargestEigenvalue_Symmetric( dRe*element.m_Ke*mal::Transposed(Re)
                                                              + mal::Transposed( dRe*element.m_Ke*mal::Transposed(Re) ),
                                                              eigen_value_dR, eigen_vector_dR, 10, 0.01f );
                    sprintf( str, "err %3.3f,%3.3f ev(Ke)=%1.2f ev(KeR)=%1.2f, ev(dR)=%1.2f", error_abs, error_rel, eigen_value_Ke, eigen_value_KeR, eigen_value_dR );
                    */

                    if( error_rel < 0.05 || error_abs < 0.05 )
                    {
                        DrawLabel2( X[nid] + dX[nid], str, sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),1) );
                        DrawVector2( X[nid] + dX[nid], mal::GColumn<0>( dR ), sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),1) );
                        DrawVector2( X[nid] + dX[nid], mal::GColumn<0>( dR_numerical ), sfr::gfx::Style(sfr::gfx::Color(0,0.5,0,1),1) );
                    }
                    else if( error_rel * error_abs < 0.01 ) //"total" error small, even if both above minimum
                    {
                        DrawLabel2( X[nid] + dX[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,1,0,1),1) );
                        DrawVector2( X[nid] + dX[nid], mal::GColumn<0>( dR ), sfr::gfx::Style(sfr::gfx::Color(1,1,0,1),1) );
                        DrawVector2( X[nid] + dX[nid], mal::GColumn<0>( dR_numerical ), sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,1),1) );
                    }
                    else //Large error
                    {
                        DrawLabel2( X[nid] + dX[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),1) );
                        DrawVector2( X[nid] + dX[nid], mal::GColumn<0>( dR ), sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),1) );
                        DrawVector2( X[nid] + dX[nid], mal::GColumn<0>( dR_numerical ), sfr::gfx::Style(sfr::gfx::Color(0.5,0.25,0.25,1),1) );
                    }
                }
                else if( !bCorrect_dR && !bCorrect_dR_numerical )
                    DrawLabel2( X[nid] + dX[nid], "dR && dR_n", sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),1) );
                else if( !bCorrect_dR )
                    DrawLabel2( X[nid] + dX[nid], "dR", sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),1) );
                else if( !bCorrect_dR_numerical )
                    DrawLabel2( X[nid] + dX[nid], "dR_n", sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),1) );
            }
        }

        // Draw |df_exact - df_numerical| along radial dx[nid] directions
        if( params.m_DrawFlags.Test(Params::eDraw_Drag_df)
            /*&& ( params.m_Method == Params::eMethod_PD_Unified
              || params.m_Method == Params::eMethod_PD_CLRM ) */ )
        {
            // Draw dR errors around drag point
            DrawPoint2( m_rAppTask.m_DragPoint, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),10) );
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            X[nid] = Vec2r( m_rAppTask.m_DragPoint );
            const unsigned int cNumDirectionalDerivatives( 8 );
            for( unsigned int it_dd = 0; it_dd < cNumDirectionalDerivatives; it_dd++ )
            {
                Real angle_rad( mal::TwoPi<Real>() * float(it_dd) / cNumDirectionalDerivatives );
                Vec2r radial( Vec2r( mal::Cos( angle_rad ), mal::Sin( angle_rad ) ) );
                Vec2r perpendicular( mal::PerpendicularCW(radial) );
                //Vec2r dX[3] = { Vec2r(0,0), Vec2r(0,0), Vec2r(0,0) };
                Vec2r dX[3] = { Vec2r(0,0),
                                params.m_DiffR_dX12_Bias_Factor * perpendicular, //params.m_DiffR_dX12_Bias_Factor * mal::Normalized(radial+perpendicular),
                                -params.m_DiffR_dX12_Bias_Factor * perpendicular }; //params.m_DiffR_dX12_Bias_Factor * mal::Normalized(radial-perpendicular) };
                dX[nid] = radial;

                // F unmodified here!
                Mat2x2r F( element.Compute_F( X, Params::eFFM_None ) );
                Mat2x2r R;
                element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R ); //\note Computes R according to selected method

                Mat2x2r dR,dR_numerical;
                /*
                bool bCorrect_dR = Compute_dR_PD_U( element, params,
                                                    X[0], X[1], X[2],
                                                    dX[0], dX[1], dX[2],
                                                    element.m_invDm, F, R, mal::Transposed(R),
                                                    dR );
                */
                bool bCorrect_dR(false);
                switch( params.m_Method )
                {
                case Params::eMethod_Id: bCorrect_dR = true; dR = Mat2x2r::Identity(); break;
                case Params::eMethod_QR_XY:
                    bCorrect_dR = Compute_dR_QR_XY( dX[0], dX[1], dX[2],
                                                    element.m_invDm, F,
                                                    dR );
                    break;
                case Params::eMethod_QR_YX: break;
                case Params::eMethod_PD: break;
                case Params::eMethod_PD_Project: break;
                case Params::eMethod_PD_Reflect: break;
                case Params::eMethod_PD_Unified:
                case Params::eMethod_DAPD_NS:
                case Params::eMethod_DAPD_SS:
                case Params::eMethod_DAPD_EX:
                    bCorrect_dR = Compute_dR_PD_U( element, params,
                                                   X[0], X[1], X[2],
                                                   dX[0], dX[1], dX[2],
                                                   element.m_invDm, F, R, mal::Transposed(R),
                                                   dR );
                    break;
                case Params::eMethod_PD_Fix: break;
                case Params::eMethod_PD_Project_Nearest: break;
                    // All SVD-based methods use the same dR from McAdams
                case Params::eMethod_PD_SVD:
                case Params::eMethod_IHFSDM:
                case Params::eMethod_ITF_LRM:
                case Params::eMethod_ECIE_CLRM:
                case Params::eMethod_ECIE_NHC0:
                case Params::eMethod_ECIE_NHC1:
                    bCorrect_dR = Compute_dR( element,
                                              dX[0], dX[1], dX[2],
                                              element.m_invDm, F, R, mal::Transposed(R),
                                              dR );
                    break;
                case Params::eMethod_PD_CLRM:
                    bCorrect_dR = Compute_dR_PD_U( element, params,
                                                   X[0], X[1], X[2],
                                                   dX[0], dX[1], dX[2],
                                                   element.m_invDm, F, R, mal::Transposed(R),
                                                   dR );
                    break;
                default: break;
                }
                bool bCorrect_dR_numerical = Compute_dR_Numerical( element, (Params::EMethod)params.m_Method,
                                                                   X[0], X[1], X[2],
                                                                   dX[0], dX[1], dX[2],
                                                                   element.m_invDm, params.m_DegenerateThresholdDetF,
                                                                   params.m_DiffR_Numerical_H,
                                                                   dR_numerical );
                if( bCorrect_dR && bCorrect_dR_numerical )
                {
                    Vec6r df( Vec6r::Zero() );
                    switch( params.m_Method )
                    {
                        // All LINEAR
                    case Params::eMethod_Id:
                    case Params::eMethod_QR_XY:
                    case Params::eMethod_QR_YX:
                    case Params::eMethod_PD:
                    case Params::eMethod_PD_Project:
                    case Params::eMethod_PD_Reflect:
                    case Params::eMethod_PD_Unified:
                    case Params::eMethod_PD_Fix:
                    case Params::eMethod_PD_Project_Nearest:
                    case Params::eMethod_PD_SVD:
                    case Params::eMethod_IHFSDM:
                    case Params::eMethod_ITF_LRM:
                        //df = element.Compute_C_df_x( X, dX, R, dR );
                        df = element.Compute_H_LRM_df_x_SymmStrain( X, dX, R, dR );
                        break;
                    case Params::eMethod_PD_CLRM:
                    case Params::eMethod_ECIE_CLRM:
                        df = element.Compute_H_CLRM_df_x( X, dX, R, dR );
                        break;
                    case Params::eMethod_ECIE_NHC0:
                    case Params::eMethod_ECIE_NHC1:
                        // NOT SUPPORTED
                        break;
                    case Params::eMethod_DAPD_NS:
                    case Params::eMethod_DAPD_SS:
                    case Params::eMethod_DAPD_EX:
                        df = element.Compute_df_Numerical( X, dX, (Params::EMethod)params.m_Method, params.m_DiffR_Numerical_H );
                    default: break;
                    }
                    Vec6r df_numerical = element.Compute_df_Numerical( X, dX, (Params::EMethod)params.m_Method, params.m_DiffR_Numerical_H );
                    Real norm_df( mal::Norm( df ) );
                    Real norm_df_n( mal::Norm( df_numerical ) );
                    Real error_abs( mal::Norm( df - df_numerical ) );
                    Real error_rel( error_abs / mal::Max( norm_df, norm_df_n ) );
                    char str[128];
                    sprintf( str, "err %3.3f,%3.3f"/*, norm %3.3f,%3.3f"*/, error_abs, error_rel );//, norm_dR, norm_dR_n );
                    if( error_rel < 0.05 || error_abs < 0.05 )
                        DrawLabel2( X[nid] + dX[nid], str, sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),1) );
                    else if( error_rel * error_abs < 0.01 ) //"total" error small, even if both above minimum
                        DrawLabel2( X[nid] + dX[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,1,0,1),1) );
                    else //Large error
                        DrawLabel2( X[nid] + dX[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),1) );
                    // Draw normalized df (shows absolute angular error and relative length error)
                    df /= mal::Max( norm_df, norm_df_n );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<0,1>( df ), sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),2) );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<2,3>( df ), sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),2) );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<4,5>( df ), sfr::gfx::Style(sfr::gfx::Color(0,0,1,1),2) );
                    df_numerical /= mal::Max( norm_df, norm_df_n );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<0,1>( df_numerical ), sfr::gfx::Style(sfr::gfx::Color(0.5,0,0,1),2) );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<2,3>( df_numerical ), sfr::gfx::Style(sfr::gfx::Color(0,0.5,0,1),2) );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<4,5>( df_numerical ), sfr::gfx::Style(sfr::gfx::Color(0,0,0.5,1),2) );
                    /*TEMP: Show also df_dR_numerical, not used in error computation
                    Vec6r df_dR_numerical;
                    if( params.m_Method == Params::eMethod_PD_Unified )
                        df_dR_numerical = element.Compute_C_df_x( X, dX, R, dR_numerical );
                    else //params.m_Method == Params::eMethod_PD_CLRM
                        df_dR_numerical = element.Compute_H_CLRM_df_x( X, dX, R, dR_numerical );
                    df_dR_numerical /= 0.5*mal::Max( norm_df, norm_df_n );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<0,1>( df_dR_numerical ), sfr::gfx::Style(sfr::gfx::Color(0.25,0,0),1) );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<2,3>( df_dR_numerical ), sfr::gfx::Style(sfr::gfx::Color(0,0.25,0),1) );
                    DrawVector2( X[nid] + dX[nid], mal::GRange<4,5>( df_dR_numerical ), sfr::gfx::Style(sfr::gfx::Color(0,0,0.25),1) );
                    */
                }
                else if( !bCorrect_dR && !bCorrect_dR_numerical )
                    DrawLabel2( X[nid] + dX[nid], "df && df_n", sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),1) );
                else if( !bCorrect_dR )
                    DrawLabel2( X[nid] + dX[nid], "df", sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),1) );
                else if( !bCorrect_dR_numerical )
                    DrawLabel2( X[nid] + dX[nid], "df_n", sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),1) );
            }
        }

        // Assemble K = Pf_Px
        if( params.m_DrawFlags.Test(Params::eDraw_Drag_df) )
        {
            Vec2r X[3] = { element.m_X[0], element.m_X[1], element.m_X[2] };
            X[nid] = Vec2r( m_rAppTask.m_DragPoint );
            Mat2x2r R;
            element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R ); //\note Computes R according to selected method
            Mat2x2r F( element.Compute_F( X, Params::eFFM_None ) );

            Mat6x6r Pf_Px( Mat6x6r::Identity() );
#ifdef __ENABLE_COMPARE_SVD1_AND_LRM
            Mat6x6r Pf_Px_SVD_exact( Mat6x6r::Identity() );
            Mat6x6r Pf_Px_LRM_exact( Mat6x6r::Identity() );
            Mat6x6r Pf_Px_SVD_trunc( Mat6x6r::Identity() );
            Mat6x6r Pf_Px_LRM_trunc( Mat6x6r::Identity() );
            /*
            Mat6x6r Pf_Px_SVD_symm( Mat6x6r::Identity() );
            Mat6x6r Pf_Px_LRM_symm( Mat6x6r::Identity() );
            */
            Mat6x6r Pf_Px_CLRM_exact( Mat6x6r::Identity() );
            Mat6x6r Pf_Px_CLRM_trunc( Mat6x6r::Identity() );
#endif
            for( int i=0; i<6; i++ )
            {
                Vec6r dx( Vec6r::Zero() );
                dx[i] = 1;
                Vec2r dX[3] = { mal::GRange<0,1>(dx),
                                mal::GRange<2,3>(dx),
                                mal::GRange<4,5>(dx) };
                Mat2x2r dR( Mat2x2r::Zero() );
                bool bCorrect_dR(false);
                {
                    switch( params.m_Method )
                    {
                    case Params::eMethod_Id: bCorrect_dR = true; dR = Mat2x2r::Identity(); break;
                    case Params::eMethod_QR_XY: break;
                    case Params::eMethod_QR_YX: break;
                    case Params::eMethod_PD: break;
                    case Params::eMethod_PD_Project: break;
                    case Params::eMethod_PD_Reflect: break;
                    case Params::eMethod_PD_Unified:
                    case Params::eMethod_DAPD_NS:
                    case Params::eMethod_DAPD_SS:
                    case Params::eMethod_DAPD_EX:
                        bCorrect_dR = Compute_dR_PD_U( element, params,
                                                       X[0], X[1], X[2],
                                                       dX[0], dX[1], dX[2],
                                                       element.m_invDm, F, R, mal::Transposed(R),
                                                       dR );
                        break;
                    case Params::eMethod_PD_Fix: break;
                    case Params::eMethod_PD_Project_Nearest: break;
                        // All SVD-based methods use the same dR from McAdams
                    case Params::eMethod_PD_SVD:
                    case Params::eMethod_IHFSDM:
                    case Params::eMethod_ITF_LRM:
                    case Params::eMethod_ECIE_CLRM:
                    case Params::eMethod_ECIE_NHC0:
                    case Params::eMethod_ECIE_NHC1:
                        bCorrect_dR = Compute_dR( element,
                                                  dX[0], dX[1], dX[2],
                                                  element.m_invDm, F, R, mal::Transposed(R),
                                                  dR );
                        break;
                    case Params::eMethod_PD_CLRM:
                        bCorrect_dR = Compute_dR_PD_U( element, params,
                                                       X[0], X[1], X[2],
                                                       dX[0], dX[1], dX[2],
                                                       element.m_invDm, F, R, mal::Transposed(R),
                                                       dR );
                        break;
                    default: break;
                    }
                }

                if( params.m_Method == Params::eMethod_DAPD_NS
                    || params.m_Method == Params::eMethod_DAPD_SS
                    || params.m_Method == Params::eMethod_DAPD_EX )
                {
                    Vec6r df = element.Compute_df_Numerical( X, dX, (Params::EMethod)params.m_Method, params.m_DiffR_Numerical_H );
                    mal::SetColumn( Pf_Px, i, df );
                    bCorrect_dR = true;
                }
                else if( bCorrect_dR )
                {
                    // Vec6r df = element.Compute_C_df_x( X, dX, R, dR );
                    //Vec6r df = element.Compute_H_LRM_df_x( X, dX, R, Mat2x2r::Zero() );//dR );
                    // Vec6r df = element.Compute_H_LRM_df_x_SymmStrain( X, dX, R, Mat2x2r::Zero() );//dR );
                    Vec6r df = element.Compute_H_LRM_df_x_SymmStrain( X, dX, R, dR );
                    mal::SetColumn( Pf_Px, i, df );
                }

#ifdef __ENABLE_COMPARE_SVD1_AND_LRM
                if( bCorrect_dR
                    && (params.m_Method == Params::eMethod_PD_SVD
                        || params.m_Method == Params::eMethod_PD_Unified) )
                {
                    Vec6r df;
                    // Exact df
                    df = element.Compute_C_df_x( X, dX, R, dR );
                    mal::SetColumn( Pf_Px_SVD_exact, i, df );
                    df = element.Compute_H_LRM_df_x( X, dX, R, dR );
                    mal::SetColumn( Pf_Px_LRM_exact, i, df );
                    // Truncated df
                    df = element.Compute_C_df_x( X, dX, R, Mat2x2r::Zero() );
                    mal::SetColumn( Pf_Px_SVD_trunc, i, df );
                    df = element.Compute_H_LRM_df_x( X, dX, R, Mat2x2r::Zero() );
                    mal::SetColumn( Pf_Px_LRM_trunc, i, df );
                    // CLRM df
                    df = element.Compute_H_CLRM_df_x( X, dX, R, dR );
                    mal::SetColumn( Pf_Px_CLRM_exact, i, df );
                    df = element.Compute_H_CLRM_df_x( X, dX, R, Mat2x2r::Zero() );
                    mal::SetColumn( Pf_Px_CLRM_trunc, i, df );
                }
#endif

            }
            // Analyze Pf_Px
            Real asymmetry( mal::NormInf( Pf_Px - mal::Transposed(Pf_Px) ) );
            Real min_max_eigenvalue_ratio( mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px, 25, 0.001f) );
            char str[128];
            sprintf( str, "Asym %3.1f m/M %1.3f", asymmetry, min_max_eigenvalue_ratio );
            DrawLabel2( X[nid], str, sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,1,1),1) );
            std::cout << "---------------------------------------------------------------" << std::endl;
            std::cout << "Asym = " << asymmetry << std::endl;
                /*
                      << "Pf_Px = " //<< std::endl
                      << Pf_Px << std::endl
                      << "|A-At|" << std::endl
                      << (Pf_Px - mal::Transposed(Pf_Px))
                      << std::endl;
                */
#ifdef __ENABLE_COMPARE_SVD1_AND_LRM
            if( params.m_Method == Params::eMethod_PD_SVD || params.m_Method == Params::eMethod_PD_Unified )
            {
                std::cout << "----------LRM-----------" << std::endl;
                std::cout << "|LRM_E| = " << mal::NormF( Pf_Px_LRM_exact ) << std::endl;
                std::cout << "|LRM_T| = " << mal::NormF( Pf_Px_LRM_trunc ) << std::endl;
                std::cout << "|SVD_E| = " << mal::NormF( Pf_Px_SVD_exact ) << std::endl;
                std::cout << "|SVD_T| = " << mal::NormF( Pf_Px_SVD_trunc ) << std::endl;
                std::cout << "|Pf_Px| = " << mal::NormF( Pf_Px ) << std::endl;
                std::cout << "|SVD-LRM|_E = " << mal::NormInf( Pf_Px_SVD_exact - Pf_Px_LRM_exact ) / mal::Max( mal::NormInf(Pf_Px_SVD_exact), mal::NormInf(Pf_Px_LRM_exact) ) << std::endl
                          << "Pf_Px_SVD_E " << "asym " << mal::NormInf( Pf_Px_SVD_exact - mal::Transposed(Pf_Px_SVD_exact) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px_SVD_exact, 25, 0.001f) << std::endl
                          << "Pf_Px_LRM_E " << "asym " << mal::NormInf( Pf_Px_LRM_exact - mal::Transposed(Pf_Px_LRM_exact) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px_LRM_exact, 25, 0.001f) << std::endl;
                std::cout << "|SVD-LRM|_T = " << mal::NormInf( Pf_Px_SVD_trunc - Pf_Px_LRM_trunc ) / mal::Max( mal::NormInf(Pf_Px_SVD_trunc), mal::NormInf(Pf_Px_LRM_trunc) ) << std::endl
                          << "Pf_Px_SVD_T " << "asym " << mal::NormInf( Pf_Px_SVD_trunc - mal::Transposed(Pf_Px_SVD_trunc) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px_SVD_trunc, 25, 0.001f) << std::endl
                          << "Pf_Px_LRM_T " << "asym " << mal::NormInf( Pf_Px_LRM_trunc - mal::Transposed(Pf_Px_LRM_trunc) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px_LRM_trunc, 25, 0.001f) << std::endl
                          << "Pf_Px       " << "asym " << mal::NormInf( Pf_Px - mal::Transposed(Pf_Px) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px, 25, 0.001f) << std::endl;

#ifdef __ENABLE_COMPUTE_APPROX_JACOBIAN_ERROR
                // std::cout << "...TESTING GSL..." << std::endl;
                //test_fucking_gsl_invert();

                std::cout << "...Computing rho..." << std::endl;
                std::cout << "LRM_E = " << Pf_Px_LRM_exact << std::endl;
                std::cout << "LRM_T = " << Pf_Px_LRM_trunc << std::endl;
                std::cout << "SVD_T = " << Pf_Px_SVD_trunc << std::endl;
                //std::cout << "rho(SVD_T) = " << ComputeSpectralRadius_Nonsymmetric( Pf_Px_SVD_trunc ) << std::endl;
                std::cout << "Inv(SVD_T) = " << ComputeInverse(Pf_Px_SVD_trunc) << std::endl;
                Real rho_LRM( ComputeSpectralRadius_Nonsymmetric( Mat6x6r::Identity() - ComputeInverse( Pf_Px_LRM_trunc ) * Pf_Px_LRM_exact ) );
                Real rho_SVD( ComputeSpectralRadius_Nonsymmetric( Mat6x6r::Identity() - ComputeInverse( Pf_Px_SVD_trunc ) * Pf_Px_SVD_exact ) );
                //Real rho_TEST( ComputeSpectralRadius_Nonsymmetric( Mat6x6r::Identity() - ComputeInverse( 0.01*Mat6x6r::Identity() + Pf_Px_LRM_exact ) * (0.01*Mat6x6r::Identity()+Pf_Px_LRM_exact) ) );
                std::cout << "|LRM_E-LRM_T| = " << rho_LRM << std::endl;
                std::cout << "|SVD_E-SVD_T| = " << rho_SVD << std::endl;
#else
                std::cout << "|LRM_E-LRM_T| = " << mal::NormInf( Pf_Px_LRM_exact - Pf_Px_LRM_trunc ) / mal::Max( mal::NormInf(Pf_Px_LRM_exact), mal::NormInf(Pf_Px_LRM_trunc) ) << std::endl;
                std::cout << "|LRM_E-SVD_T| = " << mal::NormInf( Pf_Px_LRM_exact - Pf_Px_SVD_trunc ) / mal::Max( mal::NormInf(Pf_Px_LRM_exact), mal::NormInf(Pf_Px_SVD_trunc) ) << std::endl;
#endif

                std::cout << "----------CLRM-----------" << std::endl;
                std::cout << "|CLRM_E| = " << mal::NormF( Pf_Px_CLRM_exact ) << std::endl;
                std::cout << "|CLRM_T| = " << mal::NormF( Pf_Px_CLRM_trunc ) << std::endl;
                std::cout << "Pf_Px_CLRM_E " << "asym " << mal::NormInf( Pf_Px_CLRM_exact - mal::Transposed(Pf_Px_CLRM_exact) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px_CLRM_exact, 25, 0.001f) << std::endl
                          << "Pf_Px_CLRM_T " << "asym " << mal::NormInf( Pf_Px_CLRM_trunc - mal::Transposed(Pf_Px_CLRM_trunc) ) << " ratio " << mal::GComputeMinMaxEigenvalueRatio_Symmetric(Pf_Px_CLRM_trunc, 25, 0.001f) << std::endl;

#ifdef __ENABLE_COMPUTE_APPROX_JACOBIAN_ERROR
                Real rho_CLRM( ComputeSpectralRadius_Nonsymmetric( Mat6x6r::Identity() - ComputeInverse( Pf_Px_CLRM_trunc ) * Pf_Px_CLRM_exact ) );
                std::cout << "|CLRM_E-CLRM_T| = " << rho_CLRM << std::endl;
#else
                std::cout << "|CLRM_E-CLRM_T| = " << mal::NormInf( Pf_Px_CLRM_exact - Pf_Px_CLRM_trunc ) / mal::Max( mal::NormInf(Pf_Px_CLRM_exact), mal::NormInf(Pf_Px_CLRM_trunc) ) << std::endl;
#endif
            }
#endif
        }

        return true;
    }

    bool Render3D()
    {
        const Params &params( m_rAppTask.GetParams() );
        int nid( params.m_NID );

        if( params.m_DrawFlags.Test(Params::eDraw_Axis) )
            DrawRefSys( Vec3f::Zero(), Mat3x3f::Identity(), 25.0f, sfr::gfx::Style() );

        float grid_size = 50.0f;
        float half_grid_size = 0.5f * grid_size;
        if( params.m_DrawFlags.Test(Params::eDraw_Grid) )
            DrawGrid( Vec3f(-half_grid_size,-half_grid_size,0),
                      Vec3f(1,0,0), Vec3f(0,1,0),
                      grid_size, grid_size,
                      10, 10,
                      sfr::gfx::Style(sfr::gfx::Color(0.1,0.1,0.5,1),1) );

        // Element
        Element3D &element( m_rAppTask.GetElement3D() );
        if( params.m_DrawFlags.Test(Params::eDraw_Element) )
        {
            // Rest
            DrawSegment( element.m_r[0], element.m_r[1], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment( element.m_r[0], element.m_r[2], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment( element.m_r[0], element.m_r[3], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment( element.m_r[1], element.m_r[2], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment( element.m_r[1], element.m_r[3], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawSegment( element.m_r[2], element.m_r[3], sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
            DrawPoint( element.m_r[0], sfr::gfx::Style(sfr::gfx::Color(1,0.25,0.25,1),8) );
            DrawPoint( element.m_r[1], sfr::gfx::Style(sfr::gfx::Color(0.25,1,0.25,1),8) );
            DrawPoint( element.m_r[2], sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,1,1),8) );
            DrawPoint( element.m_r[3], sfr::gfx::Style(sfr::gfx::Color(1.0,1,0.25,1),8) );
#ifdef __TODO_PORT_TO_3D
            // Rotated
            DrawSegment2( element.m_X[0], element.m_X[1], sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,5,1),1) );
            DrawSegment2( element.m_X[1], element.m_X[2], sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,5,1),1) );
            DrawSegment2( element.m_X[2], element.m_X[0], sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0,5,1),1) );
            DrawPoint2( element.m_X[0], sfr::gfx::Style(sfr::gfx::Color(0.5,0.25,0.25,1),6) );
            DrawPoint2( element.m_X[1], sfr::gfx::Style(sfr::gfx::Color(0.25,0.5,0.25,1),6) );
            DrawPoint2( element.m_X[2], sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,0.5,1),6) );
#endif
        }

        // Draw Drag stuff
        if( true )//m_rAppTask.m_bIsEnabledDrag )
        {
            Vec3r X[4] = { element.m_X[0], element.m_X[1], element.m_X[2], element.m_X[3] };
            X[nid] = Vec3r( m_rAppTask.m_DragPoint[0], m_rAppTask.m_DragPoint[1], params.m_DragPointZ * params.m_HalfSizeZ );
            DrawPoint( X[nid], sfr::gfx::Style(sfr::gfx::Color(1,0.5,0.5,1),10) );

            /*TEMP: Draw SVD details
            // Compute F
            Mat2x2r Q( mal::GMat2x2_From_Columns( X[1] - X[0], X[2] - X[0] ) );
            Mat2x2r F( Q * mal::Inverse( mal::GMat2x2_From_Columns( element.m_r[1] - element.m_r[0], element.m_r[2] - element.m_r[0] ) ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            // Draw eigenvalues
            char str[256];
            sprintf( str, "%f + %f = %f", vec_diag_F[0], vec_diag_F[1], vec_diag_F[0] + vec_diag_F[1] );
            DrawLabel2( X[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,1,1),1) );
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
            if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
            if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
            }
            Mat2x2r R( U*Vt );
            DrawRefSys2( X[nid], R, 1, sfr::gfx::Style() );
            DrawRefSys2( X[nid], U, 0.66, sfr::gfx::Style() );
            DrawRefSys2( X[nid], Vt.Transposed(), 0.33, sfr::gfx::Style() );
            */


            // Drag force/rotation
            if( params.m_DrawFlags.Test( Params::eDraw_Drag_f0
                                         | Params::eDraw_Drag_R
                                         | Params::eDraw_Drag_EE0 ) )
            {
                Mat3x3f R;
                Real elastic_energy(0);
                Vec12r f( element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R, &elastic_energy ) );
                if( params.m_DrawFlags.Test(Params::eDraw_Drag_f0) )
                {
                    // Draw elastic force
                    Vec3f f0( mal::GRange<0,2>(f) );
                    Vec3f f1( mal::GRange<3,5>(f) );
                    Vec3f f2( mal::GRange<6,8>(f) );
                    Vec3f f3( mal::GRange<9,11>(f) );
                    float max_norm_f( mal::Max( mal::Max( f0.Norm(), f1.Norm() ), mal::Max( f2.Norm(), f3.Norm() ) ) );
                    if( max_norm_f > mal::Epsilon<float>() );
                    {
                        Vec3f v0( f0 / max_norm_f );
                        Vec3f v1( f1 / max_norm_f );
                        Vec3f v2( f2 / max_norm_f );
                        Vec3f v3( f3 / max_norm_f );
                        DrawVector( X[0], v0, sfr::gfx::Style( sfr::gfx::Color( 1.0f, 0.25f, 0.25f,1 ), 2 ) );
                        DrawVector( X[1], v1, sfr::gfx::Style( sfr::gfx::Color( 0.25f, 1.0f, 0.25f,1 ), 2 ) );
                        DrawVector( X[2], v2, sfr::gfx::Style( sfr::gfx::Color( 0.25f, 0.25f, 1.0f,1 ), 2 ) );
                        DrawVector( X[3], v3, sfr::gfx::Style( sfr::gfx::Color( 1.0f, 1.0f, 0.25f,1 ), 2 ) );
                    }
                }
                if( params.m_DrawFlags.Test(Params::eDraw_Drag_R) )
                {
                    // Draw inferred rotation
                    DrawRefSys( X[nid], R, 1, sfr::gfx::Style() );
                }
                if( params.m_DrawFlags.Test(Params::eDraw_Drag_EE0) )
                {
                    char str[128];
                    sprintf( str, "EE = %f", elastic_energy );
                    DrawLabel( X[nid], str, sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),1) );
                }
            }

            // Drag trajectory
            if( params.m_DrawFlags.Test(Params::eDraw_Drag_t0) )
            {
                // Init statistics
                Real elastic_energy(0);
                Real *p_elastic_energy(0);
                unsigned int num_samples( mal::Ceil( params.m_Duration / params.m_TimeStep ) );
                if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) ) { m_rAppTask.m_TotalEnergyST.Begin( num_samples ); p_elastic_energy = &elastic_energy; }
                if( params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) ) { m_rAppTask.m_ElasticEnergyST.Begin( num_samples ); p_elastic_energy = &elastic_energy; }
                if( params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) ) m_rAppTask.m_KineticEnergyST.Begin( num_samples );
                if( params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) ) m_rAppTask.m_AngularMomentumST.Begin( num_samples );
                if( params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) ) m_rAppTask.m_DeformationRatioST.Begin( num_samples );

                // Init trajectory
                Vec3r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY, params.m_HalfSizeZ );
                Mat3x3r R;
                Vec3r prevX[4] = { X[0], X[1], X[2], X[3] };
                Vec3r velX[4] = { Vec3r::Zero(), Vec3r::Zero(), Vec3r::Zero(), Vec3r::Zero() };
                Real t(0);
                Real inv_node_mass( Real(4)/params.m_Mass );
                Real damping_coeff( params.m_DampingRatio * 2 * mal::Sqrt( params.m_YoungModulus * params.m_Mass ) );
                // Init degeneration tracking and NoC
                Real det_F( mal::Det( element.Compute_F( X, Params::eFFM_None ) ) );
                bool bDegenerate( det_F <= params.m_DegenerateThresholdDetF );
                int trajectory_noc( bDegenerate ? params.m_NoC : -1 ); //Assume NoC = 0 if starts degenerate
                do
                {
                    Vec12r f( element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R, p_elastic_energy ) );
                    // Draw elastic force
                    Vec3r f0( mal::GRange<0,2>(f) );
                    Vec3r f1( mal::GRange<3,5>(f) );
                    Vec3r f2( mal::GRange<6,8>(f) );
                    Vec3r f3( mal::GRange<9,11>(f) );
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x0) )
                    {
                        velX[0] += params.m_TimeStep * inv_node_mass * (f0 - damping_coeff * velX[0]);
                        X[0] += params.m_TimeStep * velX[0];
                        DrawSegment( prevX[0], X[0], sfr::gfx::Style( sfr::gfx::Color( 1.0f, 0.25f, 0.25f,1 ), 2 ) );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x1) )
                    {
                        velX[1] += params.m_TimeStep * inv_node_mass * (f1 - damping_coeff * velX[1]);
                        X[1] += params.m_TimeStep * velX[1];
                        DrawSegment( prevX[1], X[1], sfr::gfx::Style( sfr::gfx::Color( 0.25f, 1.0f, 0.25f, 1 ), 2 ) );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x2) )
                    {
                        velX[2] += params.m_TimeStep * inv_node_mass * (f2 - damping_coeff * velX[2]);
                        X[2] += params.m_TimeStep * velX[2];
                        DrawSegment( prevX[2], X[2], sfr::gfx::Style( sfr::gfx::Color( 0.25f, 0.25f, 1.0f, 1 ), 2 ) );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_Free_x3) )
                    {
                        velX[3] += params.m_TimeStep * inv_node_mass * (f3 - damping_coeff * velX[3]);
                        X[3] += params.m_TimeStep * velX[3];
                        DrawSegment( prevX[3], X[3], sfr::gfx::Style( sfr::gfx::Color( 1.0f, 1.0f, 0.25f, 1 ), 2 ) );
                    }
                    t += params.m_TimeStep;

                    // Compute ToC->NoC
                    det_F = mal::Det( element.Compute_F( X, Params::eFFM_None ) );
                    /*\todo
                    if( !bDegenerate && det_F <= params.m_DegenerateThresholdDetF )
                    {
                        bDegenerate = true;
                        if( params.m_bEnableNoC ) trajectory_noc = ComputeNoC( prevX, X, params.m_DegenerateThresholdDetF, element.m_Volume );
                        else trajectory_noc = params.m_NoC;
                    }
                    else if( det_F > params.m_DegenerateThresholdDetF )
                    {
                        bDegenerate = false;
                        trajectory_noc = -1;
                    }
                    element.m_NoC = trajectory_noc;
                    */

                    // Advance
                    prevX[0] = X[0];
                    prevX[1] = X[1];
                    prevX[2] = X[2];
                    prevX[3] = X[3];

                    // Gather statistics
                    if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy | Params::eTF_ElasticEnergy | Params::eTF_KineticEnergy ) )
                    {
                        Real kinetic_energy( Real(0.5) * params.m_Mass * ( mal::NormSq(velX[0]) + mal::NormSq(velX[1]) + mal::NormSq(velX[2]) + mal::NormSq(velX[2]) ) ); //Lumped KE
                        /*\todo 3D!!! Consistent KE seems more accurate than lumped M-KE... both behave similarly but the scale of consistent KE is closer to EE, which is expected...
                        Real kinetic_energy( ( params.m_Mass / Real(12) )
                                             * ( mal::NormSq(velX[0]) + mal::NormSq(velX[1]) + mal::NormSq(velX[2])
                                                 + velX[0].x()*velX[1].x() + velX[1].x()*velX[2].x() + velX[0].x()*velX[2].x()
                                                 + velX[0].y()*velX[1].y() + velX[1].y()*velX[2].y() + velX[0].y()*velX[2].y() ) );
                        */
                        if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) )
                            m_rAppTask.m_TotalEnergyST.Add( kinetic_energy + elastic_energy );
                        if( params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) )
                            m_rAppTask.m_ElasticEnergyST.Add( elastic_energy );
                        if( params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) )
                            m_rAppTask.m_KineticEnergyST.Add( kinetic_energy );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) )
                    {
                        Vec3r CoM( mal::Rcp<float>(4) * ( X[0] + X[2] + X[2] + X[3] ) );
                        m_rAppTask.m_AngularMomentumST.Add( mal::Norm( mal::Cross( velX[0], X[0] - CoM )
                                                                       + mal::Cross( velX[1], X[1] - CoM )
                                                                       + mal::Cross( velX[2], X[2] - CoM )
                                                                       + mal::Cross( velX[3], X[3] - CoM ) ) );
                    }
                    if( params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) )
                        m_rAppTask.m_DeformationRatioST.Add( det_F );
                } while( t < params.m_Duration );

                //Restore Element2D.NoC to Params.NoC, changed through GUI
                element.m_NoC = params.m_NoC;

                // End statistics and show results
                if( params.m_TrajectoryFlags.Test(Params::eTF_TotalEnergy) ) m_rAppTask.m_TotalEnergyST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_ElasticEnergy) ) m_rAppTask.m_ElasticEnergyST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_KineticEnergy) ) m_rAppTask.m_KineticEnergyST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_AngularMomentum) ) m_rAppTask.m_AngularMomentumST.End();
                if( params.m_TrajectoryFlags.Test(Params::eTF_DeformationRatio) ) m_rAppTask.m_DeformationRatioST.End();
                // Final configuration
                {
                    sfr::gfx::Color col( det_F > 0
                                         ? sfr::gfx::Color(0.25,0.75,0.25,1)
                                         : sfr::gfx::Color(0.75,0.25,0.25,1) );
                    DrawSegment( X[0], X[1], sfr::gfx::Style(col,1) );
                    DrawSegment( X[0], X[2], sfr::gfx::Style(col,1) );
                    DrawSegment( X[0], X[3], sfr::gfx::Style(col,1) );
                    DrawSegment( X[1], X[2], sfr::gfx::Style(col,1) );
                    DrawSegment( X[1], X[3], sfr::gfx::Style(col,1) );
                    DrawSegment( X[2], X[3], sfr::gfx::Style(col,1) );
                }
            }
        }

        // Sample force field
        if( params.m_DrawFlags.Test(Params::eDraw_f0_x)
            || params.m_DrawFlags.Test(Params::eDraw_f1_x)
            || params.m_DrawFlags.Test(Params::eDraw_f2_x)
            || params.m_DrawFlags.Test(Params::eDraw_f3_x) )
        {
            Vec3r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY, params.m_HalfSizeZ );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            int hcz( params.m_HalfCellsZ );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * params.m_YoungModulus );
            int idx( params.m_DrawFlags.Test(Params::eDraw_f0_x) ? 0
                     : params.m_DrawFlags.Test(Params::eDraw_f1_x) ? 1
                     : params.m_DrawFlags.Test(Params::eDraw_f2_x) ? 2
                     : 3 );
            Vec3r X[4] = { element.m_X[0], element.m_X[1], element.m_X[2], element.m_X[3] };
            for( int it_z = params.m_SliceZ0; it_z <= params.m_SliceZ1; it_z++ )
            {
                for( int it_y = params.m_SliceY0; it_y <= params.m_SliceY1; it_y++ )
                {
                    for( int it_x = params.m_SliceX0; it_x <= params.m_SliceX1; it_x++ )
                    {
                        X[nid] = Vec3r( half_sizes[0] * Real(it_x) / hcx,
                                        half_sizes[1] * Real(it_y) / hcy,
                                        half_sizes[2] * Real(it_z) / hcz );
                        Vec12r f( element.ComputeElasticForce( X, (Params::EMethod)params.m_Method ) );
                        Vec3r f_i( f[3*idx], f[3*idx+1], f[3*idx+2] );
                        float c_i( f_i.Norm() / max_norm_f );
                        Vec3r v_i( length * mal::SafeNormalized( f_i ) );
                        //(x,y,z) forces
                        DrawVector( X[nid], v_i,
                                    sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c_i,
                                                                      0.25f,
                                                                      0.25f + 0.75f*(1.0-c_i),
                                                                      1 ),
                                                     1 ) );
                        /*(x,y,c) forces
                          DrawVector( sfr::Vec3(X[nid].x(),X[nid].y(),c_i), sfr::Vec3(v_i.x(),v_i.y(),0),
                          sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c_i,
                          0.25f,
                          0.25f + 0.75f*(1.0-c_i) ),
                          1 ) );
                        */
                    }
                }
            }
        }

#ifdef __TODO_PORT_TO_3D
        if( params.m_DrawFlags.Test(Params::eDraw_P_e) )
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            float max_norm_f( half_sizes.Norm() * element.m_YoungModulus );
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    Vec2r vec_diag_F( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Vec2r vec_diag_P( - element.ComputeDiagonalP( vec_diag_F, (Params::EMethod)params.m_Method ) ); //\todo MUST NEGATE to draw... just convention
                    float c( vec_diag_P.Norm() / max_norm_f );
                    Vec2r v( length * mal::SafeNormalized( vec_diag_P ) );
                    DrawVector( sfr::Vec3(vec_diag_F.x(),vec_diag_F.y(),c), sfr::Vec3(v.x(),v.y(),0),
                                sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c,
                                                                  0.25f,
                                                                  0.25f + 0.75f*(1.0-c) ),
                                                 1 ) );
                }
            }
            DrawSegment2( sfr::Vec2(-half_sizes.x(), half_sizes.y()),
                          sfr::Vec2( half_sizes.x(),-half_sizes.y()),
                          sfr::gfx::Style( sfr::gfx::Color(0.5,0.5,0.5), 2 ) );
        }
#endif

        // Sample Rotation field
        if( params.m_DrawFlags.Test(Params::eDraw_R_x) )
        {
            Vec3r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY, params.m_HalfSizeZ );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            int hcz( params.m_HalfCellsZ );
            float length( 0.8f * mal::Min( half_sizes[0] / hcx,
                                           mal::Min( half_sizes[1] / hcy,
                                                     half_sizes[2] / hcz ) ) );
            Vec3r X[4] = { element.m_X[0], element.m_X[1], element.m_X[2], element.m_X[3] };
            for( int it_z = params.m_SliceZ0; it_z <= params.m_SliceZ1; it_z++ )
            {
                for( int it_y = params.m_SliceY0; it_y <= params.m_SliceY1; it_y++ )
                {
                    for( int it_x = params.m_SliceX0; it_x <= params.m_SliceX1; it_x++ )
                    {
                        X[nid] = Vec3r( half_sizes[0] * Real(it_x) / hcx,
                                        half_sizes[1] * Real(it_y) / hcy,
                                        half_sizes[2] * Real(it_z) / hcz );
                        Mat3x3r R;
                        element.ComputeElasticForce( X, (Params::EMethod)params.m_Method, &R );
                        DrawRefSys( X[nid], R, length, sfr::gfx::Style() );
                    }
                }
            }
        }

#ifdef __TODO_PORT_TO_3D
        // Sample Det(F)
        if( params.m_DrawFlags.Test(Params::eDraw_Det_F) )
        {
            Vec2r half_sizes( params.m_HalfSizeX, params.m_HalfSizeY );
            int hcx( params.m_HalfCellsX );
            int hcy( params.m_HalfCellsY );
            float length( 0.8f * half_sizes[0] / hcx );
            Vec2r X[3] = element.m_X;
            for( int it_y = -hcy; it_y < hcy; it_y++ )
            {
                for( int it_x = -hcx; it_x < hcx; it_x++ )
                {
                    X[nid] = Vec2r( half_sizes[0] * Real(it_x) / hcx, half_sizes[1] * Real(it_y) / hcy );
                    Mat2x2r F = element.Compute_F( X, element.m_FFM );
                    Real det_F( mal::Det(F) );
                    float c( det_F );/// max_det_F );
                    DrawVector( sfr::Vec3(X[nid].x(),X[nid].y(),0), sfr::Vec3(0,0,c),
                                sfr::gfx::Style( sfr::gfx::Color( 0.25f + 0.75f*c,
                                                                  0.25f,
                                                                  0.25f + 0.75f*(1.0-c) ),
                                                 1 ) );
                }
            }
        }
#endif

        // Draw QR-GSO critical point
        {
            Real den( element.m_invDm(0,0) + element.m_invDm(1,0) + element.m_invDm(2,0) );
            if( mal::Abs(den) > mal::Epsilon<Real>() )
            {
                Vec3r gso_critical_X0( (element.m_r[1].x() * element.m_invDm(0,0) + element.m_r[2].x() * element.m_invDm(1,0) + + element.m_r[3].x() * element.m_invDm(2,0) ) / den,
                                       (element.m_r[1].y() * element.m_invDm(0,0) + element.m_r[2].y() * element.m_invDm(1,0) + + element.m_r[3].y() * element.m_invDm(2,0) ) / den,
                                       (element.m_r[1].z() * element.m_invDm(0,0) + element.m_r[2].z() * element.m_invDm(1,0) + + element.m_r[3].z() * element.m_invDm(2,0) ) / den );
                Mat3x3r Ds( mal::GMat3x3_From_Columns( element.m_r[1] - gso_critical_X0, element.m_r[2] - gso_critical_X0, element.m_r[3] - gso_critical_X0 ) );
                Mat3x3r F( Ds * element.m_invDm );
                //TEMP: std::cout << mal::GColumn<0>(F) << std::endl;
                DrawPoint( gso_critical_X0, sfr::gfx::Style(sfr::gfx::Color(1,0,1,1),20) );
            }
        }

        return true;
    }

private:
    AppTask &m_rAppTask;
    Flags32 m_DrawFlags;

private: //hack to get stats to file ASAP
    friend class AppTask;
};

#endif //TEST_FE_APPRENDERER_H
