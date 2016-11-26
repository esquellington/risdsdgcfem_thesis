#ifndef TEST_FE_ELEMENT2D_H
#define TEST_FE_ELEMENT2D_H

#include "Config.h"
#include "Params.h"
#include <Mal/GMatDecomposition.h>

#ifdef __USE_GSL
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_permutation.h>
#  include <gsl/gsl_linalg.h>
#endif

namespace mal
{

template <typename T>
inline void GQRDecomposition( const GMat<T,2,2> &F, GMat<T,2,2> &QR_Q, GMat<T,2,2> &QR_R )
{
#ifdef __USE_GSL_QR
    /*\see http://www.gnu.org/software/gsl/manual/html_node/QR-Decomposition.html
      \note This DOES NOT GUARANTEE a right-handed rotation QR_Q
    */
    // Compute (QR,tau)
    gsl_matrix *gslQR = gsl_matrix_alloc( 2, 2 );
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            gsl_matrix_set( gslQR, i, j, F(i,j) );
    gsl_vector *gslTau = gsl_vector_alloc( 2 );
    int res = gsl_linalg_QR_decomp( gslQR, gslTau );
    if(res) printf ("error: %s\n", gsl_strerror (res));
    // Unpack (QR,tau) into (QR_Q,QR_R)
    gsl_matrix *gslQ = gsl_matrix_alloc( 2, 2 );
    gsl_matrix *gslR = gsl_matrix_alloc( 2, 2 );
    res = gsl_linalg_QR_unpack( gslQR, gslTau, gslQ, gslR );
    if(res) printf ("error: %s\n", gsl_strerror (res));
    // Save Q and R
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            QR_Q(i,j) = gsl_matrix_get( gslQ, i, j );
            QR_R(i,j) = gsl_matrix_get( gslR, i, j );
        }
    gsl_matrix_free( gslR );
    gsl_matrix_free( gslQ );
    gsl_vector_free( gslTau );
    gsl_matrix_free( gslQR );
#elif defined(__USE_GSL_QRPT)
    //\see http://www.gnu.org/software/gsl/manual/html_node/QR-Decomposition-with-Column-Pivoting.html
    // Compute (QR,tau)
    gsl_matrix *gslF = gsl_matrix_alloc( 2, 2 );
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            gsl_matrix_set( gslF, i, j, F(i,j) );
    gsl_matrix *gslQ = gsl_matrix_alloc( 2, 2 );
    gsl_matrix *gslR = gsl_matrix_alloc( 2, 2 );
    gsl_vector *gslTau = gsl_vector_alloc( 2 );
    int gslSignum(0);
    gsl_permutation *gslP = gsl_permutation_alloc( 2 );
    gsl_vector *gslNorm = gsl_vector_alloc( 2 );
    int res = gsl_linalg_QRPT_decomp2( gslF, gslQ, gslR, gslTau, gslP, &gslSignum, gslNorm );
    if(res) printf ("error: %s\n", gsl_strerror (res));
    // Save Q and R
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            QR_Q(i,j) = gsl_matrix_get( gslQ, i, j );
            QR_R(i,j) = gsl_matrix_get( gslR, i, j );
        }
    gsl_vector_free( gslNorm );
    gsl_permutation_free( gslP );
    gsl_vector_free( gslTau );
    gsl_matrix_free( gslR );
    gsl_matrix_free( gslQ );
    gsl_matrix_free( gslF );
#else
    QR_Q = GMat<T,2,2>::Identity();
    QR_R = GMat<T,2,2>::Identity();
#endif
}

template <typename T>
inline GMat<T,2,2> GRotation2x2_QR( const GMat<T,2,2> &F )
{
    GMat<T,2,2> QR_Q,QR_R;
    GQRDecomposition( F, QR_Q, QR_R ); //\todo This DOES NOT GUARANTEE a right-handed rotation QR_Q
    /*TEMP: This makes GSL QR discontinuous, BUT AT LEAST ONE CONTINUOUS PATCH matches my "geometric" QR implementation, based in NPF05...
    Real det_Q( mal::Det(QR_Q) );
    if( det_Q < 0 ) mal::GSetColumn<0>( QR_Q, -mal::GColumn<0>(QR_Q) );
    */
    return QR_Q; //funny, but QR_Q is actually R
    /*
    Real det_Q( mal::Det(QR_Q) );
    if( det_Q > 0 )
        return QR_Q; //funny, but QR_Q is actually R
    else
    {
        //
        //APP_LOG_WARNING( "QR with det %f", det_Q );
        return mal::GRotation2x2_GramSchmidtOrthonormalization_XY( F, Mat2x2r::Identity() );
        //
    }
    */
}

} //namespace mal

#ifdef __USE_NOC
#  include <Mal/GSolvePolynomialEq.h>
   int ComputeNoC( const Vec2r X0[3], const Vec2r X1[3], Real dtdf, Real area );
#else
   inline int ComputeNoC( const Vec2r X0[3], const Vec2r X1[3], Real dtdf, Real area ) { return -1; }
#endif

inline Mat2x2r RLerp( const Mat2x2r &R0, const Mat2x2r &R1, Real lambda01 )
{
    Real w1( lambda01 );
    Real w0( Real(1)-w1 );
    Mat2x2r LR( w0*R0 + w1*R1 );
    Real det_LR( mal::Det(LR) );
    return mal::GRotation2x2_PolarDecomposition( LR, det_LR ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix
    /* TEMP: non-orthonormal LR yields reasonable forcefield and
       det(LR) > 0.65 for all observed scenarios, it MAY be possible
       to reduce CPU (significantly, in 3D) cost by skipping PD(LR)
       with negligible impact in accuracy
    if( !mal::ApproxEq( det_LR, Real(1), Real(0.1) ) ) APP_LOG_WARNING("RCerp rotation severely non-orthonormal, det(LR) = %f", det_LR );
    return LR;
    */
}

inline Mat2x2r RCerp( const Mat2x2r &R0, const Mat2x2r &R1, Real lambda01 )
{
    Real w1( mal::Sq(lambda01) * (3-2*lambda01) ); //h01 polynomial from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
    Real w0( Real(1) - w1 );
    Mat2x2r CR( w0*R0 + w1*R1 );
    Real det_CR( mal::Det(CR) );
    return mal::GRotation2x2_PolarDecomposition( CR, det_CR ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix

    /* TEMP: non-orthonormal CR yields reasonable forcefield and
       det(CR) > 0.65 for all observed scenarios, it MAY be possible
       to reduce CPU (significantly, in 3D) cost by skipping PD(CR)
       with negligible impact in accuracy
    if( !mal::ApproxEq( det_CR, Real(1), Real(0.1) ) ) APP_LOG_WARNING("RCerp rotation severely non-orthonormal, det(CR) = %f", det_CR );
    return CR;
    */
}

class Element2D
{
private:
    const Params &m_rParams;
public:
    Element2D( const Params &params )
    : m_rParams(params)
    , m_LameMu(0)
    , m_LameLambda(0)
    , m_Area(0)
    , m_NoC(0)
        {}
    ~Element2D() { Clear(); }

    void Clear()
        {
        }
    void RebuildFromParams()
        {
            Clear();
            LameParameters_From_YoungAndPoisson( m_rParams.m_YoungModulus, m_rParams.m_PoissonRatio, m_LameMu, m_LameLambda );
            // Init element
            Mat2x2r R_alpha( mal::GRotation2x2_From(m_rParams.m_Alpha) );
            m_r[0] = R_alpha * Real(1) * Vec2r( mal::Cos(m_rParams.m_Angle0), mal::Sin(m_rParams.m_Angle0) );
            Vec2r dr( Real(0.5) * m_rParams.m_Dist12 * Vec2r(0,1) );
            /*TEMP: Old, it's visually easier to inspect if we rotate r0 around origin instead of r1-r2
            m_r[0] = R_alpha*Vec2r(1,0);
            Vec2r dr( Real(0.5) * dist12 * Vec2r( mal::Sin(angle0), mal::Cos(angle0) ) );
            */
            m_r[1] = R_alpha * dr;
            m_r[2] = -R_alpha * dr;
            m_invDm = mal::Inverse( mal::GMat2x2_From_Columns( m_r[1] - m_r[0], m_r[2] - m_r[0] ) );
            // Compute LinearTriangle2 Ke0
            Vec2r r0( m_r[0] );
            Vec2r r1( m_r[1] );
            Vec2r r2( m_r[2] );
            m_Area = ( Real(0.5) *( ( r1.x() - r0.x() ) * ( r2.y() - r0.y() ) - ( r2.x() - r0.x() ) * ( r1.y() - r0.y() ) ) );
            Real rcp_two_area = mal::Rcp( Real(2) * m_Area );
            Real x02( r0.x() - r2.x() ); x02 *= rcp_two_area;
            Real x10( r1.x() - r0.x() ); x10 *= rcp_two_area;
            Real x21( r2.x() - r1.x() ); x21 *= rcp_two_area;
            Real y01( r0.y() - r1.y() ); y01 *= rcp_two_area;
            Real y12( r1.y() - r2.y() ); y12 *= rcp_two_area;
            Real y20( r2.y() - r0.y() ); y20 *= rcp_two_area;
            m_B(0,0) = y12; m_B(0,1) =   0; m_B(0,2) = y20; m_B(0,3) =   0; m_B(0,4) = y01; m_B(0,5) =   0;
            m_B(1,0) =   0; m_B(1,1) = x21; m_B(1,2) =   0; m_B(1,3) = x02; m_B(1,4) =   0; m_B(1,5) = x10;
            m_B(2,0) = x21; m_B(2,1) = y12; m_B(2,2) = x02; m_B(2,3) = y20; m_B(2,4) = x10; m_B(2,5) = y01;
            Mat3x3r E;
            // Plane strain
            Real a( m_rParams.m_YoungModulus / ((1+m_rParams.m_PoissonRatio)*(1-2*m_rParams.m_PoissonRatio)) );
            Real a_times_nu( a * m_rParams.m_PoissonRatio );
            E(0,0) = a - a_times_nu; E(0,1) = a_times_nu;     E(0,2) = Real(0);
            E(1,0) = a_times_nu;     E(1,1) = a - a_times_nu; E(1,2) = Real(0);
            E(2,0) = 0;              E(2,1) = 0;              E(2,2) = Real(0.5)*a - a_times_nu;
            /*TEMP: Plane stress, tal com ho feia el Toni, es
              correspon a canvi de params E i \nu... però crec que és
              incorrecte, NO fer-ho servir! Veure IFEM Ch14, exercici 14.1 per explicació.
            Real E11( m_rParams.m_YoungModulus / (1-mal::Sq(m_rParams.m_PoissonRatio) ) );
            Real E33( 0.5 * m_rParams.m_YoungModulus / (1+m_rParams.m_PoissonRatio) );
            Real E12( m_rParams.m_PoissonRatio * E11 );
            E(0,0) = E11; E(0,1) = E12; E(0,2) = Real(0);
            E(1,0) = E12; E(1,1) = E11; E(1,2) = Real(0);
            E(2,0) = 0;   E(2,1) = 0;   E(2,2) = E33;
            TEMP*/
            m_Ke = m_Area * ( m_B.Transposed() * (E * m_B) );
            m_Ke0 = mal::GRange<0,0,1,5>( m_Ke ); //Force on node0 depending on x0,x1,x2 positions, see LinearTriangle2

            // Variable triangle node pos
            Mat2x2r R_beta( mal::GRotation2x2_From(m_rParams.m_Beta) );
            m_X[0] = R_beta * m_r[0];
            m_X[1] = m_rParams.m_Scale12 * R_beta * m_r[1];
            m_X[2] = m_rParams.m_Scale12 * R_beta * m_r[2];

            m_NoC = m_rParams.m_NoC;
        }

    Vec6r ComputeElasticForce( const Vec2r X[3], Params::EMethod method, Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            switch( method )
            {
            case Params::eMethod_Id: return ComputeElasticForce_Id( X, p_rotation, p_energy ); break;
            case Params::eMethod_QR_XY: return ComputeElasticForce_QR_XY( X, p_rotation, p_energy ); break;
            case Params::eMethod_QR_YX: return ComputeElasticForce_QR_YX( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD: return ComputeElasticForce_PD( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Project: return ComputeElasticForce_PD_Project( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Reflect: return ComputeElasticForce_PD_Reflect( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Unified: return ComputeElasticForce_PD_Unified( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Fix: return ComputeElasticForce_PD_Fix( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Project_Nearest: return ComputeElasticForce_PD_Project_Nearest( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_SVD: return ComputeElasticForce_PD_SVD( X, p_rotation, p_energy ); break;
            case Params::eMethod_IHFSDM: return ComputeElasticForce_IHFSDM( X, p_rotation, p_energy ); break;
            case Params::eMethod_ITF_LRM: return ComputeElasticForce_ITF_LRM( X, p_rotation, p_energy ); break;
            case Params::eMethod_ECIE_CLRM: return ComputeElasticForce_ECIE_CLRM( X, p_rotation, p_energy ); break;
            case Params::eMethod_ECIE_NHC0: return ComputeElasticForce_ECIE_NHC0( X, p_rotation, p_energy ); break;
            case Params::eMethod_ECIE_NHC1: return ComputeElasticForce_ECIE_NHC1( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_CLRM: return ComputeElasticForce_PD_CLRM( X, p_rotation, p_energy ); break;
            case Params::eMethod_DAPD_NS: return ComputeElasticForce_DAPD_NS( X, p_rotation, p_energy ); break;
            case Params::eMethod_DAPD_SS: return ComputeElasticForce_DAPD_SS( X, p_rotation, p_energy ); break;
            case Params::eMethod_DAPD_EX: return ComputeElasticForce_DAPD_EX( X, p_rotation, p_energy ); break;
            default: return Vec6r::Zero();
            }
        }

    Vec2r ComputeDiagonalP( const Vec2r &vec_diag_F, Params::EMethod method ) const
        {
            switch( method )
            {
            case Params::eMethod_Id: return Vec2r::Zero(); break;
            case Params::eMethod_QR_XY: return Vec2r::Zero(); break;
            case Params::eMethod_QR_YX: return Vec2r::Zero(); break;
            case Params::eMethod_PD: return Vec2r::Zero(); break;
            case Params::eMethod_PD_Project: return Vec2r::Zero(); break;
            case Params::eMethod_PD_Reflect: return Vec2r::Zero(); break;
            case Params::eMethod_PD_Unified: return Vec2r::Zero(); break;
            case Params::eMethod_PD_Fix: return Vec2r::Zero(); break;
            case Params::eMethod_PD_Project_Nearest: return Vec2r::Zero(); break;
            case Params::eMethod_PD_SVD: return Vec2r::Zero(); break;
            case Params::eMethod_IHFSDM: return Vec2r::Zero(); break;
            case Params::eMethod_ITF_LRM: return ComputeDiagonalP_ITF_LRM( vec_diag_F ); break;
            case Params::eMethod_ECIE_CLRM: return ComputeDiagonalP_ECIE_CLRM( vec_diag_F ); break;
            case Params::eMethod_ECIE_NHC0: return ComputeDiagonalP_ECIE_NHC0( vec_diag_F ); break;
            case Params::eMethod_ECIE_NHC1: return ComputeDiagonalP_ECIE_NHC1( vec_diag_F ); break;
            case Params::eMethod_PD_CLRM: return Vec2r::Zero(); break;
            case Params::eMethod_DAPD_NS: return Vec2r::Zero(); break;
            case Params::eMethod_DAPD_SS: return Vec2r::Zero(); break;
            case Params::eMethod_DAPD_EX: return Vec2r::Zero(); break;
            default: return Vec2r::Zero();
            }
        }

    Mat2x2r Compute_F( const Vec2r X[3], uint32 ffm ) const
        {
            Mat2x2r Ds( mal::GMat2x2_From_Columns( X[1] - X[0], X[2] - X[0] ) );
            Mat2x2r F( Ds * m_invDm );
            Real det_F( mal::Det(F) );
            if( det_F <= m_rParams.m_DegenerateThresholdDetF )
            {
                switch( ffm )
                {
                case Params::eFFM_Project:
                    {
                        Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                        Real dist12_sq( mal::NormSq(axis12) );
                        if( dist12_sq > 0.000001f )
                        {
                            /* Compute delta so that projected x0 yields an element with det(F) = DTDF
                               det(F) = det(Q) / det(P) = DTDF
                               => det(Q) = DTDF * det(P)
                               ( det(Q) = 2*Area(Q) = delta * dist(x1,x2) )
                               => delta * dist(x1,x2) = DTDF * det(P)
                               => delta = DTDF * det(P) / dist(x1,x2)
                               => delta = DTDF * 2*Area(P) / dist(x1,x2)
                            */

                            // Faster code that AVOIDS sqrt (using unnormalized axis12 and dist12_sq)
                            Real delta_div_dist12( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12_sq ); //== delta / dist12
                            Vec2r projected_x0( X[0] );
                            projected_x0 -= ( mal::Dot(X[0]-X[1],axis12)/dist12_sq ) * axis12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
                            projected_x0 += delta_div_dist12*axis12; //Move past collapse

                            /*\todo Old SLOW code that requires an SQRT
                              Real dist12( mal::Norm( axis12 ) );
                              axis12 /= dist12;
                              Real delta( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                              Vec2r projected_x0( X[0] - mal::Dot(X[0]-X[1],axis12)*axis12 ); //Project onto collapse line
                              projected_x0 += delta*axis12; //Move past collapse
                            */

                            Ds = mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 );
                            F = Ds * m_invDm;
                            det_F = mal::Det(F);

                            if( det_F < 0.99f*m_rParams.m_DegenerateThresholdDetF || det_F > 1.01f*m_rParams.m_DegenerateThresholdDetF )
                                APP_LOG_WARNING( "Compute_F: det(F) %f != %f should not happen if delta is correct", det_F, m_rParams.m_DegenerateThresholdDetF );
                        }
                        else
                        {
                            F = Mat2x2r::Identity();
                            det_F = 1;
                        }
                    }
                    break;
                case Params::eFFM_Project_Nearest:
                    {
                        //\todo USE sqrt-free version as in Params::eFFM_Project!!
                        // Compute all possible axis
                        Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) );
                        Real dist12( mal::Norm( axis12 ) );
                        bool bCollapsed12( dist12 < 0.0001f );

                        Vec2r axis20( mal::PerpendicularCW( X[0] - X[2] ) );
                        Real dist20( mal::Norm( axis20 ) );
                        bool bCollapsed20( dist20 < 0.0001f );

                        Vec2r axis01( mal::PerpendicularCW( X[1] - X[0] ) );
                        Real dist01( mal::Norm( axis01 ) );
                        bool bCollapsed01( dist01 < 0.0001f );
                        if( !bCollapsed12 || !bCollapsed20 || !bCollapsed01 )
                        {
                            // Compute all possible projections, discarding potentially collapsed axis
                            Real delta12(0);
                            if( !bCollapsed12 )
                            {
                                axis12 /= dist12;
                                delta12 = m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12;
                            }
                            Vec2r projected_x0( X[0] - mal::Dot(X[0]-X[1],axis12)*axis12 ); //Project onto collapse line
                            projected_x0 += delta12*axis12; //Move past collapse

                            Real delta20(0);
                            if( !bCollapsed20 )
                            {
                                axis20 /= dist20;
                                delta20 = m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist20;
                            }
                            Vec2r projected_x1( X[1] - mal::Dot(X[1]-X[0],axis20)*axis20 ); //Project onto collapse line
                            projected_x1 += delta20*axis20; //Move past collapse

                            Real delta01(0);
                            if( !bCollapsed01 )
                            {
                                axis01 /= dist01;
                                delta01 = m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist01;
                            }
                            Vec2r projected_x2( X[2] - mal::Dot(X[2]-X[0],axis01)*axis01 ); //Project onto collapse line
                            projected_x2 += delta01*axis01; //Move past collapse

                            // Compute projection distances, to pick nearest valid (uncollapsed axis) one
                            Real d0( mal::NormSq( X[0] - projected_x0 ) );
                            Real d1( mal::NormSq( X[1] - projected_x1 ) );
                            Real d2( mal::NormSq( X[2] - projected_x2 ) );
                            Real max_d( d0 + d1 + d2 );
                            if( bCollapsed12 ) d0 = max_d;
                            if( bCollapsed20 ) d1 = max_d;
                            if( bCollapsed01 ) d2 = max_d;

                            // Choose nearest projection axis
                            Vec2r x0( X[0] );
                            Vec2r x1( X[1] );
                            Vec2r x2( X[2] );
                            if( d0 <= d1 )
                            {
                                if( d0 <= d2 ) x0 = projected_x0;
                                else x2 = projected_x2;
                            }
                            else
                            {
                                if( d1 <= d2 ) x1 = projected_x1;
                                else x2 = projected_x2;
                            }
                            // Compute fixed F
                            Ds = mal::GMat2x2_From_Columns( x1 - x0, x2 - x0 );
                            F = Ds * m_invDm;
                            det_F = mal::Det(F);
                            if( det_F < 0.99f*m_rParams.m_DegenerateThresholdDetF || det_F > 1.01f*m_rParams.m_DegenerateThresholdDetF )
                                APP_LOG_WARNING( "Compute_F: det(F) %f != %f should not happen if delta is correct", det_F, m_rParams.m_DegenerateThresholdDetF );
                        }
                        else
                        {
                            F = Mat2x2r::Identity();
                            det_F = 1;
                        }
                    }
                    break;
                case Params::eFFM_Reflect:
                    {
                        if( det_F < 0 )
                        {
                            Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                            F = mal::GReflection2x2_From_Axis( mal::Normalized( axis12 ) ) * F;
                            det_F = -det_F;
                        }
                    }
                    break;
                case Params::eFFM_None:
                default:
                    break;
                }
            }
            APP_ASSERT( Params::eFFM_None == ffm || mal::Det(F) >= 0 );
            APP_ASSERT( !mal::IsNaN(F) );
            return F;
        }

    Vec6r Compute_C_df_x( const Vec2r X[3], const Vec2r dX[3], const Mat2x2r &R, const Mat2x2r &dR ) const
        {
            Vec2r ldx0( R.Transposed()*dX[0] );
            Vec2r ldx1( R.Transposed()*dX[1] );
            Vec2r ldx2( R.Transposed()*dX[2] );
            Vec6r ldx( ldx0.x(), ldx0.y(), ldx1.x(), ldx1.y(), ldx2.x(), ldx2.y() );
            Mat6x6r Re(Real(0));
            mal::GSetRange<0,0,1,1>( Re, R );
            mal::GSetRange<2,2,3,3>( Re, R );
            mal::GSetRange<4,4,5,5>( Re, R );
            // Add truncated dR main term: df = R_e * K_e * R_e^T * dx
            Vec6r df = -Re * m_Ke * ldx;
            // Term dR_e * K_e * (R_eP^T x_e - r_e) = dR_e * K_e * u_e
            Mat6x6r dRe(Real(0));
            mal::GSetRange<0,0,1,1>( dRe, dR );
            mal::GSetRange<2,2,3,3>( dRe, dR );
            mal::GSetRange<4,4,5,5>( dRe, dR );
            Vec2r u0( R.Transposed()*X[0] - m_r[0] );
            Vec2r u1( R.Transposed()*X[1] - m_r[1] );
            Vec2r u2( R.Transposed()*X[2] - m_r[2] );
            Vec6r u( u0.x(), u0.y(), u1.x(), u1.y(), u2.x(), u2.y() );
            df -= dRe * m_Ke * u;
            // Term R_e * K_e * dR_e^T * x_e
            Vec2r x0( X[0] );
            Vec2r x1( X[1] );
            Vec2r x2( X[2] );
            Vec6r x( x0.x(), x0.y(), x1.x(), x1.y(), x2.x(), x2.y() );
            df -= Re * m_Ke * mal::Transposed(dRe) * x;
            return df;
        }

    Vec6r Compute_H_LRM_df_x( const Vec2r X[3], const Vec2r dX[3], const Mat2x2r &R, const Mat2x2r &dR ) const
        {
            Mat2x2r F( Compute_F( X, Params::eFFM_None ) );
            Real det_F( mal::Det(F) );
            Mat2x2r dDs( mal::GMat2x2_From_Columns( dX[1] - dX[0], dX[2] - dX[0] ) );
            Mat2x2r dF( dDs * m_invDm );
            Mat2x2r Rt( mal::Transposed(R) );
            // Compute dP terms independent of dR, from McAdams2011
            Mat2x2r dP( 2*m_LameMu * dF
                        + m_LameLambda*mal::Trace(Rt*dF) * R );
            //TEMP
            if( mal::Det(dR) == 0 )
            {
                /*
                dP = 2*m_LameMu * dF
                     + m_LameLambda*mal::Trace(Rt*dF) * R;
                */
                /*NAN...
                dP = 2*m_LameMu * R*mal::GMatNxN_Symmetric_Part(Rt*dF)
                     + m_LameLambda*mal::Trace(Rt*dF) * R;
                */
                /* average
                dP = 2*m_LameMu * 0.5*(dF + R*mal::GMatNxN_Symmetric_Part(Rt*dF))
                     + m_LameLambda*mal::Trace(Rt*dF) * R;
                */

                //--- LERP
                const Real cDimension = 2;
                Real weight_WRP(0);

                /* Trace factor, ok
                {
                    weight_WRP = mal::Clamp01( 1 - mal::Abs( 1 - mal::Trace(F)/cDimension ) ); //larger trace => larger largest stretch => smaller weight for WRP
                }
                */

                /* BAAAD
                   {
                   Real det_F_minus_trace( mal::Det( F - (mal::Trace(F)/cDimension)*Mat2x2r::Identity() ) );
                   weight_WRP =  mal::Clamp01( 1 - mal::Abs( 1 - mal::Abs(det_F_minus_trace) ) ) ; //larger trace => larger largest stretch => smaller weight for WRP
                   }
                */

                /* SVD-Largest stretch factor
                {
                    Mat2x2r U, Vt;
                    Vec2r vec_diag_F;
                    mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
                    // Force pure rotation U,Vt by fixing potential inversion/reflection
                    //if( mal::Det(U) * mal::Det(Vt) < 0 )
                    {
                        if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                        if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
                    }
                    //\todo USING m_rParams.m_ECIE_e_threshold as stretch factor weight, or something...
                    weight_WRP = mal::Clamp01( 1 - m_rParams.m_ECIE_e_threshold*mal::Max( mal::Abs(1 - vec_diag_F[0]),
                                                                                          mal::Abs(1 - vec_diag_F[1]) ) ); //Largest POSITIVE stretch deviation from 1...
                }
                */

                /* Approx largest eigenvalue, noisy but +o- correct
                {
                    Mat2x2r S( Rt*F );
                    Real eigen_value;
                    Vec2r eigen_vector;
                    mal::GComputeLargestEigenvalue_Symmetric( S, eigen_value, eigen_vector, 10, 0.01f );
                    weight_WRP = mal::Clamp01( 1 - m_rParams.m_ECIE_e_threshold*mal::Abs(1 - eigen_value) );
                }
                */

                // Try to approximate SVD-Largest with trace
                const Real cRadiusWRP = m_rParams.m_ECIE_k_factor;
                if( cRadiusWRP > 0 )
                {
                    Mat2x2r S( Rt*F );
                    //weight_WRP = mal::Clamp01( 1 - mal::Rcp(cRadiusWRP)*mal::Abs( 1 - mal::Trace(mal::Abs(S))/cDimension ) ); //larger trace => larger largest stretch => smaller weight for WRP
                    //weight_WRP = mal::Clamp01( 1 - mal::Rcp(cRadiusWRP)*mal::Abs( 1 - mal::Trace(S)/cDimension ) ); //larger trace => larger largest stretch => smaller weight for WRP
                    //weight_WRP = mal::Clamp01( 1 - mal::Rcp(cRadiusWRP)*mal::Trace( mal::Abs( S - Mat2x2r::Identity() ) )/cDimension );

                    Real weight_LCM = mal::Clamp01( mal::Rcp(cRadiusWRP) * mal::Trace( mal::Abs( S - Mat2x2r::Identity() ) )/cDimension );
                    //Real weight_LCM = mal::Clamp01( mal::Rcp(cRadiusWRP) * mal::Trace( S - Mat2x2r::Identity() )/cDimension ); Symmetric, not centered at S=I
                    //weight_LCM = mal::Sq( weight_LCM ); //both Sq and Sqrt work better at some dist, but WORSE at stretch=1, which should be the most accurate case
                    weight_WRP = 1 - weight_LCM;
                }
                if( cRadiusWRP > 10 )
                    weight_WRP = 1;
                //

                dP = 2*m_LameMu * ((1-weight_WRP)*dF + weight_WRP*R*mal::GMatNxN_Symmetric_Part(Rt*dF) )
                     + m_LameLambda*mal::Trace(Rt*dF) * R;

                /*TEMP FAST TEST, wrong!!!
                dP = 2*m_LameMu * mal::GMatNxN_Symmetric_Part(dF) + m_LameLambda*mal::Trace(Rt*dF) * R;
                */
            }
            //TEMP

            // Add dP terms that depend on dR, from McAdams2011
            Mat2x2r S( Rt*F );
            dP += ( m_LameLambda*mal::Trace(S-Mat2x2r::Identity()) - 2*m_LameMu ) * dR;
            //dP += mal::GMatNxN_Symmetric_Part( ( m_LameLambda*mal::Trace(S-Mat2x2r::Identity()) - 2*m_LameMu ) * dR ); //TEMP: This seems to make df_dx definite, but still nonsymmetric
            /* TEMP: Warped stiffness, same Pf_Px as R*K*R^T
            Mat2x2r Rt_times_dF( Rt*dF );
            Mat2x2r dP( R * m_LameMu * ( Rt_times_dF + mal::Transposed(Rt_times_dF) ) +
                        R * m_LameLambda * mal::Trace(Rt_times_dF) );
            */
            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm );
            Mat2x2r H( - m_Area * dP * Be0.Transposed() );
            Vec2r df0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r df( df0.x(), df0.y(),
                      H(0,0), H(1,0),
                      H(0,1), H(1,1) );
            return df;
        }

    Vec6r Compute_H_LRM_df_x_SymmStrain( const Vec2r X[3], const Vec2r dX[3], const Mat2x2r &R, const Mat2x2r &dR ) const
        {
            Mat2x2r F( Compute_F( X, Params::eFFM_None ) );
            Mat2x2r dDs( mal::GMat2x2_From_Columns( dX[1] - dX[0], dX[2] - dX[0] ) );
            Mat2x2r dF( dDs * m_invDm );
            Mat2x2r Rt( mal::Transposed(R) );
            Mat2x2r S( Rt*F );
            // Use exact strain tensor where we do not assume S=S^T
            Mat2x2r E( mal::GMatNxN_Symmetric_Part(S) - Mat2x2r::Identity() );
            Mat2x2r dE( mal::GMatNxN_Symmetric_Part(mal::Transposed(dR)*F + Rt*dF) );
            // Compute dP terms independent of dR, from McAdams2011
            Mat2x2r dP( (2*m_LameMu)*(R*dE) + (m_LameLambda*mal::Trace(dE))*R );
            dP += (2*m_LameMu)*(dR*E) + (m_LameLambda*mal::Trace(E))*dR;
            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm );
            Mat2x2r H( - m_Area * dP * Be0.Transposed() );
            Vec2r df0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r df( df0.x(), df0.y(),
                      H(0,0), H(1,0),
                      H(0,1), H(1,1) );
            return df;
        }

    Vec6r Compute_H_CLRM_df_x( const Vec2r X[3], const Vec2r dX[3], const Mat2x2r &R, const Mat2x2r &dR ) const
        {
            Mat2x2r F( Compute_F( X, Params::eFFM_None ) );
            Real det_F( mal::Det(F) );
            Mat2x2r dDs( mal::GMat2x2_From_Columns( dX[1] - dX[0], dX[2] - dX[0] ) );
            Mat2x2r dF( dDs * m_invDm );
            // Compute dP terms independent of dR, from McAdams2011
            Mat2x2r dP( 2*m_LameMu * dF ); //Robust dP term
            // Compute dR-dependent terms
            dP -= 2*m_LameMu * dR;

            // Compute potentially-singular terms
            Real cEpsilonDetF( 10e-6f );
//#define __REGULARIZATION_MODE_ADD
#define __REGULARIZATION_MODE_SPLIT
//#define __REGULARIZATION_MODE_SVD
//#define __REGULARIZATION_MODE_PERTURBATION_dX
//#define __REGULARIZATION_MODE_NUMERICAL_dP
#if defined(__REGULARIZATION_MODE_ADD) //\todo This version is shit, very inaccurate
            // Compute invF-dependent dP terms only if not singular, regularize otherwise
            Real reg_det_F = det_F + mal::Sign(det_F)*cEpsilonDetF; //( mal::Abs(det_F) < cEpsilonDetF ) ? mal::Sign(det_F)*cEpsilonDetF : det_F;
            Mat2x2r invF( mal::Inverse( F, reg_det_F ) ); //\todo THIS FAILS FOR det_F = 0 => Flattened element
            Mat2x2r invFt( mal::Transposed( invF ) );
            dP += m_LameLambda * ( (2*mal::Sq(det_F) - det_F) * mal::Trace(invF*dF) * invFt
                                   - (mal::Sq(det_F) - det_F) * invFt * mal::Transposed(dF) * invFt );
#elif defined(__REGULARIZATION_MODE_SPLIT) //\todo This version is better, less inaccurate, but still problematic
            Mat2x2r atF( mal::Adjugate(F) );
            Mat2x2r atFt( mal::Transposed(atF) );
            Mat2x2r dFt( mal::Transposed(dF) );
            //TEMP: Full term dP += (m_LameLambda / reg_det_F) * ( (2*det_F-1)*mal::Trace(atF*dF)*atFt - (det_F-1)*atFt*dFt*atFt );
            dP += m_LameLambda * ( 2*mal::Trace(atF*dF)*atFt - atFt*dFt*atFt );
            /*TEMP
              Real reg_det_F = (mal::Abs(det_F) < cEpsilonDetF) ? mal::Sign(det_F)*cEpsilonDetF : det_F;
              dP += (m_LameLambda / reg_det_F) * ( atFt*dFt*atFt - mal::Trace(atF*dF)*atFt ); //\todo Truncating this term yields large error in df, but maybe reasonable
            */
            if( mal::Abs(det_F) > cEpsilonDetF )
                dP += (m_LameLambda / det_F) * ( atFt*dFt*atFt - mal::Trace(atF*dF)*atFt ); //\todo Truncating this term yields large error in df, but maybe reasonable
#elif defined(__REGULARIZATION_MODE_SVD) //\todo This is discontinuous...
            // easy term
            Mat2x2r atF( mal::AdjointTranspose(F) );
            Mat2x2r atFt( mal::Transposed(atF) );
            Mat2x2r dFt( mal::Transposed(dF) );
            dP += m_LameLambda * ( 2*mal::Trace(atF*dF)*atFt - atFt*dFt*atFt );
            // problematic term
            Mat2x2r U,Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            const Real cEpsilonSingularValue( 10e-6 );
            Vec2r vec_pseudo_inv_diag_F( ( mal::Abs(vec_diag_F[0]) < cEpsilonSingularValue ) ? cEpsilonSingularValue : mal::Rcp(vec_diag_F[0]),
                                         ( mal::Abs(vec_diag_F[1]) < cEpsilonSingularValue ) ? cEpsilonSingularValue : mal::Rcp(vec_diag_F[1]) );
            Mat2x2r pseudo_inv_diag_F( mal::GMatNxN_From_Diagonal( vec_pseudo_inv_diag_F ) );
            Mat2x2r pseudo_invF( mal::Transposed(Vt) * pseudo_inv_diag_F * mal::Transposed(U) );
            Mat2x2r pseudo_invFt( mal::Transposed( pseudo_invF ) );
            dP += m_LameLambda * ( atFt*dFt - mal::Trace(atF*dF)*Mat2x2r::Identity() ) * pseudo_invFt;

            /*\todo This is LESS ACCURATE than SPLIT version
            dP += m_LameLambda * ( (2*mal::Sq(det_F) - det_F) * mal::Trace(pseudo_invF*dF) * pseudo_invFt
                                   - (mal::Sq(det_F) - det_F) * pseudo_invFt * mal::Transposed(dF) * pseudo_invFt );
            */
#elif defined(__REGULARIZATION_MODE_PERTURBATION_dX) //BAR RESULTS...
            Mat2x2r atF( mal::AdjointTranspose(F) );
            Mat2x2r atFt( mal::Transposed(atF) );
            Mat2x2r dFt( mal::Transposed(dF) );
            dP += m_LameLambda * ( 2*mal::Trace(atF*dF)*atFt - atFt*dFt*atFt );
            if( mal::Abs(det_F) < cEpsilonDetF )
            {
                Mat2x2r perturbationF( mal::GMat2x2_From_Columns( mal::SafeNormalized(dX[1]-dX[0]),
                                                                  mal::SafeNormalized(dX[2]-dX[0]) ) );
                F += cEpsilonDetF*perturbationF;
                det_F = mal::Det(F);
                det_F = (mal::Abs(det_F) < cEpsilonDetF) ? mal::Sign(det_F)*cEpsilonDetF : det_F;
            }
            dP += (m_LameLambda / det_F) * ( atFt*dFt*atFt - mal::Trace(atF*dF)*atFt );
#elif defined(__REGULARIZATION_MODE_NUMERICAL_dP)
            if( mal::Abs(det_F) < cEpsilonDetF )
            {
                return Compute_df_Numerical( X, dX, (Params::EMethod)m_rParams.m_Method, m_rParams.m_DiffR_Numerical_H );
                /*AIXO ESTA MALAMENT...
                Real norm_dX( mal::Sqrt(mal::NormSq(dX[0])+mal::NormSq(dX[1])+mal::NormSq(dX[2])) );
                if( norm_dX == 0 ) return Vec6r::Zero();
                Mat2x2r dF_normalized( mal::Rcp(norm_dX) * dF );
                Mat2x2r leftS( mal::Transposed(R) * (F - m_rParams.m_DiffR_Numerical_H*dF_normalized) );
                Mat2x2r rightS( mal::Transposed(R) * (F + m_rParams.m_DiffR_Numerical_H*dF_normalized) );
                Real det_leftS( mal::Det(leftS) );
                Real det_rightS( mal::Det(rightS) );
                Mat2x2r leftP = R * ( 2*m_LameMu*(leftS-Mat2x2r::Identity())
                                      + m_LameLambda * (det_leftS-1) * Mat2x2r( leftS(1,1), -leftS(1,0),
                                                                                -leftS(0,1), leftS(0,0) ) );
                Mat2x2r rightP = R * ( 2*m_LameMu*(rightS-Mat2x2r::Identity())
                                       + m_LameLambda * (det_rightS-1) * Mat2x2r( rightS(1,1), -rightS(1,0),
                                                                                   -rightS(0,1), rightS(0,0) ) );
                dP = norm_dX * mal::Rcp(2*m_rParams.m_DiffR_Numerical_H) * (rightP - leftP);
                */
            }
            else
            {
                Mat2x2r invF( mal::Inverse( F, det_F ) ); //\todo THIS FAILS FOR det_F = 0 => Flattened element
                Mat2x2r invFt( mal::Transposed( invF ) );
                dP += m_LameLambda * ( (2*mal::Sq(det_F) - det_F) * mal::Trace(invF*dF) * invFt
                                       - (mal::Sq(det_F) - det_F) * invFt * mal::Transposed(dF) * invFt );
            }
#else
            // No regularization at all
            if( mal::Abs(det_F) > 10e-6f )
            {
                Mat2x2r invF( mal::Inverse( F, det_F ) ); //\todo THIS FAILS FOR det_F = 0 => Flattened element
                Mat2x2r invFt( mal::Transposed( invF ) );
                dP += m_LameLambda * ( (2*mal::Sq(det_F) - det_F) * mal::Trace(invF*dF) * invFt
                                       - (mal::Sq(det_F) - det_F) * invFt * mal::Transposed(dF) * invFt );
            }
#endif
            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm );
            Mat2x2r H( - m_Area * dP * Be0.Transposed() );
            Vec2r df0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r df( df0.x(), df0.y(),
                      H(0,0), H(1,0),
                      H(0,1), H(1,1) );
            return df;
        }

    Vec6r Compute_df_Numerical( const Vec2r X[3], const Vec2r dX[3], Params::EMethod method, Real h ) const
        {
            mal::GVec<Real,6> v;
            mal::GSetRange<0,1>( v, dX[0] );
            mal::GSetRange<2,3>( v, dX[1] );
            mal::GSetRange<4,5>( v, dX[2] );
            Real norm_sq_dx( mal::NormSq( v ) );
            if( norm_sq_dx > mal::Epsilon<Real>() )
            {
                Real norm_dx( mal::Sqrt(norm_sq_dx) );
                v /= norm_dx;
                Vec2r v0 = mal::GRange<0,1>( v );
                Vec2r v1 = mal::GRange<2,3>( v );
                Vec2r v2 = mal::GRange<4,5>( v );
                Vec2r leftX[] = { X[0] - h*v0, X[1] - h*v1, X[2] - h*v2 };
                Vec2r rightX[] = { X[0] + h*v0, X[1] + h*v1, X[2] + h*v2 };
                Vec6r left_f = ComputeElasticForce( leftX, method );
                Vec6r right_f = ComputeElasticForce( rightX, method );
                return norm_dx * mal::Rcp(2*h) * (right_f - left_f);
            }
            else
                return Vec6r::Zero();
        }

private:

    Real Compute_Effective_PoissonRatio( Real det_F ) const
        {
            if( det_F >= 0 || m_rParams.m_InvertedCompressibilityFactorDetF < 0.001 ) return m_rParams.m_PoissonRatio;
            else //det_F < 0
            {
                //Real lambda01( mal::Clamp01(1.0 + m_InvertedCompressibilityFactorDetF*det_F) ); //\todo old, not exactly what I wanted...
                Real lambda01( mal::Clamp01( Real(1) + det_F/m_rParams.m_InvertedCompressibilityFactorDetF ) ); //l01 = 1 if det_F >= 0, l01 = 0 if det_F < -ICFDF
                // Linear ipol
                //return m_rParams.m_PoissonRatio * lambda01; // + (1-lambda01)*0
                // Cubic ipol
                return m_rParams.m_PoissonRatio * mal::Sq(lambda01) * (3-2*lambda01); //h01 polynomial from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
            }
        }

    Mat6x6r Compute_KeC() const //\todo Could have params / options
        {
            // Compressible Ke => Increased young, Zero poisson
            Real compressible_young_factor( 1/*+2*m_PoissonRatio*/ ); //increase multiplier for steeper force increase wrt poisson_ratio
            Real a( compressible_young_factor * m_rParams.m_YoungModulus );
            /*Setting young_modulus to 0 in E yields
              Mat3x3r E;
              E(0,0) = a; E(0,1) = 0; E(0,2) = Real(0);
              E(1,0) = 0; E(1,1) = a; E(1,2) = Real(0);
              E(2,0) = 0; E(2,1) = 0; E(2,2) = 0.5*a;
              Mat6x6r KeC( m_Area * ( m_B.Transposed() * (E * m_B) ) );
            */
            //Assuming E(2,2) == a instead of 0.5*a actually yields better force directions, less inverted-incompressible (=> E = a*Id3)
            Mat6x6r KeC( (m_Area * a) * ( m_B.Transposed() * m_B) );
            APP_ASSERT( !mal::IsNaN(KeC) );
            return KeC;
        }

    Mat6x6r Compute_Effective_Ke( Real det_F ) const
        {
            if( det_F >= 0 ) return m_Ke;
            else //det_F < 0
            {
                if( m_rParams.m_InvertedCompressibilityFactorDetF > 0.001 )
                {
                    Real poisson_ratio( Compute_Effective_PoissonRatio(det_F) );
                    Real a( m_rParams.m_YoungModulus / ((1+poisson_ratio)*(1-2*poisson_ratio)) );
                    Real a_times_nu( a * poisson_ratio );
                    Mat3x3r E;
                    E(0,0) = a - a_times_nu; E(0,1) = a_times_nu;     E(0,2) = Real(0);
                    E(1,0) = a_times_nu;     E(1,1) = a - a_times_nu; E(1,2) = Real(0);
                    E(2,0) = 0;              E(2,1) = 0;              E(2,2) = Real(0.5)*a - a_times_nu;
                    return m_Area * ( m_B.Transposed() * (E * m_B) );
                }
                else
                    return (Real(1)-m_rParams.m_FactorDetF*det_F)*m_Ke;
            }
        }

    Vec6r Compute_Corotated_Forces( const Vec2r X[3], const Mat2x2r &R, Real det_F, Real *p_energy = 0 ) const
        {
            Vec2r u0( R.Transposed()*X[0] - m_r[0] );
            Vec2r u1( R.Transposed()*X[1] - m_r[1] );
            Vec2r u2( R.Transposed()*X[2] - m_r[2] );
            Vec6r u( u0.x(), u0.y(), u1.x(), u1.y(), u2.x(), u2.y() );
            Mat6x6r Re(Real(0));
            mal::GSetRange<0,0,1,1>( Re, R );
            mal::GSetRange<2,2,3,3>( Re, R );
            mal::GSetRange<4,4,5,5>( Re, R );
            Mat6x6r Ke( Compute_Effective_Ke(det_F) );
            if( p_energy ) *p_energy = Real(0.5)*(u * (Ke * u));
            return -Re * Ke * u;
        }

    Vec6r ComputeElasticForce_LRM_PK( bool bSymmetrizeS, const Mat2x2r &F, const Mat2x2r &R, Real *p_energy = 0 ) const
        {
            Mat2x2r S( R.Transposed() * F ); //F = R*S
            Mat2x2r E,P;
            if( bSymmetrizeS )
            {
                E = mal::GMatNxN_Symmetric_Part(S) - Mat2x2r::Identity();
                //P = (2*m_LameMu)*(R*E) + (m_LameLambda * mal::Trace(E))*R; //TEMP: Same formula as !bSymmetrizeS, WRONG?!
                P = (2*m_LameMu)*(R*E) + (m_LameLambda * mal::Trace(E))*R;
            }
            else
            {
                E = S - Mat2x2r::Identity();
                P = (2*m_LameMu)*(R*E) + (m_LameLambda * mal::Trace(E))*R;
            }
            if( p_energy ) *p_energy = m_LameMu * (E).NormSqF() + Real(0.5)*m_LameLambda*mal::Sq(mal::Trace(E));
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy
            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm );
            Mat2x2r H( - m_Area * P * Be0.Transposed() ); // Negation is present in the Siggraph2012 course notes, but not in ITF...
            Vec2r f0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r f( f0.x(), f0.y(),
                     H(0,0), H(1,0),
                     H(0,1), H(1,1) );
            return f;
        }

    Vec6r ComputeElasticForce_LRM_PK_Numerical( const Vec2r X[3], const Mat2x2r &F, const Mat2x2r &R, Real *p_energy = 0 ) const
        {
            //P = dPsi/dF = dPsi/dF_ij
            //const Real h(1e-7);
            const Real h(m_rParams.m_DiffR_Numerical_H);
            Mat2x2r P;
            for( int i=0; i<2; i++ )
            {
                for( int j=0; j<2; j++ )
                {
                    // Left
                    Mat2x2r leftF( F ); leftF(i,j) -= h;
                    // Mat2x2r leftR( mal::GRotation2x2_PolarDecomposition_From_SVD( leftF, mal::Det(leftF) ) );
                    Vec2r leftX[3] = { X[0], X[1], X[2] };
                    Mat2x2r leftDeltaDs( (leftF-F)*mal::Inverse(m_invDm) );
                    leftX[0] += -Real(0.3333333f)*(mal::GColumn<0>(leftDeltaDs) + mal::GColumn<1>(leftDeltaDs));
                    leftX[1] += Real(0.6666666f)*mal::GColumn<0>(leftDeltaDs) - Real(0.3333333f)*mal::GColumn<1>(leftDeltaDs);
                    leftX[2] += Real(0.6666666f)*mal::GColumn<1>(leftDeltaDs) - Real(0.3333333f)*mal::GColumn<0>(leftDeltaDs);
                    Mat2x2r leftR( ComputeRotation_DAPD( leftX, leftF, mal::Det(leftF) ) );
                    Mat2x2r leftS( leftR.Transposed() * leftF );
                    Mat2x2r leftE( mal::GMatNxN_Symmetric_Part(leftS) - Mat2x2r::Identity() );
                    Real leftPsi( m_LameMu * (leftE).NormSqF() + Real(0.5)*m_LameLambda*mal::Sq(mal::Trace(leftE)) );

                    // Right
                    Mat2x2r rightF( F ); rightF(i,j) += h;
                    // Mat2x2r rightR( mal::GRotation2x2_PolarDecomposition_From_SVD( rightF, mal::Det(rightF) ) );
                    Vec2r rightX[3] = { X[0], X[1], X[2] };
                    Mat2x2r rightDeltaDs( (rightF-F)*mal::Inverse(m_invDm) );
                    rightX[0] += -Real(0.3333333f)*(mal::GColumn<0>(rightDeltaDs) + mal::GColumn<1>(rightDeltaDs));
                    rightX[1] += Real(0.6666666f)*mal::GColumn<0>(rightDeltaDs) - Real(0.3333333f)*mal::GColumn<1>(rightDeltaDs);
                    rightX[2] += Real(0.6666666f)*mal::GColumn<1>(rightDeltaDs) - Real(0.3333333f)*mal::GColumn<0>(rightDeltaDs);
                    Mat2x2r rightR( ComputeRotation_DAPD( rightX, rightF, mal::Det(rightF) ) );
                    Mat2x2r rightS( rightR.Transposed() * rightF );
                    Mat2x2r rightE( mal::GMatNxN_Symmetric_Part(rightS) - Mat2x2r::Identity() );
                    Real rightPsi( m_LameMu * (rightE).NormSqF() + Real(0.5)*m_LameLambda*mal::Sq(mal::Trace(rightE)) );

                    // dPsi/dF
                    P(i,j) = mal::Rcp(2*h) * (rightPsi - leftPsi);
                }
            }
            Mat2x2r S( R.Transposed() * F ); //F = R*S
            Mat2x2r E( mal::GMatNxN_Symmetric_Part(S) - Mat2x2r::Identity() );
            if( p_energy ) *p_energy = m_LameMu * (E).NormSqF() + Real(0.5)*m_LameLambda*mal::Sq(mal::Trace(E));
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy
            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm );
            Mat2x2r H( - m_Area * P * Be0.Transposed() ); // Negation is present in the Siggraph2012 course notes, but not in ITF...
            Vec2r f0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r f( f0.x(), f0.y(),
                     H(0,0), H(1,0),
                     H(0,1), H(1,1) );
            return f;
        }

    Vec6r Compute_PiolaKirchhoff_Forces( const Mat2x2r &diag_P, const Mat2x2r &U, const Mat2x2r &Vt ) const
        {
            Mat2x2r Be0( m_invDm ); //\todo Be0 = D_m^-1 = invP from the undeformed P, deformed Q notation... confusing
            Mat2x2r H( - m_Area * U * diag_P * Vt * Be0.Transposed() ); // Negation is present in the Siggraph2012 course notes, but not in ITF...
            Vec2r f0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r f( f0.x(), f0.y(),
                     H(0,0), H(1,0),
                     H(0,1), H(1,1) );
            return f;
        }

    Vec6r ComputeElasticForce_Id( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            if( p_rotation ) *p_rotation = Mat2x2r::Identity();
            return Compute_Corotated_Forces( X, Mat2x2r::Identity(), Real(1), p_energy );
        }

    Vec6r ComputeElasticForce_QR_XY( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            /*TEMP Transposed test to match NPF05 usage on J=F^T, though PO09 uses it directly on F... but the results are not the same as transposing R afterwards!!
            Mat2x2r R( mal::GRotation2x2_GramSchmidtOrthonormalization_XY( mal::Transposed(F), Mat2x2r::Identity() ) );
            R = mal::Transposed(R);
            */
            Mat2x2r R( mal::GRotation2x2_GramSchmidtOrthonormalization_XY( F, Mat2x2r::Identity() ) );
            //TEMP:Mat2x2r R( mal::GRotation2x2_GramSchmidtOrthonormalization_XY_2( F, Mat2x2r::Identity() ) ); //\todo This does NOT guarantee right-handed R
            APP_ASSERT( !mal::IsNaN(R) );
            APP_ASSERT( mal::Det(R) > 0 );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }
    Vec6r ComputeElasticForce_QR_YX( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            //TEMP Mat2x2r R( mal::GRotation2x2_GramSchmidtOrthonormalization_YX( F, Mat2x2r::Identity() ) );
            Mat2x2r R( mal::GRotation2x2_QR( F ) );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    Vec6r ComputeElasticForce_PD( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;
            if( mal::Abs(det_F) > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, det_F ); //\todo This has a glitch at det == 0
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec6r ComputeElasticForce_PD_Reflect( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;

            /* Reflect R is biased wrt Alpha, Reflect F is not...
            if( mal::Abs(det_F) > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, det_F ); //\todo This has a glitch at det == 0
            //if( det_F < 0 ) R = mal::GReflection2x2_From_Axis( mal::Normalized( X[1] - X[2] ) ) * mal::GReflection2x2_From_Axis( mal::Normalized( mal::PerpendicularCW( X[1] - X[2] ) ) ) * R * Mat2x2r(1,0,0,-1);
            if( det_F < 0 ) R = R * mal::GReflection2x2_From_Axis( mal::Normalized( mal::PerpendicularCW( X[1] - X[2] ) ) ); //\todo Biased wrt alpha!
            */

            // Unbiased Reflect F
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else if( det_F >= 0 ) R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, det_F ); //\todo This has a glitch at det == 0
            else // det_F < 0
            {
#ifdef __USE_NOC
                APP_ASSERT( m_NoC >= 0 );
                const Vec2r &x0( X[ m_NoC ] ); //Inverted vtx
                const Vec2r &x1( X[ (m_NoC + 1) % 3 ] );
                const Vec2r &x2( X[ (m_NoC + 2) % 3 ] );
                Vec2r axis12( mal::PerpendicularCW( x2 - x1 ));
                Real dist12( mal::Norm(axis12) );
                if( dist12 > 0.001f )  //\todo use sensible epsilon
                {
                    axis12 /= dist12;
                    if( det_F <= -m_rParams.m_DegenerateThresholdDetF ) //Inverted & !Collapsed
                        R = mal::GRotation2x2_PolarDecomposition( mal::GReflection2x2_From_Axis( axis12 ) * F, -det_F ); //Safe uncollapsed simple PD
                    else if( det_F < 0 ) //-m_DegenerateThresholdDetF < det_F < 0 => Inverted & Collapsed
                        R = mal::GRotation2x2_PolarDecomposition_From_SVD( mal::GReflection2x2_From_Axis( axis12 ) * F, -det_F ); //\todo SHOULD Handle det_F = 0 properly
                    APP_ASSERT( !mal::IsNaN( R ) );
                }
#else
                F = mal::GReflection2x2_From_Axis( mal::Normalized( mal::PerpendicularCW( X[1] - X[2] ) ) ) * F;
                //TEMP: Wrong: F = F * mal::GReflection2x2_From_Axis( mal::Normalized( mal::PerpendicularCW( X[1] - X[2] ) ) );
                if( det_F < -m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, -det_F );
                else R = mal::GRotation2x2_PolarDecomposition_From_SVD( F, -det_F );
#endif
            }

            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec6r ComputeElasticForce_PD_Unified( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
                APP_ASSERT( m_NoC >= 0 );
                Vec2r tmpX[3] = { X[0], X[1], X[2] };
                Vec2r &x0( tmpX[ m_NoC ] ); //Inverted vtx
                const Vec2r &x1( tmpX[ (m_NoC + 1) % 3 ] );
                const Vec2r &x2( tmpX[ (m_NoC + 2) % 3 ] );
                Vec2r x12( x2 - x1 );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 0.000001f )
                {
                    Real dist12( mal::Sqrt(dist12_sq) );
                    Vec2r axis12( mal::PerpendicularCW( x12 ) / dist12 );
                    Real delta_h_alpha( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                    Vec2r x12_mid( 0.5 * (x1 + x2) ); //unnecessary, added to ensure dR_numerical symmetry wrt x1,x2,dx1,dx2
                    //Vec2r x12_mid( x1 ); //unnecessary, added to ensure dR_numerical symmetry wrt x1,x2,dx1,dx2
                    Real w( mal::Dot(x0-x12_mid,axis12) ); //\note REQUIRES normalized axis
                    Real s( w / delta_h_alpha );
                    Real t( Real(1) - s );
                    // Cubic spline from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
                    Real lambda_u = ( w >= 0 ) //We know that w <= \Delta h as det_F <= m_rParams.m_DegenerateThresholdDetF
                                    ? m_rParams.m_Unified_L_Factor*delta_h_alpha*( (2-t)*t*t ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                                    : m_rParams.m_Unified_L_Factor*delta_h_alpha*( t + m_rParams.m_Unified_NL_Factor*mal::Pow( -s, m_rParams.m_Unified_NL_Exponent ) ); //lambda_u^L + lambda_u^NL
                    x0 += lambda_u * axis12; //Move past collapse area, to ensure numerically safe

                    Mat2x2r Ds( mal::GMat2x2_From_Columns( tmpX[1] - tmpX[0], tmpX[2] - tmpX[0] ) );
                    Mat2x2r projF( Ds * m_invDm );
                    Real det_projF( mal::Det(projF) );
                    APP_ASSERT( det_projF >= 0 );
                    R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                }
                else
                    R = Mat2x2r::Identity();
            }
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Mat2x2r ComputeRotation_DAPD( const Vec2r X[3], const Mat2x2r& F, Real det_F ) const
        {
            Mat2x2r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else if( m_NoC < 0 ) R = Mat2x2r::Identity(); //TEMPORAL
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
                APP_ASSERT( m_NoC >= 0 );
                Vec2r tmpX[3] = { X[0], X[1], X[2] };
                Vec2r &x0( tmpX[ m_NoC ] ); //Inverted vtx
                const Vec2r &x1( tmpX[ (m_NoC + 1) % 3 ] );
                const Vec2r &x2( tmpX[ (m_NoC + 2) % 3 ] );
                Vec2r x12( x2 - x1 );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 0.000001f )
                {
                    Real dist12( mal::Sqrt(dist12_sq) );
                    Vec2r axis12( mal::PerpendicularCW( x12 ) / dist12 );
                    Real delta_h_alpha( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                    Vec2r x12_mid( 0.5 * (x1 + x2) ); //unnecessary, added to ensure dR_numerical symmetry wrt x1,x2,dx1,dx2
                    //Vec2r x12_mid( x1 ); //unnecessary, added to ensure dR_numerical symmetry wrt x1,x2,dx1,dx2
                    Real w( mal::Dot(x0-x12_mid,axis12) ); //\note REQUIRES normalized axis
                    Real s( w / delta_h_alpha );
                    Real t( Real(1) - s );
                    // Cubic spline from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
                    Real lambda_u = ( w >= 0 ) //We know that w <= \Delta h as det_F <= m_rParams.m_DegenerateThresholdDetF
                                    ? m_rParams.m_Unified_L_Factor*delta_h_alpha*( (2-t)*t*t ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                                    : m_rParams.m_Unified_L_Factor*delta_h_alpha*( t + m_rParams.m_Unified_NL_Factor*mal::Pow( -s, m_rParams.m_Unified_NL_Exponent ) ); //lambda_u^L + lambda_u^NL
                    x0 += lambda_u * axis12; //Move past collapse area, to ensure numerically safe

                    Mat2x2r Ds( mal::GMat2x2_From_Columns( tmpX[1] - tmpX[0], tmpX[2] - tmpX[0] ) );
                    Mat2x2r projF( Ds * m_invDm );
                    Real det_projF( mal::Det(projF) );
                    APP_ASSERT( det_projF >= 0 );
                    R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                }
                else
                    R = Mat2x2r::Identity();
            }
            return R;
        }

    Vec6r ComputeElasticForce_DAPD_NS( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R( ComputeRotation_DAPD( X, F, det_F ) );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return ComputeElasticForce_LRM_PK( false, F, R, p_energy );
        }

    Vec6r ComputeElasticForce_DAPD_SS( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R( ComputeRotation_DAPD( X, F, det_F ) );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return ComputeElasticForce_LRM_PK( true, F, R, p_energy );
        }

    Vec6r ComputeElasticForce_DAPD_EX( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R( ComputeRotation_DAPD( X, F, det_F ) );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return ComputeElasticForce_LRM_PK_Numerical( X, F, R, p_energy );
            // return ComputeElasticForce_DAPD_SS(X) - ComputeElasticForce_LRM_PK_Numerical( X, F, R, p_energy ); //TEMP: SS = EX + ERR => ERR = SS-EX (EX = conservative,divergence-free), but ERR MAY NOT be curl-only!!
        }

    Vec6r ComputeElasticForce_PD_Fix( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
#ifdef __USE_NOC
                APP_ASSERT( m_NoC >= 0 );
                Vec2r tmpX[3] = { X[0], X[1], X[2] };
                Vec2r &x0( tmpX[ m_NoC ] ); //Inverted vtx
                const Vec2r &x1( tmpX[ (m_NoC + 1) % 3 ] );
                const Vec2r &x2( tmpX[ (m_NoC + 2) % 3 ] );
                const Vec2r &r1( m_r[ (m_NoC + 1) % 3 ] );
                const Vec2r &r2( m_r[ (m_NoC + 2) % 3 ] );
                Vec2r x12( x2 - x1 );
                Vec2r r12( r2 - r1 );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 0.000001f )
                {
#define __ENABLE_DERIVABLE_FIX
#ifdef __ENABLE_DERIVABLE_FIX
                    /* Based in the continuity and derivability
                     analysis in DCLFEM revision report, it seems that
                     a projection length \lambda_{project} that grows
                     smoothly and derivably from \lambda = 0 when w=0
                     to some \lambda = -w + delta_h_derivable when -w
                     >= a certain maximum degeneration depth, can be
                     obtained using a cubic hermite spline that
                     interpolates the target delta_h_derivable with 0
                     derivative at both transitions on w.

                     \todo This is WIP and may be wrong, should derive
                     it on paper to ensure factors are correct, but
                     the approach is sound, I'm sure
                    */
                    Real dist12( mal::Sqrt(dist12_sq) );
                    Vec2r axis12( mal::PerpendicularCW( x12 ) / dist12 );
                    Real delta_h_alpha( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                    Real delta_h_1( 2*m_Area / dist12 );
                    Real w( mal::Dot(x0-x1,axis12) );
                    Real s01( mal::Clamp01( (delta_h_alpha - w) / (delta_h_1 - delta_h_alpha) ) ); //Degeneration fraction [0,1] between det(F)=alpha and det(F)=1
                    Real t01( mal::Sq(s01) * (3-2*s01) ); //Degeneration ipol factor between delta_h_alpha and delta_h_1
                    Real delta_h_derivable( (delta_h_alpha - w) + t01 * (delta_h_1 - delta_h_alpha) ); //h01 polynomial from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
                    Real lambda_proj( -w + delta_h_alpha ); //This is C^0 projection as in PD_Project
                    Real lambda_proj_derivable( -w + delta_h_derivable ); //This is C^2 projection (because we use a cubic polynomial)
                    x0 += lambda_proj_derivable * axis12; //Move past collapse area, to ensure numerically safe
                    Mat2x2r Ds( mal::GMat2x2_From_Columns( tmpX[1] - tmpX[0], tmpX[2] - tmpX[0] ) );
                    Mat2x2r projF( Ds * m_invDm );
                    Real det_projF( mal::Det(projF) );
                    APP_ASSERT( det_projF >= 0 );
                    R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
#else
                    // Compute rotation from r12 to x12, which is the same as global R when we rotate the resting state R=Id rigidly with x12
                    //\todo In 3D, rotR would be similarly computed from the NORMAL of the crossed face
                    Mat2x2r rotR( mal::GRotation2x2_Vec2Vec( mal::Normalized(r12) , mal::Normalized(x12) ) ); //\todo Normalized(r12) could be precomputed
                    if( det_F > m_rParams.m_ThresholdIpolDetF )
                    {
                        Vec2r axis12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                        // Faster code that AVOIDS sqrt (using unnormalized axis12 and dist12_sq)
                        Real delta_div_dist12( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12_sq ); //== delta / dist12
                        x0 -= ( mal::Dot(x0-x1,axis12) / dist12_sq ) * axis12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
                        x0 += delta_div_dist12*axis12; //Move past collapse
                        // Compute Proj(F)
                        Mat2x2r projDs( mal::GMat2x2_From_Columns( tmpX[1] - tmpX[0], tmpX[2] - tmpX[0] ) );
                        Mat2x2r projF( projDs * m_invDm );
                        Real det_projF( mal::Det(projF) );
                        APP_ASSERT( det_projF >= 0 );
                        if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                        Mat2x2r projR = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                        Real lambda01( (det_F - m_rParams.m_ThresholdIpolDetF) / (m_rParams.m_DegenerateThresholdDetF - m_rParams.m_ThresholdIpolDetF) );
                        R = RCerp( rotR, projR, lambda01 );
                    }
                    else //det_F <= 0
                        R = rotR;
#endif //__ENABLE_DERIVABLE_FIX
                }
                else
                    R = Mat2x2r::Identity();
#else
                Vec2r x12( X[2] - X[1] );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 0.000001f )
                {
                    // Compute rotation from r12 to x12, which is the same as global R when we rotate the resting state R=Id rigidly with x12
                    Vec2r r12( m_r[2]-m_r[1] );
                    //\todo In 3D, rotR would be similarly computed from the NORMAL of the crossed face
                    Mat2x2r rotR( mal::GRotation2x2_Vec2Vec( mal::Normalized(r12) , mal::Normalized(x12) ) ); //\todo Normalized(r12) could be precomputed
                    if( det_F > m_rParams.m_ThresholdIpolDetF )
                    {
                        Vec2r axis12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                        // Faster code that AVOIDS sqrt (using unnormalized axis12 and dist12_sq)
                        Real delta_div_dist12( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12_sq ); //== delta / dist12
                        Vec2r projected_x0( X[0] );
                        projected_x0 -= ( mal::Dot(X[0]-X[1],axis12)/dist12_sq ) * axis12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
                        projected_x0 += delta_div_dist12*axis12; //Move past collapse
                        // Compute Proj(F)
                        Mat2x2r projDs( mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 ) );
                        Mat2x2r projF( projDs * m_invDm );
                        Real det_projF( mal::Det(projF) );
                        APP_ASSERT( det_projF >= 0 );
                        if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                        Mat2x2r projR = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                        Real lambda01( (det_F - m_rParams.m_ThresholdIpolDetF) / (m_rParams.m_DegenerateThresholdDetF - m_rParams.m_ThresholdIpolDetF) );
                        R = RCerp( rotR, projR, lambda01 );
                    }
                    else //det_F <= 0
                        R = rotR;
                }
                else
                    R = Mat2x2r::Identity();
#endif
            }

            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec6r ComputeElasticForce_PD_Project( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
#ifdef __USE_NOC
                APP_ASSERT( m_NoC >= 0 );
                Vec2r tmpX[3] = { X[0], X[1], X[2] };
                Vec2r &x0( tmpX[ m_NoC ] ); //Inverted vtx
                const Vec2r &x1( tmpX[ (m_NoC + 1) % 3 ] );
                const Vec2r &x2( tmpX[ (m_NoC + 2) % 3 ] );
                Vec2r axis12( mal::PerpendicularCW( x2 - x1 ));
                Real dist12( mal::Norm(axis12) );
                if( dist12 > 0.001f )  //\todo use sensible epsilon
                {
                    axis12 /= dist12;
                    Real delta( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                    x0 -= mal::Dot(x0-x1,axis12) * axis12; //Project on-to axis
                    x0 += delta * axis12; //Move past collapse area, to ensure numerically safe
                    Mat2x2r Ds( mal::GMat2x2_From_Columns( tmpX[1] - tmpX[0], tmpX[2] - tmpX[0] ) );
                    Mat2x2r projF( Ds * m_invDm );
                    Real det_projF( mal::Det(projF) );
                    APP_ASSERT( det_projF >= 0 );
                    /* \todo THIS HAPPENED but now it should be fixed
                       if( det_projF < 0 )
                       {
                       APP_LOG_WARNING( "det(projF) %f < 0 , n12 (%f,%f) , (X[0]-x1)*n12 %f",
                       det_projF, axis12[0], axis12[1], mal::Dot(X[0]-X[1],axis12) );
                       }
                    */
                    if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                    R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                }
                else
                    R = Mat2x2r::Identity();
#else //__USE_NOC
                Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                Real dist12( mal::Norm( axis12 ) );
                if( dist12 > 0.0001f )
                {
                    axis12 /= dist12;
                    Real delta( m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                    Vec2r projected_x0( X[0] - mal::Dot(X[0]-X[1],axis12)*axis12 ); //Project onto collapse line
                    projected_x0 += delta*axis12; //Move past collapse, ensuring det_projF >= m_rParams.m_DegenerateThresholdDetF

                    Mat2x2r Ds( mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 ) );
                    Mat2x2r projF( Ds * m_invDm );
                    Real det_projF( mal::Det(projF) );
                    APP_ASSERT( det_projF >= 0 );
                    /* \todo THIS HAPPENED but now it should be fixed
                       if( det_projF < 0 )
                       {
                       APP_LOG_WARNING( "det(projF) %f < 0 , n12 (%f,%f) , (X[0]-x1)*n12 %f",
                       det_projF, axis12[0], axis12[1], mal::Dot(X[0]-X[1],axis12) );
                       }
                    */
                    if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                    R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                }
                else
                    R = Mat2x2r::Identity();
#endif
            }
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec6r ComputeElasticForce_PD_Project_Nearest( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
                Vec2r x0( X[0] /* * 10 */ ); //\todo Trying to scale inverted axis to collapse circular orbits into x1,x2 line, but it's a BAD IDEA
                // Compute all possible projections, discarding potentially collapsed axis
                Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) );
                Real dist12( mal::Norm( axis12 ) );
                Real delta12(0);
                bool bCollapsed12( dist12 < 0.0001f );
                if( !bCollapsed12 )
                {
                    axis12 /= dist12;
                    delta12 = m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist12;
                }
                Vec2r projected_x0( x0 - mal::Dot(x0-X[1],axis12)*axis12 ); //Project onto collapse line
                projected_x0 += delta12*axis12; //Move past collapse

                Vec2r axis20( mal::PerpendicularCW( x0 - X[2] ) );
                Real dist20( mal::Norm( axis20 ) );
                Real delta20(0);
                bool bCollapsed20( dist20 < 0.0001f );
                if( !bCollapsed20 )
                {
                    axis20 /= dist20;
                    delta20 = m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist20;
                }
                Vec2r projected_x1( X[1] - mal::Dot(X[1]-x0,axis20)*axis20 ); //Project onto collapse line
                projected_x1 += delta20*axis20; //Move past collapse

                Vec2r axis01( mal::PerpendicularCW( X[1] - x0 ) );
                Real dist01( mal::Norm( axis01 ) );
                Real delta01(0);
                bool bCollapsed01( dist01 < 0.0001f );
                if( !bCollapsed01 )
                {
                    axis01 /= dist01;
                    delta01 = m_rParams.m_DegenerateThresholdDetF * 2*m_Area / dist01;
                }
                Vec2r projected_x2( X[2] - mal::Dot(X[2]-x0,axis01)*axis01 ); //Project onto collapse line
                projected_x2 += delta01*axis01; //Move past collapse

                // Compute projection distances, to pick nearest valid (uncollapsed axis) one
                Real d0( mal::NormSq(   x0 - projected_x0 ) );
                Real d1( mal::NormSq( X[1] - projected_x1 ) );
                Real d2( mal::NormSq( X[2] - projected_x2 ) );
                Real max_d( d0 + d1 + d2 );
                if( bCollapsed12 ) d0 = max_d;
                if( bCollapsed20 ) d1 = max_d;
                if( bCollapsed01 ) d2 = max_d;

                Mat2x2r Q;
                if( d0 <= d1 )
                {
                    if( d0 <= d2 ) Q = mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 );
                    else Q = mal::GMat2x2_From_Columns( X[1] - x0, projected_x2 - x0 );
                }
                else
                {
                    if( d1 <= d2 ) Q = mal::GMat2x2_From_Columns( projected_x1 - x0, X[2] - x0 );
                    else Q = mal::GMat2x2_From_Columns( X[1] - x0, projected_x2 - x0 );
                }

                Mat2x2r projF( Q * m_invDm );
                Real det_projF( mal::Det(projF) );
                APP_ASSERT( det_projF >= 0 );
                if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF )
                    APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                R = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
            }
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec6r ComputeElasticForce_PD_SVD( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r diag_S;
            mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
//#define __TESTING_SVD_CLAMP
#ifdef __TESTING_SVD_CLAMP
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                /*
                Mat2x2r barS( diag_S[0], 0,
                              0, m_rParams.m_DegenerateThresholdDetF );
                */
                /*
                Real l0 = diag_S[0];
                Real l1 = diag_S[1];
                Real bl1 = m_rParams.m_DegenerateThresholdDetF;
                Real c( bl1 / l1 );
                Real s( mal::Sqrt( Real(1) - mal::Sq(c) ) );
                Mat2x2r barS( c*l0, s*bl1,
                              -s*l0, bl1 );
                */

                Real l0 = diag_S[0];
                Real l1 = diag_S[1];
                Real bl1 = m_rParams.m_DegenerateThresholdDetF;
                Real divisor( l0+l1 );
                if( mal::Abs(divisor) > 0.001 )
                {
                    Real c_sq( (bl1 + l0) / divisor );
                    Real c( mal::Sqrt( c_sq ) );
                    Real s_sq( Real(1) - c_sq );
                    Real s( mal::Sqrt( s_sq ) );
                    Mat2x2r barS( c_sq*l0 + s_sq*l1, c*s*l0 - c*s*l1,
                                  c*s*l0 - c*s*l1, bl1 );
                    F = U * barS * Vt;
                }
            }
            Mat2x2r R( mal::GRotation2x2_PolarDecomposition_From_SVD( F, mal::Det(F) ) );
            /*
            Mat2x2r R( F
                       * mal::Transposed(Vt)
                       * Mat2x2r( mal::Rcp(diag_S[0]), 0,
                                  0, mal::Rcp( diag_S[1] ) )
                       * Vt );
            */
#else
            // Force pure rotation U,Vt by fixing potential inversion/reflection see LeafDSH_FEM
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
            }
            Mat2x2r R( U*Vt );
#endif
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    /* Implementation of the temporally-coherent heuristic eigenvalue
       negation idea in IHFSDM, assume all inversions are detected
       across (1,0) direction.
    */
    Vec6r ComputeElasticForce_IHFSDM( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r diag_S;
            mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
            // Force pure rotation U,Vt by fixing potential inversion/reflection see LeafDSH_FEM
            if( mal::Det(F) < 0 ) //TEMP: Unnecessary
            {
                /* Invert diag_S entry whose eigenvector in V is most
                   aligned with global (1,0), which is the unambiguous
                   inversion axis in the default element setting
                   (symmetric, beta=0)
                */
                //TEMP: This assumed static (1,0) axis... if( mal::Abs(mal::GRow<0>(Vt).x()) < mal::Abs(mal::GRow<1>(Vt).x()) ) //V.Column[1] more aligned with (1,0) than V.Column[0]

                //\todo mmm, I don't really get it... which axis
                //determines alignment that allows choosing smallest
                //eigenvalue?... I think this is NOT correct, but
                //yields "reasonable" results, similar to
                //non-alpha-rotated case... MUST reread article and try to
                //clarify it...

//#define __ENABLE_HACKISH_IHFSDM
#ifdef __ENABLE_HACKISH_IHFSDM
                /*\note this works poperly for x0-inverted, and MUST use reference configuration because V is "local" to it
                // V.Column[1] more aligned with axis12 than V.Column[0]
                Vec2r axis12( mal::PerpendicularCW( m_r[2] - m_r[1] ) );
                if( mal::Abs( mal::Dot( mal::GRow<0>(Vt), axis12 ) ) < mal::Abs( mal::Dot( mal::GRow<1>(Vt), axis12 ) ) )
                */

                /*note This yields a reasonable force field, but I'm not sure why
                  U.Column[1] more aligned with axis12 than U.Column[0]
                */
                Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) );
                if( mal::Abs( mal::Dot( mal::GColumn<0>(U), axis12 ) ) < mal::Abs( mal::Dot( mal::GColumn<1>(U), axis12 ) ) )
                {
                    //if( mal::Det(U) * mal::Det(Vt) < 0 )
                    {
                        if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                        if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                    }
                }
                else
                {
                    //if( mal::Det(U) * mal::Det(Vt) < 0 )
                    {
                        if( mal::Det(U) < 0 ) mal::GSetColumn<0>( U, -mal::GColumn<0>(U) );
                        if( mal::Det(Vt) < 0 ) mal::GSetRow<0>( Vt, -mal::GRow<0>(Vt) );
                    }
                }
#else
                // Find *geometrically* shortest inversion direction in reference configuration
                //\todo NOTICE that the selected inversion direction changes with dist12
                Vec2r n_c[3];
                n_c[0] = mal::PerpendicularCW( m_r[2] - m_r[1] );
                n_c[1] = mal::PerpendicularCW( m_r[0] - m_r[2] );
                n_c[2] = mal::PerpendicularCW( m_r[1] - m_r[0] );
                Vec2r v[2];
                v[0] = mal::GRow<0>(Vt);
                v[1] = mal::GRow<1>(Vt);
                Real min_lambda( mal::Infinity<Real>() );
                int min_it_v(0);
                //Real lambda_CxV[3][2];
                for( int it_c=0;
                     it_c < 3; //it_c<1; //\todo Uncomment to see the case where c = NoC = v0, clearer discontinuity
                     it_c++ )
                {
                    for( int it_v=0; it_v<2; it_v++ )
                    {
                        Real dot_c_v = mal::Dot( v[it_v], n_c[it_c] );
                        Real lambda_c_v = mal::Abs(dot_c_v) > 0.0001f
                                          ? mal::Dot( m_r[(it_c+1)%3] - m_r[it_c], n_c[it_c] ) / dot_c_v
                                          : mal::Infinity<Real>();
                        lambda_c_v = mal::Abs(lambda_c_v);
                        //lambda_CxV[it_c][it_v] = lambda_c_v;
                        if( lambda_c_v < min_lambda )
                        {
                            min_lambda = lambda_c_v;
                            min_it_v = it_v;
                        }
                    }
                }
                // Invert selected direction
                if( 0 == min_it_v )
                {
                    //if( mal::Det(U) * mal::Det(Vt) < 0 )
                    {
                        if( mal::Det(U) < 0 ) mal::GSetColumn<0>( U, -mal::GColumn<0>(U) );
                        if( mal::Det(Vt) < 0 ) mal::GSetRow<0>( Vt, -mal::GRow<0>(Vt) );
                    }
                }
                else
                {
                    //if( mal::Det(U) * mal::Det(Vt) < 0 )
                    {
                        if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                        if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                    }
                }
#endif
            }
            Mat2x2r R( U*Vt );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    /*\note Actually, we follow ITF usage, but compute f0 as in SIGG12 course notes
      Also, we base deformation at x0, not at x2 as SIGG12, thus, final f0 is computed symmetrically

     LRM: (ECIE pg 5, also ITF pg 5)
       Energy = \mu \sum_i (\theta_i - 1)^2 + \frac{\lambda}{2} (\sum_i (\theta_i - 1) )^2
    */
    Vec2r ComputeDiagonalP_ITF_LRM( const Vec2r &vec_diag_F, Real *p_energy = 0 ) const
        {
            /* TEMP: Literal implementation from ITF, slower
            Mat2x2r diag_F( vec_diag_F[0], 0,
                            0, vec_diag_F[1] );
            Mat2x2r E( diag_F - Mat2x2r::Identity() );
            Mat2x2r diag_P( 2*m_LameMu * E
                            + m_LameLambda * mal::Trace( E ) * Mat2x2r::Identity() );
            */
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_rParams.m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_rParams.m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            if( p_energy ) *p_energy = lame_mu * ( mal::Sq(vec_diag_F[0]-Real(1)) + mal::Sq(vec_diag_F[1]-Real(1)) )
                                       + Real(0.5)*lame_lambda*mal::Sq( (vec_diag_F[0]-Real(1)) + (vec_diag_F[1]-Real(1)) );
            Real trE( vec_diag_F[0] + vec_diag_F[1] - 2 );
            return Vec2r( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * trE,
                          2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * trE );
        }
    Vec6r ComputeElasticForce_ITF_LRM( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );

            /*TEMP: Clamp F to compute U*Vt uninverted
            F = Compute_F( X, Params::eFFM_Project );
            // Compute SVD
            Vec2r dummy_vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, dummy_vec_diag_F, Vt );
            //END TEMP
            */

            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
                //\todo Vega implements this correction for V, which yields DISCONTINUOUS eigenvectors
                //  if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); }
            }

            /*TEMP: Clamp F to compute U*Vt uninverted
            F = Compute_F( X, Params::eFFM_Project );
            // Compute SVD
            Vec2r dummy_vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, dummy_vec_diag_F, Vt );
            //END TEMP
            */

            if( p_rotation ) *p_rotation = U*Vt;
            // Material = Linear Rotated Model == Corrotated (from ECIE definition)
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ITF_LRM( vec_diag_F, p_energy ) ) );
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
        }

    /* Corrotational correction (improves ITF with different energy terms)

     \todo Notice that in the ECIE Paper, a computation U * diag_P * Vt is
     suggested, but in the Supplementary Technical Document, it's not
     so clear, as they parametrize U and Vt using "Rodriges' Rotation
     Formula"...

     CLRM: (ECIE pg 5)
       Energy = \mu \sum_i (\theta_i - 1)^2 + \frac{\lambda}{2} (J - 1)^2
       J = \theta_0 \theta_1
    */
    Vec2r ComputeDiagonalP_ECIE_CLRM( const Vec2r &vec_diag_F, Real *p_energy = 0 ) const
        {
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_rParams.m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_rParams.m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real J( vec_diag_F[0] * vec_diag_F[1] );
            if( p_energy ) *p_energy = lame_mu * ( mal::Sq(vec_diag_F[0]-Real(1)) + mal::Sq(vec_diag_F[1]-Real(1)) )
                                       + Real(0.5)*lame_lambda*mal::Sq( J - Real(1) );
            return Vec2r( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * (J-1) * vec_diag_F[1],
                          2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * (J-1) * vec_diag_F[0] );
        }
    Vec6r ComputeElasticForce_ECIE_CLRM( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
            }
            if( p_rotation ) *p_rotation = U*Vt;
            // Material = CLRM
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ECIE_CLRM( vec_diag_F ,p_energy ) ) );
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
        }

    Vec6r ComputeElasticForce_PD_CLRM( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {

//#define __OLD
#ifdef __OLD //TEMP: To check equivalence between hyperelastic and PK versions
            Mat2x2r F( Compute_F( X, Params::eFFM_None ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
            }
            Mat2x2r R( U*Vt );
            if( p_rotation ) *p_rotation = R;
            Mat2x2r S( R.Transposed() * F ); //F = R*S
            Real det_S( mal::Det(S) ); //\note == det_F
#else
            // Compute PD_U R
            Mat2x2r R;
            ComputeElasticForce_PD_Unified( X, &R );
            if( p_rotation ) *p_rotation = R;
            Mat2x2r F( Compute_F( X, Params::eFFM_None ) );
            Mat2x2r S( R.Transposed() * F ); //F = R*S
            Real det_S( mal::Det(S) ); //\note == det_F
#endif
            /*TEMP The "simple" P formula without factoring R,S from F, P only fails to be continuous across det(F) = 0 (see vector field)
            Mat2x2r P( 2*m_LameMu*(F-R) );
            if( mal::Abs(det_F) > mal::Epsilon<Real>() ) P += m_LameLambda * (det_F-1) * det_F * mal::Inverse( mal::Transposed(F) );
            */

            // (S^*)^T keeping the S terms and cancelling the det_S, P is continuous across det(F) = 0
            Mat2x2r P( R * ( 2*m_LameMu*(S-Mat2x2r::Identity())
                             + m_LameLambda * (det_S-1) * Mat2x2r( S(1,1), -S(1,0),
                                                                   -S(0,1), S(0,0) ) ) );
            //

            /*\note This is ALSO correct, factoring/cancelling the det_F terms and using (F^*)^T instead
            Mat2x2r P( 2*m_LameMu*(F-R)
                       + m_LameLambda * (det_F-1) * Mat2x2r( F(1,1), -F(1,0),
                                                             -F(0,1), F(0,0) ) );
            */

            /*\note Tikhonov-regularized version... non-reciprocal
              det_F factor ALSO needs to be reg to avoid cancelling
              the whole second additive term of P This works and may
              be the solution to the dP differential problems with
              1/det(F)
            Real det_F_reg( ( mal::Abs(det_F) < mal::Epsilon<Real>() )
                            ? ( ( det_F < 0 )
                                ? det_F - mal::Epsilon<Real>()
                                : det_F + mal::Epsilon<Real>() )
                            : det_F );
            Mat2x2r P( 2*m_LameMu*(F-R)
                       + m_LameLambda * (det_F-1) * det_F_reg * mal::Rcp( det_F_reg ) * Mat2x2r( F(1,1), -F(1,0),
                                                                                                -F(0,1), F(0,0) ) );
            */

            /* \todo S^* es INCORRECTE, pq S NO ES SIMETRICA, però dona lloc a un camp de forces MOLT INTERESSANT
            Mat2x2r P( R * ( 2*m_LameMu*(S-Mat2x2r::Identity())
                             + m_LameLambda * (det_S-1) * Mat2x2r( S(1,1), -S(0,1),
                                                                   -S(1,0), S(0,0) ) ) );
            */
            /*\todo Elastic energy is tricky here... first term is
              taken from the Siggraph2012 course notes on LinearFEM
              but using "S", the second one is taken from CLRM "J"
              term using det_S instead. I think it's correct, but the
              resulting energy is not conserved, and CLRM one is.

              The results seem to be quite depending on the timestep,
              which means it's highly nonlinear... for dt=0.0001
              trajectories are much different than dt=0.001
            */
            if( p_energy ) *p_energy = m_LameMu * ( S-Mat2x2r::Identity() ).NormSqF()  //|S-I|^2_Frobenius
                                       + Real(0.5)*m_LameLambda*mal::Sq( det_S - Real(1) );
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy

            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm );
            Mat2x2r H( - m_Area * P * Be0.Transposed() ); // Negation is present in the Siggraph2012 course notes, but not in ITF...
            Vec2r f0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r f( f0.x(), f0.y(),
                     H(0,0), H(1,0),
                     H(0,1), H(1,1) );
            return f;
        }

    /* ECIE C0 extension to neo-hookean, clamping eigenvalues < ecie_e and J accordingly.
       NOT strictly demonstrated in ECIE, but suggested in other articles from the same autors
    */
    Vec2r ComputeDiagonalP_ECIE_NHC0( const Vec2r &vec_diag_F, Real *p_energy = 0 ) const
        {
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_rParams.m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_rParams.m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real ecie_e( m_rParams.m_ECIE_e_threshold );
            Vec2r vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                                      mal::Max( ecie_e, vec_diag_F[1] ) );
            Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
            if( p_energy ) *p_energy = Real(0.5)*lame_mu * ( mal::Sq(vec_clamped_diag_F[0]) + mal::Sq(vec_clamped_diag_F[1]) - Real(2) )
                                       - lame_mu * mal::Log(clamped_J)
                                       + Real(0.5)*lame_lambda*mal::Sq( mal::Log(clamped_J) );
            Real g0( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
            Real g1( lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
            /*TEMP: Trying to extract proper energy from NHC0 we see
              that it's NOT consistent... we cannot extrapolate it, if
              forces are clamped, so is energy? and the result is that
              initial TotalEnergy is not preserved...
            */
            if( p_energy )
            {
                Real delta0( mal::Min<Real>( vec_diag_F[0] - ecie_e, 0 ) );
                Real delta1( mal::Min<Real>( vec_diag_F[1] - ecie_e, 0 ) );
                Real h01( lame_lambda / clamped_J );
                *p_energy += g0 * delta0
                             + g1 * delta1
                             + h01 * delta0 * delta1;
            }
            return Vec2r( g0, g1 );
        }
    Vec6r ComputeElasticForce_ECIE_NHC0( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
            }
            if( p_rotation ) *p_rotation = U*Vt;
            // Material = C0 neohookean extension
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ECIE_NHC0( vec_diag_F, p_energy ) ) );
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
        }

    /* ECIE C1 extension to neo-hookean
       Derivatives and notation from ECIE Supplementary Technical Document
    */
    Vec2r ComputeDiagonalP_ECIE_NHC1( const Vec2r &vec_diag_F, Real *p_energy = 0 ) const
        {
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_rParams.m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_rParams.m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real ecie_e( m_rParams.m_ECIE_e_threshold );
            Real ecie_k( m_rParams.m_ECIE_k_factor * m_rParams.m_YoungModulus );
            bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
            bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
            Vec2r vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                                      mal::Max( ecie_e, vec_diag_F[1] ) );
            Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
            if( p_energy ) // clamped energy
                *p_energy = Real(0.5)*lame_mu * ( mal::Sq(vec_clamped_diag_F[0]) + mal::Sq(vec_clamped_diag_F[1]) - Real(2) )
                            - lame_mu * mal::Log(clamped_J)
                            + Real(0.5)*lame_lambda*mal::Sq( mal::Log(clamped_J) );
            Real g0( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
            Real g1( lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
            Real h01( lame_lambda / clamped_J );
            Vec2r vec_diag_P( g0, g1 );
            //\todo Cases can be generalized if delta0,delta1 are 0 for non-clamped axis
            if( bIsDegenerate0 && bIsDegenerate1 )
            {
                Real delta0( vec_diag_F[0] - ecie_e );
                Real delta1( vec_diag_F[1] - ecie_e );
                vec_diag_P[0] += h01*delta1 + 2*ecie_k*delta0;
                vec_diag_P[1] += h01*delta0 + 2*ecie_k*delta1;
                if( p_energy ) // extrapolated energy
                    *p_energy += g0 * delta0
                                 + g1 * delta1
                                 + h01 * delta0 * delta1
                                 + ecie_k * ( mal::Sq(delta0) + mal::Sq(delta1) );
            }
            else if( bIsDegenerate0 && !bIsDegenerate1 )
            {
                Real delta0( vec_diag_F[0] - ecie_e );
                vec_diag_P[0] += 2*ecie_k*delta0;
                vec_diag_P[1] += h01*delta0;
                if( p_energy ) // extrapolated energy
                    *p_energy += g0 * delta0 + ecie_k*mal::Sq(delta0);
            }
            else if( !bIsDegenerate0 && bIsDegenerate1 )
            {
                Real delta1( vec_diag_F[1] - ecie_e );
                vec_diag_P[0] += h01*delta1;
                vec_diag_P[1] += 2*ecie_k*delta1;
                if( p_energy ) // extrapolated energy
                    *p_energy += g1 * delta1 + ecie_k*mal::Sq(delta1);
            }
            return vec_diag_P;
        }
    Vec6r ComputeElasticForce_ECIE_NHC1( const Vec2r X[3], Mat2x2r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );
            // Force pure rotation U,Vt by fixing potential inversion/reflection
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) { mal::GSetColumn<1>( U, -mal::GColumn<1>(U) ); vec_diag_F[1] = -vec_diag_F[1]; }
                if( mal::Det(Vt) < 0 ) { mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) ); vec_diag_F[1] = -vec_diag_F[1]; }
            }
            if( p_rotation ) *p_rotation = U*Vt;
            // Material = Corrected Linear Rotated Model == Corrected Corrotated
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ECIE_NHC1( vec_diag_F, p_energy ) ) );
            if( p_energy ) *p_energy *= m_Area; // Convert energy density to energy
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
        }

public:
    friend class AppRenderer;
    friend class AppTask;
    // Material params
    Real m_LameMu;
    Real m_LameLambda;
    // Constant linear triangle stuff
    Vec2r m_r[3];
    Mat2x2r m_invDm;
    Mat3x6r m_B;
    Mat6x6r m_Ke;
    Mat2x6r m_Ke0;
    Real m_Area;
    // Variable triangle node pos
    Vec2r m_X[3];
    int32 m_NoC; //Current NodeOfCollapse, -1 if none
};

//TEMP: Code from LeafDSH_FEM_Solid2D, don't change it here...
#ifdef __USE_NOC
int ComputeNoC( const Vec2r X0[3], const Vec2r X1[3], Real dtdf, Real area )
{
    int noc(-1);

    // Gather node positions and displacements
    Vec2r p0( X0[0] );
    Vec2r p1( X0[1] );
    Vec2r p2( X0[2] );
    Vec2r d0( X1[0] - p0 );
    Vec2r d1( X1[1] - p1 );
    Vec2r d2( X1[2] - p2 );
    // Compute ToC, when det(F) = m_rParams.m_DegenerateThresholdDetF => det(Q) = m_rParams.m_DegenerateThresholdDetF * det(P) => det(Q) = m_rParams.m_DegenerateThresholdDetF * 2 * Area(P)
    Real D( dtdf * 2 * area );
    Real a(   d1.x() * d2.y()
              + d2.x() * d0.y()
              + d0.x() * d1.y()
              - d1.x() * d0.y()
              - d0.x() * d2.y()
              - d2.x() * d1.y() );
    Real b( (   d1.x() * p2.y() + p1.x() * d2.y() )
            + ( d2.x() * p0.y() + p2.x() * d0.y() )
            + ( d0.x() * p1.y() + p0.x() * d1.y() )
            - ( d1.x() * p0.y() + p1.x() * d0.y() )
            - ( d0.x() * p2.y() + p0.x() * d2.y() )
            - ( d2.x() * p1.y() + p2.x() * d1.y() ) );
    Real c(   p1.x() * p2.y()
              + p2.x() * p0.y()
              + p0.x() * p1.y()
              - p1.x() * p0.y()
              - p0.x() * p2.y()
              - p2.x() * p1.y()
              - D );
    Real toc1(0), toc2(1);
    if( mal::GSolvePolynomialEq2<Real>( a, b, c, toc1, toc2 ) )
    {
        //APP_LOG_WARNING( "Potentially degenerated %d with toc = %f, %f", it_e, toc1, toc2 );
        APP_ASSERT( toc1 <= toc2 );
        //\todo We knot that det(F(0)) > DTDF, and that det(F(1)) <= DTDF, so there MUST BE STRICTLY 1 crossing in toc = [0,1]
        Real toc( toc1 > 0 ? toc1 : toc2 );
        if( toc > 1 )
        {
            APP_LOG_ERROR( "toc %f > 1", toc );
        }
        /* analyze inversion interval:
           0 <= toc1 <= toc2 < 1 => Collapse->Invert->Uninvert->Uncollapse
           toc1 < 0 <= toc2 < 1 => Collapse->Invert
           toc1 <= toc2 < 0 => Extrapolated collapse outside [0,1], no actual collapse
           1 <= toc1 <= toc2 < 1 => Extrapolated collapse outside [0,1], no actual collapse
        */
        /*TEMP
          if( !(toc2 < 0 || toc1 > 1) //Inside [0,1]
          &&
          (toc1 < 0 || toc2 > 1) //2 roots inside [0,1] //\todo MAY UNINVERT ACROSS AXIS DIFERENT than collapsed one!!
          )
        */
        {
            // Collapse configuration
            Vec2r q0( p0 + toc*d0 );
            Vec2r q1( p1 + toc*d1 );
            Vec2r q2( p2 + toc*d2 );
            // Compute NoC
#ifdef __DISABLED_BY_NOW
            switch( m_Params.m_NCM )
            {
            case Params::eNCM_IHFSDM:
                {
                    // Compute current node positions and F
                    Vec2r vec_pos[3];
                    vec_pos[0] = q0;
                    vec_pos[1] = q1;
                    vec_pos[2] = q2;
                    Mat2x2 F;
                    element.ComputeF( vec_pos[0], vec_pos[1], vec_pos[2], F );
                    // Compute SVD
                    Mat2x2 U, Vt;
                    Vec2r diag_F;
                    mal::GSingularValueDecomposition_USVt( F, U, diag_F, Vt );
                    // Find *geometrically* shortest inversion direction in reference configuration
                    //\todo NOTICE that the selected inversion direction changes with dist12
                    Vec2r n_c[3];
                    n_c[0] = mal::PerpendicularCW( vec_pos[2] - vec_pos[1] );
                    n_c[1] = mal::PerpendicularCW( vec_pos[0] - vec_pos[2] );
                    n_c[2] = mal::PerpendicularCW( vec_pos[1] - vec_pos[0] );
                    Vec2r v[2];
                    v[0] = mal::GRow<0>(Vt);
                    v[1] = mal::GRow<1>(Vt);
                    Real min_lambda( mal::Infinity<Real>() );
                    int min_it_v(0);
                    //Real lambda_CxV[3][2];
                    for( int it_c=0; it_c<3; it_c++ )
                    {
                        for( int it_v=0; it_v<2; it_v++ )
                        {
                            Real dot_c_v = mal::Dot( v[it_v], n_c[it_c] );
                            Real lambda_c_v = mal::Abs(dot_c_v) > 0.0001f
                                              ? mal::Dot( vec_pos[(it_c+1)%3] - vec_pos[it_c], n_c[it_c] ) / dot_c_v
                                              : mal::Infinity<Real>();
                            lambda_c_v = mal::Abs(lambda_c_v);
                            //lambda_CxV[it_c][it_v] = lambda_c_v;
                            if( lambda_c_v < min_lambda )
                            {
                                min_lambda = lambda_c_v;
                                min_it_v = it_v;
                                noc = it_c;
                            }
                        }
                    }
                }
                break;
            case Params::eNCM_ToC:
            default: //\note Default => ToC, used ony for Viz if m_NCM != eNCM_ToC
#endif //__DISABLED_BY_NOW
                {
                    Real sql01( mal::NormSq( q1-q0 ) );
                    Real sql12( mal::NormSq( q2-q1 ) );
                    Real sql20( mal::NormSq( q0-q2 ) );
                    // Choose NoC, vertex opposite to longest edge at ToC
                    noc = (sql01 >= sql12)
                          ? (sql01 >= sql20) ? 2 : 1
                          : (sql12 >= sql20) ? 0 : 1;
                }
#ifdef __DISABLED_BY_NOW
                break;
            }
#endif
        }
        /*TEMP
          else
          {
          APP_LOG_WARNING( "No 0 or 2 degenerations in range [0,1]" );
          }
        */
    }
    else
    {
        APP_LOG_ERROR( "New Degenerated CANNOT SOLVE Eq2" );
    }
    return noc;
}
#endif

#endif //TEST_FE_ELEMENT2D_H
