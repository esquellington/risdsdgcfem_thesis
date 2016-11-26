#ifndef TEST_FE_ELEMENT3D_H
#define TEST_FE_ELEMENT3D_H

#include "Config.h"
#include "Params.h"
#include <Mal/GMatDecomposition.h>
#include <Mal/GSolvePolynomialEq.h>

#define __USE_FAST_SVD_PHYSBAM
#ifdef __USE_FAST_SVD_PHYSBAM
#  include <Mal/GMatDecomposition.h>
#else
#  define __USE_SVD_GSL
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_permutation.h>
#  include <gsl/gsl_linalg.h>
#endif

#include <gsl/gsl_poly.h>

namespace mal
{

#ifdef __USE_SVD_GSL
uint32 g_gsl_ref_count(0);
static gsl_matrix* g_gsl_U_3x3(0);
static gsl_matrix* g_gsl_V_3x3(0);
static gsl_vector* g_gsl_S_3(0);
static gsl_vector* g_gsl_Work_3(0);
static bool GSL_IsInitialized() { return g_gsl_ref_count > 0; }
static void GSL_Init()
{
    if( !GSL_IsInitialized() )
    {
        g_gsl_U_3x3 = gsl_matrix_alloc( 3, 3 );
        g_gsl_V_3x3 = gsl_matrix_alloc( 3, 3 );
        g_gsl_S_3 = gsl_vector_alloc( 3 );
        g_gsl_Work_3 = gsl_vector_alloc( 3 );
    }
    g_gsl_ref_count++;
}
static void GSL_ShutDown()
{
    if( GSL_IsInitialized() )
    {
        g_gsl_ref_count--;
        if( 0 == g_gsl_ref_count )
        {
            gsl_vector_free( g_gsl_Work_3 );
            gsl_vector_free( g_gsl_S_3 );
            gsl_matrix_free( g_gsl_V_3x3 );
            gsl_matrix_free( g_gsl_U_3x3 );
            g_gsl_U_3x3 = 0;
            g_gsl_V_3x3 = 0;
            g_gsl_S_3 = 0;
            g_gsl_Work_3 = 0;
        }
    }
}
#endif

//Solve ax^3 + bx^2 + cx + d = 0 returns x1 <= x2 <= x3
template<typename T> int GSolvePolynomialEq3( T a, T b, T c, T d, T& x1, T& x2, T& x3 )
{
    //\see http://www.gnu.org/software/gsl/manual/html_node/Cubic-Equations.html#Cubic-Equations
    if( a != T(0) ) //TEMP: Should use some epsilon....
    {
        double r1,r2,r3;
        int num_roots = gsl_poly_solve_cubic( b/a, c/a, d/a, &r1, &r2, &r3 );
        x1 = r1;
        x2 = r2;
        x3 = r3;
        return num_roots;
    }
    else
    {
        x3 = 0;
        return GSolvePolynomialEq2( b, c, d, x1, x2 );
    }
}

template <typename T>
inline void GQRDecomposition( const GMat<T,3,3> &F, GMat<T,3,3> &QR_Q, GMat<T,3,3> &QR_R )
{
#ifdef __USE_GSL_QR
    //\see http://www.gnu.org/software/gsl/manual/html_node/QR-Decomposition.html
    // Compute (QR,tau)
    gsl_matrix *gslQR = gsl_matrix_alloc( 3, 3 );
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            gsl_matrix_set( gslQR, i, j, F(i,j) );
    gsl_vector *gslTau = gsl_vector_alloc( 3 );
    int res = gsl_linalg_QR_decomp( gslQR, gslTau );
    if(res) printf ("error: %s\n", gsl_strerror (res));
    // Unpack (QR,tau) into (QR_Q,QR_R)
    gsl_matrix *gslQ = gsl_matrix_alloc( 3, 3 );
    gsl_matrix *gslR = gsl_matrix_alloc( 3, 3 );
    res = gsl_linalg_QR_unpack( gslQR, gslTau, gslQ, gslR );
    if(res) printf ("error: %s\n", gsl_strerror (res));
    // Save Q and R
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
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
    gsl_matrix *gslF = gsl_matrix_alloc( 3, 3 );
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            gsl_matrix_set( gslF, i, j, F(i,j) );
    gsl_matrix *gslQ = gsl_matrix_alloc( 3, 3 );
    gsl_matrix *gslR = gsl_matrix_alloc( 3, 3 );
    gsl_vector *gslTau = gsl_vector_alloc( 3 );
    int gslSignum(0);
    gsl_permutation *gslP = gsl_permutation_alloc( 3 );
    gsl_vector *gslNorm = gsl_vector_alloc( 3 );
    int res = gsl_linalg_QRPT_decomp2( gslF, gslQ, gslR, gslTau, gslP, &gslSignum, gslNorm );
    if(res) printf ("error: %s\n", gsl_strerror (res));
    // Save QR_Q as R
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
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
    QR_Q = GMat<T,3,3>::Identity();
    QR_R = GMat<T,3,3>::Identity();
#endif
}

template <typename T>
inline GMat<T,3,3> GRotation3x3_QR( const GMat<T,3,3> &F )
{
    GMat<T,3,3> QR_Q,QR_R;
    GQRDecomposition( F, QR_Q, QR_R );
    return QR_Q; //funny, but QR_Q is actually R
    /*
    Real det_Q( mal::Det(QR_Q) );
    if( det_Q > 0 )
        return QR_Q; //funny, but QR_Q is actually R
    else
    {
        //
        //APP_LOG_WARNING( "QR with det %f", det_Q );
        return mal::GRotation3x3_GramSchmidtOrthonormalization_XYZ( F, Mat3x3r::Identity() );
        //
    }
    */
}

#ifdef __USE_SVD_GSL
template <typename T>
inline void GSingularValueDecomposition_USVt( const GMat<T,3,3> &F,
                                              GMat<T,3,3> &U, GVec<T,3> &diag_F, GMat<T,3,3> &Vt )
{
    //\see http://www.gnu.org/software/gsl/manual/html_node/Singular-Value-Decomposition.html
    APP_ASSERT( GSL_IsInitialized() );
    // Avoids alloc/free per SVD-call
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            gsl_matrix_set( g_gsl_U_3x3, i, j, F(i,j) );
    int res = gsl_linalg_SV_decomp( g_gsl_U_3x3, g_gsl_V_3x3, g_gsl_S_3, g_gsl_Work_3 );
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            U(i,j) = gsl_matrix_get( g_gsl_U_3x3, i, j );
            Vt(i,j) = gsl_matrix_get( g_gsl_V_3x3, j, i );
        }
        diag_F[i] = gsl_vector_get( g_gsl_S_3, i );
    }

    /*\todo AVOID ALLOC/FREE per call!!
    gsl_matrix *gslU = gsl_matrix_alloc( 3, 3 );
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            gsl_matrix_set( gslU, i, j, F(i,j) );
    gsl_matrix *gslV = gsl_matrix_alloc( 3, 3 );
    gsl_vector *gslS = gsl_vector_alloc( 3 );
    gsl_vector *gslWork = gsl_vector_alloc( 3 );
    int res = gsl_linalg_SV_decomp( gslU, gslV, gslS, gslWork );
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            U(i,j) = gsl_matrix_get( gslU, i, j );
            Vt(i,j) = gsl_matrix_get( gslV, j, i );
        }
        diag_F[i] = gsl_vector_get( gslS, i );
    }
    gsl_vector_free( gslWork );
    gsl_vector_free( gslS );
    gsl_matrix_free( gslV );
    gsl_matrix_free( gslU );
    */
}
#endif

} //namespace mal

#ifdef __TODO_PORT_TO_3D
inline Mat2x2r RLerp( const Mat2x2r &R0, const Mat2x2r &R1, Real lambda01 )
{
    Real w1( lambda01 );
    Real w0( Real(1)-w1 );
    Mat2x2r LR( w0*R0 + w1*R1 );
    Real det_LR( mal::Det(LR) );
    return GRotation2x2_PolarDecomposition( LR, det_LR ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix
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
    return GRotation2x2_PolarDecomposition( CR, det_CR ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix

    /* TEMP: non-orthonormal CR yields reasonable forcefield and
       det(CR) > 0.65 for all observed scenarios, it MAY be possible
       to reduce CPU (significantly, in 3D) cost by skipping PD(CR)
       with negligible impact in accuracy
    if( !mal::ApproxEq( det_CR, Real(1), Real(0.1) ) ) APP_LOG_WARNING("RCerp rotation severely non-orthonormal, det(CR) = %f", det_CR );
    return CR;
    */
}
#endif //__TODO_PORT_TO_3D

class Element3D
{
private:
    const Params &m_rParams;
public:
    Element3D( const Params &params )
    : m_rParams(params)
    , m_LameMu(0)
    , m_LameLambda(0)
    , m_Volume(0)
    , m_NoC(0)
        {
#ifdef __USE_SVD_GSL
            mal::GSL_Init();
#endif
        }
    ~Element3D()
        {
            Clear();
#ifdef __USE_SVD_GSL
            mal::GSL_ShutDown();
#endif
        }

    void Clear()
        {
        }
    void RebuildFromParams()
        {
            Clear();
            LameParameters_From_YoungAndPoisson( m_rParams.m_YoungModulus, m_rParams.m_PoissonRatio, m_LameMu, m_LameLambda );
            // Init element
#define __APPLY_ROTATION_ALPHA_Z
#ifdef __APPLY_ROTATION_ALPHA_Z
            Mat3x3r R_alpha( mal::GRotation3x3_From( Vec3r(0,0,1), m_rParams.m_Alpha ) );
            /*
            m_r[0] = R_alpha * Real(1) * Vec2r( mal::Cos(m_rParams.m_Angle0), mal::Sin(m_rParams.m_Angle0) );
            Vec2r dr( Real(0.5) * m_rParams.m_Dist12 * Vec2r(0,1) );
            m_r[1] = R_alpha * dr;
            m_r[2] = -R_alpha * dr;
            */
#else
            Mat3x3r R_alpha( Mat3x3r::Identity() );
#endif
            m_r[0] = R_alpha * Vec3r(1,0,0);
            m_r[1] = R_alpha * ( 0.5*m_rParams.m_Dist12*Vec3r(0,1,0) );
            Real cos30deg( mal::Cos( mal::SixthPi<Real>() ) );
            Real sin30deg( mal::Sin( mal::SixthPi<Real>() ) );
            m_r[2] = R_alpha * ( 0.5*m_rParams.m_Dist12*Vec3r(0,-sin30deg,-cos30deg) );
            m_r[3] = R_alpha * ( 0.5*m_rParams.m_Dist12*Vec3r(0,-sin30deg, cos30deg) );
            Mat3x3r Dm( mal::GMat3x3_From_Columns( m_r[1] - m_r[0], m_r[2] - m_r[0], m_r[3] - m_r[0] ) );
            m_invDm = mal::Inverse( Dm );
            // Compute Ke
            Vec3r r0( m_r[0] );
            Vec3r r1( m_r[1] );
            Vec3r r2( m_r[2] );
            Vec3r r3( m_r[3] );
            m_Volume = mal::Det( Dm ) / Real(6);
            APP_ASSERT( m_Volume > Real(0) );
            //invJ from AFEM-16 p.7
            Mat4x4r invJ;
            Real x1(r0.x()),x2(r1.x()),x3(r2.x()),x4(r3.x());
            Real y1(r0.y()),y2(r1.y()),y3(r2.y()),y4(r3.y());
            Real z1(r0.z()),z2(r1.z()),z3(r2.z()),z4(r3.z());
            invJ(0,0) = x2*(y3*z4 - y4*z3) + x3*(y4*z2-y2*z4) + x4*(y2*z3 - y3*z2);
            invJ(1,0) = x1*(y4*z3 - y3*z4) + x3*(y1*z4-y4*z1) + x4*(y3*z1 - y1*z3);
            invJ(2,0) = x1*(y2*z4 - y4*z2) + x2*(y4*z1-y1*z4) + x4*(y1*z2 - y2*z1);
            invJ(3,0) = x1*(y3*z2 - y2*z3) + x2*(y1*z3-y3*z1) + x3*(y2*z1 - y1*z2);
            invJ(0,1) = (y4-y2)*(z3-z2) - (y3-y2)*(z4-z2); invJ(0,2) = (x3-x2)*(z4-z2) - (x4-x2)*(z3-z2); invJ(0,3) = (x4-x2)*(y3-y2) - (x3-x2)*(y4-y2);
            invJ(1,1) = (y3-y1)*(z4-z3) - (y3-y4)*(z1-z3); invJ(1,2) = (x4-x3)*(z3-z1) - (x1-x3)*(z3-z4); invJ(1,3) = (x3-x1)*(y4-y3) - (x3-x4)*(y1-y3);
            invJ(2,1) = (y2-y4)*(z1-z4) - (y1-y4)*(z2-z4); invJ(2,2) = (x1-x4)*(z2-z4) - (x2-x4)*(z1-z4); invJ(2,3) = (x2-x4)*(y1-y4) - (x1-x4)*(y2-y4);
            invJ(3,1) = (y1-y3)*(z2-z1) - (y1-y2)*(z3-z1); invJ(3,2) = (x2-x1)*(z1-z3) - (x3-x1)*(z1-z2); invJ(3,3) = (x1-x3)*(y2-y1) - (x1-x2)*(y3-y1);
            invJ = mal::Rcp( Real(6) * m_Volume ) * invJ;

            Real a1(invJ(0,1)), b1(invJ(0,2)), c1(invJ(0,3));
            Real a2(invJ(1,1)), b2(invJ(1,2)), c2(invJ(1,3));
            Real a3(invJ(2,1)), b3(invJ(2,2)), c3(invJ(2,3));
            Real a4(invJ(3,1)), b4(invJ(3,2)), c4(invJ(3,3));

            //B from AFEM-16 p.13
            m_B(0,0) = a1; m_B(0,1) =  0; m_B(0,2) =  0; m_B(0,3) = a2; m_B(0,4) =  0; m_B(0,5) =  0; m_B(0,6) = a3; m_B(0,7) =  0; m_B(0,8) =  0; m_B(0,9) = a4; m_B(0,10) =  0; m_B(0,11) =  0;
            m_B(1,0) =  0; m_B(1,1) = b1; m_B(1,2) =  0; m_B(1,3) =  0; m_B(1,4) = b2; m_B(1,5) =  0; m_B(1,6) =  0; m_B(1,7) = b3; m_B(1,8) =  0; m_B(1,9) =  0; m_B(1,10) = b4; m_B(1,11) =  0;
            m_B(2,0) =  0; m_B(2,1) =  0; m_B(2,2) = c1; m_B(2,3) =  0; m_B(2,4) =  0; m_B(2,5) = c2; m_B(2,6) =  0; m_B(2,7) =  0; m_B(2,8) = c3; m_B(2,9) =  0; m_B(2,10) =  0; m_B(2,11) = c4;
            m_B(3,0) = b1; m_B(3,1) = a1; m_B(3,2) =  0; m_B(3,3) = b2; m_B(3,4) = a2; m_B(3,5) =  0; m_B(3,6) = b3; m_B(3,7) = a3; m_B(3,8) =  0; m_B(3,9) = b4; m_B(3,10) = a4; m_B(3,11) =  0;
            m_B(4,0) =  0; m_B(4,1) = c1; m_B(4,2) = b1; m_B(4,3) =  0; m_B(4,4) = c2; m_B(4,5) = b2; m_B(4,6) =  0; m_B(4,7) = c3; m_B(4,8) = b3; m_B(4,9) =  0; m_B(4,10) = c4; m_B(4,11) = b4;
            m_B(5,0) = c1; m_B(5,1) =  0; m_B(5,2) = a1; m_B(5,3) = c2; m_B(5,4) =  0; m_B(5,5) = a2; m_B(5,6) = c3; m_B(5,7) =  0; m_B(5,8) = a3; m_B(5,9) = c4; m_B(5,10) =  0; m_B(5,11) = a4;

            //E from AFEM-16 p.14
            Mat6x6r E( Mat6x6r::Zero() );
            Real nu( m_rParams.m_PoissonRatio );
            Real k1( Real(1) - nu );
            Real k2( Real(0.5) - nu );
            E(0,0) = k1; E(0,1) = nu; E(0,2) = nu;
            E(1,0) = nu; E(1,1) = k1; E(1,2) = nu;
            E(2,0) = nu; E(2,1) = nu; E(2,2) = k1;
            E(3,3) = k2;
            E(4,4) = k2;
            E(5,5) = k2;
            E = ( m_rParams.m_YoungModulus / ((Real(1)+m_rParams.m_PoissonRatio)*(Real(1)-Real(2)*m_rParams.m_PoissonRatio)) ) * E;

            m_Ke = m_Volume * ( m_B.Transposed() * (E * m_B) );

            // Variable triangle node pos
#ifdef __TODO_PORT_TO_3D
            Mat2x2r R_beta( mal::GRotation2x2_From(m_rParams.m_Beta) );
#else
            Mat3x3r R_beta( Mat3x3r::Identity() );
#endif
            m_X[0] = R_beta * m_r[0];
            m_X[1] = R_beta * m_r[1];
            m_X[2] = R_beta * m_r[2];
            m_X[3] = R_beta * m_r[3];

            m_NoC = m_rParams.m_NoC;
        }

    Vec12r ComputeElasticForce( const Vec3r X[4], Params::EMethod method, Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            switch( method )
            {
            case Params::eMethod_Id: return ComputeElasticForce_Id( X, p_rotation, p_energy ); break;
            case Params::eMethod_QR_XY: return ComputeElasticForce_QR_XYZ( X, p_rotation, p_energy ); break;
            case Params::eMethod_QR_YX: return ComputeElasticForce_QR_GSL( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD: return ComputeElasticForce_PD( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Reflect: return ComputeElasticForce_PD_Reflect( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Fix: return ComputeElasticForce_PD_Fix( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Project: return ComputeElasticForce_PD_Project( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_Project_Nearest: return ComputeElasticForce_PD_Project_Nearest( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_SVD: return ComputeElasticForce_PD_SVD( X, p_rotation, p_energy ); break;
            case Params::eMethod_IHFSDM: return ComputeElasticForce_IHFSDM( X, p_rotation, p_energy ); break;
            case Params::eMethod_ITF_LRM: return ComputeElasticForce_ITF_LRM( X, p_rotation, p_energy ); break;
            case Params::eMethod_ECIE_CLRM: return ComputeElasticForce_ECIE_CLRM( X, p_rotation, p_energy ); break;
            case Params::eMethod_ECIE_NHC0: return ComputeElasticForce_ECIE_NHC0( X, p_rotation, p_energy ); break;
            case Params::eMethod_ECIE_NHC1: return ComputeElasticForce_ECIE_NHC1( X, p_rotation, p_energy ); break;
            case Params::eMethod_PD_CLRM: return ComputeElasticForce_PD_CLRM( X, p_rotation, p_energy ); break;
            default: return Vec12r::Zero();
            }
        }

    Vec3r ComputeDiagonalP( const Vec3r &vec_diag_F, Params::EMethod method ) const
        {
            switch( method )
            {
            case Params::eMethod_Id: return Vec3r::Zero(); break;
            case Params::eMethod_QR_XY: return Vec3r::Zero(); break;
            case Params::eMethod_QR_YX: return Vec3r::Zero(); break;
            case Params::eMethod_PD: return Vec3r::Zero(); break;
            case Params::eMethod_PD_Reflect: return Vec3r::Zero(); break;
            case Params::eMethod_PD_Fix: return Vec3r::Zero(); break;
            case Params::eMethod_PD_Project: return Vec3r::Zero(); break;
            case Params::eMethod_PD_Project_Nearest: return Vec3r::Zero(); break;
            case Params::eMethod_PD_SVD: return Vec3r::Zero(); break;
            case Params::eMethod_IHFSDM: return Vec3r::Zero(); break;
#ifdef __TODO_PORT_TO_3D
            case Params::eMethod_ITF_LRM: return ComputeDiagonalP_ITF_LRM( vec_diag_F ); break;
            case Params::eMethod_ECIE_CLRM: return ComputeDiagonalP_ECIE_CLRM( vec_diag_F ); break;
            case Params::eMethod_ECIE_NHC0: return ComputeDiagonalP_ECIE_NHC0( vec_diag_F ); break;
            case Params::eMethod_ECIE_NHC1: return ComputeDiagonalP_ECIE_NHC1( vec_diag_F ); break;
            case Params::eMethod_PD_CLRM: return Vec3r::Zero(); break;
#endif
            default: return Vec3r::Zero();
            }
        }

private:

    Mat3x3r Compute_F( const Vec3r X[4], uint32 ffm ) const
        {
            Mat3x3r Ds( mal::GMat3x3_From_Columns( X[1] - X[0], X[2] - X[0], X[3] - X[0] ) );
            Mat3x3r F( Ds * m_invDm );
            Real det_F( mal::Det(F) );
#ifdef __TODO_PORT_TO_3D
            if( det_F <= m_DegenerateThresholdDetF )
            {
                switch( ffm )
                {
                case eFFM_Project:
                    {
                        Vec3r axis12( mal::PerpendicularCW( X[2] - X[1] ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
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
                            Real delta_div_dist12( m_DegenerateThresholdDetF * 2*m_Area / dist12_sq ); //== delta / dist12
                            Vec3r projected_x0( X[0] );
                            projected_x0 -= ( mal::Dot(X[0]-X[1],axis12)/dist12_sq ) * axis12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
                            projected_x0 += delta_div_dist12*axis12; //Move past collapse

                            /*\todo Old SLOW code that requires an SQRT
                              Real dist12( mal::Norm( axis12 ) );
                              axis12 /= dist12;
                              Real delta( m_DegenerateThresholdDetF * 2*m_Area / dist12 );
                              Vec2r projected_x0( X[0] - mal::Dot(X[0]-X[1],axis12)*axis12 ); //Project onto collapse line
                              projected_x0 += delta*axis12; //Move past collapse
                            */

                            Ds = mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 );
                            F = Ds * m_invDm;
                            det_F = mal::Det(F);

                            if( det_F < 0.99f*m_rParams.m_DegenerateThresholdDetF || det_F > 1.01f*m_rParams.m_DegenerateThresholdDetF )
                                APP_LOG_WARNING( "Compute_F: det(F) %f != %f should not happen if delta is correct", det_F, m_DegenerateThresholdDetF );
                        }
                        else
                        {
                            F = Mat2x2r::Identity();
                            det_F = 1;
                        }
                    }
                    break;
                case eFFM_Project_Nearest:
                    {
                        F = Mat3x3r::Identity();
                        det_F = 1;
                    }
                    break;
                case eFFM_Reflect:
                    {
                        if( det_F < 0 )
                        {
                            Vec2r axis12( mal::PerpendicularCW( X[2] - X[1] ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                            F = mal::GReflection2x2_From_Axis( mal::Normalized( axis12 ) ) * F;
                            det_F = -det_F;
                        }
                    }
                    break;
                case eFFM_None:
                default:
                    break;
                }
            }
            APP_ASSERT( eFFM_None == ffm || mal::Det(F) >= 0 );
#endif
            APP_ASSERT( !mal::IsNaN(F) );
            return F;
        }

    Real Compute_Effective_PoissonRatio( Real det_F ) const
        {
            return m_rParams.m_PoissonRatio;
        }

    Mat12x12r Compute_Effective_Ke( Real det_F ) const
        {
            return m_Ke;
        }

    Vec12r Compute_Corotated_Forces( const Vec3r X[4], const Mat3x3r &R, Real det_F, Real *p_energy = 0 ) const
        {
            Vec12r u;
            mal::GSetRange<0,2>( u, R.Transposed()*X[0] - m_r[0] );
            mal::GSetRange<3,5>( u, R.Transposed()*X[1] - m_r[1] );
            mal::GSetRange<6,8>( u, R.Transposed()*X[2] - m_r[2] );
            mal::GSetRange<9,11>( u, R.Transposed()*X[3] - m_r[3] );
            Mat12x12r Re(Real(0));
            mal::GSetRange<0,0,2,2>( Re, R );
            mal::GSetRange<3,3,5,5>( Re, R );
            mal::GSetRange<6,6,8,8>( Re, R );
            mal::GSetRange<9,9,11,11>( Re, R );
            Mat12x12r Ke( Compute_Effective_Ke(det_F) );
            if( p_energy ) *p_energy = Real(0.5)*(u * (Ke * u));
            return -Re * Ke * u;
        }

#ifdef __TODO_PORT_TO_3D
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
#endif

    Vec12r ComputeElasticForce_Id( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Compute_Corotated_Forces( X, Mat3x3r::Identity(), Real(1), p_energy );
        }

    Vec12r ComputeElasticForce_QR_XYZ( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            /*TEMP Transposed test to match NPF05 usage on J=F^T, though PO09 uses it directly on F... but the results are not the same as transposing R afterwards!!
            Mat3x3r R( mal::GRotation3x3_GramSchmidtOrthonormalization_XYZ( mal::Transposed(F), Mat3x3r::Identity() ) );
            R = mal::Transposed(R);
            */
            Mat3x3r R( mal::GRotation3x3_GramSchmidtOrthonormalization_XYZ( F, Mat3x3r::Identity() ) );
            //Mat3x3r R( mal::GRotation3x3_GramSchmidtOrthonormalization_XYZ_GSO( F, Mat3x3r::Identity() ) ); //\todo This does NOT guarantee right-handed R
            APP_ASSERT( !mal::IsNaN(R) );
            APP_ASSERT( mal::Det(R) > 0 );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    Vec12r ComputeElasticForce_QR_GSL( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            Mat3x3r R( mal::GRotation3x3_QR( F ) );
            APP_ASSERT( !mal::IsNaN(R) );
            //\todo THIS FAILS FOR GSL-QR... APP_ASSERT( mal::Det(R) > 0 );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    Vec12r ComputeElasticForce_PD( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat3x3r R;
            if( mal::Abs(det_F) > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation3x3_PolarDecomposition( F, det_F );
            else R = mal::GRotation3x3_PolarDecomposition_From_SVD( F, det_F ); //\todo This has a glitch at det == 0
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec12r ComputeElasticForce_PD_Reflect( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat3x3r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation3x3_PolarDecomposition( F, det_F );
            else if( det_F >= 0 ) R = mal::GRotation3x3_PolarDecomposition_From_SVD( F, det_F ); //\todo This has a glitch at det == 0
            else // det_F < 0
            {
                Vec3r axis123( -mal::Cross( X[2] - X[1], X[3] - X[1] ) );
                F = mal::GReflection3x3_From_Axis( mal::Normalized( axis123 ) ) * F;
                if( det_F < -m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation3x3_PolarDecomposition( F, -det_F );
                else R = mal::GRotation3x3_PolarDecomposition_From_SVD( F, -det_F );
            }
            APP_LOG_ASSERT( mal::Det(R) > 0, "det(F,R) = %f, %f", det_F, mal::Det(R) ); //\todo This HAPPENS with GSL SVD, changed GRotation3x3_PolarDecomposition_From_SVD() to fix it
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec12r ComputeElasticForce_PD_Fix( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat2x2r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation2x2_PolarDecomposition( F, det_F );
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
                Vec2r x12( X[2] - X[1] );
                Real dist12_sq( mal::NormSq(x12) );
                if( dist12_sq > 0.000001f )
                {
                    // Compute rotation from r12 to x12, which is the same as global R when we rotate the resting state R=Id rigidly with x12
                    Vec2r r12( m_r[2]-m_r[1] );
                    //\todo In 3D, rotR would be similarly computed from the NORMAL of the crossed face
                    Mat2x2r rotR( mal::GRotation2x2_Vec2Vec( mal::Normalized(r12) , mal::Normalized(x12) ) ); //\todo Normalized(r12) could be precomputed
//#define __ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
#ifdef __ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
                    if( det_F > 0 )
                    {
                        Vec2r axis12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
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

                        // Compute Proj(F)
                        Mat2x2r projQ( mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 ) );
                        Mat2x2r projF( projQ * m_invDm );
                        Real det_projF( mal::Det(projF) );
                        APP_ASSERT( det_projF >= 0 );
                        if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                        Mat2x2r projR = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                        //R = RLerp( rotR, projR, det_F / m_rParams.m_DegenerateThresholdDetF );
                        R = RCerp( rotR, projR, det_F / m_rParams.m_DegenerateThresholdDetF );
                    }
                    else //det_F <= 0
                        R = rotR;
#else //__ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
                    if( det_F > m_ThresholdIpolDetF )
                    {
                        Vec2r axis12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
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

                        // Compute Proj(F)
                        Mat2x2r projQ( mal::GMat2x2_From_Columns( X[1] - projected_x0, X[2] - projected_x0 ) );
                        Mat2x2r projF( projQ * m_invDm );
                        Real det_projF( mal::Det(projF) );
                        APP_ASSERT( det_projF >= 0 );
                        if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                        Mat2x2r projR = mal::GRotation2x2_PolarDecomposition( projF, det_projF );
                        Real lambda01( (det_F - m_ThresholdIpolDetF) / (m_DegenerateThresholdDetF - m_ThresholdIpolDetF) );
                        R = RCerp( rotR, projR, lambda01 );
                    }
                    else //det_F <= 0
                        R = rotR;
#endif //__ENABLE_PD_FIX_POSITIVE_DET_F_ONLY
                }
                else
                    R = Mat2x2r::Identity();
            }

            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
        }

    Vec12r ComputeElasticForce_PD_Project( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            Real det_F( mal::Det(F) );
            Mat3x3r R;
            if( det_F > m_rParams.m_DegenerateThresholdDetF ) R = mal::GRotation3x3_PolarDecomposition( F, det_F );
            else // det_F <= m_rParams.m_DegenerateThresholdDetF
            {
                Vec3r axis123( -mal::Cross( X[2] - X[1], X[3] - X[1] ) ); //\note axis123 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
                Real base123( mal::Norm( axis123 ) ); //base123 is the area of the parallelogram, 2x area or the triangle
                if( base123 > 0.0001f )
                {
                    axis123 /= base123;
                    /* Compute delta so that projected x0 yields an element with det(F) = DTDF
                       det(F) = det(Q) / det(P) = DTDF
                       => det(Ds) = DTDF * det(Dm)
                       ( det(Ds) = 6*Volume(Ds) = delta * base(x1,x2,x3) ) //Volume of a prism with base base(x1,x2,x3) and height delta
                       => delta * base(x1,x2,x3) = DTDF * det(Dm)
                       => delta = DTDF * det(Dm) / base(x1,x2,x3)
                       => delta = DTDF * 6*Volume(Dm) / base(x1,x2,x3)
                    */
                    Real delta( m_rParams.m_DegenerateThresholdDetF * 6 * m_Volume / base123 ); //\todo It seems to be correct, the APP_LOG_WARNING below never appears...
                    Vec3r projected_x0( X[0] - mal::Dot(X[0]-X[1],axis123)*axis123 ); //Project onto collapse line
                    projected_x0 += delta*axis123; //Move past collapse, ensuring det_projF >= m_rParams.m_DegenerateThresholdDetF

                    Mat3x3r Ds( mal::GMat3x3_From_Columns( X[1] - projected_x0, X[2] - projected_x0, X[3] - projected_x0 ) );
                    Mat3x3r projF( Ds * m_invDm );
                    Real det_projF( mal::Det(projF) );
                    APP_ASSERT( det_projF >= 0 );
                    /* \todo THIS HAPPENED but now it should be fixed
                       if( det_projF < 0 )
                       {
                       APP_LOG_WARNING( "det(projF) %f < 0 , n12 (%f,%f) , (X[0]-x1)*n12 %f",
                       det_projF, axis12[0], axis12[1], mal::Dot(X[0]-X[1],axis12) );
                       }
                    */
                    if( det_projF < 0.99f*m_rParams.m_DegenerateThresholdDetF || det_projF > 1.01f*m_rParams.m_DegenerateThresholdDetF ) APP_LOG_WARNING( "det(F) %f != %f should not happen if delta is correct", det_projF, m_rParams.m_DegenerateThresholdDetF );
                    R = mal::GRotation3x3_PolarDecomposition( projF, det_projF );
                }
                else
                    R = Mat3x3r::Identity();
            }
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, det_F, p_energy );
        }

    Vec12r ComputeElasticForce_PD_Project_Nearest( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
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
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
        }

    Vec12r ComputeElasticForce_PD_SVD( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat3x3r U, Vt;
            Vec3r diag_S;
#ifdef __USE_FAST_SVD_PHYSBAM
            mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt ); //ALREADY enforces rotation, smallest singular value may be negative
#else
            mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
            // Force pure rotation U,Vt by fixing potential inversion/reflection see LeafDSH_FEM
            //if( mal::Det(U) * mal::Det(Vt) < 0 )
            {
                if( mal::Det(U) < 0 ) mal::GSetColumn<2>( U, -mal::GColumn<2>(U) );
                if( mal::Det(Vt) < 0 ) mal::GSetRow<2>( Vt, -mal::GRow<2>(Vt) );
            }
#endif
            Mat3x3r R( U*Vt );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    /* Implementation of the temporally-coherent heuristic eigenvalue
       negation idea in IHFSDM, assume all inversions are detected
       across (1,0) direction.
    */
    Vec12r ComputeElasticForce_IHFSDM( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
            Mat3x3r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat3x3r U, Vt;
            Vec3r diag_S;
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
                Vec3r n_c[4];
                n_c[0] = mal::Cross( m_r[3] - m_r[1], m_r[2] - m_r[1] );
                n_c[1] = mal::Cross( m_r[3] - m_r[0], m_r[2] - m_r[0] );
                n_c[2] = mal::Cross( m_r[3] - m_r[1], m_r[0] - m_r[1] );
                n_c[3] = mal::Cross( m_r[2] - m_r[1], m_r[0] - m_r[1] );
                Vec3r v[3];
                v[0] = mal::GRow<0>(Vt);
                v[1] = mal::GRow<1>(Vt);
                v[2] = mal::GRow<2>(Vt);
                Real min_lambda( mal::Infinity<Real>() );
                int min_it_v(0);
                //Real lambda_CxV[3][2];
                for( int it_c=0;
                     it_c < 4; //it_c<1; //\todo Uncomment to see the case where c = NoC = v0, clearer discontinuity
                     it_c++ )
                {
                    for( int it_v=0; it_v<3; it_v++ )
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
                switch( min_it_v )
                {
                case 0:
                    if( mal::Det(U) < 0 ) mal::GSetColumn<0>( U, -mal::GColumn<0>(U) );
                    if( mal::Det(Vt) < 0 ) mal::GSetRow<0>( Vt, -mal::GRow<0>(Vt) );
                    break;
                case 1:
                    if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
                    if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
                    break;
                case 2:
                    if( mal::Det(U) < 0 ) mal::GSetColumn<2>( U, -mal::GColumn<2>(U) );
                    if( mal::Det(Vt) < 0 ) mal::GSetRow<2>( Vt, -mal::GRow<2>(Vt) );
                    break;
                default:
                    APP_LOG_ERROR("IHFSDM::IMPOSSIBLE!");
                    break;
                }
#endif
            }
            Mat3x3r R( U*Vt );
            if( p_rotation ) *p_rotation = R;
            // Finally, compute forces
            return Compute_Corotated_Forces( X, R, mal::Det(F), p_energy );
        }

    /*\note Actually, we follow ITF usage, but compute f0 as in SIGG12 course notes
      Also, we base deformation at x0, not at x2 as SIGG12, thus, final f0 is computed symmetrically

     LRM: (ECIE pg 5, also ITF pg 5)
       Energy = \mu \sum_i (\theta_i - 1)^2 + \frac{\lambda}{2} (\sum_i (\theta_i - 1) )^2
    */
#ifdef __TODO_PORT_TO_3D
    Vec2r ComputeDiagonalP_ITF_LRM( const Vec2r &vec_diag_F ) const
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
            if( m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real trE( vec_diag_F[0] + vec_diag_F[1] - 2 );
            return Vec2r( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * trE,
                          2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * trE );
        }
#endif

    Vec12r ComputeElasticForce_ITF_LRM( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            // Compute SVD
            Mat2x2r U, Vt;
            Vec2r vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, vec_diag_F, Vt );

            /*TEMP: Clamp F to compute U*Vt uninverted
            F = Compute_F( X, eFFM_Project );
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
            }

            /*TEMP: Clamp F to compute U*Vt uninverted
            F = Compute_F( X, eFFM_Project );
            // Compute SVD
            Vec2r dummy_vec_diag_F;
            mal::GSingularValueDecomposition_USVt( F, U, dummy_vec_diag_F, Vt );
            //END TEMP
            */

            if( p_rotation ) *p_rotation = U*Vt;
            // Material = Linear Rotated Model == Corrotated (from ECIE definition)
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ITF_LRM(vec_diag_F) ) );
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
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
#ifdef __TODO_PORT_TO_3D
    Vec2r ComputeDiagonalP_ECIE_CLRM( const Vec2r &vec_diag_F ) const
        {
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real J( vec_diag_F[0] * vec_diag_F[1] );
            return Vec2r( 2 * lame_mu * (vec_diag_F[0] - 1) + lame_lambda * (J-1) * vec_diag_F[1],
                          2 * lame_mu * (vec_diag_F[1] - 1) + lame_lambda * (J-1) * vec_diag_F[0] );
        }
#endif

    Vec12r ComputeElasticForce_ECIE_CLRM( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
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
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ECIE_CLRM(vec_diag_F) ) );
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
        }

    Vec12r ComputeElasticForce_PD_CLRM( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
            Mat2x2r F( Compute_F( X, m_rParams.m_FFM ) );
            Mat2x2r fixedF( Compute_F( X, Element3D::eFFM_Project ) );
            //Mat2x2r fixedF( Compute_F( X, Element3D::eFFM_Reflect ) ); //TEMP: Reflect WORKS too, but requires PD_FROM_SVD in collapsed config, and has known SVD glitches
            Real det_F( mal::Det(F) );
            Real det_fixedF( mal::Det(fixedF) );
            Mat2x2r R( mal::GRotation2x2_PolarDecomposition( fixedF, det_fixedF ) );
            //Mat2x2r R( mal::GRotation2x2_PolarDecomposition_From_SVD( fixedF, det_fixedF ) ); //TEMP: Req by eFFM_Reflect, simple PD crashes here
            if( p_rotation ) *p_rotation = R;
            Mat2x2r S( R.Transposed() * F ); //F = R*S
            Real det_S( mal::Det(S) ); //\note == det_F
            Mat2x2r P( R * ( 2*m_LameMu*(S-Mat2x2r::Identity())
                             + m_LameLambda * (det_S-1) * Mat2x2r( S(1,1), -S(1,0),
                                                                   -S(0,1), S(0,0) ) ) );
            // from Compute_PiolaKirchhoff_Forces()
            Mat2x2r Be0( m_invDm ); //\todo Be0 = D_m^-1 = invP from the undeformed P, deformed Q notation... confusing
            Mat2x2r H( - m_Area * P * Be0.Transposed() ); // Negation is present in the Siggraph2012 course notes, but not in ITF...
            Vec2r f0( - mal::GColumn<0>(H) - mal::GColumn<1>(H) ); //h0 = -h1-h2;
            Vec6r f( f0.x(), f0.y(),
                     H(0,0), H(1,0),
                     H(0,1), H(1,1) );
            return f;
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
        }

    /* ECIE C0 extension to neo-hookean, clamping eigenvalues < ecie_e and J accordingly.
       NOT strictly demonstrated in ECIE, but suggested in other articles from the same autors
    */
#ifdef __TODO_PORT_TO_3D
    Vec2r ComputeDiagonalP_ECIE_NHC0( const Vec2r &vec_diag_F ) const
        {
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real ecie_e( m_ECIE_e_threshold );
            Vec2r vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                                      mal::Max( ecie_e, vec_diag_F[1] ) );
            Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
            return Vec2r( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                          lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
        }
#endif
    Vec12r ComputeElasticForce_ECIE_NHC0( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
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
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ECIE_NHC0(vec_diag_F) ) );
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
        }

    /* ECIE C1 extension to neo-hookean
       \todo Parameters inversion_threshold and k are hardcoded to
       ECIE values, but could be published...
       Derivatives and notation from ECIE Supplementary Technical Document
    */
#ifdef __TODO_PORT_TO_3D
    Vec2r ComputeDiagonalP_ECIE_NHC1( const Vec2r &vec_diag_F ) const
        {
            Real lame_mu(m_LameMu);
            Real lame_lambda(m_LameLambda);
            Real det_F( vec_diag_F[0] * vec_diag_F[1] );
            if( m_InvertedCompressibilityFactorDetF > 0 && det_F < 0 )
                LameParameters_From_YoungAndPoisson( m_YoungModulus, Compute_Effective_PoissonRatio(det_F),
                                                     lame_mu, lame_lambda );
            Real ecie_e( m_ECIE_e_threshold );
            Real ecie_k( m_ECIE_k_factor * m_YoungModulus );
            bool bIsDegenerate0( vec_diag_F[0] < ecie_e );
            bool bIsDegenerate1( vec_diag_F[1] < ecie_e );
            Vec2r vec_clamped_diag_F( mal::Max( ecie_e, vec_diag_F[0] ),
                                      mal::Max( ecie_e, vec_diag_F[1] ) );
            Real clamped_J( vec_clamped_diag_F[0] * vec_clamped_diag_F[1] );
            Vec2r vec_diag_P( lame_mu * vec_clamped_diag_F[0] + (vec_clamped_diag_F[1] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ),
                              lame_mu * vec_clamped_diag_F[1] + (vec_clamped_diag_F[0] / clamped_J) * ( lame_lambda * mal::Log(clamped_J) - lame_mu ) );
            Real h01( lame_lambda / clamped_J );
            //\todo Cases can be generalized if delta0,delta1 are 0 for non-clamped axis
            if( bIsDegenerate0 && bIsDegenerate1 )
            {
                Real delta0( vec_diag_F[0] - ecie_e );
                Real delta1( vec_diag_F[1] - ecie_e );
                vec_diag_P[0] += h01*delta1 + 2*ecie_k*delta0;
                vec_diag_P[1] += h01*delta0 + 2*ecie_k*delta1;
            }
            else if( bIsDegenerate0 && !bIsDegenerate1 )
            {
                Real delta0( vec_diag_F[0] - ecie_e );
                vec_diag_P[0] += 2*ecie_k*delta0;
                vec_diag_P[1] += h01*delta0;
            }
            else if( !bIsDegenerate0 && bIsDegenerate1 )
            {
                Real delta1( vec_diag_F[1] - ecie_e );
                vec_diag_P[0] += h01*delta1;
                vec_diag_P[1] += 2*ecie_k*delta1;
            }
            return vec_diag_P;
        }
#endif
    Vec12r ComputeElasticForce_ECIE_NHC1( const Vec3r X[4], Mat3x3r *p_rotation = 0, Real *p_energy = 0 ) const
        {
#ifdef __TODO_PORT_TO_3D
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
            Mat2x2r diag_P( mal::GMatNxN_From_Diagonal( ComputeDiagonalP_ECIE_NHC1(vec_diag_F) ) );
            // Compute forces from diag_P
            return Compute_PiolaKirchhoff_Forces( diag_P, U, Vt );
#else
            if( p_rotation ) *p_rotation = Mat3x3r::Identity();
            return Vec12r::Zero();
#endif
        }

private:
    friend class AppRenderer;
    friend class AppTask;
    // Material params
    Real m_LameMu;
    Real m_LameLambda;
    // Constant linear triangle stuff
    Vec3r m_r[4];
    Mat3x3r m_invDm;
    Mat6x12r m_B;
    Mat12x12r m_Ke;
    Real m_Volume;
    // Variable triangle node pos
    Vec3r m_X[4];
    int32 m_NoC; //Current NodeOfCollapse, -1 if none
};

struct DoC //Possible cases are VF => vit = 00,11,22,33 (vit0==vit1) and EE => vit = 01,02,03,12,13,23 (vit0<vit1)
{
    inline DoC() : m_VIT0(3), m_VIT1(0) {} //Invalid
    inline explicit DoC( int vit ) : m_VIT0(vit), m_VIT1(vit) {} //V-F
    inline DoC( int vit0, int vit1 ) : m_VIT0(vit0), m_VIT1(vit1) {} //E-E
    inline bool IsValid() const { return m_VIT0 <= m_VIT1; }
    inline bool IsInvalid() const { return m_VIT0 > m_VIT1; }
    inline bool IsVF() const { return m_VIT0 == m_VIT1; }
    inline bool IsEE() const { return m_VIT0 < m_VIT1; }
    uint8 m_VIT0, m_VIT1; //\todo Actually, 2 bits per vit would be enough (0..3)
};

/*! Compute Node-of-Collapse
  \note Internally computes ToC
  \todo Should compute Direction-of-Collapse (F1,F2) in 3D
*/
DoC ComputeDoC( const Vec3r &a0, const Vec3r &a1, const Vec3r &a2, const Vec3r &a3,
                const Vec3r &b0, const Vec3r &b1, const Vec3r &b2, const Vec3r &b3,
                Real volume, Real degenerate_threshold_det_F )
{
    DoC doc; //Invalid by default
    // Gather node displacements
    Vec3r d0( b0 - a0 );
    Vec3r d1( b1 - a1 );
    Vec3r d2( b2 - a2 );
    Vec3r d3( b3 - a3 );
    // Compute displacement differences
    Vec3r d10( d1 - d0 );
    Vec3r d20( d2 - d0 );
    Vec3r d30( d3 - d0 );
    // Compute position differences
    Vec3r a10( a1 - a0 );
    Vec3r a20( a2 - a0 );
    Vec3r a30( a3 - a0 );
    // Compute ToC, when det(F) = m_DegenerateThresholdDetF => det(Ds) = m_DegenerateThresholdDetF * det(Dm) => det(Ds) = m_DegenerateThresholdDetF * 6 * Volume(Dm)
    Real K( degenerate_threshold_det_F * Real(6) * volume );
    Real A( + d10.x() * d20.y() * d30.z()
            + d20.x() * d30.y() * d10.z()
            + d30.x() * d10.y() * d20.z()
            - d30.x() * d20.y() * d10.z()
            - d20.x() * d10.y() * d30.z()
            - d10.x() * d30.y() * d20.z() );

    Real B( + ( a10.x() * d20.y() * d30.z() + d10.x() * a20.y() * d30.z() + d10.x() * d20.y() * a30.z() )
            + ( a20.x() * d30.y() * d10.z() + d20.x() * a30.y() * d10.z() + d20.x() * d30.y() * a10.z() )
            + ( a30.x() * d10.y() * d20.z() + d30.x() * a10.y() * d20.z() + d30.x() * d10.y() * a20.z() )
            - ( a30.x() * d20.y() * d10.z() + d30.x() * a20.y() * d10.z() + d30.x() * d20.y() * a10.z() )
            - ( a20.x() * d10.y() * d30.z() + d20.x() * a10.y() * d30.z() + d20.x() * d10.y() * a30.z() )
            - ( a10.x() * d30.y() * d20.z() + d10.x() * a30.y() * d20.z() + d10.x() * d30.y() * a20.z() ) );

    Real C( + ( a10.x() * a20.y() * d30.z() + a10.x() * d20.y() * a30.z() + d10.x() * a20.y() * a30.z() )
            + ( a20.x() * a30.y() * d10.z() + a20.x() * d30.y() * a10.z() + d20.x() * a30.y() * a10.z() )
            + ( a30.x() * a10.y() * d20.z() + a30.x() * d10.y() * a20.z() + d30.x() * a10.y() * a20.z() )
            - ( a30.x() * a20.y() * d10.z() + a30.x() * d20.y() * a10.z() + d30.x() * a20.y() * a10.z() )
            - ( a20.x() * a10.y() * d30.z() + a20.x() * d10.y() * a30.z() + d20.x() * a10.y() * a30.z() )
            - ( a10.x() * a30.y() * d20.z() + a10.x() * d30.y() * a20.z() + d10.x() * a30.y() * a20.z() ) );

    Real D( + a10.x() * a20.y() * a30.z()
            + a20.x() * a30.y() * a10.z()
            + a30.x() * a10.y() * a20.z()
            - a30.x() * a20.y() * a10.z()
            - a20.x() * a10.y() * a30.z()
            - a10.x() * a30.y() * a20.z()
            - K );

    Real vec_toc[3];
    int num_roots( mal::GSolvePolynomialEq3<Real>( A, B, C, D, vec_toc[0], vec_toc[1], vec_toc[2] ) );
    if( num_roots > 0 )
    {
        Real toc(2); //no toc > 1 is accepted
        int num_roots01(0);
        for( int i=0; i<num_roots; i++ )
        {
            if( vec_toc[i] > Real(0) && vec_toc[i] <= Real(1) && vec_toc[i] < toc )
            {
                toc = vec_toc[i];
                num_roots01++;
            }
        }
        /*
        MS_LOG( "Potentially degenerated with toc = %f, %f, %f, #roots = %d, #roots01 = %d, min_toc = %f",
                vec_toc[0], vec_toc[1], vec_toc[2], num_roots, num_roots01, toc );
        */
        if( num_roots01 != 1 && num_roots01 != 3 ) return doc; //TEMPORAL TO TEST
        APP_ASSERT( num_roots01 == 1 || num_roots01 == 3 );
        APP_ASSERT( Real(0) < toc && toc <= Real(1) );
        //\todo We knot that det(F(0)) > DTDF, and that det(F(1)) <= DTDF, so there MUST BE STRICTLY 1 or 3 crossing in toc = [0,1]
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
            // Compute DoC
            Vec3r q0( a0 + toc*d0 );
            Vec3r q1( a1 + toc*d1 );
            Vec3r q2( a2 + toc*d2 );
            Vec3r q3( a3 + toc*d3 );
            // Compute collapse normal using the smallest singular-value from SVD
            /*\todo THIS SEEMS TO WORK FINE, but I'm not sure if it's
              justified... the matrix is not exaclty the
              covariance... */
            Mat3x3r Q( mal::GMat3x3_From_Columns( q1-q0, q2-q0, q3-q0 ) );
            Mat3x3r U,Vt;
            Vec3r S;
            mal::GSingularValueDecomposition_USVt( Q, U, S, Vt );
            Vec3r N = mal::GRow<2>( Vt ); //SVel from the smallest SVal
            /* \todo non-unitary normal computed by adding "some"
               non-unitary face normals SHOULD also work, and would
               avoid expensive SVD...
            for( )
            {
                if( mal::Dot(edge_ij,N) < 0 )
                    N -= cross;
                else
                    N += cross;
            }
            */
            // Compute face areas (actually, per-face unit-height prism volumes, which is the same but avoids sqrt()
            Real vec_area[4]; //index is actually the VertexInTetrahedron of the vertex opposite to the actual face
            vec_area[0] = mal::Abs( mal::Dot( N, mal::Cross( q2-q1, q3-q1 ) ) );
            vec_area[1] = mal::Abs( mal::Dot( N, mal::Cross( q2-q0, q3-q0 ) ) );
            vec_area[2] = mal::Abs( mal::Dot( N, mal::Cross( q1-q0, q3-q0 ) ) );
            vec_area[3] = mal::Abs( mal::Dot( N, mal::Cross( q1-q0, q2-q0 ) ) );
            Real sum_area( vec_area[0] + vec_area[1] + vec_area[2] + vec_area[3] );
            // Check V-F case by searching for a face with Area_i = 1/2 \sum_j Area_j
            //\todo Precision may be critical here, consider doing FIRST E-E tests and selecting the V-F case with F_i CLOSEST TO 1/2 \sum Area_i otherwise
            for( int it_fit=0; it_fit<4; it_fit++ )
                if( mal::ApproxEq( vec_area[it_fit], Real(0.5)*sum_area ) )
                    doc = DoC( it_fit );
            // If not found, we have an E-E case
            if( doc.IsInvalid() )
            {
                // Determine which outgoing edge from V_0 crosses another one on the opposite face F_0
                if( mal::Dot( mal::Cross(q1-q0, q2-q0), mal::Cross(q1-q0, q3-q0) ) < Real(0) ) //i,j,k,l = 0,1,2,3
                    doc = DoC(0,1);
                else if( mal::Dot( mal::Cross(q2-q0, q1-q0), mal::Cross(q2-q0, q3-q0) ) < Real(0) ) //i,j,k,l = 0,2,1,3
                    doc = DoC(0,2);
                else if( mal::Dot( mal::Cross(q3-q0, q1-q0), mal::Cross(q3-q0, q2-q0) ) < Real(0) ) //i,j,k,l = 0,3,1,2
                    doc = DoC(0,3);
                else
                {
                    APP_LOG_ERROR("MISSED POTENTIAL E-E crossing, forcing DoC = V-F 0");
                    doc = DoC(0);
                }
            }
            else
            {
                //MS_LOG("V-F crossing DoC = %d", doc.m_VIT0 );
            }
        }
        /*TEMP
          else
          {
          MS_LOG_WARNING( "No 0 or 2 degenerations in range [0,1]" );
          }
        */
    }
    else
    {
        APP_LOG_ERROR( "Cannot solve A*t^3 + B*t^2 + C*t + D = 0 for A = %f, B = %f, C = %f, D = %f, #roots = %d", A, B, C, D, num_roots );
    }
    return doc;
}

#endif //TEST_FE_ELEMENT3D_H
