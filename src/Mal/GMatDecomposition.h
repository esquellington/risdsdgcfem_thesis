#ifndef MAL_GMAT_DECOMPOSITION_H
#define MAL_GMAT_DECOMPOSITION_H

#include <Mal/Config.h>
#include <Mal/GMat.h>
#include <Mal/GMatUtils.h>
#include <Mal/GMatEigen.h>
#include <Mal/GVec.h>
#include <Mal/GMatDecomposition_SVD3x3_PhysBAM.h> //\note Implementation of FROM_PhysBAM::GFastSVD3x3( F, U, diag_F, Vt );

namespace mal
{

/*! Polar decomposition to extract rotation R from H = R * S \see FEM.tex
  \note We allow det_H < 0, which will result in an incorrectly oriented R with det_R < 0.
*/
template <typename T>
inline GMat<T,2,2> GRotation2x2_PolarDecomposition( const GMat<T,2,2> &H, T det_H )
{
    MAL_ASSERT( mal::Abs(det_H) > mal::Epsilon<T>() ); //should be specific epsilon for actual ops required (sqrt, rcp, etc...)
    GMat<T,2,2> H2;
    H2(0,0) = H(1,1); H2(0,1) = -H(1,0);
    H2(1,0) = -H(0,1); H2(1,1) = H(0,0);
    GMat<T,2,2> R( H + mal::Sign(det_H) * H2 );
    T inv_norm_R0( mal::Rcp( mal::Sqrt( mal::Sq(R(0,0)) + mal::Sq(R(1,0)) ) ) );
    R(0,0) *= inv_norm_R0; R(1,0) *= inv_norm_R0;
    T inv_norm_R1( mal::Rcp( mal::Sqrt( mal::Sq(R(0,1)) + mal::Sq(R(1,1)) ) ) );
    R(0,1) *= inv_norm_R1; R(1,1) *= inv_norm_R1;
    return R;
}

/*! QR decomposition => Gram-Schmidt Orthonormalization to extract rotation R from H = R * U \see FEM.tex
  \note If all columns of H are collapsed the QR decomposition is undefined, and we return the given default_R
  \note We do not use standard GSO algorithm but rather optimize it for the 2D case
  \see http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
  //\todo XY or YX This should not be relevant, but Y-first simplifies GSO/QR glitch demonstration in testSaphyre
*/
template <typename T>
inline GMat<T,2,2> GRotation2x2_GramSchmidtOrthonormalization_XY( const GMat<T,2,2> &H, const GMat<T,2,2> &default_R = GMat<T,2,2>::Identity() )
{
    GVec<T,2> e0( mal::GColumn<0>(H) );
    GVec<T,2> e1( mal::GColumn<1>(H) );
    T norm_sq_e0( mal::NormSq( e0 ) );
    T norm_sq_e1( mal::NormSq( e1 ) );
    if( norm_sq_e0 > mal::Epsilon<T>() )
    {
        e0 /= mal::Sqrt( norm_sq_e0 );
        e1 = GVec<T,2>( -e0[1], e0[0] ); //e1 = PerpendicularCW(e0)
        return GMat2x2_From_Columns( e0, e1 );
        /*\note GSO would do this:
        e1 -= mal::Dot( e0, e1 ) * e0;
        norm_sq_e1 = mal::NormSq( e1 );
        if( norm_sq_e1 > mal::Epsilon<T>() )
            return GMat2x2_From_Columns( e0, e1 / mal::Sqrt(norm_sq_e1) );
        else
            return GMat2x2_From_Columns( e0, mal::PerpendicularCW(e0) );
        */
    }
    else if( norm_sq_e1 > mal::Epsilon<T>() )
    {
        e1 /= mal::Sqrt( norm_sq_e1 );
        e0 = GVec<T,2>( e1[1], -e1[0] ); //e0 = -PerpendicularCW(e1) = PerpendicularCCW(e1)
        return GMat2x2_From_Columns( e0, e1 );
    }
    else
        return default_R;
}

/*! QR decomposition => Gram-Schmidt Orthonormalization to extract rotation R from H = R * U \see FEM.tex
  \note If all columns of H are collapsed the QR decomposition is undefined, and we return the given default_R
  \todo Direct GSO implementation, to compare with 3D version
  \see http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
*/
template <typename T>
inline GMat<T,2,2> GRotation2x2_GramSchmidtOrthonormalization_XY_GSO( const GMat<T,2,2> &H, const GMat<T,2,2> &default_R = GMat<T,2,2>::Identity() )
{
    GVec<T,2> e0( mal::GColumn<0>(H) );
    GVec<T,2> e1( mal::GColumn<1>(H) );
    T norm_sq_e0( mal::NormSq( e0 ) );
    if( norm_sq_e0 > mal::Epsilon<T>() )
    {
        e0 /= Sqrt( norm_sq_e0 );
        e1 -= mal::Dot( e0, e1 ) * e0;
        T norm_sq_e1( mal::NormSq( e1 ) );
        if( norm_sq_e1 < mal::Epsilon<T>() ) return default_R;
        e1.Normalize();
        return GMat2x2_From_Columns( e0, e1 );
    }
    else
        return default_R;
}

/*! QR decomposition => Gram-Schmidt Orthonormalization to extract rotation R from H = R * U \see FEM.tex
  \note If all columns of H are collapsed the QR decomposition is undefined, and we return the given default_R
  \note We do not use standard GSO algorithm but rather optimize it for the 2D case
  \see http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
  //\todo XY or YX This should not be relevant, but Y-first simplifies GSO/QR glitch demonstration in testSaphyre
*/
template <typename T>
inline GMat<T,2,2> GRotation2x2_GramSchmidtOrthonormalization_YX( const GMat<T,2,2> &H, const GMat<T,2,2> &default_R = GMat<T,2,2>::Identity() )
{
    GVec<T,2> e0( mal::GColumn<0>(H) );
    GVec<T,2> e1( mal::GColumn<1>(H) );
    T norm_sq_e0( mal::NormSq( e0 ) );
    T norm_sq_e1( mal::NormSq( e1 ) );
    if( norm_sq_e1 > mal::Epsilon<T>() )
    {
        e1 /= mal::Sqrt( norm_sq_e1 );
        e0 = GVec<T,2>( e1[1], -e1[0] ); //e0 = -PerpendicularCW(e1) = PerpendicularCCW(e1)
        return GMat2x2_From_Columns( e0, e1 );
        /*\note GSO would do this:
        e1 -= mal::Dot( e0, e1 ) * e0;
        norm_sq_e1 = mal::NormSq( e1 );
        if( norm_sq_e1 > mal::Epsilon<T>() )
            return GMat2x2_From_Columns( e0, e1 / mal::Sqrt(norm_sq_e1) );
        else
            return GMat2x2_From_Columns( e0, mal::PerpendicularCW(e0) );
        */
    }
    else if( norm_sq_e0 > mal::Epsilon<T>() )
    {
        e0 /= mal::Sqrt( norm_sq_e0 );
        e1 = GVec<T,2>( -e0[1], e0[0] ); //e1 = PerpendicularCW(e0)
        return GMat2x2_From_Columns( e0, e1 );
    }
    else
        return default_R;
}

/*! QR decomposition => Gram-Schmidt Orthonormalization to extract rotation R from H = R * U \see FEM.tex
  \note If all columns of H are collapsed the QR decomposition is undefined, and we return the given default_R
  \note We do not use standard GSO algorithm but rather optimize it for the 2D case
  \see http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
  //\todo XY or YX This should not be relevant, but Y-first simplifies GSO/QR glitch demonstration in testSaphyre
*/
template <typename T>
inline GMat<T,3,3> GRotation3x3_GramSchmidtOrthonormalization_XYZ( const GMat<T,3,3> &H, const GMat<T,3,3> &default_R = GMat<T,3,3>::Identity() )
{
    GVec<T,3> e0( mal::GColumn<0>(H) );
    GVec<T,3> e1( mal::GColumn<1>(H) );
    GVec<T,3> e2( mal::GColumn<2>(H) );
    T norm_sq_e0( mal::NormSq( e0 ) );
    T norm_sq_e1( mal::NormSq( e1 ) );
    T norm_sq_e2( mal::NormSq( e2 ) );
    if( norm_sq_e0 > mal::Epsilon<T>() )
    {
        e0 /= Sqrt( norm_sq_e0 );
        e2 = Cross( e0, e1 );
        if( NormSq( e2 ) > mal::Epsilon<T>() )
        {
            e2.Normalize();
            e1 = Cross( e2, e0 );
            return GMat3x3_From_Columns( e0, e1, e2 );
        }
        else
        {
            e1 = Cross( e2, e0 );
            if( NormSq( e1 ) > mal::Epsilon<T>() )
            {
                e1.Normalize();
                e2 = Cross( e0, e1 );
                return GMat3x3_From_Columns( e0, e1, e2 );
            }
            else
                return default_R;
        }
    }
    else if( norm_sq_e1 > mal::Epsilon<T>() )
    {
        e1 /= Sqrt( norm_sq_e1 );
        e0 = Cross( e1, e2 );
        if( NormSq( e0 ) > mal::Epsilon<T>() )
        {
            e0.Normalize();
            e2 = Cross( e0, e1 );
            return GMat3x3_From_Columns( e0, e1, e2 );
        }
        else
        {
            e2 = Cross( e0, e1 );
            if( NormSq( e2 ) > mal::Epsilon<T>() )
            {
                e2.Normalize();
                e0 = Cross( e1, e2 );
                return GMat3x3_From_Columns( e0, e1, e2 );
            }
            else
                return default_R;
        }
    }
    else if( norm_sq_e2 > mal::Epsilon<T>() )
    {
        e2 /= Sqrt( norm_sq_e2 );
        e1 = Cross( e2, e0 );
        if( NormSq( e1 ) > mal::Epsilon<T>() )
        {
            e1.Normalize();
            e0 = Cross( e1, e2 );
            return GMat3x3_From_Columns( e0, e1, e2 );
        }
        else
        {
            e0 = Cross( e1, e2 );
            if( NormSq( e0 ) > mal::Epsilon<T>() )
            {
                e0.Normalize();
                e1 = Cross( e2, e0 );
                return GMat3x3_From_Columns( e0, e1, e2 );
            }
            else
                return default_R;
        }
    }
    else
        return default_R;
}

/*! QR decomposition => Gram-Schmidt Orthonormalization to extract rotation R from H = R * U \see FEM.tex
  \note If all columns of H are collapsed the QR decomposition is undefined, and we return the given default_R
  \see http://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
*/
template <typename T>
inline GMat<T,3,3> GRotation3x3_GramSchmidtOrthonormalization_XYZ_GSO( const GMat<T,3,3> &H, const GMat<T,3,3> &default_R = GMat<T,3,3>::Identity() )
{
    GVec<T,3> e0( mal::GColumn<0>(H) );
    GVec<T,3> e1( mal::GColumn<1>(H) );
    GVec<T,3> e2( mal::GColumn<2>(H) );
    T norm_sq_e0( mal::NormSq( e0 ) );
    if( norm_sq_e0 > mal::Epsilon<T>() )
    {
        e0 /= Sqrt( norm_sq_e0 );
        e1 -= mal::Dot( e0, e1 ) * e0;
        T norm_sq_e1( mal::NormSq( e1 ) );
        if( norm_sq_e1 < mal::Epsilon<T>() ) return default_R;
        e1.Normalize();
        /* Standard GSO does NOT ensure right-handed R
        e2 -= mal::Dot( e0, e2 ) * e0;
        e2 -= mal::Dot( e1, e2 ) * e1;
        T norm_sq_e2( mal::NormSq( e2 ) );
        if( norm_sq_e2 < mal::Epsilon<T>() ) return default_R;
        e2.Normalize();
        */
        // Ad-hoc cross-product to ensure righ-handed R => Produces SAME R as GRotation3x3_GramSchmidtOrthonormalization_XYZ()
        e2 = mal::Normalized( mal::Cross(e0,e1) );
        return GMat3x3_From_Columns( e0, e1, e2 );
    }
    else
        return default_R;
}

/* Compute the SVD of a real-valued 2x2 matrix A = U * S * V^t
   \note Written according to ITF article to reproduce its behaviour
   \todo \post diag_S >= 0 ??
*/
template <typename T>
inline void GSingularValueDecomposition_USVt( const GMat<T,2,2> &F,
                                              GMat<T,2,2> &U, GVec<T,2> &diag_F, GMat<T,2,2> &Vt )
{
    MAL_ASSERT( !mal::IsNaN( F ) );

    // Compute sorted eigenvalues of A*A^t, l1 >= l2
    GMat<T,2,2> FtF( F.Transposed() * F );
    GVec<T,2> eigen_values = ComputeEigenvalues_Symmetric( FtF );
    if( IsNaN(eigen_values) )
    {
        //\note This has been seen to occur when F +=- Identity, and both eigen_values are -nan, while they should be (1,1)
        T det_F( mal::Det(F) );
        if( mal::Abs(det_F - T(1)) < 0.05 ) //\todo error > 0.02 observed, which is ENORMOUS...
        {
            /*\todo Silenced, but happens
            MAL_LOG_WARNING("NaN! det(F) = %f +=- 1 => eigen_values = (1,1) from F = [%f,%f,%f,%f]",
                            det_F,
                            F(0,0), F(0,1), F(1,0), F(1,1) );
            */
            U = GMat<T,2,2>::Identity();
            Vt = GMat<T,2,2>::Identity();
            diag_F = GVec<T,2>(1);
            return;
        }
        else
        {
            MAL_LOG_ERROR("NaN! det(F) = %f != 1 but eigen_values = (%f,%f)... unexpected NaN from F = [%f,%f,%f,%f]",
                          det_F,
                          eigen_values[0], eigen_values[1],
                          F(0,0), F(0,1), F(1,0), F(1,1) );
            MAL_ASSERT(false);
        }
    }
    //MAL_ASSERT( !IsNaN(eigen_values) );
    diag_F[0] = Sqrt( eigen_values[0] );
    diag_F[1] = Sqrt( eigen_values[1] );
    if( IsNaN(diag_F) )
    {
        //\note This has been seen to occur when some eigenvalue = -0.000000 ...which looks like a truncation error in ev computation that causes Sqrt(ev<0) = NaN
        MAL_LOG_WARNING("NaN! eigen_values = %f, %f, F = [%f,%f,%f,%f]",
                        eigen_values[0], eigen_values[1],
                        F(0,0), F(0,1), F(1,0), F(1,1) );
        if( eigen_values[0] < 0 ) diag_F[0] = 0;
        if( eigen_values[1] < 0 ) diag_F[1] = 0;
        MAL_ASSERT( !IsNaN(diag_F) );
    }
    //MAL_ASSERT( !IsNaN(diag_F) );
    bool bIsZeroF0( IsZero(diag_F[0]) );
    bool bIsZeroF1( IsZero(diag_F[1]) );

    /* TEMP: Negate smallest singular value if F inverted: This seems to work but does NOT fix discontinuity when 1 sval is near 0
    if( Det(F) < 0 ) diag_F[1] = -diag_F[1];
    */

    // Compute eigenvectors => U
    GMat<T,2,2> V = ComputeEigenbasisFromEigenvalues( FtF, eigen_values );
    //MAL_ASSERT( !IsNaN(V) );
    GMat<T,2,2> pseudo_inv_F( bIsZeroF0 ? 0 : Rcp(diag_F[0]), 0,
                              0, bIsZeroF1 ? 0 : Rcp(diag_F[1]) );

    //MAL_ASSERT( !IsNaN(pseudo_inv_F) );
    U = F * V * pseudo_inv_F;
    Vt = V.Transposed();

    // Consider degenerate cases
    if( bIsZeroF0 && bIsZeroF1 )
    {
        // \todo diag_F == 0 => F == 0,  U == Vt == Id may be arbitrary, consider forcing R = U*Vt = LastCorrect, maybe returning "false"
        U = GMat<T,2,2>::Identity();
        Vt = GMat<T,2,2>::Identity();
        diag_F = GVec<T,2>::Zero(); //force strict 0
    }
    else if( bIsZeroF0 && !bIsZeroF1 )
    {
        // diag_F[0] == 0 XOR diag_F[1] == 0 => Compute null eigenvector as perpendicular to the other one
        GSetRow<0>( Vt, PerpendicularCCW( GRow<1>(Vt) ) );
        GSetColumn<0>( U, PerpendicularCCW( GColumn<1>(U) ) );
        diag_F[0] = T(0); //force strict 0
    }
    else if( !bIsZeroF0 && bIsZeroF1 )
    {
        GSetRow<1>( Vt, PerpendicularCCW( GRow<0>(Vt) ) );
        GSetColumn<1>( U, PerpendicularCCW( GColumn<0>(U) ) );
        diag_F[1] = T(0); //force strict 0
    }
    //\todo diag_F[0] == diag_F[1] => Indeterminate eigenvectors?
    MAL_ASSERT( !mal::IsNaN( U ) );
    MAL_ASSERT( !mal::IsNaN( diag_F ) );
    MAL_ASSERT( !mal::IsNaN( Vt ) );
}

/* Compute the SVD of a real-valued 3x3 matrix A = U * S * V^t
   \note Written according to ITF article to reproduce its behaviour
   \todo \post diag_S >= 0 ??
*/
template <typename T>
inline void GSingularValueDecomposition_USVt( const GMat<T,3,3> &F,
                                              GMat<T,3,3> &U, GVec<T,3> &diag_F, GMat<T,3,3> &Vt )
{
    MAL_ASSERT( !mal::IsNaN( F ) );
    FROM_PhysBAM::GFastSVD3x3( F, U, diag_F, Vt );
    //\todo diag_F[0] == diag_F[1] => Indeterminate eigenvectors?
    MAL_ASSERT( !mal::IsNaN( U ) );
    MAL_ASSERT( !mal::IsNaN( diag_F ) );
    MAL_ASSERT( !mal::IsNaN( Vt ) );
}

//TEMP: Polar Decomposition rotation from SVD
template <typename T>
inline GMat<T,2,2> GRotation2x2_PolarDecomposition_From_SVD( const GMat<T,2,2> &H, T det_H )
{
    MAL_ASSERT( !mal::IsNaN( H ) );

    //TEMP: Otherwise, both eigenvectors are null, for some reason... probably eigenvalues are the same, happens when H==Id
    if( mal::ApproxEq( det_H, T(1) ) ) return GRotation2x2_PolarDecomposition(H,det_H);

    // Compute SVD of H = U * S * V^T
    GMat<T,2,2> U, Vt;
    GVec<T,2> diag_S;
    GSingularValueDecomposition_USVt( H, U, diag_S, Vt );
    //TEMP if( mal::ApproxEq( det_H, T(1) ) ) std::cout << "det(H) == 1 => diag_S =" << diag_S << std::endl;
    //std::cout << "PDSVD of" << std::endl << H << "diag_S = " << diag_S << std::endl << "U = " << std::endl << U << "Vt = " << std::endl << Vt << std::endl;
    return U*Vt; // R = U*V^T
}

template <typename T>
inline GMat<T,3,3> GRotation3x3_PolarDecomposition_From_SVD( const GMat<T,3,3> &H, T det_H )
{
    MAL_ASSERT( !IsNaN( H ) );
    //TEMP: Otherwise, eigenvectors are null, for some reason... probably eigenvalues are the same, happens when H==Id
    if( ApproxEq( det_H, T(1) ) ) return GRotation3x3_PolarDecomposition(H,det_H);
    // Compute SVD of H = U * S * V^T
    GMat<T,3,3> U, Vt;
    GVec<T,3> diag_S;
    GSingularValueDecomposition_USVt( H, U, diag_S, Vt ); //AUTOMATICALLY corrects signs
    return U*Vt; // R = U*V^T
}


/*! Polar decomposition to extract rotation R from H = R * S \see FEM.tex
  \note We allow det_H < 0, which will result in an incorrectly oriented R with det_R < 0.

  \todo Should test cPrecisionSq and cMaxIter... it seems that for
  almost-collapsed H it converges quite slowly, but for reasonable
  deformations it converges very fast. The effects of bad convergence
  are skewed and unnormalized R that result in too-large elastic
  forces. To optimize aggresively, pick cPrecisionSq and cMaxIter that
  guarantee good results for the expected deformation space.
  Some tests:
  - cMaxIter = 5 yields low prec unless DTDF is too high
  - cMaxIter = 10 is way better, no asserts with cPrecisionSq = 0.01 and DTDF = 0.06
  \todo NormSqF() MAY be quite expensive, consider doing 2x-3x atomic
  iteration blocks before cheching error instead of doing it
  per-iteration.
*/
template <typename T>
inline GMat<T,3,3> GRotation3x3_PolarDecomposition( const GMat<T,3,3> &H, T det_H )
{
    MAL_ASSERT( Abs(det_H) > Epsilon<T>() ); //should be specific epsilon for actual ops required (sqrt, rcp, etc...)
    GMat<T,3,3> H0(H);
    T error(0);
    unsigned int num_iter(0);
    const unsigned int cMaxIter(10);
    const T cPrecisionSq( 0.01f );
    do
    {
        GMat<T,3,3> H1 = T(0.5) * ( H0 + Inverse(Transposed(H0)) );
        error = (H1-H0).NormSqF();
        H0 = H1;
        num_iter++;
    } while( num_iter < cMaxIter && error > cPrecisionSq );
    /*
    if( error > cPrecisionSq ) std::cout << "PD.error " << error << std::endl;
    else std::cout << "PD.iter " << num_iter << std::endl;
    */
    return H0;
}

/*\todo This may be useful, wrote it for testFE, kept here as reference:
inline Mat2x2 RLerp( const Mat2x2 &R0, const Mat2x2 &R1, Real lambda01 )
{
    Real w1( lambda01 );
    Real w0( Real(1)-w1 );
    Mat2x2 LR( w0*R0 + w1*R1 );
    return GRotation2x2_PolarDecomposition( LR, mal::Det(LR) ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix
}

inline Mat2x2 RCerp( const Mat2x2 &R0, const Mat2x2 &R1, Real lambda01 )
{
    Real w1( mal::Sq(lambda01) * (3-2*lambda01) ); //h01 polynomial from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
    Real w0( Real(1) - w1 );
    Mat2x2 CR( w0*R0 + w1*R1 );
    return GRotation2x2_PolarDecomposition( CR, mal::Det(CR) ); //Use PD to extract "closest proper rotation" to interpolated (non-orthonormal) rotation matrix
}
*/

} // namespace mal

#endif // MAL_GMAT_UTILS_H
