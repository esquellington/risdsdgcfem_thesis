#ifndef MAL_GMAT_EIGEN_H
#define MAL_GMAT_EIGEN_H

#include <Mal/Config.h>
#include <Mal/GMat.h>
#include <Mal/GMatUtils.h>
#include <Mal/GVec.h>
#include <Mal/GRandom.h>

namespace mal
{

/*! Compute largest (eigen_value,eigen_vector) pair of a real
    symmetric matrix.

    Uses the Power Iteration method to find the largest (absolute)
    eigenvalue, which is guaranteed to be Real for a Symmetric real
    matrix.

    Returns absolute error of the last two eigen_value approximations
    computed, which will be ideally < epsilon unless max_iter is hit
    or the algorithm fails to converge.

    Returns < 0 if no significant result can be found (may be due to
    random seed vector)
*/
template <typename T, int N>
inline T GComputeLargestEigenvalue_Symmetric( const GMat<T,N,N> &m,
                                              T &eigen_value, GVec<T,N> &eigen_vector,
                                              unsigned max_iter, T epsilon )
{
    GVec<T,N> v0;
    T e0(0);
    GVec<T,N> v1( RandomV( GVec<T,N>(-1), GVec<T,N>(1) ) );
    T e1( mal::Norm(v1) );
    MAL_ASSERT( e1 > epsilon );
    T error( e1-e0 );
    unsigned int iter(0);
    do
    {
        e0 = e1;
        v0 = v1 / e1;
        v1 = m * v0;
        e1 = mal::Norm(v1);
        error = mal::Abs(e1-e0);
        iter++;
    }
    while( iter < max_iter //Limit cost
           && e1 > epsilon //Ensure non-collapsing eigenvector
           && error > epsilon ); //Ensure significant convergence
    // Compute final eigen_vector and its sign
    if( e1 > epsilon )
    {
        MAL_ASSERT( e1 > epsilon );
        eigen_vector = v1 / e1;
        // Compute eigen_value sign by checking if v0*v1 < 0 (which means it has inverted when multiplied by M in v1 = M * v0)
        eigen_value = ( Dot(v0,v1) >= 0 ) ? e1 : -e1;
        return error;
    }
    // No non-zero component found in iterated eigen_vector => no valid eigen_value
    eigen_value = 0;
    eigen_vector = GVec<T,N>::Zero();
    return T(-1);
}

/*! Same as GComputeLargestEigenvalue_Symmetric, but with eigenvalue shifting */
template <typename T, int N>
inline T GComputeLargestEigenvalue_Symmetric_Shifted( const GMat<T,N,N> &m,
                                                      T shift_eigen_value,
                                                      unsigned max_iter, T epsilon,
                                                      T &eigen_value, GVec<T,N> &eigen_vector )
{
    GVec<T,N> v0;
    T e0(0);
    GVec<T,N> v1( RandomV( GVec<T,N>(-1), GVec<T,N>(1) ) );
    T e1( mal::Norm(v1) );
    MAL_ASSERT( e1 > epsilon );
    T error( e1-e0 );
    unsigned int iter(0);
    do
    {
        e0 = e1;
        v0 = v1 / e1;
        v1 = m * v0;
        v1 -= shift_eigen_value * v0; //v1' = v1 - shift*Id*v0
        e1 = mal::Norm(v1);
        error = mal::Abs(e1-e0);
        iter++;
    }
    while( iter < max_iter //Limit cost
           && e1 > epsilon //Ensure non-collapsing eigenvector
           && error > epsilon ); //Ensure significant convergence
    // Compute final eigen_vector and its sign
    if( e1 > epsilon )
    {
        MAL_ASSERT( e1 > epsilon );
        eigen_vector = v1 / e1;
        // Compute eigen_value sign by checking if v0*v1 < 0 (which means it has inverted when multiplied by M in v1 = M * v0)
        eigen_value = ( Dot(v0,v1) >= 0 ) ? e1 : -e1;
        return error;
    }
    // No non-zero component found in iterated eigen_vector => no valid eigen_value
    eigen_value = 0;
    eigen_vector = GVec<T,N>::Zero();
    return T(-1);
}

template <typename T, int N>
inline T GComputeMinMaxEigenvalueRatio_Symmetric( const GMat<T,N,N> &m,
                                                  unsigned max_iter, T epsilon )
{
    // Find largest EV of A using power method
    GVec<T,N> eigen_vector;
    T lambda_max(0);
    GComputeLargestEigenvalue_Symmetric( m,
                                         lambda_max, eigen_vector,
                                         max_iter, epsilon );
    // Find smallest EV of A using shifted power method (http://math.stackexchange.com/questions/271864/power-iteration-smallest-eigenvalue)
    T lambda_min(0);
    GComputeLargestEigenvalue_Symmetric_Shifted( m,
                                                 lambda_max,
                                                 max_iter, epsilon,
                                                 lambda_min, eigen_vector );
    lambda_min += lambda_max;
    return lambda_min/lambda_max;
}

/*\todo Compute all (or just NV largest!??!?!?) (eigen_value,eigen_vector) pairs
template <typename T, int N>
inline T ComputeEigenvalues_Symmetric( const GMat<T,N,N> m,
                                       GVec<T,N> &vec_eigen_values, GMat<T,N> &mat_eigen_vectors,
                                       unsigned num_values, unsigned max_iter, T epsilon )
{
    return T(-1);
}
*/

////////////////////////////////////////////////////////////////
// SVD-related stuff
////////////////////////////////////////////////////////////////

/*! Computes real eigenvalues [lambda1 >= lambda2] of a Symmetric 2x2
    matrix as the solutions of the quadratic equation:
       det( M - lambda * Id ) = 0
       =>
       lambda^2 - lambda * Trace(m) + Det(m) = 0
    \note if not Symmetric, lambdas could be complex?...? or something
    \note For SVD, l1,l2 should be >=0, but l = -0.000000 has been observed, consider
          checking here insead of later in GSingularValueDecomposition_USVt
*/
template <typename T>
inline GVec<T,2> ComputeEigenvalues_Symmetric( const GMat<T,2,2> &m )
{
    //Quadratic equation: ax² + bx + c = 0 with a = 1, b = -Trace(m), c = Det(m)
    T b( -Trace(m) );
    T c( Det(m) );
    T pm( Sqrt(b*b - T(4)*c) );
    T l1( T(0.5) * (-b - pm) );
    T l2( T(0.5) * (-b + pm) );
    /* \note NO need to sort, as l2 >= l1 because pm >= 0
    if( l1 >= l2 ) return GVec<T,2>( l1, l2 );
    else return GVec<T,2>( l2, l1 );
    */
    return GVec<T,2>( l2, l1 );
}

/*! Computes the eigenvector of m associated to the given eigenvalue lambda
  M*x = lambda*x
  =>
  (M - lambda*Id)*x = 0
  =>
  E*x = 0
*/
template <typename T>
inline GVec<T,2> ComputeEigenvectorFromEigenvalue( const GMat<T,2,2> &m, T lambda )
{
    GMat<T,2,2> E( m(0,0) - lambda, m(0,1),
                   m(1,0), m(1,1) - lambda );
    if( !IsZero( GRow<0>(E) ) )
        return GVec<T,2>( E(0,1), -E(0,0) ).Normalized();
    else if( !IsZero( GRow<1>(E) ) )
        return GVec<T,2>( -E(1,1), E(1,0) ).Normalized();
    else //everything is zero
        return GVec<T,2>::Zero();
}

template <typename T>
GMat<T,2,2> ComputeEigenbasisFromEigenvalues( const GMat<T,2,2> &m, const GVec<T,2> &eigen_values )
{
    return GMat2x2_From_Columns( ComputeEigenvectorFromEigenvalue(m,eigen_values[0]),
                                 ComputeEigenvectorFromEigenvalue(m,eigen_values[1]) );
}

} // namespace mal

#endif // MAL_GMAT_EIGEN_H
