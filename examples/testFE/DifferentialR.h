#ifndef TEST_FE_DIFFERENTIAL_R_H
#define TEST_FE_DIFFERENTIAL_R_H

#include "Config.h"
#include "Element2D.h"
#include "Params.h"

inline Mat2x2r Compute_D( const Vec2r& x0, const Vec2r& x1, const Vec2r& x2 )
{
    return mal::GMat2x2_From_Columns( x1-x0, x2-x0 );
}

/*! Computes dR, from McAdams2011 tech doc
  \note Explicitly assumes S = R^T * F Symmetric
  \note Returns false if dR is 0 (optimization) or infinity (avoid div-by-zero)

  \todo This is dR of the Polar Decomposition, but is it aware of
  potential inversion?? McAdams2011 does NOT mention it... but it DOES
  WORK, tested in testFE
*/
inline bool Compute_dR( const Element2D& element,
                        const Vec2r& dx0, const Vec2r& dx1, const Vec2r& dx2,
                        const Mat2x2r& invDm, const Mat2x2r& F, const Mat2x2r& R, const Mat2x2r& Rt,
                        Mat2x2r& dR )
{
    //\todo In 2D this can be GREATLY OPTIMIZED by computing only the entries of W, S and dR strictly required
    Mat2x2r dDs = Compute_D( dx0, dx1, dx2 ); //\todo Ds can be 0 if 2 or 3 dx are equal
    Mat2x2r dF( dDs * invDm );
    Mat2x2r S( Rt * F );
    Mat2x2r W( Rt * dF );
    Real w01_minus_w10( W(0,1) - W(1,0) );
    Real s00_plus_s11( S(0,0) + S(1,1) );
    //\note Check if ANY factor is 0, s00_plus_s11 to avoid division by 0, w01_minus_w10 to early-out if dR = 0
    //\todo The singularity check may be unnecessary, see 3D version
    if( mal::Abs( w01_minus_w10 * s00_plus_s11 ) > mal::Sq( mal::Epsilon<Real>() ) )
    {
        Real skew_r( w01_minus_w10 / s00_plus_s11 );
        dR = R * Mat2x2r( 0, skew_r, -skew_r, 0 ); //Skew2x2 has just 1 dof, skew_r
        return true;
    }
    else
        return false;
}

/*\todo Compute_dR_PD_U() specific, the ONLY specific we'll implement,
 using dF_bar and dR_bar. Also, see if in COROT/CLRM we need to use
 dF or dF_bar

 \todo We ASSUME noc = 0 and that x1 != x2
*/
inline bool Compute_dR_PD_U( const Element2D& element, const Params& params,
                             const Vec2r& x0, const Vec2r& x1, const Vec2r& x2,
                             const Vec2r& dx0, const Vec2r& dx1, const Vec2r& dx2,
                             const Mat2x2r& invDm, const Mat2x2r& F, const Mat2x2r& R, const Mat2x2r& Rt,
                             Mat2x2r& dR )
{
    if( mal::Det(F) > params.m_DegenerateThresholdDetF ) return Compute_dR( element, dx0, dx1, dx2, invDm, F, R, Rt, dR );
    Vec2r x12( x2 - x1 );
    Real dist12_sq( mal::NormSq(x12) );
    if( dist12_sq > 0.000001f )
    {
        Real dist12( mal::Sqrt(dist12_sq) );

        Vec2r dir_c( mal::PerpendicularCW( x12 ) / dist12 );

        Real h( params.m_DegenerateThresholdDetF * 2 * element.m_Area / dist12 ); //\todo This is delta_h_alpha in other code, and \Delta h in notebooks/LaTeX
        Vec2r x12_mid( 0.5 * (x1 + x2) ); //unnecessary, added for symmetry with dx
        //Vec2r x12_mid( x1 ); //unnecessary, added for symmetry with dx
        Real w( mal::Dot( x0 - x12_mid, dir_c ) );
        Real s( -w / h );
        Real t( Real(1) + s );

        Real t2( t*t );
        Real t3( t*t2 );
        Real h2( h*h );
        Real h3( h*h2 );
        Real w2( w*w );
        Real w3( w*w2 );
        Real s2( s*s );
        Real s3( s*s2 );

        Real lambda = ( w < h )
                      ? ( ( w >= 0 ) ? params.m_Unified_L_Factor*h*( 2*t2 - t3 ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                          : params.m_Unified_L_Factor*h*( t + params.m_Unified_NL_Factor*mal::Pow( s, params.m_Unified_NL_Exponent ) ) ) //lambda_u^L + lambda_u^NL
                      : 0;
        Vec2r x0_bar( x0 + lambda * dir_c );

        // dR(lambda) part
        Vec2r dx12( dx2 - dx1 );
        Vec2r ddir_c( mal::PerpendicularCW( dx12*dist12_sq - x12*mal::Dot(x12,dx12) ) / (dist12*dist12_sq) );
        Vec2r dx12_mid( 0.5 * (dx1 + dx2) ); //midpoint dx, to avoid bias towards either x1 or x2
        Real dw( mal::Dot( dx0 - dx12_mid, dir_c ) + mal::Dot( x0 - x12_mid, ddir_c ) );
        Real Pl_Pw = ( w < h )
                     ? ( ( w >= 0 ) ? -params.m_Unified_L_Factor*( 4*t - 3*t2 ) //dlambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                         : -params.m_Unified_L_Factor*( 1 + params.m_Unified_NL_Factor*params.m_Unified_NL_Exponent*mal::Pow( s, params.m_Unified_NL_Exponent-1 ) ) ) //dlambda_u^L + dlambda_u^NL )
                     : 0;
        Real dlambda_w( Pl_Pw * dw );
        Vec2r dx0_bar( dx0 + dlambda_w * dir_c + lambda * ddir_c );

        Real dh( -params.m_DegenerateThresholdDetF * 2 * element.m_Area * mal::Dot(x12,dx12) / (dist12*dist12_sq) );
        Real Pl_Ph = ( w < h )
                     ? ( ( w >= 0 ) ? params.m_Unified_L_Factor * ( 1 + s2 + 2*s3  )
                         : params.m_Unified_L_Factor*( 1 - params.m_Unified_NL_Factor*(params.m_Unified_NL_Exponent-1)*mal::Pow( s, params.m_Unified_NL_Exponent ) ) )
                     : 0;
        Real dlambda_h( Pl_Ph * dh );
        dx0_bar += dlambda_h * dir_c;

        // Compute corrected \bar Ds, \bar F
        Mat2x2r Ds_bar( Compute_D( x0_bar, x1, x2 ) );
        Mat2x2r F_bar( Ds_bar * invDm );
        Mat2x2r dDs_bar( Compute_D( dx0_bar, dx1, dx2 ) );
        Mat2x2r dF_bar( dDs_bar * invDm );

        //\todo In 2D this can be GREATLY OPTIMIZED by computing only the entries of W, S and dR strictly required
        Mat2x2r S( Rt * F_bar ); //\todo THIS temporary S should be Symmetric, as Rt is computed by PD from F_bar elsewhere
        Mat2x2r W( Rt * dF_bar );
        Real w01_minus_w10( W(0,1) - W(1,0) );
        Real s00_plus_s11( S(0,0) + S(1,1) );
        //\note Check if ANY factor is 0, s00_plus_s11 to avoid division by 0, w01_minus_w10 to early-out if dR = 0
        //\todo The singularity check may be unnecessary, see 3D version
        if( true )//mal::Abs( w01_minus_w10 * s00_plus_s11 ) > mal::Sq( mal::Epsilon<Real>() ) )
        {
            Real skew_r( w01_minus_w10 / s00_plus_s11 );
            dR = R * Mat2x2r( 0, skew_r, -skew_r, 0 ); //Skew2x2 has just 1 dof, skew_r
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

/*! Compute dR by central differences.

  We compute the directional derivative numerically with \hat v
  direction and scale the result with the actual displacement dx

    d R(x,dx) = |dx| * ( R(x + h \hat v) - R(x + h \hat v) ) / 2h

  Where \hat v = Normalized(dx)

  \pre Assumes undegenerate configuration x0,x1,x2
  \todo Toni says it should have precision O(h^2)
*/
inline bool Compute_dR_Numerical( const Vec2r& x0, const Vec2r& x1, const Vec2r& x2,
                                  const Vec2r& dx0, const Vec2r& dx1, const Vec2r& dx2,
                                  const Mat2x2r& invDm, Real degenerate_threshold_det_F,
                                  Mat2x2r& dR )
{
    const Real h( 0.0001f ); //\todo With 10e-3 it ASSERTS during PD due to detF < eps..., as the perturbed left/rightF MAY BE DEGENERATE...
    mal::GVec<Real,6> v;
    mal::GSetRange<0,1>( v, dx0 );
    mal::GSetRange<2,3>( v, dx1 );
    mal::GSetRange<4,5>( v, dx2 );
    Real norm_sq_dx( mal::NormSq( v ) );
    if( norm_sq_dx > mal::Epsilon<Real>() )
    {
        Real norm_dx( mal::Sqrt(norm_sq_dx) );
        v = v / norm_dx;
        Vec2r v0 = mal::GRange<0,1>( v );
        Vec2r v1 = mal::GRange<2,3>( v );
        Vec2r v2 = mal::GRange<4,5>( v );
        Mat2x2r leftDs, leftF, leftR;
        Mat2x2r rightDs, rightF, rightR;
        leftDs = Compute_D( x0 - h*v0, x1 - h*v1, x2 - h*v2 );
        rightDs = Compute_D( x0 + h*v0, x1 + h*v1, x2 + h*v2 );
        leftF = leftDs * invDm;
        rightF = rightDs * invDm;
        Real det_LF( mal::Det(leftF) );
        Real det_RF( mal::Det(rightF) );

        /*\todo We should consider possible PD_P/R/U or SVD1 rotation
          corrections here ?!?!
          By now, we just choose bidirectional or unidirectional
          derivation according to det_LF/RF, and assume that current F
          is undegenerate.

          \todo We tried to use From_SVD version to avoid problems in
          degenerate configurations, but all methods using it started
          to diverge wildly... there must be some caveat
        */
        if( det_LF > degenerate_threshold_det_F && det_RF > degenerate_threshold_det_F )
        {
            leftR = mal::GRotation2x2_PolarDecomposition( leftF, det_LF );
            rightR = mal::GRotation2x2_PolarDecomposition( rightF, det_RF );
            dR = norm_dx * mal::Rcp(2*h) * (rightR - leftR);
            return true;
        }
        else if( det_LF > degenerate_threshold_det_F )
        {
            leftR = mal::GRotation2x2_PolarDecomposition( leftF, det_LF );
            rightDs = Compute_D( x0, x1, x2 );
            rightF = rightDs * invDm;
            rightR = mal::GRotation2x2_PolarDecomposition( rightF, mal::Det(rightF) );
            dR = norm_dx * mal::Rcp(h) * (rightR - leftR);
            return true;
        }
        else if( det_RF > degenerate_threshold_det_F )
        {
            leftDs = Compute_D( x0, x1, x2 );
            leftF = leftDs * invDm;
            leftR = mal::GRotation2x2_PolarDecomposition( leftF, mal::Det(leftF) );
            rightR = mal::GRotation2x2_PolarDecomposition( rightF, det_RF );
            dR = norm_dx * mal::Rcp(h) * (rightR - leftR);
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

/*! Compute dR by central differences.

  We compute the directional derivative numerically with \hat v
  direction and scale the result with the actual displacement dx

    d R(x,dx) = |dx| * ( R(x + h \hat v) - R(x + h \hat v) ) / 2h

  Where \hat v = Normalized(dx)

  \pre Assumes undegenerate configuration x0,x1,x2
  \todo Toni says it should have precision O(h^2)

  \todo Error is low in general position, but quite significant near X
  axis, near Y axis (collapse) and near the
  SVD-critical-point. Increasing h decreases the first two errors, but
  increases the last. Step h = 10e-2 seems the best compromise
  neighbourhood, as it keeps error reasonably low in all problematic
  areas.
*/
inline bool Compute_dR_Numerical( const Element2D& element, Params::EMethod method,
                                  const Vec2r& x0, const Vec2r& x1, const Vec2r& x2,
                                  const Vec2r& dx0, const Vec2r& dx1, const Vec2r& dx2,
                                  const Mat2x2r& invDm, Real degenerate_threshold_det_F,
                                  Real h,
                                  Mat2x2r& dR )
{
    //const Real h( 0.01f );
    mal::GVec<Real,6> v;
    mal::GSetRange<0,1>( v, dx0 );
    mal::GSetRange<2,3>( v, dx1 );
    mal::GSetRange<4,5>( v, dx2 );
    Real norm_sq_dx( mal::NormSq( v ) );
    if( norm_sq_dx > mal::Epsilon<Real>() )
    {
        Real norm_dx( mal::Sqrt(norm_sq_dx) );
        v = v / norm_dx;
        Vec2r v0 = mal::GRange<0,1>( v );
        Vec2r v1 = mal::GRange<2,3>( v );
        Vec2r v2 = mal::GRange<4,5>( v );
        Vec2r leftX[] = { x0 - h*v0, x1 - h*v1, x2 - h*v2 };
        Vec2r rightX[] = { x0 + h*v0, x1 + h*v1, x2 + h*v2 };
        Mat2x2r leftR, rightR;
        element.ComputeElasticForce( leftX, method, &leftR );
        element.ComputeElasticForce( rightX, method, &rightR );
        dR = norm_dx * mal::Rcp(2*h) * (rightR - leftR);
        return true;
    }
    else
        return false;
}

inline bool Compute_dR_QR_XY( const Vec2r& dx0, const Vec2r& dx1, const Vec2r& dx2,
                              const Mat2x2r& invDm, const Mat2x2r& F,
                              Mat2x2r& dR )
{
    Mat2x2r dDs = Compute_D( dx0, dx1, dx2 );
    Mat2x2r dF( dDs * invDm );
    Vec2r f0 = mal::GColumn<0>(F);
    Vec2r df0 = mal::GColumn<0>(dF);
    Real norm_f0 = mal::Norm( f0 );
    if( norm_f0 > 10e-3f )
    {
        Vec2r de0 = ( mal::Sq(norm_f0)*df0 - mal::Dot(f0,df0)*f0 ) / (norm_f0 * norm_f0 * norm_f0);
        dR = mal::GMat2x2_From_Columns( de0, mal::PerpendicularCW(de0) );
        return true;
    }
    else
        return false;
}

#endif //TEST_FE_DIFFERENTIAL_R_H
