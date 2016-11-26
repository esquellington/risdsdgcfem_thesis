#include "Rotation.h"
#include "TriangleElement2.h"
#include "TetrahedronElement3.h"
#include <Mal/GMatDecomposition.h>

#define __ENABLE_REGULARIZE_DAPD
#ifdef __ENABLE_REGULARIZE_DAPD
#  define REGULARIZE_DAPD_THRESHOLD 0.001f
#endif

namespace S2 {
namespace ms {
namespace fem {

////////////////////////////////////////////////////////////////
// Rotation
////////////////////////////////////////////////////////////////

Mat2x2 Compute_R_SVD( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F )
{
    //\note Undegenerate includes det(F) == 1, which FAILS for GSingularValueDecomposition_USVt() due to repeated eigenvalues, resulting in NaN
    if( detF > degenerate_threshold_det_F ) //Undegenerate
    {
        return mal::GRotation2x2_PolarDecomposition( F, detF );
    }
    else //Degenerate or Inverted
    {
        // Compute SVD
        Mat2x2 U, Vt;
        Vec2 diag_S;
        mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt );
        // Force pure rotation U,Vt by fixing potential inversion/reflection
        if( mal::Det(U) < 0 ) mal::GSetColumn<1>( U, -mal::GColumn<1>(U) );
        if( mal::Det(Vt) < 0 ) mal::GSetRow<1>( Vt, -mal::GRow<1>(Vt) );
        return U*Vt;
    }
}

Mat3x3 Compute_R_SVD( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F )
{
    /*TEMP: USING PURE SVD from GSL
    //\note Undegenerate includes det(F) == 1, which FAILS for GSingularValueDecomposition_USVt() due to repeated eigenvalues, resulting in NaN
    if( detF > m_Params.m_DegenerateThresholdDetF ) //Undegenerate
    {
    R = mal::GRotation3x3_PolarDecomposition( F, detF );
    MS_ASSERT( !mal::IsNaN( R ) );
    }
    else //Degenerate or Inverted
    */
    {
        // Compute SVD
        Mat3x3 U, Vt;
        Vec3 diag_S;
        mal::GSingularValueDecomposition_USVt( F, U, diag_S, Vt ); //Automagically Forces pure rotation U,Vt by fixing potential inversion/reflection
        return U*Vt;
    }
}

Mat2x2 Compute_R_PDP_Degenerate( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F,
                                 int32 noc,
                                 const Vec2& p0, const Vec2& p1, const Vec2& p2,
                                 Real area, const Mat2x2& inv_Dm )
{
    MS_ASSERT( noc >= 0 && noc <= 2 );
    // use noc and get x0, x1, x2 accordingly
    Vec2 localX[3] = { p0, p1, p2 };
    Vec2 &x0( localX[ noc ] ); //Inverted vtx
    const Vec2 &x1( localX[ (noc + 1) % 3 ] );
    const Vec2 &x2( localX[ (noc + 2) % 3 ] );
    Vec2 x12( x2 - x1 );
    Real dist12_sq( mal::NormSq(x12) );
    if( dist12_sq > 0.000001f )
    {
        Vec2 axis12( mal::PerpendicularCW( x12 ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
        // \note See paper DCLFEM (Detail) for delta_h algorithm that avoids normalization, using unnormalized axis12 and dist12_sq.
        Real delta_div_dist12( degenerate_threshold_det_F * 2*area / dist12_sq ); //== delta / dist12
        x0 -= ( mal::Dot(x0-x1,axis12) / dist12_sq ) * axis12; //Project onto collapse line, gathering 1/dist12 factors into a single 1/dist12_sq one
        x0 += delta_div_dist12*axis12; //Move past collapse
        // Compute Proj(F)
        Mat2x2 projDs( mal::GMat2x2_From_Columns( localX[1] - localX[0], localX[2] - localX[0] ) );
        Mat2x2 projF( projDs * inv_Dm );
        Real det_projF( mal::Det(projF) );
        MS_ASSERT( det_projF >= 0 );
        if( det_projF < 0.99f*degenerate_threshold_det_F )
            MS_LOG_WARNING( "Element det(proj(F)) %f != %f should not happen if delta is correct", det_projF, degenerate_threshold_det_F );
        Mat2x2 R( mal::GRotation2x2_PolarDecomposition( projF, det_projF ) );
        MS_ASSERT( !mal::IsNaN( R ) );
        return R;
    }
    else
    {
        //\todo Consider using last R!!
        MS_LOG_ERROR( "Degenerate element with NoC %d has dist12_sq = %f too small, collapsed? Using R=Id ", noc, dist12_sq );
        return Mat2x2::Identity();
    }
}

Mat3x3 Compute_R_PDP_Degenerate( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F,
                                 TetrahedronElement3::DoC doc,
                                 const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                                 Real volume, const Mat3x3& inv_Dm )
{
    // Get local nodes
    Vec3 localX[4] = { p0, p1, p2, p3 };
    if( doc.IsVF() )
    {
        int noc( doc.m_VIT0 );
        Vec3 &x0( localX[ noc ] ); //Inverted vtx
        const Vec3 &x1( localX[ (noc + 1) % 4 ] );
        const Vec3 &x2( localX[ (noc + 2) % 4 ] );
        const Vec3 &x3( localX[ (noc + 3) % 4 ] );
        Vec3 axis123( mal::Cross(x2-x1, x3-x1) );
        Real length123_sq( mal::NormSq(axis123) ); //\note This is the term B_{234}^2 in DCLFEM paper
        if( length123_sq > 0.000001f )
        {
            // \note See paper DCLFEM (Detail) for delta_h algorithm that avoids normalization, using unnormalized axis123 and length123_sq.
            Real delta_div_dist123( degenerate_threshold_det_F * 6 * volume / length123_sq ); //== delta / length123_sq
            //\todo ENFORCE axis123 orientation... it MAY be possible to guarantee it selecting proper x1,x2,x3 for all possible x0
            Real dot_x0_minux_x1_axis123( mal::Dot(x0-x1,axis123) );
            /*\todo As axis123 may NOT be consistently
              oriented, we flip it ONLY if dot and detF do not
              have the expected signs, that is, either x0 is
              BEHIND (x1,x2,x3), with detF < 0, and dot < 0,
              or IN FRONT, with 0 <= detF < DTDF and 0 <= dot.
            */
            if( dot_x0_minux_x1_axis123 * detF < 0 )
            {
                dot_x0_minux_x1_axis123 = -dot_x0_minux_x1_axis123;
                axis123 = -axis123;
            }
            x0 -= ( dot_x0_minux_x1_axis123 / length123_sq ) * axis123; //Project onto collapse line, gathering 1/length123_sq factors into a single 1/length123_sq one
            x0 += delta_div_dist123*axis123; //Move past collapse
            // Compute Proj(F)
            Mat3x3 projDs( mal::GMat3x3_From_Columns( localX[1] - localX[0], localX[2] - localX[0], localX[3] - localX[0] ) );
            Mat3x3 projF( projDs * inv_Dm );
            Real det_projF( mal::Det(projF) );
            MS_ASSERT( det_projF >= 0 );
            if( det_projF < 0.99f*degenerate_threshold_det_F )
                MS_LOG_WARNING( "det(proj(F)) %f != %f should not happen if delta is correct", det_projF, degenerate_threshold_det_F );
            return mal::GRotation3x3_PolarDecomposition( projF, det_projF );
        }
        else
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element with V-F DoC %d has length123_sq = %f too small, collapsed? Using R=Id ", doc.m_VIT0, length123_sq );
            return Mat3x3::Identity();
        }
    }
    else //doc.IsEE()
    {
        int vit0( doc.m_VIT0 );
        MS_ASSERT( 0 == vit0 );
        int vit1( doc.m_VIT1 );
        //\todo FIND A BETTER WAY to compute vid2,vid3...
        int vit2, vit3;
        switch( vit1 )
        {
        case 1: vit2 = 2; vit3 = 3; break;
        case 2: vit2 = 1; vit3 = 3; break;
        case 3: vit2 = 1; vit3 = 2; break;
        default: vit2 = vit3 = 1000; break;
        }
        MS_ASSERT( vit0 + vit1 + vit2 + vit3 == 0+1+2+3 ); //TEMP
        Vec3 &x0( localX[ vit0 ] );
        Vec3 &x1( localX[ vit1 ] );
        const Vec3 &x2( localX[ vit2 ] );
        const Vec3 &x3( localX[ vit3 ] );
        Vec3 axis123( mal::Cross(x1-x0, x3-x2) ); //e_ij x e_kl
        Real length123_sq( mal::NormSq(axis123) );
        if( length123_sq > 0.000001f )
        {
            // \note See paper DCLFEM (Detail) for delta_h algorithm that avoids normalization, using unnormalized axis123 and length123_sq.
            Real delta_div_dist123( degenerate_threshold_det_F * 6 * volume / length123_sq ); //== delta / length123_sq
            //\todo ENFORCE axis123 orientation... it MAY be possible to guarantee it selecting proper x2,x3 for all possible x0,x1
            Real dot_x0_minux_x2_axis123( mal::Dot(x0-x2,axis123) );
            Real dot_x1_minux_x2_axis123( mal::Dot(x1-x2,axis123) );
            MS_ASSERT( mal::ApproxEq( dot_x0_minux_x2_axis123, dot_x1_minux_x2_axis123 ) );
            /*\todo As axis123 may NOT be consistently
              oriented, we flip it ONLY if dot and detF do not
              have the expected signs, that is, either x0 is
              BEHIND (x1,x2,x3), with detF < 0, and dot < 0,
              or IN FRONT, with 0 <= detF < DTDF and 0 <= dot.
            */
            if( dot_x0_minux_x2_axis123 * detF < 0 )
            {
                dot_x0_minux_x2_axis123 = -dot_x0_minux_x2_axis123;
                dot_x1_minux_x2_axis123 = -dot_x1_minux_x2_axis123;
                axis123 = -axis123;
            }
            x0 -= ( dot_x0_minux_x2_axis123 / length123_sq ) * axis123; //Project onto collapse line, gathering 1/length123_sq factors into a single 1/length123_sq one
            x0 += delta_div_dist123*axis123; //Move past collapse
            x1 -= ( dot_x1_minux_x2_axis123 / length123_sq ) * axis123; //Project onto collapse line, gathering 1/length123_sq factors into a single 1/length123_sq one
            x1 += delta_div_dist123*axis123; //Move past collapse
            // Compute Proj(F)
            Mat3x3 projDs( mal::GMat3x3_From_Columns( localX[1] - localX[0], localX[2] - localX[0], localX[3] - localX[0] ) );
            Mat3x3 projF( projDs * inv_Dm );
            Real det_projF( mal::Det(projF) );
            MS_ASSERT( det_projF >= 0 );
            if( det_projF < 0.99f*degenerate_threshold_det_F )
                MS_LOG_WARNING( "det(proj(F)) %f != %f should not happen if delta is correct", det_projF, degenerate_threshold_det_F );
            return mal::GRotation3x3_PolarDecomposition( projF, det_projF );
        }
        else
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element with E-E DoC (%d,%d) has length123_sq = %f too small, collapsed? Using R=Id ", doc.m_VIT0, doc.m_VIT1, length123_sq );
            return Mat3x3::Identity();
        }
    }
}

Mat2x2 Compute_R_PDR_Degenerate( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F,
                                 int32 noc,
                                 const Vec2& p0, const Vec2& p1, const Vec2& p2 )
{
    // use noc and get x0, x1, x2 accordingly
    Vec2 localX[3] = { p0, p1, p2 };
    Vec2 &x0( localX[ noc ] ); //Inverted vtx
    const Vec2 &x1( localX[ (noc + 1) % 3 ] );
    const Vec2 &x2( localX[ (noc + 2) % 3 ] );
    Vec2 x12( x2 - x1 );
    Real dist12_sq( mal::NormSq(x12) );
    if( dist12_sq > 0.000001f )
    {
        Vec2 axis12( mal::Normalized( mal::PerpendicularCW( x12 ) ) ); //\note axis12 points towards uninverted x0 halfspace, as we KNOW which axis inverted a priori
        if( detF <= -degenerate_threshold_det_F ) //Inverted & !Collapsed
            return mal::GRotation2x2_PolarDecomposition( mal::GReflection2x2_From_Axis( axis12 ) * F, -detF ); //Safe uncollapsed simple PD
        else if( detF < 0 ) //-m_DegenerateThresholdDetF < det(F) < 0 => Inverted & Collapsed
            return mal::GRotation2x2_PolarDecomposition_From_SVD( mal::GReflection2x2_From_Axis( axis12 ) * F, -detF ); //\todo SHOULD Handle det(F) = 0 properly
        else // 0 <= det(F) < m_DegenerateThresholdDetF => !Inverted && Collapsed
            return mal::GRotation2x2_PolarDecomposition_From_SVD( F, detF ); //\todo SHOULD Handle det(F) = 0 properly, but yields NaN sometimes
    }
    else
    {
        //\todo Consider using last R!!
        MS_LOG_ERROR( "Degenerate element with NoC %d has dist12_sq = %f too small, collapsed? Using R=Id ", noc, dist12_sq );
        return Mat2x2::Identity();
    }
}

Mat3x3 Compute_R_PDR_Degenerate( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F,
                                 TetrahedronElement3::DoC doc,
                                 const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3 )
{
    // Get local nodes
    Vec3 localX[4] = { p0, p1, p2, p3 };
    Vec3 axis123(0,0,0);
    Real length123_sq(0);
    if( doc.IsVF() )
    {
        int noc( doc.m_VIT0 );
        //const Vec3 &x0( localX[ noc ] ); //Inverted vtx
        const Vec3 &x1( localX[ (noc + 1) % 4 ] );
        const Vec3 &x2( localX[ (noc + 2) % 4 ] );
        const Vec3 &x3( localX[ (noc + 3) % 4 ] );
        axis123 = mal::Cross(x2-x1, x3-x1);
        length123_sq = mal::NormSq(axis123);
    }
    else //doc.IsEE()
    {
        int vit0( doc.m_VIT0 );
        MS_ASSERT( 0 == vit0 );
        int vit1( doc.m_VIT1 );
        //\todo FIND A BETTER WAY to compute vid2,vid3...
        int vit2, vit3;
        switch( vit1 )
        {
        case 1: vit2 = 2; vit3 = 3; break;
        case 2: vit2 = 1; vit3 = 3; break;
        case 3: vit2 = 1; vit3 = 2; break;
        default: vit2 = vit3 = 1000; break;
        }
        MS_ASSERT( vit0 + vit1 + vit2 + vit3 == 0+1+2+3 ); //TEMP
        const Vec3 &x0( localX[ vit0 ] );
        const Vec3 &x1( localX[ vit1 ] );
        const Vec3 &x2( localX[ vit2 ] );
        const Vec3 &x3( localX[ vit3 ] );
        axis123 = mal::Cross(x1-x0, x3-x2); //e_ij x e_kl
        length123_sq = mal::NormSq(axis123);
    }
    // Use axis to compute reflection
    if( length123_sq > 0.000001f ) //\todo set proper epsilon
    {
        axis123 /= mal::Sqrt(length123_sq); //\note The orientation of axis123 does NOT matter here, as reflection along v and -v is the same.
        if( detF <= -degenerate_threshold_det_F ) //Inverted & !Collapsed
            return mal::GRotation3x3_PolarDecomposition( mal::GReflection3x3_From_Axis( axis123 ) * F, -detF ); //Safe uncollapsed simple PD
        else if( detF < 0 ) //-m_DegenerateThresholdDetF < det(F) < 0 => Inverted & Collapsed
            return mal::GRotation3x3_PolarDecomposition_From_SVD( mal::GReflection3x3_From_Axis( axis123 ) * F, -detF ); //\todo SHOULD Handle det(F) = 0 properly
        else // 0 <= det(F) < m_DegenerateThresholdDetF => !Inverted && Collapsed
            return mal::GRotation3x3_PolarDecomposition_From_SVD( F, detF ); //\todo SHOULD Handle det(F) = 0 properly, but yields NaN sometimes
    }
    else
    {
        //\todo Consider using last R!!
        MS_LOG_ERROR( "Degenerate element with DoC (%d,%d) has length123_sq = %f too small, collapsed? Using R=Id ", doc.m_VIT0, doc.m_VIT1, length123_sq );
        return Mat3x3::Identity();
    }
}

Mat2x2 Compute_R_DAPD_Degenerate( const Mat2x2& F, Real detF, Real degenerate_threshold_det_F,
                                  int32 noc,
                                  const Vec2& p0, const Vec2& p1, const Vec2& p2,
                                  Real area, const Mat2x2& inv_Dm,
                                  Real factor_L, Real factor_NL, Real exponent_NL )
{
    // use noc and get x0, x1, x2 accordingly
    Vec2 localX[3] = { p0, p1, p2 };
    Vec2 &x0( localX[ noc ] ); //Inverted vtx
    const Vec2 &x1( localX[ (noc + 1) % 3 ] );
    const Vec2 &x2( localX[ (noc + 2) % 3 ] );
    Vec2 x12( x2 - x1 );
    Real dist12_sq( mal::NormSq(x12) );
#ifdef __ENABLE_REGULARIZE_DAPD
    //dist12_sq = mal::Max( dist12_sq, mal::Sq(REGULARIZE_DAPD_THRESHOLD) );
    dist12_sq += mal::Sq(REGULARIZE_DAPD_THRESHOLD);
    if( true )
#else
    if( dist12_sq > 0.000001f )
#endif
    {
        /*\todo OPTIMIZE:
          1) Try to avoid sqrt as in Project method (w seems to need unitary axis12, but try it...)
          2) Hardcoded integer exponent 2 or 3
        */
        Real dist12( mal::Sqrt(dist12_sq) );
        Vec2 axis12( mal::PerpendicularCW( x12 ) / dist12 );
        Real delta_h_alpha( degenerate_threshold_det_F * 2 * area / dist12 );
        Vec2 x12_mid( 0.5 * (x1 + x2) ); //unnecessary, added for symmetry with dx
        Real w( mal::Dot(x0-x12_mid,axis12) ); //\note REQUIRES normalized axis
        Real s( w / delta_h_alpha );
        Real t( Real(1) - s );
        // Cubic spline from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
        Real lambda_u = ( w >= 0 ) //We know that w <= \Delta h as det_F <= degenerate_threshold_det_F
                        ? factor_L*delta_h_alpha*( (2-t)*t*t ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                        : factor_L*delta_h_alpha*( t + factor_NL*mal::Pow( -s, exponent_NL ) ); //lambda_u^L + lambda_u^NL

        x0 += lambda_u * axis12; //Move past collapse area, to ensure numerically safe
        // Compute Proj(F)
        Mat2x2 projDs( mal::GMat2x2_From_Columns( localX[1] - localX[0], localX[2] - localX[0] ) );
        Mat2x2 projF( projDs * inv_Dm );
        Real det_projF( mal::Det(projF) );
        //MS_ASSERT( det_projF >= 0 ); TEMPORALLY DISABLED
        /*\todo Actually, this CAN HAPPEN because PD_U lambda does not ensure delta_h but DOES ensure > 0
        //TEMP: Removed due to verbosity...
        if( det_projF < 0.99f*degenerate_threshold_det_F )
        MS_LOG_WARNING( "det(proj(F)) %f != %f should not happen if delta is correct (det_F was %f)", det_projF, degenerate_threshold_det_F, detF );
        */
        if( det_projF > 0.001f ) //\todo This may happen because PD_U does not guarantee minimum det(\bar F)
            return mal::GRotation2x2_PolarDecomposition( projF, det_projF );
        else
        {
            MS_LOG_WARNING( "det(proj(F)) %f <= epsilon %f", det_projF, 0.001f );
            return mal::GRotation2x2_PolarDecomposition_From_SVD( projF, det_projF );
        }
    }
    else
    {
        //\todo Consider using last R!!
        MS_LOG_ERROR( "Degenerate element with NoC %d has dist12_sq = %f too small, collapsed? Using R=Id ", noc, dist12_sq );
        return Mat2x2::Identity();
    }
}

Mat3x3 Compute_R_DAPD_Degenerate( const Mat3x3& F, Real detF, Real degenerate_threshold_det_F,
                                  TetrahedronElement3::DoC doc,
                                  const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                                  Real volume, const Mat3x3& inv_Dm,
                                  Real factor_L, Real factor_NL, Real exponent_NL )
{
    MS_ASSERT( doc.IsValid() );
    // Get local nodes
    Vec3 localX[4] = { p0, p1, p2, p3 };
    if( doc.IsVF() )
    {
        int noc( doc.m_VIT0 );
        Vec3 &x0( localX[ noc ] ); //Inverted vtx
        const Vec3 &x1( localX[ (noc + 1) % 4 ] );
        const Vec3 &x2( localX[ (noc + 2) % 4 ] );
        const Vec3 &x3( localX[ (noc + 3) % 4 ] );
        Vec3 axis123( mal::Cross(x2-x1, x3-x1) );
        Real length123_sq( mal::NormSq(axis123) ); //\note This is the term B_{234}^2 in DCLFEM paper
        if( length123_sq > 0.000001f )
        {
            /*\todo OPTIMIZE:
              1) Try to avoid sqrt as in Project method (w seems to need unitary axis12, but try it...)
              2) Hardcoded integer exponent 1, 2 or 3
            */
            Real length123( mal::Sqrt(length123_sq) );
            axis123 = axis123 / length123;
            //\todo ENFORCE axis123 orientation... it MAY be possible to guarantee it selecting proper x2,x3 for all possible x0,x1
            Real delta_h_alpha( degenerate_threshold_det_F * 6 * volume / length123 );
            Real w( mal::Dot(x0-x1,axis123) ); //\note REQUIRES normalized axis
            /*\todo As axis123 may NOT be consistently
              oriented, we flip it ONLY if w and detF do not
              have the expected signs, that is, either x0 is
              BEHIND (x1,x2,x3), with detF < 0, and dot < 0,
              or IN FRONT, with 0 <= detF < DTDF and 0 <= dot.
            */
            if( w * detF < 0 )
            {
                w = -w;
                axis123 = -axis123;
            }
            Real s( w / delta_h_alpha );
            Real t( Real(1) - s );
            // Cubic spline from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
            Real lambda_u = ( w > 0 ) //We know that w <= \Delta h as det_F <= m_rParams.m_DegenerateThresholdDetF
                            ? factor_L*delta_h_alpha*( (2 - t)*t*t ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                            : factor_L*delta_h_alpha*( t + factor_NL*mal::Pow( -s, exponent_NL ) ); //lambda_u^L + lambda_u^NL
            x0 += lambda_u * axis123; //Move past collapse area, to ensure numerically safe
            // Compute Proj(F)
            Mat3x3 projDs( mal::GMat3x3_From_Columns( localX[1] - localX[0], localX[2] - localX[0], localX[3] - localX[0] ) );
            Mat3x3 projF( projDs * inv_Dm );
            Real det_projF( mal::Det(projF) );
            MS_ASSERT( det_projF >= 0 );
            if( det_projF < 0.99f*degenerate_threshold_det_F )
                MS_LOG_WARNING( "det(proj(F)) %f != %f should not happen if delta is correct", det_projF, degenerate_threshold_det_F );
            return mal::GRotation3x3_PolarDecomposition( projF, det_projF );
        }
        else
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element with V-F DoC %d has length123_sq = %f too small, collapsed? Using R=Id ", doc.m_VIT0, length123_sq );
            return Mat3x3::Identity();
        }
    }
    else //doc.IsEE()
    {
        int vit0( doc.m_VIT0 );
        MS_ASSERT( 0 == vit0 );
        int vit1( doc.m_VIT1 );
        //\todo FIND A BETTER WAY to compute vid2,vid3...
        int vit2, vit3;
        switch( vit1 )
        {
        case 1: vit2 = 2; vit3 = 3; break;
        case 2: vit2 = 1; vit3 = 3; break;
        case 3: vit2 = 1; vit3 = 2; break;
        default: vit2 = vit3 = 1000; break;
        }
        MS_ASSERT( vit0 + vit1 + vit2 + vit3 == 0+1+2+3 ); //TEMP
        Vec3 &x0( localX[ vit0 ] );
        Vec3 &x1( localX[ vit1 ] );
        const Vec3 &x2( localX[ vit2 ] );
        const Vec3 &x3( localX[ vit3 ] );
        Vec3 axis123( mal::Cross(x1-x0, x3-x2) ); //e_ij x e_kl
        Real length123_sq( mal::NormSq(axis123) );
        if( length123_sq > 0.000001f )
        {
            /*\todo OPTIMIZE:
              1) Try to avoid sqrt as in Project method (w seems to need unitary axis12, but try it...)
              2) Hardcoded integer exponent 1, 2 or 3
            */
            Real length123( mal::Sqrt(length123_sq) );
            axis123 = axis123 / length123;
            //\todo ENFORCE axis123 orientation... it MAY be possible to guarantee it selecting proper x2,x3 for all possible x0,x1
            Real delta_h_alpha( degenerate_threshold_det_F * 6 * volume / length123 );
            Real w0( mal::Dot(x0-x2,axis123) ); //\note REQUIRES normalized axis //\note using x0 and x1 should have the SAME RESULT
            //Real w1( mal::Dot(x1-x2,axis123) ); //TEMP: unnecessary, w1 should be == w0 due to axis123 being orthogonal to x0-x1,
            /*\todo As axis123 may NOT be consistently
              oriented, we flip it ONLY if w and detF do not
              have the expected signs, that is, either x0 is
              BEHIND (x1,x2,x3), with detF < 0, and dot < 0,
              or IN FRONT, with 0 <= detF < DTDF and 0 <= dot.
            */
            if( w0 * detF < 0 )
            {
                w0 = -w0;
                axis123 = -axis123;
            }
            Real s( w0 / delta_h_alpha );
            Real t( Real(1) - s );
            // Cubic spline from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
            Real lambda_u = ( w0 > 0 ) //We know that w <= \Delta h as det_F <= m_rParams.m_DegenerateThresholdDetF
                            ? factor_L*delta_h_alpha*( (2-t)*t*t ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                            : factor_L*delta_h_alpha*( t + factor_NL*mal::Pow( -s, exponent_NL ) ); //lambda_u^L + lambda_u^NL
            x0 += lambda_u * axis123; //Move past collapse area, to ensure numerically safe
            x1 += lambda_u * axis123; //Move past collapse area, to ensure numerically safe
            // Compute Proj(F)
            Mat3x3 projDs( mal::GMat3x3_From_Columns( localX[1] - localX[0], localX[2] - localX[0], localX[3] - localX[0] ) );
            Mat3x3 projF( projDs * inv_Dm );
            Real det_projF( mal::Det(projF) );
            MS_ASSERT( det_projF >= 0 );
            if( det_projF < 0.99f*degenerate_threshold_det_F )
                MS_LOG_WARNING( "det(proj(F)) %f != %f should not happen if delta is correct", det_projF, degenerate_threshold_det_F );
            return mal::GRotation3x3_PolarDecomposition( projF, det_projF );
        }
        else
        {
            //\todo Consider using last R!!
            MS_LOG_ERROR( "Degenerate element with E-E DoC (%d,%d) has length123_sq = %f too small, collapsed? Using R=Id ", doc.m_VIT0, doc.m_VIT1, length123_sq );
            return Mat3x3::Identity();
        }
    }
}

Mat2x2 Compute_R_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                       int32 noc,
                       const Vec2& y0, const Vec2& y1, const Vec2& y2,
                       const Mat2x2& invDm, Real element_area )
{
    Mat2x2 Ds;
    TriangleElement2::Compute_D( y0, y1, y2, Ds );
    Mat2x2 F( Ds * invDm );
    Real det_F( mal::Det(F) );
    if( det_F > degenerate_threshold_det_F ) //Undegenerate
        return Compute_R_PD( F, det_F );
    else if( noc >= 0 && noc <= 2 ) //Degenerate and has proper NoC
        return Compute_R_DAPD_Degenerate( F, det_F, degenerate_threshold_det_F,
                                          noc,
                                          y0, y1, y2,
                                          element_area, invDm,
                                          factor_L, factor_NL, exponent_NL );
    else //Degenerate but no proper NoC
    {
        //\todo Consider using last R!!
        MS_LOG_ERROR( "Degenerate element has invalid NoC. Using R=Id" );
        return Mat2x2::Identity();
    }
}

Mat3x3 Compute_R_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                       TetrahedronElement3::DoC doc,
                       const Vec3& y0, const Vec3& y1, const Vec3& y2, const Vec3& y3,
                       const Mat3x3& invDm, Real element_volume )
{
    Mat3x3 Ds;
    TetrahedronElement3::Compute_D( y0, y1, y2, y3, Ds );
    Mat3x3 F( Ds * invDm );
    Real det_F( mal::Det(F) );
    if( det_F > degenerate_threshold_det_F ) //Undegenerate
        return Compute_R_PD( F, det_F );
    else if( doc.IsValid() ) //Degenerate and has proper DoC
        return Compute_R_DAPD_Degenerate( F, det_F, degenerate_threshold_det_F,
                                          doc,
                                          y0, y1, y2, y3,
                                          element_volume, invDm,
                                          factor_L, factor_NL, exponent_NL );
    else //Degenerate but no proper DoC
    {
        //\todo Consider using last R!!
        MS_LOG_ERROR( "Degenerate element has invalid DoC. Using R=Id" );
        return Mat3x3::Identity();
    }
}

////////////////////////////////////////////////////////////////
// Rotation Differential
// \note All methods returns false if dR is 0 (optimization)
//       or infinity (avoid div-by-zero)
////////////////////////////////////////////////////////////////

bool Compute_dR_QR_YX( const Vec2& dy0, const Vec2& dy1, const Vec2& dy2,
                       const Mat2x2& invDm, const Mat2x2& F,
                       Mat2x2& dR )
{
    Mat2x2 dDs;
    TriangleElement2::Compute_D( dy0, dy1, dy2, dDs );
    Mat2x2 dF( dDs * invDm );
    Vec2 f1 = mal::GColumn<1>(F);
    Vec2 df1 = mal::GColumn<1>(dF);
    Real norm_f1 = mal::Norm( f1 );
    if( norm_f1 > 10e-3f )
    {
        Vec2 de1 = ( mal::Sq(norm_f1)*df1 - mal::Dot(f1,df1)*f1 ) / (norm_f1 * norm_f1 * norm_f1);
        dR = mal::GMat2x2_From_Columns( mal::PerpendicularCCW(de1), de1 );
        return true;
    }
    else
        return false;
}

/*! Computes dR, from McAdams2011 tech doc
  \note Explicitly assumes S = R^T * F Symmetric

  \note Returns false if dR is 0 (optimization) or infinity (avoid
  div-by-zero)

  \todo This is dR of the Polar Decomposition, but is it aware of
  potential inversion?? McAdams2011 does NOT mention it... but it DOES
  WORK and results accurate with SVD1 inversion treatment, tested in
  testFE
*/
bool Compute_dR( const Vec2& dx0, const Vec2& dx1, const Vec2& dx2,
                 const Mat2x2& invDm, const Mat2x2& F, const Mat2x2& R, const Mat2x2& Rt,
                 Mat2x2& dR )
{
    //\todo In 2D this can be GREATLY OPTIMIZED by computing only the entries of W, S and dR strictly required
    Mat2x2 dDs;
    TriangleElement2::Compute_D( dx0, dx1, dx2, dDs ); //\todo Ds can be 0 if 2 or 3 dx are equal
    Mat2x2 dF( dDs * invDm );
    Mat2x2 S( Rt * F );
    Mat2x2 W( Rt * dF );
    Real w01_minus_w10( W(0,1) - W(1,0) );
    Real s00_plus_s11( S(0,0) + S(1,1) );
    //\note Check if ANY factor is 0, s00_plus_s11 to avoid division by 0, w01_minus_w10 to early-out if dR = 0
    //\todo The singularity check may be unnecessary, see 3D version
#ifdef __ENABLE_REGULARIZE_DAPD
    //s00_plus_s11 = mal::Max( s00_plus_s11, REGULARIZE_DAPD_THRESHOLD ); //\todo Consider adding, instead of maxing
    s00_plus_s11 += REGULARIZE_DAPD_THRESHOLD;
    if( true )
#else
    if( mal::Abs( w01_minus_w10 * s00_plus_s11 ) > 1e-6 )
#endif
    {
        Real skew_r( w01_minus_w10 / s00_plus_s11 );
        dR = R * Mat2x2( 0, skew_r, -skew_r, 0 ); //Skew2x2 has just 1 dof, skew_r
        MS_ASSERT( !mal::IsNaN( dR ) );
        return true;
    }
    else
        return false;
}

/*! Computes dR, from McAdams2011 tech doc
  \note Returns false if dR is 0 (optimization) or infinity (avoid div-by-zero)
*/
bool Compute_dR( const Vec3& dx0, const Vec3& dx1, const Vec3& dx2, const Vec3& dx3,
                 const Mat3x3& invDm, const Mat3x3& F, const Mat3x3& R, const Mat3x3& Rt,
                 Mat3x3& dR )
{
    //\todo This MAY BE GREATLY OPTIMIZED by computing only the entries of W, S and dR strictly required
    Mat3x3 dDs;
    TetrahedronElement3::Compute_D( dx0, dx1, dx2, dx3, dDs ); //\todo Ds can be 0 if 2, 3 or 4 dx are equal
    Mat3x3 dF( dDs * invDm );
    Mat3x3 S( Rt * F );
    Mat3x3 W( Rt * dF );
    Mat3x3 G( S - mal::Trace(S) * Mat3x3::Identity() );
    Vec3 w( W(1,2) - W(2,1),
            W(2,0) - W(0,2),
            W(0,1) - W(1,0) );
    //\todo The singularity may be unnecessary or too general, a singular G can only happen for a specific structure of S, which could be checked more efficiently or even impossible
    //\tood The 0-check in w can be done before checking if dDs is singular
    if( /*mal::NormSq(w) > mal::Epsilon<Real>()
          &&*/
        mal::Abs( mal::Det( G ) ) > mal::Epsilon<Real>() )
    {
        Vec3 r( mal::Inverse(G) * w );
        dR = R * mal::GMat3x3_From_Skew( r );
        return true;
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
bool Compute_dR_Numerical( const Vec2& x0, const Vec2& x1, const Vec2& x2,
                           const Vec2& dx0, const Vec2& dx1, const Vec2& dx2,
                           const Mat2x2& invDm, Real degenerate_threshold_det_F,
                           Real h,
                           Mat2x2& dR )
{
    mal::GVec<Real,6> v;
    mal::GSetRange<0,1>( v, dx0 );
    mal::GSetRange<2,3>( v, dx1 );
    mal::GSetRange<4,5>( v, dx2 );
    Real norm_sq_dx( mal::NormSq( v ) );
    if( norm_sq_dx > mal::Epsilon<Real>() )
    {
        Real norm_dx( mal::Sqrt(norm_sq_dx) );
        v = v / norm_dx;
        Vec2 v0 = mal::GRange<0,1>( v );
        Vec2 v1 = mal::GRange<2,3>( v );
        Vec2 v2 = mal::GRange<4,5>( v );
        Mat2x2 leftDs, leftF, leftR;
        Mat2x2 rightDs, rightF, rightR;
        TriangleElement2::Compute_D( x0 - h*v0, x1 - h*v1, x2 - h*v2, leftDs );
        TriangleElement2::Compute_D( x0 + h*v0, x1 + h*v1, x2 + h*v2, rightDs );
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
            MS_ASSERT( !mal::IsNaN( dR ) );
            return true;
        }
        else if( det_LF > degenerate_threshold_det_F )
        {
            leftR = mal::GRotation2x2_PolarDecomposition( leftF, det_LF );
            TriangleElement2::Compute_D( x0, x1, x2, rightDs );
            rightF = rightDs * invDm;
            rightR = mal::GRotation2x2_PolarDecomposition( rightF, mal::Det(rightF) );
            dR = norm_dx * mal::Rcp(h) * (rightR - leftR);
            MS_ASSERT( !mal::IsNaN( dR ) );
            return true;
        }
        else if( det_RF > degenerate_threshold_det_F )
        {
            TriangleElement2::Compute_D( x0, x1, x2, leftDs );
            leftF = leftDs * invDm;
            leftR = mal::GRotation2x2_PolarDecomposition( leftF, mal::Det(leftF) );
            rightR = mal::GRotation2x2_PolarDecomposition( rightF, det_RF );
            dR = norm_dx * mal::Rcp(h) * (rightR - leftR);
            MS_ASSERT( !mal::IsNaN( dR ) );
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

  \todo Not safe when if leftF or rightF are singular, consider using
  SVD instead of PD. Also, Compute_dR_Numerical() 2D version could be
  simplified using SVD to avoid left/right singularity treatment.
*/
bool Compute_dR_Numerical( const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3,
                           const Vec3& dx0, const Vec3& dx1, const Vec3& dx2, const Vec3& dx3,
                           const Mat3x3& invDm, Real degenerate_threshold_det_F,
                           Real h,
                           Mat3x3& dR )
{
    mal::GVec<Real,12> v;
    mal::GSetRange<0,2>( v, dx0 );
    mal::GSetRange<3,5>( v, dx1 );
    mal::GSetRange<6,8>( v, dx2 );
    mal::GSetRange<9,11>( v, dx3 );
    Real norm_sq_dx( mal::NormSq( v ) );
    if( norm_sq_dx > mal::Epsilon<Real>() )
    {
        Real norm_dx( mal::Sqrt(norm_sq_dx) );
        v = v / norm_dx;
        Vec3 v0 = mal::GRange<0,2>( v );
        Vec3 v1 = mal::GRange<3,5>( v );
        Vec3 v2 = mal::GRange<6,8>( v );
        Vec3 v3 = mal::GRange<9,11>( v );
        Mat3x3 leftDs, leftF, leftR;
        Mat3x3 rightDs, rightF, rightR;
        TetrahedronElement3::Compute_D( x0 - h*v0, x1 - h*v1, x2 - h*v2, x3 - h*v3, leftDs );
        TetrahedronElement3::Compute_D( x0 + h*v0, x1 + h*v1, x2 + h*v2, x3 + h*v3, rightDs );
        leftF = leftDs * invDm;
        rightF = rightDs * invDm;
        //\todo We should consider possible PD_P/R/U or SVD1 rotation corrections here ?!?!
        leftR = mal::GRotation3x3_PolarDecomposition( leftF, mal::Det(leftF) );
        rightR = mal::GRotation3x3_PolarDecomposition( rightF, mal::Det(rightF) );
        dR = norm_dx * mal::Rcp(2*h) * (rightR - leftR);
        return true;
    }
    else
        return false;
}

bool Compute_dR_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                      int32 noc,
                      const Vec2& y0, const Vec2& y1, const Vec2& y2,
                      const Vec2& dy0, const Vec2& dy1, const Vec2& dy2,
                      const Mat2x2& invDm, const Mat2x2& F, const Mat2x2& R, const Mat2x2& Rt, Real element_area,
                      Mat2x2& dR )
{
    if( mal::Det(F) > degenerate_threshold_det_F ) return Compute_dR( dy0, dy1, dy2, invDm, F, R, Rt, dR );
    if( noc < 0 || noc > 2 ) return false;

    // use noc and get x0, x1, x2 accordingly
    Vec2 local_X[3] = { y0, y1, y2 };
    Vec2 &x0_bar( local_X[ noc ] ); //Inverted vtx
    const Vec2 &x0( local_X[ noc ] ); //Inverted vtx
    const Vec2 &x1( local_X[ (noc + 1) % 3 ] );
    const Vec2 &x2( local_X[ (noc + 2) % 3 ] );

    Vec2 local_dX[3] = { dy0, dy1, dy2 };
    Vec2 &dx0_bar( local_dX[ noc ] ); //Inverted vtx
    const Vec2 &dx0( local_dX[ noc ] ); //Inverted vtx
    const Vec2 &dx1( local_dX[ (noc + 1) % 3 ] );
    const Vec2 &dx2( local_dX[ (noc + 2) % 3 ] );

    Vec2 x12( x2 - x1 );
    Real dist12_sq( mal::NormSq(x12) );
#ifdef __ENABLE_REGULARIZE_DAPD
    //dist12_sq = mal::Max( dist12_sq, mal::Sq(REGULARIZE_DAPD_THRESHOLD) );
    dist12_sq += mal::Sq(REGULARIZE_DAPD_THRESHOLD);
    if( true )
#else
    if( dist12_sq > 0.000001f ) //\todo CONSIDER regularization instead of an IF that induces a "hard" discontinuity
#endif
    {
        // R(lambda) part \todo REPEATED COMPUTATION from DAPD, consider saving in ElementCache. \todo THIS IS specially serious in 3D, where dir_c sign needs to be computed, at the very least cache dir_c and dist12
        Real dist12( mal::Sqrt(dist12_sq) );
        Vec2 dir_c( mal::PerpendicularCW( x12 ) / dist12 );
        Real h( degenerate_threshold_det_F * 2 * element_area / dist12 ); //\todo This is delta_h_alpha in other code, and \Delta h in notebooks/LaTeX
        //\todo CONSIDER h regularization due to dist12
        Vec2 x12_mid( 0.5 * (x1 + x2) ); //unnecessary, added for symmetry with dx
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
                      ? ( ( w >= 0 ) ? factor_L*h*( 2*t2 - t3 ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                          : factor_L*h*( t + factor_NL*mal::Pow( s, exponent_NL ) ) ) //lambda_u^L + lambda_u^NL
                      : 0;

        // dR(lambda) part
        Vec2 dx12( dx2 - dx1 );
        Vec2 ddir_c( mal::PerpendicularCW( dx12*dist12_sq - x12*mal::Dot(x12,dx12) ) / (dist12*dist12_sq) );
        Vec2 dx12_mid( 0.5 * (dx1 + dx2) ); //midpoint dx, to avoid bias towards either x1 or x2
        Real dw( mal::Dot( dx0 - dx12_mid, dir_c ) + mal::Dot( x0 - x12_mid, ddir_c ) );
        Real Pl_Pw = ( w < h )
                     ? ( ( w >= 0 ) ? -factor_L*( 4*t - 3*t2 ) //dlambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                         : -factor_L*( 1 + factor_NL*exponent_NL*mal::Pow( s, exponent_NL-1 ) ) ) //dlambda_u^L + dlambda_u^NL )
                     : 0;
        Real dlambda_w( Pl_Pw * dw );

        Real dh( -degenerate_threshold_det_F * 2 * element_area * mal::Dot(x12,dx12) / (dist12*dist12_sq) );
        Real Pl_Ph = ( w < h )
                     ? ( ( w >= 0 ) ? factor_L * ( 1 + s2 + 2*s3  )
                         : factor_L*( 1 - factor_NL*(exponent_NL-1)*mal::Pow( s, exponent_NL ) ) )
                     : 0;
        Real dlambda_h( Pl_Ph * dh );

        // Corrected x and dx
        x0_bar = x0 + lambda * dir_c;
        dx0_bar = dx0 + dlambda_w * dir_c + lambda * ddir_c + dlambda_h * dir_c;

        // Compute corrected \bar Ds, \bar F and their differentials
        Mat2x2 Ds_bar;
        TriangleElement2::Compute_D( local_X[0], local_X[1], local_X[2], Ds_bar );
        Mat2x2 F_bar( Ds_bar * invDm );
        Mat2x2 dDs_bar;
        TriangleElement2::Compute_D( local_dX[0], local_dX[1], local_dX[2], dDs_bar );
        Mat2x2 dF_bar( dDs_bar * invDm );

        // Compute \bar dR
        //\todo In 2D this can be GREATLY OPTIMIZED by computing only the entries of W, S and dR strictly required
        Mat2x2 S( Rt * F_bar ); //\todo This S is symmetric, as Rt comes from F_bar polar decomposition
        Mat2x2 W( Rt * dF_bar );
        Real w01_minus_w10( W(0,1) - W(1,0) );
        Real s00_plus_s11( S(0,0) + S(1,1) );
        //\note Check if ANY factor is 0, s00_plus_s11 to avoid division by 0, w01_minus_w10 to early-out if dR = 0
        //\todo The singularity check may be unnecessary, see 3D version
#ifdef __ENABLE_REGULARIZE_DAPD
        //s00_plus_s11 = mal::Max( s00_plus_s11, REGULARIZE_DAPD_THRESHOLD ); //\todo Consider adding, instead of maxing
        s00_plus_s11 += REGULARIZE_DAPD_THRESHOLD;
        if( true )
#else
        if( mal::Abs( w01_minus_w10 * s00_plus_s11 ) > mal::Sq( mal::Epsilon<Real>() ) ) //1e-6 )
#endif
        {
            Real skew_r( w01_minus_w10 / s00_plus_s11 );
            dR = R * Mat2x2( 0, skew_r, -skew_r, 0 ); //Skew2x2 has just 1 dof, skew_r
            MS_ASSERT( !mal::IsNaN( dR ) );
            return true;
        }
        else
            return false;
    }
    else
        return false;
}

// Local helper function
inline void Compute_lambda_dlambda_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                                         Real element_volume,
                                         Real w, Real h,
                                         const Vec3& x0, const Vec3& dx0,
                                         const Vec3& x_mid, const Vec3& dx_mid,
                                         Real norm_u, Real norm_u_sq,
                                         const Vec3& u, const Vec3& du,
                                         const Vec3& dir_c, const Vec3& ddir_c,
                                         Real& lambda, Real& dlambda )
{
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

    lambda = ( w < h )
             ? ( ( w >= 0 ) ? factor_L*h*( 2*t2 - t3 ) //lambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                 : factor_L*h*( t + factor_NL*mal::Pow( s, exponent_NL ) ) ) //lambda_u^L + lambda_u^NL
             : 0;

    Real dw( mal::Dot( dx0 - dx_mid, dir_c ) + mal::Dot( x0 - x_mid, ddir_c ) );
    Real Pl_Pw = ( w < h )
                 ? ( ( w >= 0 ) ? -factor_L*( 4*t - 3*t2 ) //dlambda_u^D, Cubic spline with P0 = 0, T0 = 0, P1 = \tau \Delta h, T1 = \tau \Delta h
                     : -factor_L*( 1 + factor_NL*exponent_NL*mal::Pow( s, exponent_NL-1 ) ) ) //dlambda_u^L + dlambda_u^NL )
                 : 0;
    Real dlambda_w( Pl_Pw * dw );

    Real dh( - degenerate_threshold_det_F * 6 * element_volume * mal::Dot(u,du) / (norm_u_sq*norm_u) );
    Real Pl_Ph = ( w < h )
                 ? ( ( w >= 0 ) ? factor_L * ( 1 + s2 + 2*s3  )
                     : factor_L*( 1 - factor_NL*(exponent_NL-1)*mal::Pow( s, exponent_NL ) ) )
                 : 0;
    Real dlambda_h( Pl_Ph * dh );
    dlambda = dlambda_w + dlambda_h;
}

bool Compute_dR_DAPD( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                      TetrahedronElement3::DoC doc,
                      const Vec3& y0, const Vec3& y1, const Vec3& y2, const Vec3& y3,
                      const Vec3& dy0, const Vec3& dy1, const Vec3& dy2, const Vec3& dy3,
                      const Mat3x3& invDm, const Mat3x3& F, const Mat3x3& R, const Mat3x3& Rt, Real element_volume,
                      Mat3x3& dR )
{
    if( mal::Det(F) > degenerate_threshold_det_F ) return Compute_dR( dy0, dy1, dy2, dy3, invDm, F, R, Rt, dR );
    if( !doc.IsValid() ) return false;

    // use noc and get x0, x1, x2 accordingly
    Vec3 local_X[4] = { y0, y1, y2, y3 };
    Vec3 local_dX[4] = { dy0, dy1, dy2, dy3 };

    if( doc.IsVF() )
    {
        int noc( doc.m_VIT0 );
        Vec3 &x0_bar( local_X[ noc ] ); //Inverted vtx
        const Vec3 &x0( local_X[ noc ] ); //Inverted vtx
        const Vec3 &x1( local_X[ (noc + 1) % 4 ] );
        const Vec3 &x2( local_X[ (noc + 2) % 4 ] );
        const Vec3 &x3( local_X[ (noc + 3) % 4 ] );

        Vec3 &dx0_bar( local_dX[ noc ] ); //Inverted vtx
        const Vec3 &dx0( local_dX[ noc ] ); //Inverted vtx
        const Vec3 &dx1( local_dX[ (noc + 1) % 4 ] );
        const Vec3 &dx2( local_dX[ (noc + 2) % 4 ] );
        const Vec3 &dx3( local_dX[ (noc + 3) % 4 ] );

        Vec3 u( mal::Cross(x2-x1, x3-x1) );
        Real norm_u_sq( mal::NormSq(u) ); //\note This is the term B_{234}^2 in DCLFEM paper
#ifdef __ENABLE_REGULARIZE_DAPD
        norm_u_sq += mal::Sq(REGULARIZE_DAPD_THRESHOLD);
        if( true )
#else
        if( norm_u_sq > 0.000001f ) //\todo CONSIDER regularization instead of an IF that induces a "hard" discontinuity
#endif
        {
           /*\todo OPTIMIZE:
              0) R(lambda) REPEATED COMPUTATION from DAPD, consider
                 saving in ElementCache. THIS IS specially serious in
                 3D, where dir_c sign needs to be computed, at the
                 very least cache dir_c and dist12
              1) Try to avoid sqrt as in Project method (w seems to
                 need unitary dir_c, but try it...)
              2) Hardcoded integer exponent 1, 2 or 3
              3) ENFORCE u orientation... it MAY be possible to
                 guarantee it selecting proper x2,x3 for all possible
                 x0,x1
            */
            Real norm_u( mal::Sqrt(norm_u_sq) );
            Vec3 dir_c( u / norm_u );
            Real h( degenerate_threshold_det_F * 6 * element_volume / norm_u );
            Vec3 x123_mid( (1.0/3.0) * (x1 + x2 + x3) ); //unnecessary, added for symmetry with dx
            Real w( mal::Dot( x0 - x123_mid, dir_c ) );
            /*\todo As u / dir_c may NOT be consistently
              oriented, we flip it ONLY if w and detF do not
              have the expected signs, that is, either x0 is
              BEHIND (x1,x2,x3), with detF < 0, and dot < 0,
              or IN FRONT, with 0 <= detF < DTDF and 0 <= dot.
            */
            if( w * mal::Det(F) < 0 )
            {
                w = -w;
                u = -u;
                dir_c = -dir_c;
            }

            // du and ddir_c
            Vec3 dx123_mid( (1.0/3.0) * (dx1 + dx2 + dx3) ); //midpoint dx, to avoid bias
            Vec3 x12( x2 - x1 );
            Vec3 x13( x3 - x1 );
            Vec3 dx12( dx2 - dx1 );
            Vec3 dx13( dx3 - dx1 );
            Vec3 du( mal::Cross(x12,dx13) + mal::Cross(dx12,x13) );

            // ddir_c \todo generic VF and EE
            Vec3 ddir_c( (du*norm_u_sq - u*mal::Dot(u,du) ) / (norm_u_sq*norm_u) );

            // Compute lambda and dlambda \todo generic VF and EE
            Real lambda, dlambda;
            Compute_lambda_dlambda_DAPD( degenerate_threshold_det_F, factor_L, factor_NL, exponent_NL,
                                         element_volume,
                                         w, h,
                                         x0, dx0,
                                         x123_mid, dx123_mid,
                                         norm_u, norm_u_sq,
                                         u, du,
                                         dir_c, ddir_c,
                                         lambda, dlambda );

            // Corrected x0 and dx0
            x0_bar = x0 + lambda * dir_c;
            dx0_bar = dx0 + dlambda * dir_c + lambda * ddir_c;
        }
        else
            return false;
    }
    else //doc.IsEE()
    {
        //\todo Probably only w and dir_c terms are specific to VF/EE,
        //consider extracting common code (such as dlambda
        //computation)

        int vit0( doc.m_VIT0 );
        MS_ASSERT( 0 == vit0 );
        int vit1( doc.m_VIT1 );
        //\todo FIND A BETTER WAY to compute vid2,vid3...
        int vit2, vit3;
        switch( vit1 )
        {
        case 1: vit2 = 2; vit3 = 3; break;
        case 2: vit2 = 1; vit3 = 3; break;
        case 3: vit2 = 1; vit3 = 2; break;
        default: vit2 = vit3 = 1000; break;
        }
        MS_ASSERT( vit0 + vit1 + vit2 + vit3 == 0+1+2+3 ); //TEMP

        Vec3 &x0_bar( local_X[ vit0 ] );
        Vec3 &x1_bar( local_X[ vit1 ] );
        const Vec3 &x0( local_X[ vit0 ] );
        const Vec3 &x1( local_X[ vit1 ] );
        const Vec3 &x2( local_X[ vit2 ] );
        const Vec3 &x3( local_X[ vit3 ] );

        Vec3 &dx0_bar( local_dX[ vit0 ] );
        Vec3 &dx1_bar( local_dX[ vit1 ] );
        const Vec3 &dx0( local_dX[ vit0 ] );
        const Vec3 &dx1( local_dX[ vit1 ] );
        const Vec3 &dx2( local_dX[ vit2 ] );
        const Vec3 &dx3( local_dX[ vit3 ] );

        Vec3 u( mal::Cross(x1-x0, x3-x2) ); //e_ij x e_kl
        Real norm_u_sq( mal::NormSq(u) ); //\note This is the term B_{234}^2 in DCLFEM paper
#ifdef __ENABLE_REGULARIZE_DAPD
        norm_u_sq += mal::Sq(REGULARIZE_DAPD_THRESHOLD);
        if( true )
#else
        if( norm_u_sq > 0.000001f ) //\todo CONSIDER regularization instead of an IF that induces a "hard" discontinuity
#endif
        {
           /*\todo OPTIMIZE:
              0) R(lambda) REPEATED COMPUTATION from DAPD, consider
                 saving in ElementCache. THIS IS specially serious in
                 3D, where dir_c sign needs to be computed, at the
                 very least cache dir_c and dist12
              1) Try to avoid sqrt as in Project method (w seems to
                 need unitary dir_c, but try it...)
              2) Hardcoded integer exponent 1, 2 or 3
              3) ENFORCE u orientation... it MAY be possible to
                 guarantee it selecting proper x2,x3 for all possible
                 x0,x1
            */
            Real norm_u( mal::Sqrt(norm_u_sq) );
            Vec3 dir_c( u / norm_u );
            Real h( degenerate_threshold_det_F * 6 * element_volume / norm_u );
            Vec3 x23_mid( (1.0/2.0) * (x2 + x3) ); //unnecessary, added for symmetry with dx
            Real w( mal::Dot( x0-x23_mid, dir_c ) ); //\note using x1 the result should be the same
            /*\todo As axis123 may NOT be consistently
              oriented, we flip it ONLY if w and detF do not
              have the expected signs, that is, either x0 is
              BEHIND (x1,x2,x3), with detF < 0, and dot < 0,
              or IN FRONT, with 0 <= detF < DTDF and 0 <= dot.
            */
            if( w * mal::Det(F) < 0 )
            {
                w = -w;
                u = -u;
                dir_c = -dir_c;
            }

            // du and ddir_c
            Vec3 dx23_mid( (1.0/2.0) * (dx2 + dx3) ); //midpoint dx, to avoid bias
            Vec3 x01( x1 - x0 );
            Vec3 x23( x3 - x2 );
            Vec3 dx01( dx1 - dx0 );
            Vec3 dx23( dx3 - dx2 );

            Vec3 du( mal::Cross(x01,dx23) + mal::Cross(dx01,x23) );

            // ddir_c \todo generic VF and EE
            Vec3 ddir_c( (du*norm_u_sq - u*mal::Dot(u,du) ) / (norm_u_sq*norm_u) );

            // Compute lambda and dlambda \todo generic VF and EE
            Real lambda, dlambda;
            Compute_lambda_dlambda_DAPD( degenerate_threshold_det_F, factor_L, factor_NL, exponent_NL,
                                         element_volume,
                                         w, h,
                                         x0, dx0,
                                         x23_mid, dx23_mid,
                                         norm_u, norm_u_sq,
                                         u, du,
                                         dir_c, ddir_c,
                                         lambda, dlambda );

            // Corrected x0,x1 and dx0,dx1
            x0_bar = x0 + lambda * dir_c;
            x1_bar = x1 + lambda * dir_c;
            dx0_bar = dx0 + dlambda * dir_c + lambda * ddir_c;
            dx1_bar = dx1 + dlambda * dir_c + lambda * ddir_c;
        }
        else
            return false;
    }

    // Compute \bar dR from modified X and dX => \bar Ds, \bar F and their differentials
    Mat3x3 Ds_bar;
    TetrahedronElement3::Compute_D( local_X[0], local_X[1], local_X[2], local_X[3], Ds_bar );
    Mat3x3 F_bar( Ds_bar * invDm );
    Mat3x3 dDs_bar;
    TetrahedronElement3::Compute_D( local_dX[0], local_dX[1], local_dX[2], local_dX[3], dDs_bar );
    Mat3x3 dF_bar( dDs_bar * invDm );

    // Compute dR
    Mat3x3 S( Rt * F_bar ); //\todo This S is symmetric, as Rt comes from F_bar polar decomposition
    Mat3x3 W( Rt * dF_bar );
    Mat3x3 G( S - mal::Trace(S) * Mat3x3::Identity() );
    Vec3 w( W(1,2) - W(2,1),
            W(2,0) - W(0,2),
            W(0,1) - W(1,0) );
    //\todo The singularity may be unnecessary or too general, a singular G can only happen for a specific structure of S, which could be checked more efficiently or even impossible
    //\tood The 0-check in w can be done before checking if dDs is singular
    if( /*mal::NormSq(w) > mal::Epsilon<Real>()
          &&*/
        mal::Abs( mal::Det( G ) ) > mal::Epsilon<Real>() )
    {
        Vec3 r( mal::Inverse(G) * w );
        dR = R * mal::GMat3x3_From_Skew( r );
        return true;
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

  //\todo Generic compute R according to selected method
*/
bool Compute_dR_DAPD_Numerical( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                                int32 noc,
                                const Vec2& x0, const Vec2& x1, const Vec2& x2,
                                const Vec2& dx0, const Vec2& dx1, const Vec2& dx2,
                                const Mat2x2& invDm, Real element_area,
                                Real h,
                                Mat2x2& dR )
{
    mal::GVec<Real,6> v;
    mal::GSetRange<0,1>( v, dx0 );
    mal::GSetRange<2,3>( v, dx1 );
    mal::GSetRange<4,5>( v, dx2 );
    Real norm_sq_dx( mal::NormSq( v ) );
#ifdef __ENABLE_REGULARIZE_DAPD
    //norm_sq_dx = mal::Max( norm_sq_dx, mal::Sq(REGULARIZE_DAPD_THRESHOLD) );
    norm_sq_dx += mal::Sq(REGULARIZE_DAPD_THRESHOLD);
    if( true )
#else
    if( norm_sq_dx > mal::Epsilon<Real>() )
#endif
    {
        Real norm_dx( mal::Sqrt(norm_sq_dx) );
        v = v / norm_dx;
        Vec2 v0 = mal::GRange<0,1>( v );
        Vec2 v1 = mal::GRange<2,3>( v );
        Vec2 v2 = mal::GRange<4,5>( v );
        Vec2 leftX[] = { x0 - h*v0, x1 - h*v1, x2 - h*v2 };
        Vec2 rightX[] = { x0 + h*v0, x1 + h*v1, x2 + h*v2 };
        Mat2x2 leftR = Compute_R_DAPD( degenerate_threshold_det_F, factor_L, factor_NL, exponent_NL,
                                       noc,
                                       leftX[0], leftX[1], leftX[2],
                                       invDm, element_area );
        Mat2x2 rightR = Compute_R_DAPD( degenerate_threshold_det_F, factor_L, factor_NL, exponent_NL,
                                        noc,
                                        rightX[0], rightX[1], rightX[2],
                                        invDm, element_area );
        dR = norm_dx * mal::Rcp(2*h) * (rightR - leftR);
        MS_ASSERT( !mal::IsNaN( dR ) );
        return true;
    }
    else
        return false;
}

bool Compute_dR_DAPD_Numerical( Real degenerate_threshold_det_F, Real factor_L, Real factor_NL, Real exponent_NL,
                                TetrahedronElement3::DoC doc,
                                const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& x3,
                                const Vec3& dx0, const Vec3& dx1, const Vec3& dx2, const Vec3& dx3,
                                const Mat3x3& invDm, Real element_volume,
                                Real h,
                                Mat3x3& dR )
{
    mal::GVec<Real,12> v;
    mal::GSetRange<0,2>( v, dx0 );
    mal::GSetRange<3,5>( v, dx1 );
    mal::GSetRange<6,8>( v, dx2 );
    mal::GSetRange<9,11>( v, dx3 );
    Real norm_sq_dx( mal::NormSq( v ) );
#ifdef __ENABLE_REGULARIZE_DAPD
    //norm_sq_dx = mal::Max( norm_sq_dx, mal::Sq(REGULARIZE_DAPD_THRESHOLD) );
    norm_sq_dx += mal::Sq(REGULARIZE_DAPD_THRESHOLD);
    if( true )
#else
    if( norm_sq_dx > mal::Epsilon<Real>() )
#endif
    {
        Real norm_dx( mal::Sqrt(norm_sq_dx) );
        v = v / norm_dx;
        Vec3 v0 = mal::GRange<0,2>( v );
        Vec3 v1 = mal::GRange<3,5>( v );
        Vec3 v2 = mal::GRange<6,8>( v );
        Vec3 v3 = mal::GRange<9,11>( v );
        Vec3 leftX[] = { x0 - h*v0, x1 - h*v1, x2 - h*v2, x3 - h*v3 };
        Vec3 rightX[] = { x0 + h*v0, x1 + h*v1, x2 + h*v2, x3 + h*v3 };
        Mat3x3 leftR = Compute_R_DAPD( degenerate_threshold_det_F, factor_L, factor_NL, exponent_NL,
                                       doc,
                                       leftX[0], leftX[1], leftX[2], leftX[3],
                                       invDm, element_volume );
        Mat3x3 rightR = Compute_R_DAPD( degenerate_threshold_det_F, factor_L, factor_NL, exponent_NL,
                                        doc,
                                        rightX[0], rightX[1], rightX[2], rightX[3],
                                        invDm, element_volume );
        dR = norm_dx * mal::Rcp(2*h) * (rightR - leftR);
        MS_ASSERT( !mal::IsNaN( dR ) );
        return true;
    }
    else
        return false;
}

#ifdef __USE_SVD_FROM_GSL_DEPRECATED
template <typename T>
inline void GSingularValueDecomposition_USVt_GSL( const GMat<T,3,3> &F,
                                                  GMat<T,3,3> &U, GVec<T,3> &diag_F, GMat<T,3,3> &Vt )
{
    //\see http://www.gnu.org/software/gsl/manual/html_node/Singular-Value-Decomposition.html
    //\note SLOW! SHOULD AVOID ALLOC/FREE per call!!
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
}
#endif

}}} //namespace S2::ms::fem
