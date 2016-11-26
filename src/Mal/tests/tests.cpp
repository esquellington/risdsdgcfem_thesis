//Test compilation with all headers included

#include <Mal/Config.h>
#include <Mal/NumberUtils.h>
#include <Mal/RealUtils.h>
#include <Mal/GFxp.h>
#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GQuatUtils.h>
#include <Mal/GPosRot.h>
#include <Mal/GMat.h>
#include <Mal/GMatUtils.h>
#include <Mal/GMatDecomposition.h>
#include <Mal/GMatEigen.h>
#include <Mal/GTransform.h>
#include <Mal/GTensor3.h>
#include <Mal/GSerialization.h>
#include <Mal/GConversion.h>
#include <Mal/GInterval.h>
#include <Mal/GRandom.h>
#include <Mal/GSRV.h>
#include "tests.h"
#include <stdio.h>
#include <iostream>

namespace mal
{

bool TestAllHeaders()
{
    return true;
}

bool TestTemplateInstantiation()
{
    typedef GVec<float,4> Vec4f;
    Vec4f vec4;

    typedef GQuat<double> Quatd;
    Quatd quatd;

    typedef GTransform<float,3> Transform3f;
    Transform3f tans3f;

    typedef GPosRot<float> PosRot3f;
    PosRot3f posrot3f;

    typedef GMat<double,3,5> Mat3x5d;
    Mat3x5d mat3x5d;
    mat3x5d(1,1) = 0;

    return true;
}

bool TestTransform()
{
    typedef GTransform<float,2> Transform2f;
    typedef GMat<float,2,2> Mat2x2f;
    typedef GVec<float,2> Vec2f;
    mal::SetRandomSeed(666);
    Vec2f a( RandomV( Vec2f(-1,-1), Vec2f(1,1) ) );
    Vec2f b( RandomV( Vec2f(-1,-1), Vec2f(1,1) ) );
    Vec2f c( RandomV( Vec2f(-1,-1), Vec2f(1,1) ) );
    Transform2f tr( Vec2f(1,1), GRotation2x2_From( 1.333f ) );
    Mat2x2f m( 1, 2, 3, 4 );
    Vec2f v1( tr * (m*(a-b) + c) );
    Vec2f v2( tr.m_Pos + tr.m_Rot*m*(a-b) + tr.m_Rot * c );
    std::cout << v1 << " =?= " << v2 << std::endl;
    return ApproxEq( NormSq( v1 - v2 ), 0.0f );
}

#define TEST_RADIX 15
bool TestQuatPrecInFxp()
{
    // do it in doubles
    typedef GQuat<double> Quatd;
    typedef GVec<double,3> Vec3d;

    Quatd quatd(0.333,0.333,0.333,0.333);
    quatd.Normalize();
    Vec3d vec(0.666,0.666,0.666);
    Vec3d vec_transformed = quatd * vec;

    // do it in Fxp
    typedef GQuat< GFxp<TEST_RADIX> > Quatr;
    typedef GVec< GFxp<TEST_RADIX> ,3> Vec3x;

    Quatr quat0;
    quat0.x = quatd.x;
    quat0.y = quatd.y;
    quat0.z = quatd.z;
    quat0.w = quatd.w;
    //quat0.Normalize();

    Vec3x vec0;
    vec0.x() = vec.x();
    vec0.y() = vec.y();
    vec0.z() = vec.z();
    Vec3x vec_transformed0 = quat0 * vec0;

    // do it in Fxp assuming 1 implicit component
    Quatr quatX, quatY, quatZ, quatW;
    quatX = quatY = quatZ = quatW = quat0;

    quatX.x = Sqrt( GFxp<TEST_RADIX>(1) - (quat0.y*quat0.y + quat0.z*quat0.z + quat0.w*quat0.w) );
    quatY.y = Sqrt( GFxp<TEST_RADIX>(1) - (quat0.x*quat0.x + quat0.z*quat0.z + quat0.w*quat0.w) );
    quatZ.z = Sqrt( GFxp<TEST_RADIX>(1) - (quat0.x*quat0.x + quat0.y*quat0.y + quat0.w*quat0.w) );
    quatW.w = Sqrt( GFxp<TEST_RADIX>(1) - (quat0.x*quat0.x + quat0.y*quat0.y + quat0.z*quat0.z) );

    Vec3x vec_transformedX = quatX * vec0;
    Vec3x vec_transformedY = quatY * vec0;
    Vec3x vec_transformedZ = quatZ * vec0;
    Vec3x vec_transformedW = quatW * vec0;

    //Check result accuracy when compared to results using double
    Vec3d vec_transformed_in_fxp0;
    vec_transformed_in_fxp0.x() = vec_transformed0.x().GetDouble();
    vec_transformed_in_fxp0.y() = vec_transformed0.y().GetDouble();
    vec_transformed_in_fxp0.z() = vec_transformed0.z().GetDouble();
    Vec3d vec_transformed_in_fxpX;
    vec_transformed_in_fxpX.x() = vec_transformedX.x().GetDouble();
    vec_transformed_in_fxpX.y() = vec_transformedX.y().GetDouble();
    vec_transformed_in_fxpX.z() = vec_transformedX.z().GetDouble();
    Vec3d vec_transformed_in_fxpY;
    vec_transformed_in_fxpY.x() = vec_transformedY.x().GetDouble();
    vec_transformed_in_fxpY.y() = vec_transformedY.y().GetDouble();
    vec_transformed_in_fxpY.z() = vec_transformedY.z().GetDouble();
    Vec3d vec_transformed_in_fxpZ;
    vec_transformed_in_fxpZ.x() = vec_transformedZ.x().GetDouble();
    vec_transformed_in_fxpZ.y() = vec_transformedZ.y().GetDouble();
    vec_transformed_in_fxpZ.z() = vec_transformedZ.z().GetDouble();
    Vec3d vec_transformed_in_fxpW;
    vec_transformed_in_fxpW.x() = vec_transformedW.x().GetDouble();
    vec_transformed_in_fxpW.y() = vec_transformedW.y().GetDouble();
    vec_transformed_in_fxpW.z() = vec_transformedW.z().GetDouble();

    vec_transformed_in_fxpX -= vec_transformed;
    vec_transformed_in_fxpY -= vec_transformed;
    vec_transformed_in_fxpZ -= vec_transformed;
    vec_transformed_in_fxpW -= vec_transformed;
    vec_transformed_in_fxp0 -= vec_transformed;

    double error0 = vec_transformed_in_fxp0.Norm();
    double errorX = vec_transformed_in_fxpX.Norm();
    double errorY = vec_transformed_in_fxpY.Norm();
    double errorZ = vec_transformed_in_fxpZ.Norm();
    double errorW = vec_transformed_in_fxpW.Norm();

    return true;
}

//Testing mid,half interval representation
template <typename T>
class GIntervalMH
{
public:
    finline GIntervalMH() : m_Mid(0), m_Half(0) {}
    finline GIntervalMH( T min, T max ) : m_Mid( T(0.5)*(min+max) ), m_Half( T(0.5)*(max-min) ) {}
    finline bool TestOverlapAbs( const GIntervalMH &other ) const
        {
            return mal::Abs( m_Mid - other.m_Mid ) <= ( m_Half + other.m_Half );
        }
private:
    T m_Mid;
    T m_Half;
};

bool TestInterval()
{
    typedef GInterval<float> interval_type;
    //typedef GIntervalMH<float> interval_type;
    unsigned int num_errors(0);
    mal::SetRandomSeed(666);

    const unsigned int cNumSamples(1000);
    interval_type vec_interval[cNumSamples];
    for( unsigned int i=0; i<cNumSamples; i++ )
    {
        float m0( mal::RandomF<float>(-1,1) );
        float M0( mal::RandomF<float>(m0,1) );
        vec_interval[i] = interval_type( m0, M0 );
    }
    unsigned int num_overlaps(0);
    for( unsigned int i=0; i<cNumSamples*cNumSamples; i++ )
        //num_overlaps += vec_interval[i%cNumSamples].TestOverlap( vec_interval[(i*i+666)%cNumSamples] ) ? 1 : 0; //fastest by almost 2x!
        num_overlaps += vec_interval[i%cNumSamples].TestOverlap( vec_interval[(i*i+666)%cNumSamples] ) ? 1 : 0; //fastest by almost 2x!
        //num_overlaps += vec_interval[i%cNumSamples].TestOverlapCmp( vec_interval[(i*i+666)%cNumSamples] ) ? 1 : 0;
        //num_overlaps += vec_interval[i%cNumSamples].TestOverlapAbs( vec_interval[(i*i+666)%cNumSamples] ) ? 1 : 0;
    printf( "TestInterval() num_overlaps = %d\n", num_overlaps );

    /* Test correctness OverlapsA, OverlapsB, OverlapsC
    for( unsigned int i=0; i<1000000; i++ )
    {
        float m0( mal::RandomF<float>(-1,1) );
        float M0( mal::RandomF<float>(m0,1) );
        float m1( mal::RandomF<float>(-1,1) );
        float M1( mal::RandomF<float>(m1,1) );
        interval_type interval0( m0, M0 );
        interval_type interval1( m1, M1 );
        bool bOverlapsA( interval0.OverlapsA(interval1) );
        bool bOverlapsB( interval0.OverlapsB(interval1) );
        bool bOverlapsC( interval0.OverlapsC(interval1) );
        //MAL_LOG( "TestInterval() intervals [%f,%f] and [%f,%f] = (%d,%d,%d)",
        //         m0, M0, m1, M1, (int)bOverlapsA, (int)bOverlapsB, (int)bOverlapsC );
        if( bOverlapsA != bOverlapsB || bOverlapsB != bOverlapsC ) num_errors++;
    }
    printf( "TestInterval() num_errors = %d\n", num_errors );
    */

    return num_errors == 0;
}

} //namespace mal
