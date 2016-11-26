#ifndef GEO_BV_GKDOP_AXIS_H
#define GEO_BV_GKDOP_AXIS_H

#include <Geo/Config.h>

namespace geo {
namespace bv {

//\todo Consider external GKDOP_Axis.h file and put ALL this boring methods there
//---- GKDOP Dot product with static axis
template<unsigned D, unsigned K, unsigned A> Real GKDOP_Dot( const mal::GVec<Real,D>& v ) { return mal::Infinity<Real>(); }
//DOP2_K4 (X,Y)
template<> finline Real GKDOP_Dot<2,4,0>( const Vec2& v ) { return v[0]; } //+X
template<> finline Real GKDOP_Dot<2,4,1>( const Vec2& v ) { return v[1]; } //+Y
//DOP2_K8 (X,Y,XY,Xy)
template<> finline Real GKDOP_Dot<2,8,0>( const Vec2& v ) { return v[0]; } //+X
template<> finline Real GKDOP_Dot<2,8,1>( const Vec2& v ) { return v[1]; } //+Y
template<> finline Real GKDOP_Dot<2,8,2>( const Vec2& v ) { return mal::RcpSqrtOfTwo<Real>()*(v[0] + v[1]); } //+X+Y
template<> finline Real GKDOP_Dot<2,8,3>( const Vec2& v ) { return mal::RcpSqrtOfTwo<Real>()*(v[0] - v[1]); } //+X-Y
//DOP3_K6 (X,Y,Z)
template<> finline Real GKDOP_Dot<3,6,0>( const Vec3& v ) { return v[0]; } //+X
template<> finline Real GKDOP_Dot<3,6,1>( const Vec3& v ) { return v[1]; } //+Y
template<> finline Real GKDOP_Dot<3,6,2>( const Vec3& v ) { return v[2]; } //+Z
//DOP3_K14 (X,Y,Z,XYZ,XYz,XyZ,Xyz)
template<> finline Real GKDOP_Dot<3,14,0>( const Vec3& v ) { return v[0]; } //+X
template<> finline Real GKDOP_Dot<3,14,1>( const Vec3& v ) { return v[1]; } //+Y
template<> finline Real GKDOP_Dot<3,14,2>( const Vec3& v ) { return v[2]; } //+Z
template<> finline Real GKDOP_Dot<3,14,3>( const Vec3& v ) { return mal::RcpSqrtOfThree<Real>()*(v[0]+v[1]+v[2]); }
template<> finline Real GKDOP_Dot<3,14,4>( const Vec3& v ) { return mal::RcpSqrtOfThree<Real>()*(v[0]+v[1]-v[2]); }
template<> finline Real GKDOP_Dot<3,14,5>( const Vec3& v ) { return mal::RcpSqrtOfThree<Real>()*(v[0]-v[1]+v[2]); }
template<> finline Real GKDOP_Dot<3,14,6>( const Vec3& v ) { return mal::RcpSqrtOfThree<Real>()*(v[0]-v[1]-v[2]); }

//---- GKDOP static axis
template<unsigned D, unsigned K, unsigned A> finline mal::GVec<Real,D> GKDOP_Axis() { return mal::GVec<Real,D>::Zero(); }
//DOP2_K4 (X,Y)
template<> finline Vec2 GKDOP_Axis<2,4,0>() { return Vec2(1,0); }
template<> finline Vec2 GKDOP_Axis<2,4,1>() { return Vec2(0,1); }
//DOP2_K8 (X,Y,XY,Xy)
template<> finline Vec2 GKDOP_Axis<2,8,0>() { return Vec2(1,0); }
template<> finline Vec2 GKDOP_Axis<2,8,1>() { return Vec2(0,1); }
template<> finline Vec2 GKDOP_Axis<2,8,2>() { return Vec2(mal::RcpSqrtOfTwo<Real>(),mal::RcpSqrtOfTwo<Real>()); }
template<> finline Vec2 GKDOP_Axis<2,8,3>() { return Vec2(mal::RcpSqrtOfTwo<Real>(),-mal::RcpSqrtOfTwo<Real>()); }
//DOP3_K6 (X,Y,Z)
template<> finline Vec3 GKDOP_Axis<3,6,0>() { return Vec3(1,0,0); }
template<> finline Vec3 GKDOP_Axis<3,6,1>() { return Vec3(0,1,0); }
template<> finline Vec3 GKDOP_Axis<3,6,2>() { return Vec3(0,0,1); }
//DOP3_K14 (X,Y,Z,XYZ,XYz,XyZ,Xyz)
template<> finline Vec3 GKDOP_Axis<3,14,0>() { return Vec3(1,0,0); }
template<> finline Vec3 GKDOP_Axis<3,14,1>() { return Vec3(0,1,0); }
template<> finline Vec3 GKDOP_Axis<3,14,2>() { return Vec3(0,0,1); }
template<> finline Vec3 GKDOP_Axis<3,14,3>() { return Vec3(mal::RcpSqrtOfThree<Real>(), mal::RcpSqrtOfThree<Real>(), mal::RcpSqrtOfThree<Real>()); }
template<> finline Vec3 GKDOP_Axis<3,14,4>() { return Vec3(mal::RcpSqrtOfThree<Real>(), mal::RcpSqrtOfThree<Real>(),-mal::RcpSqrtOfThree<Real>()); }
template<> finline Vec3 GKDOP_Axis<3,14,5>() { return Vec3(mal::RcpSqrtOfThree<Real>(),-mal::RcpSqrtOfThree<Real>(), mal::RcpSqrtOfThree<Real>()); }
template<> finline Vec3 GKDOP_Axis<3,14,6>() { return Vec3(mal::RcpSqrtOfThree<Real>(),-mal::RcpSqrtOfThree<Real>(),-mal::RcpSqrtOfThree<Real>()); }

//---- GKDOP Dot product with variable axis
template<unsigned D, unsigned K> Real GKDOP_Dot( const mal::GVec<Real,D>& v, int i ) { return mal::Infinity<Real>(); }
//DOP2_K4
template<> finline Real GKDOP_Dot<2,4>( const Vec2& v, int i )
{
    switch(i)
    {
    case 0: return GKDOP_Dot<2,4,0>(v); break;
    case 1: return GKDOP_Dot<2,4,1>(v); break;
    }
    return Real(0);
}
//DOP2_K8
template<> finline Real GKDOP_Dot<2,8>( const Vec2& v, int i )
{
    switch(i)
    {
    case 0: return GKDOP_Dot<2,8,0>(v); break;
    case 1: return GKDOP_Dot<2,8,1>(v); break;
    case 2: return GKDOP_Dot<2,8,2>(v); break;
    case 3: return GKDOP_Dot<2,8,3>(v); break;
    }
    return Real(0);
}
//DOP3_K6
template<> finline Real GKDOP_Dot<3,6>( const Vec3& v, int i )
{
    switch(i)
    {
    case 0: return GKDOP_Dot<3,6,0>(v); break;
    case 1: return GKDOP_Dot<3,6,1>(v); break;
    case 2: return GKDOP_Dot<3,6,2>(v); break;
    }
    return Real(0);
}
//DOP3_K14
template<> finline Real GKDOP_Dot<3,14>( const Vec3& v, int i )
{
    switch(i)
    {
    case 0: return GKDOP_Dot<3,14,0>(v); break;
    case 1: return GKDOP_Dot<3,14,1>(v); break;
    case 2: return GKDOP_Dot<3,14,2>(v); break;
    case 3: return GKDOP_Dot<3,14,3>(v); break;
    case 4: return GKDOP_Dot<3,14,4>(v); break;
    case 5: return GKDOP_Dot<3,14,5>(v); break;
    case 6: return GKDOP_Dot<3,14,6>(v); break;
    }
    return Real(0);
}

//---- GKDOP variable axis
template<unsigned D, unsigned K> finline mal::GVec<Real,D> GKDOP_Axis( int i ) { return mal::GVec<Real,D>::Zero(); }
// DOP2_K4
template<> finline Vec2 GKDOP_Axis<2,4>( int i )
{
    switch( i )
    {
    case 0: return GKDOP_Axis<2,4,0>(); break;
    case 1: return GKDOP_Axis<2,4,1>(); break;
    }
    return Vec2(0);
}
// DOP2_K8
template<> finline Vec2 GKDOP_Axis<2,8>( int i )
{
    switch( i )
    {
    case 0: return GKDOP_Axis<2,8,0>(); break;
    case 1: return GKDOP_Axis<2,8,1>(); break;
    case 2: return GKDOP_Axis<2,8,2>(); break;
    case 3: return GKDOP_Axis<2,8,3>(); break;
    }
    return Vec2(0);
}
// DOP3_K6
template<> finline Vec3 GKDOP_Axis<3,6>( int i )
{
    switch( i )
    {
    case 0: return GKDOP_Axis<3,6,0>(); break;
    case 1: return GKDOP_Axis<3,6,1>(); break;
    case 2: return GKDOP_Axis<3,6,2>(); break;
    }
    return Vec3(0);
}
// DOP3_K14
template<> finline Vec3 GKDOP_Axis<3,14>( int i )
{
    switch( i )
    {
    case 0: return GKDOP_Axis<3,14,0>(); break;
    case 1: return GKDOP_Axis<3,14,1>(); break;
    case 2: return GKDOP_Axis<3,14,2>(); break;
    case 3: return GKDOP_Axis<3,14,3>(); break;
    case 4: return GKDOP_Axis<3,14,4>(); break;
    case 5: return GKDOP_Axis<3,14,5>(); break;
    case 6: return GKDOP_Axis<3,14,6>(); break;
    }
    return Vec3(0);
}

}} //namespace geo::bv

#endif // GEO_BV_GKDOP_AXIS_H
