#ifndef GEO_CONFIG_H
#define GEO_CONFIG_H

#include <assert.h>
#include <pla_types.h>
#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GMat.h>
#include <Mal/GQuatUtils.h>
#include <Mal/GMatUtils.h>
#include <Mal/GTransform.h>
#include <Mal/GInterval.h>

namespace geo
{

#ifdef GEO_REAL_IS_DOUBLE
  typedef double Real;
#else
  typedef float Real;
#endif

// Standard Real types
typedef mal::GVec<Real,3> Vec3;
typedef mal::GVec<Real,3> Point3;
typedef mal::GQuat<Real> Quat;
typedef mal::GMat<Real,3,3> Mat3x3;
typedef mal::GMat<Real,2,3> Mat2x3;
typedef mal::GMat<Real,3,2> Mat3x2;
typedef mal::GMat<Real,4,4> Mat4x4;
typedef mal::GMat<Real,3,4> Mat3x4;
typedef mal::GMat<Real,4,3> Mat4x3;
typedef mal::GTransform<Real,3> Transform3;

typedef mal::GVec<Real,2> Vec2;
typedef mal::GVec<Real,2> Point2;
typedef mal::GMat<Real,2,2> Mat2x2;
typedef mal::GTransform<Real,2> Transform2;

typedef mal::GInterval<Real> Interval;

typedef mal::GVec<Real,4> Vec4;

// Global Tolerances
/* \todo NOOO they must be in geo::np::Context OR, if not limited to np, in geo::Context!!
static const Real gs_EpsilonLengthSq( 0.000001f );
static const Real gs_EpsilonLength( 0.0001f );
static const Real gs_EpsilonDir( 0.0001f );
*/
static const Real gs_EpsilonVertexLambda( 0.05f ); //[0..0.5], lambda fraction associated to Vertex instead of Edge.

} //namespace geo

#if defined(WIN32)
#  define finline _forceinline
#else
#  define finline inline
#endif

// Some defines...
// Use C++11 static_assert
#define GEO_STATIC_ASSERT(x) static_assert( (x), "x" )

// Some useful macros
#ifdef PROFILE_FINAL
#  define GEO_ASSERT( x )
#  define GEO_LOG( x, ... ) { printf("<GEO LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define GEO_LOG_WARNING( x, ... )
#  define GEO_LOG_WARNING_IF( c, x, ... )
#  define GEO_LOG_ERROR( x, ... )
#  define GEO_LOG_ERROR_IF( c, x, ... )
#  define GEO_LOG_ASSERT( c, x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define GEO_ASSERT( x ) assert(x)
#  define GEO_LOG( x, ... ) { printf("<GEO LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define GEO_LOG_WARNING( x, ... ) { printf("<GEO WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define GEO_LOG_WARNING_IF( c, x, ... ) { if(c) { printf("<GEO WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); } }
#  define GEO_LOG_ERROR( x, ... ) { printf("<GEO ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define GEO_LOG_ERROR_IF( c, x, ... ) { if(c) { printf("<GEO ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); } }
#  define GEO_LOG_ASSERT( c, x, ... )  { if(!(c)) { printf("<GEO ASSERT FAILED> "); printf( x, ##__VA_ARGS__ ); printf("\n"); assert(c); } }
#endif

// Paranoid checks can be enabled/disabled independently from regular ones
#define __GEO_ENABLE_PARANOID
#ifdef __GEO_ENABLE_PARANOID
#  include <assert.h>
#  include <stdio.h>
#  define GEO_PARANOID_ASSERT( x ) assert(x)
#else
#  define GEO_PARANOID_ASSERT( x )
#endif

/*\todo Geometric feature identification stuff... maybe should be in a
  FeatureTypes.h or in geo.h or confined to
  geo/shape/MeshTypes.h... but it's used in some NP types
  (ContactData), so it CANNOT be in shapes...  we'll see
*/
namespace geo
{

// Feature
typedef uint16 feature_count_type; //#V,#HE,#P
typedef feature_count_type feature_index_type; //VID,HEID,PID
static const feature_index_type cInvalidFeatureIndex = ~feature_index_type(0);

enum EFeatureType
{
    eFT_Invalid     = 0,
    eFT_Vertex      = 1,
    eFT_Segment     = 2,
    eFT_Triangle    = 3,
    eFT_Tetrahedron = 4
};

struct feature_id
{
    finline feature_id() : m_Type(eFT_Invalid), m_Index(cInvalidFeatureIndex) {}
    finline feature_id( EFeatureType ft, feature_index_type fid ) : m_Type(ft), m_Index(fid) {}
    finline void Set( EFeatureType ft, feature_index_type fid ) { m_Type = (uint16)ft; m_Index = fid; }
    finline bool IsValid() const { return eFT_Invalid != m_Type; }
    finline bool IsVertex() const { return eFT_Vertex == m_Type; }
    finline bool IsSegment() const { return eFT_Segment == m_Type; }
    finline bool IsTriangle() const { return eFT_Triangle == m_Type; }
    finline bool IsTetrahedron() const { return eFT_Tetrahedron == m_Type; }
    finline feature_index_type AsVertex() const { GEO_ASSERT(IsVertex()); return m_Index; }
    finline feature_index_type AsSegment() const { GEO_ASSERT(IsSegment()); return m_Index; }
    finline feature_index_type AsTriangle() const { GEO_ASSERT(IsTriangle()); return m_Index; }
    finline feature_index_type AsTetrahedron() const { GEO_ASSERT(IsTetrahedron()); return m_Index; }
    //\todo Consider redistributing bits as: uint32 m_Type : 3; uint32 m_Index : 29;
    uint16 m_Type; //EFeatureType
    feature_index_type m_Index;
};

struct PointOnFeature
{
    feature_id m_FeatureId;
    Vec4 m_BarycentricCoords; //b0,b1,b2,b3
    PointOnFeature() : m_FeatureId(), m_BarycentricCoords(0) {}
    PointOnFeature( feature_id fid, const Vec4 &bc ) : m_FeatureId(fid), m_BarycentricCoords(bc) {}
};

enum EEmbeddingMethod { eEM_Barycentric = 0,
                        eEM_MLS,
                        eEM_BNDF,
                        cNumEM };

#define __GEO_ENABLE_PARAMS
#define __GEO_ENABLE_STATS

// OPTIONAL SSE
#if __SSE4_1__
#  define __GEO_ENABLE_NP_SIMD
#endif

} //namespace geo

#endif //GEO_CONFIG_H
