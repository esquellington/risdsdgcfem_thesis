#include "TestOverlap.h"
#include "bv.h"

namespace geo {
namespace bv {

//! TEMPORAL: PLACEHOLDER IMPLEMENTATION
bool TestOverlap( const IBoundingVolume* p_bv1, const IBoundingVolume* p_bv2 )
{
    EBVType t1( p_bv1->GetType() );
    EBVType t2( p_bv2->GetType() );

    if( eBV_Void == t1 || eBV_Void == t2 ) return false; //\opt t1|t2==0
    else if( eBV_Infinite == t1 || eBV_Infinite == t2 ) return true; //\opt t1&t2==0xFF
    else
    {
        GEO_ASSERT( t1 == t2 );
        switch( t1 )
        {
        case eBV_Sphere2: return TestOverlap( p_bv1->As<Sphere2>(), p_bv2->As<Sphere2>() ); break;
        case eBV_Sphere3: return TestOverlap( p_bv1->As<Sphere3>(), p_bv2->As<Sphere3>() ); break;
        case eBV_AABB2: return TestOverlap( p_bv1->As<AABB2>(), p_bv2->As<AABB2>() ); break;
        case eBV_AABB3: return TestOverlap( p_bv1->As<AABB3>(), p_bv2->As<AABB3>() ); break;
        case eBV_LSS2: return TestOverlap( p_bv1->As<LSS2>(), p_bv2->As<LSS2>() ); break;
        case eBV_LSS3: return TestOverlap( p_bv1->As<LSS3>(), p_bv2->As<LSS3>() ); break;
        case eBV_DOP2_K4: return TestOverlap( p_bv1->As<DOP2_K4>(), p_bv2->As<DOP2_K4>() ); break;
        case eBV_DOP2_K8: return TestOverlap( p_bv1->As<DOP2_K8>(), p_bv2->As<DOP2_K8>() ); break;
        default:
            GEO_LOG_WARNING( "geo::bv::TestOverlap: BV type pair (%d,%d) not yet implemented", t1, t2 );
            return true; //BV default response is TRUE, so that an unimplemented BV does not discard a valid hit
            break;
        }
    }
}

//! TEMPORAL: PLACEHOLDER IMPLEMENTATION
bool TestOverlap( const BoundingVolume2* p_bv1, const BoundingVolume2* p_bv2 )
{
    EBVType t1( p_bv1->GetType() );
    EBVType t2( p_bv2->GetType() );

    if( eBV_Void2 == t1 || eBV_Void2 == t2 ) return false; //\opt t1|t2==0
    else if( eBV_Infinite2 == t1 || eBV_Infinite2 == t2 ) return true; //\opt t1&t2==0xFF
    else
    {
        GEO_ASSERT( t1 == t2 );
        switch( t1 )
        {
        case eBV_Sphere2: return TestOverlap( p_bv1->As<Sphere2>(), p_bv2->As<Sphere2>() ); break;
        case eBV_AABB2: return TestOverlap( p_bv1->As<AABB2>(), p_bv2->As<AABB2>() ); break;
        case eBV_LSS2: return TestOverlap( p_bv1->As<LSS2>(), p_bv2->As<LSS2>() ); break;
        case eBV_DOP2_K4: return TestOverlap( p_bv1->As<DOP2_K4>(), p_bv2->As<DOP2_K4>() ); break;
        case eBV_DOP2_K8: return TestOverlap( p_bv1->As<DOP2_K8>(), p_bv2->As<DOP2_K8>() ); break;
        default:
            GEO_LOG_WARNING( "geo::bv::TestOverlap: BV type pair (%d,%d) not yet implemented", t1, t2 );
            return true; //BV default response is TRUE, so that an unimplemented BV does not discard a valid hit
            break;
        }
    }
}

//! TEMPORAL: PLACEHOLDER IMPLEMENTATION
bool TestOverlap( const BoundingVolume3* p_bv1, const BoundingVolume3* p_bv2 )
{
    EBVType t1( p_bv1->GetType() );
    EBVType t2( p_bv2->GetType() );

    if( eBV_Void3 == t1 || eBV_Void3 == t2 ) return false; //\opt t1|t2==0
    else if( eBV_Infinite3 == t1 || eBV_Infinite3 == t2 ) return true; //\opt t1&t2==0xFF
    else
    {
        GEO_ASSERT( t1 == t2 );
        switch( t1 )
        {
        case eBV_Sphere3: return TestOverlap( p_bv1->As<Sphere3>(), p_bv2->As<Sphere3>() ); break;
        case eBV_AABB3: return TestOverlap( p_bv1->As<AABB3>(), p_bv2->As<AABB3>() ); break;
        case eBV_LSS3: return TestOverlap( p_bv1->As<LSS3>(), p_bv2->As<LSS3>() ); break;
        default:
            GEO_LOG_WARNING( "geo::bv::TestOverlap: BV type pair (%d,%d) not yet implemented", t1, t2 );
            return true; //BV default response is TRUE, so that an unimplemented BV does not discard a valid hit
            break;
        }
    }
}

}} //namespace geo::bv
