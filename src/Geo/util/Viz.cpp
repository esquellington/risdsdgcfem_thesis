#include "Viz.h"

//-- base entity types
#include <Geo/IObject.h>
#include <Geo/bv/BoundingVolume.h>
//-- shape types
#include <Geo/shape/shape.h>
//-- BV types
#include <Geo/bv/bv.h>
//-- BP types
#include <Geo/bp/QuadraticTimeBP.h>
#include <Geo/bp/SweepAndPruneBP.h>
//-- NP types
#include <Geo/np/ContactData.h>
//-- Misc
#include <Mal/GConversion.h>

#define __ENABLE_DCR3_DRAW_ERRORS
extern float g_TTS_SliverThreshold; //TEMP: Ugly hack...

namespace geo {

//----------------------------------------------------------------
// GVizBV
//\todo This function template is used in VizBVH only, but could be useful in VizBoundingVolume() too...
//----------------------------------------------------------------
template <typename BVT> void GVizBV( util::VizStream& vs, const BVT& bv, float height, const Vec4f& color ) {}
//-- 2D
template<> void GVizBV( util::VizStream& vs, const bv::Sphere2& bv, float height, const Vec4f& color )
{
    vs.BeginComplex(1,util::eType_VizDisk);
    {
        /* Geometry */
        vs.Write("pos", Vec3f( bv.GetPos()[0],bv.GetPos()[1],height ) );
        vs.Write("rot", Quatf::Identity() );
        vs.Write("radius", (float)bv.GetRadius() );
        /* Style */
        vs.BeginComplex("style",util::eType_VizStyle);
        {
            vs.Write("color",Vec4f(color));
            vs.Write("flags",Flags32(util::eVizStyle_Wire));
        }
        vs.EndComplex();
    }
    vs.EndComplex();
}
template<> void GVizBV( util::VizStream& vs, const bv::AABB2& bv, float height, const Vec4f& color )
{
    vs.BeginComplex(1,util::eType_VizRectangle);
    {
        /* Geometry */
        vs.Write("pos", Vec3f( bv.GetPos()[0],bv.GetPos()[1],height ) );
        vs.Write("rot", Quatf::Identity() );
        vs.Write("half_sizes", Vec2f(bv.GetHalfSizes()) );
        /* Style */
        vs.BeginComplex("style",util::eType_VizStyle);
        {
            vs.Write("color",Vec4f(color));
            vs.Write("flags",Flags32(util::eVizStyle_Wire));
        }
        vs.EndComplex();
    }
    vs.EndComplex();
}
template<> void GVizBV( util::VizStream& vs, const bv::DOP2_K8& bv, float height, const Vec4f& color )
{
    // Get KDOP segments
    std::pair<Vec2,Vec2> vec_segments[8];
    unsigned int num_segments = bv::DOP2_K8_Segments( bv, vec_segments );
    // Viz clipped segments, regardless of actual order
    for( unsigned int it_segment=0; it_segment < num_segments; it_segment++ )
        VIZ_SEGMENT3( vs,
                      Vec3f( vec_segments[it_segment].first[0], vec_segments[it_segment].first[1], height ),
                      Vec3f( vec_segments[it_segment].second[0], vec_segments[it_segment].second[1], height ),
                      1.0f, color );
}

//-- 3D
template<> void GVizBV( util::VizStream& vs, const bv::Sphere3& bv, float height, const Vec4f& color )
{
    vs.BeginComplex(1,util::eType_VizSphere);
    {
        /* Geometry */
        vs.Write("pos", Vec3f( bv.GetPos() ) );
        vs.Write("rot", Quatf::Identity() );
        vs.Write("radius", (float)bv.GetRadius() );
        /* Style */
        vs.BeginComplex("style",util::eType_VizStyle);
        {
            vs.Write("color",Vec4f(color));
            vs.Write("flags",Flags32(util::eVizStyle_Wire));
        }
        vs.EndComplex();
    }
    vs.EndComplex();
}
template<> void GVizBV( util::VizStream& vs, const bv::AABB3& bv, float height, const Vec4f& color )
{
    vs.BeginComplex(1,util::eType_VizBox);
    {
        /* Geometry */
        vs.Write("pos", Vec3f(bv.GetPos()) );
        vs.Write("rot", Quatf::Identity() );
        vs.Write("half_sizes", Vec3f(bv.GetHalfSizes()) );
        /* Style */
        vs.BeginComplex("style",util::eType_VizStyle);
        {
            vs.Write("color",Vec4f(color));
            vs.Write("flags",Flags32(util::eVizStyle_Wire));
        }
        vs.EndComplex();
    }
    vs.EndComplex();
}
template<> void GVizBV( util::VizStream& vs, const bv::DOP3_K6& bv, float height, const Vec4f& color )
{
    vs.BeginComplex(1,util::eType_VizBox);
    {
        /* Geometry */
        vs.Write("pos", Vec3f(bv.GetPos()) );
        vs.Write("rot", Quatf::Identity() );
        vs.Write("half_sizes", 0.5f*Vec3f(bv.GetInterval<0>().Length(),
                                          bv.GetInterval<1>().Length(),
                                          bv.GetInterval<2>().Length()) );
        /* Style */
        vs.BeginComplex("style",util::eType_VizStyle);
        {
            vs.Write("color",Vec4f(color));
            vs.Write("flags",Flags32(util::eVizStyle_Wire));
        }
        vs.EndComplex();
    }
    vs.EndComplex();
}
template<> void GVizBV( util::VizStream& vs, const bv::DOP3_K14& bv, float height, const Vec4f& color )
{
    //TEMP: SAME AS DOP3_K6 by now, axis 3..6 not shown
    vs.BeginComplex(1,util::eType_VizBox);
    {
        /* Geometry */
        vs.Write("pos", Vec3f(bv.GetPos()) );
        vs.Write("rot", Quatf::Identity() );
        vs.Write("half_sizes", 0.5f*Vec3f(bv.GetInterval<0>().Length(),
                                          bv.GetInterval<1>().Length(),
                                          bv.GetInterval<2>().Length()) );
        /* Style */
        vs.BeginComplex("style",util::eType_VizStyle);
        {
            vs.Write("color",Vec4f(color));
            vs.Write("flags",Flags32(util::eVizStyle_Wire));
        }
        vs.EndComplex();
    }
    vs.EndComplex();
}

//----------------------------------------------------------------
// VizBVH \note MUST appear before uses due to specializations,
// otherwise undefined references
// ----------------------------------------------------------------
template <typename BoundingVolumeT, typename EntryIndexT>
void GVizBVH_Internal( const bv::GBoundingVolumeHierarchy_ST_DG<BoundingVolumeT,EntryIndexT>& bvh, util::VizStream &vs, Flags32 ddf )
{
    typedef bv::GBoundingVolumeHierarchy_ST_DG<BoundingVolumeT,EntryIndexT> bvh_type;
    if( ddf.Test( eBVHDDF_Geometry ) )
    {
        if( ddf.Test( eBVHDDF_Levels ) )
            bvh.Map( [&vs]( const bvh_type& bvh, const typename bvh_type::Node& node, uint32 level )
                     {
                         float height01( float(level) / bvh.GetHeight() );
                         GVizBV( vs, node.m_Geometry.m_BV, height01, Vec4f(1,1,1-height01,0.25f) );
                         return true;
                     }
                );
        else
            for( const auto& node : bvh.m_vecNode )
                GVizBV( vs, node.m_Geometry.m_BV, 0, Vec4f(1,1,1,0.25f) );
    }
}
//----------------------------------------------------------------
// Explicit instantiations... there seems to be no way of instancing a
// tplt declared, but not defined, in Viz.h without using an
// "internal" function and explicitly instantiating them fully...
// ----------------------------------------------------------------
//2D
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::Sphere2,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::Sphere2,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB2,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB2,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP2_K8,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP2_K8,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
//3D
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::Sphere3,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::Sphere3,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB3,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::AABB3,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP3_K6,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP3_K6,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP3_K14,uint16>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }
template<> void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<bv::DOP3_K14,uint32>& bvh, util::VizStream &vs, Flags32 ddf ) { GVizBVH_Internal(bvh,vs,ddf); }

//----------------------------------------------------------------
//----------------------------------------------------------------

void VizObject( const IObject *p_go, util::VizStream &vs, Flags32 ddf )
{
    char name[256];

    const IShape *pShape( p_go->GetShapeInterface() );

    if( 0 == pShape )
    {
        GEO_LOG_WARNING( "VizObject: Object has no IShape interface... it's probably an uninitialized GObjectSS" );
        return;
    }

    if( ddf.Test( eODDF_Axis ) )
    {
        if( p_go->GetDimension() == 2 )
        {
            VIZ_TRANSFORM2( vs, static_cast<const IObject2*>(p_go)->GetTransform(), Real(1) );
        }
        else //p_go->GetDimension() == 3 )
        {
            VIZ_TRANSFORM3( vs, static_cast<const IObject3*>(p_go)->GetTransform(), Real(1) );
        }
    }

    if( ddf.Test( eODDF_BV ) )
    {
        //\todo Draw BV... which one?
    }

    if( ddf.Test( eODDF_Shape ) )
    {
        switch( pShape->GetType() )
        {
        case eShape_Plane2:
            {
                const PlaneShape2 *pPlane2( static_cast<const PlaneShape2*>(pShape) );
                Vec2 n( pPlane2->GetNormal() );
                Vec2 p(0,0);
                if( n[0] != Real(0) )
                    p[0] = -pPlane2->GetCoeffD() / n[0];
                else if( n[1] != Real(0) )
                    p[1] = -pPlane2->GetCoeffD() / n[1];

                Real size(100);
                Vec2 p1( p + size*Vec2(-n[1],n[0]) );
                Vec2 p2( p - size*Vec2(-n[1],n[0]) );

                // Plane, Normal
                VIZ_SEGMENT2( vs, p1, p2, 3.0f, Vec4f(1,0,0,0.5f) );
                VIZ_VECTOR2( vs, p, n, 1.0f, Vec4f(0,1,0,1) );
                // double-sided normal if !IsHalfSpace()
                if( !pPlane2->IsHalfSpace() ) VIZ_VECTOR2( vs, p, -n, 2.0f, Vec4f(0,0.25,0,1) );
            }
            break;
        case eShape_Plane3:
            {
                const PlaneShape3 *pPlane3( static_cast<const PlaneShape3*>(pShape) );
                Vec3 n( pPlane3->GetNormal() );
                Vec3 p(0,0,0);
                Vec3 u(0,0,0);
                if( n[0] != Real(0) )
                {
                    p[0] = -pPlane3->GetCoeffD() / n[0];
                    u = mal::Normalized( mal::Cross(n,Vec3(0,1,0) ) );
                }
                else if( n[1] != Real(0) )
                {
                    p[1] = -pPlane3->GetCoeffD() / n[1];
                    u = mal::Normalized( mal::Cross(n,Vec3(0,0,1) ) );
                }
                else if( n[2] != Real(0) )
                {
                    p[2] = -pPlane3->GetCoeffD() / n[2];
                    u = mal::Normalized( mal::Cross(n,Vec3(1,0,0) ) );
                }

                Vec3 v( mal::Cross(n,u) );

                Real size(10);
                Vec3 p1( p - size*u - size*v );
                Vec3 p2( p + size*u - size*v );
                Vec3 p3( p + size*u + size*v );
                Vec3 p4( p - size*u + size*v );

                // Plane, Normal
                VIZ_QUAD3( vs, p1, p2, p3, p4, 4.0f, Vec4f(0.5,0,0,1), util::eVizStyle_Wire );
                VIZ_VECTOR3( vs, p, n, 4.0f, Vec4f(0,1,0,1) );
                // double-sided normal if !IsHalfSpace()
                if( !pPlane3->IsHalfSpace() ) VIZ_VECTOR3( vs, p, -n, 4.0f, Vec4f(0,0.25,0,1) );
            }
            break;
        case eShape_Sphere2:
            {
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    VIZ_DISK2( vs, tr2.m_Pos, mal::GQuat_From(tr2.m_Rot),
                               static_cast<const SphereShape2*>(pShape)->GetRadius(),
                               Vec4f(1,0,0,1), util::eVizStyle_Wire );
                }
                if( ddf.Test( eSDDF_Interior ) )
                {
                    VIZ_DISK2( vs, tr2.m_Pos, mal::GQuat_From(tr2.m_Rot),
                               static_cast<const SphereShape2*>(pShape)->GetRadius(),
                               Vec4f(1,0,0,0.5f), util::eVizStyle_Solid );
                }
            }
            break;
        case eShape_Sphere3:
            {
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    Transform3f tr3( static_cast<const IObject3*>(p_go)->GetTransform() );
                    VIZ_SPHERE3( vs, tr3.m_Pos, mal::GQuat_From(tr3.m_Rot),
                                 static_cast<const SphereShape3*>(pShape)->GetRadius(),
                                 Vec4f(1,0,0,1), util::eVizStyle_Wire );
                }
                if( ddf.Test( eSDDF_Interior ) )
                {
                    Transform3f tr3( static_cast<const IObject3*>(p_go)->GetTransform() );
                    VIZ_SPHERE3( vs, tr3.m_Pos, mal::GQuat_From(tr3.m_Rot),
                                 static_cast<const SphereShape3*>(pShape)->GetRadius(),
                                 Vec4f(1,0,0,0.5f), util::eVizStyle_Solid );
                }
            }
            break;
        case eShape_Capsule2:
            {
                const CapsuleShape2 *pCapsule( static_cast<const CapsuleShape2*>(pShape) );
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    VIZ_DISK2( vs,
                               tr2.m_Pos - mal::Column(1,tr2.m_Rot) * pCapsule->GetHalfHeight(),
                               mal::GQuat_From(tr2.m_Rot),
                               pCapsule->GetRadius(),
                               Vec4f(1,0,1,1), util::eVizStyle_Wire );
                    VIZ_DISK2( vs,
                               tr2.m_Pos + mal::Column(1,tr2.m_Rot) * pCapsule->GetHalfHeight(),
                               mal::GQuat_From(tr2.m_Rot),
                               pCapsule->GetRadius(),
                               Vec4f(1,0,1,1), util::eVizStyle_Wire );
                    VIZ_RECT2( vs, tr2.m_Pos, mal::GQuat_From(tr2.m_Rot),
                               Vec2( pCapsule->GetRadius(), pCapsule->GetHalfHeight() ),
                               Vec4f(1,0,0,1), util::eVizStyle_Wire );
                }
                if( ddf.Test( eSDDF_Interior ) )
                {
                    VIZ_DISK2( vs,
                               tr2.m_Pos - mal::Column(1,tr2.m_Rot) * pCapsule->GetHalfHeight(),
                               mal::GQuat_From(tr2.m_Rot),
                               pCapsule->GetRadius(),
                               Vec4f(1,0,1,0.5f), util::eVizStyle_Solid );
                    VIZ_DISK2( vs,
                               tr2.m_Pos + mal::Column(1,tr2.m_Rot) * pCapsule->GetHalfHeight(),
                               mal::GQuat_From(tr2.m_Rot),
                               pCapsule->GetRadius(),
                               Vec4f(1,0,1,0.5f), util::eVizStyle_Solid );
                    VIZ_RECT2( vs, tr2.m_Pos, mal::GQuat_From(tr2.m_Rot),
                               Vec2( pCapsule->GetRadius(), pCapsule->GetHalfHeight() ),
                               Vec4f(1,0,0,0.5f), util::eVizStyle_Solid );
                }
            }
            break;
        case eShape_Box2:
            {
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    VIZ_RECT2( vs, tr2.m_Pos, mal::GQuat_From(tr2.m_Rot),
                               static_cast<const BoxShape2*>(pShape)->GetHalfSizes(),
                               Vec4f(1,0,0,1), util::eVizStyle_Wire );
                }
                if( ddf.Test( eSDDF_Interior ) )
                {
                    VIZ_RECT2( vs, tr2.m_Pos, mal::GQuat_From(tr2.m_Rot),
                               static_cast<const BoxShape2*>(pShape)->GetHalfSizes(),
                               Vec4f(1,0,0,0.5f), util::eVizStyle_Solid );
                }
            }
            break;
        case eShape_Box3:
            {
                Transform3f tr3( static_cast<const IObject3*>(p_go)->GetTransform() );
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    VIZ_BOX3( vs, tr3.m_Pos, mal::GQuat_From(tr3.m_Rot),
                              static_cast<const BoxShape3*>(pShape)->GetHalfSizes(),
                              Vec4f(1,0,0,1), util::eVizStyle_Wire );
                }
                if( ddf.Test( eSDDF_Interior ) )
                {
                    VIZ_BOX3( vs, tr3.m_Pos, mal::GQuat_From(tr3.m_Rot),
                              static_cast<const BoxShape3*>(pShape)->GetHalfSizes(),
                              Vec4f(1,0,0,0.5f), util::eVizStyle_Solid );
                }
            }
            break;
        case eShape_SphereSet2:
            {
                const SphereSetShape2* pSSS( static_cast<const SphereSetShape2*>(pShape) );
                const GObjectSDOF<2,Vec2> *pSSSGO( static_cast<const GObjectSDOF<2,Vec2> *>(p_go) );
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                if( ddf.Test( eSDDF_Boundary ) )
                    for( unsigned int it_sphere=0; it_sphere < pSSS->GetCount(); it_sphere++ )
                    {
                        VIZ_DISK2_NO_ROT( vs, pSSSGO->GetSDOF(it_sphere), pSSS->GetRadius(), Vec4f(1,0,0,1), util::eVizStyle_Wire );
                    }
                if( ddf.Test( eSDDF_Interior ) )
                    for( unsigned int it_sphere=0; it_sphere < pSSS->GetCount(); it_sphere++ )
                    {
                        VIZ_DISK2_NO_ROT( vs, pSSSGO->GetSDOF(it_sphere), pSSS->GetRadius(), Vec4f(1,0,0,0.5f), util::eVizStyle_Solid );
                    }
            }
            break;
        case eShape_SphereSet3:
            //\todo!!
            GEO_ASSERT(false);
            break;
        case eShape_Polygonal2:
            {
                const PolygonalShape2* pPS2( static_cast<const PolygonalShape2*>(pShape) );
                const GObjectSDOF<2,Vec2> *pPSGO( static_cast<const GObjectSDOF<2,Vec2> *>(p_go) );
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                int num_edges( pPS2->IsClosed() ? pPS2->GetNumVertices() : pPS2->GetNumVertices()-1 );
                const Vec2 *actual_sdof( pPSGO->GetVecSDOF() );
                if( 0 == actual_sdof )
                {
                    GEO_LOG_WARNING("Object has no SDOF, using default");
                    actual_sdof = pPS2->GetVecDefaultSDOF();
                }

                // Vertices
                if( ddf.Test( eSDDF_Vertices ) )
                    for( unsigned int it_v=0; it_v < pPS2->GetNumVertices(); it_v++ )
                    {
                        Vec2 p( tr2 * pPS2->V_Pos( it_v, actual_sdof ) );
                        if( ddf.Test( eSDDF_FeatureId ) )
                        {
                            sprintf( name, "V%u", it_v );
                            VIZ_POINT2_NAMED( vs, name, p, 1.0f, Vec4f(1,1,0,0.5f) );
                        }
                        else
                            VIZ_POINT2( vs, p, 5.0f, Vec4f(1,0.5,0.5,1) );
                    }

                // Boundary | Edges
                if( ddf.Test( eSDDF_Boundary ) )
                    for( int it_edge=0; it_edge < num_edges; it_edge++ )
                    {
                        Vec2 p0( tr2 * pPS2->V_Pos( it_edge, actual_sdof ) );
                        Vec2 p1( tr2 * pPS2->V_Pos( (it_edge+1) % pPS2->GetNumVertices(), actual_sdof ) );
                        VIZ_SEGMENT2( vs, p0, p1, 3.0f, Vec4f(1,0,0,1) );
                    }

                // \todo Interior, Topology....

                // Normals
                if( ddf.Test( eSDDF_Normals ) )
                {
                    if( pPS2->HasVDL( eVDL_Normals ) )
                        for( unsigned it_vtx=0; it_vtx < pPS2->GetNumVertices(); it_vtx++ )
                        {
                            Vec2 p( tr2 * pPS2->V_Pos( it_vtx, actual_sdof ) );
                            Vec2 n( tr2.m_Rot * pPS2->GetNormal( it_vtx ) ); //\todo Normal is not recomputed on actual sdof?!
                            VIZ_SEGMENT2( vs, p, p+n, 1.0f, Vec4f(0,1,0,0.5f) );
                        }
                }
                // Tangents
                if( ddf.Test( eSDDF_Tangents ) )
                {
                   if( pPS2->HasVDL( eVDL_Tangents ) )
                       for( unsigned int it_vtx=0; it_vtx < pPS2->GetNumVertices(); it_vtx++ )
                       {
                           Vec2 p( tr2 * pPS2->V_Pos( it_vtx, actual_sdof ) );
                           Vec2 t( tr2.m_Rot * pPS2->GetTangent( it_vtx ) ); //\todo Tangent is not recomputed on actual sdof?!
                           VIZ_SEGMENT2( vs, p, p+t, 1.0f, Vec4f(0,0,1,0.5f) );
                       }
                }

                /* Lambdas
                   if( pPS2->HasVDL( eVDL_Lambdas ) )
                   for( int it_vtx=0; it_vtx < pPS2->GetNumVertices(); it_vtx++ )
                   {
                   Vec2 p( tr2 * pPS2->V_Pos( it_vtx, actual_sdof ) );
                   VIZ_DISK2_NO_ROT( vs, p, pPS2->GetLambda(it_vtx), Vec4f(1,1,1,0.5f), util::eVizStyle_Wire );  //\todo Lambdas is not recomputed on actual sdof?!
                   }
                */
                /*
                //\todo eVDL_Flags
                //Barycenter
                Vec2 p( tr2 * pPS2->GetBarycenter( actual_sdof ) );
                VIZ_DISK2_NO_ROT( vs, p, 0.1, Vec4f(1,1,0,1), util::eVizStyle_Solid );
                */

#ifdef __GEO_ENABLE_IMPLICIT_RBF
                // Implicit
                const IImplicit2 *pImplicit( pPS2->GetImplicit2() );
                bv::IBoundingVolume bv(bv::eBV_AABB2);
                pPS2->ComputeBV( bv, tr2, 0 );
                bv::AABB2 &aabb2( bv.As<bv::AABB2>() );
                aabb2.SetPosHalfSizes( aabb2.GetPos(), Real(2)*aabb2.GetHalfSizes() ); //x2
                Vec2 cStep( aabb2.GetHalfSizes() / 40 );
                Transform2f inv_tr2( tr2.Inverse() );
                for( Real x=-aabb2.GetHalfSizes().x(); x<aabb2.GetHalfSizes().x(); x+=cStep.x() )
                    for( Real y=-aabb2.GetHalfSizes().y(); y<aabb2.GetHalfSizes().y(); y+=cStep.y() )
                    {
                        Vec2 p( aabb2.GetPos() + Vec2(x,y) );
                        Real f( pImplicit->F( inv_tr2 * p ) ); //Eval implicit in local coords
                        if( f < 0 )
                        {
                            //\todo Enable to viz Outside region
                            //VIZ_POINT2( vs, p, 1.0f, Vec4f(1,1,0,1) );

                            //VIZ_DISK2_NO_ROT( vs, p, 0.05, Vec4f(1,1,0,1), util::eVizStyle_Solid );
                            //VIZ_SEGMENT3( vs, Vec3f(p.x(),p.y(),0), Vec3f(p.x(),p.y(),f), 1.0f, Vec4f(1,1,0,1) );
                        }
                        else
                        {
                            VIZ_POINT2( vs, p, 3.0f, Vec4f(0,1,1,1) );
                            //VIZ_DISK2_NO_ROT( vs, p, 0.1, Vec4f(0,1,1,1), util::eVizStyle_Solid );
                            //VIZ_SEGMENT3( vs, Vec3f(p.x(),p.y(),0), Vec3f(p.x(),p.y(),f), 1.0f, Vec4f(0,1,1,1) );
                        }
                    }
#endif
            }
            break;
        case eShape_Polygonal3:
            {
                GEO_ASSERT(false);
                /*
                const PolygonalShape3* pPS3( static_cast<const PolygonalShape3*>(pShape) );
                const GObjectSDOF<3,Vec3> *pPSGO( static_cast<const GObjectSDOF<3,Vec3> *>(p_go) );
                Transform3f tr3( static_cast<const IObject3*>(p_go)->GetTransform() );
                int num_edges( pPS3->IsClosed() ? pPS3->GetNumVertices() : pPS3->GetNumVertices()-1 );
                const Vec3 *actual_sdof( pPSGO->GetVecSDOF() );
                if( 0 == actual_sdof )
                {
                    GEO_LOG_WARNING("Object has no SDOF, using default");
                    actual_sdof = pPS3->GetVecDefaultSDOF();
                }
                for( int it_edge=0; it_edge < num_edges; it_edge++ )
                {
                    Vec3f p0( tr3 * pPS3->V_Pos( it_edge, actual_sdof ) );
                    Vec3f p1( tr3 * pPS3->V_Pos( (it_edge+1) % pPS3->GetNumVertices(), actual_sdof ) );
                    VIZ_SEGMENT3( vs, p0, p1, 1.0f, Vec4f(1,0,0,0.5f) );
                }
                */
            }
            break;
        case eShape_Path2:
            {
                const PathShape2* pPS2( static_cast<const PathShape2*>(pShape) );
                const GObjectSDOF<2,Vec2> *pPSGO( static_cast<const GObjectSDOF<2,Vec2> *>(p_go) );
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                const Vec2 *actual_sdof( pPSGO->GetVecSDOF() );
                if( 0 == actual_sdof )
                {
                    GEO_LOG_WARNING("Object has no SDOF, using default");
                    actual_sdof = pPS2->GetVecDefaultSDOF();
                }
                // Boundary | Edges
                for( unsigned int it_edge=0; it_edge < pPS2->GetNumE(); it_edge++ )
                {
                    Vec2 p0( tr2 * pPS2->V_Pos( pPS2->E_OriginVID(it_edge), actual_sdof ) );
                    Vec2 p1( tr2 * pPS2->V_Pos( pPS2->E_FinalVID(it_edge), actual_sdof ) );
                    if( ddf.Test( eSDDF_Vertices ) )
                    {
                        VIZ_POINT2( vs, p0, 5.0f, Vec4f(1,0.25,0.25,1) );
                        VIZ_POINT2( vs, p1, 5.0f, Vec4f(1,0.25,0.25,1) );
                    }
                    switch( pPS2->E_CurveType(it_edge) )
                    {
                    case PathShape2::edge_type::eCT_Line:
                        if( ddf.Test( eSDDF_Boundary ) )
                        {
                            VIZ_SEGMENT2( vs, p0, p1, 3.0f, Vec4f(1,0,0,1) );
                        }
                        break;
                    case PathShape2::edge_type::eCT_Bezier3:
                        {
                            // Bezier curve
                            Vec2 bezier3_a( tr2 * pPS2->E_Data(it_edge).m_ParamA ); //\todo If params are also SDOF, use actual_sdof here!!
                            Vec2 bezier3_b( tr2 * pPS2->E_Data(it_edge).m_ParamB );
                            if( ddf.Test( eSDDF_Boundary ) )
                            {
                                Vec2 b0( p0 );
                                const unsigned int cNumSamples(10);
                                for( unsigned int i=1; i<cNumSamples; i++ )
                                {
                                    Vec2 b1 = Eval_Bezier3( p0, p1, bezier3_a, bezier3_b, Real(i) / (cNumSamples-1) );
                                    VIZ_SEGMENT2( vs, b0, b1, 3.0f, Vec4f(1,0,0,1) );
                                    b0 = b1;
                                }
                            }
                            // Bezier points
                            if( ddf.Test( eSDDF_Control ) )
                            {
                                VIZ_SEGMENT2( vs, p0, bezier3_a, 1.0f, Vec4f(1,0,1,1) );
                                VIZ_SEGMENT2( vs, bezier3_a, bezier3_b, 1.0f, Vec4f(1,0,1,1) );
                                VIZ_SEGMENT2( vs, bezier3_b, p1, 1.0f, Vec4f(1,0,1,1) );
                            }
                        }
                        break;
                    default: break;
                    }
                }
                /* Orientation \todo There's better ways to do this...
                Vec2 p0( tr2 * pPS2->V_Pos( pPS2->E_OriginVID(0), actual_sdof ) );
                Vec2 p1( tr2 * pPS2->V_Pos( pPS2->E_FinalVID(0), actual_sdof ) );
                VIZ_DISK2_NO_ROT( vs, p0, 5.0f, Vec4f(1,1,1,0.5f), util::eVizStyle_Wire );
                VIZ_DISK2_NO_ROT( vs, p1, 1.0f, Vec4f(1,1,1,0.5f), util::eVizStyle_Wire );
                */
            }
            break;
        case eShape_MeshSolid2:
            {
                const MeshSolidShape2* pMSS2( static_cast<const MeshSolidShape2*>(pShape) );
                const GObjectSDOF<2,Vec2> *pMSSGO( static_cast<const GObjectSDOF<2,Vec2> *>(p_go) );
                Transform2f tr2( static_cast<const IObject2*>(p_go)->GetTransform() );
                const Vec2 *actual_sdof( pMSSGO->GetVecSDOF() );
                if( 0 == actual_sdof )
                {
                    GEO_LOG_WARNING("Object has no SDOF, using default");
                    actual_sdof = pMSS2->GetVecDefaultSDOF();
                }

                // Boundary | Edges
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    for( unsigned int it_bp=0; it_bp < pMSS2->GetNumBoundaryP(); it_bp++ )
                    {
                        // Barycenter
                        Vec2 barycenter( tr2 * pMSS2->BP_Barycenter( it_bp, actual_sdof ) );
                        // HE, expanded from barycenter to differentiate from symmetrics
                        const Real cContractionCoeff( 1.01f );
                        unsigned int it_he( pMSS2->BP_FirstHEID(it_bp) );
                        int num_neighbours(0);
                        Vec2 p0( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) );
                        if( ddf.Test( eSDDF_Vertices ) )
                        {
                            VIZ_POINT2( vs, p0, 5.0f, Vec4f(1,0.25,0.25,1) );
                        }
                        do
                        {
                            if( pMSS2->HE_Sym(it_he) != cInvalidFeatureIndex ) num_neighbours++;
                            it_he = pMSS2->HE_Next(it_he);
                            Vec2 p1( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) );
                            VIZ_SEGMENT2( vs,
                                          barycenter + cContractionCoeff*(p0-barycenter),
                                          barycenter + cContractionCoeff*(p1-barycenter),
                                          1.0f, Vec4f(1,0,0,1) );

                            if( ddf.Test( eSDDF_Vertices ) )
                            {
                                VIZ_POINT2( vs, p1, 5.0f, Vec4f(1,0.25,0.25,1) );
                            }

                            /*TEMP
                              sprintf( name, "HE%d,S%d", it_he, pMSS2->HE_Sym(it_he) );
                              VIZ_POINT2_NAMED( vs, name,
                              barycenter + 0.5f*0.8f*( cContractionCoeff*(p0-barycenter) + cContractionCoeff*(p1-barycenter) ),
                              2.0f, Vec4f(1,0,1,0.5f) );
                              TEMP*/

                            p0 = p1;
                        } while ( it_he != pMSS2->BP_FirstHEID(it_bp) );

                        /*TEMP: Num neighbours
                          sprintf( name, "P%d,#N%d", it_bp, num_neighbours );
                          VIZ_POINT2_NAMED( vs, name, barycenter, 5.0f, Vec4f(1,1,0,0.5f) );
                        */
                    }
                }

                // Interior (Polygons)
                if( ddf.Test( eSDDF_Interior ) )
                {
                    for( unsigned int it_p=0; it_p < pMSS2->GetNumP(); it_p++ )
                    {
                        // Barycenter
                        Vec2 barycenter( tr2 * pMSS2->P_Barycenter( it_p, actual_sdof ) );
                        //VIZ_DISK2_NO_ROT( vs, barycenter, 0.01, Vec4f(1,1,0,1), util::eVizStyle_Solid );

                        // HE, contracted towards barycenter to differentiate from symmetrics
                        const Real cContractionCoeff( 0.95f );
                        unsigned int it_he( pMSS2->P_FirstHEID(it_p) );
                        int num_neighbours(0);
                        Vec2 p0( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) );
                        do
                        {
                            /*TEMP
                            sprintf( name, "HE%d,S%d,L%d,R%d", it_he, pMSS2->HE_Sym(it_he), pMSS2->HE_LeftPID(it_he), pMSS2->HE_RightPID(it_he) );
                            //TEMP*/
                            if( pMSS2->HE_Sym(it_he) != cInvalidFeatureIndex ) num_neighbours++;
                            it_he = pMSS2->HE_Next(it_he);
                            Vec2 p1( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) );
                            VIZ_SEGMENT2( vs,
                                          barycenter + cContractionCoeff*(p0-barycenter),
                                          barycenter + cContractionCoeff*(p1-barycenter),
                                          1.0f, Vec4f(1,0,0,0.5) );
                            /*TEMP
                            VIZ_POINT2_NAMED( vs, name,
                                              barycenter + 0.5f*0.8f*( cContractionCoeff*(p0-barycenter) + cContractionCoeff*(p1-barycenter) ),
                                              2.0f, Vec4f(1,0,1,0.5f) );
                            //TEMP*/

                            p0 = p1;
                        } while ( it_he != pMSS2->P_FirstHEID(it_p) );

                        /*TEMP: PID
                          sprintf( name, "P%d", it_p );
                          VIZ_POINT2_NAMED( vs, name, barycenter, 1.0f, Vec4f(1,1,0,0.5f) );
                        */

                        /*TEMP: Num neighbours
                          sprintf( name, "P%d,#N%d", it_p, num_neighbours );
                          VIZ_POINT2_NAMED( vs, name, barycenter, 5.0f, Vec4f(1,1,0,0.5f) );
                        */
                    }

                    // Interior Vertices
                    for( unsigned int it_v=0; it_v < pMSS2->GetNumV(); it_v++ )
                    {
                        Vec2 p( tr2 * pMSS2->V_Pos( it_v, actual_sdof ) );
                        VIZ_POINT2( vs, p, 5.0f, Vec4f(1,0.25,0.25,0.5) );
                        /*TEMP: Named
                          sprintf( name, "V[%d]=%d", it_v, pMSS2->V_NumEdges(it_v) );
                          VIZ_POINT2_NAMED( vs, name, p, 5.0f, Vec4f(1,1,1,0.5f) );
                        */

                        /* iterator_around_vertex_ccw
                        //TEMP: debug it_ccw
                        int num_ccw(0);
                        MeshSolidShape2::iterator_around_vertex_ccw it_ccw( pMSS2->GetIteratorAroundVertexCCW(it_v) );
                        while( it_ccw.IsValid() )
                        {
                        ++num_ccw;
                        ++it_ccw;
                        }
                        sprintf( name, "V%d,HE%d,#P%d,O%d", it_v, pMSS2->V_OutHEID(it_v), num_ccw, (int)it_ccw.IsOpen() );
                        VIZ_POINT2_NAMED( vs, name, p, 5.0f, Vec4f(1,1,1,0.5f) );
                        */
                        /*TEMP: debug it_cw
                          int num_cw(0);
                          MeshSolidShape2::iterator_around_vertex_ccw it_ccw( pMSS2->GetIteratorAroundVertexCCW(it_v) );
                          while( it_ccw.IsValid() )
                          {
                          ++num_cw;
                          --it_ccw;
                          }
                          sprintf( name, "V%d,HE%d,#P%d,O%d", it_v, pMSS2->V_OutHEID(it_v), num_cw, (int)it_ccw.IsOpen() );
                          VIZ_POINT2_NAMED( vs, name, p, 5.0f, Vec4f(1,1,1,0.5f) );
                        */
                    }

                    /* Layers \todo Consider specific DDF
                    for( unsigned int it_l=0; it_l < pMSS2->GetNumL(); it_l++ )
                    {
                        Vec4f color( 0.25f + 0.5f*(1-(float)it_l/pMSS2->GetNumL()), 0, 0, 1 );
                        unsigned int it_p = pMSS2->L_FirstPID(it_l);
                        for( unsigned int it_pil=0; it_pil < pMSS2->L_NumP(it_l); it_pil++, it_p++ )
                        {
                            if( pMSS2->P_NumEdges(it_p) == 3 )
                            {
                                unsigned int it_he = pMSS2->P_FirstHEID(it_p);
                                Vec2 p0( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS2->HE_Next(it_he);
                                Vec2 p1( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS2->HE_Next(it_he);
                                Vec2 p2( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS2->HE_Next(it_he);
                                VIZ_TRIANGLE2( vs, p0, p1, p2, color, util::eVizStyle_Solid );
                            }
                            else
                            {
                                // Draw layer-id
                                Vec2 barycenter( tr2 * pMSS2->P_Barycenter( it_p, actual_sdof ) );
                                sprintf( name, "L%d", it_l );
                                VIZ_POINT2_NAMED( vs, name, barycenter, 1.0f, Vec4f(1,1,0,0.5f) );
                            }
                        }
                    }
                    */
                }

                // Topology (Triangles)
                if( ddf.Test( eSDDF_Topology ) )
                {
                    for( unsigned int it_p=0; it_p < pMSS2->GetNumP(); it_p++ )
                    {
                        // Barycenter
                        Vec2 barycenter1( tr2 * pMSS2->P_Barycenter( it_p, actual_sdof ) );
                        unsigned int it_he( pMSS2->P_FirstHEID(it_p) );
                        do
                        {
                            feature_index_type npid( pMSS2->HE_RightPID(it_he) );
                            it_he = pMSS2->HE_Next(it_he);
                            if( npid != cInvalidFeatureIndex )
                            {
                                Vec2 barycenter2( tr2 * pMSS2->P_Barycenter( npid, actual_sdof ) );
                                VIZ_SEGMENT2( vs, barycenter1, barycenter2, 1.0f, Vec4f(1,1,0,0.5) );
                            }
                        } while ( it_he != pMSS2->P_FirstHEID(it_p) );
                    }
                }
                // Annotations::DCR
                if( ddf.Test( eSDDF_DCR )
                    && 0 != pMSS2->GetDCR() )
                    VizDCR( pMSS2->GetDCR(), pMSS2, tr2, p_go->GetVecDOF(), vs, ddf );
                // Annotations::BVH
                if( ddf.Test( eSDDF_BVH )
                    && 0 != pMSS2->GetBVH() )
                    GVizBVH( *pMSS2->GetBVH(), vs, ddf );
            }
            break;
        case eShape_TetSolid3:
            {
                const TetSolidShape3* pTSS3( static_cast<const TetSolidShape3*>(pShape) );
                const GObjectSDOF<3,Vec3> *pTSSGO( static_cast<const GObjectSDOF<3,Vec3> *>(p_go) );
                Transform3f tr3( static_cast<const IObject3*>(p_go)->GetTransform() );
                const Vec3 *actual_sdof( pTSSGO->GetVecSDOF() );
                if( 0 == actual_sdof )
                {
                    GEO_LOG_WARNING("Object has no SDOF, using default");
                    actual_sdof = pTSS3->GetVecDefaultSDOF();
                }
                // Vertices
                // Vertices
                if( ddf.Test( eSDDF_Vertices ) )
                {
                    if( ddf.Test( eSDDF_FeatureId ) )
                        for( unsigned int it_v=0; it_v < pTSS3->GetNumV(); it_v++ )
                        {
                            sprintf( name, "%u", it_v );
                            VIZ_POINT3_NAMED( vs, name, tr3 * pTSS3->V_Pos( it_v, actual_sdof ), 4.0f, Vec4f(1,0.5,0.5,1) );
                        }
                    else
                        for( unsigned int it_v=0; it_v < pTSS3->GetNumV(); it_v++ )
                            VIZ_POINT3( vs, tr3 * pTSS3->V_Pos( it_v, actual_sdof ), 5.0f, Vec4f(1,0.25,0.25,1) );
                }
                // Boundary | Triangles
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    for( unsigned int it_bf=0; it_bf < pTSS3->GetNumBF(); it_bf++ )
                    {
                        // Barycenter
                        Vec3 barycenter( tr3 * pTSS3->BF_Barycenter( it_bf, actual_sdof ) );
                        // HE, expanded from barycenter to differentiate from symmetrics
                        Vec3 p0( tr3 * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,0), actual_sdof ) );
                        Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,1), actual_sdof ) );
                        Vec3 p2( tr3 * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,2), actual_sdof ) );
                        // Vertices
                        if( ddf.Test( eSDDF_Vertices ) )
                        {
                            if( ddf.Test( eSDDF_FeatureId ) )
                            {
                                sprintf( name, "%d", pTSS3->BF_VID(it_bf,0) );
                                VIZ_POINT3_NAMED( vs, name, p0, 4.0f, Vec4f(1,0.5,0.5,1) );
                                sprintf( name, "%d", pTSS3->BF_VID(it_bf,1) );
                                VIZ_POINT3_NAMED( vs, name, p1, 4.0f, Vec4f(1,0.5,0.5,1) );
                                sprintf( name, "%d", pTSS3->BF_VID(it_bf,2) );
                                VIZ_POINT3_NAMED( vs, name, p2, 4.0f, Vec4f(1,0.5,0.5,1) );
                            }
                            else
                            {
                                VIZ_POINT3( vs, p0, 5.0f, Vec4f(1,0.25,0.25,1) );
                                VIZ_POINT3( vs, p1, 5.0f, Vec4f(1,0.25,0.25,1) );
                                VIZ_POINT3( vs, p2, 5.0f, Vec4f(1,0.25,0.25,1) );
                            }
                        }
                        // Wire
                        VIZ_TRIANGLE3( vs, p0, p1, p2, 1, Vec4f(1,0,0,1), util::eVizStyle_Wire );
                        // Solid (contracted)
                        const Real cContractionCoeff( 0.85f );
                        p0 = barycenter + cContractionCoeff*(p0-barycenter);
                        p1 = barycenter + cContractionCoeff*(p1-barycenter);
                        p2 = barycenter + cContractionCoeff*(p2-barycenter);
                        VIZ_TRIANGLE3( vs, p0, p1, p2, 1, Vec4f(0.5,0,0,1), util::eVizStyle_Solid );

                        if( ddf.Test( eSDDF_Normals ) )
                            VIZ_SEGMENT3( vs, barycenter, barycenter+0.05f*mal::Normalized(mal::Cross(p1-p0,p2-p0)), 1.0f, Vec4f(1,1,1,1) );

                        /*\todo if( ddf.Test( eSDDF_FeatureId ) )
                        sprintf( name, "(%d,%d,%d)", pTSS3->BF_VID(it_bf,0), pTSS3->BF_VID(it_bf,1), pTSS3->BF_VID(it_bf,2) );
                        VIZ_POINT3_NAMED( vs, name, barycenter, 4.0f, Vec4f(1,0.5,0.5,1) );
                        */
                    }
                }
                // Interior (Tetrahedrons)
                if( ddf.Test( eSDDF_Interior ) )
                {
                    for( unsigned int it_tet=0; it_tet < pTSS3->GetNumT(); it_tet++ )
                    {
                        // Gather VID
                        Vec3 barycenter( tr3 * pTSS3->T_Barycenter( it_tet, actual_sdof ) );
                        Vec3 p0( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,0), actual_sdof ) );
                        Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,1), actual_sdof ) );
                        Vec3 p2( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,2), actual_sdof ) );
                        Vec3 p3( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,3), actual_sdof ) );
                        // Vertices
                        if( ddf.Test( eSDDF_Vertices ) )
                        {
                            if( ddf.Test( eSDDF_FeatureId ) )
                            {
                                sprintf( name, "%d", pTSS3->T_VID(it_tet,0) );
                                VIZ_POINT3_NAMED( vs, name, p0, 4.0f, Vec4f(0.5,1,0.5,1) );
                                sprintf( name, "%d", pTSS3->T_VID(it_tet,1) );
                                VIZ_POINT3_NAMED( vs, name, p1, 4.0f, Vec4f(0.5,1,0.5,1) );
                                sprintf( name, "%d", pTSS3->T_VID(it_tet,2) );
                                VIZ_POINT3_NAMED( vs, name, p2, 4.0f, Vec4f(0.5,1,0.5,1) );
                                sprintf( name, "%d", pTSS3->T_VID(it_tet,3) );
                                VIZ_POINT3_NAMED( vs, name, p3, 4.0f, Vec4f(0.5,1,0.5,1) );
                            }
                            else
                            {
                                VIZ_POINT3( vs, p0, 5.0f, Vec4f(1,0.25,0.25,1) );
                                VIZ_POINT3( vs, p1, 5.0f, Vec4f(1,0.25,0.25,1) );
                                VIZ_POINT3( vs, p2, 5.0f, Vec4f(1,0.25,0.25,1) );
                                VIZ_POINT3( vs, p3, 5.0f, Vec4f(1,0.25,0.25,1) );
                            }
                        }
                        const Real cContractionCoeff( 0.95f );
                        VIZ_TETRAHEDRON3( vs,
                                          barycenter + cContractionCoeff*(p0-barycenter),
                                          barycenter + cContractionCoeff*(p1-barycenter),
                                          barycenter + cContractionCoeff*(p2-barycenter),
                                          barycenter + cContractionCoeff*(p3-barycenter),
                                          //Vec4f(1,0,0,0.25), util::eVizStyle_Solid );
                                          Vec4f(1,0,0,0.25), util::eVizStyle_Wire );
                        if( ddf.Test( eSDDF_FeatureId ) )
                        {
                            sprintf( name, "%u", it_tet );
                            VIZ_POINT3_NAMED( vs, name, barycenter, 2.0f, Vec4f(1,0.5,0.5,1) );
                        }

                        /*\todo if( ddf.Test( eSDDF_FeatureId ) )
                        sprintf( name, "%d = (%d,%d,%d,%d)", it_tet,
                                 pTSS3->T_NTID(it_tet,0)==cInvalidFeatureIndex ? -1:pTSS3->T_NTID(it_tet,0),
                                 pTSS3->T_NTID(it_tet,1)==cInvalidFeatureIndex ? -1:pTSS3->T_NTID(it_tet,1),
                                 pTSS3->T_NTID(it_tet,2)==cInvalidFeatureIndex ? -1:pTSS3->T_NTID(it_tet,2),
                                 pTSS3->T_NTID(it_tet,3)==cInvalidFeatureIndex ? -1:pTSS3->T_NTID(it_tet,3) );
                        VIZ_POINT3_NAMED( vs, name, barycenter, 4.0f, Vec4f(1,0.5,0.5,1) );
                        //TEMP*/
                    }

                    /* Layers \todo Consider specific DDF
                    for( unsigned int it_l=0; it_l < pTSS3->GetNumL(); it_l++ )
                    {
                        Vec4f color( 0.25f + 0.75f*(1-(float)it_l/pTSS3->GetNumL()), 0, 0, 1 );
                        unsigned int it_tet = pTSS3->L_FirstTID(it_l);
                        for( unsigned int it_til=0; it_til < pTSS3->L_NumT(it_l); it_til++, it_tet++ )
                        {
                            Vec3 barycenter( tr3 * pTSS3->T_Barycenter( it_tet, actual_sdof ) );
                            VIZ_SPHERE3( vs, barycenter, Quat::Identity(), 0.01f * (1+(float)it_l/pTSS3->GetNumL()), color, util::eVizStyle_Solid );
                        }
                    }
                    */
                }
                // Topology (Tetrahedrons)
                if( ddf.Test( eSDDF_Topology ) )
                {
                    for( unsigned int it_tet=0; it_tet < pTSS3->GetNumT(); it_tet++ )
                    {
                        // Gather VID
                        Vec3 barycenter1( tr3 * pTSS3->T_Barycenter( it_tet, actual_sdof ) );
                        for( unsigned int it_ntit=0; it_ntit<4; it_ntit++ )
                        {
                            uint32 ntid( pTSS3->T_NTID(it_tet,it_ntit) );
                            if( ntid != cInvalidFeatureIndex
                                && ntid < pTSS3->GetNumT() ) //TEMPORAL: REQUIRED because some baked/loaded TetSolidShape3 do not have compiled topology!!
                            {
                                Vec3 barycenter2( tr3 * pTSS3->T_Barycenter( ntid, actual_sdof ) );
                                VIZ_SEGMENT3( vs, barycenter1, barycenter2, 1.0f, Vec4f(1,1,0,0.5) );
                            }
                            else
                            {
                                Vec3 p0( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,(it_ntit+1)%4), actual_sdof ) );
                                Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,(it_ntit+2)%4), actual_sdof ) );
                                Vec3 p2( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tet,(it_ntit+3)%4), actual_sdof ) );
                                //VIZ_SEGMENT3( vs, barycenter1, mal::Rcp<Real>(3)*(p0+p1+p2), 5.0f, Vec4f(1,1,1,1) );
                                //VIZ_TRIANGLE3( vs, p0, p1, p2, 1, Vec4f(1,1,1,1), util::eVizStyle_Wire );
                                VIZ_SEGMENT3( vs, p0, p1, 1.0f, Vec4f(1,1,1,0.5) );
                                VIZ_SEGMENT3( vs, p1, p2, 1.0f, Vec4f(1,1,1,0.5) );
                                VIZ_SEGMENT3( vs, p2, p0, 1.0f, Vec4f(1,1,1,0.5) );
                            }
                        }
                    }
                }
                // Annotations::DCR
                if( ddf.Test( eSDDF_DCR )
                    && 0 != pTSS3->GetDCR() )
                    VizDCR( pTSS3->GetDCR(), pTSS3, tr3, p_go->GetVecDOF(), vs, ddf );
                // Annotations::BVH
                if( ddf.Test( eSDDF_BVH )
                    && 0 != pTSS3->GetBVH() )
                    GVizBVH( *pTSS3->GetBVH(), vs, ddf );
            }
            break;
        case eShape_TriSurface3:
            {
                const TriSurfaceShape3* pTSS3( static_cast<const TriSurfaceShape3*>(pShape) );
                const GObjectSDOF<3,Vec3> *pTSSGO( static_cast<const GObjectSDOF<3,Vec3> *>(p_go) );
                Transform3f tr3( static_cast<const IObject3*>(p_go)->GetTransform() );
                const Vec3 *actual_sdof( pTSSGO->GetVecSDOF() );
                if( 0 == actual_sdof )
                {
                    GEO_LOG_WARNING("Object has no SDOF, using default");
                    actual_sdof = pTSS3->GetVecDefaultSDOF();
                }
                // Vertices
                if( ddf.Test( eSDDF_Vertices ) )
                {
                    if( ddf.Test( eSDDF_FeatureId ) )
                        for( unsigned int it_v=0; it_v < pTSS3->GetNumV(); it_v++ )
                        {
                            sprintf( name, "%u", it_v );
                            VIZ_POINT3_NAMED( vs, name, tr3 * pTSS3->V_Pos( it_v, actual_sdof ), 4.0f, Vec4f(0.5,1,0.5,1) );
                        }
                    else
                        for( unsigned int it_v=0; it_v < pTSS3->GetNumV(); it_v++ )
                            VIZ_POINT3( vs, tr3 * pTSS3->V_Pos( it_v, actual_sdof ), 5.0f, Vec4f(0.25,1,0.25,1) );
                }
                // Triangles
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    for( unsigned int it_tri=0; it_tri < pTSS3->GetNumT(); it_tri++ )
                    {
                        Vec3 barycenter( tr3 * pTSS3->T_Barycenter( it_tri, actual_sdof ) );
                        Vec3 p0( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,0), actual_sdof ) );
                        Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,1), actual_sdof ) );
                        Vec3 p2( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,2), actual_sdof ) );
                        // Vertices
                        if( ddf.Test( eSDDF_Vertices ) )
                        {
                            if( ddf.Test( eSDDF_FeatureId ) )
                            {
                                sprintf( name, "%d", pTSS3->T_VID(it_tri,0) );
                                VIZ_POINT3_NAMED( vs, name, p0, 4.0f, Vec4f(0.5,1,0.5,1) );
                                sprintf( name, "%d", pTSS3->T_VID(it_tri,1) );
                                VIZ_POINT3_NAMED( vs, name, p1, 4.0f, Vec4f(0.5,1,0.5,1) );
                                sprintf( name, "%d", pTSS3->T_VID(it_tri,2) );
                                VIZ_POINT3_NAMED( vs, name, p2, 4.0f, Vec4f(0.5,1,0.5,1) );
                            }
                            else
                            {
                                VIZ_POINT3( vs, p0, 5.0f, Vec4f(0.25,1,0.25,1) );
                                VIZ_POINT3( vs, p1, 5.0f, Vec4f(0.25,1,0.25,1) );
                                VIZ_POINT3( vs, p2, 5.0f, Vec4f(0.25,1,0.25,1) );
                            }
                        }
                        // Wire
                        VIZ_TRIANGLE3( vs, p0, p1, p2, 1, Vec4f(0,1,0,1), util::eVizStyle_Wire );
                        // Solid (contracted)
                        const Real cContractionCoeff( 0.95f );
                        p0 = barycenter + cContractionCoeff*(p0-barycenter);
                        p1 = barycenter + cContractionCoeff*(p1-barycenter);
                        p2 = barycenter + cContractionCoeff*(p2-barycenter);
                        VIZ_TRIANGLE3( vs, p0, p1, p2, 1, Vec4f(0,0.5,0,1), util::eVizStyle_Solid );
                        // Normals
                        if( ddf.Test( eSDDF_Normals ) )
                            VIZ_SEGMENT3( vs, barycenter, barycenter+0.05f*mal::Normalized(mal::Cross(p1-p0,p2-p0)), 1.0f, Vec4f(1,1,1,1) );
                    }
                }
                // Topology
                if( ddf.Test( eSDDF_Topology ) )
                {
                    for( unsigned int it_tri=0; it_tri < pTSS3->GetNumT(); it_tri++ )
                    {
                        Vec3 barycenter1( tr3 * pTSS3->T_Barycenter( it_tri, actual_sdof ) );
                        for( unsigned int it_ntit=0; it_ntit < 3; it_ntit++ )
                        {
                            uint32 ntid( pTSS3->T_NTID(it_tri,it_ntit) );
                            if( ntid != geo::cTriSurface3_InvalidFeatureIndex )
                            {
                                Vec3 barycenter2( tr3 * pTSS3->T_Barycenter( ntid, actual_sdof ) );
                                VIZ_SEGMENT3( vs, barycenter1, barycenter2, 1.0f, Vec4f(1,1,0,1) );
                            }
                            else
                            {
                                Vec3 p0( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,it_ntit), actual_sdof ) );
                                Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,(it_ntit+1)%3), actual_sdof ) );
                                VIZ_POINT3( vs, p0, 6.0f, Vec4f(1,0.5,0.5,1) );
                                VIZ_POINT3( vs, p1, 6.0f, Vec4f(1,0.5,0.5,1) );
                                VIZ_SEGMENT3( vs, p0, p1, 3.0f, Vec4f(1,1,1,1) );
                            }
                        }
                        /*\todo if( ddf.Test( eSDDF_FeatureId ) )
                        sprintf( name, "%d", it_tri );
                        VIZ_POINT3_NAMED( vs, name, barycenter1, 1.0f, Vec4f(1,0,1,0.5f) );
                        //TEMP
                        */
                    }
                }
                // Slivers
                if( ddf.Test( eSDDF_Interior ) ) //TEMP
                {
                    for( unsigned int it_tri=0; it_tri < pTSS3->GetNumT(); it_tri++ )
                    {
                        Vec3 p0( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,0), actual_sdof ) );
                        Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,1), actual_sdof ) );
                        Vec3 p2( tr3 * pTSS3->V_Pos( pTSS3->T_VID(it_tri,2), actual_sdof ) );
                        Vec3 d01( p1-p0 );
                        Vec3 d02( p2-p0 );
                        Vec3 d12( p2-p1 );
                        Real a( mal::Norm(d01) );
                        Real b( mal::Norm(d02) );
                        Real c( mal::Norm(d12) );
                        Real circumradius( a*b*c / mal::Sqrt( (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c) ) ); //\see http://mathworld.wolfram.com/Circumradius.html
                        if( circumradius / mal::Min(a,mal::Min(b,c)) > g_TTS_SliverThreshold ) //larger ratio => worse sliver (min value = 0.5)
                            VIZ_TRIANGLE3( vs, p0, p1, p2, 4, Vec4f(1,0,1,1), util::eVizStyle_Wire );
                    }
                }
                // \todo SUPPORT trisurface3_feature_index_type in BVH viz Annotations::BVH
                if( ddf.Test( eSDDF_BVH )
                    && 0 != pTSS3->GetBVH() )
                    GVizBVH( *pTSS3->GetBVH(), vs, ddf );
            }
            break;
        default:
            GEO_LOG_WARNING( "geo::VizObject: viz shape type %d not yet implemented", pShape->GetType() );
            break;
        }
    }
}

#ifdef __BY_NOW_VIZSHAPE_IS_DISABLED
void VizShape2( const IShape2* p_shape, const Transform2& tr, const Real* vec_dof, util::VizStream& vs, Flags32 ddf )
{
    if( 0 == p_shape )
    {
        GEO_LOG_WARNING( "VizStream: Null p_shape!" );
        return;
    }
    if( 0 == vec_dof )
    {
        GEO_LOG_WARNING( "VizStream: Null vec_dof!" );
        return;
    }

    if( ddf.Test( eODDF_Axis ) )
    {
        VIZ_TRANSFORM2( vs, tr, Real(1) );
    }

    if( ddf.Test( eODDF_BV ) )
    {
        //\todo Draw BV... which one?
    }

    if( ddf.Test( eODDF_Shape ) )
    {
        switch( p_shape->GetType() )
        {
        case eShape_MeshSolid2:
            {
                const MeshSolidShape2* pMSS2( static_cast<const MeshSolidShape2*>(p_shape) );
                const Vec2* vec_sdof( reinterpret_cast<const Vec2*>( vec_dof ) );
                // Boundary | Edges
                if( ddf.Test( eSDDF_Boundary ) )
                {
                    for( unsigned int it_bp=0; it_bp < pMSS2->GetNumBoundaryP(); it_bp++ )
                    {
                        // Barycenter
                        Vec2 barycenter( tr * pMSS2->BP_Barycenter( it_bp, vec_sdof ) );
                        // HE, expanded from barycenter to differentiate from symmetrics
                        const Real cContractionCoeff( 1.01f );
                        unsigned int it_he( pMSS2->BP_FirstHEID(it_bp) );
                        int num_neighbours(0);
                        Vec2 p0( tr * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), vec_sdof ) );
                        if( ddf.Test( eSDDF_Vertices ) )
                        {
                            VIZ_POINT2( vs, p0, 5.0f, Vec4f(1,0.25,0.25,1) );
                        }
                        do
                        {
                            if( pMSS2->HE_Sym(it_he) != cInvalidFeatureIndex ) num_neighbours++;

                            /*TEMP
                              sprintf( name, "HE%d,S%d", it_he, pMSS2->HE_Sym(it_he) );
                              TEMP*/

                            it_he = pMSS2->HE_Next(it_he);
                            Vec2 p1( tr * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), vec_sdof ) );
                            VIZ_SEGMENT2( vs,
                                          barycenter + cContractionCoeff*(p0-barycenter),
                                          barycenter + cContractionCoeff*(p1-barycenter),
                                          1.0f, Vec4f(1,0,0,1) );

                            if( ddf.Test( eSDDF_Vertices ) )
                            {
                                VIZ_POINT2( vs, p1, 5.0f, Vec4f(1,0.25,0.25,1) );
                            }

                            /*TEMP
                              VIZ_POINT2_NAMED( vs, name,
                              barycenter + 0.5f*0.8f*( cContractionCoeff*(p0-barycenter) + cContractionCoeff*(p1-barycenter) ),
                              2.0f, Vec4f(1,0,1,0.5f) );
                              TEMP*/

                            p0 = p1;
                        } while ( it_he != pMSS2->BP_FirstHEID(it_bp) );

                        /*TEMP: Num neighbours
                          sprintf( name, "P%d,#N%d", it_bp, num_neighbours );
                          VIZ_POINT2_NAMED( vs, name, barycenter, 5.0f, Vec4f(1,1,0,0.5f) );
                        */
                    }
                }

                // Interior (Polygons)
                if( ddf.Test( eSDDF_Interior ) )
                {
                    for( unsigned int it_p=0; it_p < pMSS2->GetNumP(); it_p++ )
                    {
                        // Barycenter
                        Vec2 barycenter( tr * pMSS2->P_Barycenter( it_p, vec_sdof ) );
                        //VIZ_DISK2_NO_ROT( vs, barycenter, 0.01, Vec4f(1,1,0,1), util::eVizStyle_Solid );

                        // HE, contracted towards barycenter to differentiate from symmetrics
                        const Real cContractionCoeff( 0.95f );
                        unsigned int it_he( pMSS2->P_FirstHEID(it_p) );
                        int num_neighbours(0);
                        Vec2 p0( tr * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), vec_sdof ) );
                        do
                        {
                            if( pMSS2->HE_Sym(it_he) != cInvalidFeatureIndex ) num_neighbours++;

                            /*TEMP
                              sprintf( name, "HE%d,S%d", it_he, pMSS2->HE_Sym(it_he) );
                              TEMP*/

                            it_he = pMSS2->HE_Next(it_he);
                            Vec2 p1( tr * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), vec_sdof ) );
                            VIZ_SEGMENT2( vs,
                                          barycenter + cContractionCoeff*(p0-barycenter),
                                          barycenter + cContractionCoeff*(p1-barycenter),
                                          1.0f, Vec4f(1,0,0,0.5) );
                            /*TEMP
                              VIZ_POINT2_NAMED( vs, name,
                              barycenter + 0.5f*0.8f*( cContractionCoeff*(p0-barycenter) + cContractionCoeff*(p1-barycenter) ),
                              2.0f, Vec4f(1,0,1,0.5f) );
                              TEMP*/

                            p0 = p1;
                        } while ( it_he != pMSS2->P_FirstHEID(it_p) );

                        /*TEMP: Num neighbours
                          sprintf( name, "P%d,#N%d", it_p, num_neighbours );
                          VIZ_POINT2_NAMED( vs, name, barycenter, 5.0f, Vec4f(1,1,0,0.5f) );
                        */
                    }

                    // Interior Vertices
                    for( unsigned int it_v=0; it_v < pMSS2->GetNumV(); it_v++ )
                    {
                        Vec2 p( tr * pMSS2->V_Pos( it_v, vec_sdof ) );
                        VIZ_POINT2( vs, p, 5.0f, Vec4f(1,0.25,0.25,0.5) );
                        /*TEMP: Named
                          sprintf( name, "V%d", it_v );
                          VIZ_POINT2_NAMED( vs, name, p, 5.0f, Vec4f(1,1,1,0.5f) );
                        */

                        /* iterator_around_vertex_ccw
                        //TEMP: debug it_ccw
                        int num_ccw(0);
                        MeshSolidShape2::iterator_around_vertex_ccw it_ccw( pMSS2->GetIteratorAroundVertexCCW(it_v) );
                        while( it_ccw.IsValid() )
                        {
                        ++num_ccw;
                        ++it_ccw;
                        }
                        sprintf( name, "V%d,HE%d,#P%d,O%d", it_v, pMSS2->V_OutHEID(it_v), num_ccw, (int)it_ccw.IsOpen() );
                        VIZ_POINT2_NAMED( vs, name, p, 5.0f, Vec4f(1,1,1,0.5f) );
                        */
                        /*TEMP: debug it_cw
                          int num_cw(0);
                          MeshSolidShape2::iterator_around_vertex_ccw it_ccw( pMSS2->GetIteratorAroundVertexCCW(it_v) );
                          while( it_ccw.IsValid() )
                          {
                          ++num_cw;
                          --it_ccw;
                          }
                          sprintf( name, "V%d,HE%d,#P%d,O%d", it_v, pMSS2->V_OutHEID(it_v), num_cw, (int)it_ccw.IsOpen() );
                          VIZ_POINT2_NAMED( vs, name, p, 5.0f, Vec4f(1,1,1,0.5f) );
                        */
                    }
                }
            }
        default:
            GEO_LOG_WARNING( "geo::VizShape2: viz shape type %d not yet implemented", p_shape->GetType() );
            break;
        }
    }
}
#endif //__BY_NOW_VIZSHAPE_IS_DISABLED

void VizBoundingVolume( const bv::IBoundingVolume* p_bv, util::VizStream& vs, Flags32 ddf )
{
    Vec4f color(1,1,0,0.25f);
    switch( p_bv->GetType() )
    {
    case bv::eBV_Sphere2:
        {
            const bv::Sphere2 &sphere2( p_bv->As<bv::Sphere2>() );
            VIZ_DISK2_NO_ROT( vs, sphere2.GetPos(), sphere2.GetRadius(), color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_Sphere3:
        {
            const bv::Sphere3 &sphere3( p_bv->As<bv::Sphere3>() );
            VIZ_SPHERE3_NO_ROT( vs, sphere3.GetPos(), sphere3.GetRadius(), color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_AABB2:
        {
            const bv::AABB2 &aabb2( p_bv->As<bv::AABB2>() );
            VIZ_RECT2_NO_ROT( vs, aabb2.GetPos(), aabb2.GetHalfSizes(), color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_AABB3:
        {
            const bv::AABB3 &aabb3( p_bv->As<bv::AABB3>() );
            VIZ_BOX3_NO_ROT( vs, aabb3.GetPos(), aabb3.GetHalfSizes(), color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_LSS2:
        {
            const bv::LSS2 &lss2( p_bv->As<bv::LSS2>() );
            VIZ_DISK2_NO_ROT( vs, lss2.GetPos0(), lss2.GetRadius(), color, util::eVizStyle_Wire );
            VIZ_SEGMENT2( vs, lss2.GetPos0(), lss2.GetPos1(), 5.0f, color );
            VIZ_DISK2_NO_ROT( vs, lss2.GetPos1(), lss2.GetRadius(), color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_LSS3:
        {
            const bv::LSS3 &lss3( p_bv->As<bv::LSS3>() );
            VIZ_SPHERE3_NO_ROT( vs, lss3.GetPos0(), lss3.GetRadius(), color, util::eVizStyle_Wire );
            VIZ_SEGMENT3( vs, lss3.GetPos0(), lss3.GetPos1(), 5.0f, color );
            VIZ_SPHERE3_NO_ROT( vs, lss3.GetPos1(), lss3.GetRadius(), color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_DOP2_K4:
        {
            // \note Special case for AABB
            const bv::DOP2_K4 &dop2k4( p_bv->As<bv::DOP2_K4>() );
            VIZ_RECT2_NO_ROT( vs,
                              Vec2( dop2k4.GetInterval(0).Mid(),
                                    dop2k4.GetInterval(1).Mid() ),
                              Vec2( Real(0.5)*dop2k4.GetInterval(0).Length(),
                                    Real(0.5)*dop2k4.GetInterval(1).Length() ),
                              color, util::eVizStyle_Wire );
        }
        break;
    case bv::eBV_DOP2_K8:
        {
            // Get KDOP segments
            const bv::DOP2_K8& dop2k8( p_bv->As<bv::DOP2_K8>() );
            std::pair<Vec2,Vec2> vec_segments[8];
            unsigned int num_segments = bv::DOP2_K8_Segments( dop2k8, vec_segments );
            // Viz clipped segments, regardless of actual order
            for( unsigned int it_segment=0; it_segment < num_segments; it_segment++ )
                VIZ_SEGMENT2( vs, vec_segments[it_segment].first,
                              vec_segments[it_segment].second,
                              1.0f, color );
        }
        break;
    case bv::eBV_Void:
    case bv::eBV_Infinite:
    default:
        GEO_LOG_WARNING( "geo::VizBoundingVolume: viz BV type %d not yet implemented", p_bv->GetType() );
        break;
    }
}

void VizBroadPhase( const bp::IBroadPhase *p_bp, util::VizStream &vs, Flags32 ddf )
{
    switch( p_bp->GetType() )
    {
    case bp::eBP_Basic:
    case bp::eBP_QuadraticTime:
        for( const bp::Proxy *pProxy = p_bp->FirstProxy(); 0 != pProxy; pProxy = p_bp->NextProxy() )
            VizBoundingVolume( pProxy->m_pBV, vs );
        break;
        //\todo case bp::eBP_SweepAndPrune:
    default:
        GEO_LOG_WARNING( "geo::VizBroadPhase: viz BP type %d not yet implemented", p_bp->GetType() );
        break;
    }
}

void VizContactData( const geo::np::ContactData2 &cd, util::VizStream &vs, Flags32 ddf )
{
    if( cd.Size() > 0 )
    {
        Vec2 avg_pos(0,0);
        for( unsigned int it_cp=0; it_cp<cd.Size(); it_cp++ )
        {
            // Pos
            if( ddf.Test( eContactDataDDF_Points ) )
            {
                VIZ_POINT2( vs, cd.GetCP(it_cp).m_Pos1, 5.0f, Vec4f(0,1,1,1) );
                VIZ_POINT2( vs, cd.GetCP(it_cp).m_Pos2, 5.0f, Vec4f(1,0,1,1) );
            }
            // Normal
            if( ddf.Test( eContactDataDDF_Normals ) )
                VIZ_VECTOR2( vs, cd.GetCP(it_cp).m_Pos1, cd.GetCP(it_cp).m_Normal * cd.GetCP(it_cp).m_Depth, 2.0f, Vec4f(0,1,0,1) );
            avg_pos += cd.GetCP(it_cp).m_Pos1 + cd.GetCP(it_cp).m_Pos2;
        }
        // Avg Normal*Depth
        if( ddf.Test( eContactDataDDF_AvgNormal ) )
            VIZ_VECTOR2( vs, avg_pos / cd.Size(), cd.m_AvgNormal * cd.m_AvgDepth, 2.0f, Vec4f(0.5,1,0.5,1) );
    }
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
    if( ddf.Test( eContactDataDDF_Stochastic_IP ) )
        for( unsigned int it_ip=0; it_ip<cd.m_VD.m_vecIP.size(); it_ip++ )
            VIZ_POINT2( vs, 0.5f*(cd.m_VD.m_vecIP[it_ip].first + cd.m_VD.m_vecIP[it_ip].second), 8, Vec4f(1,0,0,1) );
    if( ddf.Test( eContactDataDDF_Stochastic_CP ) )
        for( unsigned int it_cp=0; it_cp<cd.m_VD.m_vecCP.size(); it_cp++ )
        {
            VIZ_POINT2( vs, cd.m_VD.m_vecCP[it_cp].first, 6, Vec4f(1,0.5,0,1) );
            VIZ_POINT2( vs, cd.m_VD.m_vecCP[it_cp].second, 6, Vec4f(1,0.5,0,1) );
        }
    if( ddf.Test( eContactDataDDF_Stochastic_NP ) )
    {
        for( unsigned int it_np=0; it_np<cd.m_VD.m_vecNP.size(); it_np++ )
        {
            VIZ_POINT2( vs, cd.m_VD.m_vecNP[it_np].first, 4, Vec4f(0,1,0,1) );
            VIZ_POINT2( vs, cd.m_VD.m_vecNP[it_np].second, 4, Vec4f(0,1,0,1) );
        }
        for( unsigned int it_fp=0; it_fp<cd.m_VD.m_vecFP.size(); it_fp++ )
        {
            VIZ_POINT2( vs, cd.m_VD.m_vecFP[it_fp].first, 3, Vec4f(0.25,0.25,1,1) );
            VIZ_POINT2( vs, cd.m_VD.m_vecFP[it_fp].second, 3, Vec4f(0.25,0.25,1,1) );
        }
    }
    if( ddf.Test( eContactDataDDF_Stochastic_RNP ) )
        for( unsigned int it_rnp=0; it_rnp<cd.m_VD.m_vecRNP.size(); it_rnp++ )
        {
            VIZ_POINT2( vs, cd.m_VD.m_vecRNP[it_rnp].first, 2, Vec4f(0.5,0.5,0.5,1) );
            VIZ_POINT2( vs, cd.m_VD.m_vecRNP[it_rnp].second, 2, Vec4f(1,1,1,1) );
        }
    if( ddf.Test( eContactDataDDF_Stochastic_PCA ) )
        for( unsigned int it_pca=0; it_pca<cd.m_VD.m_vecPCA.size(); it_pca++ )
            VIZ_SEGMENT2( vs,
                          cd.m_VD.m_vecPCA[it_pca].first,
                          cd.m_VD.m_vecPCA[it_pca].first + cd.m_VD.m_vecPCA[it_pca].second,
                          2, Vec4f(1,1,0,1) );
    /*IM
    for( unsigned int it_im_p=0; it_im_p<cd.m_VD.m_vecIM_XPoint.size(); it_im_p++ )
        VIZ_POINT2( vs, cd.m_VD.m_vecIM_XPoint[it_im_p], 8, Vec4f(0.5,0.5,0.5,1) );
    for( unsigned int it_im_s=0; it_im_s<cd.m_VD.m_vecIM_XSegment.size(); it_im_s++ )
        VIZ_SEGMENT2( vs,
                      cd.m_VD.m_vecIM_XSegment[it_im_s].first,
                      cd.m_VD.m_vecIM_XSegment[it_im_s].second,
                      4, Vec4f(1,1,0,1) );
    */
#endif
}

void VizContactData( const geo::np::ContactData3 &cd, util::VizStream &vs, Flags32 ddf )
{
    if( cd.Size() > 0 )
    {
        Vec3 avg_pos( Vec3::Zero() );
        for( unsigned int it_cp=0; it_cp<cd.Size(); it_cp++ )
        {
            // Pos
            if( ddf.Test( eContactDataDDF_Points ) )
            {
                VIZ_POINT3( vs, cd.GetCP(it_cp).m_Pos1, 5.0f, Vec4f(0,1,1,1) );
                VIZ_POINT3( vs, cd.GetCP(it_cp).m_Pos2, 5.0f, Vec4f(1,0,1,1) );
            }
            // Normal
            if( ddf.Test( eContactDataDDF_Normals ) )
                VIZ_VECTOR3( vs, cd.GetCP(it_cp).m_Pos1, cd.GetCP(it_cp).m_Normal * cd.GetCP(it_cp).m_Depth, 2.0f, Vec4f(0,1,0,1) );
            avg_pos += cd.GetCP(it_cp).m_Pos1 + cd.GetCP(it_cp).m_Pos2;
        }
        // Avg Normal*Depth
        if( ddf.Test( eContactDataDDF_AvgNormal ) )
            VIZ_VECTOR3( vs, avg_pos / cd.Size(), cd.m_AvgNormal * cd.m_AvgDepth, 2.0f, Vec4f(0.5,1,0.5,1) );
    }
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
    if( ddf.Test( eContactDataDDF_Stochastic_IP ) )
        for( unsigned int it_ip=0; it_ip<cd.m_VD.m_vecIP.size(); it_ip++ )
            VIZ_POINT3( vs, 0.5f*(cd.m_VD.m_vecIP[it_ip].first + cd.m_VD.m_vecIP[it_ip].second), 8, Vec4f(1,0,0,1) );
    if( ddf.Test( eContactDataDDF_Stochastic_CP ) )
        for( unsigned int it_cp=0; it_cp<cd.m_VD.m_vecCP.size(); it_cp++ )
        {
            VIZ_POINT3( vs, cd.m_VD.m_vecCP[it_cp].first, 6, Vec4f(1,0.5,0,1) );
            VIZ_POINT3( vs, cd.m_VD.m_vecCP[it_cp].second, 6, Vec4f(1,0.5,0,1) );
        }
    if( ddf.Test( eContactDataDDF_Stochastic_NP ) )
    {
        for( unsigned int it_np=0; it_np<cd.m_VD.m_vecNP.size(); it_np++ )
        {
            VIZ_POINT3( vs, cd.m_VD.m_vecNP[it_np].first, 4, Vec4f(0,1,0,1) );
            VIZ_POINT3( vs, cd.m_VD.m_vecNP[it_np].second, 4, Vec4f(0,1,0,1) );
        }
        for( unsigned int it_fp=0; it_fp<cd.m_VD.m_vecFP.size(); it_fp++ )
        {
            VIZ_POINT3( vs, cd.m_VD.m_vecFP[it_fp].first, 3, Vec4f(0.25,0.25,1,1) );
            VIZ_POINT3( vs, cd.m_VD.m_vecFP[it_fp].second, 3, Vec4f(0.25,0.25,1,1) );
        }
    }
    if( ddf.Test( eContactDataDDF_Stochastic_RNP ) )
        for( unsigned int it_rnp=0; it_rnp<cd.m_VD.m_vecRNP.size(); it_rnp++ )
        {
            VIZ_POINT3( vs, cd.m_VD.m_vecRNP[it_rnp].first, 2, Vec4f(0.5,0.5,0.5,1) );
            VIZ_POINT3( vs, cd.m_VD.m_vecRNP[it_rnp].second, 2, Vec4f(1,1,1,1) );
        }
    if( ddf.Test( eContactDataDDF_Stochastic_PCA ) )
        for( unsigned int it_pca=0; it_pca<cd.m_VD.m_vecPCA.size(); it_pca++ )
            VIZ_SEGMENT3( vs,
                          cd.m_VD.m_vecPCA[it_pca].first,
                          cd.m_VD.m_vecPCA[it_pca].first + cd.m_VD.m_vecPCA[it_pca].second,
                          2, Vec4f(1,1,0,1) );
    /*IM
    for( unsigned int it_im_p=0; it_im_p<cd.m_VD.m_vecIM_XPoint.size(); it_im_p++ )
        VIZ_POINT3( vs, cd.m_VD.m_vecIM_XPoint[it_im_p], 8, Vec4f(0.5,0.5,0.5,1) );
    for( unsigned int it_im_s=0; it_im_s<cd.m_VD.m_vecIM_XSegment.size(); it_im_s++ )
        VIZ_SEGMENT3( vs,
                      cd.m_VD.m_vecIM_XSegment[it_im_s].first,
                      cd.m_VD.m_vecIM_XSegment[it_im_s].second,
                      4, Vec4f(1,1,0,1) );
    */
#endif
}

void VizDCR( const DCR_MeshSolidShape2* p_dcr,
             const MeshSolidShape2* p_mss, const Transform2& tr, const Real* vec_dof,
             util::VizStream &vs, Flags32 ddf )
{
//#define __ENABLE_DCR_DRAW_BV
//#define __ENABLE_DCR_COMPARE_BV //\todo THIS SHOULD NOT BE IN VIZ...
#ifdef __ENABLE_DCR_COMPARE_BV
    Real total_area_aabb_E(0);
    Real total_area_kdop_E(0);
    Real total_area_BSlab(0);
    Real total_area_aabb_BSlab(0);
    Real total_area_kdop_BSlab(0);
    Real total_area_aabb_BDOP(0);
    Real total_area_kdop_BDOP(0);
#endif
    std::vector<Vec2> alloc_v(p_dcr->m_NumVertices); //\todo Scope-allocation, consider optimizing with a Viz-specific scratchpad!
    Vec2* vec_v = &alloc_v[0]; //\todo p_context->m_ScratchPad.NewPOD<Vec2>(pDCR1->m_NumVertices); //POD => Don't call ctor/dtor
    const Vec2* default_sdof( p_mss->GetVecDefaultSDOF() );
    const Vec2* actual_sdof( reinterpret_cast<const Vec2*>(vec_dof) );
    if( 0 == actual_sdof ) actual_sdof = default_sdof;
    // Apply baricentric transform and tr to all V
    for( unsigned int it_e=0; it_e<p_dcr->m_NumElements; it_e++ )
    {
        const DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[it_e] );
        uint32 vec_nid[3];
        p_mss->P_VecVID( it_e, vec_nid, 3 );
        Mat2x2 B( mal::GMat2x2_From_Columns(actual_sdof[vec_nid[1]]-actual_sdof[vec_nid[0]],
                                            actual_sdof[vec_nid[2]]-actual_sdof[vec_nid[0]]) );
        Mat2x2 B0( mal::GMat2x2_From_Columns(default_sdof[vec_nid[1]]-default_sdof[vec_nid[0]],
                                             default_sdof[vec_nid[2]]-default_sdof[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat2x2 B_InvB0( B * mal::Inverse(B0) );
        /*\note Direct but slower code
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
            vec_v[ed.m_FirstVID + it_vie] = tr * ( B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]])
                                                    + actual_sdof[vec_nid[0]] );
        */
        // Optimized version
        Transform2 tr_B_InvB0( tr.m_Pos, tr.m_Rot * B_InvB0 );
        Vec2 p0_0( tr.m_Rot * actual_sdof[vec_nid[0]] );
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
            vec_v[ed.m_FirstVID + it_vie] = tr_B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]]) + p0_0;

        //TEMP: Draw BDOP and AABB(BDOP)
        if( !ed.m_BDOP[0].IsEmpty()
            && !ed.m_BDOP[1].IsEmpty()
            && !ed.m_BDOP[2].IsEmpty() )
        {
            bv::AABB2 aabb_E, aabb_BSlab, aabb_BDOP;
            GEBV_MeshSolidShape2_E( p_mss, tr, actual_sdof, feature_index_type(it_e), aabb_E );
            GEBV_MeshSolidShape2_BSlab_Safe( p_mss, tr, actual_sdof, feature_index_type(it_e), aabb_BSlab );
            GEBV_MeshSolidShape2_BDOP_Safe( p_mss, tr, actual_sdof, feature_index_type(it_e), aabb_BDOP );
            bv::DOP2_K8 dop2k8_E, dop2k8_BSlab, dop2k8_BDOP;
            GEBV_MeshSolidShape2_E( p_mss, tr, actual_sdof, feature_index_type(it_e), dop2k8_E );
            GEBV_MeshSolidShape2_BSlab_Safe( p_mss, tr, actual_sdof, feature_index_type(it_e), dop2k8_BSlab );
            GEBV_MeshSolidShape2_BDOP_Safe( p_mss, tr, actual_sdof, feature_index_type(it_e), dop2k8_BDOP );
#ifdef __ENABLE_DCR_COMPARE_BV
            total_area_aabb_E += bv::ComputeVolume(aabb_E);
            total_area_kdop_E += bv::ComputeVolume(dop2k8_E);
            total_area_aabb_BSlab += bv::ComputeVolume(aabb_BSlab);
            total_area_kdop_BSlab += bv::ComputeVolume(dop2k8_BSlab);
            total_area_aabb_BDOP += bv::ComputeVolume(aabb_BDOP);
            total_area_kdop_BDOP += bv::ComputeVolume(dop2k8_BDOP);
#endif
            Vec2 bslab_pos[4];
            {
                // compute global node pos
                Vec2 node_pos[3] = { tr*actual_sdof[vec_nid[0]], tr*actual_sdof[vec_nid[1]], tr*actual_sdof[vec_nid[2]] };
                int i( ed.m_BDOP_BestSlabIdx );
                bslab_pos[0] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+1)%3];
                bslab_pos[1] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+2)%3];
                bslab_pos[2] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+1)%3];
                bslab_pos[3] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+2)%3];
#ifdef __ENABLE_DCR_COMPARE_BV
                // Trapezoid area = h*(min(l1,l2) + 1/2*abs(l1-l2))
                Real h = mal::Abs( mal::Dot( bslab_pos[0]-bslab_pos[2],
                                             mal::Normalized( mal::PerpendicularCCW(bslab_pos[0]-bslab_pos[1]) ) ) );
                total_area_BSlab = h * ( mal::Min( mal::Norm(bslab_pos[0]-bslab_pos[1]),
                                                   mal::Norm(bslab_pos[2]-bslab_pos[3]) )
                                         + 0.5f*mal::Abs( mal::Norm(bslab_pos[0]-bslab_pos[1]) - mal::Norm(bslab_pos[2]-bslab_pos[3]) ) );
#endif
                Vec4f bslab_color = Vec4f(1,0.5,0.5,1);
                bslab_color[i] = 1;
#ifdef __ENABLE_DCR_DRAW_BV
                VIZ_SEGMENT2( vs, bslab_pos[0], bslab_pos[1], 1.0f, bslab_color );
                VIZ_SEGMENT2( vs, bslab_pos[2], bslab_pos[3], 2.0f, bslab_color );
#endif
            }
#ifdef __ENABLE_DCR_DRAW_BV
            /*TEMP
            VizBoundingVolume( &aabb_E, vs );
            VizBoundingVolume( &dop2k8_E, vs );
            VizBoundingVolume( &aabb_BSlab, vs );
            VizBoundingVolume( &dop2k8_BSlab, vs );
            */
            VizBoundingVolume( &aabb_BDOP, vs );
            VizBoundingVolume( &dop2k8_BDOP, vs );

            /*
            Vec2 barycenter( 0.3333f*(actual_sdof[vec_nid[0]]+actual_sdof[vec_nid[1]]+actual_sdof[vec_nid[2]]) );
            VIZ_POINT2( vs, tr*barycenter, 5.0f, Vec4f(1,0,0,1) );
            VIZ_POINT2_NAMED( vs, "n0", tr*actual_sdof[vec_nid[0]], 1.0f, Vec4f(1,1,1,1) );
            VIZ_POINT2_NAMED( vs, "n1", tr*actual_sdof[vec_nid[1]], 1.0f, Vec4f(1,1,1,1) );
            VIZ_POINT2_NAMED( vs, "n2", tr*actual_sdof[vec_nid[2]], 1.0f, Vec4f(1,1,1,1) );
            */
#endif
        }
    }
#ifdef __ENABLE_DCR_COMPARE_BV
    GEO_LOG_WARNING("Areas:\n\tAABB(E)     %f\n\tKDOP(E)     %f\n\tBSlab       %f\n\tAABB(BSlab) %f\n\tKDOP(BSlab) %f\n\tKDOP(BDOP)  %f\n\tKDOP(BDOP)  %f",
                    total_area_aabb_E, total_area_kdop_E,
                    total_area_BSlab, total_area_aabb_BSlab, total_area_kdop_BSlab,
                    total_area_aabb_BDOP, total_area_kdop_BDOP );
#endif

    // Per-patch segments
    mal::SetRandomSeed(666);
    for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
    {
        Vec4f color( mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) );
        const geo::DCR_MeshSolidShape2::Patch& pd( p_dcr->m_vecP[it_patch] );
        for( unsigned int it_sid=0; it_sid < pd.m_NumSegments; it_sid++ )
        {
            unsigned int sid( pd.m_FirstSID + it_sid );
            VIZ_SEGMENT2( vs,
                          vec_v[p_dcr->m_vecS[sid].GetVID(0)],
                          vec_v[p_dcr->m_vecS[sid].GetVID(1)],
                          3.0f, color );
            //TEMP: show segment direction
            // Vec2 n( mal::Normalized( mal::PerpendicularCW(vec_v[p_dcr->m_vecS[sid].GetVID(1)]-vec_v[p_dcr->m_vecS[sid].GetVID(0)] ) ) );
            // VIZ_TRIANGLE2( vs,
            //                vec_v[p_dcr->m_vecS[sid].GetVID(0)]-n*0.01,
            //                vec_v[p_dcr->m_vecS[sid].GetVID(0)]+n*0.01,
            //                vec_v[p_dcr->m_vecS[sid].GetVID(1)],
            //                1, color, util::eVizStyle_Solid );
        }
    }
    /* Vertices
       for( unsigned int it_v=0; it_v < m_pDCR->m_NumVertices; it_v++ )
       VIZ_POINT2( m_VizIS, translation + m_pDCR->m_vecV[it_v], 5.0f, Vec4f(1,1,1,1) );
    */
    /* Topology
       for( unsigned int it_segment=0; it_segment < m_pDCR->m_NumSegments; it_segment++ )
       {
       const geo::DCR_MeshSolidShape2::patch_segment_topology& sd( m_pDCR->m_vecS[it_segment] );
       geo::Vec2 lp0 = m_pDCR->m_vecV[ m_pDCR->m_vecS[sd.m_vecNeighbourSID[0]].GetVID(1) ];
       geo::Vec2 lD = m_pDCR->m_vecV[ sd.GetVID(1) ] - m_pDCR->m_vecV[ sd.GetVID(0) ];
       geo::Vec2 lN = mal::PerpendicularCCW( lD );
       geo::Vec2 lp = 0.5 * ( m_pDCR->m_vecV[ sd.GetVID(0) ] + m_pDCR->m_vecV[ sd.GetVID(1) ] ) + 0.1f * lN;
       geo::Vec2 lp1 = m_pDCR->m_vecV[ m_pDCR->m_vecS[sd.m_vecNeighbourSID[1]].GetVID(0) ];
       VIZ_SEGMENT2( m_VizIS, translation + lp0, translation + lp, 2.0f, Vec4f(1,1,0,1) );
       VIZ_SEGMENT2( m_VizIS, translation + lp, translation + lp1, 2.0f, Vec4f(0,1,1,1) );
       }
    */

#ifdef __ENABLE_MSS_DCR_PATCH_BDT
    // Draw DCR BST
    if( true )//ddf.Test( eODDF_BV ) )
    {
        mal::SetRandomSeed(666);
        for( unsigned int it_p=0; it_p < p_dcr->m_NumPatches; it_p++ )
        {
            Vec4f color( mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) );
            const geo::DCR_MeshSolidShape2::Patch& pd( p_dcr->m_vecP[it_p] );
            bv::BDOP3 root_bdop( p_dcr->m_vecE[pd.m_EID].m_BDOP );
            typedef std::pair< DCR_MeshSolidShape2::Patch::bdt_node_sip_range,uint32> stack_entry_type;
            std::vector< stack_entry_type > stackBDTN;
            stackBDTN.push_back( stack_entry_type(DCR_MeshSolidShape2::Patch::bdt_node_sip_range(0,pd.m_NumSegments),0) );
            // GEO_LOG("DCR.P[%u].BDT",it_p);
            while( !stackBDTN.empty() )
            {
                // Pop PS subarray
                stack_entry_type se( stackBDTN.back() );
                stackBDTN.pop_back();
                DCR_MeshSolidShape2::Patch::bdt_node_sip_range node_sr( se.first );
                // Get node
                const DCR_MeshSolidShape2::Segment& bdtn( p_dcr->m_vecS[pd.m_FirstSID+node_sr.first] );
                bv::BDOP3 bdop( bdtn.BDTN_BDOPq() );

                //TEMP
                // {
                //     GEO_LOG("BDTN level %u, volume %f", se.second, bdop[0].Length()*bdop[1].Length()*bdop[2].Length() );
                //     GEO_LOG("- BDOP [%f,%f]x[%f,%f]x[%f,%f]", bdop[0].Min(), bdop[0].Max(), bdop[1].Min(), bdop[1].Max(), bdop[2].Min(), bdop[2].Max() );
                // }

// #define __ENABLE_DCR_DRAW_BDT
#ifdef __ENABLE_DCR_DRAW_BDT
                /* Draw BDOP[0] slab
                   \todo This is pretty useful but it'd be better
                   if we actually (1) computed height01 as a global
                   height fraction (2) computed exact BDOP segments as
                   we do with DOP2_K8
                */
                //if( se.second == 0 )
                {
                    uint32 vec_nid[3];
                    p_mss->P_VecVID( pd.m_EID, vec_nid, 3 );
                    Vec2 node_pos[3] = { tr*actual_sdof[vec_nid[0]], tr*actual_sdof[vec_nid[1]], tr*actual_sdof[vec_nid[2]] };
                    for( int i=0; i<3; i++ )
                    {
                        Vec2 bslab_pos[4];
                        bslab_pos[0] = bdop[i].Min() * node_pos[i] + (1-bdop[i].Min()) * node_pos[(i+1)%3];
                        bslab_pos[1] = bdop[i].Min() * node_pos[i] + (1-bdop[i].Min()) * node_pos[(i+2)%3];
                        bslab_pos[2] = bdop[i].Max() * node_pos[i] + (1-bdop[i].Max()) * node_pos[(i+1)%3];
                        bslab_pos[3] = bdop[i].Max() * node_pos[i] + (1-bdop[i].Max()) * node_pos[(i+2)%3];
                        // VIZ_SEGMENT2( vs, bslab_pos[0], bslab_pos[1], 1.0f, color );
                        // VIZ_SEGMENT2( vs, bslab_pos[2], bslab_pos[3], 2.0f, color );
                        Real height( 0.1f * se.second );
                        VIZ_SEGMENT3( vs, mal::Concat(bslab_pos[0],height), mal::Concat(bslab_pos[1],height), 1.0f, color );
                        VIZ_SEGMENT3( vs, mal::Concat(bslab_pos[2],height), mal::Concat(bslab_pos[3],height), 2.0f, color );
                    }
                }
#endif

                // Recurse
                int length_sr( node_sr.second - node_sr.first );
                int remaining_sr( length_sr - 1 );
                if( remaining_sr > 1 )
                {
                    stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1, node_sr.first+1+remaining_sr/2 ), se.second+1 ) );
                    stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1+remaining_sr/2, node_sr.second ), se.second+1) );
                }
                else if( remaining_sr == 1 )
                    stackBDTN.push_back( stack_entry_type( DCR_MeshSolidShape2::Patch::bdt_node_sip_range( node_sr.first+1, node_sr.first+2 ), se.second+1 ) );
            }
        }
    }
#endif
}

void VizDCR( const DCR_TetSolidShape3* p_dcr,
             const TetSolidShape3* p_tss, const Transform3& tr, const Real* vec_dof,
             util::VizStream &vs, Flags32 ddf )
{
    std::vector<Vec3> alloc_v(p_dcr->m_NumVertices); //\todo Scope-allocation, consider optimizing with a Viz-specific scratchpad!
    Vec3* vec_v = &alloc_v[0]; //\todo p_context->m_ScratchPad.NewPOD<Vec3>(pDCR1->m_NumVertices); //POD => Don't call ctor/dtor
    const Vec3* default_sdof( p_tss->GetVecDefaultSDOF() );
    const Vec3* actual_sdof( reinterpret_cast<const Vec3*>(vec_dof) );
    if( 0 == actual_sdof ) actual_sdof = default_sdof;
    // Apply baricentric transform and tr to all V
    for( unsigned int it_e=0; it_e<p_dcr->m_NumElements; it_e++ )
    {
        const DCR_TetSolidShape3::Element& ed( p_dcr->m_vecE[it_e] );
        uint32 vec_nid[4] = { p_tss->T_VID(it_e,0), p_tss->T_VID(it_e,1), p_tss->T_VID(it_e,2), p_tss->T_VID(it_e,3) };
        Mat3x3 B( mal::GMat3x3_From_Columns(actual_sdof[vec_nid[1]]-actual_sdof[vec_nid[0]],
                                            actual_sdof[vec_nid[2]]-actual_sdof[vec_nid[0]],
                                            actual_sdof[vec_nid[3]]-actual_sdof[vec_nid[0]] ) );
        Mat3x3 B0( mal::GMat3x3_From_Columns(default_sdof[vec_nid[1]]-default_sdof[vec_nid[0]],
                                             default_sdof[vec_nid[2]]-default_sdof[vec_nid[0]],
                                             default_sdof[vec_nid[3]]-default_sdof[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
        Mat3x3 B_InvB0( B * mal::Inverse(B0) );
        /*\note Direct but slower code
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
            vec_v[ed.m_FirstVID + it_vie] = tr * ( B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]])
                                                    + actual_sdof[vec_nid[0]] );
        */
        // Optimized version
        Transform3 tr_B_InvB0( tr.m_Pos, tr.m_Rot * B_InvB0 );
        Vec3 p0_0( tr.m_Rot * actual_sdof[vec_nid[0]] );
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
            vec_v[ed.m_FirstVID + it_vie] = tr_B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]]) + p0_0;

//#define __ENABLE_DCR3_DRAW_BV
#ifdef __ENABLE_DCR3_DRAW_BV
        // Draw BDOP
        if( ddf.Test( eODDF_BV )
            && !ed.m_BDOP[0].IsEmpty()
            && !ed.m_BDOP[1].IsEmpty()
            && !ed.m_BDOP[2].IsEmpty()
            && !ed.m_BDOP[3].IsEmpty() )
        {
            Vec3 bslab_pos[6];
            // compute global node pos
            Vec3 node_pos[4] = { tr*actual_sdof[vec_nid[0]], tr*actual_sdof[vec_nid[1]], tr*actual_sdof[vec_nid[2]], tr*actual_sdof[vec_nid[3]] };
            int i( ed.m_BDOP_BestSlabIdx );
            bslab_pos[0] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+1)%4];
            bslab_pos[1] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+2)%4];
            bslab_pos[2] = ed.m_BDOP[i].Min() * node_pos[i] + (1-ed.m_BDOP[i].Min()) * node_pos[(i+3)%4];
            bslab_pos[3] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+1)%4];
            bslab_pos[4] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+2)%4];
            bslab_pos[5] = ed.m_BDOP[i].Max() * node_pos[i] + (1-ed.m_BDOP[i].Max()) * node_pos[(i+3)%4];
            VIZ_TRIANGLE3( vs, bslab_pos[0], bslab_pos[1], bslab_pos[2], 1, Vec4(1,0.5,0.5,1), util::eVizStyle_Solid );
            VIZ_TRIANGLE3( vs, bslab_pos[3], bslab_pos[4], bslab_pos[5], 1, Vec4(0.5,1,0.5,1), util::eVizStyle_Solid );
        }
#endif
    }

    // Per-patch triangles
    mal::SetRandomSeed(666);
    for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
    {
        Vec4f color( mal::RandomV( Vec4f(0,0,0,0.25f), Vec4f(1,1,1,0.25f) ) );
        const geo::DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_patch] );
        for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
        {
            unsigned int tid( pd.m_FirstTID + it_tid );
            const DCR_TetSolidShape3::Triangle& tri( p_dcr->m_vecT[tid] );
            VIZ_TRIANGLE3( vs,
                           vec_v[tri.GetVID(0)],
                           vec_v[tri.GetVID(1)],
                           vec_v[tri.GetVID(2)],
                           1, color, util::eVizStyle_Wire );
#ifdef __ENABLE_DCR3_DRAW_ERRORS
            Vec3 barycenter( Real(0.333333)*(vec_v[tri.GetVID(0)]+
                                             vec_v[tri.GetVID(1)]+
                                             vec_v[tri.GetVID(2)]) );
            // Topology
            for( unsigned int it_ntit=0; it_ntit < 3; it_ntit++ )
            {
                unsigned int ntid( tri.GetNTID(it_ntit) );
                if( ntid != DCR_TetSolidShape3::Triangle::cInvalidTID )
                {
                    // const DCR_TetSolidShape3::Triangle& ntri( p_dcr->m_vecT[ntid] );
                    // Vec3 nbarycenter( Real(0.333333)*(vec_v[ntri.GetVID(0)]+
                    //                                   vec_v[ntri.GetVID(1)]+
                    //                                   vec_v[ntri.GetVID(2)]) );
                    // VIZ_SEGMENT3( vs, barycenter, nbarycenter, 1.0f, Vec4f(1,1,0,1) );//color );
                }
                else //draw open edges
                {
                    VIZ_SEGMENT3( vs,
                                  vec_v[tri.GetVID(it_ntit)], vec_v[tri.GetVID((it_ntit+1)%3)],
                                  5.0f, Vec4f(1,1,1,1 ) );
                }
            }
            // Sliver triangles
            {
                Vec3 d01( vec_v[tri.GetVID(1)]-vec_v[tri.GetVID(0)] );
                Vec3 d02( vec_v[tri.GetVID(2)]-vec_v[tri.GetVID(0)] );
                Vec3 d12( vec_v[tri.GetVID(2)]-vec_v[tri.GetVID(1)] );
                Real a( mal::Norm(d01) );
                Real b( mal::Norm(d02) );
                Real c( mal::Norm(d12) );
                Real circumradius( a*b*c / mal::Sqrt( (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c) ) ); //\see http://mathworld.wolfram.com/Circumradius.html
                if( circumradius / mal::Min(a,mal::Min(b,c)) > g_TTS_SliverThreshold ) //larger ratio => worse sliver (min value = 0.5)
                {
                    VIZ_TRIANGLE3( vs, vec_v[tri.GetVID(0)], vec_v[tri.GetVID(1)], vec_v[tri.GetVID(2)], 4.0f, Vec4f(1,0,1,1), util::eVizStyle_Wire );
                    // for( unsigned int it_vit=0; it_vit < 3; it_vit++ )
                    // {
                    //     VIZ_SEGMENT3( vs,
                    //                   Vec3(0,0,0), vec_v[tri.GetVID(it_vit)],
                    //                   5.0f, Vec4f(1,0,0,1) );
                    // }
                }
            }
#endif
        }
    }

#ifdef __ENABLE_DCR3_DRAW_ERRORS
    for( unsigned int it_e=0; it_e<p_dcr->m_NumElements; it_e++ )
    {
        const DCR_TetSolidShape3::Element& ed( p_dcr->m_vecE[it_e] );
        for( unsigned int it_p=ed.m_FirstPID; it_p<ed.m_FirstPID+ed.m_NumPatches; it_p++ )
        {
            const DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_p] );
            /*\todo ENABLE when NPID become available
            if( pd.GetNPID(0) < 0 || pd.GetNPID(0) >= p_dcr->m_NumPatches
                || pd.GetNPID(1) < 0 || pd.GetNPID(1) >= p_dcr->m_NumPatches
                || pd.GetNPID(2) < 0 || pd.GetNPID(2) >= p_dcr->m_NumPatches )
                GEO_LOG_ERROR( "P_DCR->E[%u].P[%u] = P_DCR->P[%u] has NPID (%u,%u,%u,%u) outside [%u...%u)",
                               it_e, it_p-ed.m_FirstPID, it_p,
                               pd.GetNPID(0), pd.GetNPID(1), pd.GetNPID(2), pd.GetNPID(3),
                               0, p_dcr->m_NumPatches );
            */
            for( unsigned int it_t=pd.m_FirstTID; it_t<pd.m_FirstTID+pd.m_NumTriangles; it_t++ )
            {
                const DCR_TetSolidShape3::Triangle& tri( p_dcr->m_vecT[it_t] );
                if( tri.GetVID(0) < ed.m_FirstVID || tri.GetVID(0) >= ed.m_FirstVID+ed.m_NumVertices
                    || tri.GetVID(1) < ed.m_FirstVID || tri.GetVID(1) >= ed.m_FirstVID+ed.m_NumVertices
                    || tri.GetVID(2) < ed.m_FirstVID || tri.GetVID(2) >= ed.m_FirstVID+ed.m_NumVertices )
                {
                    Vec3 p0( vec_v[tri.GetVID(0)] );
                    Vec3 p1( vec_v[tri.GetVID(1)] );
                    Vec3 p2( vec_v[tri.GetVID(2)] );
                    Vec3 bary( Real(0.333333) * (p0+p1+p2) );
                    Vec3 n( mal::SafeNormalized( mal::Cross( p1 - p0, p2 - p0 ) ) );
                    VIZ_SEGMENT3( vs,
                                  bary,
                                  bary + 0.1*n,
                                  3.0f, Vec4f(0.5,0,0,1) );
                    VIZ_TRIANGLE3( vs,
                                   p0,// + n*0.001,
                                   p1,// + n*0.001,
                                   p2,// + n*0.001,
                                   1, Vec4f(1,0,0,1), util::eVizStyle_Solid );
                }
                if( tri.GetNTID(0) < 0 || tri.GetNTID(0) >= p_dcr->m_NumTriangles
                    || tri.GetNTID(1) < 0 || tri.GetNTID(1) >= p_dcr->m_NumTriangles
                    || tri.GetNTID(2) < 0 || tri.GetNTID(2) >= p_dcr->m_NumTriangles )
                {
                    VIZ_TRIANGLE3( vs,
                                   vec_v[tri.GetVID(0)],
                                   vec_v[tri.GetVID(1)],
                                   vec_v[tri.GetVID(2)],
                                   1, Vec4f(0,1,0,1), util::eVizStyle_Solid );
                }
            }
        }
    }
#endif
}

} // namespace geo
