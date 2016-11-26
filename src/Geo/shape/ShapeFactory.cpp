#include "ShapeFactory.h"
//-- shape types
#include "shape.h"

namespace geo {

ShapeFactory::~ShapeFactory()
{
    for( MapId2Shape::const_iterator it_shape = m_mapSharedShapes.begin();
         it_shape != m_mapSharedShapes.end();
         ++it_shape )
        delete it_shape->second;
}

IShape *ShapeFactory::AllocShape( const ShapeDef &shape_def )
{
    if( eType_Property_Group == shape_def.GetType() ) { GEO_LOG_WARNING( "shape_def is of eType_Property_Group, assuming correct... (loaded from TXT?)" ); }
    else GEO_ASSERT( eType_ShapeDef == shape_def.GetType() );
    util::ItemStream::ItemIt params_it( shape_def.GetSubItem() );
    switch( (EShapeType)params_it.Find("shape_type").Get<uint32>() )
    {
    case eShape_Plane2: return new PlaneShape2(); break;
    case eShape_Plane3: return new PlaneShape3(); break;
    case eShape_Sphere2: return new SphereShape2(); break;
    case eShape_Sphere3: return new SphereShape3(); break;
    case eShape_Capsule2: return new CapsuleShape2(); break;
    case eShape_Capsule3: return new CapsuleShape3(); break;
    case eShape_Box2: return new BoxShape2(); break;
    case eShape_Box3: return new BoxShape3(); break;
    case eShape_SphereSet2: return new SphereSetShape2(); break;
    case eShape_SphereSet3: return new SphereSetShape3(); break;
    case eShape_Polygonal2: return new PolygonalShape2(); break;
    case eShape_Polygonal3: return new PolygonalShape3(); break;
    case eShape_Path2: return new PathShape2(); break;
    case eShape_MeshSolid2: return new MeshSolidShape2(); break;
    case eShape_TetSolid3: return new TetSolidShape3(); break;
    case eShape_TriSurface3: return new TriSurfaceShape3(); break;
    default: GEO_ASSERT(false); return 0; break;
    }
    return 0;
}

IShape *ShapeFactory::CreateExclusive( const ShapeDef &shape_def )
{
    if( eType_Property_Group == shape_def.GetType() ) { GEO_LOG_WARNING( "shape_def is of eType_Property_Group, assuming correct... (loaded from TXT?)" ); }
    else GEO_ASSERT( eType_ShapeDef == shape_def.GetType() );
    // Alloc & Init as non-shared
    IShape *pShape( AllocShape( shape_def ) );
    InitShape( pShape, shape_def, false );
    return pShape;
}

IShape *ShapeFactory::CreateExclusive( ShapeID sid )
{
    ShapeDef shape_def( m_ShapeLib.Lookup(sid) );
    GEO_ASSERT( shape_def.IsValid() );
    return CreateExclusive(shape_def);
}

IShape *ShapeFactory::CreateShared( ShapeID sid )
{
    //find or create
    MapId2Shape::const_iterator it_shape = m_mapSharedShapes.find(sid);
    if( it_shape == m_mapSharedShapes.end() )
    {
        ShapeDef shape_def( m_ShapeLib.Lookup(sid) );
        GEO_ASSERT( shape_def.IsValid() );
        // Alloc & Init as shared
        IShape *pShape( AllocShape( shape_def ) );
        InitShape( pShape, shape_def, true );
        m_mapSharedShapes[sid] = pShape;
        return pShape;
    }
    else
        return it_shape->second;
}

bool InitShape( IShape *p_shape, const ShapeDef &shape_def, bool b_shared )
{
    if( eType_Property_Group == shape_def.GetType() ) { GEO_LOG_WARNING( "shape_def is of eType_Property_Group, assuming correct... (loaded from TXT?)" ); }
    else GEO_ASSERT( eType_ShapeDef == shape_def.GetType() );
    util::ItemStream::ItemIt params_it( shape_def.GetSubItem() );
    GEO_ASSERT( (uint32)p_shape->GetType() == params_it.Find("shape_type").Get<uint32>() );
    switch( p_shape->GetType() )
    {
    case eShape_Plane2:
        reinterpret_cast<PlaneShape2*>(p_shape)->Init( params_it.Find("normal").Get<Vec2>(),
                                                       params_it.Find("coeff_d").Get<Real>(),
                                                       params_it.Find("is_half_space").Get<bool>() );
        break;
    case eShape_Plane3:
        reinterpret_cast<PlaneShape3*>(p_shape)->Init( params_it.Find("normal").Get<Vec3>(),
                                                       params_it.Find("coeff_d").Get<Real>(),
                                                       params_it.Find("is_half_space").Get<bool>() );
        break;
    case eShape_Sphere2: reinterpret_cast<SphereShape2*>(p_shape)->SetRadius(params_it.Find("radius").Get<Real>()); break;
    case eShape_Sphere3: reinterpret_cast<SphereShape3*>(p_shape)->SetRadius(params_it.Find("radius").Get<Real>()); break;
    case eShape_Capsule2: reinterpret_cast<CapsuleShape2*>(p_shape)->SetRadiusAndHalfHeight( params_it.Find("radius").Get<Real>(),
                                                                                             params_it.Find("half_height").Get<Real>() ); break;
    case eShape_Capsule3: reinterpret_cast<CapsuleShape3*>(p_shape)->SetRadiusAndHalfHeight( params_it.Find("radius").Get<Real>(),
                                                                                             params_it.Find("half_height").Get<Real>() ); break;
    case eShape_Box2: reinterpret_cast<BoxShape2*>(p_shape)->SetHalfSizes(params_it.Find("half_sizes").Get<Vec2>()); break;
    case eShape_Box3: reinterpret_cast<BoxShape3*>(p_shape)->SetHalfSizes(params_it.Find("half_sizes").Get<Vec3>()); break;
    case eShape_SphereSet2:
        reinterpret_cast<SphereSetShape2*>(p_shape)->SetCount( params_it.Find("count").Get<uint32>() );
        reinterpret_cast<SphereSetShape2*>(p_shape)->SetRadius( params_it.Find("radius").Get<Real>() );
        break;
    case eShape_SphereSet3:
        reinterpret_cast<SphereSetShape3*>(p_shape)->SetCount( params_it.Find("count").Get<uint32>() );
        reinterpret_cast<SphereSetShape3*>(p_shape)->SetRadius( params_it.Find("radius").Get<Real>() );
        break;
    case eShape_Polygonal2:
        {
            PolygonalShape2 *pPS2( reinterpret_cast<PolygonalShape2*>(p_shape) );
            //\todo GEO_ASSERT( params_it.Find("version").Get<uint32>(0) == PolygonalShape2::cVersion ); //check matching version > 0
            pPS2->Init( params_it.Find("num_vertices").Get<uint32>(),
                        params_it.Find("is_closed").Get<bool>(),
                        params_it.Find("radius").Get<Real>(),
                        reinterpret_cast<const Vec2*>( params_it.Find("vec_points").GetArrayPtr<uint8>() ),
                        reinterpret_cast<const Vec2*>( params_it.Find("vec_normals").SafeGetArrayPtr<uint8>(0) ),
                        reinterpret_cast<const Vec2*>( params_it.Find("vec_tangents").SafeGetArrayPtr<uint8>(0) ),
                        reinterpret_cast<const Real*>( params_it.Find("vec_lambdas").SafeGetArrayPtr<uint8>(0) ),
                        reinterpret_cast<const Flags32*>( params_it.Find("vec_flags").SafeGetArrayPtr<uint8>(0) ) );
        }
        break;
    case eShape_Polygonal3:
        {
            PolygonalShape3 *pPS3( reinterpret_cast<PolygonalShape3*>(p_shape) );
            //\todo GEO_ASSERT( params_it.Find("version").Get<uint32>(0) == PolygonalShape3::cVersion ); //check matching version > 0
            pPS3->Init( params_it.Find("num_vertices").Get<uint32>(),
                        params_it.Find("is_closed").Get<bool>(),
                        params_it.Find("radius").Get<Real>(),
                        reinterpret_cast<const Vec3*>( params_it.Find("vec_points").GetArrayPtr<uint8>() ),
                        reinterpret_cast<const Vec3*>( params_it.Find("vec_normals").SafeGetArrayPtr<uint8>(0) ),
                        reinterpret_cast<const Vec3*>( params_it.Find("vec_tangents").SafeGetArrayPtr<uint8>(0) ),
                        reinterpret_cast<const Real*>( params_it.Find("vec_lambdas").SafeGetArrayPtr<uint8>(0) ),
                        reinterpret_cast<const Flags32*>( params_it.Find("vec_flags").SafeGetArrayPtr<uint8>(0) ) );
        }
        break;
    case eShape_Path2:
        {
            PathShape2 *pPS2( reinterpret_cast<PathShape2*>(p_shape) );
            //\todo GEO_ASSERT( params_it.Find("version").Get<uint32>(0) == PathShape2::cVersion ); //check matching version > 0
            pPS2->SetBakedData( b_shared,
                                params_it.Find("num_v").Get<uint32>(),
                                params_it.Find("num_e").Get<uint32>(),
                                params_it.Find("is_closed").Get<bool>(),
                                reinterpret_cast<const Vec2*>( params_it.Find("vec_points").GetArrayPtr<uint8>() ),
                                reinterpret_cast<const PathShape2::edge_type*>( params_it.Find("vec_e").GetArrayPtr<uint8>() ) );
        }
        break;
    case eShape_MeshSolid2:
        {
            MeshSolidShape2 *pMSS2( reinterpret_cast<MeshSolidShape2*>(p_shape) );
            GEO_ASSERT( params_it.Find("version").SafeGet<uint32>(0) == MeshSolidShape2::cVersion ); //check matching version > 0
            pMSS2->SetBakedData( b_shared,
                                 params_it.Find("num_v").Get<uint32>(),
                                 params_it.Find("num_p").Get<uint32>(),
                                 params_it.Find("num_he").Get<uint32>(),
                                 params_it.Find("num_boundary_p").Get<uint32>(),
                                 params_it.Find("num_boundary_he").Get<uint32>(),
                                 params_it.Find("num_l").Get<uint32>(),
                                 reinterpret_cast<const Vec2*>( params_it.Find("vec_points").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const MeshSolidShape2::vertex_type*>( params_it.Find("vec_v").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const MeshSolidShape2::polygon_type*>( params_it.Find("vec_p").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const MeshSolidShape2::half_edge_type*>( params_it.Find("vec_he").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const MeshSolidShape2::polygon_layer_type*>( params_it.Find("vec_l").GetArrayPtr<uint8>() ) );
            // Optional annotations
            util::ItemStream::ItemIt dcr_it( params_it.Find("dcr") );
            if( dcr_it.IsValid() )
            {
                dcr_it = dcr_it.GetSubItem();
                GEO_ASSERT( dcr_it.Find("version").SafeGet<uint32>(0) == DCR_MeshSolidShape2::cVersion ); //check matching version > 0
                DCR_MeshSolidShape2* pDCR = new DCR_MeshSolidShape2();
                pDCR->SetBakedData( b_shared,
                                    dcr_it.Find("num_e").Get<uint32>(),
                                    dcr_it.Find("num_p").Get<uint32>(),
                                    dcr_it.Find("num_s").Get<uint32>(),
                                    dcr_it.Find("num_v").Get<uint32>(),
                                    reinterpret_cast<const DCR_MeshSolidShape2::Element*>( dcr_it.Find("vec_e").GetArrayPtr<uint8>() ),
                                    reinterpret_cast<const DCR_MeshSolidShape2::Patch*>( dcr_it.Find("vec_p").GetArrayPtr<uint8>() ),
                                    reinterpret_cast<const DCR_MeshSolidShape2::Segment*>( dcr_it.Find("vec_s").GetArrayPtr<uint8>() ),
                                    reinterpret_cast<const Vec2*>( dcr_it.Find("vec_v").GetArrayPtr<uint8>() ) );
                pMSS2->SetBakedDCR_StrictlyNonshared_UglyHack( pDCR );
            }
        }
        break;
    case eShape_TetSolid3:
        {
            TetSolidShape3 *pTSS3( reinterpret_cast<TetSolidShape3*>(p_shape) );
            GEO_ASSERT( params_it.Find("version").SafeGet<uint32>(0) == TetSolidShape3::cVersion ); //check matching version > 0
            pTSS3->SetBakedData( b_shared,
                                 params_it.Find("num_v").Get<uint32>(),
                                 params_it.Find("num_t").Get<uint32>(),
                                 params_it.Find("num_bf").Get<uint32>(),
                                 params_it.Find("num_l").Get<uint32>(),
                                 reinterpret_cast<const Vec3*>( params_it.Find("vec_points").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const TetSolidShape3::vertex_type*>( params_it.Find("vec_v").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const TetSolidShape3::tetrahedron_type*>( params_it.Find("vec_t").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const TetSolidShape3::boundary_face_type*>( params_it.Find("vec_bf").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const TetSolidShape3::tetrahedron_layer_type*>( params_it.Find("vec_l").GetArrayPtr<uint8>() ) );
            // Optional annotations
            util::ItemStream::ItemIt dcr_it( params_it.Find("dcr") );
            if( dcr_it.IsValid() )
            {
                dcr_it = dcr_it.GetSubItem();
                GEO_ASSERT( dcr_it.Find("version").SafeGet<uint32>(0) == DCR_TetSolidShape3::cVersion ); //check matching version > 0
                DCR_TetSolidShape3* pDCR = new DCR_TetSolidShape3();
                pDCR->SetBakedData( b_shared,
                                    dcr_it.Find("num_e").Get<uint32>(),
                                    dcr_it.Find("num_p").Get<uint32>(),
                                    dcr_it.Find("num_t").Get<uint32>(),
                                    dcr_it.Find("num_v").Get<uint32>(),
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
                                    dcr_it.Find("num_npid").Get<uint32>(),
#endif
                                    reinterpret_cast<const DCR_TetSolidShape3::Element*>( dcr_it.Find("vec_e").GetArrayPtr<uint8>() ),
                                    reinterpret_cast<const DCR_TetSolidShape3::Patch*>( dcr_it.Find("vec_p").GetArrayPtr<uint8>() ),
                                    reinterpret_cast<const DCR_TetSolidShape3::Triangle*>( dcr_it.Find("vec_t").GetArrayPtr<uint8>() ),
                                    reinterpret_cast<const Vec3*>( dcr_it.Find("vec_v").GetArrayPtr<uint8>() )
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
                                    ,reinterpret_cast<const uint32*>( dcr_it.Find("vec_npid").GetArrayPtr<uint8>() )
#endif
                    );
                pTSS3->SetBakedDCR_StrictlyNonshared_UglyHack( pDCR );
            }
        }
        break;
    case eShape_TriSurface3:
        {
            TriSurfaceShape3 *pTSS3( reinterpret_cast<TriSurfaceShape3*>(p_shape) );
            GEO_ASSERT( params_it.Find("version").SafeGet<uint32>(0) == TriSurfaceShape3::cVersion ); //check matching version > 0
            pTSS3->SetBakedData( b_shared,
                                 params_it.Find("flags").Get<int32>(), //\todo Support GFlags<T> serialization!?
                                 params_it.Find("num_v").Get<uint32>(),
                                 params_it.Find("num_t").Get<uint32>(),
                                 reinterpret_cast<const Vec3*>( params_it.Find("vec_points").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const TriSurfaceShape3::vertex_type*>( params_it.Find("vec_v").GetArrayPtr<uint8>() ),
                                 reinterpret_cast<const TriSurfaceShape3::triangle_type*>( params_it.Find("vec_t").GetArrayPtr<uint8>() ) );
        }
        break;
    default: return false;
    }
    return true;
}

} //namespace geo
