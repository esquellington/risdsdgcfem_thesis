#include "ObjectFactory.h"
#include <Geo/shape/shape.h> //shape types

namespace geo {

IObject* ObjectFactory::STATIC_CreateES( const ShapeDef& shape_def )
{
    if( eType_Property_Group == shape_def.GetType() ) { GEO_LOG_WARNING( "shape_def is of eType_Property_Group, assuming correct... (loaded from TXT?)" ); }
    else GEO_ASSERT( eType_ShapeDef == shape_def.GetType() );
    util::ItemStream::ItemIt params_it( shape_def.GetSubItem() );
    // Alloc
    IObject *pObject(0);
    switch( (EShapeType)params_it.Find("shape_type").Get<uint32>() )
    {
    case eShape_Plane2: pObject = new GObjectES<PlaneShape2>(); break;
    case eShape_Plane3: pObject = new GObjectES<PlaneShape3>(); break;
    case eShape_Sphere2: pObject = new GObjectES<SphereShape2>(); break;
    case eShape_Sphere3: pObject = new GObjectES<SphereShape3>(); break;
    case eShape_Capsule2: pObject = new GObjectES<CapsuleShape2>(); break;
    case eShape_Capsule3: pObject = new GObjectES<CapsuleShape3>(); break;
    case eShape_Box2: pObject = new GObjectES<BoxShape2>(); break;
    case eShape_Box3: pObject = new GObjectES<BoxShape3>(); break;
    case eShape_SphereSet2: GEO_ASSERT(false); return 0; break; //\todo
    case eShape_SphereSet3: GEO_ASSERT(false); return 0; break; //\todo
    case eShape_Polygonal2: pObject = new GObjectES<PolygonalShape2>(); break;
    case eShape_Polygonal3: pObject = new GObjectES<PolygonalShape3>(); break;
    case eShape_Path2: pObject = new GObjectES<PathShape2>(); break;
    case eShape_MeshSolid2: pObject = new GObjectES<MeshSolidShape2>(); break;
    case eShape_TetSolid3: pObject = new GObjectES<TetSolidShape3>(); break;
    case eShape_TriSurface3: pObject = new GObjectES<TriSurfaceShape3>(); break;
    default: GEO_ASSERT(false); return 0; break;
    }
    // Init as non-shared
    InitShape( pObject->GetShapeInterface(), shape_def, false );
    pObject->ResetDOF();
    return pObject;
}

IObject* ObjectFactory::CreateES( const ShapeDef& shape_def )
{
    IObject* pObject = STATIC_CreateES(shape_def);
    //\todo add pObject to managed object map or whatever...
    return pObject;
}

IObject* ObjectFactory::CreateSS( const IShape* p_shape )
{
    IObject* pObject(0);
    switch( p_shape->GetType() )
    {
    case eShape_Plane2:
        {
            GObjectSS<PlaneShape2> *pOSS( new GObjectSS<PlaneShape2>() );
            pOSS->SetShape( static_cast<const PlaneShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Plane3:
        {
            GObjectSS<PlaneShape3> *pOSS( new GObjectSS<PlaneShape3>() );
            pOSS->SetShape( static_cast<const PlaneShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Sphere2:
        {
            GObjectSS<SphereShape2> *pOSS( new GObjectSS<SphereShape2>() );
            pOSS->SetShape( static_cast<const SphereShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Sphere3:
        {
            GObjectSS<SphereShape3> *pOSS( new GObjectSS<SphereShape3>() );
            pOSS->SetShape( static_cast<const SphereShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Capsule2:
        {
            GObjectSS<CapsuleShape2> *pOSS( new GObjectSS<CapsuleShape2>() );
            pOSS->SetShape( static_cast<const CapsuleShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Capsule3:
        {
            GObjectSS<CapsuleShape3> *pOSS( new GObjectSS<CapsuleShape3>() );
            pOSS->SetShape( static_cast<const CapsuleShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Box2:
        {
            GObjectSS<BoxShape2> *pOSS( new GObjectSS<BoxShape2>() );
            pOSS->SetShape( static_cast<const BoxShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Box3:
        {
            GObjectSS<BoxShape3> *pOSS( new GObjectSS<BoxShape3>() );
            pOSS->SetShape( static_cast<const BoxShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
//\todo GSphereSetShape
    case eShape_Polygonal2:
        {
            GObjectSS<PolygonalShape2> *pOSS( new GObjectSS<PolygonalShape2>() );
            pOSS->SetShape( static_cast<const PolygonalShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Polygonal3:
        {
            GObjectSS<PolygonalShape3> *pOSS( new GObjectSS<PolygonalShape3>() );
            pOSS->SetShape( static_cast<const PolygonalShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_Path2:
        {
            GObjectSS<PathShape2> *pOSS( new GObjectSS<PathShape2>() );
            pOSS->SetShape( static_cast<const PathShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_MeshSolid2:
        {
            GObjectSS<MeshSolidShape2> *pOSS( new GObjectSS<MeshSolidShape2>() );
            pOSS->SetShape( static_cast<const MeshSolidShape2*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_TetSolid3:
        {
            GObjectSS<TetSolidShape3> *pOSS( new GObjectSS<TetSolidShape3>() );
            pOSS->SetShape( static_cast<const TetSolidShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
    case eShape_TriSurface3:
        {
            GObjectSS<TriSurfaceShape3> *pOSS( new GObjectSS<TriSurfaceShape3>() );
            pOSS->SetShape( static_cast<const TriSurfaceShape3*>(p_shape) );
            pObject = pOSS;
        }
        break;
    default: GEO_LOG_ERROR("ObjectFactory Cannot create shape type %d", (int)p_shape->GetType() ); GEO_ASSERT(false); return 0; break;
    }
    pObject->ResetDOF();
    return pObject;
}

IObject* ObjectFactory::CreateES( ShapeID sid )
{
    ShapeDef shape_def( m_pShapeFactory->GetLibrary().Lookup(sid) );
    GEO_ASSERT( shape_def.IsValid() );
    return CreateES(shape_def);
}

IObject* ObjectFactory::CreateSS( ShapeID sid )
{
    IShape* pShape = m_pShapeFactory->CreateShared(sid);
    GEO_ASSERT(0 != pShape);
    return CreateSS(pShape);
}

void ObjectFactory::Release( IObject* p_object )
{
    //\todo Remove from factory, when tracked
    delete p_object;
}

} //namespace geo
