#ifndef GEO_OBJECT_FACTORY_H
#define GEO_OBJECT_FACTORY_H

#include <Geo/Config.h>
#include <Geo/IObject.h>
#include <Geo/shape/IShape.h>
#include <Geo/shape/ShapeFactory.h>

namespace geo {

/*! Creates object instances for given shape definitions.

  - ShapeID instantiation: Requires a valid ShapeFactory to be specified on creation.
  - IShape and shape_def instantiation: Does NOT NEED a
    ShapeFactory. Can be done with STATIC methods that do not store
    any local copies of IObjects.
  \todo Created IObjects should be stored internally as ShapeFactory does

  Objects with exclusive simple shapes try to allocate them inside the
  object (GObjectES<>) itself, while objects with shared shapes use
  GObjectSS<>
*/
class ObjectFactory
{
public:
    ObjectFactory( ShapeFactory* p_shape_factory = 0 ) : m_pShapeFactory(p_shape_factory) {}
    ~ObjectFactory() {}

    static IObject* STATIC_CreateES( const ShapeDef& shape_def ); //!< Create GObjectES with embedded shape TEMPORAL

    IObject* CreateES( const ShapeDef& shape_def ); //!< Create GObjectES with embedded shape
    IObject* CreateSS( const IShape* p_shape );     //!< Create GObjectSS with given shared shape

    IObject* CreateES( ShapeID sid );   //!< Create GObjectES with exclusive shape from a ShapeLibrary
    IObject* CreateSS( ShapeID sid );   //!< Create GObjectSS with shared shape from a ShapeFactory

    void Release( IObject* p_object );

private:
    ShapeFactory* m_pShapeFactory;
    //\todo could store all created objects here, use a custom allocator, etc... the factory owns the IObjects.
};

} //namespace geo

#endif // GEO_OBJECT_FACTORY_H
