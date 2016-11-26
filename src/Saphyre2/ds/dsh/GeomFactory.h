#ifndef S2_DS_DSH_GEOM_FACTORY_H
#define S2_DS_DSH_GEOM_FACTORY_H

#include "IGeom.h"
#include <Geo/IObject.h>
#include <Geo/ObjectFactory.h>

namespace S2 {
namespace ds {

//! Geom factory to simplify creation in MoP
class GeomFactory
{
public:
    GeomFactory()
    : m_ShapeFactory(m_ShapeLib)
    , m_GeoObjFactory(&m_ShapeFactory)
    {
        m_ShapeLib.Reserve( 1<<15 );
    }
    ~GeomFactory() {}

    IGeom *Create( machine_uint_type uid, IDynamicSystemHierarchy *parent, const geo::ShapeDef &shape_def );
    IGeom *Create( machine_uint_type uid, IDynamicSystemHierarchy *parent, geo::ShapeID shape_id );
    
    geo::ShapeLibrary &GetShapeLibrary() { return m_ShapeLib; }
    const geo::ShapeLibrary &GetShapeLibrary() const { return m_ShapeLib; }
    geo::ObjectFactory &GetGeoObjFactory() { return m_GeoObjFactory; }
    
private:
    geo::ShapeLibrary m_ShapeLib;
    geo::ShapeFactory m_ShapeFactory;
    geo::ObjectFactory m_GeoObjFactory;
};
    
}} // namespace S2::ds

#endif // S2_DS_DSH_GEOM_FACTORY_H
