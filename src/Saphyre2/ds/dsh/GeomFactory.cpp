#include "GeomFactory.h"
#include <Saphyre2/ds/DSG.h>

namespace S2 {
namespace ds {

//---- GeomFactory implementation
IGeom *GeomFactory::Create( machine_uint_type uid, IDynamicSystemHierarchy *parent, const geo::ShapeDef &shape_def )
{
    geo::IObject* pGeoObj = m_GeoObjFactory.CreateES(shape_def);
    if( 2 == pGeoObj->GetDimension() ) return new Geom2(uid,parent,static_cast< geo::IObject2* >(pGeoObj));
    else if( 3 == pGeoObj->GetDimension() ) return new Geom3(uid,parent,static_cast< geo::IObject3* >(pGeoObj));
    else return 0;    
}

IGeom *GeomFactory::Create( machine_uint_type uid, IDynamicSystemHierarchy *parent, geo::ShapeID shape_id )
{
    geo::IObject *pGeoObj = m_GeoObjFactory.CreateSS(shape_id);
    if( 2 == pGeoObj->GetDimension() ) return new Geom2(uid,parent,static_cast< geo::IObject2 * >(pGeoObj));
    else if( 3 == pGeoObj->GetDimension() ) return new Geom3(uid,parent,static_cast< geo::IObject3 * >(pGeoObj));
    else return 0;
}

}} // namespace S2::ds
