#ifndef GEO_SHAPE_SHAPE_FACTORY_H
#define GEO_SHAPE_SHAPE_FACTORY_H

#include <Geo/Config.h>
#include <Geo/shape/IShape.h>
#include <Geo/shape/ShapeLibrary.h>
#include <map>

namespace geo {

/*! Creates instances of shape definitions. Allows exclusive or shared
  instance creation.
*/
class ShapeFactory
{
public:
    ShapeFactory( const ShapeLibrary &shape_lib ) : m_ShapeLib(shape_lib) {}
    ~ShapeFactory();

    IShape *CreateExclusive( const ShapeDef &shape_def ); //!< Instantiate new shape and return
    IShape *CreateExclusive( ShapeID sid );      //!< Instantiate new shape and return
    
    IShape *CreateShared( ShapeID sid );  //!< Instantiate if non-existent and return
    
    const ShapeLibrary &GetLibrary() const { return m_ShapeLib; }
    
private:
    IShape *AllocShape( const ShapeDef &shape_def );
    
private:
    const ShapeLibrary &m_ShapeLib;
    typedef std::map<ShapeID,IShape*> MapId2Shape;
    MapId2Shape m_mapSharedShapes;
};

//! Free helper function to polymorphically initialize IShape params from a ShapeDef
bool InitShape( IShape *p_shape, const ShapeDef &shape_def, bool b_shared );

} //namespace geo

#endif // GEO_SHAPE_SHAPE_FACTORY_H
