#ifndef GEO_SHAPE_SHAPE_LIBRARY_H
#define GEO_SHAPE_SHAPE_LIBRARY_H

#include <Geo/Config.h>
#include <Geo/ShapeTypes.h>
#include <Geo/shape/IShape.h>
#include <util/ItemStream.h>

enum EPlaGeoTypes { eType_ShapeDef = cNumPlaTypes+1 };

namespace geo {

typedef util::ItemStream ShapeDefStream;
typedef util::ItemStream::ItemIt ShapeDef;

/*! Holds a collection of shape definitions indexed by ShapeID, and
  offers methods to Register new ones.

  ShapeDef may be stored in an externally or internally (using
  Reserve()) ShapeDefStream.

  IMPORTANT: No realloc supported by now, which because it would
  invalidate all pointers to the SDS.
*/
class ShapeLibrary
{
public:
    ShapeLibrary();
    // ShapeLibrary( ShapeDefStream *p_sds ); \todo Disabled to avoid realloc conundrums
    ~ShapeLibrary();

    bool Reserve( uint32 sds_size, Flags32 realloc_flags = util::ItemStream::eRealloc_Default );

    //! \name Generic shape definition
    //@{
    ShapeID Register( const IShape &shape, const char* name = 0 );
    ShapeID Register( const IShape *p_shape, const char* name = 0 );
    ShapeID Register( const ShapeDef &shape_def, const char* name = 0 );
    //@}

    //!\name Simple shape definition
    //@{
    ShapeID Register_Plane2( const Vec2 &normal, Real coeff_d, bool b_half_space, const char* name = 0 );
    ShapeID Register_Plane3( const Vec3 &normal, Real coeff_d, bool b_half_space, const char* name = 0 );
    ShapeID Register_Sphere2( Real radius, const char* name = 0 );
    ShapeID Register_Sphere3( Real radius, const char* name = 0 );
    ShapeID Register_Capsule2( Real radius, Real half_height, const char* name = 0 );
    ShapeID Register_Capsule3( Real radius, Real half_height, const char* name = 0 );
    ShapeID Register_Box2( const Vec2 &half_sizes, const char* name = 0 );
    ShapeID Register_Box3( const Vec3 &half_sizes, const char* name = 0 );
    // \todo SphereSetN (var alloc!)
    ShapeID Register_Polygonal2( const PolygonalShape2* p_ps2, const char* name = 0 );
    ShapeID Register_Polygonal3( const PolygonalShape3* p_ps3, const char* name = 0 );
    ShapeID Register_PathShape2( const PathShape2 *p_ps2, const char* name = 0 );
    ShapeID Register_MeshSolidShape2( const MeshSolidShape2 *p_mss2, const char* name = 0 );
    ShapeID Register_TetSolidShape3( const TetSolidShape3 *p_tss3, const char* name = 0 );
    ShapeID Register_TriSurfaceShape3( const TriSurfaceShape3 *p_tss3, const char* name = 0 );
    //...
    //@}

    ShapeDef Lookup( ShapeID sid ) const;
    ShapeID GetShapeIdByName( const char *name ) const;
    ShapeDef GetFirstShapeDef() const;
    void Clear();
    bool IsEmpty() const { return 0 == m_pSDS || m_pSDS->IsEmpty(); }

    //!\name Serialization TEMPORAL, may change API or even do it externally...
    //@{
    bool Save( const char *file_name, bool b_binary = false ) const;
    bool Load( const char *file_name, bool b_binary = false );
    //@}

private:
    ShapeDefStream *m_pSDS;
    bool m_bIsExternalSDS;
    ShapeID m_LastShapeId;
};

} //namespace geo

#endif // GEO_SHAPE_SHAPE_LIBRARY_H
