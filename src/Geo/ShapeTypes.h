#ifndef GEO_SHAPE_TYPES_H
#define GEO_SHAPE_TYPES_H

#include "Config.h"

namespace geo {

//! ShapeID
typedef int32 ShapeID; //\note MUST be the same as util::ItemStream Item id, signed int to allow named shapes

//! BS constants
enum EConstants {
    cInvalidShapeId = 0
};

//! All Shape Types
enum EShapeType {
    eShape_External = 0,

    //---- 2D
    //-- Primitives
    eShape_Sphere2,
    eShape_Capsule2,
    eShape_Box2, //Superseeded by ConvexPolygon2
    eShape_Triangle2, //Superseeded by ConvexPolygon2
    //eShape_Quadric2, //Implicit: Circle,Ellipse (+infinite Parabola,Cone,Hyperbole... unless tappered)
    eShape_Segment2, //superseeded by Path2
    //eShape_Spline2, //superseeded by Path2
    //-- Complexes
    //eShape_ConvexPolygon2, //\todo superseedes Triangle2 and Box2, but more expensive...
    eShape_Polygonal2, //\todo deprecated/superseeded by Path2
    eShape_Path2,
    eShape_MeshSolid2,
    //\todo eShape_TriSolid2: specialization of MeshSolid2 with faster ops based in triangular faces
    eShape_SphereSet2,
    //-- Infinite
    eShape_Plane2,

    //---- 3D
    //-- Primitives
    eShape_Sphere3,
    eShape_Capsule3,
    eShape_Cylinder3,
    eShape_Box3,
    eShape_Triangle3,
    eShape_Tetrahedron3,
    //eShape_Quadric3, //Implicit: Sphere, Ellipsoid (+infinite Paraboloid Cone,Hyperboloid... unless tappered)
    eShape_Segment3, //superseeded by Path3
    //eShape_Spline3, //superseeded by Path3
    //-- Complexes
    //eShape_ConvexPolygon3, //\todo superseedes Triangle3 but more expensive...
    //eShape_ConvexPolyhedron3, //\todo superseedes Tetrahedron3 and Box3, but potentially inefficient
    eShape_Polygonal3, //\todo deprecated/superseeded by Path3
    eShape_Path3,
    eShape_TriSurface3, //!< \note Better AVOID MeshSurface3 with N-sided faces, too complex and not directly renderable
    eShape_TetSolid3, //!< \note Better AVOID MeshSolid3 with N-faced solid elements, non tetra volume elements may be HARD to do right...
    eShape_SphereSet3,
    //-- Infinite
    eShape_Plane3,

    cNumShapeTypes
};


//! \name Fwd-declare all shape classes
//@{
//---- Generic D
template<unsigned D> class GSphereShape; typedef GSphereShape<2> SphereShape2; typedef GSphereShape<3> SphereShape3;
template<unsigned D> class GCapsuleShape; typedef GCapsuleShape<2> CapsuleShape2; typedef GCapsuleShape<3> CapsuleShape3;
template<unsigned D> class GBoxShape; typedef GBoxShape<2> BoxShape2; typedef GBoxShape<3> BoxShape3;
template<unsigned D> class GSphereSetShape; typedef GSphereSetShape<2> SphereSetShape2; typedef GSphereSetShape<3> SphereSetShape3;
template<unsigned D> class GPlaneShape; typedef GPlaneShape<2> PlaneShape2; typedef GPlaneShape<3> PlaneShape3;

//---- Pure 2D
class SegmentShape2;
class PolygonalShape2;
class PathShape2;
//class SplineShape2;

class TriangleShape2;
class MeshSolidShape2;
//class TriSolidShape2;

//---- Pure 3D
class SegmentShape3;
class PolygonalShape3;
//class SplineShape3;
//class SphereShape3;
class CylinderShape3;

class TriangleShape3;
class TetrahedronShape3;
class TriSurfaceShape3;
class TetSolidShape3;
//@}

//! \name Editable shapes
//@{
class EditablePolygonalShape2;
class EditablePolygonalShape3;
//class EditableSplineShape2;
//class EditableSplineShape3;
class EditablePathShape2;
//....
//@}

//! \name ClassToEnum Type-Trait for all shape classes
//@{
template <class ShapeT> inline EShapeType shape_type_of() { return eShape_External; }
//-- Primitives
template<> inline EShapeType shape_type_of<SphereShape2>() { return eShape_Sphere2; }
template<> inline EShapeType shape_type_of<SphereShape3>() { return eShape_Sphere3; }
template<> inline EShapeType shape_type_of<CapsuleShape2>() { return eShape_Capsule2; }
template<> inline EShapeType shape_type_of<CapsuleShape3>() { return eShape_Capsule3; }
template<> inline EShapeType shape_type_of<BoxShape2>() { return eShape_Box2; }
template<> inline EShapeType shape_type_of<BoxShape3>() { return eShape_Box3; }
//-- Complexes
template<> inline EShapeType shape_type_of<SphereSetShape2>() { return eShape_SphereSet2; }
template<> inline EShapeType shape_type_of<SphereSetShape3>() { return eShape_SphereSet3; }
//-- Infinite
template<> inline EShapeType shape_type_of<PlaneShape2>() { return eShape_Plane2; }
template<> inline EShapeType shape_type_of<PlaneShape3>() { return eShape_Plane3; }

//---- 2D
template<> inline EShapeType shape_type_of<SegmentShape2>() { return eShape_Segment2; }
template<> inline EShapeType shape_type_of<PolygonalShape2>() { return eShape_Polygonal2; }
//template<> inline EShapeType shape_type_of<SplineShape2>() { return eShape_Spline2; }

template<> inline EShapeType shape_type_of<TriangleShape2>() { return eShape_Triangle2; }
template<> inline EShapeType shape_type_of<MeshSolidShape2>() { return eShape_MeshSolid2; }
//template<> inline EShapeType shape_type_of<TriSolidShape2>() { return eShape_TriSolid2; }
//---- 3D
template<> inline EShapeType shape_type_of<SegmentShape3>() { return eShape_Segment3; }
template<> inline EShapeType shape_type_of<PolygonalShape3>() { return eShape_Polygonal3; }
//template<> inline EShapeType shape_type_of<SplineShape3>() { return eShape_Spline3; }
template<> inline EShapeType shape_type_of<CylinderShape3>() { return eShape_Cylinder3; }

template<> inline EShapeType shape_type_of<TriangleShape3>() { return eShape_Triangle3; }
template<> inline EShapeType shape_type_of<TetrahedronShape3>() { return eShape_Tetrahedron3; }
template<> inline EShapeType shape_type_of<TriSurfaceShape3>() { return eShape_TriSurface3; }
template<> inline EShapeType shape_type_of<TetSolidShape3>() { return eShape_TetSolid3; }
//@}

}

#endif //GEO_SHAPE_TYPES_H
