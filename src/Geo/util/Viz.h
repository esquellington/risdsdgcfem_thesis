/*! \file Viz.h
  Visualization helpers for Geo entities.

  \note ALL viz functionality for Geo is here, completely decoupled
  from Geo entities implementation.
*/

#ifndef GEO_UTIL_VIZ_H
#define GEO_UTIL_VIZ_H

#include <Geo/Config.h>
#include <util/VizStream.h>
//#include <Geo/shape/IShape.h> //\todo ONLY required because IShape2/3 cannot be fwd declared... shitty

//TEMP: We need all this because we use bv::BVH_ST_DG_AABB2 here, and
//a fwd declaration is not easy because it's a typedef template...
// => Consider separate bv/Viz.h for BV/BVH and try again to use fwd declaration...
//@{
#include <Geo/bv/GSphere.h>
#include <Geo/bv/GAABB.h>
#include <Geo/bv/GKDOP.h>
#include <Geo/bv/GBoundingVolumeHierarchy.h>
//@}

//
#include <Geo/np/ContactData.h>

namespace geo {


/* Draw Debug Flags
   Disjoint flags for Object, Shape, BVH allow merging them into a single mask:
   - Object: bits 0..3
   - Shape: bits 4..16
   - BVH: bits 29..31
*/
enum EObjectDDF { eODDF_Nothing    = 0,
                  eODDF_Axis       = 1<<0,
                  eODDF_Shape      = 1<<1,
                  eODDF_Embedding  = 1<<2,
                  eODDF_BV         = 1<<3,
                  eODDF_Everything = 0xF0000000,
                  eODDF_Default    = eODDF_Shape };

enum EShapeDDF { eSDDF_Nothing    = 0,
                 // Aspects
                 eSDDF_Boundary   = 1<<4, //Basic shape
                 eSDDF_Interior   = 1<<5, //Internal structure, if any (Ex: MeshSolidShape2 elements)
                 eSDDF_Topology   = 1<<6,
                 // Features
                 eSDDF_Vertices   = 1<<7,
                 eSDDF_Edges      = 1<<8,
                 eSDDF_Faces      = 1<<9,
                 eSDDF_Volumes    = 1<<10,
                 eSDDF_Control    = 1<<11, //Control cage for Path and other spline/nurbs shapes
                 eSDDF_Normals    = 1<<12,
                 eSDDF_Tangents   = 1<<13,
                 // Misc
                 eSDDF_FeatureId  = 1<<14,
                 // Annotations
                 eSDDF_DCR        = 1<<15,
                 eSDDF_BVH        = 1<<16,
                 // \todo Other... eSDDF_Coplanar, eSDDF_Curvature..., eSDDF_Barycenter...
                 eSDDF_Everything = 0x0FFFFFFF,  //\note Upper 4 bits reserved for IObject DDF
                 eSDDF_Default    = eSDDF_Boundary };

enum EBVHDDF { eBVHDDF_Nothing    = 0,
               // Aspects
               eBVHDDF_Geometry   = 1<<29,
               eBVHDDF_Topology   = 1<<30,
               eBVHDDF_Levels     = 1<<31,
               // misc
               eBVHDDF_Everything = 0x0FFFFFFF,
               eBVHDDF_Default    = eBVHDDF_Geometry };

class IObject;
namespace bv { class IBoundingVolume; }
namespace bp { class IBroadPhase; }

void VizObject( const IObject* p_go, util::VizStream& vs, Flags32 ddf = eODDF_Default | eSDDF_Default );
//void VizShape2( const IShape2* p_shape, const Transform2& tr, const Real* vec_dof, util::VizStream& vs, Flags32 ddf = eODDF_Default | eSDDF_Default );
void VizBoundingVolume( const bv::IBoundingVolume* p_bv, util::VizStream& vs, Flags32 ddf = eODDF_Default | eSDDF_Default );
void VizBroadPhase( const bp::IBroadPhase* p_bp, util::VizStream& vs, Flags32 ddf = eODDF_Default | eSDDF_Default );

//\todo Consider internal GVizContactData<D> implementation, as they're IDENTICAL
enum EContactDataDDF { eContactDataDDF_Nothing        = 0,
                       eContactDataDDF_Points         = 1<<0,
                       eContactDataDDF_Normals        = 1<<1,
                       eContactDataDDF_AvgNormal      = 1<<2,
                       eContactDataDDF_Features       = 1<<3,
                       // Stochastic
                       eContactDataDDF_Stochastic_IP  = 1<<4,
                       eContactDataDDF_Stochastic_CP  = 1<<5,
                       eContactDataDDF_Stochastic_NP  = 1<<6,
                       eContactDataDDF_Stochastic_RNP = 1<<7,
                       eContactDataDDF_Stochastic_PCA = 1<<8,
                       // Stochastic Mapping
                       eContactDataDDF_Stochastic_IM  = 1<<9,
                       // \todo Bruteforce, GJK, SAT, specific flags...
                       eContactDataDDF_Everything     = 0xFFFFFFFF,
                       eContactDataDDF_Default        = eContactDataDDF_Points | eContactDataDDF_Normals };
void VizContactData( const geo::np::ContactData2& cd, util::VizStream& vs, Flags32 ddf = eContactDataDDF_Default );
void VizContactData( const geo::np::ContactData3& cd, util::VizStream& vs, Flags32 ddf = eContactDataDDF_Default );

//\todo Consider moving this somewhere else
// TEMP
class MeshSolidShape2;
class DCR_MeshSolidShape2;
void VizDCR( const DCR_MeshSolidShape2* p_dcr,
             const MeshSolidShape2* p_mss, const Transform2& tr, const Real* vec_dof,
             util::VizStream &vs, Flags32 ddf = 0xFFFFFFFF ); //\todo ENUM DDF?!?!!
class TetSolidShape3;
class DCR_TetSolidShape3;
void VizDCR( const DCR_TetSolidShape3* p_dcr,
             const TetSolidShape3* p_tss, const Transform3& tr, const Real* vec_dof,
             util::VizStream &vs, Flags32 ddf = 0xFFFFFFFF ); //\todo ENUM DDF?!?!!
/*
void VizBVH( const bv::BVH_ST_DG_Sphere2& bvh, util::VizStream &vs, Flags32 ddf = eBVHDDF_Default );
void VizBVH( const bv::BVH_ST_DG_AABB2& bvh, util::VizStream &vs, Flags32 ddf = eBVHDDF_Default );
void VizBVH( const bv::BVH_ST_DG_DOP2_K8& bvh, util::VizStream &vs, Flags32 ddf = eBVHDDF_Default );

void VizBVH( const bv::BVH_ST_DG_Sphere3& bvh, util::VizStream &vs, Flags32 ddf = eBVHDDF_Default );
void VizBVH( const bv::BVH_ST_DG_AABB3& bvh, util::VizStream &vs, Flags32 ddf = eBVHDDF_Default );
*/

template <typename BoundingVolumeT, typename EntryIndexT>
void GVizBVH( const bv::GBoundingVolumeHierarchy_ST_DG<BoundingVolumeT,EntryIndexT>& bvh, util::VizStream &vs, Flags32 ddf = eBVHDDF_Default );

}

#endif //GEO_UTIL_VIZ_H
