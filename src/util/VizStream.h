#ifndef PLA_UTIL_VIZ_STREAM_H
#define PLA_UTIL_VIZ_STREAM_H

#include <pla_types.h>
#include "ItemStream.h"

namespace util {

typedef ItemStream VizStream;
typedef ItemStream::ItemIt VizItem;

/*! Enumerated Type and attributes of standard viz types

  All viz types are Complex, and contain mandatory geometry-specific
  attributes and an optional Style attribute.

  \name Named VizItems will display a Label at a standard positiion
  within them.
*/
enum EVizTypes {
    eVizTypesBase = cNumPlaTypes,

    // Viz types
    eType_VizSet,      //!< Named set/collection of viz data (eg: S2::eVizCollision....)

    eType_VizPoint,    //!< {pos=Vec3f}
    eType_VizSegment,  //!< {pos1=Vec3f,pos2=Vec3f}
    eType_VizVec,      //!< {pos=Vec3f,vec=Vec3f}
    eType_VizRefSys,   //!< {pos=Vec3f,rot=Quatf}

    eType_VizTriangle, //!< {p1,p2,p3=Vec3f}
    eType_VizQuad,     //!< {p1,p2,p3,p4=Vec3f}
    eType_VizRectangle,//!< {pos=Vec3f,rot=Quatf,half_sizes=Vec2f}
    eType_VizDisk,     //!< {pos=Vec3f,rot=Quatf,radius=f}

    eType_VizBox,      //!< {pos=Vec3f,rot=Quatf,half_sizes=Vec3f}
    eType_VizSphere,   //!< {pos=Vec3f,rot=Quatf,radius=float}
    eType_VizCapsule,  //!< {pos=Vec3f,axis=Vec3f,Radius=f,Height=f}
    eType_VizCylinder, //!< {pos=Vec3f,axis=Vec3f,Radius=f,Height=f}
    eType_VizTetrahedron, //!< {p0,p1,p2,p3=Vec3f}

    eType_VizAABB,     //!< {pos_min=Vec3f,pos_max=Vec3f}
    eType_VizGrid2D,   //!< {pos=Vec3f,dir1=Vec3f,dir2=Vec3f,size1=uint,size2=uint}

    eType_VizTransform2, //!< Transform2f
    eType_VizTransform3, //!< Transform3f

    eType_VizStyle,    //!< {color=Vec4f,pen_size=float,flags=Flags32}
};

enum EVizStyleFlags {
    eVizStyle_Nothing = 0,
    eVizStyle_Solid   = 1,
    eVizStyle_Wire    = 2
};

} // namespace util

//! Ugly viz macros are defined here
#define __PLA_UTIL_ENABLE_VIZ_MACROS
#include "VizMacros.h"

#endif // PLA_UTIL_VIZ_STREAM_H
