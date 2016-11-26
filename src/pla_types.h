#ifndef PLA_TYPES_H
#define PLA_TYPES_H

enum EPlaBasicTypes {
    eType_Unknown = 0,

    eType_VoidPointer, //!< void* SIMPLE type. Size is platform-dependent, 32 or 64b, CANNOT be serialized.

    eType_Int64,
    eType_UInt64,
    eType_Int32,
    eType_UInt32,
    eType_Int16,
    eType_UInt16,
    eType_Int8,
    eType_UInt8,

    eType_Bool32,

    eType_Float32,
    eType_Float64,

    eType_Flags8,
    eType_Flags16,
    eType_Flags32,

    eType_String16,
    eType_String32,
    eType_String64,

    cLastPlaBasicType = eType_String64
};

enum EPlaMathTypes {
    eType_Intervalf = cLastPlaBasicType + 1,
    eType_Vec2f,
    eType_Vec3f,
    eType_Vec4f,
    eType_Mat2x2f,
    eType_Mat3x3f,
    eType_Transform2f,
    eType_Transform3f,
    eType_Quatf,

    cLastPlaMathType = eType_Quatf
};

//! Enumeration of standard types
enum EPlaPropertyTypes {
    eType_Property = cLastPlaMathType + 1,

    eType_Property_Object,
    eType_Property_Group,

    eType_Property_NIR,
    eType_Property_NIR_Int32,
    eType_Property_NIR_UInt32,
    eType_Property_NIR_Float32,
    eType_Property_NIR_Float64,

    eType_Property_Enum,
    eType_Property_Flags,

    eType_Pair_Float32_Float32,

    /*\todo
    eType_Property_Position2,
    eType_Property_Position3,
    eType_Property_Direction2,
    eType_Property_Direction3,
    eType_Property_Angle,
    eType_Property_Quaternion,
    eType_Property_Scale2,
    eType_Property_Scale3
    */

    cNumPlaTypes, //!< ALL user-types should be >
};

// pLa type-system
typedef unsigned short int PlaTypeID;
//! type-trait for automatic identification of simple types (\note MUST BE IN GLOBAL NAMESPACE!!)
template <typename T> struct pla_type_id { static const PlaTypeID value = eType_Unknown; };

#include "pla_basic_types.h"
#include "pla_math_types.h"
#include "pla_property_types.h"

#endif //PLA_TYPES_H
