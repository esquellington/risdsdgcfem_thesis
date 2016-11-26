#ifndef PLA_TYPES_H
#error "You must only include pla_types.h!!"
#endif

#ifndef PLA_MATH_TYPES_H
#define PLA_MATH_TYPES_H

#include <Mal/GInterval.h>
#include <Mal/GVec.h>
#include <Mal/GMat.h>
#include <Mal/GTransform.h>
#include <Mal/GQuat.h>

//!\name Global math types
//@{
typedef mal::GInterval<float> Intervalf;
typedef mal::GVec<float,2> Vec2f;
typedef mal::GVec<float,3> Vec3f;
typedef mal::GVec<float,4> Vec4f;
typedef mal::GMat<float,2,2> Mat2x2f;
typedef mal::GMat<float,3,3> Mat3x3f;
typedef mal::GTransform<float,2> Transform2f;
typedef mal::GTransform<float,3> Transform3f;
typedef mal::GQuat<float> Quatf;
//@}

//!\name Traits of math types
//@{
template <> struct pla_type_id< Intervalf > { static const uint16 value = eType_Intervalf; };
template <> struct pla_type_id< Vec2f > { static const uint16 value = eType_Vec2f; };
template <> struct pla_type_id< Vec3f > { static const uint16 value = eType_Vec3f; };
template <> struct pla_type_id< Vec4f > { static const uint16 value = eType_Vec4f; };
template <> struct pla_type_id< Mat2x2f > { static const uint16 value = eType_Mat2x2f; };
template <> struct pla_type_id< Mat3x3f > { static const uint16 value = eType_Mat3x3f; };
template <> struct pla_type_id< Transform2f > { static const uint16 value = eType_Transform2f; };
template <> struct pla_type_id< Transform3f > { static const uint16 value = eType_Transform3f; };
template <> struct pla_type_id< Quatf > { static const uint16 value = eType_Quatf; };
//@}

#endif //PLA_MATH_TYPES_H
