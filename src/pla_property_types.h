#ifndef PLA_TYPES_H
#error "You must only include pla_types.h!!"
#endif

#ifndef PLA_PROPERTY_TYPES_H
#define PLA_PROPERTY_TYPES_H

//! Number in a range property
template <typename T>
struct GProperty_NumberInRange
{
    GProperty_NumberInRange() : m_Value(0), m_Min(0), m_Max(0) {}
    GProperty_NumberInRange( T value, T min, T max ) : m_Value(value), m_Min(min), m_Max(max) {}
    T m_Value;
    T m_Min, m_Max;
};
typedef GProperty_NumberInRange<int32> Property_NIR_int32;
typedef GProperty_NumberInRange<uint32> Property_NIR_uint32;
typedef GProperty_NumberInRange<float> Property_NIR_float32;
typedef GProperty_NumberInRange<double> Property_NIR_float64;

//!\name Traits of property types
//@{
template <> struct pla_type_id< Property_NIR_int32 > { static const uint16 value = eType_Property_NIR_Int32; };
template <> struct pla_type_id< Property_NIR_uint32 > { static const uint16 value = eType_Property_NIR_UInt32; };
template <> struct pla_type_id< Property_NIR_float32 > { static const uint16 value = eType_Property_NIR_Float32; };
template <> struct pla_type_id< Property_NIR_float64 > { static const uint16 value = eType_Property_NIR_Float64; };
//@}

#endif //PLA_PROPERTY_TYPES_H
