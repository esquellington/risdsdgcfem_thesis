#ifndef PLA_TYPES_H
#error "You must only include pla_types.h!!"
#endif

#ifndef PLA_BASIC_TYPES_H
#define PLA_BASIC_TYPES_H

//------------------------------------------------------------------------------------
//!\name Basic types \todo USE c++11 basic types uint64_t, etc...
//@{
typedef void * void_pointer_type;
typedef long long int int64;
typedef long long unsigned int uint64;
/*\todo??
typedef long int long_int_type;
typedef long unsigned int long_uint_type;
*/
typedef int int32;
typedef unsigned int uint32;
typedef short int16;
typedef unsigned short uint16;
typedef char int8;
typedef unsigned char uint8;
typedef float float32;
typedef double float64;
//@}

#ifdef __LP64__
typedef int64 machine_int_type;
typedef uint64 machine_uint_type;
#else
typedef int32 machine_int_type;
typedef uint32 machine_uint_type;
#endif

//! Simple Flags template to make code more readable :-P
template <typename BitsT>
struct GFlags
{
    inline GFlags() : m_Bits(0) {}
    inline GFlags( BitsT flags ) : m_Bits(flags) {}
    inline bool Test( BitsT flags ) const { return (0 != (m_Bits & flags)); }
    inline bool Includes( BitsT flags ) const { return (flags == (m_Bits & flags)); }
    inline bool Any() const { return (0 != m_Bits); }
    inline bool None() const { return (0 == m_Bits); }
    inline GFlags& Enable( BitsT flags ) { m_Bits |= flags; return *this; }
    inline GFlags& Disable( BitsT flags ) { m_Bits &= ~flags; return *this; }
    inline GFlags& Toggle( BitsT flags ) { m_Bits = m_Bits ^ flags; return *this; }
    inline GFlags& Set( BitsT flags ) { m_Bits = flags; return *this; }
    inline GFlags& Reset() { m_Bits = BitsT(0); return *this; }
    inline operator BitsT() const { return m_Bits; }
    BitsT m_Bits;
};
typedef GFlags<int32> Flags32;
typedef GFlags<int16> Flags16;
typedef GFlags<int8> Flags8;

/*! String with fixed memory allocation of L bytes, and maximum length
  L-2 (+1 for null-termination +1 for m_Length)
*/
template <unsigned L>
struct GString
{
    inline GString() { Clear(); }
    inline GString( const char *str ) { Set(str); }
    inline GString( const char *str, unsigned int length ) { Set(str,length); }
    inline void Clear() { m_Length = 0; m_vecChar[0] = 0; }
    inline bool Set( const char *str ) //!< 0-ended string
        {
            m_Length = 0;
            for( unsigned i=0; i<L-2 && str[i] != 0; i++, m_Length++ ) m_vecChar[i] = str[i];
            m_vecChar[m_Length] = 0;
            return str[m_Length] == 0;
        }
    inline bool Set( const char *str, unsigned int length ) //!< non-0-ended string
        {
            m_Length = 0;
            for( unsigned int i=0; i<L-2 && i<length; i++, m_Length++ ) m_vecChar[i] = str[i];
            m_vecChar[m_Length] = 0;
            return str[m_Length] == 0;
        }
    inline uint32 GetLength() const { return m_Length; }
    inline uint32 GetMaxLength() const { return L-2; }
    inline const char *GetStr() const { return reinterpret_cast<const char *>(&m_vecChar); }
    inline char *GetStr() { return reinterpret_cast<char *>(&m_vecChar); }
    inline operator const char *() const { return reinterpret_cast<const char *>(&m_vecChar); }
    inline operator char *() { return reinterpret_cast<char *>(&m_vecChar); }
    char m_vecChar[L-1];
    uint8 m_Length;
};

typedef GString<16> String16;
typedef GString<32> String32;
typedef GString<64> String64;

//\todo UGLY GPair type to export time-value samples, monitor and plot them... will probably change
template <typename T1, typename T2>
struct GPair
{
    inline GPair() : m_First(), m_Second() {}
    inline GPair( const T1 &first, const T2 &second ) : m_First(first), m_Second(second) {}
    T1 m_First;
    T2 m_Second;
};
typedef GPair<float32,float32> Pair_f32_f32;

//!\name Traits of basic types
//@{
template <> struct pla_type_id<void_pointer_type> { static const uint16 value = eType_VoidPointer; };

template <> struct pla_type_id<int64> { static const uint16 value = eType_Int64; };
template <> struct pla_type_id<uint64> { static const uint16 value = eType_UInt64; };
template <> struct pla_type_id<int32> { static const uint16 value = eType_Int32; };
template <> struct pla_type_id<uint32> { static const uint16 value = eType_UInt32; };
template <> struct pla_type_id<int16> { static const uint16 value = eType_Int16; };
template <> struct pla_type_id<uint16> { static const uint16 value = eType_UInt16; };
template <> struct pla_type_id<int8> { static const uint16 value = eType_Int8; };
template <> struct pla_type_id<uint8> { static const uint16 value = eType_UInt8; };

template <> struct pla_type_id<bool> { static const uint16 value = eType_Bool32; };

template <> struct pla_type_id<float32> { static const uint16 value = eType_Float32; };
template <> struct pla_type_id<float64> { static const uint16 value = eType_Float64; };

template <> struct pla_type_id< Flags8 > { static const uint16 value = eType_Flags8; };
template <> struct pla_type_id< Flags16 > { static const uint16 value = eType_Flags16; };
template <> struct pla_type_id< Flags32 > { static const uint16 value = eType_Flags32; };

template <> struct pla_type_id< String16 > { static const uint16 value = eType_String16; };
template <> struct pla_type_id< String32 > { static const uint16 value = eType_String32; };
template <> struct pla_type_id< String64 > { static const uint16 value = eType_String64; };

template <> struct pla_type_id< Pair_f32_f32 > { static const uint16 value = eType_Pair_Float32_Float32; };
//@}

#endif //PLA_BASIC_TYPES
