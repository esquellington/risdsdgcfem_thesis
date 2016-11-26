#ifndef UTIL_ITEM_STREAM_H
#define UTIL_ITEM_STREAM_H

#include <util/Config.h>
#include <util/MemBlock.h>

namespace util {

#define __ENABLE_ITEM_STREAM_REALLOC

#define UTIL_ITEM_STREAM_MAGICNUMBER uint32( (uint8('U')<<24) | (uint8('I')<<16) | (uint8('S')<<8) | uint8('M') )
#define UTIL_ITEM_STREAM_VERSION uint32(0)

/*! Stream/Sequence of generic Items
Supports:
- Sequential addition of Identified and type-checked, Simple, Array
  and (nested) Complex data types
- Sequential and hierarchical iteration and item retrieval with bounds
  and type-checking
- Named items using either positive integer identifiers or strings (using a
  fixed-size string pool)
- 4-byte aligned memory
- Automatic realloc, if b_can_realloc = true

\todo Missing functionality and optimizations:
- Guarantee proper alignment for Arrays/Structs
- Find Recursiu (complex items)
- Find filtrat (per type <T>)
- Iterator que filtra per tipus <T>
- Optimize GetStringOffset()
  - Using some sort of string-hash map (or better: Id32/64?)
- Map() a functor onto all items of the stream
\todo Consider adding a RawMem entry type with a DIFFERENT header for
      LARGE blocks (uint32) and limit size_type to uint16 to avoid
      enlarging all headers...
*/
class ItemStream
{
public:
    enum EReallocFlags {
        eRealloc_Nothing     = 0,
        eRealloc_Data        = 1,
        eRealloc_Identifiers = 2,
        eRealloc_Everything  = (eRealloc_Data | eRealloc_Identifiers),
        eRealloc_Default     = (eRealloc_Data | eRealloc_Identifiers)
    };

public:

    typedef uint32 size_type; //\todo This enlarges ALL HEADERS, but is required to support large binarrays (uint8)

    // First 4 bits of type descriptor are reserved for {Simple,Pointer,Array,Complex} categories
    static const uint16 cItemType_Simple   = 0x1000;
    static const uint16 cItemType_Pointer  = 0x2000;
    static const uint16 cItemType_Array    = 0x4000;
    static const uint16 cItemType_Complex  = 0x8000;
    // 5th upper bit is reserved for IsTouched
    static const uint16 cItemFlag_Touched  = 0x0800;
    // Masks and limits
    static const uint16 cMax_ItemType      = 0x07FF;
    static const uint16 cMask_ItemType     = 0x07FF;
    static const uint16 cMask_ItemMetaType = 0xF000;

    static const int cMaxNestedComplex = 32; //!< Completely arbitrary

    //! Simple POD and Complex elements header
    struct ItemHeader
    {
        inline ItemHeader( int32 id, uint16 type, size_type size ) : m_Id(id), m_Type(type), m_Size(size) {}
        int32 m_Id;    //!< id>0 => id=identifier, id<0 => id=string map key, id=0 => Invalid
        uint16 m_Type; //!< [15=>IsArray, 14=>IsComplex, 13..0=>EItemType]
        size_type m_Size; //!< Item size in bytes ( including {Item|Array}Header )
    };

    //! Array of POD elements header
    struct ArrayHeader: public ItemHeader
    {
        inline ArrayHeader( int32 id, uint16 type, size_type size, uint32 count )
        : ItemHeader(id,type,size), m_Count(count) {}
        uint32 m_Count;
    };

    //! Item iterator
    class ItemIt
    {
    public:
        inline ItemIt() : m_pIS(0), m_BeginOffset(0), m_EndOffset(0) {}
        inline ItemIt( const ItemStream &item_stream, uint32 begin_offset = 0, uint32 end_offset = 0 )
        : m_pIS(&item_stream), m_BeginOffset(begin_offset), m_EndOffset(end_offset) {}
        inline ~ItemIt() {}

        //! \name Structure consultors
        //@{
        inline bool IsSimple() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Simple & GetItemHeader().m_Type); }
        inline bool IsPointer() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Pointer & GetItemHeader().m_Type); }
        inline bool IsArray() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Array & GetItemHeader().m_Type); }
        inline bool IsComplex() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Complex & GetItemHeader().m_Type); }
        //@}

        //! \name Header consultors
        //@{
        inline int32 GetId() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Id; }
        inline uint16 GetType() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Type & cMask_ItemType; }
        inline uint16 GetMetaType() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Type & cMask_ItemMetaType; }
        inline uint32 GetFullType() const
        {
            UTIL_ASSERT(IsValid());
            return (IsArray()) ? (uint32(GetItemHeader().m_Type) << 16) | GetArrayCount()
                               : uint32(GetItemHeader().m_Type);
        }
        inline size_type GetSize() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Size; }
        inline bool IsNamed() const { UTIL_ASSERT(IsValid()); return GetId() < 0; }
        inline const char *GetName() const { UTIL_ASSERT(IsValid()); UTIL_ASSERT(GetId()<0); return m_pIS->GetString(-GetId()); }
        //@}

        //! \name Validity checking
        //@{
        inline bool IsValid() const { return m_BeginOffset < m_EndOffset; }
        //inline operator bool() const { return IsValid(); }
        //@}

        //! \name Touching
        //@{
        inline bool IsTouched() const { UTIL_ASSERT(IsValid()); return 0 != (cItemFlag_Touched & GetItemHeader().m_Type); }
        //@}

        //! Get simple type element (with type-checking)
        template <typename T> inline const T &Get() const
        {
            UTIL_ASSERT( IsValid() && IsSimple() && pla_type_id<T>::value == GetType() );
            return *reinterpret_cast<const T*>(GetDataPtr());
        }
        //! Safe get for simple type element with "default" value support
        template <typename T> inline T SafeGet( const T &default_value ) const
        {
            if( IsValid() ) return Get<T>();
            else return default_value;
        }
        //! Get as type T, casting from internal type if possible
        template <typename T> inline T CastGet( const T &default_value = T() ) const
            {
                UTIL_ASSERT( IsSimple() );
                switch( GetType() )
                {
                case eType_Int64: return T( Get<int64>() ); break;
                case eType_UInt64: return T( Get<uint64>() ); break;
                case eType_Int32: return T( Get<int32>() ); break;
                case eType_UInt32: return T( Get<uint32>() ); break;
                case eType_Int16: return T( Get<int16>() ); break;
                case eType_UInt16: return T( Get<uint16>() ); break;
                case eType_Int8: return T( Get<int8>() ); break;
                case eType_UInt8: return T( Get<uint8>() ); break;
                case eType_Float32: return T( Get<float32>() ); break;
                case eType_Float64: return T( Get<float64>() ); break;
                default: return default_value; break;
                }
            }
        //! Get as number type T, even if it's a different number type T1 or a GNumberInRange<T> type
        template <typename T> inline T GetNumber( const T &default_value = T() ) const
            {
                UTIL_ASSERT( IsSimple() );
                if( IsValid() )
                {
                    if( GetType() == pla_type_id<T>::value )
                        return Get<T>();
                    else if( GetType() == pla_type_id< GProperty_NumberInRange<T> >::value )
                        return Get< GProperty_NumberInRange<T> >().m_Value;
                    else
                        return CastGet<T>();
                }
                else return default_value;
            }
        //! Get Pointer type element (with type-checking)
        template <typename T> inline const T *GetPtr() const
        {
            UTIL_ASSERT( IsValid() && IsPointer() && pla_type_id<T>::value == GetType() );
            return reinterpret_cast<const T*>(GetDataPtr());
        }
        //! Get simple type array (with type-checking)
        template <typename T> inline const T *GetArrayPtr() const
        {
            UTIL_ASSERT( IsValid() && IsArray() && pla_type_id<T>::value == GetType() );
            return reinterpret_cast<const T*>(GetDataPtr()+4);
        }
        //! Get array with default value if !IsValid
        template <typename T> inline const T *SafeGetArrayPtr( const T *default_ptr ) const
        {
            if(IsValid()) return GetArrayPtr<T>();
            else return default_ptr;
        }
        inline uint32 GetArrayCount() const
        {
            UTIL_ASSERT( IsValid() && IsArray() );
            return GetArrayHeader().m_Count;
        }
        inline const char *GetString() const { return GetArrayPtr<char>(); }
        inline const char *SafeGetString( const char *default_str ) const { return SafeGetArrayPtr<char>(default_str); }

        //! Get complex item iterator pointing at the first SubItem
        inline ItemIt GetSubItem() const
        {
            UTIL_ASSERT( IsValid() && IsComplex() );
            return ItemIt( *m_pIS, m_BeginOffset+sizeof(ItemHeader), m_BeginOffset+GetSize() );
        }

        //! Advance to next item (in the same level)
        inline ItemIt &operator ++()
        {
            UTIL_ASSERT( IsValid() );
            m_BeginOffset += GetSize();
            return *this;
        }
        //! Get next, without advancing this
        inline ItemIt Next() const
        {
            UTIL_ASSERT( IsValid() );
            return ItemIt( *m_pIS, m_BeginOffset+GetSize(), m_EndOffset );
        }

        //! Find an Item by ID (in the same level) after the current one
        ItemIt Find( int32 id ) const;
        //! Find an Item by name (in the same level) after the current one
        ItemIt Find( const char *name ) const;

        //! Raw-copy data to buffer
        void ToBuffer( int8 *p_buffer ) const;

        //! Comparison, works on !IsValid() iterators too
        inline bool operator==( const ItemIt &other ) const
        {
            return m_pIS == other.m_pIS
                && m_BeginOffset == other.m_BeginOffset
                && m_EndOffset == other.m_EndOffset;
        }
        inline bool operator!=( const ItemIt &other ) const { return !(*this == other); }

    public:
        inline const int8 *GetDataPtr() const { UTIL_ASSERT(IsValid()); return &(*m_pIS)[m_BeginOffset+sizeof(ItemHeader)]; }

    private:
        inline const ItemHeader &GetItemHeader() const { UTIL_ASSERT(IsValid()); return *reinterpret_cast<const ItemHeader*>(&(*m_pIS)[m_BeginOffset]); }
        inline const ArrayHeader &GetArrayHeader() const { UTIL_ASSERT( IsArray() ); return *reinterpret_cast<const ArrayHeader*>(&(*m_pIS)[m_BeginOffset]); }

        friend class ItemStream;

    private:
        const ItemStream *m_pIS;
        uint32 m_BeginOffset;
        uint32 m_EndOffset;
    };

    //! Item iterator
    class ItemItRW
    {
    public:
        inline ItemItRW() : m_pIS(0), m_BeginOffset(0), m_EndOffset(0) {}
        inline ItemItRW( ItemStream &item_stream, uint32 begin_offset = 0, uint32 end_offset = 0 )
        : m_pIS(&item_stream), m_BeginOffset(begin_offset), m_EndOffset(end_offset) {}
        inline ~ItemItRW() {}

        //! \name Structure consultors
        //@{
        inline bool IsSimple() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Simple & GetItemHeader().m_Type); }
        inline bool IsPointer() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Pointer & GetItemHeader().m_Type); }
        inline bool IsArray() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Array & GetItemHeader().m_Type); }
        inline bool IsComplex() const { UTIL_ASSERT(IsValid()); return 0 != (cItemType_Complex & GetItemHeader().m_Type); }
        //@}

        //! \name Header consultors
        //@{
        inline int32 GetId() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Id; }
        inline uint16 GetType() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Type & cMask_ItemType; }
        inline uint16 GetMetaType() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Type & cMask_ItemMetaType; }
        inline uint32 GetFullType() const
        {
            UTIL_ASSERT(IsValid());
            return (IsArray()) ? (uint32(GetItemHeader().m_Type) << 16) | GetArrayCount()
                               : uint32(GetItemHeader().m_Type);
        }
        inline size_type GetSize() const { UTIL_ASSERT(IsValid()); return GetItemHeader().m_Size; }
        inline bool IsNamed() const { UTIL_ASSERT(IsValid()); return GetId() < 0; }
        inline const char *GetName() const { UTIL_ASSERT(GetId()<0); return m_pIS->GetString(-GetId()); }
        //@}

        //! \name Validity checking
        //@{
        inline bool IsValid() const { return m_BeginOffset < m_EndOffset; }
        //inline operator bool() const { return IsValid(); }
        //@}

        //! \name Touching
        //@{
        inline bool IsTouched() const { UTIL_ASSERT(IsValid()); return 0 != (cItemFlag_Touched & GetItemHeader().m_Type); }
        inline void Touch() { UTIL_ASSERT(IsValid()); GetItemHeader().m_Type = uint16(GetItemHeader().m_Type | cItemFlag_Touched); }
        inline void Untouch() { UTIL_ASSERT(IsValid()); GetItemHeader().m_Type = uint16(GetItemHeader().m_Type & ~cItemFlag_Touched); }
        void TouchRecursive();
        void UntouchRecursive();
        //@}

        //! Get simple type element (with type-checking)
        template <typename T> inline void Set( T value )
        {
            UTIL_ASSERT( IsValid() && IsSimple() && pla_type_id<T>::value == GetType() );
            *reinterpret_cast<T*>(GetDataPtr()) = value;
            Touch();
        }
        template <typename T> inline T &Get() //Does NOT automatically Touch() the element if changed
        {
            UTIL_ASSERT( IsValid() && IsSimple() && pla_type_id<T>::value == GetType() );
            return *reinterpret_cast<T*>(GetDataPtr());
        }
        template <typename T> inline const T &Get() const
        {
            UTIL_ASSERT( IsValid() && IsSimple() && pla_type_id<T>::value == GetType() );
            return *reinterpret_cast<const T*>(GetDataPtr());
        }
        //! Safe get for simple type element with "default" value support
        template <typename T> inline T SafeGet( const T &default_value ) const
        {
            if(IsValid()) return Get<T>();
            else return default_value;
        }
        //! Get as type T, casting from internal type if possible
        template <typename T> inline T CastGet( const T &default_value = T() ) const
            {
                UTIL_ASSERT( IsSimple() );
                switch( GetType() )
                {
                case eType_Int64: return T( Get<int64>() ); break;
                case eType_UInt64: return T( Get<uint64>() ); break;
                case eType_Int32: return T( Get<int32>() ); break;
                case eType_UInt32: return T( Get<uint32>() ); break;
                case eType_Int16: return T( Get<int16>() ); break;
                case eType_UInt16: return T( Get<uint16>() ); break;
                case eType_Int8: return T( Get<int8>() ); break;
                case eType_UInt8: return T( Get<uint8>() ); break;
                case eType_Float32: return T( Get<float32>() ); break;
                case eType_Float64: return T( Get<float64>() ); break;
                    //\todo All other numeric types
                default: return default_value; break;
                }
            }
        //! Get as number type T, even if it's a different number type T1 or a GNumberInRange<T> type
        template <typename T> inline T GetNumber( const T &default_value = T() ) const
            {
                UTIL_ASSERT( IsSimple() );
                if( IsValid() )
                {
                    if( GetType() == pla_type_id<T>::value )
                        return Get<T>();
                    else if( GetType() == pla_type_id< GProperty_NumberInRange<T> >::value )
                        return Get< GProperty_NumberInRange<T> >().m_Value;
                    else
                        return CastGet<T>();
                }
                else return default_value;
            }
        //! Get Pointer type element (with type-checking)
        template <typename T> inline T *GetPtr()
        {
            UTIL_ASSERT( IsValid() && IsPointer() && pla_type_id<T>::value == GetType() );
            return reinterpret_cast<T*>(GetDataPtr());
        }
        //! Get simple type array (with type-checking)
        template <typename T> inline T *GetArrayPtr()
        {
            UTIL_ASSERT( IsValid() && IsArray() && pla_type_id<T>::value == GetType() );
            return reinterpret_cast<T*>(GetDataPtr()+4);
        }
        //! Get simple type array (with type-checking)
        template <typename T> inline const T *GetArrayPtr() const
        {
            UTIL_ASSERT( IsValid() && IsArray() && pla_type_id<T>::value == GetType() );
            return reinterpret_cast<T*>(GetDataPtr()+4);
        }
        //! Get array with default value if !IsValid
        template <typename T> inline const T *SafeGetArrayPtr( const T *default_ptr ) const
        {
            if(IsValid()) return GetArrayPtr<T>();
            else return default_ptr;
        }
        inline uint32 GetArrayCount() const
        {
            UTIL_ASSERT( IsValid() && IsArray() );
            return GetArrayHeader().m_Count;
        }
        inline char *GetString() { return GetArrayPtr<char>(); }
        inline const char *SafeGetString( const char *default_str ) const { return SafeGetArrayPtr<char>(default_str); }

        //! Get complex item iterator pointing at the first SubItem
        inline ItemItRW GetSubItem()
        {
            UTIL_ASSERT( IsValid() && IsComplex() );
            return ItemItRW( *m_pIS, m_BeginOffset+sizeof(ItemHeader), m_BeginOffset+GetSize() );
        }

        //! Advance to next item (in the same level)
        inline ItemItRW &operator ++()
        {
            UTIL_ASSERT( IsValid() );
            m_BeginOffset += GetSize();
            return *this;
        }
        //! Get next, without advancing this
        inline ItemItRW Next() const
        {
            UTIL_ASSERT( IsValid() );
            return ItemItRW( *m_pIS, m_BeginOffset+GetSize(), m_EndOffset );
        }

        //! Find an Item by ID (in the same level) after the current one
        ItemItRW Find( int32 id ) const;
        //! Find an Item by name (in the same level) after the current one
        ItemItRW Find( const char *name ) const;

        //! Raw-copy data to buffer
        void ToBuffer( int8 *p_buffer ) const;

        //! Comparison, works on !IsValid iterators too
        inline bool operator==( const ItemItRW &other ) const
        {
            return m_pIS == other.m_pIS
                && m_BeginOffset == other.m_BeginOffset
                && m_EndOffset == other.m_EndOffset;
        }
        inline bool operator!=( const ItemItRW &other ) const { return !(*this == other); }

        // Cast to const itemit
        inline operator ItemIt() const { return ItemIt( *m_pIS, m_BeginOffset, m_EndOffset ); }

    public:
        inline int8 *GetDataPtr() const { UTIL_ASSERT(IsValid()); return &(*m_pIS)[m_BeginOffset+sizeof(ItemHeader)]; }

    private:
        inline ItemHeader &GetItemHeader() const { UTIL_ASSERT(IsValid()); return *reinterpret_cast<ItemHeader*>(&(*m_pIS)[m_BeginOffset]); }
        inline ArrayHeader &GetArrayHeader() const { UTIL_ASSERT( IsArray() ); return *reinterpret_cast<ArrayHeader*>(&(*m_pIS)[m_BeginOffset]); }

        friend class ItemStream;

    private:
        ItemStream *m_pIS;
        uint32 m_BeginOffset;
        uint32 m_EndOffset;
    };

public:
    //Constructor: By default, no memory is reserved in MemBlocks
    inline ItemStream( unsigned int max_size_in_bytes = 0,
                       unsigned int max_identifier_mem_size = 0,
                       Flags32 realloc_flags = eRealloc_Default )
    : m_DataMB(max_size_in_bytes)
    , m_StringMB(max_identifier_mem_size)
#ifdef __ENABLE_ITEM_STREAM_REALLOC
    , m_ReallocFlags(realloc_flags)
#endif
    , m_BeginOffset(0), m_EndOffset(0)
    , m_StringEndOffset(1) //Starts at 1
    , m_NestedComplexLevel(0)
    {}
    inline ~ItemStream() {}

    inline void Init( unsigned int max_size_in_bytes = 0,
                      unsigned int max_identifier_mem_size = 0,
                      Flags32 realloc_flags = eRealloc_Default )
    {
        Clear();
        m_DataMB.Init(max_size_in_bytes);
        m_StringMB.Init(max_identifier_mem_size);
#ifdef __ENABLE_ITEM_STREAM_REALLOC
        m_ReallocFlags = realloc_flags;
#endif
    }

    inline void SetReallocFlags( Flags32 realloc_flags ) { m_ReallocFlags = realloc_flags; }
    inline Flags32 GetReallocFlags() const { return m_ReallocFlags; }

    inline void Clear()
    {
        m_BeginOffset = 0;
        m_EndOffset = 0;
        m_StringEndOffset = 1; //Starts at 1
        m_NestedComplexLevel = 0;
    }

    inline bool IsEmpty() const { return 0 == m_EndOffset; }
// #ifdef __ENABLE_ITEM_STREAM_REALLOC
//     inline bool CanRealloc() const { return m_bCanRealloc; }
// #else
//     inline bool CanRealloc() const { return false; }
// #endif

    //! \name Item access, search and count
    //@{
    inline ItemIt Begin() const { return ItemIt( *this, m_BeginOffset, m_EndOffset ); }
    inline ItemIt Find( int32 id ) const { return Begin().Find( id ); }
    inline ItemIt Find( const char *name ) const { return Begin().Find( name ); }

    inline ItemItRW BeginRW() { return ItemItRW( *this, m_BeginOffset, m_EndOffset ); }
    inline ItemItRW FindRW( int32 id ) { return BeginRW().Find( id ); }
    inline ItemItRW FindRW( const char *name ) { return BeginRW().Find( name ); }

    //! Linear-time count (first level only)
    inline uint32 Count() const { uint32 c=0; for( ItemIt it=Begin(); it.IsValid(); ++it, ++c ); return c; }
    //@}

    //!\name Sequential write
    //@{
    // Simple types (req pla_type_id<T> trait and T=POD)
    template <typename T> inline void Write( int32 id, const T& data )
    {
        UTIL_ASSERT( id > 0 );
        Write_Internal( id, data );
    }
    template <typename T> inline void Write( const char *name, const T& data )
    {
        uint32 id = GetStringOffset( name );
        Write_Internal( -id, data );
    }
    // Deep-Copy an Item, possibly from another ItemStream
    inline ItemItRW WriteItem( int32 id, const ItemIt &it )
    {
        UTIL_ASSERT( id > 0 );
        return WriteItem_Internal( id, it );
    }
    inline ItemItRW WriteItem( const char *name, const ItemIt &it )
    {
        uint32 id = GetStringOffset( name );
        return WriteItem_Internal( -id, it );
    }

    // Pointer types (req pla_type_id<T> trait and T=POD)
    template <typename T> inline void WritePtr( int32 id, const T* ptr )
    {
        UTIL_ASSERT( id > 0 );
        WritePtr_Internal( id, ptr );
    }
    template <typename T> inline void WritePtr( const char *name, const T* ptr )
    {
        uint32 id = GetStringOffset( name );
        WritePtr_Internal( -id, ptr );
    }
    // Simple type arrays (req pla_type_id<T> trait and T=POD)
    template <typename T> inline void WriteArray( int32 id, const T* data, unsigned int num_elems )
    {
        UTIL_ASSERT( id > 0 );
        WriteArray_Internal( id, data, num_elems );
    }
    template <typename T> inline void WriteArray( const char *name,
                                                  const T* data, unsigned int num_elems )
    {
        uint32 id = GetStringOffset( name );
        WriteArray_Internal( -id, data, num_elems );
    }
    // Simple type array allocation (req pla_type_id<T> trait and T=POD)
    // Returns insertion pointer for allocated elements
    template <typename T> inline T* AllocArray( int32 id, unsigned int num_elems )
    {
        UTIL_ASSERT( id > 0 );
        return AllocArray_Internal<T>( id, num_elems );
    }
    template <typename T> inline T* AllocArray( const char *name, unsigned int num_elems )
    {
        uint32 id = GetStringOffset( name );
        return AllocArray_Internal<T>( -id, num_elems );
    }

    // String type, helper on Array of Char
    void WriteString( int32 id, const char* str, int length = -1 );
    void WriteString( const char *name, const char* str, int length = -1 );

    // Complex types (item-type must be explicitly given)
    void BeginComplex( int32 id, uint16 type );
    void BeginComplex( const char *name, uint16 type );
    /*! Ends last begun complex (computes size), and returns an Iterator to the Complex Item begining */
    ItemItRW EndComplex();

    //! Parse Item definitions from a string, return first non-consumed character
    const char *Parse( const char *str );
    //@}

    //! \name Quick-and-dirty Load/Save to file
    //@{
    bool LoadBin( const char *file_name );
    bool SaveBin( const char *file_name ) const;
    bool LoadTxt( const char *file_name );
    bool SaveTxt( const char *file_name ) const;
    //@}

private:
    friend class ItemIt;
    friend class ItemItRW;

    bool ReallocDataMB( size_type required_capacity )
    {
        size_type capacity( m_DataMB.GetSize() );
#ifdef __ENABLE_ITEM_STREAM_REALLOC
        if( required_capacity >= capacity )
        {
            if( m_ReallocFlags.Test(eRealloc_Data) )
            {
                size_type new_capacity = (required_capacity < 2*capacity) ? 2*capacity : required_capacity;
                return m_DataMB.Realloc( new_capacity );
            }
            else
                return false;
        }
        else
            return true;
#else
        return required_capacity < capacity;
#endif
    }

    bool ReallocStringMB( size_type required_capacity )
    {
        size_type capacity( m_StringMB.GetSize() );
#ifdef __ENABLE_ITEM_STREAM_REALLOC
        if( required_capacity >= capacity )
        {
            if( m_ReallocFlags.Test(eRealloc_Identifiers) )
            {
                size_type new_capacity = (required_capacity < 2*capacity) ? 2*capacity : required_capacity;
                return m_StringMB.Realloc( new_capacity );
            }
            else
                return false;
        }
        else
            return true;
#else
        return required_capacity < capacity;
#endif
    }

    //! \name Internal data writing methods
    //@{
    //! Simple data
    template <typename T> inline void Write_Internal( int32 id, const T& data )
    {
        UTIL_ASSERT( pla_type_id<T>::value < cMax_ItemType );

        //Last term computes 4-aligned sizeof
        machine_uint_type required_size = sizeof(ItemHeader) + 4*((sizeof(T)+3)/4);
        size_type size( required_size );
        UTIL_ASSERT( required_size == size );

        bool bReallocOk = ReallocDataMB( m_EndOffset + size );
        UTIL_ASSERT( bReallocOk );

        ItemHeader *p_header = (ItemHeader*)&m_DataMB[m_EndOffset];
        T* p_data = (T*)(&m_DataMB[m_EndOffset] + sizeof(ItemHeader));
        *p_header = ItemHeader( id, pla_type_id<T>::value | cItemType_Simple, size);
        *p_data = data;
        m_EndOffset += size;
    }

    // Deep-Copy an Item, possibly from another ItemStream
    ItemItRW WriteItem_Internal( int32 id, const ItemIt &it );

    //! Pointer data
    template <typename T> inline void WritePtr_Internal( int32 id, const T* ptr )
    {
        UTIL_ASSERT( pla_type_id<T>::value < cMax_ItemType );

        //Last term computes 4-aligned sizeof
        machine_uint_type required_size = sizeof(ItemHeader) + 4*((sizeof(T)+3)/4);
        size_type size( required_size );
        UTIL_ASSERT( required_size == size );

        bool bReallocOk = ReallocDataMB( m_EndOffset + size );
        UTIL_ASSERT( bReallocOk );

        ItemHeader *p_header = (ItemHeader*)&m_DataMB[m_EndOffset];
        T** p_data = (T**)(&m_DataMB[m_EndOffset] + sizeof(ItemHeader));
        *p_header = ItemHeader( id, pla_type_id<T>::value | cItemType_Pointer, size);
        *p_data = (T*)ptr; //Must cast to non-const T* here for some stupid reason...
        m_EndOffset += size;
    }
    //! Simple type arrays
    template <typename T> inline void WriteArray_Internal( int32 id,
                                                           const T* data, unsigned int num_elems )
    {
        UTIL_ASSERT( pla_type_id<T>::value < cMax_ItemType );

        //Last term computes 4-aligned sizeof
        machine_uint_type required_size = sizeof(ArrayHeader) + 4*((num_elems*sizeof(T)+3)/4);
        size_type size( required_size );
        UTIL_ASSERT( required_size == size );

        bool bReallocOk = ReallocDataMB( m_EndOffset + size );
        UTIL_ASSERT( bReallocOk );

        ArrayHeader *p_header = (ArrayHeader*)&m_DataMB[m_EndOffset];
        T* p_data = (T*)(&m_DataMB[m_EndOffset] + sizeof(ArrayHeader));
        *p_header = ArrayHeader( id, pla_type_id<T>::value | cItemType_Array, size, num_elems );
        for( unsigned int i=0; i<num_elems; i++ ) p_data[i] = data[i];
        m_EndOffset += size;
    }
    // Simple type array allocation that returns insertion pointer
    template <typename T> inline T* AllocArray_Internal( int32 id, unsigned int num_elems )
    {
        UTIL_ASSERT( pla_type_id<T>::value < cMax_ItemType );

        //Last term computes 4-aligned sizeof
        machine_uint_type required_size = sizeof(ArrayHeader) + 4*((num_elems*sizeof(T)+3)/4);
        size_type size( required_size );
        UTIL_ASSERT( required_size == size );

        bool bReallocOk = ReallocDataMB( m_EndOffset + size );
        UTIL_ASSERT( bReallocOk );

        ArrayHeader *p_header = (ArrayHeader*)&m_DataMB[m_EndOffset];
        T* p_data = (T*)(&m_DataMB[m_EndOffset] + sizeof(ArrayHeader));
        *p_header = ArrayHeader( id, pla_type_id<T>::value | cItemType_Array, size, num_elems );
        m_EndOffset += size;
        return p_data;
    }
    //! Complex types
    void BeginComplex_Internal( int32 id, uint16 type );
    //@}

private:
    //! \name low-level byte access
    //@{
    inline int8 &operator[]( uint32 offset ) { return m_DataMB[offset]; }
    inline const int8 &operator[]( uint32 offset ) const { return m_DataMB[offset]; }
    //@}

    //! \name low-level identifier allocation and retrieval
    //@{
    inline const char *GetString( uint32 offset ) const { return (const char*)&m_StringMB[offset]; }
    uint32 GetStringOffset( const char *name );   //!< Create if non-existent
    uint32 QueryStringOffset( const char *name ); //!< Returns 0 if non-existent
    //@}

    //!\name Low-level parsing methods
    //@{
    const char *ParseItemList( const char *str, const char *str_stop_chars );
    const char *ParseComment( const char *str );
    const char *ParseItem( const char *str );
    const char *ParseKeyword( const char *name, const char *str );
    const char *ParseNumber( const char *name, const char *str );
    const char *ParseString( const char *name, const char *str );
    const char *ParseStringN( const char *name, const char *str );
    const char *ParseComplex( const char *name, const char *str );
    const char *ParseArray( const char *name, const char *str );
    const char *ParseHexBuffer( const char *name, const char *str );
    const char *ParseNiR( const char *name, const char *str );
    //@}

private:
    MemBlock m_DataMB;   //!< Holds actual data
    MemBlock m_StringMB; //!< Holds item identifiers

#ifdef __ENABLE_ITEM_STREAM_REALLOC
    Flags32 m_ReallocFlags;
#endif
    uint32 m_BeginOffset;
    uint32 m_EndOffset;   //!< Offset behind last inserted data in m_DataMB
    uint32 m_StringEndOffset; //!< Offset behind last inserted string in m_StringMB

    uint32 m_NestedComplexOffsetStack[cMaxNestedComplex]; //!< begin_offset of the complex being defined
    int m_NestedComplexLevel;
};

} //util

#endif //UTIL_ITEM_STREAM_H
