#ifndef PLA_UTIL_ARCHETYPE_H
#define PLA_UTIL_ARCHETYPE_H

#include "Config.h"
#include "ItemStream.h"

#define __ENABLE_ARCHETYPE_PROPERTY_GROUP

namespace util {


#ifdef __DRAFTING //\todo Non-intrusive static Archetypes that do NOT change memory footprint/layout of "introspectable" classes
class Something
{
    //\todo Should be optional
    bool Rebuild();
};
// fwd declaration of the archetype class
class GArchetype<Something>; //
template <typename T>
GArchetype
{
    static bool Rebuild( T *p_instance ); //\todo Could be empty
    static void InitArchetype( util::ArchetypeLibrary &al );
    /* \todo Instance Naming can be separated from everything
       else... It's only used in LoadInstance(), and could be done
       externally, REMOVE it from IArchetypeInstance virtual
       API. Named instances are possible, but not mandatory!
    void SetName( const char *name ) {}
    const char *GetName() const { return "geo::np::Stats"; }
    */
};

//template class specialization
class GArchetype<Something>
{
    ..
    ...
}
#endif

class IArchetypeInstance
{
public:
    typedef bool (*notify_touched_property_function_type)( IArchetypeInstance * );

public:
    IArchetypeInstance() {}
    virtual ~IArchetypeInstance() {}
    virtual bool Rebuild() { return true; }
    virtual void SetName( const char *name ) {}
    virtual const char *GetName() const = 0;
    static bool NTPF_Ignore( IArchetypeInstance *p_instance ) { return false; }
    static bool NTPF_Rebuild( IArchetypeInstance *p_instance ) { return true; }
};

/* Archetype is a static property list, which must be explicitly
 initialized. Property archetypes know the actual property value
 offset in class instances, and thus can read/write them.

 \todo ONLY InitArchetype() needs to be written for each archetype,
 the other methods are standard and could be in an SArchetype base
 class, together with archetype property map.

 For implementation purposes, it's a lot easier if the item types allowed are clearly predefined:
 - Fundamental: float,int,uint,bool
 - Simple: Vec,Quat,Flags
 - GProperty_NumberInRange
 - Property_Enumerated
 - Property_FlagMask
 - String

 Everything is a lot easier if archetypes that are not equal to values
 were complex items, instead of simple items. Ex:
 GProperty_NumberInRange<T> could be serialized as a complex with
 sub-items of type <T> "min", "max" and "value". Property_Enumerated
 is an int "value" and two arrays "enum" and "names" of int and
 StringX types. Similarly with Property_FlagMask. This way, and
 Archetype ALWAYS has a "value" sub-item with exactly the same
 type as an Instance property's value.

 ...all in all, Property_XXXX classes are a bit unnecessary, if we can
 store them as Complex in an ItemStream.

 Thus:
 - Archetype property: Complex( prop_name, arch_type ) + Simple( "value", prop_type, value ) +  etc...
 - Instance property: Simple( prop_name, prop_type, value )

 \todo ALL property values MUST have a FIXED maximum size in memory to
       be usable in POD-instance types. Thus, we will NOT ALLOW
       var-sized strings or arrays.
*/
class ArchetypeLibrary
{
public:
    ArchetypeLibrary();
    ~ArchetypeLibrary();

public:
    bool IsEmpty() const { return m_ALIS.IsEmpty(); }

    void BeginArchetype( const char *name );
    void EndArchetype();

    void BeginProperty_Group( const char *name );
    void EndProperty_Group();

    //! Simple property (floatX, intX, uintX, bool, StringN)
    template <typename T>
    void AddProperty( const char *name, const T& value, machine_uint_type offset,
                      IArchetypeInstance::notify_touched_property_function_type ntpf = IArchetypeInstance::NTPF_Rebuild )
        {
            m_ALIS.BeginComplex( name, eType_Property );
            {
                m_ALIS.Write( "value", value );
                m_ALIS.Write( "offset", offset );
                m_ALIS.Write( "size", static_cast<machine_uint_type>( sizeof(T) ) );
                m_ALIS.Write( "ntpf", ntpf );
            }
            m_ALIS.EndComplex();
        }

    //! Complex property: NiR
    template <typename T>
    void AddProperty_NIR( const char *name, const T& value, const T& min, const T& max, machine_uint_type offset,
                          IArchetypeInstance::notify_touched_property_function_type ntpf )//= IArchetypeInstance::NTPF_Rebuild )
        {
            UTIL_ASSERT( min <= max );
            UTIL_LOG_ERROR_IF( value < min || value > max, "NIR<T> Property '%s' with default value out of range [min,max]. Clamping automatically", name );
            m_ALIS.BeginComplex( name, eType_Property );
            {
                m_ALIS.BeginComplex( "value", eType_Property_NIR );
                {
                    m_ALIS.Write( "value", mal::Clamp(value,min,max) );
                    m_ALIS.Write( "min", min );
                    m_ALIS.Write( "max", max );
                }
                m_ALIS.EndComplex();
                m_ALIS.Write( "offset", offset );
                m_ALIS.Write( "size", static_cast<machine_uint_type>( sizeof(T) ) );
                m_ALIS.Write( "ntpf", ntpf );
            }
            m_ALIS.EndComplex();
        }

    //! Complex property Enumerated
    void AddProperty_Enum32( const char *name, uint32 value,
                             unsigned int count, const char** vec_name, const uint32* vec_value,
                             machine_uint_type offset,
                             IArchetypeInstance::notify_touched_property_function_type ntpf );
    //! Complex property Flags
    void AddProperty_Flags32( const char *name, Flags32 value,
                              unsigned int count, const char** vec_name, const int32* vec_value,
                              machine_uint_type offset,
                              IArchetypeInstance::notify_touched_property_function_type ntpf );
    //! \todo Property Array<T>

    //! \name Instance serialization and tweaking
    //@{
    bool LoadInstance( const char *archetype_name, IArchetypeInstance *p_instance, ItemStream::ItemIt instance );
    bool SaveInstance( const char *archetype_name, const IArchetypeInstance *p_instance, ItemStream &is );
    bool ExportInstance( const char *archetype_name, const IArchetypeInstance *p_instance, ItemStream &is );
    bool SyncInstance( const char *archetype_name, IArchetypeInstance *p_instance, ItemStream::ItemItRW instance );
    //@}

private:

#ifdef __ENABLE_ARCHETYPE_PROPERTY_GROUP
    bool LoadInstance_Internal( ItemStream::ItemIt it_archetype, IArchetypeInstance *p_instance, ItemStream::ItemIt it_instance );
    bool SaveInstance_Internal( ItemStream::ItemIt it_archetype, const IArchetypeInstance *p_instance, ItemStream &is );
    void ExportInstance_Internal( ItemStream::ItemIt it_archetype, const IArchetypeInstance *p_instance, ItemStream &is );
    unsigned int SyncInstance_Internal( ItemStream::ItemIt it_archetype, IArchetypeInstance *p_instance, ItemStream::ItemItRW it_instance );
#endif

    bool IsCompatible( ItemStream::ItemIt ap, ItemStream::ItemIt ip ) const;
    void CopyPropertyValueToInstance( ItemStream::ItemIt value,
                                      void *p_instance_data, machine_uint_type offset, machine_uint_type size ) const;
    void CopyInstanceValueToProperty( const void *p_instance_data, machine_uint_type offset, machine_uint_type size,
                                      ItemStream::ItemItRW ip ) const;
    bool AreEqualInstanceAndProperty( const void *p_instance_data,
                                      machine_uint_type offset, machine_uint_type size,
                                      ItemStream::ItemItRW ip ) const;
    bool ValidatePropertyValueInInstance( ItemStream::ItemIt archetype_property,
                                          void *p_instance_data, machine_uint_type offset, machine_uint_type size ) const;

protected:
    ItemStream m_ALIS; //Archetype description
};

#define archetype_offset_of( instance, property ) \
    ( reinterpret_cast<const uint8*>(&instance.property) - reinterpret_cast<const uint8*>(&instance) )

#define __ENABLE_TEST_UTIL_ARCHETYPE
#ifdef __ENABLE_TEST_UTIL_ARCHETYPE
void Run_Test_Archetype();
#endif

} //namespace util

#endif //PLA_UTIL_ARCHETYPE_H
