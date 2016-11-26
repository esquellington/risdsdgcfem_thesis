#include "Archetype.h"
#include "memory.h"
#include "ItemStreamSerialization.h"

//#define __ENABLE_TRACE_ARCHETYPE_SYNC
//#define __ENABLE_TRACE_ARCHETYPE_GROUP

namespace util
{

ArchetypeLibrary::ArchetypeLibrary() : m_ALIS(1<<15,1<<15) {}

ArchetypeLibrary::~ArchetypeLibrary() {}

void ArchetypeLibrary::BeginArchetype( const char *name )
{
    m_ALIS.BeginComplex( name, 0 );//\todo eType_Archetype );
}

void ArchetypeLibrary::EndArchetype()
{
    m_ALIS.EndComplex();
}

void ArchetypeLibrary::BeginProperty_Group( const char *name )
{
#ifdef __ENABLE_ARCHETYPE_PROPERTY_GROUP
    m_ALIS.BeginComplex( name, eType_Property );
    m_ALIS.Write( "offset", machine_uint_type(0) );
    m_ALIS.Write( "size", machine_uint_type(0) );
    m_ALIS.Write( "ntpf", (IArchetypeInstance::notify_touched_property_function_type) IArchetypeInstance::NTPF_Ignore );
    m_ALIS.BeginComplex( "value", eType_Property_Group );
#endif
}

void ArchetypeLibrary::EndProperty_Group()
{
#ifdef __ENABLE_ARCHETYPE_PROPERTY_GROUP
    m_ALIS.EndComplex();
    m_ALIS.EndComplex();
#endif
}


void ArchetypeLibrary::AddProperty_Enum32( const char *name, uint32 value,
                                           unsigned int count, const char** vec_name, const uint32* vec_value,
                                           machine_uint_type offset,
                                           IArchetypeInstance::notify_touched_property_function_type ntpf )
{
    UTIL_ASSERT( count > 0 );
    //\todo static-assert( IsInteger<T> );
    bool bValid(false);
    for( unsigned int i=0; i<count && !bValid; i++ ) bValid = value == vec_value[i];
    UTIL_LOG_ERROR_IF( !bValid, "Enum32 Property '%s' with invalid enumerated value %d. Defaulting to %d", name, value, vec_value[0] );
    m_ALIS.BeginComplex( name, eType_Property );
    {
        m_ALIS.BeginComplex( "value", eType_Property_Enum );
        {
            m_ALIS.Write( "value", (bValid) ? value : vec_value[0] );
            uint32 *array_value = m_ALIS.AllocArray<uint32>( "vec_value", count );
            String32 *array_name = m_ALIS.AllocArray<String32>( "vec_name", count );
            for( unsigned int i=0; i<count; i++ )
            {
                array_value[i] = vec_value[i];
                array_name[i].Set( vec_name[i] );
            }
        }
        m_ALIS.EndComplex();
        m_ALIS.Write( "offset", offset );
        m_ALIS.Write( "size", static_cast<machine_uint_type>( sizeof(uint32) ) );
        m_ALIS.Write( "ntpf", ntpf );
    }
    m_ALIS.EndComplex();
}

//! Complex property Flags
void ArchetypeLibrary::AddProperty_Flags32( const char *name, Flags32 value,
                                            unsigned int count, const char** vec_name, const int32* vec_value,
                                            machine_uint_type offset,
                                            IArchetypeInstance::notify_touched_property_function_type ntpf )
{
    UTIL_ASSERT( count > 0 );
    Flags32 tmp_value( value );
    for( unsigned int i=0; i<count; i++ ) tmp_value.Disable( vec_value[i] );
    bool bValid( 0 == int32(tmp_value) ); //Valid if no bits set after disabling all valid flags
    UTIL_LOG_ERROR_IF( !bValid, "Flags32 Property '%s' with invalid value %x (unknown flags %x). Defaulting to %x",
                       name, int32(value), int32(tmp_value), vec_value[0] );
    m_ALIS.BeginComplex( name, eType_Property );
    {
        m_ALIS.BeginComplex( "value", eType_Property_Flags );
        {
            m_ALIS.Write( "value", (bValid) ? value : Flags32(vec_value[0]) );
            int32 *array_value = m_ALIS.AllocArray<int32>( "vec_value", count );
            String32 *array_name = m_ALIS.AllocArray<String32>( "vec_name", count );
            for( unsigned int i=0; i<count; i++ )
            {
                array_value[i] = vec_value[i];
                array_name[i].Set( vec_name[i] );
            }
        }
        m_ALIS.EndComplex();
        m_ALIS.Write( "offset", offset );
        m_ALIS.Write( "size", static_cast<machine_uint_type>( sizeof(Flags32) ) );
        m_ALIS.Write( "ntpf", ntpf );
    }
    m_ALIS.EndComplex();
}

#ifdef __ENABLE_ARCHETYPE_PROPERTY_GROUP //----------------------------------------------------------------
bool ArchetypeLibrary::LoadInstance_Internal( ItemStream::ItemIt it_archetype, IArchetypeInstance *p_instance, ItemStream::ItemIt it_instance )
{
    /* Use the Archetype properties descriptions (and offsets)
       to initialize p_instance Values read from
       instance. Any Value not present in instance will
       have the default value from the Archetype.
    */
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Copy from instance if exists and matches type...or set to archetype default otherwise
        const char *name = it_ap.GetName();
        ItemStream::ItemIt it_ap_value( it_ap.GetSubItem().Find("value") );
        bool bIsPropertyGroup( eType_Property_Group == it_ap_value.GetType() );
        if( bIsPropertyGroup )
        {
            ItemStream::ItemIt it_ip = it_instance.Find(name);
#ifdef __ENABLE_TRACE_ARCHETYPE_GROUP
            UTIL_LOG("ArchetypeLibrary::LoadInstance() Group %s", name );
#endif
            if( it_ip.IsValid() && eType_Property_Group == it_ip.GetType() ) LoadInstance_Internal( it_ap_value, p_instance, it_ip.GetSubItem() );
            else
            {
#ifdef __ENABLE_TRACE_ARCHETYPE_GROUP
                UTIL_LOG("ArchetypeLibrary::LoadInstance() Group %s not found, setting default", name );
#endif
                /* Call with invalid it_instance, so that all basic
                   props get copied recursively (no it_instance.Find()
                   will return valid it_ip) */
                LoadInstance_Internal( it_ap_value, p_instance, ItemStream::ItemIt() );
            }
        }
        else
        {
            ItemStream::ItemIt it_ip = it_instance.Find(name);
            machine_uint_type offset( it_ap.GetSubItem().Find("offset").Get<machine_uint_type>() );
            machine_uint_type size( it_ap.GetSubItem().Find("size").Get<machine_uint_type>() );
            if( IsCompatible( it_ap_value, it_ip ) ) //fails if !it_ip.IsValid()
            {
                // Set instance value
                CopyPropertyValueToInstance( it_ip, p_instance, offset, size ); //it_ip == value
            }
            else
            {
                // Set archetype default value
                CopyPropertyValueToInstance( it_ap_value, p_instance, offset, size );
                if( it_ip.IsValid() ) UTIL_LOG_WARNING( "ArchetypeLibrary::LoadInstance() Incompatible type in property %s, using default", name );
            }
            // After copying instance value, we validate according to archetype (min,max,enum...)"
            bool bValid = ValidatePropertyValueInInstance( it_ap, p_instance, offset, size );
            //UTIL_LOG_ERROR_IF( !bValid, "ArchetypeLibrary: Invalid value in property %s", name );
        }
    }
    return true;
}

bool ArchetypeLibrary::LoadInstance( const char *archetype_name, IArchetypeInstance *p_instance, ItemStream::ItemIt instance )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::LoadInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    ItemStream::ItemIt first_ip( instance.GetSubItem() );
    LoadInstance_Internal( it_archetype, p_instance, first_ip );
    p_instance->SetName( instance.GetName() ); //\todo THIS COULD BE DONE OUTSIDE and avoid requiring GetName() method
    p_instance->Rebuild();
    return true;
}


bool ArchetypeLibrary::SaveInstance_Internal( ItemStream::ItemIt it_archetype, const IArchetypeInstance *p_instance, ItemStream &is )
{
    /* Use the archetype to iterate over properties, retrieve their
       p_instance values using the offsets and write them to the
       output is.
    */
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Write archetype property value to output is
        ItemStream::ItemIt it_ap_value( it_ap.GetSubItem().Find("value") );
        bool bIsPropertyGroup( eType_Property_Group == it_ap_value.GetType() );
        if( bIsPropertyGroup )
        {
#ifdef __ENABLE_TRACE_ARCHETYPE_GROUP
            UTIL_LOG("ArchetypeLibrary::SaveInstance() Group %s", it_ap.GetName() );
#endif
            is.BeginComplex( it_ap.GetName(), eType_Property_Group );
                SaveInstance_Internal( it_ap_value, p_instance, is );
            is.EndComplex();
        }
        else
        {
            ItemStream::ItemItRW it_ip_value;
            if( it_ap_value.IsSimple() )
                it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value );
            else if( it_ap_value.IsComplex() )
                it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value.GetSubItem().Find("value") );
            else
            {
                UTIL_LOG_ERROR("ArchetypeLibrary::SaveInstance() Non-Simple or Complex properties not supported yet");
                return false;
            }
            // Overwrite output property value with instance value
            CopyInstanceValueToProperty( p_instance,
                                         it_ap.GetSubItem().Find("offset").Get<machine_uint_type>(),
                                         it_ap.GetSubItem().Find("size").Get<machine_uint_type>(),
                                         it_ip_value );
        }
    }
    return true;
}

bool ArchetypeLibrary::SaveInstance( const char *archetype_name, const IArchetypeInstance *p_instance, ItemStream &is )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::SaveInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    //\todo NAME COULD BE PASSED AS A PARAMETER and avoid requiring GetName() method
    is.BeginComplex( p_instance->GetName(), 0 ); //eType_Instance );
        SaveInstance_Internal( it_archetype, p_instance, is );
    is.EndComplex();
    return true;
}

void ArchetypeLibrary::ExportInstance_Internal( ItemStream::ItemIt it_archetype, const IArchetypeInstance *p_instance, ItemStream &is )
{
    /* Use the archetype to iterate over properties and copy their "value" entry to the output is.
       \note Archetype bookkeeping info such as "offset" and "size" is not exported.
    */
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Write archetype property value to output is
        ItemStream::ItemIt it_ap_value( it_ap.GetSubItem().Find("value") );
        bool bIsPropertyGroup( eType_Property_Group == it_ap_value.GetType() );
        if( bIsPropertyGroup )
        {
#ifdef __ENABLE_TRACE_ARCHETYPE_GROUP
            UTIL_LOG("ArchetypeLibrary::SaveInstance() Group %s", it_ap.GetName() );
#endif
            is.BeginComplex( it_ap.GetName(), eType_Property_Group );
                ExportInstance_Internal( it_ap_value, p_instance, is );
            is.EndComplex();
        }
        else
        {
            ItemStream::ItemItRW it_ip_value;
            if( it_ap_value.IsSimple() )
                it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value ); //Simple value def
            else if( it_ap_value.IsComplex() )
                it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value ); //Complex value def (NIR, Enum...)
            else
            {
                UTIL_LOG_ERROR("ArchetypeLibrary::ExportInstance() Non-Simple or Complex properties not supported yet");
                return;
            }
            /*TEMP: Debugging, remove when fixed...
            UTIL_LOG("CopyInstanceValueToProperty() AP"); std::cout << it_ap << std::endl;
            UTIL_LOG("CopyInstanceValueToProperty() IP"); std::cout << it_ip_value << std::endl;
            */
            // Overwrite output property value with instance value
            CopyInstanceValueToProperty( p_instance,
                                         it_ap.GetSubItem().Find("offset").Get<machine_uint_type>(),
                                         it_ap.GetSubItem().Find("size").Get<machine_uint_type>(),
                                         it_ip_value );
        }
    }
}

bool ArchetypeLibrary::ExportInstance( const char *archetype_name, const IArchetypeInstance *p_instance, ItemStream &is )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::ExportInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    is.BeginComplex( p_instance->GetName(), eType_Property_Group ); //eType_Instance );
        ExportInstance_Internal( it_archetype, p_instance, is );
    is.EndComplex();
    return true;
}

unsigned int ArchetypeLibrary::SyncInstance_Internal( ItemStream::ItemIt it_archetype, IArchetypeInstance *p_instance, ItemStream::ItemItRW it_instance )
{
    unsigned int num_touched(0);
    /* Use the archetype to iterate over properties, find them in the
       instance and copy the value one way or another depending on
       which one changed, calling NTPF if required and available.
    */
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Copy from instance if exists and matches type...or set to archetype default otherwise
        const char *name = it_ap.GetName();
        ItemStream::ItemIt it_ap_value( it_ap.GetSubItem().Find("value") );
        ItemStream::ItemItRW it_ip = it_instance.Find(name);

#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
        if( !it_ip.IsValid() )
        {
            //TEMP: This should never happen if Export/Sync are properly matched!!
            UTIL_LOG_ERROR("ArchetypeLibrary::SyncInstance() unmatched property %s", name );
            UTIL_LOG("Archetype = "); std::cout << it_archetype << std::endl;
            UTIL_LOG("Instance = "); std::cout << instance << std::endl;
        }
#endif

        bool bIsPropertyGroup( eType_Property_Group == it_ap_value.GetType() );
        if( bIsPropertyGroup )
        {
#ifdef __ENABLE_TRACE_ARCHETYPE_GROUP
            UTIL_LOG("ArchetypeLibrary::SyncInstance() Group %s", name );
#endif
            //TEMP: This should never raise if Export/Sync are properly matched!!
            UTIL_LOG_ASSERT( eType_Property_Group == it_ip.GetType(), "ArchetypeLibrary::SyncInstance() unmatched Group type %s", name );
            num_touched += SyncInstance_Internal( it_ap_value, p_instance, it_ip.GetSubItem() );
        }
        else
        {
            machine_uint_type offset( it_ap.GetSubItem().Find("offset").Get<machine_uint_type>() );
            machine_uint_type size( it_ap.GetSubItem().Find("size").Get<machine_uint_type>() );
            if( it_ip.IsTouched() )
            {
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
                UTIL_LOG( "Archetype::SyncInstance() touched property %s", name );
#endif
                CopyPropertyValueToInstance( it_ip, p_instance, offset, size );
                it_ip.Untouch();
                // notify single property touched
                ItemStream::ItemIt it_ntpf( it_ap.GetSubItem().Find("ntpf") );
                UTIL_ASSERT( it_ntpf.IsValid() );
                if( it_ntpf.Get<IArchetypeInstance::notify_touched_property_function_type>()( p_instance ) )
                    num_touched++;
            }
            else if( !AreEqualInstanceAndProperty( p_instance, offset, size, it_ip ) )
            {
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
                UTIL_LOG( "Archetype::SyncInstance() touched instance %s", name );
#endif
                CopyInstanceValueToProperty( p_instance, offset, size, it_ip );
                it_ip.Touch();
            }
            /*
            // After copying instance value, we validate according to archetype (min,max,enum...)"
            bool bValid = ValidatePropertyValueInInstance( it_ap, p_instance, offset, size );
            UTIL_LOG_ERROR_IF( !bValid, "ArchetypeLibrary::SyncInstance() Invalid value in property %s", name );
            */
            //\todo REDUNDANT: This should be unnecessary if any property change respects its archetype limits, which MUST HAPPEN, otherwise it's an error
            UTIL_LOG_ASSERT( ValidatePropertyValueInInstance( it_ap, p_instance, offset, size ),
                             "ArchetypeLibrary::SyncInstance() Invalid value in property %s",
                             name );
        }
    }
    return num_touched;
}

bool ArchetypeLibrary::SyncInstance( const char *archetype_name, IArchetypeInstance *p_instance, ItemStream::ItemItRW instance )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::SyncInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    ItemStream::ItemItRW first_ip( instance.GetSubItem() );
    unsigned int num_touched = SyncInstance_Internal( it_archetype, p_instance, first_ip );
    // Rebuild instance if touched
    if( num_touched > 0 )
    {
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
        UTIL_LOG( "ArchetypeLibrary::SyncInstance() Rebuilding... (num_touched = %d)", num_touched );
#endif
        p_instance->Rebuild();
    }
    return true;
}

#else //__ENABLE_ARCHETYPE_PROPERTY_GROUP ----------------------------------------------------------------

bool ArchetypeLibrary::LoadInstance( const char *archetype_name, IArchetypeInstance *p_instance, ItemStream::ItemIt instance )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::LoadInstance() Unknown Archetype %s", archetype_name );
        return false;
    }

    /* Use the Archetype properties descriptions (and offsets)
       to initialize p_instance Values read from
       instance. Any Value not present in instance will
       have the default value from the Archetype.
    */
    ItemStream::ItemIt first_ip( instance.GetSubItem() );
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Copy from instance if exists and matches type...or set to archetype default otherwise
        const char *name = it_ap.GetName();
        ItemStream::ItemIt it_ip = first_ip.Find(name);
        machine_uint_type offset( it_ap.GetSubItem().Find("offset").Get<machine_uint_type>() );
        machine_uint_type size( it_ap.GetSubItem().Find("size").Get<machine_uint_type>() );
        if( IsCompatible( it_ap.GetSubItem().Find("value"), it_ip ) )
        {
            // Set instance value
            CopyPropertyValueToInstance( it_ip, p_instance, offset, size ); //it_ip == value
        }
        else
        {
            // Set archetype default value
            CopyPropertyValueToInstance( it_ap.GetSubItem().Find("value"), p_instance, offset, size );
            if( it_ip.IsValid() ) UTIL_LOG_WARNING( "ArchetypeLibrary::LoadInstance() Incompatible type in property %s, using default", name );
        }
        // After copying instance value, we validate according to archetype (min,max,enum...)"
        bool bValid = ValidatePropertyValueInInstance( it_ap, p_instance, offset, size );
        //UTIL_LOG_ERROR_IF( !bValid, "ArchetypeLibrary: Invalid value in property %s", name );
    }
    p_instance->SetName( instance.GetName() );
    p_instance->Rebuild();
    return true;
}

bool ArchetypeLibrary::SaveInstance( const char *archetype_name, const IArchetypeInstance *p_instance, ItemStream &is )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::SaveInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    /* Use the archetype to iterate over properties, retrieve
       their p_instance using the offsets and write them to
       the output is.
    */
    is.BeginComplex( p_instance->GetName(), 0 ); //eType_Instance );
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Write archetype property value to output is
        ItemStream::ItemIt it_ap_value( it_ap.GetSubItem().Find("value") );
        ItemStream::ItemItRW it_ip_value;
        if( it_ap_value.IsSimple() )
            it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value );
        else if( it_ap_value.IsComplex() )
            it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value.GetSubItem().Find("value") );
        else
            UTIL_LOG_ERROR("ArchetypeLibrary::SaveInstance() Non-Simple or Complex properties not supported yet");
        // Overwrite output property value with instance value
        CopyInstanceValueToProperty( p_instance,
                                     it_ap.GetSubItem().Find("offset").Get<machine_uint_type>(),
                                     it_ap.GetSubItem().Find("size").Get<machine_uint_type>(),
                                     it_ip_value );
    }
    is.EndComplex();
    return true;
}

bool ArchetypeLibrary::ExportInstance( const char *archetype_name, const IArchetypeInstance *p_instance, ItemStream &is )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::ExportInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    /* Use the archetype to iterate over properties and copy their "value" entry to the output is.
       \note Archetype bookkeeping info such as "offset" and "size" is not exported.
    */
    is.BeginComplex( p_instance->GetName(), eType_Property_Group ); //eType_Instance );
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Write archetype property value to output is
        ItemStream::ItemIt it_ap_value( it_ap.GetSubItem().Find("value") );
        ItemStream::ItemItRW it_ip_value;
        if( it_ap_value.IsSimple() )
            it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value ); //Simple value def
        else if( it_ap_value.IsComplex() )
            it_ip_value = is.WriteItem( it_ap.GetName(), it_ap_value ); //Complex value def (NIR, Enum...)
        else
            UTIL_LOG_ERROR("ArchetypeLibrary::ExportInstance() Non-Simple or Complex properties not supported yet");
        // Overwrite output property value with instance value
        CopyInstanceValueToProperty( p_instance,
                                     it_ap.GetSubItem().Find("offset").Get<machine_uint_type>(),
                                     it_ap.GetSubItem().Find("size").Get<machine_uint_type>(),
                                     it_ip_value );
    }
    is.EndComplex();
    return true;
}

bool ArchetypeLibrary::SyncInstance( const char *archetype_name, IArchetypeInstance *p_instance, ItemStream::ItemItRW instance )
{
    ItemStream::ItemIt it_archetype( m_ALIS.Begin().Find(archetype_name) );
    if( !it_archetype.IsValid() )
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::SyncInstance() Unknown Archetype %s", archetype_name );
        return false;
    }
    /* Use the archetype to iterate over properties, find them in the
       instance and copy the value one way or another depending on
       which one changed, calling NTPF if required and available.
    */
    unsigned int num_touched(0);
    ItemStream::ItemItRW first_ip( instance.GetSubItem() );
    for( ItemStream::ItemIt it_ap = it_archetype.GetSubItem(); it_ap.IsValid(); ++it_ap )
    {
        // Copy from instance if exists and matches type...or set to archetype default otherwise
        const char *name = it_ap.GetName();
        ItemStream::ItemItRW it_ip = first_ip.Find(name);
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
        if( !it_ip.IsValid() )
        {
            UTIL_LOG_ERROR("Archetype::SyncInstance() unmatched property %s", name );
            UTIL_LOG("Archetype = "); std::cout << it_archetype << std::endl;
            UTIL_LOG("Instance = "); std::cout << instance << std::endl;
        }
#endif
        machine_uint_type offset( it_ap.GetSubItem().Find("offset").Get<machine_uint_type>() );
        machine_uint_type size( it_ap.GetSubItem().Find("size").Get<machine_uint_type>() );
        if( it_ip.IsTouched() )
        {
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
            UTIL_LOG( "Archetype::SyncInstance() touched property %s", name );
#endif
            CopyPropertyValueToInstance( it_ip, p_instance, offset, size );
            it_ip.Untouch();
            // notify single property touched
            ItemStream::ItemIt it_ntpf( it_ap.GetSubItem().Find("ntpf") );
            UTIL_ASSERT( it_ntpf.IsValid() );
            if ( it_ntpf.Get<IArchetypeInstance::notify_touched_property_function_type>()( p_instance ) )
                num_touched++;
        }
        else if( !AreEqualInstanceAndProperty( p_instance, offset, size, it_ip ) )
        {
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
            UTIL_LOG( "Archetype::SyncInstance() touched instance %s", name );
#endif
            CopyInstanceValueToProperty( p_instance, offset, size, it_ip );
            it_ip.Touch();
        }
        /*
        // After copying instance value, we validate according to archetype (min,max,enum...)"
        bool bValid = ValidatePropertyValueInInstance( it_ap, p_instance, offset, size );
        UTIL_LOG_ERROR_IF( !bValid, "ArchetypeLibrary::SyncInstance() Invalid value in property %s", name );
        */
        //\todo REDUNDANT: This should be unnecessary if any property change respects its archetype limits, which MUST HAPPEN, otherwise it's an error
        UTIL_LOG_ASSERT( ValidatePropertyValueInInstance( it_ap, p_instance, offset, size ),
                         "ArchetypeLibrary::SyncInstance() Invalid value in property %s",
                         name );
    }
    // Rebuild instance if touched
    if( num_touched > 0 )
    {
#ifdef __ENABLE_TRACE_ARCHETYPE_SYNC
        UTIL_LOG_WARNING( "ArchetypeLibrary::SyncInstance() Rebuilding... (num_touched = %d)", num_touched );
#endif
        p_instance->Rebuild();
    }
    return true;
}
#endif //__ENABLE_ARCHETYPE_PROPERTY_GROUP ----------------------------------------------------------------

/* Check compatibility between an archetype_property (ap) and an instance_property (ip)
   Compatible properties include:
   - Same full type
   - GProperty_NumberInRange<T> and T
   - Property_Enumerated and enumerated_value_type (int32?)
   - Property_FlagMask and flagmask_value_type (int32? or Flags32...)
   - Same-size arrays of fundamental type T
   - Shorter or equal Strings
*/
bool ArchetypeLibrary::IsCompatible( ItemStream::ItemIt av, ItemStream::ItemIt ip ) const
{
    if( !ip.IsValid() ) return false;
    return ( av.GetFullType() == ip.GetFullType()
             // NiR<T> and <T>
             ||
             ( av.GetType() == eType_Property_NIR
               && av.GetSubItem().Find("value").GetFullType() == ip.GetFullType() )
             // Enumerated \todo Consider supporting String32 value_name, too...
             ||
             ( av.GetType() == eType_Property_Enum
               && av.GetSubItem().Find("value").GetType() == pla_type_id<uint32>::value )
             // Flags32
             ||
             ( av.GetType() == eType_Property_Flags
               && av.GetSubItem().Find("value").GetType() == pla_type_id<Flags32>::value )
             // Same size string ( //\todo accept GString<N> types in ip too)
             ||
             ( av.GetType() == pla_type_id<char>::value
               && av.GetType() == ip.GetType()
               && av.IsArray() && ip.IsArray() )
             // Same size array
             ||
             ( av.GetType() == ip.GetType()
               && av.IsArray() && ip.IsArray()
               && av.GetArrayCount() == ip.GetArrayCount() ) );
}

void ArchetypeLibrary::CopyPropertyValueToInstance( ItemStream::ItemIt value,
                                                    void *p_instance_data,
                                                    machine_uint_type offset, machine_uint_type size ) const
{
    UTIL_ASSERT( value.IsValid() );
    if( value.IsSimple() ) // Simple property: floatX, intX, uintX, FlagsX, bool, VecN, QuatX)
        memcpy( reinterpret_cast<uint8*>(p_instance_data)+offset,
                reinterpret_cast<const uint8*>( value.GetDataPtr() ),
                size );
    else if( value.IsComplex() ) // Complex property: NiR, Enumerated, FlagMask
        memcpy( reinterpret_cast<uint8*>(p_instance_data)+offset,
                reinterpret_cast<const uint8*>( value.GetSubItem().Find("value").GetDataPtr() ), //Only copy actual value
                size );
    else //\todo if( value.IsArray() ) //support fixed-size arrays
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::CopyPropertyValueToInstance() cannot copy non-Simple or non-Complex Property to Instance value" );
    }
}

bool ArchetypeLibrary::AreEqualInstanceAndProperty( const void *p_instance_data,
                                                    machine_uint_type offset, machine_uint_type size,
                                                    ItemStream::ItemItRW ip ) const
{
    UTIL_ASSERT( ip.IsValid() );
    if( ip.IsSimple() )
        return 0 == memcmp( reinterpret_cast<const uint8*>( ip.GetDataPtr() ),
                            reinterpret_cast<const uint8*>(p_instance_data)+offset,
                            size );
    else if( ip.IsComplex() )
        return 0 == memcmp( reinterpret_cast<const uint8*>( ip.GetSubItem().Find("value").GetDataPtr() ),
                            reinterpret_cast<const uint8*>(p_instance_data)+offset,
                            size );
    else //\todo if( value.IsArray() ) //support fixed-size arrays
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::CopyInstanceValueToProperty() cannot copy Instance value to non-Simple or non-Complex Property" );
        return true;
    }
}

void ArchetypeLibrary::CopyInstanceValueToProperty( const void *p_instance_data,
                                                    machine_uint_type offset, machine_uint_type size,
                                                    ItemStream::ItemItRW ip ) const
{
    UTIL_ASSERT( ip.IsValid() );
    if( ip.IsSimple() )
        memcpy( reinterpret_cast<uint8*>( ip.GetDataPtr() ),
                reinterpret_cast<const uint8*>(p_instance_data)+offset,
                size );
    else if( ip.IsComplex() )
    {
        memcpy( reinterpret_cast<uint8*>( ip.GetSubItem().Find("value").GetDataPtr() ), //Copy on actual value
                reinterpret_cast<const uint8*>(p_instance_data)+offset,
                size );
    }
    else //\todo if( value.IsArray() ) //support fixed-size arrays
    {
        UTIL_LOG_ERROR( "ArchetypeLibrary::CopyInstanceValueToProperty() cannot copy Instance value to non-Simple or non-Complex Property" );
    }
}

bool ArchetypeLibrary::ValidatePropertyValueInInstance( ItemStream::ItemIt archetype_property,
                                                        void *p_instance_data, machine_uint_type offset, machine_uint_type size ) const
{
    UTIL_ASSERT( archetype_property.IsValid() );
    ItemStream::ItemIt av( archetype_property.GetSubItem().Find("value") );
    // Complex property: NiR, Enumerated, FlagMask
    if( av.IsComplex() )
    {
        // NiR<T>
        if( av.GetType() == eType_Property_NIR )
        {
            bool bInvalid(false);
            switch( av.GetSubItem().Find("value").GetType() )
            {
            case eType_Int32:
                {
                    int32 *ptr( reinterpret_cast<int32*>( reinterpret_cast<uint8*>(p_instance_data)+offset ) );
                    int32 min( av.GetSubItem().Find("min").Get<int32>() );
                    int32 max( av.GetSubItem().Find("max").Get<int32>() );
                    bInvalid = ( *ptr < min || *ptr > max );
                    UTIL_LOG_ERROR_IF( bInvalid,
                                       "ArchetypeLibrary: Invalid value %d in property %s [%d..%d]",
                                       *ptr, archetype_property.GetName(), min, max );
                    if( bInvalid ) *ptr = mal::Clamp(*ptr,min,max);
                }
                break;
            case eType_UInt32:
                {
                    uint32 *ptr( reinterpret_cast<uint32*>( reinterpret_cast<uint8*>(p_instance_data)+offset ) );
                    uint32 min( av.GetSubItem().Find("min").Get<uint32>() );
                    uint32 max( av.GetSubItem().Find("max").Get<uint32>() );
                    bInvalid = ( *ptr < min || *ptr > max );
                    UTIL_LOG_ERROR_IF( bInvalid,
                                       "ArchetypeLibrary: Invalid value %u in property %s [%u..%u]",
                                       *ptr, archetype_property.GetName(), min, max );
                    if( bInvalid ) *ptr = mal::Clamp(*ptr,min,max);
                }
                break;
            case eType_Float32:
                {
                    float32 *ptr( reinterpret_cast<float32*>( reinterpret_cast<uint8*>(p_instance_data)+offset ) );
                    float32 min( av.GetSubItem().Find("min").Get<float32>() );
                    float32 max( av.GetSubItem().Find("max").Get<float32>() );
                    bInvalid = ( *ptr < min || *ptr > max );
                    UTIL_LOG_ERROR_IF( bInvalid,
                                       "ArchetypeLibrary: Invalid value %f in property %s [%f..%f]",
                                       *ptr, archetype_property.GetName(), min, max );
                    if( bInvalid ) *ptr = mal::Clamp(*ptr,min,max);
                }
                break;
            case eType_Float64:
                {
                    float64 *ptr( reinterpret_cast<float64*>( reinterpret_cast<uint8*>(p_instance_data)+offset ) );
                    float64 min( av.GetSubItem().Find("min").Get<float64>() );
                    float64 max( av.GetSubItem().Find("max").Get<float64>() );
                    bInvalid = ( *ptr < min || *ptr > max );
                    UTIL_LOG_ERROR_IF( bInvalid,
                                       "ArchetypeLibrary: Invalid value %lf in property %s [%lf..%lf]",
                                       *ptr, archetype_property.GetName(), min, max );
                    if( bInvalid ) *ptr = mal::Clamp(*ptr,min,max);
                }
                break;
            default:
                UTIL_LOG_WARNING("ArchetypeLibrary: Cannot validate eType_Property_NIR type %d", av.GetSubItem().Find("value").GetType() );
                break;
            }
            return !bInvalid;
        }
        else if( av.GetType() == eType_Property_Enum ) //Enumerated
        {
            uint32 *ptr( reinterpret_cast<uint32*>( reinterpret_cast<uint8*>(p_instance_data)+offset ) );
            uint32 count( av.GetSubItem().Find("vec_value").GetArrayCount() );
            const uint32 *array_value( av.GetSubItem().Find("vec_value").GetArrayPtr<uint32>() );
            bool bValid(false);
            for( unsigned int i=0; i<count && !bValid; i++ ) bValid = *ptr == array_value[i];
            UTIL_LOG_ERROR_IF( !bValid,
                               "ArchetypeLibrary: Invalid value %d in enumerated property %s",
                               *ptr, archetype_property.GetName() );
            if( !bValid ) *ptr = av.GetSubItem().Find("value").Get<uint32>();
            return bValid;
        }
        else if( av.GetType() == eType_Property_Flags ) //Flags
        {
            Flags32 *ptr( reinterpret_cast<Flags32*>( reinterpret_cast<uint8*>(p_instance_data)+offset ) );
            uint32 count( av.GetSubItem().Find("vec_value").GetArrayCount() );
            const int32 *array_value( av.GetSubItem().Find("vec_value").GetArrayPtr<int32>() );

            Flags32 tmp_value( *ptr );
            for( unsigned int i=0; i<count; i++ ) tmp_value.Disable( array_value[i] );
            bool bValid( 0 == int32(tmp_value) ); //Valid if no bits set after disabling all valid flags
            UTIL_LOG_ERROR_IF( !bValid,
                               "ArchetypeLibrary: Invalid value %x in Flags property %s (unknown flags %x).",
                               int32(*ptr), archetype_property.GetName(), int32(tmp_value) );
            if( !bValid ) *ptr = av.GetSubItem().Find("value").Get<Flags32>();
            return bValid;
        }
        else
        {
            UTIL_LOG_WARNING("ArchetypeLibrary: Cannot validate property type %d", av.GetType() );
            return true;
        }
    }
    else // Simple properties have no restritions
        return true;
}

} //namespace util


#ifdef __ENABLE_TEST_UTIL_ARCHETYPE
namespace util
{

class TestArchetype: public IArchetypeInstance
{
public:
    TestArchetype() : m_Int32(0), m_UInt32(0), m_Float32(0) {}
    ~TestArchetype() {}

    //!\name IArchetypeInstance implementation
    //@{
    bool Rebuild() { return true; }
    void SetName( const char *name ) { m_Name.Set(name); }
    const char *GetName() const { return m_Name; }
    //@}

public:
    static void InitArchetype( ArchetypeLibrary &al )
        {
            al.BeginArchetype( "TestArchetype" );
            {
                TestArchetype ti;
                al.AddProperty_NIR<int32>( "int32", 1, -10, 10, archetype_offset_of(ti,m_Int32), NTPF_Rebuild );
                al.AddProperty<uint32>( "uint32", uint32(21232), archetype_offset_of(ti,m_UInt32), NTPF_Rebuild );
                al.AddProperty_NIR<float>( "float32", 0.1f, 0.0f, 1.0f, archetype_offset_of(ti,m_Float32), NTPF_Rebuild );
                al.AddProperty<String32>( "string32", String32("jajaja"), archetype_offset_of(ti,m_String32), NTPF_Rebuild );

                enum EEnum { eEnumValue0 = 0, eEnumValue1 = 10, eEnumValue2 = 20, cNumEnumValues = 3 };
                const char *vec_names[] = { "v0", "v1", "v2" };
                uint32 vec_values[] = { (uint32)eEnumValue0, (uint32)eEnumValue1, (uint32)eEnumValue2 };
                al.AddProperty_Enum32( "enum32",
                                       (uint32)eEnumValue1,
                                       (uint32)cNumEnumValues, vec_names, vec_values, //\todo Cannot use { "v1", "v2" }, and {v1,v2} here, could in C++11
                                       archetype_offset_of(ti,m_Enum32), NTPF_Rebuild );
            }
            al.EndArchetype();
        }
public:
    uint32 m_Int32;
    uint32 m_UInt32;
    float32 m_Float32;
    String32 m_String32;
    uint32 m_Enum32;

    String32 m_Name;
};

void Run_Test_Archetype()
{
    // Init Archetype
    ArchetypeLibrary al;
    TestArchetype::InitArchetype(al);

    // Fill instance stream
    ItemStream is( 1024, 1024 );
    is.BeginComplex( "instance1", 0 );
    {
        is.Write( "int32", int32(-1) );
        is.Write( "uint32", uint32(11) );
        //is.Write( "float32", float32(111.0f) );
        is.Write( "enum32", uint32(1) );
    }
    is.EndComplex();
    is.BeginComplex( "instance3", 0 );
    {
        is.Write( "int32", int32(-2) );
        //is.Write( "uint32", uint32(22) );
        is.Write( "float32", float32(222.0f) );
        is.Write( "enum32", uint32(20) );
    }
    is.EndComplex();
    is.BeginComplex( "instance2", 0 );
    {
        //is.Write( "int32", int32(-3) );
        is.Write( "uint32", uint32(33) );
        is.Write( "float32", float32(333.0f) );
        is.Write( "string32", String32("kkkkkkkk") );
        is.Write( "enum32", uint32(666) );
    }
    is.EndComplex();

    // Save instance stream
    is.SaveTxt( "util_Archetype_Run_Test_Archetype_INPUT.txt" );

    // Load instance stream
    is.Clear();
    is.LoadTxt( "util_Archetype_Run_Test_Archetype_INPUT.txt" );

    //TEMPORAL: std::cout << is << std::endl;

    // Load instances
    //\todo Library could internally associate <ti> type to "archetype_name" to avoid having to specify it on instance load/save.
    TestArchetype vec_ti[3];
    unsigned int count(0);
    for( ItemStream::ItemIt it_instance=is.Begin(); it_instance.IsValid(); ++it_instance )
        al.LoadInstance( "TestArchetype", &vec_ti[count++], it_instance );

    // Save instances
    is.Clear();
    for( unsigned int i=0; i<count; i++ )
        al.SaveInstance( "TestArchetype", &vec_ti[i], is );
    is.SaveTxt( "util_Archetype_Run_Test_Archetype_OUTPUT.txt" );
    //UTIL_LOG_WARNING( "Instance values = [%d, %d, %f]", ti.m_Int32, ti.m_UInt32, ti.m_Float32 );
}

} //namespace util

#endif //__ENABLE_TEST_UTIL_ARCHETYPE
