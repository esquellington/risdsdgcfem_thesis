#include <util/ItemStream.h>
#include <util/StringUtils.h>
#include <string.h>
#include <stdio.h> //TEMPORAL: For fopen/fread/fwrite
#include <util/ItemStreamSerialization.h>

namespace util {


//----------------------------------------------------------------
// ItemStream implementation
//----------------------------------------------------------------

void ItemStream::WriteString( int32 id, const char* str, int length )
{
    UTIL_ASSERT( id > 0 );
    unsigned int l(0);
    if( length < 0 ) while( str[l] != 0 ) l++;
    else l = length;
    char *ptr = AllocArray<char>( id, l+1 ); //l includes '\0' ending mark
    for( unsigned i=0; i<l; i++ ) ptr[i] = str[i];
    ptr[l] = '\0';
}

void ItemStream::WriteString( const char *name, const char* str, int length )
{
    unsigned int l(0);
    if( length < 0 ) while( str[l] != 0 ) l++;
    else l = length;
    char *ptr = AllocArray<char>( name, l+1 ); //l includes '\0' ending mark
    for( unsigned i=0; i<l; i++ ) ptr[i] = str[i];
    ptr[l] = '\0';
}

void ItemStream::BeginComplex( int32 id, uint16 type )
{
    UTIL_ASSERT( id > 0 );
    UTIL_ASSERT( m_NestedComplexLevel < cMaxNestedComplex );
    BeginComplex_Internal( id, type );
}

void ItemStream::BeginComplex( const char *name, uint16 type )
{
    uint32 id = GetStringOffset( name );
    BeginComplex_Internal( -id, type );
}

ItemStream::ItemItRW ItemStream::EndComplex()
{
    UTIL_ASSERT( m_NestedComplexLevel > 0 );
    m_NestedComplexLevel--;

    /*IMPORTANT: THIS CODE FAILED if there was a ReallocDataMB()
      between BeginComplex_Internal() and EndComplex(), as
      m_NestedComplexStack saves POINTERS, not OFFSETS
    m_NestedComplexStack[m_NestedComplexLevel]->m_Size =
        &m_DataMB[m_EndOffset] - (int8*)m_NestedComplexStack[m_NestedComplexLevel];
    // returning Iterator to Complex element just finished...
    uint32 item_begin_offset( (int8*)m_NestedComplexStack[m_NestedComplexLevel] - &m_DataMB[m_BeginOffset] );
    return ItemItRW( *this, item_begin_offset, item_begin_offset + m_NestedComplexStack[m_NestedComplexLevel]->m_Size );
    */

    // Get just-ended Complex header from its stored offset, set its total size and return its ItemItRW
    uint32 header_begin_offset( m_NestedComplexOffsetStack[m_NestedComplexLevel] );
    ItemHeader *p_header = (ItemHeader*)&m_DataMB[ header_begin_offset ];
    p_header->m_Size = m_EndOffset - header_begin_offset;

    return ItemItRW( *this, header_begin_offset, m_EndOffset );
}

void ItemStream::BeginComplex_Internal( int32 id, uint16 type )
{
    UTIL_ASSERT( type < cMax_ItemType );

    // we cannot check the whole complex size because it will be defined later
    size_type size = sizeof(ItemHeader);
    bool bReallocOk = ReallocDataMB( m_EndOffset + size );
    UTIL_ASSERT( bReallocOk );

    /*IMPORTANT: THIS CODE FAILED if there was a ReallocDataMB()
      between BeginComplex_Internal() and EndComplex(), as
      m_NestedComplexStack saves POINTERS, not OFFSETS
    ItemHeader *p_header = (ItemHeader*)&m_DataMB[m_EndOffset];
    *p_header = ItemHeader( id, type | cItemType_Complex, sizeof(ItemHeader) );
    m_EndOffset += sizeof(ItemHeader);
    m_NestedComplexStack[ m_NestedComplexLevel ] = p_header;
    */

    // Save just-begun Complex header in the stack to access it on EndComplex()
    ItemHeader *p_header = (ItemHeader*)&m_DataMB[m_EndOffset];
    *p_header = ItemHeader( id, type | cItemType_Complex, sizeof(ItemHeader) );
    m_NestedComplexOffsetStack[ m_NestedComplexLevel ] = m_EndOffset;
    m_EndOffset += sizeof(ItemHeader);

    m_NestedComplexLevel++;
}

/*! Returns string offset, allocating it if it's not present.

  \note String offsets start at 1 because they are stored as NEGATIVE
  Identifiers in in item headers so that they can be discriminated
  from actually numeric user-given IDs (which must be > 0). An offset
  of 0 cannot be made negative and this is discarded.
*/
uint32 ItemStream::GetStringOffset( const char *name )
{
    //find...
    uint32 offset = 1;
    while( offset < m_StringEndOffset )
    {
        if( 0 == strcmp( name, &m_StringMB[offset] ) )
            return offset;
        else
            offset += strlen(&m_StringMB[offset]) + 1;
    }

    //...or add
    uint32 size = strlen(name)+1;
    bool bReallocOK = ReallocStringMB( m_StringEndOffset + size );
    UTIL_ASSERT( bReallocOK );

    strcpy( &m_StringMB[m_StringEndOffset], name );
    uint32 last_offset = m_StringEndOffset;
    m_StringEndOffset += size;
    return last_offset;
}

uint32 ItemStream::QueryStringOffset( const char *name )
{
    //find...
    uint32 offset = 1;
    while( offset < m_StringEndOffset )
    {
        if( 0 == strcmp( name, &m_StringMB[offset] ) )
            return offset;
        else
            offset += strlen(&m_StringMB[offset]) + 1;
    }
    // return 0 if not found
    return 0;
}

//----------------------------------------------------------------
// ItemIt implementation
//----------------------------------------------------------------
ItemStream::ItemIt ItemStream::ItemIt::Find( int32 id ) const
{
    UTIL_ASSERT( 0 != id );
    // find...
    for( ItemStream::ItemIt it=*this; it.IsValid(); ++it )
        //\todo ALLOW SIGNED to find named items by id if( it.GetId() > 0 && id == it.GetId() )
        if( id == it.GetId() )
            return it;
    //...or return invalid it
    return ItemStream::ItemIt();
}

//! Find an Item by name (in the same level) after the current one
ItemStream::ItemIt ItemStream::ItemIt::Find( const char *name ) const
{
    // find...
    for( ItemStream::ItemIt it=*this; it.IsValid(); ++it )
        if( it.GetId() < 0 && 0 == strcmp(it.GetName(),name) )
            return it;
    //...or return invalid it
    return ItemStream::ItemIt();

    /*! \todo Alternative implementation that compares ids instead of strings...
    int32 id = QueryStringOffset(name);
    if( id != 0 ) return Find(-id);
    else return ItemIt::ItemIt();
    */
}

void ItemStream::ItemIt::ToBuffer( int8 *p_buffer ) const
{
    UTIL_ASSERT( !IsArray() );
    const int8 *p_data = GetDataPtr();
    for( unsigned i=0; i<GetItemHeader().m_Size-sizeof(ItemHeader); i++ ) p_buffer[i] = p_data[i];
}

//----------------------------------------------------------------
// ItemItRW implementation
//----------------------------------------------------------------
ItemStream::ItemItRW ItemStream::ItemItRW::Find( int32 id ) const
{
    UTIL_ASSERT( 0 != id );
    // find...
    for( ItemStream::ItemItRW it=*this; it.IsValid(); ++it )
        //\todo ALLOW SIGNED to find named items by id if( it.GetId() > 0 && id == it.GetId() )
        if( id == it.GetId() )
            return it;
    //...or return invalid it
    return ItemStream::ItemItRW();
}

//! Find an Item by name (in the same level) after the current one
ItemStream::ItemItRW ItemStream::ItemItRW::Find( const char *name ) const
{
    // find...
    for( ItemStream::ItemItRW it=*this; it.IsValid(); ++it )
        if( it.GetId() < 0 && 0 == strcmp(it.GetName(),name) )
            return it;
    //...or return invalid it
    return ItemStream::ItemItRW();

    /*! \todo Alternative implementation that compares ids instead of strings...
    int32 id = QueryStringOffset(name);
    if( id != 0 ) return Find(-id);
    else return ItemIt::ItemIt();
    */
}

void ItemStream::ItemItRW::ToBuffer( int8 *p_buffer ) const
{
    UTIL_ASSERT( !IsArray() );
    const int8 *p_data = GetDataPtr();
    for( unsigned i=0; i<GetItemHeader().m_Size-sizeof(ItemHeader); i++ ) p_buffer[i] = p_data[i];
}

void ItemStream::ItemItRW::TouchRecursive()
{
    Touch();
    if( IsComplex() )
        for( ItemStream::ItemItRW sub_it = GetSubItem(); sub_it.IsValid(); ++sub_it )
            sub_it.TouchRecursive();
}
void ItemStream::ItemItRW::UntouchRecursive()
{
    Untouch();
    if( IsComplex() )
        for( ItemStream::ItemItRW sub_it = GetSubItem(); sub_it.IsValid(); ++sub_it )
            sub_it.UntouchRecursive();
}

ItemStream::ItemItRW ItemStream::WriteItem_Internal( int32 id, const ItemIt &it )
{
    UTIL_ASSERT(it.IsValid());
    uint32 item_begin_offset = m_EndOffset;
    if( it.IsComplex() )
    {
        // Realloc ONCE for the whole Complex size, instead of per-sub-WriteItem_Internal
        //\todo THIS CODE IS UNDER TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        size_type size_in_bytes = it.GetItemHeader().m_Size;
        UTIL_ASSERT( 0 == size_in_bytes % 4 ); //4-aligned
        bool bReallocOk = ReallocDataMB( m_EndOffset + size_in_bytes );
        UTIL_ASSERT( bReallocOk );

        /*! Deep-copy requires recursive calls in order to add non-existent
          strings to the pool.
          \todo Of only ID were used (or a unique global string pool),
          the complex item memory sub-array could be copied at once,
          avoiding recursion.
          \todo This would ALSO solve the Touched flag copy, which
          must be explicitly set in the new item by now
        */
        BeginComplex_Internal( id, it.GetType() );
        {
            ItemIt sub_it = it.GetSubItem();
            while( sub_it.IsValid() )
            {
                if( sub_it.IsNamed() )
                    WriteItem_Internal( -GetStringOffset(sub_it.GetName()), sub_it );
                else
                    WriteItem_Internal( sub_it.GetId(), sub_it );
                ++sub_it;
            }
        }
        ItemStream::ItemItRW new_it = EndComplex();

        //TEMPORAL: Touched items must be copied as touched... it's required in S2 params stuff... ugly?
        if( it.IsTouched() )
        {
            //UTIL_LOG_WARNING( "Touching a copied item" );
            new_it.Touch();
        }
        return new_it;
    }
    else
    {
        size_type size_in_bytes = it.GetItemHeader().m_Size;
        UTIL_ASSERT( 0 == size_in_bytes % 4 ); //4-aligned

        bool bReallocOk = ReallocDataMB( m_EndOffset + size_in_bytes );
        UTIL_ASSERT( bReallocOk );

        // copy header and data in dwords
        const uint32 *p_input = reinterpret_cast<const uint32*>( &it.GetItemHeader() );
        uint32 *p_data = reinterpret_cast<uint32*>( &m_DataMB[m_EndOffset] );
        // copy in 4-byte packets
        for( unsigned int i=0; i<(size_in_bytes/4); i++ ) p_data[i] = p_input[i];

        // Set header Id to id
        ItemHeader *p_header = reinterpret_cast<ItemHeader*>( &m_DataMB[m_EndOffset] );
        p_header->m_Id = id;

        m_EndOffset += size_in_bytes;

        // Return ItemItRW
        return ItemItRW( *this, item_begin_offset, m_EndOffset );
    }
}

//---- Parse Item definitions from a string
const char *ItemStream::Parse( const char *str )
{
    /* Grammar:
       item_list -> item | item + ',' + item_list
       item -> name + '=' + value
       name -> Identifier
       value -> int32 | uint32 | float | double
                | nir_int32 | nir_uint32 | nir_float | nir_double
                | vec2f | vec3f | vec4f | quatf
                | string
                | array
                | complex
       int32 -> Int_Literal (+ 'i')
       float -> Float_Literal (+ 'f')
       nir_float -> Float_Literal + 'fir' + '(' + Float_Literal + ',' + Float_Literal + ',' + Float_Literal + ')'
       array -> '[' + (iufd) + value_list + ']'
       value_list -> value
                     | value + ',' + value
    */

    // Parse until EOF
    return ParseItemList( str, "\0" );
}

const char *ItemStream::ParseItemList( const char *str, const char *stop_chars )
{
    const char *p_begin = str;
    while( 0 != p_begin && !util::IsCharInSet(*p_begin,stop_chars) )
    {
        // Skip ALL comments
        do { p_begin = ParseComment(p_begin); } while ( 0 != p_begin && *p_begin == '/' );
        // Parse item, if any
        if( 0 != p_begin && !util::IsCharInSet(*p_begin,stop_chars) )
        {
            //std::cout << "After comment: " << *p_begin << std::endl;
            const char *p_end = ParseItem(p_begin);
            //UTIL_LOG( "Parsed %s, +remaining %s", p_begin, p_end );
            if( 0 == p_end )
            {
                //UTIL_LOG( "ParseItemList: finished after '%s'", p_begin );
                return 0;
            }
            p_begin = util::Skip(p_end,", \n"); //comma is NOT required, but if present, we just ignore it
            //\post At this point, p_begin may point to: next item | str_end | stop_char
        }
    }
    //if( p_begin ) UTIL_LOG( "ParseItemList: remaining %s", p_begin );
    return p_begin;
}

const char *ItemStream::ParseComment( const char *str )
{
    // Skip spaces, find '/', then skip whole line if found, and
    // return p_end AFTER the comment-line
    const char *p_begin = util::Skip(str,' ');
    if( *str == '/' )
    {
        const char *p_end = util::Find(p_begin,"\n");

        /* Just ignore comments...
        unsigned int length = p_end - p_begin;
        if( length < 128 )
        {
            char tmp_str[128];
            strncpy( tmp_str, p_begin, length );
            tmp_str[length] = '\0';
            std::cout << "[COMMENT] '" << tmp_str << "'" << std::endl;
        }
        */

        if( p_end ) p_end = util::Skip(p_end,"\n \t");
        return p_end;
    }
    else return p_begin;
}

const char *ItemStream::ParseItem( const char *str )
{
    // Extract param name
    const char *p_name_begin = util::Skip(str," ");
    const char *p_name_end = util::Find(p_name_begin,"\n =");
    if( 0 == p_name_end )
    {
        UTIL_LOG_WARNING( "ItemStream::ParseItem() Unnamed items not supported yet" );
        return 0; //\todo Allow nameless params by setting a unique numeric ID instead
    }

    // Check if its a valid token
    String32 param_name;
    unsigned int length = p_name_end - p_name_begin;
    if( IsLetter(*p_name_begin)
        || *p_name_begin=='#' ) //TEMP
    {
        if( length < param_name.GetMaxLength() )
        {
            strncpy( param_name.GetStr(), p_name_begin, length );
            param_name.GetStr()[length] = '\0';
            //std::cout << "Parsed param name '" << param_name << "'" << std::endl;
        }
        else
        {
            UTIL_LOG_ERROR( "ItemStream::ParseItem() Item name too long > %d", param_name.GetMaxLength() );
            return 0;
        }
    }
    else
    {
        UTIL_LOG_WARNING( "ItemStream::ParseItem() Item name should begin with a Letter, begins with %c", *p_name_begin );
        return 0;
    }

    // Extract value and store with proper name and type
    const char *p_value_begin = util::Skip(p_name_end,"\n =");
    const char *p_value_end = 0;
    switch( *p_value_begin )
    {
    case '_': p_value_end = ParseKeyword(param_name,p_value_begin); break;
    case '[': p_value_end = ParseArray(param_name,p_value_begin); break;
    case '^': p_value_end = ParseHexBuffer(param_name,p_value_begin); break;
    case '<': p_value_end = ParseNiR(param_name,p_value_begin); break;
    case '"': p_value_end = ParseString(param_name,p_value_begin); break;
    case '\'': p_value_end = ParseStringN(param_name,p_value_begin); break;
        //case '$': break; //p_value_end = ParseVariable(param_name,p_value_begin) break;
    case '{': p_value_end = ParseComplex(param_name,p_value_begin); break;
    default:
        if( IsDigit(*p_value_begin) || IsNumericSymbol(*p_value_begin) )
            p_value_end = ParseNumber(param_name,p_value_begin);
        break;
    }
    return p_value_end;
}

const char *ItemStream::ParseKeyword( const char *name, const char *str )
{
    const char *p_begin = util::Skip(str,'_');
    const char *p_end = util::Find(p_begin," =(\n");
    String64 keyword_str;
    unsigned int length = p_end - p_begin;
    UTIL_ASSERT( length < keyword_str.GetMaxLength() );
    strncpy( keyword_str, p_begin, length );
    keyword_str.GetStr()[length] = '\0';
    if( 0 == strcmp( keyword_str.GetStr(), "true" ) ) Write<bool>(name,true);
    else if( 0 == strcmp( keyword_str.GetStr(), "false" ) ) Write<bool>(name,false);
    else UTIL_LOG_ERROR( "Item %s uses unknown keyword _%s", name, keyword_str.GetStr() );
    return p_end;
}

const char *ItemStream::ParseNumber( const char *name, const char *str )
{
    // Extract number substring
    const char *p_begin = util::Skip(str," ");
    const char *p_end = util::SkipNotNL(p_begin,"0123456789+-.");
    String64 number_str;
    unsigned int length = p_end - p_begin;
    UTIL_ASSERT( length < number_str.GetMaxLength() );
    strncpy( number_str, p_begin, length );
    number_str.GetStr()[length] = '\0';
    // Convert to number and store it
    switch( (p_end) ? *p_end : 0 )
    {
        // Explicit types (skip suffix afterwards)
    case 'i': Write<int32>( name, (int32)atoi(number_str) );  p_end++; break;
    case 'u': Write<uint32>( name, (uint32)atoi(number_str) );  p_end++; break;
    case 'f': Write<float>( name, (float)atof(number_str) );  p_end++; break;
    case 'd': Write<double>( name, (double)atof(number_str) );  p_end++; break;
        // Implicit double (no suffix)
    default: Write<double>( name, (double)atof(number_str) );  break;
    }
    return p_end;
}

const char *ItemStream::ParseString( const char *name, const char *str )
{
    // Extract substring
    const char *p_begin = util::Skip(str,'"');
    if( 0 == p_begin ) return 0;
    const char *p_end = util::Find(p_begin,'"');
    if( 0 == p_end ) return 0;
    unsigned int length = p_end - p_begin;
    WriteString( name, p_begin, length );
    return p_end+1; //jump closing '"'
}
const char *ItemStream::ParseStringN( const char *name, const char *str )
{
    // Extract substring
    const char *p_begin = util::Skip(str,'\'');
    if( 0 == p_begin ) return 0;
    const char *p_end = util::Find(p_begin,'\'');
    if( 0 == p_end ) return 0;
    unsigned int length = p_end - p_begin;

    // Extract length substring
    const char *p_begin_number = util::Skip(p_end,"'");
    const char *p_end_number = util::SkipNotNL(p_begin_number,"0123456789");
    String64 number_str;
    unsigned int length_number = p_end_number - p_begin_number;
    UTIL_ASSERT( length_number < number_str.GetMaxLength() );
    strncpy( number_str, p_begin_number, length_number );
    number_str.GetStr()[length_number] = '\0';
    int size = atoi(number_str);
    // Write appropiate GString<N>
    if( size <= 16 ) Write( name, String16(p_begin,length) );
    else if( size <= 32 ) Write( name, String32(p_begin,length) );
    else if( size <= 64 ) Write( name, String64(p_begin,length) );
    else UTIL_LOG_ERROR( "ParseStringN: Unsupported String<%d> length requested. Use {16,32,64}", size );
    return p_end_number;
}

const char *ItemStream::ParseComplex( const char *name, const char *str )
{
    const char *p_begin = util::Find(str,'{');
    if( !p_begin )
    {
        UTIL_LOG_ERROR( "ParseComplex: Missing starting '{'" );
        return 0;
    }
    //\todo Generic type by now... Ideally, any specific type Complex should also accept an eType_Property_Group item
    BeginComplex( name, eType_Property_Group );
    const char *p_end = ParseItemList( p_begin+1, "}" );
    if( p_end && *p_end == '}' )
    {
        EndComplex();
        p_end++;
    }
    else
    {
        UTIL_LOG_ERROR( "ParseComplex: Unbalanced '{' with no ending '}' " );
        //UTIL_LOG_ERROR( "chars = %c, %c == %d, %d", *p_begin, p_end?*p_end:'Ñ', *p_begin, p_end?*p_end:'Ñ' );
        UTIL_LOG_ERROR( "chars = %c, %c == %d, %d", *p_begin, p_end?*p_end:'@', *p_begin, p_end?*p_end:'@' );
        return 0;
    }
    return p_end;
}

const char *ItemStream::ParseArray( const char *name, const char *str )
{
    const char *p_begin = str;
    UTIL_ASSERT( *p_begin == '[' );
    p_begin++;
    const char *p_end = util::Find( p_begin, ']' );
    if( !p_end ) return 0;
    // Array type (double by default)
    char array_type = 'd';
    if( util::IsCharInSet( *p_begin, "iufd" ) )
    {
        array_type = *p_begin;
        p_begin++;
    }
    // Array count
    const int cMaxArrayCount = 128;
    const int cMaxArrayItemSize = 8;
    int8 array_data[ cMaxArrayCount * cMaxArrayItemSize ];
    void *p_array_data( &array_data[0] );
    // Add values to array
    int array_count(0);
    while( p_begin && p_begin < p_end && array_count < cMaxArrayCount && *p_begin != ']' )
    {
        // Extract number substring
        p_begin = util::Skip(p_begin,", ");
        if( *p_begin != ']' )
        {
            const char *p_end_number = util::SkipNotNL(p_begin,"0123456789+-.");
            String64 number_str;
            unsigned int length = p_end_number - p_begin;
            UTIL_ASSERT( length < number_str.GetMaxLength() );
            strncpy( number_str, p_begin, length );
            number_str.GetStr()[length] = '\0';
            // Convert to number and store it
            switch( array_type )
            {
            case 'i': reinterpret_cast<int32*>(p_array_data)[array_count++] = (int32)atoi(number_str); break;
            case 'u': reinterpret_cast<uint32*>(p_array_data)[array_count++] = (uint32)atoi(number_str); break;
            case 'f': reinterpret_cast<float32*>(p_array_data)[array_count++] = (float32)atof(number_str); break;
            case 'd': reinterpret_cast<float64*>(p_array_data)[array_count++] = (float64)atof(number_str); break;
            default: break;
            }
            if( p_end_number && util::IsCharInSet( *p_end_number, "iufd" ) ) p_end_number++;
            p_begin = p_end_number;
        }
    }
    UTIL_ASSERT( array_count < cMaxArrayCount );
    switch( array_type )
    {
    case 'i': WriteArray( name, reinterpret_cast<int32*>(p_array_data), array_count ); break;
    case 'u': WriteArray( name, reinterpret_cast<uint32*>(p_array_data), array_count ); break;
    case 'f': WriteArray( name, reinterpret_cast<float32*>(p_array_data), array_count ); break;
    case 'd': WriteArray( name, reinterpret_cast<float64*>(p_array_data), array_count ); break;
    default: break;
    }
    return p_end+1; // Skip last ']'
}

static uint8 MapHexCharTo4bits( char c )
{
    if( c >= '0' && c <= '9' ) return uint8( c - '0' );
    else if( c >= 'A' && c <= 'F' ) return uint8( 10 + c - 'A' );
    else UTIL_LOG_ERROR( "MapHexCharTo4bits() with char %c out of hex range [0..F]", c );
    return 0;
}

const char *ItemStream::ParseHexBuffer( const char *name, const char *str )
{
    UTIL_ASSERT( *str == '^' );
    const char *p_begin = util::Skip(str,"^ ");
    if( 0 == p_begin ) return 0;
    const char *p_end = util::SkipNotNL(p_begin,"0123456789ABCDEF");
    unsigned int length( p_end - p_begin );
    UTIL_ASSERT( length % 2 == 0 );
    unsigned int buffer_size( length / 2 );
    uint8 *p_hex_buffer = AllocArray<uint8>( name, buffer_size );
    for( unsigned int i=0; i<buffer_size; i++, p_begin += 2 )
    {
        uint8 b47( MapHexCharTo4bits(*p_begin) );
        uint8 b03( MapHexCharTo4bits(*(p_begin+1)) );
        p_hex_buffer[i] = uint8( (b47 << 4) | b03 );
    }
    return p_end;
}

const char *ItemStream::ParseNiR( const char *name, const char *str )
{
    const char *p_begin = str;
    UTIL_ASSERT( *p_begin == '<' );
    p_begin++;
    const char *p_end = util::Find( p_begin, '>' );
    if( !p_end ) return 0;
    // Array type (double by default)
    char array_type = 'd';
    if( util::IsCharInSet( *p_begin, "iufd" ) )
    {
        array_type = *p_begin;
        p_begin++;
    }
    // Array count
    const int cMaxArrayCount = 3;
    const int cMaxArrayItemSize = 8;
    int8 array_data[ cMaxArrayCount * cMaxArrayItemSize ];
    void *p_array_data( &array_data[0] );
    // Add values to array
    int array_count(0);
    while( p_begin && p_begin < p_end && array_count < cMaxArrayCount && *p_begin != '>' )
    {
        // Extract number substring
        p_begin = util::Skip(p_begin,", ");
        if( *p_begin != '>' )
        {
            const char *p_end_number = util::SkipNotNL(p_begin,"0123456789+-.");
            String64 number_str;
            unsigned int length = p_end_number - p_begin;
            UTIL_ASSERT( length < number_str.GetMaxLength() );
            strncpy( number_str, p_begin, length );
            number_str.GetStr()[length] = '\0';
            // Convert to number and store it
            switch( array_type )
            {
            case 'i': reinterpret_cast<int32*>(p_array_data)[array_count++] = (int32)atoi(number_str); break;
            case 'u': reinterpret_cast<uint32*>(p_array_data)[array_count++] = (uint32)atoi(number_str); break;
            case 'f': reinterpret_cast<float32*>(p_array_data)[array_count++] = (float32)atof(number_str); break;
            case 'd': reinterpret_cast<float64*>(p_array_data)[array_count++] = (float64)atof(number_str); break;
            default: break;
            }
            if( p_end_number && util::IsCharInSet( *p_end_number, "iufd" ) ) p_end_number++;
            p_begin = p_end_number;
        }
    }
    UTIL_ASSERT( array_count <= cMaxArrayCount );
    switch( array_type )
    {
    case 'i': Write< Property_NIR_int32 >( name,
                                           Property_NIR_int32( reinterpret_cast<int32*>(p_array_data)[0],
                                                               reinterpret_cast<int32*>(p_array_data)[1],
                                                               reinterpret_cast<int32*>(p_array_data)[2] ) ); break;
    case 'u': Write< Property_NIR_uint32 >( name,
                                            Property_NIR_uint32( reinterpret_cast<uint32*>(p_array_data)[0],
                                                                 reinterpret_cast<uint32*>(p_array_data)[1],
                                                                 reinterpret_cast<uint32*>(p_array_data)[2] ) ); break;
    case 'f': Write< Property_NIR_float32 >( name,
                                             Property_NIR_float32( reinterpret_cast<float32*>(p_array_data)[0],
                                                                   reinterpret_cast<float32*>(p_array_data)[1],
                                                                   reinterpret_cast<float32*>(p_array_data)[2] ) ); break;
    case 'd': Write< Property_NIR_float64 >( name,
                                             Property_NIR_float64( reinterpret_cast<float64*>(p_array_data)[0],
                                                                   reinterpret_cast<float64*>(p_array_data)[1],
                                                                   reinterpret_cast<float64*>(p_array_data)[2] ) ); break;
    default: break;
    }
    return p_end+1; // Skip last '>'
}

//---- Quick-and-dirty Load/Save to file
bool ItemStream::LoadBin( const char *file_name )
{
    // fopen
    FILE *pFile = fopen( file_name, "rb" );
    if( !pFile ) return false;
    // Check magic number, which ALSO checks that file endianness == machine endianness
    uint32 tmp;
    fread( &tmp, 1, sizeof(uint32), pFile );
    UTIL_ASSERT( tmp == UTIL_ITEM_STREAM_MAGICNUMBER );
    /* TEMPORAL: Unnecessary, magicnumber already checks this!
    // Check that file and machine endianness are the same
    uint8 *p_file_endian( reinterpret_cast<uint8*>(&tmp) );
    uint32 machine_magicnumber( UTIL_ITEM_STREAM_MAGICNUMBER );
    uint8 *p_machine_endian( reinterpret_cast<uint8*>(&machine_magicnumber) );
    UTIL_ASSERT( p_file_endian[0] == p_machine_endian[0] &&
                 p_file_endian[1] == p_machine_endian[1] &&
                 p_file_endian[2] == p_machine_endian[2] &&
                 p_file_endian[3] == p_machine_endian[3] );
    */
    // Check version
    fread( &tmp, sizeof(uint32), 1, pFile );
    UTIL_ASSERT( tmp == UTIL_ITEM_STREAM_VERSION );
    // Load Header
    fread( &m_BeginOffset, sizeof(uint32), 1, pFile );
    fread( &m_EndOffset, sizeof(uint32), 1, pFile );
    fread( &m_StringEndOffset, sizeof(uint32), 1, pFile );
    // Load DataMB and DataMB sizes
    uint32 size_datamb;
    fread( &size_datamb, sizeof(uint32), 1, pFile );
    uint32 size_stringmb;
    fread( &size_stringmb, sizeof(uint32), 1, pFile );
    // Alloc and Load DataMB and StringMB data
    m_DataMB.Init(size_datamb);
    m_StringMB.Init(size_stringmb);
    if( size_datamb > 0 ) fread( m_DataMB.GetData(), 1, size_datamb, pFile );
    if( size_stringmb > 0 ) fread( m_StringMB.GetData(), 1, size_stringmb, pFile );
    m_NestedComplexLevel = 0;
    fclose( pFile );
    return true;
}

bool ItemStream::SaveBin( const char *file_name ) const
{
    UTIL_ASSERT( m_NestedComplexLevel == 0 ); //IS must be "closed" to save it...
    // fopen
    FILE *pFile = fopen( file_name, "wb" );
    if( !pFile ) return false;
    // Write magic number
    uint32 tmp;
    tmp = UTIL_ITEM_STREAM_MAGICNUMBER;
    fwrite( &tmp, sizeof(uint32), 1, pFile );
    // Write version
    tmp = UTIL_ITEM_STREAM_VERSION;
    fwrite( &tmp, sizeof(uint32), 1, pFile );
    // Write Header
    fwrite( &m_BeginOffset, sizeof(uint32), 1, pFile );
    fwrite( &m_EndOffset, sizeof(uint32), 1, pFile );
    fwrite( &m_StringEndOffset, sizeof(uint32), 1, pFile );
#define __USE_ITEMSTREAM_SAVE_EXACT_SIZE //\todo THIS should be default when sufficiently tested
#ifdef __USE_ITEMSTREAM_SAVE_EXACT_SIZE
    uint32 required_size_data = 4*((m_EndOffset+3)/4);
    uint32 required_size_string = 4*((m_StringEndOffset+3)/4);
    // Write DataMB and StringMB exact 4-aligned sizes
    fwrite( &required_size_data, sizeof(uint32), 1, pFile  );
    fwrite( &required_size_string, sizeof(uint32), 1, pFile );
    // Write DataMB and StringMB data
    if( required_size_data > 0 ) fwrite( m_DataMB.GetData(), 1, required_size_data, pFile );
    if( required_size_string > 0 ) fwrite( m_StringMB.GetData(), 1, required_size_string, pFile );
#else
    // Write DataMB and StringMB sizes
    tmp = m_DataMB.GetSize();
    fwrite( &tmp, sizeof(uint32), 1, pFile  );
    tmp = m_StringMB.GetSize();
    fwrite( &tmp, sizeof(uint32), 1, pFile );
    // Write DataMB and StringMB data
    if( m_DataMB.GetSize() > 0 ) fwrite( m_DataMB.GetData(), 1, m_DataMB.GetSize(), pFile );
    if( m_StringMB.GetSize() > 0 ) fwrite( m_StringMB.GetData(), 1, m_StringMB.GetSize(), pFile );
#endif
    fclose( pFile );
    return true;
}

bool ItemStream::LoadTxt( const char *file_name )
{
    FILE *pFile = fopen( file_name, "rt" );
    if( !pFile ) return false;
    bool bResult = util::ItemStreamFromFileTxt( *this, pFile );
    fclose( pFile );
    return bResult;
}

bool ItemStream::SaveTxt( const char *file_name ) const
{
    FILE *pFile = fopen( file_name, "wt" );
    if( !pFile ) return false;
    bool bResult = util::ItemStreamToFileTxt( *this, pFile );
    fclose( pFile );
    return bResult;
}

} //util
