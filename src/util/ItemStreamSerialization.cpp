#include <util/ItemStreamSerialization.h>
#include <stdio.h>

namespace util
{

static char g_Map_4bits_to_hex[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };

inline void TabulationToFileTxt( FILE *p_file, int num_spaces )
{
    for( int i=0; i<num_spaces; i++ ) fputs( " ", p_file );
}

static void ItemToFileTxt( ItemStream::ItemIt it, FILE *p_file, int level, int tabsize )
{
    if( it.IsSimple() )
    {
        TabulationToFileTxt( p_file, level*tabsize );
        if( it.IsNamed() ) fprintf( p_file, "%s ", it.GetName() );
        else fprintf( p_file, "#%d ", (int)it.GetId() );

        switch( it.GetType() )
        {
            /*
        case eType_Int8: fprintf( p_file, "%di\n", it.Get<int32>() ); break;
        case eType_UInt8: fprintf( p_file, "%du\n", it.Get<uint32>() ); break;
        case eType_Int16: fprintf( p_file, "%di\n", it.Get<int32>() ); break;
        case eType_UInt16: fprintf( p_file, "%du\n", it.Get<uint32>() ); break;
            */
        case eType_Int32: fprintf( p_file, "%di\n", it.Get<int32>() ); break;
        case eType_UInt32: fprintf( p_file, "%du\n", it.Get<uint32>() ); break;
        case eType_Float32: fprintf( p_file, "%ff\n", it.Get<float>() ); break;
        case eType_Float64: fprintf( p_file, "%lfd\n", it.Get<double>() ); break;
        case eType_Bool32: fprintf( p_file, "%s\n", it.Get<bool>() ? "_true" : "_false" ); break;
            // NiR
        case eType_Property_NIR_Int32:
            {
                const Property_NIR_int32 lNIR( it.Get<Property_NIR_int32>() );
                fprintf( p_file, "<i %d %d %d >\n", lNIR.m_Value, lNIR.m_Min, lNIR.m_Max );
            }
            break;
        case eType_Property_NIR_UInt32:
            {
                const Property_NIR_uint32 lNIR( it.Get<Property_NIR_uint32>() );
                fprintf( p_file, "<u %d %d %d >\n", lNIR.m_Value, lNIR.m_Min, lNIR.m_Max );
            }
            break;
        case eType_Property_NIR_Float32:
            {
                const Property_NIR_float32 lNIR( it.Get<Property_NIR_float32>() );
                fprintf( p_file, "<f %f %f %f >\n", lNIR.m_Value, lNIR.m_Min, lNIR.m_Max );
            }
            break;
        case eType_Property_NIR_Float64:
            {
                const Property_NIR_float64 lNIR( it.Get<Property_NIR_float64>() );
                fprintf( p_file, "<d %lf %lf %lf >\n", lNIR.m_Value, lNIR.m_Min, lNIR.m_Max );
            }
            break;
            // GString<N>
        case eType_String16: fprintf( p_file, "'%s'16 \n", it.Get<String16>().GetStr() ); break;
        case eType_String32: fprintf( p_file, "'%s'32 \n", it.Get<String32>().GetStr() ); break;
        case eType_String64: fprintf( p_file, "'%s'64 \n", it.Get<String64>().GetStr() ); break;
            //\todo GFlags<N>
            //case eType_Flags32: fprintf( p_file, "'%s'64 \n", it.Get<String64>().GetStr() ); break;
        default:
            UTIL_LOG_ERROR("ItemToFileTxt() Cannot serialize <Simple> item with pla_type_id = %d", it.GetType() );
            fprintf( p_file, "{}\n" );
            break;
        }
    }
    else if( it.IsArray() )
    {
        TabulationToFileTxt( p_file, level*tabsize );
        if( it.IsNamed() ) fprintf( p_file, "%s ", it.GetName() );
        else fprintf( p_file, "#%d ", (int)it.GetId() );

        switch( it.GetType() )
        {
        case eType_UInt8:
            {
                const uint8 *p_data = it.GetArrayPtr<uint8>();
                fputc( '^', p_file );
                for( unsigned int i=0; i < it.GetArrayCount(); i++ )
                {
                    fputc( g_Map_4bits_to_hex[ 0x0F & (p_data[i]>>4) ], p_file );
                    fputc( g_Map_4bits_to_hex[ 0x0F & p_data[i] ], p_file );
                }
                fputc( '\n', p_file );
            }
            break;
        case eType_Int8:
            fprintf( p_file, "\"%s\" \n", it.GetString() );
            break;
        case eType_Int32:
            {
                const int32 *p_data = it.GetArrayPtr<int32>();
                fprintf( p_file, "[i " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "%d ", p_data[i] );
                fputs( "]\n", p_file );
            }
            break;
        case eType_UInt32:
            {
                const uint32 *p_data = it.GetArrayPtr<uint32>();
                fprintf( p_file, "[u " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "%d ", p_data[i] );
                fputs( "]\n", p_file );
            }
            break;
        case eType_Float32:
            {
                const float32 *p_data = it.GetArrayPtr<float32>();
                fprintf( p_file, "[f " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "%f ", p_data[i] );
                fputs( "]\n", p_file );
            }
            break;
        case eType_Float64:
            {
                const float64 *p_data = it.GetArrayPtr<float64>();
                fprintf( p_file, "[d " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "%lf ", p_data[i] );
                fputs( "]\n", p_file );
            }
            break;
        case eType_String16:
            {
                const String16 *p_data = it.GetArrayPtr<String16>();
                fprintf( p_file, "[''16 " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "'%s' ", p_data[i].GetStr() );
                fputs( "]\n", p_file );
            }
            break;
        case eType_String32:
            {
                const String32 *p_data = it.GetArrayPtr<String32>();
                fprintf( p_file, "[''32 " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "'%s' ", p_data[i].GetStr() );
                fputs( "]\n", p_file );
            }
            break;
        case eType_String64:
            {
                const String64 *p_data = it.GetArrayPtr<String64>();
                fprintf( p_file, "[''64 " );
                for( unsigned int i=0; i< it.GetArrayCount(); i++ ) fprintf( p_file, "'%s' ", p_data[i].GetStr() );
                fputs( "]\n", p_file );
            }
            break;
        default:
            UTIL_LOG_ERROR("ItemToFileTxt() Cannot serialize <Array> item with pla_type_id = %d", it.GetType() );
            fprintf( p_file, "{}\n" );
            break;
        }
    }
    else if( it.IsComplex() )
    {
        //\todo check type == Object or Group or whatever??
        TabulationToFileTxt( p_file, level*tabsize );
        if( it.IsNamed() ) fprintf( p_file, "%s\n", it.GetName() );
        else fprintf( p_file, "#%d\n", (int)it.GetId() );

        // print {} with tabulation level
        TabulationToFileTxt( p_file, level*tabsize );
        fputs( "{\n", p_file );
        for( ItemStream::ItemIt sub_it=it.GetSubItem(); sub_it.IsValid(); ++sub_it )
            ItemToFileTxt( sub_it, p_file, level+1, tabsize );
        TabulationToFileTxt( p_file, level*tabsize );
        fputs( "}\n", p_file );
    }
    else
        UTIL_LOG_ERROR("Cannot serialize %s item that is !Simple && !Complex && !Array", it.GetName() );
}

bool ItemStreamToFileTxt( const ItemStream &is, FILE *p_file )
{
    UTIL_ASSERT( p_file );
    for( ItemStream::ItemIt it=is.Begin(); it.IsValid(); ++it )
        ItemToFileTxt( it, p_file, 0, 4 );
    return true;
}

bool ItemStreamFromFileTxt( ItemStream &is, FILE *p_file )
{
    UTIL_ASSERT( p_file );
    // Get filesize
    fseek( p_file , 0 , SEEK_END );
    uint32 file_size = ftell( p_file );
    rewind( p_file );
    // Read into memory buffer
    char *file_buffer = new char[file_size+1];
    if( 0 == file_buffer )
    {
        UTIL_LOG_ERROR( "ItemStreamFromFileTxt() OUT OF MEMORY allocating file buffer of size %d", file_size );
        fclose( p_file );
        return false;
    }
    fread( file_buffer, 1, file_size, p_file );
    file_buffer[file_size] = '\0';
    // Parse
#ifndef __ENABLE_ITEM_STREAM_REALLOC
    //Binary size not known, must grow dynamically
    UTIL_LOG_WARNING( "ItemStreamFromFileTxt() does NOT grow alloc memory, may crash. Use LoadBin() if possible" );
#endif
    is.Clear();
    is.Parse( file_buffer );
    // Cleanup
    delete [] file_buffer;
    return true;
}

} //namespace util
