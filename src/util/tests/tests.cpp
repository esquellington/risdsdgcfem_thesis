//Test whole Util

#include "tests.h"
#include <util/Chronos.h>
#include <util/LogStream.h>
#include <util/ItemStream.h>
#include <util/ItemStreamSerialization.h>
#include <util/SimpleTracer.h>

#include <Mal/GSerialization.h>
#include <iostream>
#include <string.h>

namespace util
{

void WriteAtLevel( ItemStream &its, int level )
{
    for( int l=0; l<level; l++ )
        its.BeginComplex(l+1, 666);

    its.Write(555, 2.0);

    for( int l=0; l<level; l++ )
        its.EndComplex();
}

bool TestItemStream()
{
    ItemStream item_stream(1<<20, 1<<10);

    // Basic types
    item_stream.Write( 1, int32(-1) );
    item_stream.Write( 2, int16(-2) );
    item_stream.Write( 3, int8(-3) );
    item_stream.Write( 4, uint32(1) );
    item_stream.Write( 5, uint16(2) );
    item_stream.Write( 6, uint8(3) );
    item_stream.Write( 7, float(1.235f) );
    item_stream.Write( 8, double(1.235f) );

    // Basic type arrays
    int vec_int[] = { 1, 2, 3, 4, 5 };
    char vec_char[] = "estriiiiing";
    double vec_double[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7 };
    item_stream.WriteArray( 12, vec_int, 5 );
    item_stream.WriteArray( 13, vec_char, strlen(vec_char) );
    item_stream.WriteArray( 14, vec_double, 7 );

    // Complex
    item_stream.BeginComplex( 15, 666 );
      item_stream.Write( 2, -1.2345f );
      item_stream.WriteArray( 3, vec_char, strlen(vec_char) );
      item_stream.Write( 4, 123 );
      // Complex L1
      item_stream.BeginComplex( 123, 666 );
        item_stream.Write( 5, short(444) );
        item_stream.WriteArray( 3, vec_int, 5 );
        item_stream.Write( 4, 555 );
      item_stream.EndComplex();
    item_stream.EndComplex();

    // User-types
    item_stream.Write( 16, Vec2f(0.1,0.2) );
    item_stream.Write( 17, Vec3f(0.1,0.2,0.5) );
    item_stream.Write( 18, Quatf(0,0,0,1) );

    // Write at very deep level
    WriteAtLevel( item_stream, 5 );

    item_stream.WriteArray( "PRE-Last", "the end", 7 );
    item_stream.Write( "LAST", -1 );

    ItemStream::ItemIt it = item_stream.Begin();

    bool bResult = true;
    bResult = bResult && (-1 == it.Get<int32>());
    ++it;

    // Dump ItemStream
    std::cout << item_stream << std::endl;
    std::cout << "Count = " << item_stream.Count() << std::endl;

    // Search things
    std::cout << "Finding Item[1]..." << std::endl << item_stream.Find(1) << std::endl;
    std::cout << "Finding Item[LAST]..." << std::endl << item_stream.Find("LAST") << std::endl;

    // Clear and refill with named items
    item_stream.Clear();
    /*
    item_stream.Write( "int32", int(-1) );
    item_stream.Write( "uint32", uint32(1) );
    item_stream.Write( "float", float(1.235f) );
    item_stream.Write( "double", double(1.235f) );
    item_stream.WriteArray( "vec_int32", vec_int, 5 );
    item_stream.WriteArray( "vec_char", vec_char, strlen(vec_char) );
    item_stream.WriteArray( "vec_double", vec_double, 7 );
    */

    item_stream.Write( "shape_type", uint32(11) );
    item_stream.Write( "num_v", uint32(11) );
    item_stream.Write( "num_p", uint32(12) );
    item_stream.Write( "num_he", uint32(36) );
    item_stream.Write( "num_boundary_p", uint32(1) );
    item_stream.Write( "num_boundary_he", uint32(8) );

    // Save
    item_stream.SaveTxt( "TestItemStream_Save.txt" );
    item_stream.SaveBin( "TestItemStream_Save.bin" );

    item_stream.LoadBin( "TestItemStream_Load.bin" );
    item_stream.SaveTxt( "TestItemStream_Save.bin.txt" );

    // Clear and Dump again
    item_stream.Clear();
    std::cout << "Cleared..." << std::endl;
    std::cout << item_stream << std::endl;
    std::cout << "Count = " << item_stream.Count() << std::endl;

    // Create empty IS and add elements, forcing a Realloc()
    unsigned int cNumAddedElements( 10000 );
    std::cout << "Adding " << cNumAddedElements << " elements to an empty ItemStream...";
    Chronos chronos;
    chronos.ResetTime();
    ItemStream item_stream2(0,0);
    for( unsigned int i=0; i<cNumAddedElements; i++ )
    {
        // String Ids are VERY SLOW due to linear searches in ItemStream::GetStringOffset()
        char str[16];
        sprintf( str, "%d", i );
        item_stream2.Write( str, i );
        /* Non-string Ids are A LOT FASTER
           item_stream2.Write( 1, i );
        */
    }
    std::cout << "...took " << chronos.GetTime() << " seconds." << std::endl;

    return bResult;
}

} // namespace util

namespace  util
{

bool TestChronos()
{
    // Chronos (gettimeofday) has precicion = 1 microsecond
    Chronos chronos;
    chronos.ResetTime();
    float tmp = 0;
    for( int i=0; i<100000; i++ )
        tmp += float(i);
    double elapsed = chronos.GetTime();
    std::cout << "Chronos: Elapsed = " << elapsed << " Value = " << tmp << std::endl;

    // clock() has a LOT less precision
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    tmp = 0;
    for( int i=0; i<100000; i++ )
        tmp += float(i);
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "Chronos: Elapsed = " << elapsed << " Value = " << tmp << std::endl;

    return true;
}

void Scope1( float& value )
{
    UTIL_TRACE_BEGIN_SCOPE("Scope1");
    float tmp = value;
    UTIL_TRACE_LOCAL( "local_tmp", tmp );
    value *= value;
    UTIL_TRACE_GLOBAL( "global_value", value );
    UTIL_TRACE_END_SCOPE();
}

bool TestTracer()
{
    unsigned int cNumValues(5);
    float value(0);
    float vec_values[cNumValues];
    for( unsigned int i=0; i<cNumValues; i++ ) vec_values[i] = 0;
    UTIL_TRACE_INIT();
    for( unsigned int i=0; i<cNumValues; i++ )
    {
        UTIL_TRACE_BEGIN_SLICE( i );
        {
            UTIL_TRACE_BEGIN_SCOPE("Scope0");
            {
                UTIL_TRACE_LOCAL( "local_i", float(i) );
                value += float(i);
                UTIL_TRACE_GLOBAL( "global_value", value );
            }
            UTIL_TRACE_END_SCOPE();
            Scope1( value );
            vec_values[i] = value;
            UTIL_TRACE_GLOBAL_ARRAY( "global_vec_values", vec_values, cNumValues );
        }
        UTIL_TRACE_END_SLICE(true);
    }
    UTIL_TRACE_SHUTDOWN();
    return true;
}

} // namespace util
