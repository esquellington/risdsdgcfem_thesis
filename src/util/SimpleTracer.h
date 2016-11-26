#ifndef UTIL_SIMPLE_TRACER_H
#define UTIL_SIMPLE_TRACER_H

#include <util/Config.h>
#include <util/ItemStream.h>

namespace util
{
class SimpleTracer
{
public:
    enum EConstants { cMaxNestedScope = 32 };
public:
    SimpleTracer() {}
    ~SimpleTracer() {}
    bool Init();
    bool ShutDown();
    bool BeginSlice( float slice_time );
    bool EndSlice( bool b_flush = false );
    bool BeginScope( const char* name );
    bool EndScope();
    bool TraceLocal( const char* name, const float& value );
    bool TraceLocal( const char* name, const float* vec_values, uint32 num_values );
    bool TraceGlobal( const char* name, const float& value );
    bool TraceGlobal( const char* name, const float* vec_values, uint32 num_values );

private:
    // Slice
    uint32 m_SliceId;
    float m_SliceTime;
    //ItemIt m_SliceItemIt; //Current Slice Complex
    // Scope
    uint32 m_NestedScopeLevel;
    String64 m_NestedScopeNameStack[cMaxNestedScope];
    // Slices itemstream
};
extern SimpleTracer* g_pTracer;
}

#define UTIL_TRACE_INIT() { if( 0 == util::g_pTracer ) { util::g_pTracer = new util::SimpleTracer; util::g_pTracer->Init(); } }
#define UTIL_TRACE_SHUTDOWN() { if( 0 != util::g_pTracer ) { util::g_pTracer->ShutDown(); delete util::g_pTracer; util::g_pTracer = 0; } }
#define UTIL_TRACE_BEGIN_SLICE(x) { if( 0 != util::g_pTracer ) { util::g_pTracer->BeginSlice( x ); } }
#define UTIL_TRACE_END_SLICE(b_flush) { if( 0 != util::g_pTracer ) { util::g_pTracer->EndSlice( b_flush ); } }
#define UTIL_TRACE_BEGIN_SCOPE(name) { if( 0 != util::g_pTracer ) { util::g_pTracer->BeginScope( name ); } }
#define UTIL_TRACE_END_SCOPE() { if( 0 != util::g_pTracer ) { util::g_pTracer->EndScope(); } }
#define UTIL_TRACE_LOCAL(name,value) { if( 0 != util::g_pTracer ) { util::g_pTracer->TraceLocal( name, value ); } }
#define UTIL_TRACE_LOCAL_ARRAY(name,vec_value,num_values) { if( 0 != util::g_pTracer ) { util::g_pTracer->TraceLocal( name, vec_value, num_values ); } }
#define UTIL_TRACE_GLOBAL(name,value) { if( 0 != util::g_pTracer ) { util::g_pTracer->TraceGlobal( name, value ); } }
#define UTIL_TRACE_GLOBAL_ARRAY(name,vec_value,num_values) { if( 0 != util::g_pTracer ) { util::g_pTracer->TraceGlobal( name, vec_value, num_values ); } }

#endif //UTIL_SIMPLE_TRACER_H
