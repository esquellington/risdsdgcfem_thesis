#ifndef UTIL_CONFIG_H
#define UTIL_CONFIG_H

#include <pla_types.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define UTIL_ASSERT( x )
#  define UTIL_LOG( x, ... )
#  define UTIL_LOG_WARNING( x, ... )
#  define UTIL_LOG_ERROR( x, ... )
#  define UTIL_LOG_ERROR_IF( c, x, ... )
#  define UTIL_LOG_ASSERT( c, x, ... )
#else
#  include <assert.h>
#  define UTIL_ASSERT( x ) assert(x)
#  define UTIL_LOG( x, ... ) { printf("<UTIL LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define UTIL_LOG_WARNING( x, ... ) { printf("<UTIL WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define UTIL_LOG_ERROR( x, ... ) { printf("<UTIL ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define UTIL_LOG_ERROR_IF( c, x, ... ) { if(c) { printf("<UTIL ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); } }
#  define UTIL_LOG_ASSERT( c, x, ... )  { if(!(c)) { printf("<UTIL ASSERT FAILED> "); printf( x, ##__VA_ARGS__ ); printf("\n"); assert(c); } }
#endif

#endif //UTIL_CONFIG_H
