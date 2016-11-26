#ifndef S2_NS_CONFIG_H
#define S2_NS_CONFIG_H

#include <pla_types.h>
#include <Saphyre2/Config.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define NS_ASSERT( x )
#  define NS_LOG( x, ... )
#  define NS_LOG_WARNING( x, ... )
#  define NS_LOG_ERROR( x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define NS_ASSERT( x ) assert(x)
#  define NS_LOG( x, ... ) { printf("<NS LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define NS_LOG_WARNING( x, ... ) { printf("<NS WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define NS_LOG_ERROR( x, ... ) { printf("<NS ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#endif

// Enable debug stuff
#define __S2_NS_ENABLE_STATS
#define __S2_NS_ENABLE_PARAMS

// OPTIONAL SSE
#if __SSE4_1__
#  define __S2_NS_ENABLE_SIMD
#endif

#endif // S2_NS_CONFIG_H
