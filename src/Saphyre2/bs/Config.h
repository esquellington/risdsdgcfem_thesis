#ifndef S2_BS_CONFIG_H
#define S2_BS_CONFIG_H

#include <Saphyre2/Config.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define BS_ASSERT( x )
#  define BS_LOG( x, ... )
#  define BS_LOG_WARNING( x, ... )
#  define BS_LOG_ERROR( x, ... )
#  define BS_LOG_ASSERT( c, x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define BS_ASSERT( x ) assert(x)
#  define BS_LOG( x, ... ) { printf("<BS LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define BS_LOG_WARNING( x, ... ) { printf("<BS WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define BS_LOG_ERROR( x, ... ) { printf("<BS ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define BS_LOG_ASSERT( c, x, ... )  { if(!(c)) { printf("<BS ASSERT FAILED> "); printf( x, ##__VA_ARGS__ ); printf("\n"); assert(c); } }
#endif

// Global Defines
static const unsigned cDefaultMaxObjInUniverse = 100;

// Enable debug stuff
#define __S2_BS_ENABLE_DS_STATS
#define __S2_BS_ENABLE_DS_PARAMS

#endif // S2_BS_CONFIG_H
