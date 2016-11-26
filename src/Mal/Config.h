#ifndef MAL_CONFIG_H
#define MAL_CONFIG_H

#if defined(USE_BOOST)
#  include <boost/static_assert.hpp>
#  define MAL_STATIC_ASSERT BOOST_STATIC_ASSERT
#else
#  define MAL_STATIC_ASSERT(x) //nothing
#endif

// Some useful macros
#ifdef PROFILE_FINAL
#  define MAL_ASSERT( x )
#  define MAL_LOG( x, ... )
#  define MAL_LOG_WARNING( x, ... )
#  define MAL_LOG_ERROR( x, ... )
#  define MAL_LOG_ERROR_IF( c, x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define MAL_ASSERT( x ) assert(x)
#  define MAL_LOG( x, ... ) { printf("<MAL LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define MAL_LOG_WARNING( x, ... ) { printf("<MAL WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define MAL_LOG_ERROR( x, ... ) { printf("<MAL ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define MAL_LOG_ERROR_IF( c, x, ... ) { if(c) { printf("<MAL ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); } }
#endif

#if defined(WIN32)
#  define finline _forceinline
#else
#  define finline inline
#endif

#endif //MAL_CONFIG_H
