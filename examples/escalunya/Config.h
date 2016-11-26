#ifndef ESCALUNYA_CONFIG_H
#define ESCALUNYA_CONFIG_H

// Some useful macros
#ifdef PROFILE_FINAL
#  define APP_ASSERT( x )
#  define APP_LOG( x, ... )
#  define APP_LOG_WARNING( x, ... )
#  define APP_LOG_ERROR( x, ... )
#  define APP_LOG_ASSERT( x, y, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define APP_ASSERT( x ) assert(x)
#  define APP_LOG( x, ... ) { printf("<APP LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_WARNING( x, ... ) { printf("<APP WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_ERROR( x, ... ) { printf("<APP ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_ASSERT( x, y, ... ) { if(!(x)) { printf("<APP ASSERT> "); printf( y, ##__VA_ARGS__ ); printf("\n"); assert(x); } }
#endif

#endif //ESCALUNYA_CONFIG_H
