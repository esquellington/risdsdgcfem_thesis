#ifndef GEO_GEO_H
#define GEO_GEO_H

/*TEMP: Unnecessary.... but maybe geo should include "all root functionality" to offer single header include to users
#include <Geo/Config.h>
#include <Geo/shape/IShape.h>
#include <Geo/IObject.h>
*/

namespace geo {

bool Init();
bool ShutDown();
bool IsInialized();

}

#endif
