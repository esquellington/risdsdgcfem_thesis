#ifndef GEO_GEO_PARAMS_H
#define GEO_GEO_PARAMS_H

#include <util/ItemStream.h>

namespace geo {

util::ItemStream::ItemItRW QueryParams( util::ItemStream &params_is );
void SyncParams( util::ItemStream::ItemItRW params_it );

}

#endif //GEO_GEO_PARAMS_H
