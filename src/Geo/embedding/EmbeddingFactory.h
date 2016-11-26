#ifndef GEO_EMBEDDING_EMBEDDING_FACTORY_H
#define GEO_EMBEDDING_EMBEDDING_FACTORY_H

#include <Geo/Config.h>
#include <Geo/embedding/IEmbedding.h>

namespace geo {

/* Default embedding creation
   \note Both objects are assumed to have up to date global Transform and local DOF
   \todo could accept an EEmbeddingType param (barycentric, harmonic, etc...)
   \todo Could be moved to an EmbeddingFactory class, even have an EmbeddingLibrary as Shapes... whatever
*/
IEmbedding* EmbeddingFactory_Create( const IObject* p_deformer, const IObject* p_embedded, EEmbeddingMethod em );

} //namespace geo

#endif //GEO_EMBEDDING_EMBEDDING_FACTORY_H
