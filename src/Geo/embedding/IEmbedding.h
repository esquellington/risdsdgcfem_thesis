#ifndef GEO_EMBEDDING_IEMBEDDING_H
#define GEO_EMBEDDING_IEMBEDDING_H

#include <Geo/Config.h>

namespace geo {

class IObject;

/*\name feature indes/count types
  - In principle, the embedding aspect of an IObject can have
    different feature counts than general shapes Polygonal, MeshSolid,
    TriSS, TetSS, etc...
  - In practice, if smaller than their specific feature_count_type,
    problems may arise... so keep them at uint32 unless specific
    optimizations are required.
*/
//@{
typedef uint32 deformer_feature_index_type;
typedef uint32 deformer_feature_count_type;
typedef uint32 embedded_feature_index_type;
typedef uint32 embedded_feature_count_type;
static const uint32 cDeformer_InvalidFeatureIndex = 0xFFFFFFFF;
static const uint32 cEmbedded_InvalidFeatureIndex = 0xFFFFFFFF;
//@}

class IEmbedding
{
public:
    IEmbedding() {}
    virtual ~IEmbedding() {}
    virtual void Apply( const IObject* p_deformer, IObject* p_embedded ) const = 0;
    //\todo Maybe virtual void Apply( const IObject* p_deformer, IObject* p_embedded, embedded dof range ) const = 0; //\todo Apply only to a subset of vertices
    //\todo Maybe virtual void Apply( const IObject* p_deformer, IObject* p_embedded, deformer element ) const = 0; //\todo Apply only to a triangle/tetrahedron
    //\todo Add Jacobian() related stuff, ideally localized to a specific embedded vertex => 3/4 deformer vertices
};

} //namespace geo

#endif //GEO_EMBEDDING_IEMBEDDING_H
