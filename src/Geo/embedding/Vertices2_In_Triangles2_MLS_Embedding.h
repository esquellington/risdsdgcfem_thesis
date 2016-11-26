#ifndef GEO_EMBEDDING_V2INT2MLSE_H
#define GEO_EMBEDDING_V2INT2MLSE_H

#include <Geo/embedding/IEmbedding.h>
#include <Geo/IObject.h>

namespace geo {

/*
  \todo IObject::Embed( IObject, Method )
  - Methods = { Barycentric, DeformationField, MLS_{All,Cutoff,Topological}_{Gaussinan,Spline} }
  Configurable MLS:
  - Influences (options)
    - All influences per embedded vertex
    - Fixed K influences per embedded vertex (largest weight cutoff K)
    - Fixed K topological influences per embedded vertex
  - Weight function (options)
    - Truncated Gaussian
    - Spline
*/
class Vertices2_In_Triangles2_MLS_Embedding: public IEmbedding
{
public:
    enum EConstants { cMaxInfluencesPerVertex = 16 };

    /* Baked array of influences-per-vertex.
       - Weights could be quantized.
       - Global weight array could be used, with BEV storing first,count (specific sub-array), probably in 32b total (24:8)
    */
    struct BakedEmbeddedVertex
    {
        unsigned int m_NumInfluences;
        feature_index_type m_vecNID[cMaxInfluencesPerVertex];
        Real m_vecWeight[cMaxInfluencesPerVertex];
    };

public:
    Vertices2_In_Triangles2_MLS_Embedding() : m_NumEV(0), m_vecEV(0) {}
    ~Vertices2_In_Triangles2_MLS_Embedding() { ClearBakedData(); }

    void SetBakedData( bool b_shared,
                       const Transform2& transform_e2d,
                       unsigned int num_ev, const BakedEmbeddedVertex* vec_ev );

    void ClearBakedData();
    void Apply( const IObject* p_deformer, IObject* p_embedded ) const;

protected:
    Transform2 m_TransformE2D;

    unsigned int m_NumEV;
    const BakedEmbeddedVertex* m_vecEV;

    /*\todo Global influence array
    unsigned int m_NumEmbeddingInfluences;
    const Influence* m_vecEmbeddingInfluences;
    */
};

class Editable_Vertices2_In_Triangles2_MLS_Embedding: public Vertices2_In_Triangles2_MLS_Embedding
{
public:
    Editable_Vertices2_In_Triangles2_MLS_Embedding() : m_allocEV(0) {}
    ~Editable_Vertices2_In_Triangles2_MLS_Embedding() { Clear(); }
    void Init( const IObject2* p_deformer, const IObject2* p_embedded );
    void Clear();

private:
    void ClearEditData();

private:
    BakedEmbeddedVertex* m_allocEV;
};

} //namespace geo

#endif //GEO_EMBEDDING_V2INT2MLSE_H
