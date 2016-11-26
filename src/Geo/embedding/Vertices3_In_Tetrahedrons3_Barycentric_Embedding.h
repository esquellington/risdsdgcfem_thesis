#ifndef GEO_EMBEDDING_V3INT3BE_H
#define GEO_EMBEDDING_V3INT3BE_H

#include <Geo/embedding/IEmbedding.h>
#include <Geo/IObject.h>

namespace geo {

/* A collection of 3D Vertices embedded in a collection of 3D Tetrahedrons
   - TriSurfaceShape3
   - Etc...
   Requisites:
   - Both deformer and embedded SDOF are vertex positions stored as Vec3 arrays
*/
class Vertices3_In_Tetrahedrons3_Barycentric_Embedding: public IEmbedding
{
public:
    struct BakedTetrahedronEmbeddedVertexSubset //\todo Change this name or pull my eyes out...
    {
        deformer_feature_index_type m_TetrahedronID; //\todo Unnecessary if m_vecTetrahedronVID is available, or if the crust elements are always 0..N-1 in the global array, with their index being implicit
        deformer_feature_index_type m_vecTetrahedronVID[4]; //\todo This is implicit in TetrahedronId, but requires accessing topology, it's cheaper storing ids here and processing only the vec_cage_dof/geometry
        Mat4x4 m_InvTetrahedronBCM; //\todo THIS IS IN GLOBAL COORDS!! \todo This is redundant and can could recomputed from VID and vec_default_cage_dof during embedding eval, if required, but if the ratio embedded_vertices/cage_eleements is large enough, it's irrelevant
        unsigned int m_FirstEmbeddedVertexIndex;
        unsigned int m_NumEmbeddedVertices;
    };

public:
    Vertices3_In_Tetrahedrons3_Barycentric_Embedding() : m_NumTEVS(0), m_NumEmbeddedVID(0), m_vecTEVS(0), m_vecEmbeddedVID(0) {}
    ~Vertices3_In_Tetrahedrons3_Barycentric_Embedding() { ClearBakedData(); }

    void SetBakedData( bool b_shared,
                       const Transform3& transform_e2d,
                       unsigned int num_tevs, unsigned num_embedded_vid,
                       const BakedTetrahedronEmbeddedVertexSubset* vec_tevs, const embedded_feature_index_type* vec_embedded_vid );
    void ClearBakedData();
    void Apply( const IObject* p_deformer, IObject* p_embedded ) const;

protected:
    Transform3 m_TransformE2D;
    unsigned int m_NumTEVS;
    unsigned int m_NumEmbeddedVID;
    const BakedTetrahedronEmbeddedVertexSubset* m_vecTEVS; //\todo Ideally, sorted by m_TetrahedronID
    const embedded_feature_index_type* m_vecEmbeddedVID;
};

class Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding: public Vertices3_In_Tetrahedrons3_Barycentric_Embedding
{
public:
    Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding() : m_allocTEVS(0), m_allocEmbeddedVID(0) {}
    ~Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding() { Clear(); }

    /* \todo See same comment in 2D version */
    void Init( const IObject3* p_deformer, const IObject3* p_embedded );
    void Clear();

private:
    void ClearEditData();

private:
    BakedTetrahedronEmbeddedVertexSubset* m_allocTEVS; //\todo Ideally, sorted by m_TetrahedronID
    embedded_feature_index_type* m_allocEmbeddedVID;
};

} //namespace geo

#endif //GEO_EMBEDDING_V3INT3BE_H
