#ifndef GEO_EMBEDDING_V2INT2DFE_H
#define GEO_EMBEDDING_V2INT2DFE_H

#include <Geo/embedding/IEmbedding.h>
#include <Geo/IObject.h>

namespace geo {

/*
  TEMP: Essentially the SAME AS Vertices2_In_Triangles2_Barycentric_Embedding, but changing Apply() to compute and use the deformation field in deformer nodes
*/
class Vertices2_In_Triangles2_DeformationField_Embedding: public IEmbedding
{
public:
    struct BakedTriangleEmbeddedVertexSubset //\todo Change this name or pull my eyes out...
    {
        feature_index_type m_TriangleID; //\todo Unnecessary if m_vecTriangleVID is available, or if the crust elements are always 0..N-1 in the global array, with their index being implicit
        feature_index_type m_vecTriangleVID[3]; //\todo This is implicit in TriangleId, but requires accessing topology, it's cheaper storing ids here and processing only the vec_cage_dof/geometry
        Mat3x3 m_InvTriangleBCM; //\todo THIS IS IN GLOBAL COORDS!! \todo This is redundant and can could recomputed from VID and vec_default_cage_dof during embedding eval, if required, but if the ratio embedded_vertices/cage_eleements is large enough, it's irrelevant
        unsigned int m_FirstEmbeddedVertexIndex;
        unsigned int m_NumEmbeddedVertices;
    };

public:
    Vertices2_In_Triangles2_DeformationField_Embedding() : m_NumTEVS(0), m_NumEmbeddedVID(0), m_vecTEVS(0), m_vecEmbeddedVID(0) {}
    ~Vertices2_In_Triangles2_DeformationField_Embedding() { ClearBakedData(); }

    void SetBakedData( bool b_shared,
                       const Transform2& transform_e2d,
                       unsigned int num_tevs, unsigned num_embedded_vid,
                       const BakedTriangleEmbeddedVertexSubset* vec_tevs, const feature_index_type* vec_embedded_vid );
    void ClearBakedData();
    void Apply( const IObject* p_deformer, IObject* p_embedded ) const;

protected:
    Transform2 m_TransformE2D;
    unsigned int m_NumTEVS;
    unsigned int m_NumEmbeddedVID;
    const BakedTriangleEmbeddedVertexSubset* m_vecTEVS; //\todo Ideally, sorted by m_TriangleID
    const feature_index_type* m_vecEmbeddedVID;
};

class Editable_Vertices2_In_Triangles2_DeformationField_Embedding: public Vertices2_In_Triangles2_DeformationField_Embedding
{
public:
    Editable_Vertices2_In_Triangles2_DeformationField_Embedding() : m_allocTEVS(0), m_allocEmbeddedVID(0) {}
    ~Editable_Vertices2_In_Triangles2_DeformationField_Embedding() { Clear(); }
    void Init( const IObject2* p_deformer, const IObject2* p_embedded );
    void Clear();

private:
    void ClearEditData();

private:
    BakedTriangleEmbeddedVertexSubset *m_allocTEVS; //\todo Ideally, sorted by m_TriangleID
    feature_index_type *m_allocEmbeddedVID;
};

} //namespace geo

#endif //GEO_EMBEDDING_V2INT2DFE_H
