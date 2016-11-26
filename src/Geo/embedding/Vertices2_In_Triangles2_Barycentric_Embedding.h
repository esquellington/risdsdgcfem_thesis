#ifndef GEO_EMBEDDING_V2INT2BE_H
#define GEO_EMBEDDING_V2INT2BE_H

#include <Geo/embedding/IEmbedding.h>
#include <Geo/IObject.h>

namespace geo {

/* A collection of 2D Vertices embedded in a collection of 2D Triangles.
   - Polygonal2_In_MeshSolid2
   - Etc...
   Requisites:
   - Both deformer and embedded SDOF are vertex positions stored as Vec2 arrays
*/
class Vertices2_In_Triangles2_Barycentric_Embedding: public IEmbedding
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
    Vertices2_In_Triangles2_Barycentric_Embedding() : m_NumTEVS(0), m_NumEmbeddedVID(0), m_vecTEVS(0), m_vecEmbeddedVID(0) {}
    ~Vertices2_In_Triangles2_Barycentric_Embedding() { ClearBakedData(); }

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

class Editable_Vertices2_In_Triangles2_Barycentric_Embedding: public Vertices2_In_Triangles2_Barycentric_Embedding
{
public:
    Editable_Vertices2_In_Triangles2_Barycentric_Embedding() : m_allocTEVS(0), m_allocEmbeddedVID(0) {}
    ~Editable_Vertices2_In_Triangles2_Barycentric_Embedding() { Clear(); }

    /* \todo Chose among several options for Init( IObject p_deformer, IObject p_embedded ):
       a) Embed <embedded_default_dof> into <deformer_default_dof>, assuming both transforms = Identity (thus ignoring them)
         => No extra memory required
       b) Embed <embedded_default_dof,embedded_transform> into <deformer_default_dof,deformer_transform>, using IObject transforms at embedding-time, but default_dof arrays
         => Needs to store relative transform from embedded to deformer captured at embedding-time
       c) Embed <embedded_default_dof,embedded_transform> into <deformer_dof,deformer_transform>, using IObject transforms at embedding-time and deformer dof, not default ones
         => Needs to store p_deformer DOF at embedding-time somehow, requires 1.5xDOF memory
       d) Embed <embedded_dof,embedded_transform> into <deformer_dof,deformer_transform>, using IObject transforms and vec_dof at embedding-time
         => Needs to store both p_deformer and p_embedded DOF at embedding-time as an embedding-pose, requires 2xDOF memory

       By now, we choose (b), as it allows independent transforms for
       deformer and embedded at embedding-time with little extra
       memory, but avoids the overhead of additional dof arrays as (c)
       and (d)
     */
    void Init( const IObject2* p_deformer, const IObject2* p_embedded );
    void Clear();

private:
    void ClearEditData();

private:
    BakedTriangleEmbeddedVertexSubset *m_allocTEVS; //\todo Ideally, sorted by m_TriangleID
    feature_index_type *m_allocEmbeddedVID;
};

} //namespace geo

#endif //GEO_EMBEDDING_V2INT2BE_H
