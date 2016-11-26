#include "EmbeddingFactory.h"
#include <Geo/IObject.h>
// 2D
#include "Vertices2_In_Triangles2_Barycentric_Embedding.h"
#include "Vertices2_In_Triangles2_DeformationField_Embedding.h"
#include "Vertices2_In_Triangles2_MLS_Embedding.h"
#include "Vertices2_In_Triangles2_BEDF_Embedding.h"
// 3D
#include "Vertices3_In_Tetrahedrons3_Barycentric_Embedding.h"
#include "Vertices3_In_Tetrahedrons3_MLS_Embedding.h"

namespace geo {

/*TEMPORAL: maybe...
//template< DeformerT = MeshSolid2, EmbeddedT = Polygonal, EmbbeddingT = Barycentric >
class GIEmbedding
{

};
*/

////////////////////////////////////////////////////////////////
// EmbeddingFactory implementation
////////////////////////////////////////////////////////////////

// This could be double-dispatched, factorized and whatnot, but by now we only support a few combinations
IEmbedding* EmbeddingFactory_Create( const IObject* p_deformer, const IObject* p_embedded, EEmbeddingMethod em )
{
    EShapeType deformer_st = p_deformer->GetShapeInterface()->GetType();
    EShapeType embedded_st = p_embedded->GetShapeInterface()->GetType();
    IEmbedding* result(0);
    if( eShape_Polygonal2 == embedded_st
        && eShape_MeshSolid2 == deformer_st )
    {
        switch( em )
        {
        case eEM_Barycentric:
            {
                Editable_Vertices2_In_Triangles2_Barycentric_Embedding* pEmbedding = new Editable_Vertices2_In_Triangles2_Barycentric_Embedding();
                pEmbedding->Init( static_cast<const IObject2*>(p_deformer), static_cast<const IObject2*>(p_embedded) );
                result = pEmbedding;
            }
            break;
        case eEM_MLS:
            {
                Editable_Vertices2_In_Triangles2_MLS_Embedding* pEmbedding = new Editable_Vertices2_In_Triangles2_MLS_Embedding();
                pEmbedding->Init( static_cast<const IObject2*>(p_deformer), static_cast<const IObject2*>(p_embedded) );
                result = pEmbedding;
            }
            break;
        case eEM_BNDF:
            {
                Editable_Vertices2_In_Triangles2_DeformationField_Embedding* pEmbedding = new Editable_Vertices2_In_Triangles2_DeformationField_Embedding();
                pEmbedding->Init( static_cast<const IObject2*>(p_deformer), static_cast<const IObject2*>(p_embedded) );
                result = pEmbedding;
            }
            break;
            /*case eEM_BEDF:
              {
              Editable_Vertices2_In_Triangles2_BEDF_Embedding* pEmbedding = new Editable_Vertices2_In_Triangles2_BEDF_Embedding();
              pEmbedding->Init( static_cast<const IObject2*>(p_deformer), static_cast<const IObject2*>(p_embedded) );
              result = pEmbedding;
              }
              break;
            */
        default: break;
        }

        if( !result )
            GEO_LOG_ERROR("IObjectBase::Embed() Unsupported embedding type %d", (int)em );
        /*
        else
            GEO_LOG("IObjectBase::Embed() Successful embedding type %d of %d into %d", (int)em, (int)embedded_st, (int)deformer_st );
        */
        return result;
    }
    else if( eShape_TriSurface3 == embedded_st
             && eShape_TetSolid3 == deformer_st )
    {
        switch( em )
        {
        case eEM_Barycentric:
            {
                Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding* pEmbedding = new Editable_Vertices3_In_Tetrahedrons3_Barycentric_Embedding();
                pEmbedding->Init( static_cast<const IObject3*>(p_deformer), static_cast<const IObject3*>(p_embedded) );
                result = pEmbedding;
            }
            break;
        case eEM_MLS:
            {
                Editable_Vertices3_In_Tetrahedrons3_MLS_Embedding* pEmbedding = new Editable_Vertices3_In_Tetrahedrons3_MLS_Embedding();
                pEmbedding->Init( static_cast<const IObject3*>(p_deformer), static_cast<const IObject3*>(p_embedded) );
                result = pEmbedding;
            }
            break;
        default: break;
        }
        if( !result )
            GEO_LOG_ERROR("IObjectBase::Embed() Unsupported embedding type %d", (int)em );
        /*
        else
            GEO_LOG("IObjectBase::Embed() Successful embedding type %d of %d into %d", (int)em, (int)embedded_st, (int)deformer_st );
        */
        return result;
    }
    else
    {
        GEO_LOG_ERROR("IObjectBase::Embed() Unsupported embedding of %d into %d", (int)embedded_st, (int)deformer_st );
        return 0;
    }
}

} //namespace geo
