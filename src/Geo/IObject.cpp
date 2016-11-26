#include "IObject.h"
#include <Geo/embedding/IEmbedding.h>
#include <Geo/embedding/EmbeddingFactory.h>

namespace geo {

/*! Lightweight Embedding annotation class
  \todo Consider alloc from a Pool if it's uniform
  \todo Consider Shared and Exclusive embedding
  \note m_pEmbedded NEEDS to be saved, and NEEDS to be non-const to
        use it in SyncEmbedding() const in lazy evaluation methods that ARE
        const (GetTransform(), GetVecDOF)
*/
struct EmbeddingData
{
    const IObjectBase* m_pDeformer;
    IObjectBase* m_pEmbedded;
    IEmbedding* m_pEmbedding;
    uint32 m_KnownDeformerTimeStamp;
    finline EmbeddingData() : m_pDeformer(0), m_pEmbedding(0), m_KnownDeformerTimeStamp(0) {}
    ~EmbeddingData() { if( 0 != m_pEmbedding ) delete m_pEmbedding; } //\todo Exclusive IEmbedding ownership by now
};

IObjectBase::~IObjectBase()
{
    if( 0 != m_pEmbeddingData ) delete m_pEmbeddingData;
}

// \todo Explicit default embedding, create IEmbedding from current config
void IObjectBase::Embed( const IObject* p_deformer, EEmbeddingMethod em )
{
    // Remove current embedding
    if( 0 != m_pEmbeddingData )
    {
        SyncEmbedding(); //We sync before unembedding to preserve last up-to-date embedded configuration
        delete m_pEmbeddingData;
        m_pEmbeddingData = 0;
        //\note The embedded object transform is in absolute coords and therefore there is no need to change it here despite having been updated by the deformer so far.
    }
    if( 0 == p_deformer ) return;
    // Create new embedding
    IEmbedding *pEmbedding = EmbeddingFactory_Create( p_deformer, this, em );
    if( 0 != pEmbedding  )
    {
        m_pEmbeddingData = new EmbeddingData();
        m_pEmbeddingData->m_pDeformer = static_cast<const IObjectBase*>(p_deformer);
        m_pEmbeddingData->m_pEmbedded = static_cast<IObjectBase*>(this);
        m_pEmbeddingData->m_pEmbedding = pEmbedding;
        m_pEmbeddingData->m_KnownDeformerTimeStamp = m_pEmbeddingData->m_pDeformer->GetTimeStamp()-1; //\todo DIFFERENT timestamp at init to force sync asap
    }
}

bool IObjectBase::SyncEmbedding() const
{
    if( 0 != m_pEmbeddingData
        && m_pEmbeddingData->m_KnownDeformerTimeStamp != m_pEmbeddingData->m_pDeformer->GetTimeStamp() )
    {
        //GEO_LOG( "IObjectBase::SyncEmbedding()" );
        // Apply embedding to Transform and DOF
        m_pEmbeddingData->m_pEmbedding->Apply( m_pEmbeddingData->m_pDeformer, m_pEmbeddingData->m_pEmbedded );
        m_pEmbeddingData->m_pEmbedded->Touch();
        // Set known deformer timestamp as up to date
        m_pEmbeddingData->m_KnownDeformerTimeStamp = m_pEmbeddingData->m_pDeformer->GetTimeStamp();
        return true;
    }
    else
        return false;
}

} //namespace geo
