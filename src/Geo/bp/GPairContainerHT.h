#ifndef GEO_BP_G_PAIR_CONTAINER_HT_H
#define GEO_BP_G_PAIR_CONTAINER_HT_H

#include <Geo/Config.h>
#include <Geo/bp/Pair.h>

//#include <util/GPoolSA.h>
#include <util/GPoolDA.h>
#define GPool GPoolDA

namespace geo { namespace bp
{

/*! Persistent Pair allocator and container with specialized ops:
  - New is Expected O(1), Worst O(2^HashBits)
  - Remove is O(1) due to HashTable
  - Find is Expected O(1), Worst O(NumPairs)
  - Persist is O(1)
  - Clear is O(aMaxPairs + 2^HashBits)
  - Update is O(NumPairs)

  Fast Add/Find/Remove is implemented using a very simple Hash Table
  (HT):  
  - Fixed 2^N dimension that allows fast "X mod 2^N"  
  - The hash function is designed to be fast and symmetric on
    (p_proxy1,p_proxy2) and to remove lower bits that are constant due
    to pointer (p_proxy1,p_proxy2) alignment.    
  - The HT is only rebuild when a new pair cannot be Added.
  - Rebuilding takes O(2^N) time.
    - HashCells have some additional slots (cNumSafetyEIC) reserved to
      minimize full rebuilds.      
  - Param HashBits is a tradeoff between Space and Hash Collision
    probability. As Rebuild takes O(2^HashBits), small values [8..10]
    should be used in general.
*/
template <unsigned HashBits, unsigned AlignmentBits, unsigned NumSafetySlots, unsigned MaxHashCollisions>
class GPairContainerHT
{
public:
    typedef util::GPool<Pair> CPoolPair;
    typedef CPoolPair::iterator iterator;
    
public:
    GPairContainerHT( unsigned int aMaxPairs )
    : m_TimeStamp(0), m_PoolPairs(aMaxPairs)
    {
        GEO_ASSERT( (aMaxPairs) < (1<<16) );
        m_BufferEntries = new Pair*[aMaxPairs + (cNumHashCells*cNumSafetyEIC)];
        HT_Rebuild();
    }
    ~GPairContainerHT()
    {
        delete m_BufferEntries;
    }
    
    //\name Non-persistent Interface
    //@{
    /*! Allocate a New non-persistent pair efficiently (no map)
      Non-Persistent pairs are ONLY removed with Clear()
    */
    inline Pair *NewNP( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        Pair *pPair = m_PoolPairs.New();
        pPair->m_pProxy1 = p_proxy1;
        pPair->m_pProxy2 = p_proxy2;        
        pPair->m_TimeStamp = m_TimeStamp;
        pPair->m_State = Pair::eNew;
        return pPair;
    }
    //@}
    
    //\name Persistent Interface
    //@{
    inline Pair *New( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        Pair *pPair = m_PoolPairs.New();
        pPair->m_pProxy1 = p_proxy1;
        pPair->m_pProxy2 = p_proxy2;
        pPair->m_TimeStamp = m_TimeStamp;
        pPair->m_State = Pair::eNew;
        if( !HT_Add(pPair,p_proxy1,p_proxy2) ) HT_Rebuild(); //< If Add fails, HT_Rebuild will add ALL pairs in the pool
        return pPair;
    }
    inline Pair *Find( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        return HT_Find(p_proxy1,p_proxy2);
    }
    inline Pair *FindOrNew( Proxy *p_proxy1, Proxy *p_proxy2, bool &aFound )
    {
        Pair *pPair = HT_Find(p_proxy1,p_proxy2);
        if( 0 != pPair ) { aFound = true; return pPair; }
        else { aFound = false; return New(p_proxy1,p_proxy2); }
    }    
    inline void Remove( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        Pair *pPair = HT_Remove(p_proxy1,p_proxy2);
        m_PoolPairs.Delete(pPair);
    }
    //! Called after an Update(), flags the pair as persistent (so that next update may vanish it, but NOT remove it
    inline void Persist( Pair *p_pair ) const { p_pair->m_TimeStamp = m_TimeStamp; p_pair->m_State = Pair::ePersistent; }
    //! Clear all pairs. It's expensive, so early-out if there are no pairs.
    inline void Clear() { m_TimeStamp = 0; if( m_PoolPairs.Size() > 0 ) { m_PoolPairs.Clear(); HT_Clear(); } }

    inline uint32 GetTimeStamp() const { return m_TimeStamp; }
    
    //! Purges old pairs and flags the rest as Vanished
    void Update()
    {        
        iterator it_pair=m_PoolPairs.Begin();
        while( it_pair.IsValid() )
        {
            if( it_pair->m_TimeStamp != m_TimeStamp )
            {
                iterator it_delete(it_pair);
                ++it_pair;
                HT_Remove(it_delete->m_pProxy1,it_delete->m_pProxy2);
                m_PoolPairs.Delete( it_delete.GetPtr() );
            }
            else
            {
                it_pair->m_State = Pair::eVanished;
                ++it_pair;
            }
        }
        m_TimeStamp++;
    }
    //@}
    
    inline unsigned int Size() const { return m_PoolPairs.Size(); }
    inline iterator Begin() const { return m_PoolPairs.Begin(); }
    
private:
    inline uint32 HT_Hash( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        //Very simple "X mod 2^N" hash of lower N=cHashBits from the subset of significative bits [32..cAlignmentBits]
#ifdef __LP64__
        return ( (reinterpret_cast<uint64>(p_proxy1) ^ reinterpret_cast<uint64>(p_proxy2)) >> cAlignmentBits ) & cHashMask;
#else
        return ( (reinterpret_cast<uint32>(p_proxy1) ^ reinterpret_cast<uint32>(p_proxy2)) >> cAlignmentBits ) & cHashMask;
#endif
        
    }

    //! Add if possible, return false if not (which will spawn a RebuildHT)
    bool HT_Add( Pair *p_pair, Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        GEO_ASSERT( 0 == HT_Find(p_proxy1,p_proxy2) );
        uint32 hash( HT_Hash(p_proxy1,p_proxy2) );
        HashCell &hc( m_HashCells[hash] );
        if( hc.m_NumEntries < hc.m_MaxEntries ) { m_BufferEntries[ hc.m_FirstEntryIdx + hc.m_NumEntries++ ] = p_pair; return true; }
        else return false;
    }
    Pair *HT_Find( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        uint32 hash( HT_Hash(p_proxy1,p_proxy2) );
        HashCell &hc( m_HashCells[hash] );
        for( int it_eic=0; it_eic<hc.m_NumEntries; it_eic++ )
        {
            Pair *pPair( m_BufferEntries[hc.m_FirstEntryIdx+it_eic] );            
            if( (pPair->m_pProxy1 == p_proxy1 && pPair->m_pProxy2 == p_proxy2)
                || (pPair->m_pProxy1 == p_proxy2 && pPair->m_pProxy2 == p_proxy1) ) //\todo USE 0 == (a1-b1) | (a2-b2) trick?
                return pPair;
            /* \todo USE 0 == (a1-b1) | (a2-b2) trick?
            if( (0 == (int32(pPair->m_pProxy1)-int32(p_proxy1)) | (int32(pPair->m_pProxy2)-int32(p_proxy2)))
                ||
                (0 == (int32(pPair->m_pProxy1)-int32(p_proxy2)) | (int32(pPair->m_pProxy2)-int32(p_proxy1))) )
                return pPair;                
            */
        }
        return 0;
    }
    Pair *HT_Remove( Proxy *p_proxy1, Proxy *p_proxy2 )
    {
        uint32 hash( HT_Hash(p_proxy1,p_proxy2) );
        HashCell &hc( m_HashCells[hash] );
        for( int it_eic=0; it_eic<hc.m_NumEntries; it_eic++ )
        {
            Pair *pPair( m_BufferEntries[hc.m_FirstEntryIdx+it_eic] );
            if( (pPair->m_pProxy1 == p_proxy1 && pPair->m_pProxy2 == p_proxy2)
                || (pPair->m_pProxy1 == p_proxy2 && pPair->m_pProxy2 == p_proxy1) )
            {
                m_BufferEntries[hc.m_FirstEntryIdx + it_eic] = m_BufferEntries[hc.m_FirstEntryIdx + --hc.m_NumEntries];
                return pPair;
            }
            /* \todo USE 0 == (a1-b1) | (a2-b2) trick?
            if( (0 == (int32(pPair->m_pProxy1)-int32(p_proxy1)) | (int32(pPair->m_pProxy2)-int32(p_proxy2)))
                ||
                (0 == (int32(pPair->m_pProxy1)-int32(p_proxy2)) | (int32(pPair->m_pProxy2)-int32(p_proxy1))) )
            {
                m_BufferEntries[hc.m_FirstEntryIdx + it_eic] = m_BufferEntries[hc.m_FirstEntryIdx + --hc.m_NumEntries];
                return pPair;
            }
            */            
        }
        // Error if not found
        GEO_ASSERT( false );
        return 0;
    }
    void HT_Clear()
    {
        // Clear cells        
        for( int it_hc=0; it_hc<cNumHashCells; it_hc++ ) m_HashCells[it_hc].Clear(); //\todo Use MemZero()??
        //sal::MemorySet( &m_HashCells[0], 0, cNumHashCells*sizeof(HashCell) );
    }
    
    /*! Rebuilds the whole HT. \todo Split into a temptative rebuild
        method and an external control loop that retries if the
        distribution is bad (ex: Ensure that size(HC_i) < cMaxEIC << 255)
    */
    void HT_Rebuild()
    {        
        HT_Clear();
        
        // Count entries per hash cell
        for( iterator it_pair=m_PoolPairs.Begin(); it_pair.IsValid(); ++it_pair )
        {
            uint32 hash( HT_Hash(it_pair->m_pProxy1,it_pair->m_pProxy2) );
            GEO_ASSERT( m_HashCells[hash].m_NumEntries < (255-cNumSafetyEIC) ); //\todo uint8 entries per cell!!
            m_HashCells[hash].m_NumEntries++;
        }
        // Allocate per HC subarray
        int last_eic(0);
        for( int it_hc=0; it_hc < cNumHashCells; it_hc++ )
        {
            HashCell &hc( m_HashCells[it_hc] );
            hc.m_FirstEntryIdx = static_cast<uint16>(last_eic);
            hc.m_MaxEntries = static_cast<uint8>(hc.m_NumEntries + cNumSafetyEIC);
            hc.m_NumEntries = 0;
            last_eic += hc.m_MaxEntries;
        }
        // Add entries to HC
        for( iterator it_pair=m_PoolPairs.Begin(); it_pair.IsValid(); ++it_pair )
        {
            uint32 hash( HT_Hash(it_pair->m_pProxy1,it_pair->m_pProxy2) );
            HashCell &hc( m_HashCells[hash] );
            m_BufferEntries[ hc.m_FirstEntryIdx + hc.m_NumEntries++ ] = it_pair.GetPtr();
        }
    }
    
private:
    uint32 m_TimeStamp;
    CPoolPair m_PoolPairs;

    //! \name Hash-Table stuff
    //@{
    enum EConstants {
        cHashBits         = HashBits,
        cNumHashCells     = (1<<HashBits),
        cHashMask         = (1<<HashBits)-1,
        cAlignmentBits    = AlignmentBits,
        cNumSafetyEIC     = NumSafetySlots,
        cMaxEIC           = MaxHashCollisions
    };
    Pair **m_BufferEntries;
    struct HashCell
    {
        uint16 m_FirstEntryIdx;
        uint8 m_NumEntries;
        uint8 m_MaxEntries;
        inline void Clear() { *reinterpret_cast<uint32*>(this) = 0; }
    };
    HashCell m_HashCells[cNumHashCells];
    //@}
};

/*! Default GPairContainerHT<> specialization:
  - 2^8 hash cells
  - 2-bit alignment (comes from dword-aligned struct GPool::Entry)
  - 2 safety slots per cell
  - 10 hash collisions allowed before rehashing
*/
typedef GPairContainerHT<8,2,2,10> PairContainer;

}} //namespace geo::bp

#endif // GEO_BP_G_PAIR_CONTAINER_HT_H
