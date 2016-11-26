#ifndef UTIL_GPOOL_DA_H
#define UTIL_GPOOL_DA_H

#include <util/Config.h>
#include <vector>

#define __ENABLE_ISVALID_TRACKING //!< Enable when debugging the pool

namespace util
{

/*! Iterable pool of T items with dynamic allocation in cBucketSize
  buckets and O(1) runtime operations.
  
  Pool functionality:
  - Creation is O(N)
  - New is O(1)
  - Delete is O(1)
  - Size() is O(1)

  Iterator functionality:
  - Begin() is O(1)
  - *Iterator is O(1)
  - ++Iterator is O(1)
*/
template < typename T, unsigned BucketSizePowerOfTwo = 5 >
class GPoolDA
{
public:
    class iterator;
    
    enum EConstants { cInvalidEID = uint16(0xFFFF),
                      cBucketSizePowerOfTwo = BucketSizePowerOfTwo,
                      cBucketSize = (1<<cBucketSizePowerOfTwo) };
    
private:
    struct EntryId;
    struct Entry;
    
public:
    GPoolDA( unsigned int hint_max_size )
    : m_HintMaxSize(uint16(hint_max_size))
    , m_Size(0)
    {
        UTIL_ASSERT( hint_max_size > 1 && hint_max_size <= (1<<16) );
        Clear();
    }
    ~GPoolDA() { DeallocBuckets(); }

    inline void Clear()
    {
        // Dealloc buckets
        DeallocBuckets();
        m_vecBuckets.clear();
        // Init first bucket
        uint16 bucket_id = NewBucket();
        // Init free and valid list
        m_FirstFreeId = EntryId(bucket_id,0);
        m_FirstValidId = EntryId(cInvalidEID);
        m_Size = 0;
    }
    
    //! \name Object creation/destruction
    //@{
    T *New()
    {
        /*
        LOG_WARNING_IF( m_Size > m_HintMaxSize,
                        "GPoolDA: hint_max_size = %d surpassed, a GPoolSA would have crashed here",
                        m_HintMaxSize );
        */
        UTIL_ASSERT(m_Size < (1<<16) );
        
        EntryId id = m_FirstFreeId;
        Entry &entry( GetEntry(id) );
        
#ifdef __ENABLE_ISVALID_TRACKING
        entry.m_IsValid = true;
#endif
        
        //---- Remove from free list
        if( entry.m_NextId == id )
        {
            // If no more free entries, alloc new bucket
            uint16 bucket_id = NewBucket();
            m_FirstFreeId = EntryId(bucket_id,0);
        }
        else
        {
            // Recompute first free and relink
            m_FirstFreeId = entry.m_NextId;
            GetEntry( entry.m_PrevId ).m_NextId = m_FirstFreeId;
            GetEntry( m_FirstFreeId ).m_PrevId = entry.m_PrevId;
        }

        //---- Add to valid list
        if( m_FirstValidId == EntryId(cInvalidEID) )
        {
            // Set as single valid
            m_FirstValidId = id;
            entry.m_PrevId = id;
            entry.m_NextId = id;
        }
        else
        {
            // Link to neighbours
            entry.m_PrevId = GetEntry( m_FirstValidId ).m_PrevId;
            entry.m_NextId = m_FirstValidId;
            // Link from neighbours
            GetEntry( entry.m_PrevId ).m_NextId = id;
            GetEntry( entry.m_NextId ).m_PrevId = id;
        }
        m_Size++;
        return &(entry.m_Item);
    }
    
    inline void Delete( T *p_item )
    {
        UTIL_ASSERT( IsValid(p_item) );
        Entry &entry( *(Entry*)p_item );        
        EntryId id( entry.m_Id );
        
#ifdef __ENABLE_ISVALID_TRACKING
        entry.m_IsValid = false;
#endif

        //---- Remove from valid list
        if( entry.m_NextId == id )
            m_FirstValidId = EntryId(cInvalidEID);
        else
        {
            // Recompute first valid and relink
            m_FirstValidId = entry.m_NextId;
            GetEntry( entry.m_PrevId ).m_NextId = m_FirstValidId;
            GetEntry( m_FirstValidId ).m_PrevId = entry.m_PrevId;
        }
        
        //---- Add to free list
        if( m_FirstFreeId == EntryId(cInvalidEID) )
        {
            // Set as single free
            m_FirstFreeId = id;
            entry.m_PrevId = id;
            entry.m_NextId = id;
        }
        else
        {
            // Link to neighbours
            entry.m_PrevId = GetEntry(m_FirstFreeId).m_PrevId;
            entry.m_NextId = m_FirstFreeId;
            // Link from neighbours
            GetEntry( entry.m_PrevId ).m_NextId = id;
            GetEntry( entry.m_NextId ).m_PrevId = id;
        }        
        m_Size--;
    }
    
    inline bool IsValid( const T* p_item ) const
    {        
        const Entry *entry = reinterpret_cast<const Entry*>(p_item);
#ifdef __ENABLE_ISVALID_TRACKING
        return ( entry != 0                       
                 && entry->m_Id.GetBucketId() < m_vecBuckets.size()
                 && entry->m_IsValid );
#else
        return ( entry != 0                       
                 && entry->m_Id.GetBucketId() < m_vecBuckets.size() );
#endif
    }
    //@}
    
    //! \name Object iteration (valid objects only)
    //@{
    inline bool IsEmpty() const { return 0 == m_Size; }
    inline unsigned int Size() const { return m_Size; }
    inline unsigned int Capacity() const { return (1<<16)-1; } //!< Limited by uint16 only
    inline iterator Begin() const { return iterator(*this); }
    inline iterator GetIterator( const T *p_item ) const { return iterator(*this,p_item); }
    //@}
    
    inline uint32 GetEntrySize() const { return sizeof(Entry); }
    
    //!\pre Requires T::operator<()
    void Sort()
        {
            if( m_Size < 2 ) return;
            // Sort whole list and reposition first_id (may have changed)
            m_FirstValidId = Sort_IPMS( m_FirstValidId, m_Size );
        }

    /*! In-place sublist sort.
      \return updated begin iterator.
      \note We use entry-id's directly, faster than "user" iterators.
    */
    EntryId Sort_IPMS( EntryId bid, unsigned int count )
        {
            // Split sublist in two halves
            unsigned int lc(count/2);
            unsigned int rc(count - lc);
            EntryId lid( bid );
            EntryId rid( bid );
            for( unsigned int i=0; i<lc; i++ ) rid = GetEntry(rid).m_NextId;
            // Sort independentrly
            if( lc > 1 ) lid = Sort_IPMS( lid, lc );
            if( rc > 1 ) rid = Sort_IPMS( rid, rc );
            // Merge inplace
            return Merge_IPMS( lid, lc, rid, rc );
        }

    /*! In-place sorted sublist merge.
      \pre  [left,right] are sorted
      \return updated begin iterator.

      \note We use entry-id's directly, faster than "user" iterators.      
        
      Idea: Iterate over left side, checking order. Whenever an
      element is bigger than first right element, copy the right one
      to the current left slot and insert this to the right side,
      reusing the newly freed "slot".
    */
    EntryId Merge_IPMS( EntryId lid, unsigned int lc, EntryId rid, unsigned int rc )
        {
            // Recompute begin_id
            EntryId begin_id = ( GetEntry(rid).m_Item < GetEntry(lid).m_Item ) ? rid : lid;
            // Merge
            EntryId left_id( lid );
            EntryId right_id( rid );
            unsigned int ri(0);
            unsigned int li(0);                        
            while( li < lc && ri < rc )
            {
                Entry &left_entry( GetEntry(left_id) );
                Entry &right_entry( GetEntry(right_id) );
                if( right_entry.m_Item < left_entry.m_Item )
                {
                    // Extract right entry from list
                    GetEntry(right_entry.m_PrevId).m_NextId = right_entry.m_NextId;
                    GetEntry(right_entry.m_NextId).m_PrevId = right_entry.m_PrevId;

                    // Reconnect right entry before current left entry
                    right_id = right_entry.m_NextId;
                    right_entry.m_PrevId = left_entry.m_PrevId;
                    right_entry.m_NextId = left_entry.m_Id;
                    
                    // Reconnect left and previous entries to newly inserted right entry
                    GetEntry(right_entry.m_PrevId).m_NextId = right_entry.m_Id;
                    left_entry.m_PrevId = right_entry.m_Id;
                    
                    // Advance
                    ri++;
                }
                else
                {
                    left_id = left_entry.m_NextId;
                    li++;
                }
            }
            return begin_id;
        }

    
public:
    //! Simple forward iterator
    class iterator
    {
    private:
        iterator( const GPoolDA &pool ) : m_pPool(&pool), m_Id(pool.m_FirstValidId) {} //!< Sets to first valid id (may be cInvalidEID)
        iterator( const GPoolDA &pool, const T* p_item ) //!< Create iterator pointing to a specific item
        : m_pPool(&pool)
        {
            UTIL_ASSERT( pool.IsValid(p_item) );
            m_Id = reinterpret_cast<const Entry*>(p_item)->m_Id;
        }
        
        friend class GPoolDA;
    public:
        iterator() : m_pPool(0), m_Id(cInvalidEID) {}
        iterator( const iterator &it ) : m_pPool(it.m_pPool), m_Id(it.m_Id) {}
        
        inline const T* operator->() const { UTIL_ASSERT(IsValid()); return &m_pPool->GetEntry(m_Id).m_Item; }
        inline const T& operator*() const { UTIL_ASSERT(IsValid()); return m_pPool->GetEntry(m_Id).m_Item; }
        
        inline T* operator->() { UTIL_ASSERT(IsValid()); return &m_pPool->GetEntry(m_Id).m_Item; }
        inline T& operator*() { UTIL_ASSERT(IsValid()); return m_pPool->GetEntry(m_Id).m_Item; }

        inline T* GetPtr() { UTIL_ASSERT(IsValid()); return &m_pPool->GetEntry(m_Id).m_Item; }
        inline const T* GetPtr() const { UTIL_ASSERT(IsValid()); return &m_pPool->GetEntry(m_Id).m_Item; }

        inline const GPoolDA *GetPool() const { return m_pPool; }
        
        inline void operator++()
        {
            UTIL_ASSERT(IsValid());
            
            // If its the last element the iterator becomes invalid, otherwise advance
            if( m_pPool->GetEntry(m_Id).m_NextId == m_pPool->m_FirstValidId )
                m_Id = EntryId(cInvalidEID);
            else
                m_Id = m_pPool->GetEntry(m_Id).m_NextId;
        }
        inline void operator+=( unsigned int inc )
        {
            UTIL_ASSERT(IsValid());
            unsigned int count(0);
            while( count++ < inc ) operator++();
        }
        inline bool IsValid() const { return m_Id != EntryId(cInvalidEID); }
        inline EntryId GetEID() const { return m_Id; }
        
    private:
        const GPoolDA *m_pPool;
        EntryId m_Id;
    };

private:
    finline Entry &GetEntry( EntryId eid ) const { return m_vecBuckets[eid.GetBucketId()]->GetEntry(eid.GetIndex()); }
    
    /*! Allocates a new bucket and initializes and links its entries
      as free. Returns new bucket Id.
    */
    inline uint16 NewBucket()
    {
        uint16 bucket_id( (uint16)m_vecBuckets.size() );
        
        // Reserve first bucket
        m_vecBuckets.push_back( new Bucket );
        
        // First
        Entry &lFirst( GetEntry( EntryId(bucket_id,0) ) );
#ifdef __ENABLE_ISVALID_TRACKING
        lFirst.m_IsValid = false;
#endif
        lFirst.m_Id = EntryId( bucket_id, 0 );
        lFirst.m_PrevId = EntryId( bucket_id, cBucketSize-1 );
        lFirst.m_NextId = EntryId( bucket_id, 1 );
        
        // Inner
        for( uint16 i=1; i<cBucketSize-1; i++ )
        {
            Entry &lInner( GetEntry( EntryId( bucket_id, i ) ) );
#ifdef __ENABLE_ISVALID_TRACKING                    
            lInner.m_IsValid = false;
#endif
            lInner.m_Id = EntryId( bucket_id, i );
            lInner.m_PrevId = EntryId( bucket_id, i-1 );
            lInner.m_NextId = EntryId( bucket_id, i+1 );
        }
        // Last
        Entry &lLast( GetEntry( EntryId( bucket_id, cBucketSize-1 ) ) );
#ifdef __ENABLE_ISVALID_TRACKING
        lLast.m_IsValid = false;
#endif
        lLast.m_Id = EntryId( bucket_id, cBucketSize-1 );
        lLast.m_PrevId = EntryId( bucket_id, cBucketSize-2 );
        lLast.m_NextId = EntryId( bucket_id, 0 );

        return bucket_id;
    }

    inline void DeallocBuckets()
    {
        for( uint16 it_bucket=0; it_bucket<m_vecBuckets.size(); it_bucket++ )
            delete m_vecBuckets[it_bucket];
    }

private:
    /*! Entry with pool annotations */
    struct Entry
    {
        T m_Item;
#ifdef __ENABLE_ISVALID_TRACKING
        bool m_IsValid; //!< \todo Actually unnecessary, as Validity is implicit in free/valid lists. ONLY left for debugging purposes
#endif
        EntryId m_Id; // required to delete it in O(1), otherwise we'd need to find bucket ptr in m_vecBuckets
        EntryId m_PrevId, m_NextId;
    };
    
    /*! Fixed-size bucket */
    struct Bucket
    {        
        Entry m_vecEntry[cBucketSize];
        finline Entry &GetEntry( uint16 aIndex ) { return m_vecEntry[aIndex]; }
        finline const Entry &GetEntry( uint16 aIndex ) const { return m_vecEntry[aIndex]; }
    };

private:
    //! 16bit Id [Bucket:Entry]
    struct EntryId
    {
        finline EntryId() : m_Bits( cInvalidEID ) {}
        finline EntryId( uint16 bits ) : m_Bits( bits ) {}
        finline EntryId( uint16 bucked_id, uint16 index ) : m_Bits( uint16( (bucked_id<<cBucketSizePowerOfTwo) | index) ) {}
        finline uint16 GetBucketId() const { return uint16(m_Bits >> cBucketSizePowerOfTwo); }
        finline uint16 GetIndex() const { return uint16( m_Bits & ((1<<cBucketSizePowerOfTwo) - 1) ); }
        finline bool operator==( const EntryId &entry_id ) const { return m_Bits == entry_id.m_Bits; }
        finline bool operator!=( const EntryId &entry_id ) const { return m_Bits != entry_id.m_Bits; }
        uint16 m_Bits;
    };

private:
    std::vector<Bucket*> m_vecBuckets;

    uint16 m_HintMaxSize;
    uint16 m_Size;
    EntryId m_FirstValidId;
    EntryId m_FirstFreeId;

    friend class iterator;
};

} // namespace util

#endif // UTIL_GPOOL_DA_H
