#ifndef GEO_UTIL_GSIMPLESPATIALHASH_H
#define GEO_UTIL_GSIMPLESPATIALHASH_H

#include "../Config.h"
#include <vector>
#include <memory.h> //TEMPORAL req by memset()

namespace geo
{

template <typename EntryIndexT, typename EntryT>
class GSimpleSpatialHash
{
public:
    typedef EntryIndexT entry_index_type;
    typedef EntryT entry_type;

private:
    struct CellCoords
    {
        uint8 i,j,k;
        inline CellCoords() {}
        inline CellCoords( uint8 a, uint8 b, uint8 c ) : i(a), j(b), k(c) {}
        inline CellCoords( uint32 a, uint32 b, uint32 c ) : i(uint8(a)), j(uint8(b)), k(uint8(c)) {}
        inline uint8 &operator[](int index) { return reinterpret_cast<uint8*>(this)[index]; }
        inline const uint8 &operator[](int index) const { return reinterpret_cast<const uint8*>(this)[index]; }
        inline bool operator==( const CellCoords &coords ) const { return 0 == ( (i-coords.i)|(j-coords.j)|(k-coords.k) ); }
    };
    struct CellData
    {
        uint32 m_FirstEntryIdx;
        uint32 m_NumEntries;
    };

    // Compute discrete cell coords, clamps outside points
    inline CellCoords Discretize( const Vec3& pos ) const
    {
        CellCoords coords;
        Vec3 lPosRel( pos - m_PosMin );
        Vec3 lNormalizedCoords( m_RcpCellSizes[0] * lPosRel[0],
                                m_RcpCellSizes[1] * lPosRel[1],
                                m_RcpCellSizes[2] * lPosRel[2] );
        coords.i = mal::Clamp<uint8>( mal::IntPart<uint8>( lNormalizedCoords[0] ), 0, m_Dimensions[0]-1 );
        coords.j = mal::Clamp<uint8>( mal::IntPart<uint8>( lNormalizedCoords[1] ), 0, m_Dimensions[1]-1 );
        coords.k = mal::Clamp<uint8>( mal::IntPart<uint8>( lNormalizedCoords[2] ), 0, m_Dimensions[2]-1 );
        return coords;
    }

    // Hash valid cell coords into the flat CellData array
    inline uint32 Hash( const CellCoords& coords ) const
    {
        GEO_ASSERT( coords.i < m_Dimensions[0] && coords.j < m_Dimensions[1] && coords.k < m_Dimensions[2] );
        return uint32( coords.i*m_Dimensions[1]*m_Dimensions[2] + coords.j*m_Dimensions[2] + coords.k );
    }

public:
    inline GSimpleSpatialHash( entry_index_type max_entries, uint32 dimensions[3], const Vec3& pos_min, const Vec3& pos_max )
    : m_MaxEntries(max_entries), m_NumEntries(0)
    , m_NumCells( dimensions[0] * dimensions[1] * dimensions[2] )
    , m_PosMin(pos_min), m_PosMax(pos_max)
    {
        GEO_ASSERT( m_MaxEntries < entry_index_type(-1) );
        GEO_ASSERT( m_MaxEntries > 0 );
        GEO_ASSERT( m_NumCells > 0 );
        GEO_ASSERT( dimensions[0] < 256 && dimensions[1] < 256 && dimensions[2] < 256 );
        m_Dimensions[0] = dimensions[0];
        m_Dimensions[1] = dimensions[1];
        m_Dimensions[2] = dimensions[2];
        m_vecCells = new CellData[m_NumCells];
        m_vecEntries = new entry_type[m_MaxEntries];
        m_RcpCellSizes[0] = Real(dimensions[0]) / (m_PosMax[0] - m_PosMin[0]);
        m_RcpCellSizes[1] = Real(dimensions[1]) / (m_PosMax[1] - m_PosMin[1]);
        m_RcpCellSizes[2] = Real(dimensions[2]) / (m_PosMax[2] - m_PosMin[2]);
    }

    inline ~GSimpleSpatialHash()
    {
        if( m_vecCells ) delete[] m_vecEntries;
        if( m_vecEntries ) delete[] m_vecCells;
    }

    inline void Clear()
    {
        m_NumEntries = 0;
        memset( &m_vecCells[0], 0, m_NumCells*sizeof(CellData) );
    }

    //\name First-pass add to compute per-cell entry count
    //@{
    inline void BeginClassify()
    {
        Clear();
    }
    // Allocates a slot for the entry in its containing cell, but does *NOT* add it yet
    inline void Classify( const entry_type& entry, const Vec3& point )
    {
        GEO_ASSERT( m_NumEntries <= m_MaxEntries );
        // Count entries per cell
        CellData& cell( m_vecCells[Hash( Discretize( point ) )] );
        cell.m_NumEntries++;
        m_NumEntries++;
    }
    inline void EndClassify()
    {
        // Set per-cell subarray first indices and reset per-cell entry count to 0
        uint32 acc_num_entries(0);
        for( uint32 it_cell=0; it_cell < m_NumCells; it_cell++ )
        {
            CellData& cell(m_vecCells[it_cell]);
            cell.m_FirstEntryIdx = acc_num_entries;
            acc_num_entries += cell.m_NumEntries;
            cell.m_NumEntries = 0; //Reset numentries, as they are not yet added!
        }
        GEO_ASSERT( acc_num_entries == m_NumEntries );
    }
    //@}

    //\name Second-pass add that effectively adds to cells and returns potentially overlapping entries
    //@{
    inline void BeginTestAdd()
    {
    }
    /* Add entry to previously reserved cell in Classify(), no testing.
       - Useful to add "passive" entries that are only "overlapped-with" (eg: boundary geometry)
    */
    inline void Add( const entry_type& entry, const Vec3& point, Real radius )
    {
        CellData& cell( m_vecCells[ Hash( Discretize( point ) ) ] );
        m_vecEntries[ cell.m_FirstEntryIdx + cell.m_NumEntries ] = entry;
        cell.m_NumEntries++;
    }
    /* Test entry against SH contents first, and add it to its
       containing cell afterwards. This avoids self-overlap reporting.
       - Useful to add "active" entries that need all their overlaps to be reported (eg: actors)
    */
    inline bool TestAdd( const entry_type& entry, const Vec3& point, Real radius,
                         std::vector<entry_type>& vec_potential_overlaps )
    {
        bool bResult = TestNoAdd(point,radius,vec_potential_overlaps);
        Add(entry,point,radius);
        return bResult;
    }
    inline void EndTestAdd()
    {
    }
    //@}

    /* Test but do not add, used to query the SH without modifying it
       - Useful to query the GSSH with an arbitrary sphere after filling it
    */
    inline bool TestNoAdd( const Vec3& point, Real radius,
                           std::vector<entry_type>& vec_potential_overlaps )
    {
        // Gather potential overlaps
        CellCoords coords_min = Discretize( point - Vec3(radius) );
        CellCoords coords_max = Discretize( point + Vec3(radius) );
        vec_potential_overlaps.clear();
        // Fast-check if only 1 cell
        if( 0 == ( (coords_min.i - coords_max.i) | (coords_min.j - coords_max.j) | (coords_min.k - coords_max.k) ) )
            GetCellEntries( Hash(coords_min), vec_potential_overlaps );
        else // Fill cells from coords_min to coords_max in all 3 dimensions
            for( uint32 i=coords_min.i; i <= coords_max.i; i++ )
                for( uint32 j=coords_min.j; j <= coords_max.j; j++ )
                    for( uint32 k=coords_min.k; k <= coords_max.k; k++ )
                        GetCellEntries( Hash(CellCoords(i,j,k)), vec_potential_overlaps );
        return vec_potential_overlaps.size() > 0;
    }

private:
    // inline const entry_type& GetEntryByIdx( entry_index_type index ) const { GEO_ASSERT( index < m_NumEntries ); return m_vecEntries[index]; }
    inline void GetCellEntries( uint32 index, std::vector<entry_type>& vec_entries ) const
    {
        GEO_ASSERT( index < m_NumCells );
        const CellData& cell( m_vecCells[index] );
        for( uint32 it_eic=0; it_eic<cell.m_NumEntries; it_eic++ )
            vec_entries.push_back( m_vecEntries[cell.m_FirstEntryIdx + it_eic] );
    }
    //non-copiable
    GSimpleSpatialHash( const GSimpleSpatialHash& ssh ) {}

private:
    uint32 m_MaxEntries;
    uint32 m_NumEntries;
    uint32 m_Dimensions[3];
    uint32 m_NumCells;
    Vec3 m_PosMin;
    Vec3 m_PosMax;
    Vec3 m_RcpCellSizes;

    CellData* m_vecCells;
    entry_type* m_vecEntries;
};

} //namespace geo

#endif //GEO_UTIL_GSIMPLESPATIALHASH_H
