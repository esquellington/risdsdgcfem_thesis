#ifndef S2_MS_SPH_GRID_2D_H
#define S2_MS_SPH_GRID_2D_H

#include "Config.h"
#include <vector>

namespace S2 { namespace ms
{

//! Uniform 2D grid spatial particle classification class.
/*! This class implements a 2D uniform grid for fast spatial
  classification of 2D spheres.

  A sphere is classified in a SINGLE cells that contains its center.
  The sphere-cell intersction test is:
  Inside(center,cell) || Min( dist(center,cell_corner_i) ) <= radius

  Where dist(p,cell_corner) can be seen as the inclusion test of a point to
  the minkowski sum of a rectangle (cell) and a circle).
  
  \todo Convert it into an Infinite Hash Grid instead of a Finite Grid
  \todo Reserve fixed size cell-buckets represented by 2 pointers and an integer
  struct Bucket { T* m_Data; Bucket *m_Next; int size };
  (better GBucket<T,N> with N max size as a global type)

  Grid2D(x,y) => Array[j][i] mapping:
  
       Y=i      (pos_max)
       ^ .......
       | .......
       | .......
       | .......
       |-------> X=j
  (pos_min)  
*/
class SPH_Grid2D
{
public:
    typedef uint16 ParticleID;
    struct CellID { int i; int j; CellID() {} CellID(int ci, int cj):i(ci),j(cj){} };
    
public:
  
    SPH_Grid2D(): m_NumCells(0), m_NumParticles(0) {}
    ~SPH_Grid2D() { Clear(); }
    
    //! \name Construction and Modification
    //@{
    void Init( const Point2 &pos_min, const Point2 &pos_max,
               unsigned int dim_x, unsigned int dim_y )
    {
        m_PosMin = pos_min;
        m_PosMax = pos_max;
        m_DimX = dim_x;
        m_DimY = dim_y;
        m_NumCells = m_DimX*m_DimY;

        m_CellSizes = m_PosMax-m_PosMin;
        m_CellSizes.x() /= m_DimX;
        m_CellSizes.y() /= m_DimY;
        m_CellHalfSizes = 0.5*m_CellSizes;
        
        m_MinCellSize = mal::Min( m_CellSizes.x(), m_CellSizes.y() );
        m_vecCells.resize( m_NumCells ); //!< Use resize(), NOT reserve()        
    }
    
    void Clear()
    {
        for(unsigned it_cells=0; it_cells<m_NumCells; it_cells++)
            m_vecCells[it_cells].clear();
        m_NumParticles = 0;
    }

    //radius <= size/2
    CellID AddParticle( ParticleID pid, const Point2 &pos )
    {
        // Find cell that contains the sphere center
        CellID center_cell_id = GetCell( pos );
        if( !IsValid(center_cell_id) )
            return center_cell_id;
        
        AddToCell( center_cell_id, pid );

        m_NumParticles++;
        return center_cell_id;
    }
    //@}
    
    //!\name Grid level consultors
    //@{    
    unsigned int GetNeighbours( const Point2 &pos, Real radius,
                                const Point2 *vec_pos,
                                ParticleID *vec_neighbours ) const
    {
        MS_ASSERT(radius <= 0.5f*m_MinCellSize);
        
        unsigned int num_neighbours = 0;
        
        // Find cell that contains the sphere center
        CellID center_cell_id = GetCell( pos );
        if( !IsValid(center_cell_id) )
            return 0;

        Real sq_radius = radius*radius;
        num_neighbours += GetNeighboursInCell( center_cell_id, pos, sq_radius, vec_pos, vec_neighbours );
        
        // Check direction and distance from cell center
        Point2 cell_pos = GetPos( center_cell_id );
        Point2 cell_offset = pos - cell_pos;
        // X-axis
        if( cell_offset.x() > (m_CellHalfSizes.x()-radius)
            && center_cell_id.j+1 < m_DimX )
        {
            //check i,j+1 cell and i+,j+ corner
            num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i,center_cell_id.j+1),
                                                   pos, sq_radius, vec_pos,
                                                   &vec_neighbours[num_neighbours] );
            if( center_cell_id.i+1 < m_DimY
                && (pos - GetCornerPos<1,1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i+1,center_cell_id.j+1),
                                                       pos, sq_radius, vec_pos,
                                                       &vec_neighbours[num_neighbours] );
        }
        else if( cell_offset.x() < -(m_CellHalfSizes.x()-radius)
                 && center_cell_id.j-1 >= 0 )
        {
            //check i,j-1 cell and i-,j- corner
            num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i,center_cell_id.j-1),
                                                   pos, sq_radius, vec_pos,
                                                   &vec_neighbours[num_neighbours] );
            if( center_cell_id.i-1 >= 0
                && (pos - GetCornerPos<-1,-1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i-1,center_cell_id.j-1),
                                                       pos, sq_radius, vec_pos,
                                                       &vec_neighbours[num_neighbours] );
        }        
        // Y-axis
        if( cell_offset.y() > (m_CellHalfSizes.y()-radius)
            && center_cell_id.i+1 < m_DimY )
        {
            //check i+1,j cell and i+,j- corner
            num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i+1,center_cell_id.j),
                                                   pos, sq_radius, vec_pos,
                                                   &vec_neighbours[num_neighbours] );
            if( center_cell_id.j-1 >= 0
                && (pos - GetCornerPos<1,-1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i+1,center_cell_id.j-1),
                                                       pos, sq_radius, vec_pos,
                                                       &vec_neighbours[num_neighbours] );
        }        
        else if( cell_offset.y() < -(m_CellHalfSizes.y()-radius)
                 && center_cell_id.i-1 >= 0 )
        {
            //check i-1,j cell and i-,j+ corner
            num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i-1,center_cell_id.j),
                                                   pos, sq_radius, vec_pos,
                                                   &vec_neighbours[num_neighbours] );
            if( center_cell_id.j+1 < m_DimX
                && (pos - GetCornerPos<-1,1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell( CellID(center_cell_id.i-1,center_cell_id.j+1),
                                                       pos, sq_radius, vec_pos,
                                                       &vec_neighbours[num_neighbours] );
        }
        return num_neighbours;
    }

    unsigned int GetNeighbours2( const Point2 &pos, Real radius,
                                 const Point2 *vec_pos,
                                 ParticleID *vec_neighbours,
                                 float *vec_distances ) const
    {
        MS_ASSERT(radius <= 0.5f*m_MinCellSize);
        
        unsigned int num_neighbours = 0;
        
        // Find cell that contains the sphere center
        CellID center_cell_id = GetCell( pos );
        if( !IsValid(center_cell_id) )
            return 0;

        Real sq_radius = radius*radius;
        num_neighbours += GetNeighboursInCell2( center_cell_id, pos, sq_radius, vec_pos,
                                                vec_neighbours,  vec_distances );
        
        // Check direction and distance from cell center
        Point2 cell_pos = GetPos( center_cell_id );
        Point2 cell_offset = pos - cell_pos;
        // X-axis
        if( cell_offset.x() > (m_CellHalfSizes.x()-radius)
            && center_cell_id.j+1 < m_DimX )
        {
            //check i,j+1 cell and i+,j+ corner
            num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i,center_cell_id.j+1),
                                                    pos, sq_radius, vec_pos,
                                                    &vec_neighbours[num_neighbours],
                                                    &vec_distances[num_neighbours] );
            if( center_cell_id.i+1 < m_DimY
                && (pos - GetCornerPos<1,1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i+1,center_cell_id.j+1),
                                                        pos, sq_radius, vec_pos,
                                                        &vec_neighbours[num_neighbours],
                                                        &vec_distances[num_neighbours] );
        }
        else if( cell_offset.x() < -(m_CellHalfSizes.x()-radius)
                 && center_cell_id.j-1 >= 0 )
        {
            //check i,j-1 cell and i-,j- corner
            num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i,center_cell_id.j-1),
                                                    pos, sq_radius, vec_pos,
                                                    &vec_neighbours[num_neighbours],
                                                    &vec_distances[num_neighbours] );
            if( center_cell_id.i-1 >= 0
                && (pos - GetCornerPos<-1,-1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i-1,center_cell_id.j-1),
                                                        pos, sq_radius, vec_pos,
                                                        &vec_neighbours[num_neighbours],
                                                        &vec_distances[num_neighbours] );
        }        
        // Y-axis
        if( cell_offset.y() > (m_CellHalfSizes.y()-radius)
            && center_cell_id.i+1 < m_DimY )
        {
            //check i+1,j cell and i+,j- corner
            num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i+1,center_cell_id.j),
                                                    pos, sq_radius, vec_pos,
                                                    &vec_neighbours[num_neighbours],
                                                    &vec_distances[num_neighbours] );
            if( center_cell_id.j-1 >= 0
                && (pos - GetCornerPos<1,-1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i+1,center_cell_id.j-1),
                                                        pos, sq_radius, vec_pos,
                                                        &vec_neighbours[num_neighbours],
                                                        &vec_distances[num_neighbours] );
        }        
        else if( cell_offset.y() < -(m_CellHalfSizes.y()-radius)
                 && center_cell_id.i-1 >= 0 )
        {
            //check i-1,j cell and i-,j+ corner
            num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i-1,center_cell_id.j),
                                                    pos, sq_radius, vec_pos,
                                                    &vec_neighbours[num_neighbours],
                                                    &vec_distances[num_neighbours] );
            if( center_cell_id.j+1 < m_DimX
                && (pos - GetCornerPos<-1,1>(center_cell_id)).NormSq() <= sq_radius )
                num_neighbours += GetNeighboursInCell2( CellID(center_cell_id.i-1,center_cell_id.j+1),
                                                        pos, sq_radius, vec_pos,
                                                        &vec_neighbours[num_neighbours],
                                                        &vec_distances[num_neighbours] );
        }
        return num_neighbours;
    }


    /*! Gathers "all" (up to max_count) particles inside a sphere (pos,radius). */
    int GetParticlesInSphere( const Point2 &pos, Real radius,
                              const Point2 *vec_pos,
                              ParticleID *vec_pid, unsigned int max_count ) const
    {        
        // Find cell that contains the sphere center
        // We get the minimum and maximun indices of the cells in each
        // dimension.        
        CellID min_cell_id = GetCell( pos - Vec2(radius,radius) );
        CellID max_cell_id = GetCell( pos + Vec2(radius,radius) );

        // Clamp indices / Clip rectangle against grid bounds
        min_cell_id.i = mal::Clamp( min_cell_id.i, 0, m_DimY-1 );
        min_cell_id.j = mal::Clamp( min_cell_id.j, 0, m_DimX-1 );
        max_cell_id.i = mal::Clamp( max_cell_id.i, 0, m_DimY-1 );
        max_cell_id.j = mal::Clamp( max_cell_id.j, 0, m_DimX-1 );

        // Iterate over all cells in rectangle collecting particles
        // inside sphere
        unsigned int num_neighbours = 0;
        for( int i = min_cell_id.i; i <= max_cell_id.i; i++ )
            for( int j = min_cell_id.j; j <= max_cell_id.j; j++ )
            {
                num_neighbours += GetNeighboursInCell_MaxCount( CellID(i,j),
                                                                pos, mal::Sq(radius),
                                                                vec_pos,
                                                                &vec_pid[num_neighbours], max_count-num_neighbours );
                // Early exit if max_count reached
                if( num_neighbours == max_count ) return max_count;
            }
        return num_neighbours;
    }
    
    unsigned int GetNeighboursInCell( CellID cell_id,
                                      const Point2 &pos, Real sq_radius,
                                      const Point2 *vec_pos,
                                      ParticleID *vec_neighbours ) const
    {
        unsigned int num_neighbours = 0;
        const std::vector<ParticleID> &vec_particles_in_cell = m_vecCells[Id2Index(cell_id)];
        for( unsigned int it_particle=0; it_particle < vec_particles_in_cell.size(); it_particle++ )
            if( (pos-vec_pos[ vec_particles_in_cell[it_particle] ]).NormSq() < sq_radius )
                vec_neighbours[num_neighbours++] = vec_particles_in_cell[it_particle];
        return num_neighbours;
    }

    unsigned int GetNeighboursInCell_MaxCount( CellID cell_id,
                                               const Point2 &pos, Real sq_radius,
                                               const Point2 *vec_pos,
                                               ParticleID *vec_neighbours, unsigned int max_count ) const
    {
        unsigned int num_neighbours = 0;
        const std::vector<ParticleID> &vec_particles_in_cell = m_vecCells[Id2Index(cell_id)];
        for( unsigned int it_particle=0; it_particle < vec_particles_in_cell.size(); it_particle++ )
        {
            if( (pos-vec_pos[ vec_particles_in_cell[it_particle] ]).NormSq() < sq_radius )
                vec_neighbours[num_neighbours++] = vec_particles_in_cell[it_particle];
            // Early exit if max_count reached
            if( num_neighbours == max_count ) return max_count;
        }
        return num_neighbours;
    }
    
    unsigned int GetNeighboursInCell2( CellID cell_id,
                                       const Point2 &pos, Real sq_radius,
                                       const Point2 *vec_pos,
                                       ParticleID *vec_neighbours,
                                       float *vec_distances ) const
    {
        unsigned int num_neighbours = 0;
        const std::vector<ParticleID> &vec_particles_in_cell = m_vecCells[Id2Index(cell_id)];
        for( unsigned int it_particle=0; it_particle < vec_particles_in_cell.size(); it_particle++ )
        {
            float sq_dist = (pos-vec_pos[ vec_particles_in_cell[it_particle] ]).NormSq();
            if( sq_dist < sq_radius )
            {
                vec_neighbours[num_neighbours] = vec_particles_in_cell[it_particle];
                vec_distances[num_neighbours] = mal::Sqrt(sq_dist);
                num_neighbours++;
            }
        }
        return num_neighbours;
    }    
    
    unsigned int GetNumCells() const { return m_NumCells; }
    unsigned int GetNumParticles() const { return m_NumParticles; }
    void GetDimensions( unsigned int &dim_x, unsigned int &dim_y ) const { dim_x = m_DimX; dim_y = m_DimY; }
        
    bool IsValid( int index ) const { return (index >= 0 && index < (int)m_NumCells); }
    bool IsValid( CellID id ) const { return (id.i >= 0 && id.i < m_DimY && id.j >= 0 && id.j < m_DimX); }
    
    unsigned int GetNonEmptyCells( std::vector<int> &vec_cells )
    {
        vec_cells.clear();
        for(unsigned it_cells=0; it_cells<m_NumCells; it_cells++)
            if( m_vecCells[it_cells].size() > 0 )
                vec_cells.push_back( it_cells );
        return vec_cells.size();
    }

    unsigned int GetNumNonEmptyCells()
    {
        unsigned int num_nec = 0;
        for(unsigned it_cells=0; it_cells<m_NumCells; it_cells++)
            if( m_vecCells[it_cells].size() > 0) 
                num_nec++;
        return num_nec;
    }
    
    void GetAABB( Point2 &pos_min, Point2 &pos_max ) const { pos_min = m_PosMin; pos_max = m_PosMax; }

    /*! Computes statistics:
      \param max_iic: Maximum particles in any cell
      \param avg_iic: Avg particles per non-empty cell
    */
    void ComputeStatistics( unsigned int &max_iic, float &avg_iic, float &non_emtpy_cell_fraction )
    {
        unsigned int acc_iic = 0;
        max_iic = 0;
        for(unsigned it_cells=0; it_cells<m_NumCells; it_cells++)
        {
            unsigned int iic = m_vecCells[it_cells].size();
            acc_iic += iic;
            max_iic = mal::Max(max_iic,iic);
        }
        unsigned int num_nec = GetNumNonEmptyCells();        
        avg_iic = (num_nec > 0) ? float(acc_iic)/num_nec : 0.0f;
        non_emtpy_cell_fraction = float(num_nec)/m_NumCells;
    }    
    //@}
    
    //!\name Cell level consultors
    //@{
    //! Cell center
    Point2 GetPos( CellID id ) const { return Point2( m_PosMin.x() + m_CellSizes.x()*(float(id.j) + 0.5f),
                                                      m_PosMin.y() + m_CellSizes.y()*(float(id.i) + 0.5f) ); }
    Point2 GetPos( int index ) const { return GetPos( Index2Id(index)); }

    //! Cell corners <i,j> = <+-1,+-1>
    template <int DI, int DJ>
    Point2 GetCornerPos( CellID id ) const
    {
        return GetPos(id) + Vec2( m_CellSizes.x()*DJ, m_CellSizes.y()*DI );
    }
    
    Vec2 GetSizes( CellID id ) const { return m_CellSizes; }

    /*! If pos is on a cell boundary, the CellID returned will be the
      one given by floor() (be careful with the upper/right boundary,
      as the returned CellID will be out of bounds!! */    
    CellID GetCell( const Point2 &pos ) const
    {
        Point2 offset_pos_min = pos - m_PosMin;
        offset_pos_min.x() /= m_CellSizes.x();
        offset_pos_min.y() /= m_CellSizes.y();
        return CellID( (int)mal::Floor(offset_pos_min.y()), (int)mal::Floor(offset_pos_min.x()) );    
    }
    
    const std::vector<ParticleID> &GetParticles( int index ) const { return m_vecCells[index]; }
    const std::vector<ParticleID> &GetParticles( CellID id ) const { return GetParticles( Id2Index(id) ); }
    //@}
    
    //!\name Coordinates and Index conversions between the grid and the world.
    int Id2Index( CellID id ) const { return id.i*m_DimX + id.j; }
    CellID Index2Id( int index ) const { return CellID( index / m_DimX, index % m_DimX ); }
    //@}

private:
    void AddToCell( CellID id, ParticleID pid )
    {
        m_vecCells[ Id2Index(id) ].push_back(pid);
        //NO: m_NumParticles++;
    }
    
private:
    Point2 m_PosMin, m_PosMax;
    int m_DimX, m_DimY;
    Vec2 m_CellSizes;
    Vec2 m_CellHalfSizes;
    Real m_MinCellSize;
    
    unsigned int m_NumCells;
    unsigned int m_NumParticles;
    
    std::vector< std::vector<ParticleID> > m_vecCells;
};

} } //namespace S2::ms

#endif // S2_MS_SPH_GRID_2D_H
