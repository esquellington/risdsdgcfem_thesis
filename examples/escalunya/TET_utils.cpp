#include "TET_utils.h"
#include <Geo/np/Overlap.h>
#include <Geo/np/Clip.h>

#define __USE_OVERLAP_TRITET
#define __USE_OVERLAP_TRITET_BVH //3-4x faster even with counting BVH BottomUp construction

namespace geo
{

//----------------------------------------------------------------
// Tetrahedron stuff
// \todo This is LINEAR geom, not bv nor shape... it's more
//       fundamental. GTriangle, GLine, GRay, GPlane, GHalfSpace could be
//       together in some sub-dir
//----------------------------------------------------------------
class Tetrahedron //\todo GSimplex<D>??
{
public:
    enum EConstants { cDimension = 3 };
    typedef Vec3 vec_type;
public:
    finline Tetrahedron() {}
    finline Tetrahedron( const vec_type& r0, const vec_type& r1, const vec_type& r2, const vec_type& r3 )
        { m_vec_r[0]=r0; m_vec_r[1]=r1; m_vec_r[2]=r2; m_vec_r[3]=r3; }
    finline const vec_type& r( unsigned int i ) const { return m_vec_r[i]; }
public:
    vec_type m_vec_r[4];
};

finline Real ComputeVolume( const Tetrahedron& tet )
{
    return mal::Det( mal::GMat3x3_From_Columns( tet.r(1)-tet.r(0), tet.r(2)-tet.r(0), tet.r(3)-tet.r(0) ) );
}

/* Circumsphere that passes through all 4 vertices
  \see http://www.geometrictools.com/Source/Containment.html
  \see http://www.geometrictools.com/Documentation/CentersOfSimplex.pdf
  \todo Consider degenerate cases (tet/tri/segment/point) and return
        circumsphere/circumcircle/midpoint/point accordingly
*/
bv::Sphere3 ComputeCircumsphere( const Tetrahedron& tet )
{
    Real volume( ComputeVolume(tet) );
    Real rcp_volume_12( mal::Rcp<Real>(12*volume) );
    GEO_ASSERT( !mal::IsNaN(rcp_volume_12) );
    Vec3 e1( tet.r(1)-tet.r(0) );
    Vec3 e2( tet.r(2)-tet.r(0) );
    Vec3 e3( tet.r(3)-tet.r(0) );
    Vec3 L2( mal::NormSq(e1), mal::NormSq(e2), mal::NormSq(e3) );
    Vec3 c( tet.r(0).x() + rcp_volume_12*( +(e2.y()*e3.z()-e3.y()*e2.z())*L2.x() - (e1.y()*e3.z()-e3.y()*e1.z())*L2.y() + (e1.y()*e2.z()-e2.y()*e1.z())*L2.z() ),
            tet.r(0).y() + rcp_volume_12*( -(e2.x()*e3.z()-e3.x()*e2.z())*L2.x() + (e1.x()*e3.z()-e3.x()*e1.z())*L2.y() - (e1.x()*e2.z()-e2.x()*e1.z())*L2.z() ),
            tet.r(0).z() + rcp_volume_12*( +(e2.x()*e3.y()-e3.x()*e2.y())*L2.x() - (e1.x()*e3.y()-e3.x()*e1.y())*L2.y() + (e1.x()*e2.y()-e2.x()*e1.y())*L2.z() ) );
    return bv::Sphere3( c, mal::Norm(tet.r(0)-c) );
    /*\todo BY NOW WE COMPUTE bary/radius, but for correct ODT WE SHOULD USE http://mathworld.wolfram.com/Circumsphere.html, req Det(4x4)...
    Vec3 barycenter( Real(0.25) * (tet.r(0)+tet.r(1)+tet.r(2)+tet.r(3)) );
    return bv::Sphere3( barycenter, mal::Norm(tet.r(0)-barycenter) );
    */
}

/* Tetrahedralize an AABB using a simple regular lattice + minimal hexahedron tetrahedralization.
   \todo Find the "goodness" of this tets (2 different fundamental shapes, I think)
*/
bool Make_TetSolidShape3_From_AABB3_Minimal( const geo::Vec3& pos_min, const geo::Vec3& pos_max,
                                             geo::Real cell_size,
                                             geo::EditableTetSolidShape3& solid )
{
    GEO_ASSERT( !mal::IsNaN(pos_min) && !mal::IsNaN(pos_max) && !mal::IsNaN(cell_size) );
    // Compute grid resolution for the largest extent
    const uint32 cMaxDimension = 32; //\todo publish param
    Vec3 sizes( pos_max-pos_min );
    uint32 dimensions[3] = {3,2,2}; //\note dimension[d] = #vertices[d]
    for( int i=0; i<3; i++ )
        dimensions[i] = mal::Clamp<uint32>( mal::IntPart<uint32>( mal::Ceil( sizes[i] / cell_size ) ) + 1, 2, cMaxDimension );
    // GEO_LOG( "Make_TetSolidShape3_From_AABB3_Simple() lattice dim = (%d,%d,%d)", dimensions[0], dimensions[1], dimensions[2] );
    solid.Clear();
    geo::Make_TetSolidShape3_Box( solid, Real(0.5)*sizes, dimensions[0], dimensions[1], dimensions[2] );
    return true;
}

/* Tetrahedralize an AABB using BCC lattice
   \see "Isosurface Stuffing: Fast Tetrahedral Meshes with Good Dihedral Angles"
   \todo Find the "goodness" of this tets (a single fundamental shape, I think)
*/
bool Make_TetSolidShape3_From_AABB3_BCC( const geo::Vec3& pos_min, const geo::Vec3& pos_max,
                                         geo::Real cell_size,
                                         //\todo Other params...?
                                         geo::EditableTetSolidShape3& solid )
{
    //\todo
    return false;
}

/* Fit a given TetSolid to a contained TriSurface:
   - Remove external tet
   - Adjust boundary tet
   - Adjust internal tet
*/
bool Fit_TetSolidShape3_To_TriSurfaceShape3_OLD( const geo::TriSurfaceShape3& surface,
                                                 //\todo float cell_ratio, float cell_size,
                                                 geo::EditableTetSolidShape3& solid )
{
#ifdef __USE_OVERLAP_TRITET_BVH
    const BVH_TriSurfaceShape3* pSurfaceBVH( surface.GetBVH() ); //\todo THIS is assumed up-to-date with surface object( transform, dof )
    GEO_ASSERT( 0 != pSurfaceBVH );
#endif

    //---- Remove completely external tet
    /* Classify vertices as in/out
       \todo OPTIMIZATION:
       - Compute surface BVH
       - "Rasterize" using axis-aligned raycasts between aabb-surface nodes, tracking in/out status
    */
    std::vector<bool> vec_outside_vid;
    vec_outside_vid.reserve( solid.GetNumV() );
    for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
    {
#ifdef __USE_OVERLAP_TRITET_BVH_DISABLED_BECAUSE_REJECTS_COMPLETELY_INTERIOR_TET //TEMP ACTUALLY this accepts only near-surface vertices... just testing
        std::vector< BVH_TriSurfaceShape3::entry_index_type > vecOverlaps;
        if( pSurfaceBVH->Test( BVH_TriSurfaceShape3::bv_type(solid.V_Pos_0(it_v)), vecOverlaps ) )
            vec_outside_vid.push_back( geo::np::Overlap_Point3_TriSurface3( solid.V_Pos_0(it_v),
                                                                            &surface, geo::Transform3::Identity(), surface.GetVecDefaultDOF() ) );
        else
            vec_outside_vid.push_back(false);
#else
        vec_outside_vid.push_back( geo::np::Overlap_Point3_TriSurface3( solid.V_Pos_0(it_v),
                                                                        &surface, geo::Transform3::Identity(), surface.GetVecDefaultDOF() ) );
#endif
    }
    // Gather outside tets
    std::vector<uint32> vec_outside_tid;
    for( unsigned int it_tet=0; it_tet<solid.GetNumT(); it_tet++ )
    {
        bool bInside(false);
        //-- Exclude tets with any inside vid
        for( unsigned int it_vit=0; !bInside && it_vit<4 ; it_vit++ )
            bInside = vec_outside_vid[ solid.T_VID(it_tet,it_vit) ];
        // if( bInside ) GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() tet[%d] has interior vid", it_tet );

#ifdef __USE_OVERLAP_TRITET
        /*\todo OPTIMIZATION:
          - Instead of testing tet and tri edges explicitly:
            - Test whole tets against surface BVH
            - If overlap, test Overlap_Tri_Tet (including all edge tests and potential vertex containment)
            \todo TriVsTet will ALSO be useful in DCR computation, to gather all tri that must be partitioned into an element
        */
#  ifdef __USE_OVERLAP_TRITET_BVH
        std::vector< BVH_TriSurfaceShape3::entry_index_type > vecOverlaps;
        geo::Vec3 vec_tet_pos[4] = { solid.V_Pos_0( solid.T_VID(it_tet,0) ),
                                     solid.V_Pos_0( solid.T_VID(it_tet,1) ),
                                     solid.V_Pos_0( solid.T_VID(it_tet,2) ),
                                     solid.V_Pos_0( solid.T_VID(it_tet,3) ) };
        auto tetBV = BVH_TriSurfaceShape3::bv_type(vec_tet_pos[0]).Merge(vec_tet_pos[1]).Merge(vec_tet_pos[2]).Merge(vec_tet_pos[3]);
        if( pSurfaceBVH->Test( tetBV, vecOverlaps ) )
        {
            for( unsigned int it_o=0; !bInside && it_o<vecOverlaps.size(); it_o++ )
            {
                uint32 tid = vecOverlaps[it_o];
                bInside = geo::np::Overlap_Triangle3_Tetrahedron3( surface.V_Pos_0( surface.T_VID(tid,0) ),
                                                                   surface.V_Pos_0( surface.T_VID(tid,1) ),
                                                                   surface.V_Pos_0( surface.T_VID(tid,2) ),
                                                                   vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] );
            }
        }
#  else
        for( unsigned int it_tri=0; !bInside && it_tri<surface.GetNumT(); it_tri++ )
            bInside = geo::np::Overlap_Triangle3_Tetrahedron3( surface.V_Pos_0( surface.T_VID(it_tri,0) ),
                                                               surface.V_Pos_0( surface.T_VID(it_tri,1) ),
                                                               surface.V_Pos_0( surface.T_VID(it_tri,2) ),
                                                               solid.V_Pos_0( solid.T_VID(it_tet,0) ),
                                                               solid.V_Pos_0( solid.T_VID(it_tet,1) ),
                                                               solid.V_Pos_0( solid.T_VID(it_tet,2) ),
                                                               solid.V_Pos_0( solid.T_VID(it_tet,3) ) );
#  endif
#else //TEMP: DEPRECATED CODE, kept as reference form a while...
        //-- Exclude tets any edge crossing the surface \todo EACH edge is tested multiple times... also consider BVH!!
        for( unsigned int it_eit=0; !bInside && it_eit<6 ; it_eit++ )
        {
            geo::Vec3 p0,p1;
            solid.T_Edge_0( it_tet, it_eit, p0, p1 ); //\todo Enumerate 6 edges (0->1,0->2,0->3) + (1->2->3->1)in a tet, get endpoints
            geo::np::RayHit3 rh;
            //todo bOutside = !geo::np::RayCast_TriSurfaceShape3( surface, p0, p1-p0, geo::Interval(0,1), rh );
            for( unsigned int it_tri=0; !bInside && it_tri<surface.GetNumT(); it_tri++ )
                bInside = geo::np::RayCast_Triangle3_DoubleSided( p0, p1-p0, Interval(0,1),
                                                                  surface.V_Pos_0( surface.T_VID(it_tri,0) ),
                                                                  surface.V_Pos_0( surface.T_VID(it_tri,1) ),
                                                                  surface.V_Pos_0( surface.T_VID(it_tri,2) ),
                                                                  rh );
            // if( bInside ) GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() tet[%d] edges pierces the surface", it_tet );
        }
        //-- Exclude tets that contain any surface vertices or are have any crossing surface edges
        if( !bInside )
        {
            geo::Vec3 vec_node_pos[4] = { };
            // for each face in the tetrahedron
            for( unsigned int it_fit=0; !bInside && it_fit<3 ; it_fit++ )
            {
                // for each triangle in the surface
                for( unsigned int it_tri=0; !bInside && it_tri<surface.GetNumT(); it_tri++ )
                {
                    // for each edge in the triangle
                    for( unsigned int it_eit=0; !bInside && it_eit<3 ; it_eit++ )
                    {
                        geo::Vec3 p0,p1;
                        surface.T_Edge_0( it_tri, it_eit, p0, p1 );
                        // Test p0 inside tet
                        bInside = geo::np::Overlap_Point3_Tetrahedron3( p0,
                                                                        solid.V_Pos_0( solid.T_VID(it_tet,0) ),
                                                                        solid.V_Pos_0( solid.T_VID(it_tet,1) ),
                                                                        solid.V_Pos_0( solid.T_VID(it_tet,2) ),
                                                                        solid.V_Pos_0( solid.T_VID(it_tet,3) ) );
                        // if( !!bInside ) GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() tet[%d] contains surface vertices", it_tet );
                        // Test tri-edge against tet-face
                        if( !bInside )
                        {
                            geo::np::RayHit3 rh;
                            bInside = geo::np::RayCast_Triangle3_DoubleSided( p0, p1-p0, Interval(0,1),
                                                                              solid.V_Pos_0( solid.T_VID(it_tet,(it_fit+1)%4) ),
                                                                              solid.V_Pos_0( solid.T_VID(it_tet,(it_fit+2)%4) ),
                                                                              solid.V_Pos_0( solid.T_VID(it_tet,(it_fit+3)%4) ),
                                                                              rh );
                        }
                    }
                }
                // if( bInside ) GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() tet[%d] is pierced by surface edges", it_tet );
            }
        }
#endif
        //-- Remaining POT area actually outside
        if( !bInside )
            vec_outside_tid.push_back( it_tet );
    }
    // Remove tets (\todo MUST ALSO REMOVE ISOLATED nodes, AUTOMATICALLY, in rebuildtopology)
    if( !vec_outside_tid.empty() )
    {
        solid.RemoveTetrahedrons( vec_outside_tid.size(), &vec_outside_tid[0] );
        GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() discarded %d exterior tet, remaining %d interior|boundary", (int)vec_outside_tid.size(), solid.GetNumT() );
    }
    //\todo adjust boundary and internal tets according to cell_size/cell_ratio/etc...
    return true;
}

/* Fit a given TetSolid to a contained TriSurface:
   Classify:
   - External: all out vertices, no piercing edges
   - Internal: all in vertices, no piercing edges
   - Boundary: both in && out nodes or piercing edges
   Fit:
   - Remove External
   - Adjust Boundary and Internal
*/
bool Fit_TetSolidShape3_To_TriSurfaceShape3( const geo::TriSurfaceShape3& surface,
                                             float odt_relaxation_coeff, float lpc_relaxation_coeff, unsigned int max_iter, bool b_fix_nmf,
                                             geo::EditableTetSolidShape3& solid )
{
    const BVH_TriSurfaceShape3* pSurfaceBVH( surface.GetBVH() ); //\todo THIS is assumed up-to-date with surface object( transform, dof )
    GEO_ASSERT( 0 != pSurfaceBVH );

    /* Classify vertices as in/out
       \todo OPTIMIZATION:
       - "Rasterize" using axis-aligned raycasts (with BVH) between aabb-surface nodes, tracking in/out status
    */
    std::vector<bool> vec_is_inside_vid;
    vec_is_inside_vid.reserve( solid.GetNumV() );
    for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        vec_is_inside_vid.push_back( geo::np::Overlap_Point3_TriSurface3( solid.V_Pos_0(it_v),
                                                                          &surface, geo::Transform3::Identity(), surface.GetVecDefaultDOF() ) );
    //---- Classify tets
    std::vector<uint32> vec_external_tid;
    std::vector<uint32> vec_internal_tid;
    std::vector<uint32> vec_boundary_tid;
    for( unsigned int it_tet=0; it_tet<solid.GetNumT(); it_tet++ )
    {
        //-- Test tet vid inside surface
        bool bInternal(true);
        for( unsigned int it_vit=0; bInternal && it_vit<4; it_vit++ )
            bInternal = vec_is_inside_vid[ solid.T_VID(it_tet,it_vit) ];

        // Test tet against surface BVH and, if overlap, test against all tri, which considers all piercing-edge and contained-tri cases.
        bool bPiercing(false);
        std::vector< BVH_TriSurfaceShape3::entry_index_type > vecOverlaps;
        geo::Vec3 vec_tet_pos[4] = { solid.V_Pos_0( solid.T_VID(it_tet,0) ),
                                     solid.V_Pos_0( solid.T_VID(it_tet,1) ),
                                     solid.V_Pos_0( solid.T_VID(it_tet,2) ),
                                     solid.V_Pos_0( solid.T_VID(it_tet,3) ) };
        auto tetBV = BVH_TriSurfaceShape3::bv_type(vec_tet_pos[0]).Merge(vec_tet_pos[1]).Merge(vec_tet_pos[2]).Merge(vec_tet_pos[3]);
        tetBV.Extend(0.01f);
        if( pSurfaceBVH->Test( tetBV, vecOverlaps ) )
        {
            for( unsigned int it_o=0; !bPiercing && it_o<vecOverlaps.size(); it_o++ )
            {
                uint32 tid = vecOverlaps[it_o];
                bPiercing = geo::np::Overlap_Triangle3_Tetrahedron3( surface.V_Pos_0( surface.T_VID(tid,0) ),
                                                                     surface.V_Pos_0( surface.T_VID(tid,1) ),
                                                                     surface.V_Pos_0( surface.T_VID(tid,2) ),
                                                                     vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] );
            }
        }
        //-- Classify tet
        if( bPiercing )
            vec_boundary_tid.push_back( it_tet );
        else if( bInternal )
            vec_internal_tid.push_back( it_tet );
        else
            vec_external_tid.push_back( it_tet );
    }
    // Remove tets (and any resulting isolated nodes)
    if( !vec_external_tid.empty() )
    {
        solid.RemoveTetrahedrons( vec_external_tid.size(), &vec_external_tid[0] );
        GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() discarded %d external tet, remaining %d tet (%d interior and %d boundary)",
                (int)vec_external_tid.size(), solid.GetNumT(), (int)vec_internal_tid.size(), (int)vec_boundary_tid.size() );
    }

    /*Remove topology connections in solid that are not present in surface:
      - Non-manifold Vertices: Split external V shared by several non-adjacent boundary tet (split 1V into nV, 1 per adjacent connected set of tet)
        - Better, ANY non-manifold V should be split, not just external
      - Non-manifold Edges: Split external E shared by non-adjacent boundary tet (split 1E into nE, 1 per adjacent connected set of tet)
        - Better, ANY non-manifold E should be split, not just external.
        - \note A non-manifold E CANNOT pierce the surface, as in such case ALL its original adjacent tet would be present and it would not be non-manifold.
      - Faces: Split external F shared by adjacent boundary tet. (split 1F into 2F, one par adjacent tet, IMPORTANT: BUT NOT ALWAYS V need to be split...)
        - Actually, we should gather the whole connected surface (potentially non-manifold) formed by adjacent excess Faces and split it as a whole...
      - Elements: Split elements T that contain disconnected parts of the surface (=> 1 virtual tet for connected component)
      IMPORTANT: We MUST AVOID reconnecting split V/E/F/T in RebuildTopology()
      - Vertices are merged *geometrically*, thus we should MOVE them apart > epsilon
      - Elements are connected topologically (common face must have matching topology, VID), therefore won't be connected if V are not merged
      - Alternatively, add/pass an explicit list of do-not-merge VID pairs and consider it in RebuildTopology()
      \todo THIS TOPOLOGY-correction would be an improvement over "Automatic Construction of Coarse, High-Quality
            Tetrahedralizations that Enclose and Approximate Surfaces for Animation", eg, Homer's legs would be separated!
    */
    if( b_fix_nmf )
    {
        GEO_LOG("************ Fit_TetSolidShape3_To_TriSurfaceShape3() BEFORE FixNonManifoldFeatures() #V=%u #T=%u", solid.GetNumV(), solid.GetNumT() );
        // Gather split F array
        vec_is_inside_vid.clear();
        vec_is_inside_vid.reserve( solid.GetNumV() );
        for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
            vec_is_inside_vid.push_back( geo::np::Overlap_Point3_TriSurface3( solid.V_Pos_0(it_v),
                                                                              &surface, geo::Transform3::Identity(), surface.GetVecDefaultDOF() ) );
        std::vector< geo::tetsolid3_face_id_type > vec_split_fid; //todo consider scope alloc
        // Gather all non-boundary faces that are completely external, that is, have all V outside AND do not pierce the surface
        for( unsigned int it_tet=0; it_tet<solid.GetNumT(); it_tet++ )
        {
            for( unsigned int it_ntit=0; it_ntit<4; it_ntit++ )
            {
                uint32 ntid( solid.T_NTID(it_tet,it_ntit) );
                if( ntid != cTetSolid3_InvalidFeatureIndex
                    && it_tet < ntid )
                {
                    uint32 vec_fvid[3] = { solid.T_VID(it_tet,(it_ntit+1)%4),
                                           solid.T_VID(it_tet,(it_ntit+2)%4),
                                           solid.T_VID(it_tet,(it_ntit+3)%4) };
                    if( !vec_is_inside_vid[vec_fvid[0]]
                        && !vec_is_inside_vid[vec_fvid[1]]
                        && !vec_is_inside_vid[vec_fvid[2]]
                        && !geo::np::Overlap_Triangle3_TriSurface3( solid.V_Pos_0(vec_fvid[0]), solid.V_Pos_0(vec_fvid[1]), solid.V_Pos_0(vec_fvid[2]),
                                                                    &surface, Transform3::Identity(), surface.GetVecDefaultDOF() ) )
                    {
                        vec_split_fid.push_back( tetsolid3_face_id_type(it_tet,it_ntit) );
                        // Find opposite face-in-tet and add opposite fid
                        int ofit(-1);
                        for( unsigned int it_ofit=0; it_ofit<4; it_ofit++ )
                            if( solid.T_VID(ntid,it_ofit) != vec_fvid[0]
                                && solid.T_VID(ntid,it_ofit) != vec_fvid[1]
                                && solid.T_VID(ntid,it_ofit) != vec_fvid[2] )
                                ofit = it_ofit;
                        GEO_ASSERT(ofit != -1);
                        vec_split_fid.push_back( tetsolid3_face_id_type(ntid,ofit) );
                    }
                }
            }
        }
        GEO_LOG("************ Fit_TetSolidShape3_To_TriSurfaceShape3() found %d two-sided split-F", (int)vec_split_fid.size()/2 );
        solid.FixNonManifoldFeatures( vec_split_fid ); //Splits V,E and explicitly given F
        GEO_LOG("************ Fit_TetSolidShape3_To_TriSurfaceShape3() AFTER FixNonManifoldFeatures() #V=%u #T=%u", solid.GetNumV(), solid.GetNumT() );
    }

#define __ENABLE_FIT_TETSS_TRISS_ODT
// #define __TRACE_ODT
// #define __TRACE_ODT_VERBOSE
#ifdef __ENABLE_FIT_TETSS_TRISS_ODT

    /* IMPORTANT: Reclassify vertices as in/out, as their indices have
       changed from the original solid due to
       solid.RemoveTetrahedrons()
       \todo OPTIMIZATION:
       - "Rasterize" using axis-aligned raycasts (with BVH) between aabb-surface nodes, tracking in/out status
    */
    vec_is_inside_vid.clear();
    vec_is_inside_vid.reserve( solid.GetNumV() );
    for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        vec_is_inside_vid.push_back( geo::np::Overlap_Point3_TriSurface3( solid.V_Pos_0(it_v),
                                                                          &surface, geo::Transform3::Identity(), surface.GetVecDefaultDOF() ) );
//#define __TRACE_IS_INSIDE
#ifdef __TRACE_IS_INSIDE
    for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        GEO_LOG("V[%d] inside = %d",it_v,(int)vec_is_inside_vid[it_v]);
#endif

    //\todo adjust boundary and internal tets according to cell_size/cell_ratio/etc...
    // Gather V-tet and V-tri adjacency
    std::vector< std::vector<tetsolid3_feature_index_type> > vec_vtet( solid.GetNumV() );
    std::vector< std::vector<std::pair<tetsolid3_feature_index_type,tetsolid3_feature_index_type> > > vec_vtri( solid.GetNumV() );
    for( unsigned int it_tet=0; it_tet<solid.GetNumT(); it_tet++ )
    {
        for( unsigned int it_vit=0; it_vit<4; it_vit++ )
        {
            uint32 vid( solid.T_VID(it_tet,it_vit) );
            vec_vtet[vid].push_back(it_tet);
            // Gather pairs of other-VID that form fully internal/external triangles with a given VID in vec_vtri[vid]
            uint32 vec_other_vid[3] = { solid.T_VID(it_tet,(it_vit+1)%4),
                                        solid.T_VID(it_tet,(it_vit+2)%4),
                                        solid.T_VID(it_tet,(it_vit+3)%4) };
            for( unsigned int it_tit=0; it_tit<3; it_tit++ )
            {
                uint32 vid1( vec_other_vid[it_tit] );
                uint32 vid2( vec_other_vid[(it_tit+1)%3] );
                if( (vec_is_inside_vid[vid] && vec_is_inside_vid[vid1] && vec_is_inside_vid[vid2]) //internal \todo will appear TWICE, as they are shared by 2 tets...
                    || (!vec_is_inside_vid[vid] && !vec_is_inside_vid[vid1] && !vec_is_inside_vid[vid2]) ) //external
                {
                    /* Even for fully internal/external VID, the
                       actual triangle MAY intersect the surface (eg
                       1,2,7 in 1x1x1 cube test) therefore we MUST
                       test the whole tit against surface
                    */
                    if( !geo::np::Overlap_Triangle3_TriSurface3( solid.V_Pos_0(vid), solid.V_Pos_0(vid1), solid.V_Pos_0(vid2),
                                                                 &surface, Transform3::Identity(), surface.GetVecDefaultDOF() ) )
                        vec_vtri[vid].push_back( std::make_pair(vid1,vid2) );
                }
            }
        }
    }
#ifdef __TRACE_ODT_VERBOSE
    for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
    {
        GEO_LOG("VTRI[%u]----",it_v);
        for( auto vtri : vec_vtri[it_v] )
            GEO_LOG("->(%d,%d)",vtri.first,vtri.second);
    }
#endif
    // Gather vertices to be modified (layer0+layer1 | all)
    std::vector<Vec3> vec_pos0( solid.GetNumV() );
    for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        vec_pos0[it_v] = solid.V_Pos_0(it_v);
    const Real cSafetyMarginFactor(0.99f);
    unsigned int num_iter(0);
    bool bConverged(false);
    GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() starting layer[0] fitting with %d tet", solid.L_NumT(0));
    while( num_iter < max_iter && !bConverged )
    {
        GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() iteration %d", num_iter );
        /* Compute Optimal Delaunay Triangulation objective node positions for Layer[0] nodes
           \see "Automatic Construction of Coarse, High-Quality
           Tetrahedralizations that Enclose and Approximate Surfaces for
           Animation"
        */
        std::vector<Vec3> vec_pos1( solid.GetNumV(), Vec3::Zero() );
        std::vector<Real> vec_vol( solid.GetNumV(), 0 );
        for( unsigned int it_tet=solid.L_FirstTID(0); it_tet<solid.L_FirstTID(0)+solid.L_NumT(0); it_tet++ )
            //TEMP: consider fitting all if layer0-only result in degenerate tet
            //for( unsigned int it_tet=0; it_tet<solid.GetNumT(); it_tet++ )
        {
            Tetrahedron tetrahedron( vec_pos0[solid.T_VID(it_tet,0)],
                                     vec_pos0[solid.T_VID(it_tet,1)],
                                     vec_pos0[solid.T_VID(it_tet,2)],
                                     vec_pos0[solid.T_VID(it_tet,3)] );
            Real volume( geo::ComputeVolume(tetrahedron) );
            bv::Sphere3 circumsphere( geo::ComputeCircumsphere(tetrahedron) );
            for( unsigned int it_vit=0; it_vit<4 ; it_vit++ )
            {
                vec_pos1[ solid.T_VID(it_tet,it_vit) ] += volume*circumsphere.GetPos();
                vec_vol[ solid.T_VID(it_tet,it_vit) ] += volume;
            }
        }
        // Compute ODT final pos or copy previous (for non layer0 tet vertices)
        for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        {
            if( vec_vol[it_v] > 0 )
                vec_pos1[it_v] /= vec_vol[it_v];
            else
                vec_pos1[it_v] = vec_pos0[it_v];
        }

        // Simple Laplacian smoothing
        std::vector<Vec3> vec_smoothing_target( solid.GetNumV(), Vec3::Zero() );
        std::vector<uint32> vec_smoothing_count( solid.GetNumV(), 0 );
        for( unsigned int it_bf=0; it_bf<solid.GetNumBF(); it_bf++ )
        {
            // Acc BF effect on each V
            uint32 vid0( solid.BF_VID(it_bf,0) );
            uint32 vid1( solid.BF_VID(it_bf,1) );
            uint32 vid2( solid.BF_VID(it_bf,2) );
            vec_smoothing_target[vid0] += vec_pos0[vid1] + vec_pos0[vid2];
            vec_smoothing_target[vid1] += vec_pos0[vid0] + vec_pos0[vid2];
            vec_smoothing_target[vid2] += vec_pos0[vid0] + vec_pos0[vid1];
            vec_smoothing_count[vid0]++;
            vec_smoothing_count[vid1]++;
            vec_smoothing_count[vid2]++;
        }
        // divide by 2x size of 1-ring V set, because we've counted each neighbour-V twice in previous loop
        for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        {
            if( mal::NormSq(vec_smoothing_target[it_v]) > 0 )
                vec_smoothing_target[it_v] /= 2*vec_smoothing_count[it_v]; //\todo THIS IS WRONG -> 2*vec_vtri.size();
            else
                vec_smoothing_target[it_v] = vec_pos0[it_v];
        }

        // Move vertices to their ODT+smoothing target pos but stopping before any collision/collapse
        unsigned int num_unfulfilled = 0;
        for( unsigned int it_v=0; it_v<solid.GetNumV(); it_v++ )
        {
            // Compute maximum displacement
            //\todo Constrain d with per-vertex cached constraint plane normals and zero them (if still active will be rediscovered)
            Vec3 d( odt_relaxation_coeff * (vec_pos1[it_v] - vec_pos0[it_v]) );
            // if( num_iter % 2 == 0 )
                d += lpc_relaxation_coeff * (vec_smoothing_target[it_v] - vec_pos0[it_v]);
            // else //\todo Taubin smoothing, using \mu = -\lambda
            //     d -= 0.5*lpc_relaxation_coeff * (vec_smoothing_target[it_v] - vec_pos0[it_v]);
            if( mal::NormSq(d) == 0 ) continue;

#ifdef __TRACE_ODT
            GEO_LOG("V[%d]=(%f,%f,%f), d0=(%f,%f,%f), |d_0|=%f",
                    it_v,
                    vec_pos0[it_v][0], vec_pos0[it_v][1], vec_pos0[it_v][2],
                    d[0], d[1], d[2],
                    mal::Norm(d) );
#endif

            // Compute V-tri swept BV
            auto bv = BVH_TriSurfaceShape3::bv_type(vec_pos0[it_v]);
            bv.Merge( vec_pos0[it_v] + d );
            for( auto vtri : vec_vtri[it_v] )
            {
                bv.Merge( vec_pos0[vtri.first] );
                bv.Merge( vec_pos0[vtri.second] );
            }
            // Gather potential overlapping surface.tri for the BV
            std::vector< BVH_TriSurfaceShape3::entry_index_type > vecOverlaps;
            if( pSurfaceBVH->Test( bv, vecOverlaps ) )
            {
#ifdef __TRACE_ODT
                GEO_LOG("=> %d potential overlaps", (int)vecOverlaps.size() );
#endif
                // Perform ccd-toi-bisection iterations, testing each
                // V-tri swept tetrahedron against potential
                // surface.tri, and select left bracket if ANY actual
                // overlap
                Vec3 d0( Vec3::Zero() );
                Vec3 d1( 2*d ); //\note This makes original sampling mid-point to be actually the target d
                unsigned int num_iter_bisection(0);
                while( mal::NormSq(d1-d0) > 1e-6 && num_iter_bisection < 5 )
                {
                    // Test midpoint
                    d = 0.5*(d0+d1);
                    bool bPiercing(false);
                    for( unsigned int it_po=0; it_po<vecOverlaps.size() && !bPiercing; it_po++ )
                    {
                        trisurface3_feature_index_type tid( vecOverlaps[it_po] );
                        Vec3 vec_tri_pos[3] = { surface.V_Pos_0( surface.T_VID(tid,0) ),
                                                surface.V_Pos_0( surface.T_VID(tid,1) ),
                                                surface.V_Pos_0( surface.T_VID(tid,2) ) };
                        for( unsigned int it_vtri=0; it_vtri<vec_vtri[it_v].size() && !bPiercing; it_vtri++ )
                        {
                            Vec3 vec_tet_pos[4] = { vec_pos0[it_v],
                                                    vec_pos0[it_v]+d,
                                                    vec_pos0[vec_vtri[it_v][it_vtri].first],
                                                    vec_pos0[vec_vtri[it_v][it_vtri].second] };
                            // Test if non-degenerate tet (swept triangle), test it, otherwise test only final pos triangle.
                            if( mal::Abs(mal::Det(mal::GMat3x3_From_Columns( vec_tet_pos[1]-vec_tet_pos[0],
                                                                             vec_tet_pos[1]-vec_tet_pos[1],
                                                                             vec_tet_pos[1]-vec_tet_pos[2] ))) > 1e-6 )
                                bPiercing = geo::np::Overlap_Triangle3_Tetrahedron3( vec_tri_pos[0], vec_tri_pos[1], vec_tri_pos[2],
                                                                                     vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] );
                            else
                                bPiercing = geo::np::Overlap_Triangle3_Triangle3( vec_tri_pos[0], vec_tri_pos[1], vec_tri_pos[2],
                                                                                  vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3] );
                        }
                    }
                    // Detect collapse/inversion
                    bool bCollapse(false);
                    for( unsigned int it_tiv=0; it_tiv<vec_vtet[it_v].size() && !bCollapse; it_tiv++ )
                    {
                        uint32 tid( vec_vtet[it_v][it_tiv] );
                        Vec3 vec_tet_pos[4] = { vec_pos0[solid.T_VID(tid,0)],
                                                vec_pos0[solid.T_VID(tid,1)],
                                                vec_pos0[solid.T_VID(tid,2)],
                                                vec_pos0[solid.T_VID(tid,3)] };
                        // Compute volume before and after displacement d
                        Real vol0( mal::Det(mal::GMat3x3_From_Columns( vec_tet_pos[1] - vec_tet_pos[0],
                                                                       vec_tet_pos[2] - vec_tet_pos[0],
                                                                       vec_tet_pos[3] - vec_tet_pos[0] )) );
                        for( unsigned int it_vit=0; it_vit<4; it_vit++ )
                            if( it_v == solid.T_VID(tid,it_vit) )
                                vec_tet_pos[it_vit] += d;
                        Real vol1( mal::Det(mal::GMat3x3_From_Columns( vec_tet_pos[1] - vec_tet_pos[0],
                                                                       vec_tet_pos[2] - vec_tet_pos[0],
                                                                       vec_tet_pos[3] - vec_tet_pos[0] )) );
                        // Detect inversion using volume change fraction threshold (captures collapse and inversion/sign-change)
                        bCollapse = vol1/vol0 < 0.1;
                    }
                    // Choose bracket
                    if( bPiercing || bCollapse ) d1 = d;
                    else d0 = d;
                    num_iter_bisection++;
#ifdef __TRACE_ODT
                    GEO_LOG("=> |d|=%f", mal::Norm(d) );
#endif
                }
                d = d0;
                if(num_iter_bisection > 1) num_unfulfilled++;
#ifdef __TRACE_ODT
                GEO_LOG("=> %d bisection iterations", num_iter_bisection );
#endif
            }
            //Advance each V independently-in Gauss-Seidel fashion...
            vec_pos0[it_v] += d;
#ifdef __TRACE_ODT
            GEO_LOG("=> |d_1| = %f", mal::Norm(d) );
#endif
        }
        bConverged = num_unfulfilled == 0;
        num_iter++;
    }
    GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() fitting converged?=%s after %d iter", bConverged?"true":"false", num_iter );
    //\todo Change solid vertices to fitted positions
    solid.BeginEdition();
    {
        for( unsigned int it_v=0; it_v<vec_pos0.size(); it_v++ )
            solid.SetVertex( it_v, vec_pos0[it_v] );
    }
    solid.EndEdition();

    // Final report
    Real vol_solid( ComputeVolume(solid,Transform3::Identity(),0) );
    Real vol_surface( ComputeVolume(surface,Transform3::Identity(),0) );
    GEO_LOG("Fit_TetSolidShape3_To_TriSurfaceShape3() Vol(solid)=%f, Vol(surface)=%f, ratio=%f", vol_solid, vol_surface, vol_solid/vol_surface );
#endif //__ENABLE_FIT_TETSS_TRISS_ODT

    return true;
}

/* Clip a TriSurfaceShape3 to a TetSolidShape3 */
bool Clip_TriSurfaceShape3_TetSolidShape3( const geo::TriSurfaceShape3& surface,
                                           const geo::EditableTetSolidShape3& solid, const BVH_TetSolidShape3* pSolidBVH, //TEMP pass BVH as param by now
                                           geo::EditableTriSurfaceShape3& clipped_surface )
{
    // const BVH_TetSolidShape3* pSolidBVH( solid.GetBVH() ); //\todo THIS is assumed up-to-date with surface object( transform, dof )
//    GEO_ASSERT( 0 != pSolidBVH );

    GEO_LOG("Clip_TriSurfaceShape3_TetSolidShape3() surface with %d tri against solid with %d tet", surface.GetNumT(), solid.GetNumT() );
    unsigned int total_clipped_triangles(0);
    clipped_surface.Clear();
    clipped_surface.BeginEdition();
    for( unsigned int it_tri=0; it_tri<surface.GetNumT(); it_tri++ )
    {
        // Test tri against solid BVH and, if overlap, test against all tet
        geo::Vec3 vec_tri_pos[3] = { surface.V_Pos_0( surface.T_VID(it_tri,0) ),
                                     surface.V_Pos_0( surface.T_VID(it_tri,1) ),
                                     surface.V_Pos_0( surface.T_VID(it_tri,2) ) };
        std::vector< BVH_TetSolidShape3::entry_index_type > vecOverlaps;
        auto triBV = BVH_TetSolidShape3::bv_type(vec_tri_pos[0]).Merge(vec_tri_pos[1]).Merge(vec_tri_pos[2]);
        triBV.Extend(0.01f);
        bool bClipped(false);
        //TEMP: Hack to support no-BVH
        if( pSolidBVH ) pSolidBVH->Test( triBV, vecOverlaps );
        else for( uint32 it_tet=0; it_tet<solid.GetNumT(); it_tet++ ) vecOverlaps.push_back(it_tet);
        //END TEMP
        if( !vecOverlaps.empty() )
        {
            for( unsigned int it_o=0; it_o<vecOverlaps.size(); it_o++ )
            {
                uint32 tid = vecOverlaps[it_o];
                geo::Vec3 vec_clipped_points[6];
                uint32 num_clipped_v = geo::np::Clip_Triangle3_Tetrahedron3( vec_tri_pos[0], vec_tri_pos[1], vec_tri_pos[2],
                                                                             solid.V_Pos_0( solid.T_VID(tid,0) ),
                                                                             solid.V_Pos_0( solid.T_VID(tid,1) ),
                                                                             solid.V_Pos_0( solid.T_VID(tid,2) ),
                                                                             solid.V_Pos_0( solid.T_VID(tid,3) ),
                                                                             vec_clipped_points );
                // Clipped may be internal or piercing
                if( num_clipped_v > 0 )
                {
                    // Add all clipped points
                    uint32 vid0 = clipped_surface.AddVertex( vec_clipped_points[0] );
                    for( uint32 i=1; i<num_clipped_v; i++ )
                        clipped_surface.AddVertex( vec_clipped_points[i] );
                    // Add triangulated clipped polygon
                    for( uint32 i=1; i<num_clipped_v-1; i++ )
                        clipped_surface.AddTriangle( vid0, vid0+i, vid0+i+1 );
                    // Clipped/External/Internal
                    bClipped = true;
                }
            }
        }
        if( !bClipped )
        {
            // If not clipped, the triangle is COMPLETELY EXTERNAL... this should be an error!
            uint32 vid0 = clipped_surface.AddVertex( vec_tri_pos[0] );
            clipped_surface.AddVertex( vec_tri_pos[1] );
            clipped_surface.AddVertex( vec_tri_pos[2] );
            clipped_surface.AddTriangle( vid0, vid0+1, vid0+2 );
        }
        else total_clipped_triangles++;
    }
    clipped_surface.EndEdition();
    // GEO_LOG("Clip_TriSurfaceShape3_TetSolidShape3() finished with %d clipped and %d external triangles",
    //         total_clipped_triangles, surface.GetNumT() - total_clipped_triangles );
    return true;
}

} //namespace geo
