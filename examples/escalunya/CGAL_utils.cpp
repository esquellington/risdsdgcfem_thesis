#include "CGAL_utils.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//\todo This requires linking with libboost_system
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
// For features()
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>

//\note For surface mesh simplification
//#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

//#define __ENABLE_CGAL_UTILS_TRACE
#ifdef __ENABLE_CGAL_UTILS_TRACE
#  include <iostream>
#endif

namespace geo
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGAL_KERNEL;
//typedef CGAL::Simple_cartesian<double> CGAL_KERNEL; //\todo this seems faster but FAILS in refine_mesh_3() for bunny.off

// Polyhedron (all params after the kernel are required for the "_with_features_" version... for reasons...
typedef CGAL::Polyhedron_3< CGAL_KERNEL,
                            CGAL::Mesh_3::Mesh_polyhedron_items<int>,
                            CGAL::HalfedgeDS_default,
                            std::allocator<int> > Polyhedron;

/* Ugly "modifier" class required by Polyhedron.delegate() to actually generate it from a TSS3...
*/
template <class HDS>
class Build_Polyhedron3_From_geoTriSurface3_CGAL_MODIFIER : public CGAL::Modifier_base<HDS>
{
public:
    Build_Polyhedron3_From_geoTriSurface3_CGAL_MODIFIER( const geo::TriSurfaceShape3& tss3 )
    : m_rTSS3(tss3) {}
    void operator()( HDS& hds )
        {
            // Postcondition: hds is a valid polyhedral surface.
            CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, false );//true );
            typedef typename HDS::Vertex Vertex;
            typedef typename Vertex::Point Point;
            B.begin_surface( m_rTSS3.GetNumV(), m_rTSS3.GetNumT() );
            {
                // std::cout << "adding vertices..." << std::endl;
                for( uint32 it_v=0; it_v<m_rTSS3.GetNumV(); it_v++ )
                {
                    Vec3 v( m_rTSS3.V_Pos_0(it_v) );
                    B.add_vertex( Point(v.x(),v.y(),v.z()) );
                }
                // std::cout << "adding facets..." << std::endl;
                for( uint32 it_t=0; it_t<m_rTSS3.GetNumT(); it_t++ )
                {
                    B.begin_facet();
                    B.add_vertex_to_facet( m_rTSS3.T_VID(it_t,0) );
                    B.add_vertex_to_facet( m_rTSS3.T_VID(it_t,1) );
                    B.add_vertex_to_facet( m_rTSS3.T_VID(it_t,2) );
                    B.end_facet();
                }
            }
            // std::cout << "ending surface..." << std::endl;
            B.end_surface();
        }
private:
    const geo::TriSurfaceShape3& m_rTSS3;
};

bool Make_TriSurfaceShape3_From_TriSurfaceShape3_Detail( const geo::TriSurfaceShape3& src_surface,
                                                         float refinement_size,
                                                         geo::EditableTriSurfaceShape3& dst_surface )
{
    //\todo Refine or Simplify SRC to the required detail
    //\todo Consider CGAL http://doc.cgal.org/latest/Surface_mesh_simplification/
    dst_surface.Set( src_surface );
    return true;
}

/* 3D mesh simplification
   - See http://doc.cgal.org/latest/Surface_mesh_simplification/
*/
bool Make_TriSurfaceShape3_From_TriSurfaceShape3_Simplify( const geo::TriSurfaceShape3& src_surface,
                                                           float refinement_size,
                                                           geo::EditableTriSurfaceShape3& dst_surface )
{
    namespace SMS = CGAL::Surface_mesh_simplification;

    //\todo Create Polyhedron from geo::TriSurfaceShape3
    Polyhedron polyhedron;
    {
        Build_Polyhedron3_From_geoTriSurface3_CGAL_MODIFIER<typename Polyhedron::HalfedgeDS> modifier( src_surface );
        polyhedron.delegate( modifier );
#ifdef __ENABLE_CGAL_UTILS_TRACE
        std::cout << "Polyhedron created with #V = " << polyhedron.size_of_vertices()
                  << " #T = " << polyhedron.size_of_facets()
                  << " #HE = " << polyhedron.size_of_halfedges()
                  << " ?Closed = " << (polyhedron.is_closed()?"true":"false")
                  << std::endl;
#endif
    }

    // This is a stop predicate (defines when the algorithm terminates).
    // In this example, the simplification stops when the number of undirected edges
    // left in the surface mesh drops below the specified number (1000)
    SMS::Count_stop_predicate<Polyhedron> stop( refinement_size * polyhedron.size_of_halfedges()/2 );//1000 ); TEMPORAL HACK, use refinement_size as a fraction 0..1!!
    // This the actual call to the simplification algorithm.
    // The surface mesh and stop conditions are mandatory arguments.
    // The index maps are needed because the vertices and edges
    // of this surface mesh lack an "id()" field.
    int r = SMS::edge_collapse( polyhedron
                                ,stop
                                ,CGAL::vertex_index_map( get(CGAL::vertex_external_index,polyhedron) )
                                .halfedge_index_map( get(CGAL::halfedge_external_index,polyhedron) )
                                .get_cost( SMS::Edge_length_cost<Polyhedron>() )
                                .get_placement( SMS::Midpoint_placement<Polyhedron>() ) );
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "\nFinished...\n" << r << " edges removed and "
              << (polyhedron.size_of_halfedges()/2) << " final edges"
              << " is_pure_triangle = " << (polyhedron.is_pure_triangle() ? "true" : "false")
              << std::endl;
#endif
    // Abort if not pure triangle mesh
    if( !polyhedron.is_pure_triangle() )
        return false;

    // Convert polyhedron to tss3
    dst_surface.Clear();
    dst_surface.BeginEdition();
    {
        typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
        typedef typename Polyhedron::Vertex::Point Point_3;
        typedef typename Polyhedron::Vertex_const_iterator Vertex_const_iterator;
        typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;

        // Map Vertex_handle->vertex index
        std::map<Vertex_const_handle, int> map_vh_to_vid;
        uint32 it_v(0);
        // Vertices
        for( Vertex_const_iterator vit = polyhedron.vertices_begin(); vit != polyhedron.vertices_end(); ++vit )
        {
            map_vh_to_vid[vit] = it_v++;
            dst_surface.AddVertex( Vec3(vit->point().x(),vit->point().y(),vit->point().z()) );
        }
        // Faces
        for( Facet_const_iterator fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit )
        {
            Polyhedron::Halfedge_const_handle he( fit->halfedge() );
            Vertex_const_handle vh0( he->vertex() ); he = he->next();
            Vertex_const_handle vh1( he->vertex() ); he = he->next();
            Vertex_const_handle vh2( he->vertex() ); he = he->next();
            dst_surface.AddTriangle( map_vh_to_vid[vh0], map_vh_to_vid[vh1], map_vh_to_vid[vh2] );
        }
    }
    dst_surface.EndEdition();
    return true;
}

bool Make_TriSurfaceShape3_From_TriSurfaceShape3_Offset( const geo::TriSurfaceShape3& src_surface,
                                                         float offset,
                                                         geo::EditableTriSurfaceShape3& dst_surface )
{
    // Compute vertex normals
    Vec3* vec_vertex_normals = new Vec3[src_surface.GetNumV()];
    memset( vec_vertex_normals, 0, sizeof(Vec3)*src_surface.GetNumV() );
    for( uint32 it_tri=0; it_tri<src_surface.GetNumT(); it_tri++ )
    {
        uint32 vid0 = src_surface.T_VID(it_tri,0);
        uint32 vid1 = src_surface.T_VID(it_tri,1);
        uint32 vid2 = src_surface.T_VID(it_tri,2);
        Vec3f weighted_normal = mal::Cross( src_surface.V_Pos_0(vid1) - src_surface.V_Pos_0(vid0),
                                            src_surface.V_Pos_0(vid2) - src_surface.V_Pos_0(vid0) );
        vec_vertex_normals[vid0] += weighted_normal;
        vec_vertex_normals[vid1] += weighted_normal;
        vec_vertex_normals[vid2] += weighted_normal;
    }
    for( uint32 it_v=0; it_v<src_surface.GetNumV(); it_v++ )
        if( mal::NormSq( vec_vertex_normals[it_v] ) > mal::Epsilon<Real>() )
            vec_vertex_normals[it_v] = mal::Normalized( vec_vertex_normals[it_v] );
    // Create offset surface
    dst_surface.Clear();
    dst_surface.BeginEdition();
    {
        for( uint32 it_v=0; it_v<src_surface.GetNumV(); it_v++ )
            dst_surface.AddVertex( src_surface.V_Pos_0(it_v) + offset*vec_vertex_normals[it_v] );
        for( uint32 it_tri=0; it_tri<src_surface.GetNumT(); it_tri++ )
            dst_surface.AddTriangle( src_surface.T_VID(it_tri,0),
                                     src_surface.T_VID(it_tri,1),
                                     src_surface.T_VID(it_tri,2) );
    }
    dst_surface.EndEdition();
    delete[] vec_vertex_normals;
    return true;
}


/* 3D Mesh generation
  1) Convert geo::TriSurfaceShape3 to Polyhedral_mesh_domain_3
  2) Use CGAL 3D mesh generation from a http://doc.cgal.org/latest/Mesh_3/index.html to generate a Mesh_3
     - See Mesh_3/mesh_polyhedral_domain.cpp CGAL example
  3) Convert Mesh_3 to geo::TetSolidShape3
*/
bool Make_TetSolidShape3_From_TriSurfaceShape3_CDT( const geo::TriSurfaceShape3& surface,
                                                    float cdt_facet_angle, float cdt_facet_size, float cdt_facet_distance,
                                                    float cdt_cell_ratio, float cdt_cell_size,
                                                    bool cdt_lloyd, bool cdt_odt, bool cdt_perturb, bool cdt_exude,
                                                    geo::EditableTetSolidShape3& solid )
{
    //typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, CGAL_KERNEL> Mesh_domain;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<CGAL_KERNEL> Mesh_domain;
    // Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
    typedef CGAL::Mesh_triangulation_3<
        Mesh_domain,
        CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
        CGAL::Parallel_tag                        // Tag to activate parallelism
        >::type Tr;
#else
    typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif
    //typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    //typedef CGAL::Mesh_criteria_with_features_3<Tr> Mesh_criteria;

    //\todo Create Polyhedron from geo::TriSurfaceShape3
    Polyhedron polyhedron;
    {
        Build_Polyhedron3_From_geoTriSurface3_CGAL_MODIFIER<typename Polyhedron::HalfedgeDS> modifier( surface );
        polyhedron.delegate( modifier );
#ifdef __ENABLE_CGAL_UTILS_TRACE
        std::cout << "Polyhedron created with #V = " << polyhedron.size_of_vertices()
                  << " #T = " << polyhedron.size_of_facets()
                  << " #HE = " << polyhedron.size_of_halfedges()
                  << " ?Closed = " << (polyhedron.is_closed()?"true":"false")
                  << std::endl;
#endif
    }

    // Create domain
    Mesh_domain domain(polyhedron);

    // TEMP
    domain.detect_features(120); //max dihedral angle to consider edge as feature (180 => ALL) \todo publish param
    //\todo We could add_features for all undirected edges in the Polyhedron... this would be the same as the "constraints" added in 2D mesh gen, I guess...

    // Meshing criteria \todo _with_features?!?!?!
    Mesh_criteria criteria( CGAL::parameters::edge_size = cdt_facet_size, //0.025 //\todo
                            CGAL::parameters::facet_angle = cdt_facet_angle, //15,
                            CGAL::parameters::facet_size = cdt_facet_size, //0.15,
                            CGAL::parameters::facet_distance = cdt_facet_distance, //0.008,
                            CGAL::parameters::cell_size = cdt_cell_size, //0.1,
                            CGAL::parameters::cell_radius_edge_ratio = cdt_cell_ratio ); //3 );

    // Mesh generation
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Make_TetSolidShape3_From_TriSurfaceShape3_CDT() Making Mesh_complex_3_in_triangulation_3..." << std::endl;
#endif
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>( domain, criteria,
                                         CGAL::parameters::features(),
                                         cdt_lloyd ? CGAL::parameters::lloyd() : CGAL::parameters::no_lloyd(),
                                         cdt_odt ? CGAL::parameters::odt() : CGAL::parameters::no_odt(),
                                         cdt_perturb ? CGAL::parameters::perturb() : CGAL::parameters::no_perturb(),
                                         cdt_exude ? CGAL::parameters::exude() : CGAL::parameters::no_exude() );

#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "...Made Mesh_complex_3_in_triangulation_3 with #V = " << c3t3.triangulation().number_of_vertices() //\note Intuitively it would be c3t3.number_of_vertices_in_complex() but it's always 0 !?
              << " #F = " << c3t3.number_of_facets_in_complex()
              << " #T = " << c3t3.number_of_cells_in_complex()
              << std::endl;
#endif

    /* Mesh refinement TEMPORALLY disabled, it from the doc it seems
       that make_mesh_3 ALREADY does everything with the SAME
       available criteria, this is just an incremental optimization.

    // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets //\todo receive params!
    Mesh_criteria new_criteria(CGAL::parameters::cell_radius_edge_ratio=3, CGAL::parameters::cell_size=0.1);//TEMP 0.03);
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Make_TetSolidShape3_From_TriSurfaceShape3_CDT() Refining Mesh3..." << std::endl;
#endif
    CGAL::refine_mesh_3(c3t3, domain, new_criteria);
    */

    // Create geo::TetSolidShape3 from C3t3
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Make_TetSolidShape3_From_TriSurfaceShape3_CDT() Creating TetSolidShape3... " << std::endl;
#endif
    solid.Clear();
    solid.BeginEdition();
    {
        //\note This is HELL, Mesh_3 access taken from C3T3 output_to_medit() function in File_medit.h
        typedef typename C3t3::Triangulation Tr;
        typedef typename C3t3::Facets_in_complex_iterator Facet_iterator;
        typedef typename C3t3::Cells_in_complex_iterator Cell_iterator;
        typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef typename Tr::Vertex_handle Vertex_handle;
        typedef typename Tr::Point Point_3;
        const Tr& tr = c3t3.triangulation();
        // Map Vertex_handle->vertex index
        std::map<Vertex_handle, int> map_vh_to_vid;
        uint32 it_v(0);
        // Vertices
        for( Finite_vertices_iterator vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit )
        {
            map_vh_to_vid[vit] = it_v++;
            solid.AddVertex( Vec3(vit->point().x(),vit->point().y(),vit->point().z()) );
        }
        // Tetrahedra
        for( Cell_iterator cit = c3t3.cells_in_complex_begin(); cit != c3t3.cells_in_complex_end(); ++cit )
            solid.AddTetrahedron( map_vh_to_vid[cit->vertex(0)], map_vh_to_vid[cit->vertex(1)], map_vh_to_vid[cit->vertex(2)], map_vh_to_vid[cit->vertex(3)] );
        // Boundary Faces
        for( Facet_iterator fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit )
        {
            uint32 vec_bf_vid[3];
            Vec3 vec_bf_pos[3];
            uint32 num_bv_vid(0);
            /*This code is insane... it seems that the fit stores as
              "second" the index of a vertex that is NOT on the
              facet... probably the out-of-plane vtx of the
              tetrahedron that generated this facet... whatever */
            for( int i=0; i<4; i++ )
                if( i != fit->second )
                {
                    vec_bf_vid[ num_bv_vid ] = map_vh_to_vid[ (*fit).first->vertex(i) ];
                    vec_bf_pos[ num_bv_vid ] = Vec3( (*fit).first->vertex(i)->point().x(), (*fit).first->vertex(i)->point().y(), (*fit).first->vertex(i)->point().z() );
                    num_bv_vid++;
                }

            /*TEMP: The following block is deprecated because BF are
              now automatically computed in EndEdition() from Tet
              topology, but we'll keep it to remember CGAL's strange
              behaviour, just in case...
            */
#ifdef __CGAL_FIX_FACET_ORIENTATION
            /*\todo C3t3 does NOT seem to guarantee CCW facet
              orientation, therefore we must enforce... however, the
              following code DOES NOT WORK... maybe fit->second is NOT
              what I think it is... or maybe it cannot be retrieved
              with first->vertex( second )... all in all I'll compute
              BF *automatically* using topology in geo::TetSolidShape3
            */
            Vec3 interior_p( (*fit).first->vertex( fit->second )->point().x(), (*fit).first->vertex( fit->second )->point().y(), (*fit).first->vertex( fit->second )->point().z() );
            if( mal::Dot( interior_p - vec_bf_pos[0], mal::Cross(vec_bf_pos[1]-vec_bf_pos[0],vec_bf_pos[2]-vec_bf_pos[0]) ) <= 0 )
                solid.AddBoundaryFace( vec_bf_vid[0], vec_bf_vid[1], vec_bf_vid[2] );
            else
                solid.AddBoundaryFace( vec_bf_vid[0], vec_bf_vid[2], vec_bf_vid[1] );
#else
            // Just add it and EndEdition() will correct them...
            //solid.AddBoundaryFace( vec_bf_vid[0], vec_bf_vid[1], vec_bf_vid[2] );
#endif
        }
    }
    solid.EndEdition();
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Make_TetSolidShape3_From_TriSurfaceShape3_CDT() Created TetSolidShape3 with #V = " << solid.GetNumV()
              << " #T = " << solid.GetNumT()
              << " #BF = " << solid.GetNumBF()
              << std::endl;
#endif

    return true;
}

} //namespace geo
