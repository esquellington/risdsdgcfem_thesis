#include "CGAL_utils.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/create_offset_polygons_2.h>

//#define __ENABLE_CGAL_UTILS_TRACE
#ifdef __ENABLE_CGAL_UTILS_TRACE
#  include <iostream>
#endif

#define __ENABLE_BEZIER_TARGET_EDGE_LENGTH

namespace geo
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGAL_KERNEL;

void Make_PolygonalShape2_From_PathShape2( const PathShape2& path,
                                           float refinement_size, bool b_subdivide_line_ct,
                                           EditablePolygonalShape2& polygonal )
{
    polygonal.Clear();
    if( !polygonal.IsAlloc() ) polygonal.Alloc( 10000, geo::eVDL_Points, geo::eVDL_Points ); //TEMPORAL
    polygonal.BeginEdition();
    {
        // Init v0
        geo::Vec2 p0 = path.V_Pos_0(0);
        /*IMPORTANT: We do NOT add the 1st vtx to avoid
          repetition, CGAL::Polygon_2 does NOT duplicate
          first-last vtx.
          polygon.push_back( K::Point_2( p0.x(), p0.y() ) );
        */
        for( unsigned int it_edge=0; it_edge<path.GetNumE(); it_edge++ )
        {
            // Add v1 stuff
            geo::Vec2 p1 = path.V_Pos_0( path.E_FinalVID(it_edge) );
            // Add vertices and constraints between v0 and v1
            switch( path.E_CurveType(it_edge) )
            {
            case geo::PathShape2::edge_type::eCT_Line:
                if( b_subdivide_line_ct )
                {
                    float total_arc_length( mal::Norm(p1-p0) );
                    // Compute num segments required
                    unsigned int num_segments = mal::Ceil( total_arc_length / refinement_size );
                    APP_ASSERT( num_segments > 0 );
                    // Discretize using arc-length parametrization
                    geo::Vec2 b0( p0 );
                    // Adaptative refinement internal segments
                    for( unsigned int i=1; i<num_segments; i++ )
                    {
                        // Compute next vertex's arclength
                        float lambda( float(i) / num_segments );
                        // Lerp and add
                        geo::Vec2 b1 = (1-lambda)*p0 + lambda* p1;
                        polygonal.AddPoint( b1 );
                        b0 = b1;
                    }
                }
                polygonal.AddPoint( p1 );
                break;
            case geo::PathShape2::edge_type::eCT_Bezier3:
                {
                    // Bezier curve
                    geo::Vec2 bezier3_a( path.E_Data(it_edge).m_ParamA ); //\todo If params are also SDOF, use actual_sdof here!!
                    geo::Vec2 bezier3_b( path.E_Data(it_edge).m_ParamB );
#  ifdef __ENABLE_BEZIER_TARGET_EDGE_LENGTH
                    //-- Sampled curve with target edge length
                    // Compute curve length
                    float total_arc_length(0);
                    geo::Vec2 b0( p0 );
                    const unsigned int cNumArcLengthSamples(1000); // With so many samples there's no need to Lerp() during refinement (*)
                    float vec_arc_length[cNumArcLengthSamples];
                    for( unsigned int i=0; i<cNumArcLengthSamples; i++ )
                    {
                        vec_arc_length[i] = total_arc_length;
                        float lambda( float(i) / (cNumArcLengthSamples-1) );
                        geo::Vec2 b1 = geo::Eval_Bezier3( p0, p1, bezier3_a, bezier3_b, lambda );
                        total_arc_length += mal::Norm( b0 - b1 );
                        b0 = b1;
                    }
                    // Compute num segments required
                    unsigned int num_segments = mal::Ceil( total_arc_length / refinement_size );
                    APP_ASSERT( num_segments > 0 );
                    // Discretize using arc-length parametrization
                    b0 = p0;
                    // Adaptative refinement internal segments
                    unsigned int lambda_index(0);
                    for( unsigned int i=1; i<num_segments; i++ )
                    {
                        // Compute next vertex's arclength
                        float arc_length( float(i) * (total_arc_length / num_segments) );
                        // Find lambda at desired arclength
                        while( vec_arc_length[lambda_index+1] < arc_length ) lambda_index++;
                        // (*) Could Lerp(l,l+1) with proper weights, but with sufficient cNumArcLengthSamples it's unnecessary
                        float lambda( float(lambda_index) / (cNumArcLengthSamples-1) );
                        // Eval bezier and add
                        geo::Vec2 b1 = geo::Eval_Bezier3( p0, p1, bezier3_a, bezier3_b, lambda );
                        polygonal.AddPoint( b1 );
                        b0 = b1;
                    }
#  else //__ENABLE_BEZIER_TARGET_EDGE_LENGTH
                    //-- Sampled curve with FIXED segments
                    geo::Vec2 b0( p0 );
                    // Add edge refinement internal segments
                    const unsigned int cNumSegments(10);
                    for( unsigned int i=1; i<cNumSegments; i++ )
                    {
                        geo::Vec2 b1 = geo::Eval_Bezier3( p0, p1, bezier3_a, bezier3_b, geo::Real(i) / cNumSegments );
                        polygonal.AddPoint( b1 );
                        b0 = b1;
                        bh0 = bh1;
                    }
#  endif //__ENABLE_BEZIER_TARGET_EDGE_LENGTH
                    polygonal.AddPoint( p1 );
                }
                break;
            default:
                APP_LOG_ERROR("RebuildCDT: Unknown geo::PathShape2::ECurveType %d!!", (int)path.E_CurveType(it_edge) );
                break;
            }
            // Advance
            p0 = p1;
        }
    }
    polygonal.SetClosed( true );
    polygonal.EndEdition();
}

bool Make_PolygonShape2_From_PolygonShape2_Offset( const geo::PolygonalShape2& src_polygonal,
                                                   float offset,
                                                   geo::EditablePolygonalShape2& dst_polygonal )
{
    // Fill cgal_polygon from src_polygonal
    CGAL::Polygon_2<CGAL_KERNEL> cgal_polygon;
    for( unsigned int it_v=0; it_v<src_polygonal.GetNumVertices(); it_v++ )
    {
        geo::Vec2 p( src_polygonal.V_Pos_0( it_v ) );
        cgal_polygon.push_back( CGAL_KERNEL::Point_2( p.x(), p.y() ) );
    }

    // Compute polygon offset
    std::vector< boost::shared_ptr< CGAL::Polygon_2<CGAL_KERNEL> > > vec_offset_polygon = CGAL::create_exterior_skeleton_and_offset_polygons_2( mal::Max<double>(1e-6,offset), cgal_polygon );
    /*IMPORTANT: For SOME reason, the interesting polygon is
      always [1], while [0] is something like an AABB... see
      "Exterior Skeletons and Exterior Offset Contours" in
      \see http://doc.cgal.org/latest/Straight_skeleton_2/index.html#Chapter_2D_Straight_Skeleton_and_Polygon_Offsetting
    */
    APP_ASSERT( vec_offset_polygon.size() > 1 );
    const CGAL::Polygon_2<CGAL_KERNEL>& offset_polygon( *vec_offset_polygon[1] );

    // Fill dst_polygonal from cgal_polygon
    dst_polygonal.Clear();
    dst_polygonal.BeginEdition();
    for( CGAL::Polygon_2<CGAL_KERNEL>::Vertex_const_iterator it_vtx = offset_polygon.vertices_begin();
         it_vtx != offset_polygon.vertices_end();
         ++it_vtx )
        dst_polygonal.AddPoint( geo::Vec2( it_vtx->x(), it_vtx->y() ) );
    dst_polygonal.SetClosed( true );
    dst_polygonal.EndEdition();

    return true;
}

/* Constrained Delaunay Triangulation
  1) *Build a CGAL::Polygon_2 Poly0 from source geo::PathShape2 EPS2, conveniently refined \see http://doc.cgal.org/latest/Polygon/classCGAL_1_1Polygon__2.html
  2) *Compute CGAL::Polygon_2 Poly1 from CGAL::create_exterior_skeleton_and_offset_polygons_2( offset, Poly0 ) \see http://doc.cgal.org/latest/Straight_skeleton_2/index.html#Chapter_2D_Straight_Skeleton_and_Polygon_Offsetting
  3) *Compute CGAL::Constrained_Delaunay_triangulation_2 CDT from Poly1 \see http://doc.cgal.org/latest/Mesh_2/index.html#Chapter_2D_Conforming_Triangulations_and_Meshes
  4) *Build geo::MeshSolidShape2 from CDT
*/
bool Make_MeshSolidShape2_From_PolygonalShape2_CDT( const geo::PolygonalShape2& polygonal,
                                                    float cdt_ratio, float cdt_size,
                                                    geo::EditableMeshSolidShape2& mesh )
{
    typedef CGAL::Triangulation_vertex_base_2<CGAL_KERNEL> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<CGAL_KERNEL> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<CGAL_KERNEL,Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

    // CDT: insert vertices/constraints
    CDT cdt;
    geo::Vec2 p = polygonal.V_Pos_0( 0 );
    CDT::Vertex_handle vh_first = cdt.insert( CDT::Point( p.x(), p.y() ) );
    CDT::Vertex_handle vh0 = vh_first;
    unsigned int num_constraints(0);
    for( unsigned int it_v=1; it_v<polygonal.GetNumVertices(); it_v++ )
    {
        p = polygonal.V_Pos_0( it_v );
        CDT::Vertex_handle vh1 = cdt.insert( CDT::Point( p.x(), p.y() ) );
        cdt.insert_constraint( vh0, vh1 );
        vh0 = vh1;
        num_constraints++;
    }
    cdt.insert_constraint( vh0, vh_first );
    num_constraints++;
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Inserted constraints = " << num_constraints << std::endl;
#endif
    // Compute CDT from polygon
    std::list<CDT::Point> list_of_seeds;
    list_of_seeds.push_back(CDT::Point(1000, 1000)); //\todo OUTSIDE point seed to mark the exterior of the mesh
    CGAL::refine_Delaunay_mesh_2( cdt, list_of_seeds.begin(), list_of_seeds.end(),
                                  Criteria( cdt_ratio, cdt_size ) );
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Number of vertices before: "
              << cdt.number_of_vertices() << std::endl;
#endif
    // make it conforming Delaunay
    CGAL::make_conforming_Delaunay_2(cdt);
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Number of vertices after make_conforming_Delaunay_2: "
              << cdt.number_of_vertices() << std::endl;
#endif
    // then make it conforming Gabriel
    CGAL::make_conforming_Gabriel_2(cdt);
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Number of vertices after make_conforming_Gabriel_2: "
              << cdt.number_of_vertices() << std::endl;
    // count faces
    std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
#endif
    int mesh_faces_counter = 0;
    for( CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
         fit != cdt.finite_faces_end();
         ++fit )
        if( fit->is_in_domain() )
            ++mesh_faces_counter;
#ifdef __ENABLE_CGAL_UTILS_TRACE
    std::cout << "Number of faces in the mesh domain: " << mesh_faces_counter << std::endl;
#endif
    APP_ASSERT( mesh_faces_counter > 0 );

    // Build MeshSolidShape2 from CDT
    mesh.Clear();
    mesh.BeginEdition();
    {
        // Add vertices while building map from CDT::Vertex_handle to vertex index
        std::map< CDT::Vertex_handle, geo::feature_index_type > map_vh_to_vid;
        for( CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit )
        {
            geo::feature_index_type vid = mesh.AddVertex( geo::Vec2( vit->point().x(), vit->point().y() ) );
            map_vh_to_vid[ vit ] = vid;
        }
        // Add faces using map_vh_to_vid
        for( CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit )
            if( fit->is_in_domain() )
            {
                CDT::Vertex_handle vh0( fit->vertex(0) );
                CDT::Vertex_handle vh1( fit->vertex(1) );
                CDT::Vertex_handle vh2( fit->vertex(2) );
                mesh.AddPolygon3( map_vh_to_vid[vh0], map_vh_to_vid[vh1], map_vh_to_vid[vh2] );
            }
    }
    mesh.EndEdition();
    return true;
}

} //namespace geo
