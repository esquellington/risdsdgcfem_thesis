#include <iostream>
#include <fstream>
#include <string>

#include <Geo/geo.h>
#include <Geo/IObject.h>
#include <Geo/shape/shape.h>
#include <Geo/util/Viz.h>
#include <Geo/shape/ShapeLibrary.h>
#include <Geo/shape/ShapeFactory.h>
#include <Geo/bv/BoundingVolume.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define APP_ASSERT( x )
#  define APP_LOG( x, ... )
#  define APP_LOG_WARNING( x, ... )
#  define APP_LOG_ERROR( x, ... )
#  define APP_LOG_ASSERT( x, y, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define APP_ASSERT( x ) assert(x)
#  define APP_LOG( x, ... ) { printf("<APP LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_WARNING( x, ... ) { printf("<APP WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_ERROR( x, ... ) { printf("<APP ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define APP_LOG_ASSERT( x, y, ... ) { if(!(x)) { printf("<APP ASSERT> "); printf( y, ##__VA_ARGS__ ); printf("\n"); assert(x); } }
#endif

int main( int argc, const char *argv[] )
{
    geo::Init();

    // Parse params
    if( argc > 4
        && std::string(argv[1]) == "-f"   //input tetgen file(s)
        && std::string(argv[3]) == "-o" ) //output
    {
        //Open
        std::string tetgen_file_base_name( argv[2] );
        std::string node_file( tetgen_file_base_name + ".node" );
        std::string ele_file( tetgen_file_base_name + ".ele" );
        std::string face_file( tetgen_file_base_name + ".face" );

        //force normalized scale
        float scale_factor(1);
        if( argc > 6
            && std::string(argv[5]) == "-s" )
        {
            float natural_scale = atof( argv[6] );
            // Compute AABB
            geo::bv::AABB3 aabb;
            std::ifstream fs;
            int base_vid(0); //Base vertex id IS NOT NECESSARILY 1
            int num_nodes(0);
            int num_tet(0);
            int num_tri(0);
            //-- AddVertex
            {
                int dimension, num_attrib, num_boundary;
                fs.open( node_file.c_str() );
                fs >> num_nodes;
                fs >> dimension;
                APP_ASSERT( dimension == 3 );
                fs >> num_attrib;
                fs >> num_boundary;
                for( int it_node=0; it_node<num_nodes; it_node++ )
                {
                    int vid;
                    fs >> vid;
                    if( 0 == it_node ) base_vid = vid;
                    Vec3f pos;
                    fs >> pos[0]; fs >> pos[1]; fs >> pos[2];
                    if( 0 == it_node ) aabb.SetMinMax( pos, pos );
                    else aabb.Merge( pos );
                }
                fs.close();
            }
        }

        geo::EditableTetSolidShape3 etss3;
        etss3.BeginEdition();
        {
            std::ifstream fs;
            int base_vid(0); //Base vertex id IS NOT NECESSARILY 1
            int num_nodes(0);
            int num_tet(0);
            int num_tri(0);
            geo::bv::AABB3 aabb;
            //-- AddVertex
            {
                int dimension, num_attrib, num_boundary;
                fs.open( node_file.c_str() );
                fs >> num_nodes;
                fs >> dimension;
                APP_ASSERT( dimension == 3 );
                fs >> num_attrib;
                fs >> num_boundary;
                for( int it_node=0; it_node<num_nodes; it_node++ )
                {
                    int vid;
                    fs >> vid;
                    if( 0 == it_node ) base_vid = vid;
                    Vec3f pos;
                    fs >> pos[0]; fs >> pos[1]; fs >> pos[2];
                    // Skip attrib and boundary markers
                    for( int it_attrib=0; it_attrib<num_attrib; it_attrib++ ) { float dummy; fs >> dummy; }
                    for( int it_bm=0; it_bm<num_boundary; it_bm++ ) { float dummy; fs >> dummy; }
                    //APP_LOG_WARNING("Adding vtx[%d] (%f,%f,%f)",vid,pos[0],pos[1],pos[2]);
                    etss3.AddVertex( pos );
                    // Update AABB
                    if( 0 == it_node ) aabb.SetMinMax( pos, pos );
                    else aabb.Merge( pos );
                }
                fs.close();
            }
            //-- AddTetrahedron
            {
                int num_nodes_per_tet, num_attrib;
                fs.open( ele_file.c_str() );
                fs >> num_tet;
                fs >> num_nodes_per_tet;
                APP_ASSERT( num_nodes_per_tet == 4 );
                fs >> num_attrib;
                for( int it_tet=0; it_tet<num_tet; it_tet++ )
                {
                    int tid;
                    fs >> tid;
                    // read base-1 indices
                    int vid0,vid1,vid2,vid3;
                    fs >> vid0; fs >> vid1; fs >> vid2; fs >> vid3;
                    // convert to base-0 indices
                    vid0 -= base_vid;
                    vid1 -= base_vid;
                    vid2 -= base_vid;
                    vid3 -= base_vid;
                    // Skip attrib
                    for( int it_attrib=0; it_attrib<num_attrib; it_attrib++ ) { float dummy; fs >> dummy; }
                    //APP_LOG_WARNING("Adding tet[%d] (%d,%d,%d,%d)",tid,vid0+1,vid1+1,vid2+1,vid3+1);
                    etss3.AddTetrahedron( vid0, vid1, vid2, vid3 );
                }
                fs.close();
            }
            //-- AddBoundaryFace
            {
                int num_boundary;
                fs.open( face_file.c_str() );
                fs >> num_tri;
                fs >> num_boundary;
                for( int it_tri=0; it_tri<num_tri; it_tri++ )
                {
                    int bfid;
                    fs >> bfid;
                    // read base-1 indices
                    int vid0,vid1,vid2;
                    fs >> vid0; fs >> vid1; fs >> vid2;
                    // convert to base-0 indices
                    vid0 -= base_vid;
                    vid1 -= base_vid;
                    vid2 -= base_vid;
                    // Skip boundary markers
                    for( int it_bm=0; it_bm<num_boundary; it_bm++ ) { float dummy; fs >> dummy; }
                    //APP_LOG_WARNING("Adding bf[%d] (%d,%d,%d)",bfid,vid0+1,vid1+1,vid2+1);
                    //DEPRECATED! etss3.AddBoundaryFace( vid0, vid2, vid1 ); //2-1 swapped because TetGen Boundary Faces are CW, and Geo expects CCW
                }
                fs.close();
            }
            APP_LOG("Converted %s with %d nodes, %d tets, %d tris (base_vid = %d)",
                    tetgen_file_base_name.c_str(), num_nodes, num_tet, num_tri, base_vid );

            //-- Force natural scale and center, optionally
            if( argc > 6
                && std::string(argv[5]) == "-nsac" )
            {
                float natural_scale = atof( argv[6] );
                float scale_factor = natural_scale / mal::Max(2*aabb.GetHalfSizes());
                // Center
                etss3.Transform( Transform3f( -aabb.GetPos(),
                                              Mat3x3f::Identity() ) );
                // Scale at center
                etss3.Transform( Transform3f( Vec3f::Zero(),
                                              mal::GMatNxN_From_Diagonal( Vec3f(scale_factor,scale_factor,scale_factor) ) ) );
                APP_LOG("Applied natural scaling factor %f and centering (%f,%f,%f) to Original AABB Center (%f,%f,%f) HalfSizes (%f,%f,%f)",
                        scale_factor, -aabb.GetPos()[0], -aabb.GetPos()[1], -aabb.GetPos()[2],
                        aabb.GetPos()[0], aabb.GetPos()[1], aabb.GetPos()[2],
                        aabb.GetHalfSizes()[0], aabb.GetHalfSizes()[1], aabb.GetHalfSizes()[2] );
            }

        }
        etss3.EndEdition();

        // Save
        {
            std::string output_file_name( argv[4] );
            // Create SL
            geo::ShapeLibrary sl;
            sl.Reserve( 1<<24 ); //2^24 = 16mb
            // Register TSS
            sl.Register( etss3 );
            // Save SL
            sl.Save( output_file_name.c_str(), false );
        }
    }
    else
    {
        APP_LOG_WARNING( "Usage: teggen_to_geo -f tetgen_file_base_name -o geo_file_name" );
    }
    geo::ShutDown();
    return 0;
}
