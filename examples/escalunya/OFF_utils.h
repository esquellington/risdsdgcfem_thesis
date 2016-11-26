#ifndef ESCALUNYA_OFF_UTILS_H
#define ESCALUNYA_OFF_UTILS_H

#include "Config.h"
#include <string>
#include <fstream>
#include <Geo/shape/shape.h>

/*! Load an OFF mesh */
bool OFF_Load_ETSS3( geo::EditableTriSurfaceShape3& etss3, const std::string& file_name, bool b_invert_normals = false )
{
    std::ifstream ifs( file_name );
    std::string magic;
    ifs >> magic;
    if( magic != "OFF" ) return false;

    uint32 num_vertices, num_triangles, num_edges;
    ifs >> num_vertices;
    ifs >> num_triangles;
    ifs >> num_edges;

    etss3.BeginEdition();
    {
        for( uint32 it_v=0; it_v<num_vertices; it_v++ )
        {
            geo::Vec3 v;
            ifs >> v[0] >> v[1] >> v[2];
            etss3.AddVertex( v );
        }
        for( uint32 it_t=0; it_t<num_triangles; it_t++ )
        {
            uint32 num_vid;
            ifs >> num_vid;
            if( num_vid != 3 ) { etss3.EndEdition(); return false; }
            uint32 vec_vid[3];
            ifs >> vec_vid[0] >> vec_vid[1] >> vec_vid[2];
            if( b_invert_normals )
                etss3.AddTriangle( vec_vid[0], vec_vid[2], vec_vid[1] );
            else
                etss3.AddTriangle( vec_vid[0], vec_vid[1], vec_vid[2] ); //\todo OFF format does NOT guarantee CW (bunny) or CCW (elephant, squirrel...)
        }
    }
    etss3.EndEdition();
    return true;
}

#endif //ESCALUNYA_OFF_UTILS_H
