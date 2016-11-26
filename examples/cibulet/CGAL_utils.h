#ifndef CIBULET_CGAL_UTILS_H
#define CIBULET_CGAL_UTILS_H

#include "Config.h"
#include <Geo/shape/shape.h>

namespace geo
{

/*\todo This function does not actually use CGAL, consider moving elsewhere */
void Make_PolygonalShape2_From_PathShape2( const geo::PathShape2& ps2,
                                           float refinement_size, bool b_subdivide_line_ct,
                                           geo::EditablePolygonalShape2& eps2 );

bool Make_PolygonShape2_From_PolygonShape2_Offset( const geo::PolygonalShape2& src_polygonal,
                                                   float offset,
                                                   geo::EditablePolygonalShape2& dst_polygonal );

bool Make_MeshSolidShape2_From_PolygonalShape2_CDT( const geo::PolygonalShape2& polygonal,
                                                    float cdt_ratio, float cdt_size,
                                                    geo::EditableMeshSolidShape2& mesh );

}

#endif //CIBULET_CGAL_UTILS_H
