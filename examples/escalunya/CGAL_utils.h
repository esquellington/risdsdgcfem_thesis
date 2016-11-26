#ifndef ESCALUNYA_CGAL_UTILS_H
#define ESCALUNYA_CGAL_UTILS_H

#include "Config.h"
#include <Geo/shape/shape.h>

namespace geo
{

bool Make_TriSurfaceShape3_From_TriSurfaceShape3_Detail( const geo::TriSurfaceShape3& src_surface,
                                                         float refinement_size,
                                                         geo::EditableTriSurfaceShape3& dst_surface );

bool Make_TriSurfaceShape3_From_TriSurfaceShape3_Simplify( const geo::TriSurfaceShape3& src_surface,
                                                           float refinement_size,
                                                           geo::EditableTriSurfaceShape3& dst_surface );

bool Make_TriSurfaceShape3_From_TriSurfaceShape3_Offset( const geo::TriSurfaceShape3& src_surface,
                                                         float offset,
                                                         geo::EditableTriSurfaceShape3& dst_surface );

bool Make_TetSolidShape3_From_TriSurfaceShape3_CDT( const geo::TriSurfaceShape3& surface,
                                                    float cdt_facet_angle, float cdt_facet_size, float cdt_facet_distance,
                                                    float cdt_cell_ratio, float cdt_cell_size,
                                                    bool cdt_lloyd, bool cdt_odt, bool cdt_perturb, bool cdt_exude,
                                                    geo::EditableTetSolidShape3& solid );

}

#endif //ESCALUNYA_CGAL_UTILS_H
