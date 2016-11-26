#ifndef ESCALUNYA_TET_UTILS_H
#define ESCALUNYA_TET_UTILS_H

#include "Config.h"
#include <Geo/shape/shape.h>

namespace geo
{

/* Tetrahedralize an AABB using minimal hexaedral grid tetrahedralization. */
bool Make_TetSolidShape3_From_AABB3_Minimal( const geo::Vec3& pos_min, const geo::Vec3& pos_max,
                                             geo::Real cell_size,
                                             geo::EditableTetSolidShape3& solid );

/* Tetrahedralize an AABB using BCC lattice
   \see "Isosurface Stuffing: Fast Tetrahedral Meshes with Good Dihedral Angles"
*/
bool Make_TetSolidShape3_From_AABB3_BCC( const geo::Vec3& pos_min, const geo::Vec3& pos_max,
                                         geo::Real cell_size,
                                         //\todo Other params...?
                                         geo::EditableTetSolidShape3& solid );

/* Fit a given TetSolid to a contained TriSurface:
   - Remove external tet
   - Adjust boundary tet
   - Adjust internal tet
*/
bool Fit_TetSolidShape3_To_TriSurfaceShape3( const geo::TriSurfaceShape3& surface,
                                             float odt_relaxation_coeff, float lpc_relaxation_coeff, unsigned int max_iter, bool b_fix_nmf,
                                             geo::EditableTetSolidShape3& solid );

/* Clip a TriSurfaceShape3 to a TetSolidShape3 */
bool Clip_TriSurfaceShape3_TetSolidShape3( const geo::TriSurfaceShape3& surface,
                                           const geo::EditableTetSolidShape3& solid, const BVH_TetSolidShape3* pSolidBVH, //TEMP pass BVH as param by now
                                           geo::EditableTriSurfaceShape3& clipped_surface );

}

#endif //ESCALUNYA_TET_UTILS_H
