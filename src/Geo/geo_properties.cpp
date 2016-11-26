#include "geo_properties.h"
#include <Geo/util/Viz.h>

namespace geo {

// Add a DDF property to an archetype being edited in an AL
void ArchetypeLibrary_AddProperty_DDF( util::ArchetypeLibrary& al, const char* name,
                                       Flags32 value, machine_uint_type offset,
                                       util::IArchetypeInstance::notify_touched_property_function_type ntpf )
{
    const char* vec_names[] = { // Object
                                "O.Frm",
                                "O.Shp",
                                "O.Emb",
                                "O.BV",
                                //Shape
                                "S.Bnd",
                                "S.Int",
                                "S.Tpl",
                                "S.V",
                                "S.E",
                                "S.F",
                                "S.Vol",
                                "S.Ctrl",
                                "S.Norm",
                                "S.Tan",
                                "S.Id",
                                "S.DCR",
                                "S.BVH",
                                // BVH
                                "BVH.Geo",
                                "BVH.Tpl",
                                "BVH.Lev" };
    int32 vec_values[] = { //Object
                           (int32)geo::eODDF_Axis,
                           (int32)geo::eODDF_Shape,
                           (int32)geo::eODDF_Embedding,
                           (int32)geo::eODDF_BV,
                           //Shape
                           (int32)geo::eSDDF_Boundary,
                           (int32)geo::eSDDF_Interior,
                           (int32)geo::eSDDF_Topology,
                           (int32)geo::eSDDF_Vertices,
                           (int32)geo::eSDDF_Edges,
                           (int32)geo::eSDDF_Faces,
                           (int32)geo::eSDDF_Volumes,
                           (int32)geo::eSDDF_Control,
                           (int32)geo::eSDDF_Normals,
                           (int32)geo::eSDDF_Tangents,
                           (int32)geo::eSDDF_FeatureId,
                           (int32)geo::eSDDF_DCR,
                           (int32)geo::eSDDF_BVH,
                           //BVH
                           (int32)geo::eBVHDDF_Geometry,
                           (int32)geo::eBVHDDF_Topology,
                           (int32)geo::eBVHDDF_Levels };
    al.AddProperty_Flags32( name, value,
                            20, vec_names, vec_values,
                            offset,
                            ntpf );
}

} //namespace geo
