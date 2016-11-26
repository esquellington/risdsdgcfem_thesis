#ifndef GEO_GEO_PROPERTIES_H
#define GEO_GEO_PROPERTIES_H

#include <util/Archetype.h>

//\todo This file may be split into per-namespace properties as happens with stats and params... bv/properties, util/properties, np/properties, etc... by now it's just a draft

namespace geo {

// Add a DDF property to an archetype being edited in an AL
void ArchetypeLibrary_AddProperty_DDF( util::ArchetypeLibrary& al, const char* name,
                                       Flags32 value, machine_uint_type offset,
                                       util::IArchetypeInstance::notify_touched_property_function_type ntpf = util::IArchetypeInstance::NTPF_Rebuild );
} //namespace geo

#endif //GEO_GEO_PROPERTIES_H
