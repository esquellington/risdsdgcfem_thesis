#ifndef PLA_UTIL_VIZ_MACROS_H
#define PLA_UTIL_VIZ_MACROS_H

#ifdef __PLA_UTIL_ENABLE_VIZ_MACROS

//----- POINT2, POINT3
#define VIZ_POINT2( vs, pos2d, pen_size, style_color ) \
{                                                                  \
    vs.BeginComplex(1,util::eType_VizPoint);                       \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("pos", mal::CastDimension<3,float,2>( pos2d ) );  \
        /* Style */                                                \
        vs.BeginComplex("style",util::eType_VizStyle);             \
        {                                                          \
            vs.Write("color",Vec4f(style_color));                  \
            vs.Write("pen_size",float(pen_size));                  \
        }                                                          \
        vs.EndComplex();                                           \
     }                                                             \
     vs.EndComplex();                                              \
}

#define VIZ_POINT3( vs, pos3d, pen_size, style_color ) \
{                                                                  \
    vs.BeginComplex(1,util::eType_VizPoint);                       \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("pos", pos3d );                                   \
        /* Style */                                                \
        vs.BeginComplex("style",util::eType_VizStyle);             \
        {                                                          \
            vs.Write("color",Vec4f(style_color));                  \
            vs.Write("pen_size",float(pen_size));                  \
        }                                                          \
        vs.EndComplex();                                           \
     }                                                             \
     vs.EndComplex();                                              \
}

//---- TRIANGLE2
#define VIZ_TRIANGLE2( vs, p1, p2, p3, pen_size, style_color, style_flags ) \
{ \
    vs.BeginComplex(1,util::eType_VizTriangle); \
    { \
        /* Geometry */ \
        vs.Write("p1", mal::CastDimension<3,float,2>( p1 ) ); \
        vs.Write("p2", mal::CastDimension<3,float,2>( p2 ) ); \
        vs.Write("p3", mal::CastDimension<3,float,2>( p3 ) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}

//---- TRIANGLE2
#define VIZ_TRIANGLE3( vs, p1, p2, p3, pen_size, style_color, style_flags ) \
{ \
    vs.BeginComplex(1,util::eType_VizTriangle); \
    { \
        /* Geometry */ \
        vs.Write("p1", Vec3f( p1 ) ); \
        vs.Write("p2", Vec3f( p2 ) ); \
        vs.Write("p3", Vec3f( p3 ) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}

//---- DISK2
#define VIZ_DISK2( vs, pos2d, quat, radius, style_color, style_flags )  \
{                                                                  \
    vs.BeginComplex(1,util::eType_VizDisk);                         \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("pos", mal::CastDimension<3,float,2>( pos2d ) );  \
        vs.Write("rot", Quatf(quat) );                             \
        vs.Write("radius", (float)radius );                        \
        /* Style */                                                \
        vs.BeginComplex("style",util::eType_VizStyle);              \
        {                                                          \
            vs.Write("color",Vec4f(style_color));                  \
            vs.Write("flags",Flags32(style_flags));                \
        }                                                          \
        vs.EndComplex();                                           \
     }                                                             \
     vs.EndComplex();                                              \
}

#define VIZ_DISK2_NO_ROT( vs, pos2d, radius, style_color, style_flags ) \
    VIZ_DISK2(vs,pos2d,Quatf::Identity(),radius,style_color,style_flags)

//---- SPHERE3
#define VIZ_SPHERE3( vs, pos3d, quat, radius, style_color, style_flags ) \
{                                                                  \
    vs.BeginComplex(1,util::eType_VizSphere);                      \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("pos", Vec3f(pos3d) );                            \
        vs.Write("rot", Quatf(quat) ); \
        vs.Write("radius", (float)radius );                        \
        /* Style */                                                \
        vs.BeginComplex("style",util::eType_VizStyle);             \
        {                                                          \
            vs.Write("color",Vec4f(style_color));                  \
            vs.Write("flags",Flags32(style_flags));                \
        }                                                          \
        vs.EndComplex();                                           \
     }                                                             \
     vs.EndComplex();                                              \
}

#define VIZ_SPHERE3_NO_ROT( vs, pos3d, radius, style_color, style_flags ) \
    VIZ_SPHERE3(vs,pos3d,Quatf::Identity(),radius,style_color,style_flags)

//---- RECT2
#define VIZ_RECT2( vs, pos2d, quat, half_sizes2d, style_color, style_flags )  \
{ \
    vs.BeginComplex(1,util::eType_VizRectangle); \
    { \
        /* Geometry */ \
        vs.Write("pos", mal::CastDimension<3,float,2>( pos2d ) ); \
        vs.Write("rot", Quatf(quat) ); \
        vs.Write("half_sizes", Vec2f(half_sizes2d) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}

#define VIZ_RECT2_NO_ROT( vs, pos2d, half_sizes2d, style_color, style_flags ) \
    VIZ_RECT2(vs,pos2d,Quatf::Identity(),half_sizes2d,style_color,style_flags)

//---- BOX3
#define VIZ_BOX3( vs, pos3d, quat, half_sizes3d, style_color, style_flags )  \
{ \
    vs.BeginComplex(1,util::eType_VizBox); \
    { \
        /* Geometry */ \
        vs.Write("pos", Vec3f(pos3d) ); \
        vs.Write("rot", Quatf(quat) ); \
        vs.Write("half_sizes", Vec3f(half_sizes3d) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}

#define VIZ_BOX3_NO_ROT( vs, pos3d, half_sizes3d, style_color, style_flags ) \
    VIZ_BOX3(vs,pos3d,Quatf::Identity(),half_sizes3d,style_color,style_flags)

//---- TETRAHEDRON3
#define VIZ_TETRAHEDRON3( vs, p0, p1, p2, p3, style_color, style_flags ) \
{ \
    vs.BeginComplex(1,util::eType_VizTetrahedron); \
    { \
        /* Geometry */ \
        vs.Write("p0", Vec3f( p0 ) ); \
        vs.Write("p1", Vec3f( p1 ) ); \
        vs.Write("p2", Vec3f( p2 ) ); \
        vs.Write("p3", Vec3f( p3 ) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}

#define VIZ_QUAD3( vs, p1, p2, p3, p4, pen_size, style_color, style_flags ) \
{ \
    vs.BeginComplex(1,util::eType_VizQuad); \
    { \
        /* Geometry */ \
        vs.Write("p1", Vec3f( p1 )); \
        vs.Write("p2", Vec3f( p2 )); \
        vs.Write("p3", Vec3f( p3 )); \
        vs.Write("p4", Vec3f( p4 )); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}

/*
#define VIZ_PLANE3( vs, normal, coeff_d, style_color, style_flags ) \
{ \
    vs.BeginComplex(1,util::eType_VizPlane); \
    { \
        // Geometry                            \
        vs.Write("normal", Vec3f( normal ) ); \
        vs.Write("coeff_d", float( coeff_d ) ); \
        // Style                                       \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("flags",Flags32(style_flags)); \
        } \
        vs.EndComplex(); \
     } \
     vs.EndComplex(); \
}
*/

//---- SEGMENT2, SEGMENT3
#define VIZ_SEGMENT2( vs, p0, p1, pen_size, style_color ) \
{ \
    vs.BeginComplex(1,util::eType_VizSegment); \
    { \
        /* Geometry */ \
        vs.Write("pos1", mal::CastDimension<3,float,2>(p0) ); \
        vs.Write("pos2", mal::CastDimension<3,float,2>(p1) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
    } \
    vs.EndComplex(); \
}

#define VIZ_SEGMENT3( vs, p0, p1, pen_size, style_color ) \
{ \
    vs.BeginComplex(1,util::eType_VizSegment); \
    { \
        /* Geometry */ \
        vs.Write("pos1", Vec3f(p0) );                         \
        vs.Write("pos2", Vec3f(p1) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
    } \
    vs.EndComplex(); \
}

//---- VECTOR2, VECTOR3
#define VIZ_VECTOR2( vs, pos2d, vec2d, pen_size, style_color ) \
{ \
    vs.BeginComplex(1,util::eType_VizVec); \
    { \
        /* Geometry */ \
        vs.Write("pos", mal::CastDimension<3,float,2>(pos2d) ); \
        vs.Write("vec", mal::CastDimension<3,float,2>(vec2d) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
    } \
    vs.EndComplex(); \
}

#define VIZ_VECTOR3( vs, pos3d, vec3d, pen_size, style_color ) \
{ \
    vs.BeginComplex(1,util::eType_VizVec); \
    { \
        /* Geometry */ \
        vs.Write("pos", Vec3f(pos3d) ); \
        vs.Write("vec", Vec3f(vec3d) ); \
        /* Style */ \
        vs.BeginComplex("style",util::eType_VizStyle); \
        { \
            vs.Write("color",Vec4f(style_color)); \
            vs.Write("pen_size",float(pen_size)); \
        } \
        vs.EndComplex(); \
    } \
    vs.EndComplex(); \
}

#define VIZ_TRANSFORM2( vs, transform, scale )  \
{                                                                  \
    vs.BeginComplex(1,util::eType_VizTransform2);                  \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("transform", Transform2f(transform) );            \
        vs.Write("scale", float(scale) );                          \
    }                                                              \
    vs.EndComplex();                                               \
}

#define VIZ_TRANSFORM3( vs, transform, scale )  \
{                                                                  \
    vs.BeginComplex(1,util::eType_VizTransform3);                  \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("transform", Transform3f(transform) );            \
        vs.Write("scale", float(scale) );                          \
    }                                                              \
    vs.EndComplex();                                               \
}

#define VIZ_POINT2_NAMED( vs, name, pos2d, pen_size, style_color ) \
{                                                                  \
    vs.BeginComplex(name,util::eType_VizPoint);                    \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("pos", mal::CastDimension<3,float,2>( pos2d ) );  \
        /* Style */                                                \
        vs.BeginComplex("style",util::eType_VizStyle);             \
        {                                                          \
            vs.Write("color",Vec4f(style_color));                  \
            vs.Write("pen_size",float(pen_size));                  \
        }                                                          \
        vs.EndComplex();                                           \
     }                                                             \
     vs.EndComplex();                                              \
}

#define VIZ_POINT3_NAMED( vs, name, pos, pen_size, style_color )   \
{                                                                  \
    vs.BeginComplex(name,util::eType_VizPoint);                    \
    {                                                              \
        /* Geometry */                                             \
        vs.Write("pos", pos );                                     \
        /* Style */                                                \
        vs.BeginComplex("style",util::eType_VizStyle);             \
        {                                                          \
            vs.Write("color",Vec4f(style_color));                  \
            vs.Write("pen_size",float(pen_size));                  \
        }                                                          \
        vs.EndComplex();                                           \
     }                                                             \
     vs.EndComplex();                                              \
}

#else //__PLA_UTIL_ENABLE_VIZ_MACROS

#endif //__PLA_UTIL_ENABLE_VIZ_MACROS

#endif // PLA_UTIL_VIZ_MACROS_H
