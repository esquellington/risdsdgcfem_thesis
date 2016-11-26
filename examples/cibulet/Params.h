#ifndef CIBULET_PARAMS_H
#define CIBULET_PARAMS_H

#include <util/Archetype.h>
#include <Geo/IObject.h>

static util::ArchetypeLibrary s_Cibulet_ArchetypeLibrary;
static util::ItemStream s_ParamsIS( 4096, 4096 );

enum EDirtyFlags { eDF_None  = 0,
                   eDF_SRC   = 1<<0,
                   eDF_VIZ   = 1<<1,
                   eDF_EXT   = 1<<2,
                   eDF_INT   = 1<<3,
                   eDF_SIM   = 1<<4,
                   eDF_DCR   = 1<<5,
                   eDF_BVH   = 1<<6,
                   eDF_EMB   = 1<<7,
                   eDF_All   = 0xFFFFFFFF };

struct Params: public util::IArchetypeInstance
{
    float32 m_SRC_NaturalSize;

    float32 m_VIZ_Detail;
    float32 m_VIZ_DistortionScale;
    float32 m_VIZ_DistortionFreq;

    float32 m_EXT_Detail;
    float32 m_EXT_Offset;

    float32 m_SIM_CDT_Size;
    float32 m_SIM_CDT_Ratio;

    //\todo ????m_DCR_Method;

    enum EBVHMethod { eBVHM_BV_E = 0,
                      eBVHM_BV_BSLAB,
                      eBVHM_BV_BDOP,
                      cNumBVHM } m_BVH_Method;

    geo::EEmbeddingMethod m_EMB_Method;
    uint32 m_EMB_Subd_Iter;
    float32 m_EMB_DistortionScale;

    Flags32 m_DirtyFlags;

    inline Params()
    : m_SRC_NaturalSize(1)
    , m_VIZ_Detail(0.1f)
    , m_VIZ_DistortionScale(0)//0.015)
    , m_VIZ_DistortionFreq(20)
    , m_EXT_Detail(0.2f)
    , m_EXT_Offset(0)
    , m_SIM_CDT_Size(0.25)
    , m_SIM_CDT_Ratio(0.125)
    , m_BVH_Method( eBVHM_BV_E )
    , m_EMB_Method(geo::eEM_Barycentric)
    , m_EMB_Subd_Iter(0)
    , m_EMB_DistortionScale(0)
    , m_DirtyFlags( eDF_None )
        {}

    void SetName( const char *name ) {}
    const char *GetName() const { return "cibulet::Params"; }
    static void InitArchetype( util::ArchetypeLibrary &al );
    static bool NCPF_SRC( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_SRC ); return true; }
    static bool NCPF_VIZ( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_VIZ ); return true; }
    static bool NCPF_EXT( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_EXT ); return true; }
    static bool NCPF_SIM( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_SIM ); return true; }
    static bool NCPF_DCR( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_DCR ); return true; }
    static bool NCPF_BVH( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_BVH ); return true; }
    static bool NCPF_EMB( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_EMB ); return true; }
};

void Params::InitArchetype( util::ArchetypeLibrary &al )
{
    Params params;
    al.BeginArchetype( "Archetype_Cibulet_Params" );
    {
        al.BeginProperty_Group( "<SRC>" );
        {
            al.AddProperty_NIR( "natural_scale", params.m_SRC_NaturalSize, 0.1f, 10.0f, archetype_offset_of(params,m_SRC_NaturalSize), Params::NCPF_SRC );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<VIZ>" );
        {
            al.AddProperty_NIR( "detail", params.m_VIZ_Detail, 0.001f, 0.5f, archetype_offset_of(params,m_VIZ_Detail), Params::NCPF_VIZ );
            al.AddProperty_NIR( "distortion.Scale", params.m_VIZ_DistortionScale, 0.0f, 0.1f, archetype_offset_of(params,m_VIZ_DistortionScale), Params::NCPF_VIZ );
            al.AddProperty_NIR( "distortion.Freq", params.m_VIZ_Detail, 0.0f, 100.0f, archetype_offset_of(params,m_VIZ_DistortionFreq), Params::NCPF_VIZ );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<EXT>" );
        {
            al.AddProperty_NIR( "detail", params.m_EXT_Detail, 0.001f, 0.5f, archetype_offset_of(params,m_EXT_Detail), Params::NCPF_EXT );
            al.AddProperty_NIR( "offset", params.m_EXT_Offset, 0.0f, 0.25f, archetype_offset_of(params,m_EXT_Offset), Params::NCPF_EXT );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<SIM>" );
        {
            al.AddProperty_NIR( "cdt_size", params.m_SIM_CDT_Size, 0.05f, 1.0f, archetype_offset_of(params,m_SIM_CDT_Size), Params::NCPF_SIM );
            al.AddProperty_NIR( "cdt_ratio", params.m_SIM_CDT_Ratio, 0.0f, 0.5f, archetype_offset_of(params,m_SIM_CDT_Ratio), Params::NCPF_SIM );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<DCR>" );
        {
            //\todo No params by now...
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<BVH>" );
        {
            const char *vec_names_bvhm[] = { "BV(E)", "BV(BSlab)", "BV(BDOP)" };
            uint32 vec_values_bvhm[] = { (uint32)eBVHM_BV_E,
                                         (uint32)eBVHM_BV_BSLAB,
                                         (uint32)eBVHM_BV_BDOP };
            al.AddProperty_Enum32( "method", (uint32)params.m_BVH_Method,
                                   cNumBVHM, vec_names_bvhm, vec_values_bvhm,
                                   archetype_offset_of(params,m_BVH_Method), Params::NCPF_BVH );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<EMB>" );
        {
            const char *vec_names_em[] = { "B", "MLS", "BNDF" };
            uint32 vec_values_em[] = { (uint32)geo::eEM_Barycentric,
                                       (uint32)geo::eEM_MLS,
                                       (uint32)geo::eEM_BNDF };
            al.AddProperty_Enum32( "method", (uint32)params.m_EMB_Method,
                                   geo::cNumEM, vec_names_em, vec_values_em,
                                   archetype_offset_of(params,m_EMB_Method), Params::NCPF_EMB );
            al.AddProperty_NIR<uint32>( "#subd", params.m_EMB_Subd_Iter, 0, 10, archetype_offset_of(params,m_EMB_Subd_Iter), Params::NCPF_EMB );
            al.AddProperty_NIR( "distortion", params.m_EMB_DistortionScale, 0.0f, 0.1f, archetype_offset_of(params,m_EMB_DistortionScale), Params::NCPF_EMB );
        }
        al.EndProperty_Group();
    }
    al.EndArchetype();
}

#endif //CIBULET_PARAMS_H
