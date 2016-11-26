#ifndef ESCALUNYA_PARAMS_H
#define ESCALUNYA_PARAMS_H

#include <util/Archetype.h>
#include <Geo/IObject.h>
#include <Geo/util/Viz.h>
#include <Geo/geo_properties.h>

static util::ArchetypeLibrary s_Escalunya_ArchetypeLibrary;
static util::ItemStream s_ParamsIS( 4096, 4096 );

enum EDirtyFlags { eDF_None   = 0,
                   eDF_SRC    = 1<<0,
                   eDF_VIZ    = 1<<1,
                   eDF_EXT    = 1<<2,
                   eDF_INT    = 1<<3,
                   eDF_SIM    = 1<<4,
                   eDF_CLP    = 1<<5,
                   eDF_DCR    = 1<<6,
                   eDF_BVH    = 1<<7,
                   eDF_EMB    = 1<<8,
                   eDF_Redraw = 1<<9,
                   eDF_All    = 0xFFFFFFFF };

struct Params: public util::IArchetypeInstance
{
    float32 m_SRC_NaturalSize;
    Flags32 m_SRC_DDF;

    float32 m_VIZ_Detail;
    Flags32 m_VIZ_DDF;

    float32 m_EXT_Detail;
    float32 m_EXT_Offset;
    Flags32 m_EXT_DDF;

    // CDT params
    float32 m_SIM_CDT_Facet_Angle;
    float32 m_SIM_CDT_Facet_Size;
    float32 m_SIM_CDT_Facet_Distance;
    float32 m_SIM_CDT_Cell_Size;
    float32 m_SIM_CDT_Cell_Ratio;
    bool m_SIM_CDT_Lloyd;
    bool m_SIM_CDT_Odt;
    bool m_SIM_CDT_Perturb;
    bool m_SIM_CDT_Exude;

    // Fast Tetrahedralization params
    float32 m_SIM_Fit_ODT_RelaxationCoeff;
    float32 m_SIM_Fit_Lpc_RelaxationCoeff;
    uint32 m_SIM_Fit_MaxIter;
    bool m_SIM_Fit_FixNMF;
    // Other
    Flags32 m_SIM_DDF;

    uint32 m_CLP_MaxIter_CloseHoles;
    float32 m_CLP_SliverThreshold_Log10;
    float32 m_CLP_EpsilonLength;
    Flags32 m_CLP_DDF;

    enum EDCRSource { eDCRS_VIZ = 0,
                      eDCRS_CLP,
                      cNumDCRS } m_DCR_Source;

    enum EBVHMethod { eBVHM_BV_E = 0,
                      eBVHM_BV_BSLAB,
                      eBVHM_BV_BDOP,
                      cNumBVHM } m_BVH_Method;

    enum EEMBSource { eEMBS_VIZ = 0,
                      eEMBS_CLP,
                      cNumEMBS } m_EMB_Source;
    geo::EEmbeddingMethod m_EMB_Method;
    float32 m_EMB_DistortionScale;
    Flags32 m_EMB_DDF;

    Flags32 m_DirtyFlags;

    inline Params()
    : m_SRC_NaturalSize(1)
    , m_SRC_DDF( geo::eODDF_Default | geo::eSDDF_Default )
    , m_VIZ_Detail(0.1f)
    , m_VIZ_DDF( geo::eODDF_Default | geo::eSDDF_Default )
    , m_EXT_Detail(1)
    , m_EXT_Offset(0)
    , m_EXT_DDF( geo::eODDF_Default | geo::eSDDF_Default )
    , m_SIM_CDT_Facet_Angle(25) //degrees?!
    , m_SIM_CDT_Facet_Size(0.25)
    , m_SIM_CDT_Facet_Distance(0.008)
    , m_SIM_CDT_Cell_Size(0.25)
    , m_SIM_CDT_Cell_Ratio(3) //??
    , m_SIM_CDT_Lloyd(false), m_SIM_CDT_Odt(false), m_SIM_CDT_Perturb(false), m_SIM_CDT_Exude(false)
    , m_SIM_Fit_ODT_RelaxationCoeff(1), m_SIM_Fit_Lpc_RelaxationCoeff(0), m_SIM_Fit_MaxIter(0), m_SIM_Fit_FixNMF(true)
    , m_SIM_DDF(geo::eODDF_Shape|geo::eSDDF_Boundary)
    , m_CLP_MaxIter_CloseHoles(1000)
    , m_CLP_SliverThreshold_Log10(mal::Log10(25.0f))
    , m_CLP_EpsilonLength(0)
    , m_CLP_DDF(geo::eODDF_Default|geo::eSDDF_Default|geo::eSDDF_Interior) //\todo interior to see slivers
    , m_DCR_Source(eDCRS_CLP)
    , m_BVH_Method( eBVHM_BV_E )
    , m_EMB_Source(eEMBS_VIZ)
    , m_EMB_Method(geo::eEM_Barycentric)
    , m_EMB_DistortionScale(0)
    , m_EMB_DDF( geo::eODDF_Default | geo::eSDDF_Default )
    , m_DirtyFlags( eDF_None )
        {}

    void SetName( const char *name ) {}
    const char *GetName() const { return "escalunya::Params"; }
    static void InitArchetype( util::ArchetypeLibrary &al );
    static bool NCPF_SRC( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_SRC ); return true; }
    static bool NCPF_VIZ( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_VIZ ); return true; }
    static bool NCPF_EXT( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_EXT ); return true; }
    static bool NCPF_SIM( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_SIM ); return true; }
    static bool NCPF_CLP( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_CLP ); return true; }
    static bool NCPF_DCR( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_DCR ); return true; }
    static bool NCPF_BVH( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_BVH ); return true; }
    static bool NCPF_EMB( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_EMB ); return true; }
    static bool NCPF_Redraw( IArchetypeInstance *p_instance ) { static_cast<Params*>(p_instance)->m_DirtyFlags.Enable( eDF_Redraw ); return true; }
};

void Params::InitArchetype( util::ArchetypeLibrary &al )
{
    Params params;
    al.BeginArchetype( "Archetype_Escalunya_Params" );
    {
        al.BeginProperty_Group( "<SRC>" );
        {
            al.AddProperty_NIR( "natural_scale", params.m_SRC_NaturalSize, 0.1f, 5.0f, archetype_offset_of(params,m_SRC_NaturalSize), Params::NCPF_SRC );
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_SRC_DDF, archetype_offset_of(params,m_SRC_DDF), Params::NCPF_Redraw );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<VIZ>" );
        {
            al.AddProperty_NIR( "detail", params.m_VIZ_Detail, 0.001f, 1.0f, archetype_offset_of(params,m_VIZ_Detail), Params::NCPF_VIZ );
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_VIZ_DDF, archetype_offset_of(params,m_VIZ_DDF), Params::NCPF_Redraw );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<EXT>" );
        {
            al.AddProperty_NIR( "detail", params.m_EXT_Detail, 0.001f, 1.0f, archetype_offset_of(params,m_EXT_Detail), Params::NCPF_EXT );
            al.AddProperty_NIR( "offset", params.m_EXT_Offset, 0.0f, 0.25f, archetype_offset_of(params,m_EXT_Offset), Params::NCPF_EXT );
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_EXT_DDF, archetype_offset_of(params,m_EXT_DDF), Params::NCPF_Redraw );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<SIM>" );
        {
            /*TEMP: disabled, no longer used... consider removing all related code or keep as reference of an unsuccessful approach...
            al.BeginProperty_Group( "<CDT>" );
            {
                al.AddProperty_NIR( "facet_angle", params.m_SIM_CDT_Facet_Angle, 0.0f, 90.0f, archetype_offset_of(params,m_SIM_CDT_Facet_Angle), Params::NCPF_SIM );
                al.AddProperty_NIR( "facet_size", params.m_SIM_CDT_Facet_Size, 0.0f, 1.0f, archetype_offset_of(params,m_SIM_CDT_Facet_Size), Params::NCPF_SIM );
                al.AddProperty_NIR( "facet_distance", params.m_SIM_CDT_Facet_Distance, 0.0f, 0.1f, archetype_offset_of(params,m_SIM_CDT_Facet_Distance), Params::NCPF_SIM );
                al.AddProperty_NIR( "cell_size", params.m_SIM_CDT_Cell_Size, 0.05f, 1.0f, archetype_offset_of(params,m_SIM_CDT_Cell_Size), Params::NCPF_SIM );
                al.AddProperty_NIR( "cell_ratio", params.m_SIM_CDT_Cell_Ratio, 0.0f, 5.0f, archetype_offset_of(params,m_SIM_CDT_Cell_Ratio), Params::NCPF_SIM );
                al.AddProperty( "?Lloyd", params.m_SIM_CDT_Lloyd, archetype_offset_of(params,m_SIM_CDT_Lloyd), Params::NCPF_SIM );
                al.AddProperty( "?odt", params.m_SIM_CDT_Odt, archetype_offset_of(params,m_SIM_CDT_Odt), Params::NCPF_SIM );
                al.AddProperty( "?perturb", params.m_SIM_CDT_Perturb, archetype_offset_of(params,m_SIM_CDT_Perturb), Params::NCPF_SIM );
                al.AddProperty( "?exude", params.m_SIM_CDT_Exude, archetype_offset_of(params,m_SIM_CDT_Exude), Params::NCPF_SIM );
            }
            al.EndProperty_Group();
            */
            al.BeginProperty_Group( "<FIT>" );
            {
                al.AddProperty_NIR( "cell_size", params.m_SIM_CDT_Cell_Size, 0.05f, 1.0f, archetype_offset_of(params,m_SIM_CDT_Cell_Size), Params::NCPF_SIM );
                al.AddProperty_NIR( "ODT coeff", params.m_SIM_Fit_ODT_RelaxationCoeff, 0.0f, 1.0f, archetype_offset_of(params,m_SIM_Fit_ODT_RelaxationCoeff), Params::NCPF_SIM );
                al.AddProperty_NIR( "Lpc coeff", params.m_SIM_Fit_Lpc_RelaxationCoeff, 0.0f, 1.0f, archetype_offset_of(params,m_SIM_Fit_Lpc_RelaxationCoeff), Params::NCPF_SIM );
                al.AddProperty_NIR( "max_iter", params.m_SIM_Fit_MaxIter, (uint32)0, (uint32)20, archetype_offset_of(params,m_SIM_Fit_MaxIter), Params::NCPF_SIM );
                al.AddProperty( "Fix NMF", params.m_SIM_Fit_FixNMF, archetype_offset_of(params,m_SIM_Fit_FixNMF), Params::NCPF_SIM );
            }
            al.EndProperty_Group();
            //\todo Consider moving to geo/viz or elsewhere...
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_SIM_DDF, archetype_offset_of(params,m_SIM_DDF), Params::NCPF_Redraw );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<CLP>" );
        {
            al.AddProperty_NIR( "max_iter CH", params.m_CLP_MaxIter_CloseHoles, (uint32)0, (uint32)1000, archetype_offset_of(params,m_CLP_MaxIter_CloseHoles), Params::NCPF_CLP );
            al.AddProperty_NIR( "Log10(sliver_ratio)", params.m_CLP_SliverThreshold_Log10, 1.0f, 3.0f, archetype_offset_of(params,m_CLP_SliverThreshold_Log10), Params::NCPF_CLP );
            al.AddProperty_NIR( "eps.length CH", params.m_CLP_EpsilonLength, 0.0f, 0.1f, archetype_offset_of(params,m_CLP_EpsilonLength), Params::NCPF_CLP );
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_CLP_DDF, archetype_offset_of(params,m_CLP_DDF), Params::NCPF_Redraw );
        }
        al.EndProperty_Group();
        al.BeginProperty_Group( "<DCR>" );
        {
            const char *vec_names_dcrs[] = { "VIZ", "CLP" };
            uint32 vec_values_dcrs[] = { (uint32)eDCRS_VIZ,
                                       (uint32)eDCRS_CLP  };
            al.AddProperty_Enum32( "source", (uint32)params.m_DCR_Source,
                                   cNumDCRS, vec_names_dcrs, vec_values_dcrs,
                                   archetype_offset_of(params,m_DCR_Source), Params::NCPF_DCR );
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
            const char *vec_names_es[] = { "VIZ", "CLP" };
            uint32 vec_values_es[] = { (uint32)eEMBS_VIZ,
                                       (uint32)eEMBS_CLP  };
            al.AddProperty_Enum32( "source", (uint32)params.m_EMB_Source,
                                   cNumEMBS, vec_names_es, vec_values_es,
                                   archetype_offset_of(params,m_EMB_Source), Params::NCPF_EMB );
            const char *vec_names_em[] = { "B", "MLS", "BNDF" };
            uint32 vec_values_em[] = { (uint32)geo::eEM_Barycentric,
                                       (uint32)geo::eEM_MLS,
                                       (uint32)geo::eEM_BNDF };
            al.AddProperty_Enum32( "method", (uint32)params.m_EMB_Method,
                                   geo::cNumEM, vec_names_em, vec_values_em,
                                   archetype_offset_of(params,m_EMB_Method), Params::NCPF_EMB );
            al.AddProperty_NIR( "distortion", params.m_EMB_DistortionScale, 0.0f, 0.1f, archetype_offset_of(params,m_EMB_DistortionScale), Params::NCPF_Redraw );
            geo::ArchetypeLibrary_AddProperty_DDF( al, "DDF", params.m_EMB_DDF, archetype_offset_of(params,m_EMB_DDF), Params::NCPF_Redraw );
        }
        al.EndProperty_Group();
    }
    al.EndArchetype();
}

#endif //ESCALUNYA_PARAMS_H
