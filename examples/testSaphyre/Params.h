#ifndef TEST_SAPHYRE_PARAMS_H
#define TEST_SAPHYRE_PARAMS_H

#include "Config.h"
#include <util/Archetype.h>

class Params: public util::IArchetypeInstance
{
public:

    struct SceneParams
    {
        bool m_bLog = true;
        bool m_bViz = true;
        bool m_bProf = false;
    } scene;

    enum ERenderMethod {
        eRM_None  = 0,
        eRM_Wire  = 1,
        eRM_Solid = 2,
        eRM_Flat  = 3,
        cNumRenderMethods,
        eRM_Default = eRM_Wire
    };

    struct SceneRender
    {
        enum EDrawFlags {
            eDF_None   = 0,
            eDF_Axis   = (1<<0),
            eDF_GridX  = (1<<1),
            eDF_GridY  = (1<<2),
            eDF_GridZ  = (1<<3),
            eDF_Coords = (1<<4),
            eDF_DCR    = (1<<5),
            eDF_BVH    = (1<<6),
            eDF_Default = (eDF_DCR),
            eDF_Grids  = (eDF_GridX | eDF_GridY | eDF_GridZ),
            eDF_All    = 0xFFFFFFFF,
            cNumDrawFlags = 7
        };
        Flags32 m_DF = eDF_Default;
        ERenderMethod m_GroundRM = eRM_Flat;
        ERenderMethod m_KineRM = eRM_Wire;
        ERenderMethod m_DynRM = eRM_Wire;
        ERenderMethod m_EmbRM = eRM_Flat;
        ERenderMethod m_KineDCRRM = eRM_Flat;
        ERenderMethod m_DynDCRRM = eRM_None;
        bool m_DCR_UsePatchColor = true;
        bool m_DCR_BoundaryOnly = true;
    } scene_render;

    struct ContactRender
    {
        enum EDrawFlags {
            eDF_None                      = 0,
            eDF_SelectedObjects           = 1<<0,
            eDF_ContactPoints             = 1<<1,
            eDF_ContactNormal             = 1<<2,
            eDF_DCR2Primitive_IC          = 1<<3,
            eDF_DCR2DCR_IC                = 1<<4,
            eDF_DCR2DCR_IC_UnconnectedCTP = 1<<5,
            eDF_DCR2DCR_IC_SeedHE         = 1<<6,
            eDF_DCR2DCR_IB_FloodT         = 1<<7,
            eDF_DCR2DCR_IB_FloodP         = 1<<8,
            eDF_DCR2DCR_IB_PosNormal      = 1<<9,
            eDF_DCR2DCR_IB_CE             = 1<<10,
            eDF_DCR2DCR_IB_CT             = 1<<11,
            cNumDrawFlags                 = 12,
            eDF_Default = ( eDF_SelectedObjects
                            | eDF_ContactPoints | eDF_ContactNormal
                            | eDF_DCR2Primitive_IC
                            | eDF_DCR2DCR_IC
                            | eDF_DCR2DCR_IC_UnconnectedCTP
                            | eDF_DCR2DCR_IC_SeedHE
                            | eDF_DCR2DCR_IB_FloodT
                            | eDF_DCR2DCR_IB_FloodP
                            | eDF_DCR2DCR_IB_PosNormal
                            | eDF_DCR2DCR_IB_CE
                            | eDF_DCR2DCR_IB_CT )
        };
        Flags32 m_DF = eDF_Default;
        ERenderMethod m_FloodRM = eRM_Wire;
    } contact_render;

public:
    Params() {}
    ~Params() {}

    //!\name IArchetypeInstance implementation
    //@{
    inline bool Rebuild() { return true; }
    inline void SetName( const char *name ) {}
    inline const char *GetName() const { return "AppParams"; }
    static void InitArchetype( util::ArchetypeLibrary &al )
        {
            al.BeginArchetype( "Archetype_AP" );
            {
                Params ap;
                al.BeginProperty_Group( "<Scene>" );
                {
                    al.AddProperty( "Log?", ap.scene.m_bLog, archetype_offset_of(ap,scene.m_bLog) );
                    al.AddProperty( "Viz?", ap.scene.m_bViz, archetype_offset_of(ap,scene.m_bViz) );
                    al.AddProperty( "Prof?", ap.scene.m_bProf, archetype_offset_of(ap,scene.m_bProf) );
                }
                al.EndProperty_Group();
                al.BeginProperty_Group( "<Scn Render>" );
                {
                    const char *vec_names_df[] = { "Axis", "G.X", "G.Y", "G.Z", "Cords", "DCR", "BVH" };
                    int32 vec_values_df[] = { (int32)SceneRender::eDF_Axis,
                                              (int32)SceneRender::eDF_GridX,
                                              (int32)SceneRender::eDF_GridY,
                                              (int32)SceneRender::eDF_GridZ,
                                              (int32)SceneRender::eDF_Coords,
                                              (int32)SceneRender::eDF_DCR,
                                              (int32)SceneRender::eDF_BVH };
                    al.AddProperty_Flags32( "DF", ap.scene_render.m_DF,
                                            SceneRender::cNumDrawFlags, vec_names_df, vec_values_df,
                                            archetype_offset_of(ap,scene_render.m_DF),
                                            IArchetypeInstance::NTPF_Ignore );
                    const char *vec_names_rm[] = { "None", "Wire", "Solid", "Flat" };
                    uint32 vec_values_rm[] = { eRM_None, eRM_Wire, eRM_Solid, eRM_Flat };
                    al.AddProperty_Enum32( "RM.Ground", (uint32)ap.scene_render.m_GroundRM,
                                           cNumRenderMethods, vec_names_rm, vec_values_rm,
                                           archetype_offset_of(ap,scene_render.m_GroundRM),
                                           IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_Enum32( "RM.Kine", (uint32)ap.scene_render.m_KineRM,
                                           cNumRenderMethods, vec_names_rm, vec_values_rm,
                                           archetype_offset_of(ap,scene_render.m_KineRM),
                                           IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_Enum32( "RM.Dyn", (uint32)ap.scene_render.m_DynRM,
                                           cNumRenderMethods, vec_names_rm, vec_values_rm,
                                           archetype_offset_of(ap,scene_render.m_DynRM),
                                           IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_Enum32( "RM.Emb", (uint32)ap.scene_render.m_EmbRM,
                                           cNumRenderMethods, vec_names_rm, vec_values_rm,
                                           archetype_offset_of(ap,scene_render.m_EmbRM),
                                           IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_Enum32( "RM.KineDCR", (uint32)ap.scene_render.m_KineDCRRM,
                                           cNumRenderMethods, vec_names_rm, vec_values_rm,
                                           archetype_offset_of(ap,scene_render.m_KineDCRRM),
                                           IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty_Enum32( "RM.DynDCR", (uint32)ap.scene_render.m_DynDCRRM,
                                           cNumRenderMethods, vec_names_rm, vec_values_rm,
                                           archetype_offset_of(ap,scene_render.m_DynDCRRM),
                                           IArchetypeInstance::NTPF_Ignore );
                    al.AddProperty( "DCR.PatchColors?", ap.scene_render.m_DCR_UsePatchColor, archetype_offset_of(ap,scene_render.m_DCR_UsePatchColor) );
                    al.AddProperty( "DCR.BoundaryOnly?", ap.scene_render.m_DCR_BoundaryOnly, archetype_offset_of(ap,scene_render.m_DCR_BoundaryOnly) );
                }
                al.EndProperty_Group();

                al.BeginProperty_Group( "<CD Render>" );
                {
                    const char *vec_names_df[] = { "Obj", "CP", "CN", "PIC", "IC", "UCP", "SHE", "FT", "FP", "IBPN", "CE", "CT" };
                    int32 vec_values_df[] = {  (int32)ContactRender::eDF_SelectedObjects,
                                               (int32)ContactRender::eDF_ContactPoints,
                                               (int32)ContactRender::eDF_ContactNormal,
                                               (int32)ContactRender::eDF_DCR2Primitive_IC,
                                               (int32)ContactRender::eDF_DCR2DCR_IC,
                                               (int32)ContactRender::eDF_DCR2DCR_IC_UnconnectedCTP,
                                               (int32)ContactRender::eDF_DCR2DCR_IC_SeedHE,
                                               (int32)ContactRender::eDF_DCR2DCR_IB_FloodT,
                                               (int32)ContactRender::eDF_DCR2DCR_IB_FloodP,
                                               (int32)ContactRender::eDF_DCR2DCR_IB_PosNormal,
                                               (int32)ContactRender::eDF_DCR2DCR_IB_CE,
                                               (int32)ContactRender::eDF_DCR2DCR_IB_CT };
                    al.AddProperty_Flags32( "DF", ap.contact_render.m_DF,
                                            ContactRender::cNumDrawFlags, vec_names_df, vec_values_df,
                                            archetype_offset_of(ap,contact_render.m_DF),
                                            IArchetypeInstance::NTPF_Ignore );
                    const char *vec_names_rm_flood[] = { "None", "Wire", "Solid", "Flat" };
                    uint32 vec_values_rm_flood[] = { eRM_None, eRM_Wire, eRM_Solid, eRM_Flat };
                    al.AddProperty_Enum32( "RM.Flood", (uint32)ap.contact_render.m_FloodRM,
                                           cNumRenderMethods, vec_names_rm_flood, vec_values_rm_flood,
                                           archetype_offset_of(ap,contact_render.m_FloodRM),
                                           IArchetypeInstance::NTPF_Ignore );
                }
                al.EndProperty_Group();
            }
            al.EndArchetype();
        }
    //@}
};

#endif //TEST_SAPHYRE_PARAMS_H
