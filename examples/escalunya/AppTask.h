#ifndef ESCALUNYA_APP_TASK_H
#define ESCALUNYA_APP_TASK_H

#include "Params.h"
#include <Safra/Safra.h>
#include <Safra/gui/CameraControllerML.h>
#include <Safra/gfx/VizRenderer.h>

#include <boost/bind.hpp>

#include <iostream>
#include <fstream>
#include <string>

#include "OFF_utils.h"
#include "CGAL_utils.h"
#include "TET_utils.h"

#include <Geo/IObject.h>
#include <Geo/shape/shape.h>
#include <Geo/util/Viz.h>
#include <Geo/shape/ShapeLibrary.h>
#include <Geo/shape/ShapeFactory.h>

#include <util/Archetype.h>
#include <util/VizMacros.h>

//#define __ENABLE_BVH_TEST
#ifdef __ENABLE_BVH_TEST
#  define __USE_BVH_AABB
//\todo #  define __USE_BVH_KDOP
//#  define __USE_BVH_SPHERE
#  include <Geo/bv/GBoundingVolumeHierarchy.h> //TEMPORAL!!
#endif

#include <Geo/geo.h>
//#define __ENABLE_CLIP_TEST
#ifdef __ENABLE_CLIP_TEST
#  include <Geo/np/Clip.h>
#endif

extern float g_TTS_SliverThreshold;

class AppTask: public sfr::IUpdateTarget
{
public:
    Params m_Params;

#ifdef __ENABLE_BVH_TEST
#  ifdef __USE_BVH_AABB
    typedef geo::bv::BVH_ST_DG_AABB3 bvh_type;
#  elif defined(__USE_BVH_KDOP)
\todo    typedef geo::bv::BVH_ST_DG_DOP2_K8 bvh_type;
#  else
    typedef geo::bv::BVH_ST_DG_Sphere3 bvh_type;
#  endif
#endif

public:
    AppTask()
    : m_VizIS(1<<20,1<<15)
    , m_TweakIS(1024,1024)
    , m_pAppView(0)
      // Args
    , m_FileNameSRC(""), m_OutputName( "escalunya_output" )
      // State
    , m_bIsActive(true)
    , m_DirtyFlags( eDF_None )
        {
            // Init libraries
            geo::Init();
#ifdef __ENABLE_CLIP_TEST
            mal::SetRandomSeed(666); //TEMP
#endif
        }
    ~AppTask() {}

    inline void Init( const std::string& src_file_name, const std::string& output_name, const Params& params )
        {
            // Save args
            m_FileNameSRC = src_file_name;
            m_OutputName = output_name;
            if( src_file_name == "" )
            {
                geo::Make_TriSurfaceShape3_Box( m_SRC_Shape, geo::Vec3(2,0.5,6) ); //1,1,1
                m_SRC_Shape.BeginEdition();
                m_SRC_Shape.Transform( geo::Transform3( geo::Vec3::Zero(), mal::GRotation3x3_From( mal::Normalized(geo::Vec3(1,1,1)), geo::Real(0.1) ) ) );
                m_SRC_Shape.EndEdition();
            }
            else
            {
                // Load EPS2
                APP_LOG( "Loading file %s", src_file_name.c_str() );
                OFF_Load_ETSS3( m_SRC_Shape, m_FileNameSRC );
            }
            m_DirtyFlags.Enable( eDF_SRC );
            // Init params
            m_Params = params;
            m_Params.m_DirtyFlags.Enable( eDF_SRC ); //TEMP: params.flags do OVERWRITE flags at OnSync_Params
        }

    void RebuildAllPhasesAndSave()
        {
            Rebuild_SRC();
            Rebuild_VIZ();
            Rebuild_EXT();
            Rebuild_INT();
            Rebuild_SIM();
            Rebuild_CLP();
            Rebuild_DCR();
            Rebuild_BVH();
            Rebuild_EMB();
            OnSaveCall(0);
        }

    void OnSync_Params( sfr::gui::PropertyIt pit )
        {
            s_Escalunya_ArchetypeLibrary.SyncInstance( "Archetype_Escalunya_Params", &m_Params, pit );
            //TEMP Sync DF with params.DF
            m_DirtyFlags = m_Params.m_DirtyFlags;
            m_Params.m_DirtyFlags.Set( eDF_None );
        }

    const util::ItemStream &GetVizIS() const { return m_VizIS; }

    void SetView( sfr::IView *p_view )
        {
            m_pAppView = p_view;

            if( s_Escalunya_ArchetypeLibrary.IsEmpty() )
            {
                Params::InitArchetype( s_Escalunya_ArchetypeLibrary );
            }
            s_ParamsIS.BeginComplex( "AppParams", eType_Property_Group );
            {
                s_Escalunya_ArchetypeLibrary.ExportInstance( "Archetype_Escalunya_Params", &m_Params, s_ParamsIS );
            }
            s_ParamsIS.EndComplex();
            util::ItemStream::ItemItRW sit = s_ParamsIS.BeginRW().GetSubItem();
            sfr::gui::IWidget *pPTW
                = m_pAppView->GetDesktop()->CreateWidget_PropertyTree( "Escalunya",
                                                                       sfr::gui::eDWF_Default,
                                                                       sfr::Vec2(10,10),
                                                                       sit,
                                                                       boost::bind( &AppTask::OnSync_Params, this, _1 ) );
            pPTW->GetBasicAPI()->Minimize();
            pPTW->GetBasicAPI()->Unroll();
            pPTW->GetBasicAPI()->AddFootButton( "(Load)", boost::bind( &AppTask::OnLoadCall, this, _1 ), 0 );
            pPTW->GetBasicAPI()->AddFootButton( "(Save)", boost::bind( &AppTask::AppTask::OnSaveCall, this, _1 ), 0 );
            /*
            pPTW->GetBasicAPI()->AddFootButton( "(Subd)", boost::bind( &AppTask::AppTask::OnSubdCall, this, _1 ), 0 );
            pPTW->GetBasicAPI()->AddFootButton( "(CropSimExt)", boost::bind( &AppTask::AppTask::OnCropSimExtCall, this, _1 ), 0 );
            */
        }
    inline sfr::IView *GetAppView() const { return m_pAppView; }

    bool Update( double dt )
        {
            //APP_LOG( "Updating AppTask" );
            if( !m_bIsActive ) return false;
            else
            {
                bool bRedraw = m_DirtyFlags.Test( eDF_All );
                if( m_DirtyFlags.Test( eDF_SRC ) ) { Rebuild_SRC(); m_DirtyFlags.Enable( eDF_VIZ | eDF_EXT | eDF_INT ); }
                if( m_DirtyFlags.Test( eDF_VIZ ) ) { Rebuild_VIZ(); m_DirtyFlags.Enable( eDF_SIM ); /*m_DirtyFlags.Enable( eDF_CLP); m_DirtyFlags.Enable( eDF_DCR ); m_DirtyFlags.Enable( eDF_EMB );*/ }
                if( m_DirtyFlags.Test( eDF_EXT ) ) { Rebuild_EXT(); m_DirtyFlags.Enable( eDF_SIM ); }
                if( m_DirtyFlags.Test( eDF_INT ) ) { Rebuild_INT(); m_DirtyFlags.Enable( eDF_SIM ); }
                if( m_DirtyFlags.Test( eDF_SIM ) ) { Rebuild_SIM(); m_DirtyFlags.Enable( eDF_CLP); m_DirtyFlags.Enable( eDF_DCR ); m_DirtyFlags.Enable( eDF_EMB ); }
                if( m_DirtyFlags.Test( eDF_CLP ) ) { Rebuild_CLP(); }
                if( m_DirtyFlags.Test( eDF_DCR ) ) { Rebuild_DCR(); m_DirtyFlags.Enable( eDF_BVH ); }
                if( m_DirtyFlags.Test( eDF_BVH ) ) { Rebuild_BVH(); }
                if( m_DirtyFlags.Test( eDF_EMB ) ) { Rebuild_EMB(); }
                if( bRedraw || m_DirtyFlags.Test( eDF_Redraw ) ) RedrawAll();
                m_DirtyFlags.Set( eDF_None );
                return true;
            }
        }

    // SRC( natural_scale )
    inline void Rebuild_SRC()
        {
            APP_LOG( "Rebuild_SRC()" );
            // Apply natural scale
            geo::bv::AABB3 aabb;
            m_SRC_Shape.ComputeBVD( aabb, geo::Transform3::Identity(), 0 );
            Vec3f aabb_sizes( 2.0f * aabb.GetHalfSizes() );
            m_SRC_Shape.BeginEdition();
            m_SRC_Shape.Transform( geo::Transform3( -aabb.GetPos(), geo::Mat3x3::Identity() ) ); //center
            m_SRC_Shape.Scale( m_Params.m_SRC_NaturalSize / mal::Max(aabb_sizes) );
            m_SRC_Shape.EndEdition();
            //Try to close AFTER scaling (tolerances change a lot!)
            if( !m_SRC_Shape.IsClosed() )
            {
                APP_LOG_WARNING("Rebuild_SRC() SRC is open, trying to close holes...");
                unsigned int num_iter(0);
                do
                {
                    APP_LOG("TryToCloseHoles() iter %u", num_iter);
                    m_SRC_Shape.BeginEdition();
                    m_SRC_Shape.TryToCloseHoles(m_Params.m_CLP_EpsilonLength);
                    m_SRC_Shape.EndEdition();
                    num_iter++;
                }
                while( !m_SRC_Shape.IsClosed() && num_iter < m_Params.m_CLP_MaxIter_CloseHoles );
                if( !m_SRC_Shape.IsClosed() ) { APP_LOG_WARNING( "Rebuild_SRC() cannot close, SRC remains OPEN" ); }
                else { APP_LOG( "Rebuild_SRC() successfully CLOSED" ); }
            }
            // Reset path and redraw
            m_SRC_Object.ResetShape( &m_SRC_Shape );
            m_SRC_Object.SetTransform( geo::Transform3::Identity() );
        }

    // VIZ( SRC, viz.criteria_size )
    inline void Rebuild_VIZ()
        {
            APP_LOG( "Rebuild_VIZ()" );
            // Refine or Simplify SRC to the required VIZ detail
            geo::Make_TriSurfaceShape3_From_TriSurfaceShape3_Detail( m_SRC_Shape, m_Params.m_VIZ_Detail, m_VIZ_Shape );
            m_VIZ_Shape.AddBVH(); //\todo ASSUME VIZ == CD
            // BVH stats
            {
                float volume(0);
                float h_volume(0);
                float b_volume(0);
                m_VIZ_Shape.GetBVH()->Map( [&volume,&h_volume,&b_volume]
                                           ( const geo::BVH_TriSurfaceShape3& bvh, const geo::BVH_TriSurfaceShape3::Node& node, uint32 level )
                                           {
                                               float v = geo::bv::ComputeVolume( node.m_Geometry.m_BV );
                                               volume += v;
                                               h_volume += v / (level+1);
                                               b_volume += v / mal::Pow<float>(2,level); //\note 1/2^l seems reasonable, as overlap probability /2 at each level...
                                               return true;
                                           } );
                APP_LOG( "BVH(VIZ) Volume = %f, HierarchyVolume = %f, BinaryVolume = %f", volume, h_volume, b_volume );
            }
            m_VIZ_Object.ResetShape( &m_VIZ_Shape );
            m_VIZ_Object.SetTransform( geo::Transform3( geo::Vec3( -1.1 * m_Params.m_SRC_NaturalSize, 0, 0 ), geo::Mat3x3::Identity() ) );

#ifdef __ENABLE_BVH_TEST
            /*\todo If we do this, then Fit_TetSolidShape3_To_TriSurfaceShape3 fails because BVH(VIZ) and SIM do not match...
            m_VIZ_Shape.GetBVH()->Refit( boost::bind<void>( &geo::GEBV_TriSolidShape3_E<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                            &m_VIZ_Shape, m_VIZ_Object.GetTransform(), m_VIZ_Object.GetVecSDOF(), _1, _2) );
            */
#endif
        }

    // +Offset( Discretize( SRC ) )
    // \todo Avoid intersection with VIZ while minimizing detail and distance to VIZ
    inline void Rebuild_EXT()
        {
            APP_LOG( "Rebuild_EXT() with %f and %f", m_Params.m_EXT_Detail, m_Params.m_EXT_Offset );
#define __ESCALUNYA_OFFSET_AFTER_SIMPLIFY //\todo This seems to work better but causes irregular expansion
#ifdef __ESCALUNYA_OFFSET_AFTER_SIMPLIFY
            if( m_Params.m_EXT_Detail < 1.0f || m_Params.m_EXT_Offset > 0.0f )
            {
                // TMP = Simplify( SRC )
                geo::EditableTriSurfaceShape3 TMP_Shape;
                geo::Make_TriSurfaceShape3_From_TriSurfaceShape3_Simplify( m_SRC_Shape, m_Params.m_EXT_Detail, TMP_Shape );
                APP_LOG( "Rebuild_EXT() EXT simplified with #V = %u, T = %u", TMP_Shape.GetNumV(), TMP_Shape.GetNumT() );
                // EXT = +Offset( TMP )
                geo::Make_TriSurfaceShape3_From_TriSurfaceShape3_Offset( TMP_Shape, m_Params.m_SRC_NaturalSize*m_Params.m_EXT_Offset, m_EXT_Shape );
                APP_LOG( "Rebuild_EXT() EXT offset with %d vertices", m_EXT_Shape.GetNumV() );
            }
            else
                m_EXT_Shape.Set( m_SRC_Shape );
#else
            // TMP = +Offset( SRC )
            geo::EditableTriSurfaceShape3 TMP_Shape;
            geo::Make_TriSurfaceShape3_From_TriSurfaceShape3_Offset( m_SRC_Shape, m_Params.m_SRC_NaturalSize*m_Params.m_EXT_Offset, TMP_Shape );
            APP_LOG( "Rebuild_EXT() EXT offset with %d vertices", TMP_Shape.GetNumV() );
            // EXT = Simplify( TMP )
            geo::Make_TriSurfaceShape3_From_TriSurfaceShape3_Simplify( TMP_Shape, m_Params.m_EXT_Detail, m_EXT_Shape );
            APP_LOG( "Rebuild_EXT() EXT discretized with %d vertices", m_EXT_Shape.GetNumV() );
#endif
            m_EXT_Object.ResetShape( &m_EXT_Shape );
            m_EXT_Object.SetTransform( geo::Transform3( geo::Vec3( 1.1*m_Params.m_SRC_NaturalSize, 0, 0 ), geo::Mat3x3::Identity() ) );
        }

    // Offset( -Discretize( SRC ) )
    // \todo Avoid intersection with VIZ while minimizing detail and distance to VIZ
    inline void Rebuild_INT()
        {
            //APP_LOG( "Rebuild_INT()" );
        }

    // CDT( EXT, INT )
    // \todo Embedding( CD, SIM )
    inline void Rebuild_SIM()
        {
            APP_LOG( "Rebuild_SIM()" );
#define __USE_FAST_TETRAHEDRALIZATION
#ifdef __USE_FAST_TETRAHEDRALIZATION
            geo::bv::AABB3 aabb;
            m_VIZ_Shape.ComputeBVD( aabb, geo::Transform3::Identity(), 0 );
            //aabb.Extend( 0.01 * mal::Max(aabb.GetHalfSizes()) ); //\note Grow 1% to avoid boundary coincidences in Clip_TriSurfaceShape3_TetSolidShape3()
            aabb.Extend( 0.1 * mal::Max(aabb.GetHalfSizes()) ); //\note Grow 1% to avoid boundary coincidences in Clip_TriSurfaceShape3_TetSolidShape3()
            geo::Make_TetSolidShape3_From_AABB3_Minimal( aabb.GetMin(), aabb.GetMax(), m_Params.m_SIM_CDT_Cell_Size, m_SIM_Shape );
            geo::Fit_TetSolidShape3_To_TriSurfaceShape3( m_VIZ_Shape,
                                                         m_Params.m_SIM_Fit_ODT_RelaxationCoeff, m_Params.m_SIM_Fit_Lpc_RelaxationCoeff, m_Params.m_SIM_Fit_MaxIter, m_Params.m_SIM_Fit_FixNMF,
                                                         m_SIM_Shape );
            GEO_LOG_WARNING("SIM stats: #V = %d, #T = %d, #BF = %d, #L = %d", m_SIM_Shape.GetNumV(), m_SIM_Shape.GetNumT(), m_SIM_Shape.GetNumBF(), m_SIM_Shape.GetNumL() );
#else
            geo::Make_TetSolidShape3_From_TriSurfaceShape3_CDT( m_EXT_Shape,
                                                                m_Params.m_SIM_CDT_Facet_Angle, m_Params.m_SIM_CDT_Facet_Size, m_Params.m_SIM_CDT_Facet_Distance,
                                                                m_Params.m_SIM_CDT_Cell_Ratio, m_Params.m_SIM_CDT_Cell_Size,
                                                                m_Params.m_SIM_CDT_Lloyd, m_Params.m_SIM_CDT_Odt, m_Params.m_SIM_CDT_Perturb, m_Params.m_SIM_CDT_Exude,
                                                                m_SIM_Shape );
#endif
            m_SIM_Object.ResetShape( &m_SIM_Shape );
            m_SIM_Object.SetTransform( geo::Transform3( geo::Vec3( 2.2 * m_Params.m_SRC_NaturalSize, 0, 0 ), geo::Mat3x3::Identity() ) );

            // TEMP: Add BVH to SIM... ideally should be used in all overlap tests automatically
            m_SIM_Shape.AddBVH();
            // BVH stats
            {
                float volume(0);
                float h_volume(0);
                float b_volume(0);
                m_SIM_Shape.GetBVH()->Map( [&volume,&h_volume,&b_volume]
                                           ( const geo::BVH_TetSolidShape3& bvh, const geo::BVH_TetSolidShape3::Node& node, uint32 level )
                                           {
                                               float v = geo::bv::ComputeVolume( node.m_Geometry.m_BV );
                                               volume += v;
                                               h_volume += v / (level+1);
                                               b_volume += v / mal::Pow<float>(2,level); //\note 1/2^l seems reasonable, as overlap probability /2 at each level...
                                               return true;
                                           } );
                APP_LOG( "BVH(SIM) Volume = %f, HierarchyVolume = %f, BinaryVolume = %f", volume, h_volume, b_volume );
            }
        }

    inline void Rebuild_BVH()
        {
#ifdef __ENABLE_BVH_TEST
            switch( m_Params.m_BVH_Method )
            {
            case Params::eBVHM_BV_E: m_BVH.Rebuild_TopDown( m_SIM_Shape.GetNumT(),
                                                            boost::bind<void>( &geo::GEBV_TetSolidShape3_E<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                               &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
#  ifdef __TODO_ESCALUNYA
                //\todo!!
            case Params::eBVHM_BV_BSLAB: m_BVH.Rebuild_BottomUp( m_SIM_Shape.GetNumP(),
                                                                 boost::bind<void>( &geo::GEBV_MeshSolidShape2_BSlab_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                                    &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
            case Params::eBVHM_BV_BDOP: m_BVH.Rebuild_BottomUp( m_SIM_Shape.GetNumP(),
                                                                boost::bind<void>( &geo::GEBV_MeshSolidShape2_BDOP_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                                   &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
#  endif
            default: break;
            }
            // BVH stats
            {
                float volume(0);
                float h_volume(0);
                float b_volume(0);
                m_BVH.Map( [&volume,&h_volume,&b_volume]
                           ( const bvh_type& bvh, const bvh_type::Node& node, uint32 level )
                           {
                               float v = geo::bv::ComputeVolume( node.m_Geometry.m_BV );
                               volume += v;
                               h_volume += v / (level+1);
                               b_volume += v / mal::Pow<float>(2,level); //\note 1/2^l seems reasonable, as overlap probability /2 at each level...
                               return true;
                           }
                    );
                APP_LOG( "BVH(SIM) Volume = %f, HierarchyVolume = %f, BinaryVolume = %f", volume, h_volume, b_volume );
            }
#endif
        }

    // Clip( VIZ, SIM )
    inline void Rebuild_CLP()
        {
            APP_LOG( "Rebuild_CLP()" );

            g_TTS_SliverThreshold = mal::Exp10( m_Params.m_CLP_SliverThreshold_Log10 );

#define __USE_BVH_FOR_CLIP
#ifdef __USE_BVH_FOR_CLIP
            //\todo We rebuild BVH at origin to match VIZ SHAPE, not Object... should do this differently to avoid moving stuff to match shapes...
            m_SIM_Shape.GetBVH()->Rebuild_TopDown( m_SIM_Shape.GetNumT(),
                                                   boost::bind<void>( &geo::GEBV_TetSolidShape3_E<geo::BVH_TetSolidShape3::entry_index_type,geo::BVH_TetSolidShape3::bv_type>,
                                                   &m_SIM_Shape, geo::Transform3::Identity(), m_SIM_Object.GetVecSDOF(), _1, _2) );
            geo::Clip_TriSurfaceShape3_TetSolidShape3( m_VIZ_Shape, m_SIM_Shape, m_SIM_Shape.GetBVH(), m_CLP_Shape );
#else
            geo::Clip_TriSurfaceShape3_TetSolidShape3( m_VIZ_Shape, m_SIM_Shape, 0, m_CLP_Shape );
#endif
            APP_LOG_WARNING("Rebuild_CLP() VIZ.IsClosed() = %s with  #V = %d, #T = %d => CLP.IsClosed() = %s with #V = %d, #T = %d",
                            m_VIZ_Shape.IsClosed()?"true":"false",
                            m_VIZ_Shape.GetNumV(), m_VIZ_Shape.GetNumT(),
                            m_CLP_Shape.IsClosed()?"true":"false",
                            m_CLP_Shape.GetNumV(), m_CLP_Shape.GetNumT() );
            if( m_VIZ_Shape.IsClosed() && !m_CLP_Shape.IsClosed() )
            {
                APP_LOG_WARNING("Rebuild_CLP() Trying to close holes in CLP that did not exist in VIZ...");
                unsigned int num_iter(0);
                do
                {
                    //\todo With 0, collapse only shortest, CLP mesh
                    //is minimally changed, however, it takes a LOT of
                    //iterations to fix complex meshes (some of them
                    //are NOT fixed in 20 iter (armadillo_HD?). Using
                    //default epsilon (infinity) collapses ALL open
                    //edges in a single pass, which fixes models very
                    //fast, BUT distorts CLP and causes non-manifold
                    //meshes... consider alternatives or a more
                    //adaptive
                    //behaviour, such as collapse-shortest-and-all-below-epsilon
                    APP_LOG("TryToCloseHoles() iter %u", num_iter);
                    m_CLP_Shape.BeginEdition();
                    m_CLP_Shape.TryToCloseHoles(m_Params.m_CLP_EpsilonLength);
                    m_CLP_Shape.EndEdition();
                    num_iter++;
                }
                while( !m_CLP_Shape.IsClosed() && num_iter < m_Params.m_CLP_MaxIter_CloseHoles );
                APP_LOG_WARNING("Rebuild_CLP() VIZ.IsClosed() = %s with  #V = %d, #T = %d => CLP.IsClosed() = %s with #V = %d, #T = %d",
                                m_VIZ_Shape.IsClosed()?"true":"false",
                                m_VIZ_Shape.GetNumV(), m_VIZ_Shape.GetNumT(),
                                m_CLP_Shape.IsClosed()?"true":"false",
                                m_CLP_Shape.GetNumV(), m_CLP_Shape.GetNumT() );
            }
            m_CLP_Object.ResetShape( &m_CLP_Shape );
            m_CLP_Object.SetTransform( geo::Transform3( geo::Vec3( 3.3 * m_Params.m_SRC_NaturalSize, 0, 0 ), geo::Mat3x3::Identity() ) );
        }

    // Embedding( VIZ, SIM )
    inline void Rebuild_EMB()
        {
            APP_LOG( "Rebuild_EMB()" );
            // Set EMB to VIZ (\todo Consider CLP instead)
            switch( m_Params.m_EMB_Source )
            {
            case Params::eEMBS_VIZ: m_EMB_Shape.Set( m_VIZ_Shape ); break;
            case Params::eEMBS_CLP: m_EMB_Shape.Set( m_CLP_Shape ); break;
            default: break;
            }
            // Reset embedding
            m_EMB_Object.Embed(0);
            m_EMB_Object.ResetShape( &m_EMB_Shape );
            m_EMB_Object.SetTransform( geo::Transform3::Identity() );
            // Recompute embedding (\todo with Tr = Identity, by now)
            geo::Transform3 deformer_tr = m_SIM_Object.GetTransform();
            m_SIM_Object.SetTransform( geo::Transform3::Identity() );
            m_EMB_Object.Embed( &m_SIM_Object, m_Params.m_EMB_Method );
            m_SIM_Object.SetTransform( deformer_tr );
        }

    inline void Rebuild_DCR()
        {
            APP_LOG( "Rebuild_DCR()" );
            // Create DCR
            switch( m_Params.m_DCR_Source )
            {
            case Params::eDCRS_VIZ: m_SIM_Shape.AddDCR( &m_VIZ_Shape, geo::Transform3::Identity() ); break; //\note tr is identity because tr_s2tss is relative and they are both in globals here
            case Params::eDCRS_CLP: m_SIM_Shape.AddDCR( &m_CLP_Shape, geo::Transform3::Identity() ); break; //\note tr is identity because tr_s2tss is relative and they are both in globals here
            default: break;
            }
            // Log
            const geo::DCR_TetSolidShape3* pDCR( m_SIM_Shape.GetDCR() );
            if( pDCR )
            {
                APP_LOG( "DCR(SIM,VIZ) with %u elements, %u patches, %u triangles and %u vertices",
                         pDCR->m_NumElements, pDCR->m_NumPatches, pDCR->m_NumTriangles, pDCR->m_NumVertices );
            }
            /*TEMP: Verbose!
            APP_LOG( "----------Elements----------" );
            for( unsigned int it_element=0; it_element < m_pDCR->m_NumElements; it_element++ )
            {
                const geo::DCR_MeshSolidShape2::ElementData& ed( m_pDCR->m_vecED[it_element] );
                APP_LOG( "Element[%d]: FirstGPID %d, NumGP %d, FirstVID %d, NumV %d", it_element, ed.m_FirstGPID, ed.m_NumPatches, ed.m_FirstVID, ed.m_NumVertices );
            }
            APP_LOG( "----------Patches----------" );
            for( unsigned int it_patch=0; it_patch < m_pDCR->m_NumPatches; it_patch++ )
            {
                const geo::DCR_MeshSolidShape2::PatchData& pd( m_pDCR->m_vecPD[it_patch] );
                APP_LOG( "Patch[%d]: FirstSID %d, NumS %d", it_patch, pd.m_FirstSID, pd.m_NumSegments );
            }
            APP_LOG( "----------Segments----------" );
            for( unsigned int it_segment=0; it_segment < m_pDCR->m_NumSegments; it_segment++ )
            {
                const geo::DCR_MeshSolidShape2::patch_segment_topology& sd( m_pDCR->m_vecS[it_segment] );
                APP_LOG( "Segment[%d]: VID = %d,%d, NSID = %d,%d", it_segment, sd.m_vecVID[0], sd.m_vecVID[1], sd.m_vecNeighbourSID[0], sd.m_vecNeighbourSID[1] );
            }
            */
        }

    void RedrawAll()
        {
            m_VizIS.Clear();
#ifdef __ENABLE_CLIP_TEST
            {
                geo::Vec3 vec_polygon_clipped_pos[6];

                /* This params generate 7 cp that are reduced to 5, keep them as a coincicences test
                geo::Real height(0.5);
                geo::Real length(0.5);
                geo::Vec3 offset_tri(0.25,0,0.25);
                geo::Vec3 vec_tri_pos[3] = { offset_tri + geo::Vec3(length,height,length),
                                             offset_tri + geo::Vec3(0,height,-length),
                                             offset_tri + geo::Vec3(-length,height,0) };
                geo::Vec3 vec_tet_pos[4] = { geo::Vec3(0,0,0), geo::Vec3(1,0,0), geo::Vec3(0,1,0), geo::Vec3(0,0,1) };
                */
                geo::Real height(0.5);
                geo::Real length(0.5);
                geo::Vec3 offset_tri(0.25,0,0.25);
                geo::Vec3 vec_tri_pos[3] = { offset_tri + geo::Vec3(length,height,length) + mal::RandomV(geo::Vec3(-length,-height,-length),geo::Vec3(length,height,length)),
                                             offset_tri + geo::Vec3(0,height,-length) + mal::RandomV(geo::Vec3(-length,-height,-length),geo::Vec3(length,height,length)),
                                             offset_tri + geo::Vec3(-length,height,0) + mal::RandomV(geo::Vec3(-length,-height,-length),geo::Vec3(length,height,length)) };
                geo::Vec3 vec_tet_pos[4] = { geo::Vec3(0,0,0), geo::Vec3(1,0,0), geo::Vec3(0,1,0), geo::Vec3(0,0,1) };
                unsigned int num_clipped = geo::np::Clip_Triangle3_Tetrahedron3( vec_tri_pos[0], vec_tri_pos[1], vec_tri_pos[2],
                                                                                 vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3],
                                                                                 vec_polygon_clipped_pos );
                APP_LOG("#Clipped = %u",num_clipped);
                for( unsigned int it_v=0; it_v<num_clipped; it_v++ )
                {
                    VIZ_SEGMENT3( m_VizIS, vec_polygon_clipped_pos[it_v], vec_polygon_clipped_pos[(it_v+1)%num_clipped], 2, Vec4f(0,0,1,1) );
                    char str[32];
                    sprintf(str,"%u",it_v);
                    VIZ_POINT3_NAMED( m_VizIS, str, vec_polygon_clipped_pos[it_v], 2, Vec4f(1,1,1,1) );
                }
                VIZ_TRIANGLE3( m_VizIS, vec_tri_pos[0], vec_tri_pos[1], vec_tri_pos[2], Vec4f(1,0,0,1), util::eVizStyle_Wire );
                VIZ_TETRAHEDRON3( m_VizIS, vec_tet_pos[0], vec_tet_pos[1], vec_tet_pos[2], vec_tet_pos[3], Vec4f(0,1,0,1), util::eVizStyle_Wire );
            }
#else
            geo::VizObject( &m_SRC_Object, m_VizIS, m_Params.m_SRC_DDF );
            geo::VizObject( &m_VIZ_Object, m_VizIS, m_Params.m_VIZ_DDF );
            geo::VizObject( &m_EXT_Object, m_VizIS, m_Params.m_EXT_DDF );
            geo::VizObject( &m_SIM_Object, m_VizIS, m_Params.m_SIM_DDF );
            geo::VizObject( &m_CLP_Object, m_VizIS, m_Params.m_CLP_DDF );
            Redraw_DCR();
            geo::VizObject( &m_EMB_Object, m_VizIS, m_Params.m_EMB_DDF );
            Redraw_EMB(); //\todo BEFORE RedrawBVH() because it REFITS IT!!

            //Labels
            {
                char str[128];
                snprintf(str,128,"SRC(%u,T%u)", m_SRC_Shape.GetNumV(), m_SRC_Shape.GetNumT() );
                VIZ_POINT3_NAMED(m_VizIS,str,m_SRC_Object.GetTransform().Pos()+Vec3f(0,m_Params.m_SRC_NaturalSize,0),5,Vec4f(1,1,1,1));
                snprintf(str,128,"VIZ(%u,%u)", m_VIZ_Shape.GetNumV(), m_VIZ_Shape.GetNumT() );
                VIZ_POINT3_NAMED(m_VizIS,str,m_VIZ_Object.GetTransform().Pos()+Vec3f(0,m_Params.m_SRC_NaturalSize,0),5,Vec4f(1,1,1,1));
                snprintf(str,128,"EXT(%u,%u)", m_EXT_Shape.GetNumV(), m_EXT_Shape.GetNumT() );
                VIZ_POINT3_NAMED(m_VizIS,str,m_EXT_Object.GetTransform().Pos()+Vec3f(0,m_Params.m_SRC_NaturalSize,0),5,Vec4f(1,1,1,1));
                snprintf(str,128,"SIM(%u,%u)", m_SIM_Shape.GetNumV(), m_SIM_Shape.GetNumT() );
                VIZ_POINT3_NAMED(m_VizIS,str,m_SIM_Object.GetTransform().Pos()+Vec3f(0,m_Params.m_SRC_NaturalSize,0),5,Vec4f(1,1,1,1));
                snprintf(str,128,"DCR(%u,%u,%u,%u)", m_SIM_Shape.GetDCR()->m_NumElements, m_SIM_Shape.GetDCR()->m_NumPatches, m_SIM_Shape.GetDCR()->m_NumVertices, m_SIM_Shape.GetDCR()->m_NumTriangles );
                VIZ_POINT3_NAMED(m_VizIS,str,m_SIM_Object.GetTransform().Pos()+Vec3f(0,m_Params.m_SRC_NaturalSize+0.25,0),5,Vec4f(1,1,1,1));
                snprintf(str,128,"CLP(%u,%u)", m_CLP_Shape.GetNumV(), m_CLP_Shape.GetNumT() );
                VIZ_POINT3_NAMED(m_VizIS,str,m_CLP_Object.GetTransform().Pos()+Vec3f(0,m_Params.m_SRC_NaturalSize,0),5,Vec4f(1,1,1,1));
            }

#endif
#ifdef __ENABLE_BVH_TEST
            geo::GVizBVH( m_BVH, m_VizIS, geo::eBVHDDF_Default | geo::eBVHDDF_Levels );
#endif
        }

    void Redraw_EMB()
        {
            // Distort Deformer, redraw embedded (automatically synced) and restore Deformer afterwards
            if( 0 != m_SIM_Object.GetShape() && 0 != m_EMB_Object.GetShape() )
            {
                // Save deformer SDOF and transform
                uint32 num_deformer_sdof = m_SIM_Object.GetShape()->GetNumSDOF();
                geo::Vec3* vec_saved_deformer_sdof = new geo::Vec3[ num_deformer_sdof ];
                memcpy( vec_saved_deformer_sdof, m_SIM_Object.GetVecSDOF(), sizeof(geo::Vec3)*num_deformer_sdof );
                geo::Transform3 deformer_tr = m_SIM_Object.GetTransform();
                // Distort deformer SDOF
                mal::SetRandomSeed(666);
                for( uint32 it_sdof=0; it_sdof < num_deformer_sdof; it_sdof++ )
                    m_SIM_Object.GetVecSDOF()[it_sdof] += m_Params.m_SRC_NaturalSize * m_Params.m_EMB_DistortionScale * mal::RandomUnitVec<geo::Real,3>();
                m_SIM_Object.SetTransform( geo::Transform3( geo::Vec3( 4.4 * m_Params.m_SRC_NaturalSize, 0, 0 ), geo::Mat3x3::Identity() ) );
                // Draw distorted m_SIM_Object and restore SDOF afterwards
                geo::VizObject( &m_SIM_Object, m_VizIS, geo::eODDF_Shape | geo::eSDDF_Interior );
                geo::VizObject( &m_EMB_Object, m_VizIS, geo::eODDF_Shape | geo::eSDDF_Boundary );
// #ifdef __ENABLE_BVH_TEST
//                 switch( m_Params.m_BVH_Method )
//                 {
//                 case Params::eBVHM_BV_E: m_BVH.Refit( boost::bind<void>( &geo::GEBV_MeshSolidShape2_E<bvh_type::entry_index_type,bvh_type::bv_type>,
//                                                                          &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
//                 case Params::eBVHM_BV_BSLAB: m_BVH.Refit( boost::bind<void>( &geo::GEBV_MeshSolidShape2_BSlab_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
//                                                                              &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
//                 case Params::eBVHM_BV_BDOP: m_BVH.Refit( boost::bind<void>( &geo::GEBV_MeshSolidShape2_BDOP_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
//                                                                             &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
//                 default: break;
//                 }
// #endif
                // Restore undistorted SDOF
                memcpy( m_SIM_Object.GetVecSDOF(), vec_saved_deformer_sdof, sizeof(geo::Vec3)*num_deformer_sdof );
                m_SIM_Object.SetTransform( deformer_tr );
                //\note m_EMB_Object transform is OUT OF SYNC here, but should be recomputed when accessed
                //m_EMB_Object.SyncEmbedding();
                // Dealloc
                delete[] vec_saved_deformer_sdof;
            }
        }

    void Redraw_DCR()
        {
            /*TEMP: This is done automatically in VizObject( m_SIM_Object )
            if( m_pDCR == 0 ) return;
            geo::VizDCR( m_pDCR,
                         &m_SIM_Shape,
                         geo::Transform2( geo::Vec2( 2.2 * m_Params.m_SRC_NaturalSize, 0 ),
                                          geo::Mat2x2::Identity() ),
                         m_SIM_Object.GetVecDOF(),
                         m_VizIS, 0xFFFFFFFF );
            */
        }

    void OnSaveCall( void *p_user_data )
        {
            //\todo TRY to save Params ItemIt to the output file, maybe as a comment...
            APP_LOG( "Saving %s + { .txt, .bin, .params }", m_OutputName.c_str() );
            // Create SL
            geo::ShapeLibrary sl;
            sl.Reserve( 1<<20, util::ItemStream::eRealloc_Data | util::ItemStream::eRealloc_Identifiers ); //Allow full realloc, no risk of shape pointer invalidation
            sl.Register( m_SRC_Shape, "SRC" );
            sl.Register( m_VIZ_Shape, "VIZ" );
            sl.Register( m_SIM_Shape, "SIM" );
            // Register embedded shape //\todo This IGNORES EMBEDDING, as it's in the Object, not Shape, but saves the detailed/refined TriSurface shape
            sl.Register( m_EMB_Shape, "EMB" );
            // Save SL
            sl.Save( (m_OutputName + ".txt").c_str(), false );
            sl.Save( (m_OutputName + ".bin").c_str(), true ); //\todo redundant, keep while testing txt/bin formats
            // Save Params
            {
                //\todo InitArchetype may be neceessary if called from command-line, where no SetView() happens
                if( s_Escalunya_ArchetypeLibrary.IsEmpty() )
                {
                    Params::InitArchetype( s_Escalunya_ArchetypeLibrary );
                }
                util::ItemStream params_is;
                params_is.BeginComplex( "AppParams", eType_Property_Group );
                {
                    s_Escalunya_ArchetypeLibrary.ExportInstance( "Archetype_Escalunya_Params", &m_Params, params_is );
                }
                params_is.EndComplex();
                params_is.SaveTxt( (m_OutputName + ".params").c_str() );
            }
        }

    void OnLoadCall( void *p_user_data )
        {
#ifdef __TODO_ESCALUNYA
            geo::ShapeLibrary sl;
            sl.Reserve( 1<<20, util::ItemStream::eRealloc_Data | util::ItemStream::eRealloc_Identifiers ); //Allow full realloc, no risk of shape pointer invalidation

            sl.Load( (m_OutputName + ".txt").c_str(), false );
            //sl.Load( "shape_library.bin", true );
            /* testing bin save
            sl.Save( "shape_library.save.bin", true );
            sl.Save( "shape_library.save.txt", false );
            */
            geo::ShapeFactory sf( sl );

            // Load EPS
            geo::IShape *pShape;
            geo::ShapeID sid1 = sl.GetShapeIdByName("#1");
            APP_LOG( "Creating shape '#1' with sid %d", sid1 );
            pShape = sf.CreateExclusive( sid1 );
            APP_ASSERT( pShape->GetType() == geo::eShape_Path2 );
            m_SRC_Shape.Set( *static_cast<const geo::PathShape2*>( pShape ) );
            m_SRC_Object.ResetShape( &m_SRC_Shape );

            // Load MSS
            //pShape = sf.CreateExclusive( 1 );
            geo::ShapeID sid2 = sl.GetShapeIdByName("#2");
            APP_LOG( "Creating shape '#2' with sid %d", sid2 );
            pShape = sf.CreateExclusive( sid2 );
            APP_ASSERT( pShape->GetType() == geo::eShape_MeshSolid2 );
            m_SIM_Shape.Set( *static_cast<const geo::MeshSolidShape2*>( pShape ) );
            m_SIM_Object.ResetShape( &m_SIM_Shape );

            m_Params.m_DirtyFlags.Enable( eDF_EMB );
            m_DirtyFlags.Enable( eDF_EMB );
#endif
        }

    void OnSubdCall( void *p_user_data )
        {
#ifdef __TODO_ESCALUNYA
            m_SIM_Shape.Subdivide();
            m_SIM_Object.ResetShape( &m_SIM_Shape );
            m_Params.m_DirtyFlags.Enable( eDF_EMB );
            m_DirtyFlags.Enable( eDF_EMB );
#endif
        }

    void OnCropSimExtCall( void *p_user_data )
        {
#ifdef __TODO_ESCALUNYA
#endif //__TODO_ESCALUNYA
        }

#ifdef __ENABLE_BVH_TEST
    bool TestBVH( const Vec2f& p ) const
    {
#  ifdef __ENABLE_BVH_TEST_TODO
        bvh_type::bv_type bv( p );
        std::vector< bvh_type::entry_index_type > vec_entry;
        bool bHit = m_BVH.Test( bv, vec_entry );
        if( bHit ) APP_LOG("TestBVH hit!");
        return bHit;
#  else
        return false;
#  endif
    }
#endif

private:
    util::ItemStream m_VizIS;
    util::ItemStream m_TweakIS;
    sfr::IView *m_pAppView;

    // Args
    std::string m_FileNameSRC;
    std::string m_OutputName;

    // State
    bool m_bIsActive;
    Flags32 m_DirtyFlags;

    // SRC
    geo::EditableTriSurfaceShape3 m_SRC_Shape;
    geo::GObjectSS<geo::EditableTriSurfaceShape3> m_SRC_Object;

    // VIZ( SRC )
    geo::EditableTriSurfaceShape3 m_VIZ_Shape;
    geo::GObjectSS<geo::EditableTriSurfaceShape3> m_VIZ_Object;

    // EXT( SRC )
    geo::EditableTriSurfaceShape3 m_EXT_Shape;
    geo::GObjectSS<geo::EditableTriSurfaceShape3> m_EXT_Object;

    // SIM( EXT, INT, \todo VIZ => CD )
    geo::EditableTetSolidShape3 m_SIM_Shape;
    geo::GObjectSS<geo::EditableTetSolidShape3> m_SIM_Object;

    // CLP( VIZ, SIM )
    geo::EditableTriSurfaceShape3 m_CLP_Shape;
    geo::GObjectSS<geo::EditableTriSurfaceShape3> m_CLP_Object;

    // EMB( VIZ, SIM )
    geo::EditableTriSurfaceShape3 m_EMB_Shape;
    geo::GObjectSS<geo::EditableTriSurfaceShape3> m_EMB_Object;

#ifdef __ENABLE_BVH_TEST
    // BVH( CD )
    bvh_type m_BVH;
#endif

};

#endif //ESCALUNYA_APP_TASK_H
