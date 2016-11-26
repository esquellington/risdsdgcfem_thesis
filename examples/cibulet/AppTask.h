#ifndef CIBULET_APP_TASK_H
#define CIBULET_APP_TASK_H

#include "Params.h"
#include <Safra/Safra.h>
#include <Safra/gui/CameraControllerML.h>
#include <Safra/gfx/VizRenderer.h>

#include <boost/bind.hpp>

#include <iostream>
#include <fstream>
#include <string>

#include "SVG_utils.h"
#include "CGAL_utils.h"

#include <Geo/IObject.h>
#include <Geo/shape/shape.h>
#include <Geo/util/Viz.h>
#include <Geo/shape/ShapeLibrary.h>
#include <Geo/shape/ShapeFactory.h>

#include <util/Archetype.h>

#define __ENABLE_EMB_SubdSIM //TEMP, testing, should be default

#define __ENABLE_CROP_SIM_EXT
#ifdef __ENABLE_CROP_SIM_EXT
#include <Geo/np/Overlap.h> //For Overlap_Point2_Polygonal2()
#include <Geo/np/Intersection.h> //For Intersection_Segment2_Segment2()
#endif

//#define __ENABLE_BVH_TEST //TEMP: disabled while developing __ENABLE_EMB_SubdSIM
#ifdef __ENABLE_BVH_TEST
//#  define __USE_BVH_AABB
#  define __USE_BVH_KDOP
//#  define __USE_BVH_SPHERE
#  include <Geo/bv/GBoundingVolumeHierarchy.h> //TEMPORAL!!
#endif

class AppTask: public sfr::IUpdateTarget
{
public:
    Params m_Params;

#ifdef __ENABLE_BVH_TEST
#  ifdef __USE_BVH_AABB
    typedef geo::bv::BVH_ST_DG_AABB2 bvh_type;
#  elif defined(__USE_BVH_KDOP)
    typedef geo::bv::BVH_ST_DG_DOP2_K8 bvh_type;
#  else
    typedef geo::bv::BVH_ST_DG_Sphere2 bvh_type;
#  endif
#endif

public:
    AppTask()
    : m_VizIS(1<<20,1<<15)
    , m_TweakIS(1024,1024)
    , m_pAppView(0)
      // Args
    , m_FileNameSVG( "s2s/shapes2d.svg" ), m_PathElementNameSVG( "Patito" ), m_OutputName( "cibulet_output" )
      // State
    , m_bIsActive(true)
    , m_DirtyFlags( eDF_None )
    , m_pDCR(0)
        {
        }
    ~AppTask() {}

    inline void Init( const char *svg_file_name, const char *svg_path_element_name, const char *output_name, const Params& params )
        {
            // Save args
            m_FileNameSVG = svg_file_name;
            m_PathElementNameSVG = svg_path_element_name;
            m_OutputName = output_name;
            APP_LOG( "Loading file %s, path element %s", svg_file_name, svg_path_element_name );
            // Load EPS2
            SVG_Load_EPS2( m_SRC_Shape, m_FileNameSVG, m_PathElementNameSVG );
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
            Rebuild_DCR();
            Rebuild_BVH();
            Rebuild_EMB();
            OnSaveCall(0);
        }

    void OnSync_Params( sfr::gui::PropertyIt pit )
        {
            s_Cibulet_ArchetypeLibrary.SyncInstance( "Archetype_Cibulet_Params", &m_Params, pit );
            //TEMP Sync DF with params.DF
            m_DirtyFlags = m_Params.m_DirtyFlags;
            m_Params.m_DirtyFlags.Set( eDF_None );
        }

    const util::ItemStream &GetVizIS() const { return m_VizIS; }

    void SetView( sfr::IView *p_view )
        {
            m_pAppView = p_view;

            if( s_Cibulet_ArchetypeLibrary.IsEmpty() )
            {
                Params::InitArchetype( s_Cibulet_ArchetypeLibrary );
            }
            s_ParamsIS.BeginComplex( "AppParams", eType_Property_Group );
            {
                s_Cibulet_ArchetypeLibrary.ExportInstance( "Archetype_Cibulet_Params", &m_Params, s_ParamsIS );
            }
            s_ParamsIS.EndComplex();
            util::ItemStream::ItemItRW sit = s_ParamsIS.BeginRW().GetSubItem();
            sfr::gui::IWidget *pPTW
                = m_pAppView->GetDesktop()->CreateWidget_PropertyTree( "Cibulet",
                                                                       sfr::gui::eDWF_Default,
                                                                       sfr::Vec2(10,10),
                                                                       sit,
                                                                       boost::bind( &AppTask::OnSync_Params, this, _1 ) );
            pPTW->GetBasicAPI()->Minimize();
            pPTW->GetBasicAPI()->Unroll();
            pPTW->GetBasicAPI()->AddFootButton( "(Load)", boost::bind( &AppTask::OnLoadCall, this, _1 ), 0 );
            pPTW->GetBasicAPI()->AddFootButton( "(Save)", boost::bind( &AppTask::AppTask::OnSaveCall, this, _1 ), 0 );
            pPTW->GetBasicAPI()->AddFootButton( "(Subd)", boost::bind( &AppTask::AppTask::OnSubdCall, this, _1 ), 0 );
            pPTW->GetBasicAPI()->AddFootButton( "(CropSimExt)", boost::bind( &AppTask::AppTask::OnCropSimExtCall, this, _1 ), 0 );
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
                if( m_DirtyFlags.Test( eDF_VIZ ) ) { Rebuild_VIZ(); m_DirtyFlags.Enable( eDF_DCR ); m_DirtyFlags.Enable( eDF_EMB ); }
                if( m_DirtyFlags.Test( eDF_EXT ) ) { Rebuild_EXT(); m_DirtyFlags.Enable( eDF_SIM ); }
                if( m_DirtyFlags.Test( eDF_INT ) ) { Rebuild_INT(); m_DirtyFlags.Enable( eDF_SIM ); }
                if( m_DirtyFlags.Test( eDF_SIM ) ) { Rebuild_SIM(); m_DirtyFlags.Enable( eDF_DCR ); m_DirtyFlags.Enable( eDF_EMB ); }
                if( m_DirtyFlags.Test( eDF_DCR ) ) { Rebuild_DCR(); m_DirtyFlags.Enable( eDF_BVH ); }
                if( m_DirtyFlags.Test( eDF_BVH ) ) { Rebuild_BVH(); }
                if( m_DirtyFlags.Test( eDF_EMB ) ) { Rebuild_EMB(); }
                if( bRedraw ) RedrawAll();
                m_DirtyFlags.Set( eDF_None );
                return true;
            }
        }

    // SRC( natural_scale )
    inline void Rebuild_SRC()
        {
            APP_LOG( "Rebuild_SRC()" );
            // Force Closed
            if( !m_SRC_Shape.IsClosed() )
            {
                APP_LOG_WARNING( "Rebuild_SRC() closing open shape" );
                m_SRC_Shape.BeginEdition();
                m_SRC_Shape.Close();
                m_SRC_Shape.EndEdition();
                m_SRC_Object.ResetShape( &m_SRC_Shape );
            }
            // Apply natural scale
            geo::bv::AABB2 aabb;
            m_SRC_Shape.ComputeBVD( aabb, geo::Transform2::Identity(), 0 );
            Vec2f aabb_sizes( 2.0f * aabb.GetHalfSizes() );
            m_SRC_Shape.BeginEdition();
            //m_SRC_Shape.Transform( geo::Transform2( -aabb.GetPos(), geo::Mat2x2::Identity() ) ); //center?
            m_SRC_Shape.Scale( m_Params.m_SRC_NaturalSize / mal::Max(aabb_sizes) );
            m_SRC_Shape.EndEdition();
            // Reset path and redraw
            m_SRC_Object.ResetShape( &m_SRC_Shape );
            m_SRC_Object.SetTransform( geo::Transform2::Identity() );
        }

    // VIZ( SRC, viz.criteria_size )
    inline void Rebuild_VIZ()
        {
            APP_LOG( "Rebuild_VIZ()" );
            // Discretize SRC to the required VIZ detail
            geo::Make_PolygonalShape2_From_PathShape2( m_SRC_Shape, m_Params.m_VIZ_Detail, true, m_VIZ_Shape );
            m_VIZ_Shape.Distort( m_Params.m_VIZ_DistortionScale * m_Params.m_SRC_NaturalSize, m_Params.m_VIZ_DistortionFreq * mal::Pi<geo::Real>() / m_Params.m_SRC_NaturalSize );
            m_VIZ_Object.ResetShape( &m_VIZ_Shape );
            m_VIZ_Object.SetTransform( geo::Transform2( geo::Vec2( -m_Params.m_SRC_NaturalSize, 0 ), geo::Mat2x2::Identity() ) );
        }

    // +Offset( Discretize( SRC ) )
    // \todo Avoid intersection with VIZ while minimizing detail and distance to VIZ
    inline void Rebuild_EXT()
        {
            APP_LOG( "Rebuild_EXT() with %f and %f", m_Params.m_EXT_Detail, m_Params.m_EXT_Offset );
            if( !m_EXT_Shape.IsAlloc() ) m_EXT_Shape.Alloc( 10000, geo::eVDL_Points, geo::eVDL_Points ); //TEMPORAL
            // TMP = Discretize( SRC )
            geo::EditablePolygonalShape2 discretized_polygonal;
            discretized_polygonal.Alloc( 10000, geo::eVDL_Points, geo::eVDL_Points ); //TEMPORAL
            geo::Make_PolygonalShape2_From_PathShape2( m_SRC_Shape, m_Params.m_EXT_Detail, true, discretized_polygonal );
            APP_LOG( "Rebuild_EXT() EXT discretized with %d vertices", discretized_polygonal.GetNumVertices() );
            // EXT = +Offset( TMP )
            geo::Make_PolygonShape2_From_PolygonShape2_Offset( discretized_polygonal, m_Params.m_SRC_NaturalSize*m_Params.m_EXT_Offset, m_EXT_Shape );
            APP_LOG( "Rebuild_EXT() EXT offset with %d vertices", m_EXT_Shape.GetNumVertices() );
            m_EXT_Object.ResetShape( &m_EXT_Shape );
            m_EXT_Object.SetTransform( geo::Transform2( geo::Vec2( m_Params.m_SRC_NaturalSize, 0 ), geo::Mat2x2::Identity() ) );
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
            geo::Make_MeshSolidShape2_From_PolygonalShape2_CDT( m_EXT_Shape, m_Params.m_SIM_CDT_Ratio, m_Params.m_SIM_CDT_Size, m_SIM_Shape );
            m_SIM_Object.ResetShape( &m_SIM_Shape );
            m_SIM_Object.SetTransform( geo::Transform2( geo::Vec2( 2 * m_Params.m_SRC_NaturalSize, 0 ), geo::Mat2x2::Identity() ) );
        }

    inline void Rebuild_BVH()
        {
#ifdef __ENABLE_BVH_TEST
            APP_LOG( "Rebuild_BVH()" );
            switch( m_Params.m_BVH_Method )
            {
            case Params::eBVHM_BV_E: m_BVH.Rebuild_TopDown( m_SIM_Shape.GetNumP(),
                                                            boost::bind<void>( &geo::GEBV_MeshSolidShape2_E<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                               &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
            case Params::eBVHM_BV_BSLAB: m_BVH.Rebuild_TopDown( m_SIM_Shape.GetNumP(),
                                                                boost::bind<void>( &geo::GEBV_MeshSolidShape2_BSlab_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                                   &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
            case Params::eBVHM_BV_BDOP: m_BVH.Rebuild_TopDown( m_SIM_Shape.GetNumP(),
                                                               boost::bind<void>( &geo::GEBV_MeshSolidShape2_BDOP_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                                  &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
            default: break;
            }
#  if __cplusplus > 199711L //C++11 FTW!
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
            APP_LOG( "BVH Volume = %f, HierarchyVolume = %f, BinaryVolume = %f", volume, h_volume, b_volume );
#  endif
#endif
        }

    // Embedding( VIZ, SIM )
    inline void Rebuild_EMB()
        {
            APP_LOG( "Rebuild_EMB()" );

            // (5) Embed PathShape2 into MeshSolidShape2
            // Compute refined polygonal to be embedded
            // \todo PolygonalShape2 is archaic, and does not grow automatically, requires explicit Alloc(), consider modernizing it or using embedded PathShape2 instead
            geo::Make_PolygonalShape2_From_PathShape2( m_SRC_Shape, m_Params.m_VIZ_Detail, true, m_EMB_Shape );
            // Reset embedded polygonal
            m_EMB_Object.Embed(0);
            m_EMB_Object.ResetShape( &m_EMB_Shape );
            m_EMB_Object.SetTransform( geo::Transform2::Identity() );
#ifdef __ENABLE_EMB_SubdSIM
            // Recompute embedding (\todo with Tr = Identity, by now)
            // with Subdivided SIM, but keeping its original nodes unchanged
            m_EMB_SubdSIM_Shape.Set( m_SIM_Shape );
            for( uint32 it_subd=0; it_subd<m_Params.m_EMB_Subd_Iter; it_subd++ )
                m_EMB_SubdSIM_Shape.Subdivide( m_SIM_Shape.GetNumV() );
            m_EMB_SubdSIM_Object.ResetShape(&m_EMB_SubdSIM_Shape);
            m_EMB_SubdSIM_Object.SetTransform( geo::Transform2::Identity() );
            m_EMB_Object.Embed( &m_EMB_SubdSIM_Object, m_Params.m_EMB_Method );
            //\todo here we could DISCARD EMB_SubdSIM, we'l re-subdivide on Redraw_EMB()
#else
            // Recompute embedding (\todo with Tr = Identity, by now)
            geo::Transform2 deformer_tr = m_SIM_Object.GetTransform();
            m_SIM_Object.SetTransform( geo::Transform2::Identity() );
            m_EMB_Object.Embed( &m_SIM_Object, m_Params.m_EMB_Method );
            m_SIM_Object.SetTransform( deformer_tr );
#endif
        }

    inline void Rebuild_DCR()
        {
            APP_LOG( "Rebuild_DCR()" );
            //\todo Add CD info to SIM
            //\todo No... pq NO VULL UNA POLYGONAL, vull dades INDEPENDENTS per cada element de SIM! \todo COULD BE USEFUL to generate partitioned VIZ, though...
            //       geo::Make_PolygonalShape2_From_PolygonalShape2_Clipped_With_MeshSolidShape2( m_CD_Shape, m_VIZ_Shape, m_SIM_Shape );
            //m_SIM_Shape.AddDCR( m_VIZ_Shape ); //\todo Choose eDCR_Partition, eDCR_DistanceField...
            //\todo VizDCR( ) //Draw DCR patches in different colors to show element transitions...

            if( m_SIM_Shape.AddDCR( &m_VIZ_Shape, geo::Transform2::Identity() ) ) //\note tr is identity because tr_s2mss is relative and they are both in globals here
                m_pDCR = m_SIM_Shape.GetDCR();
            // Log
            if( m_pDCR )
            {
                APP_LOG( "VIZ with %d V", m_VIZ_Shape.GetNumVertices() );
                APP_LOG( "SIM with %d V, %d HE, %d P", m_SIM_Shape.GetNumV(), m_SIM_Shape.GetNumHE(), m_SIM_Shape.GetNumP() );
                APP_LOG( "DCR(SIM,VIZ) with %d elements, %d patches, %d segments and %d vertices",
                         m_pDCR->m_NumElements, m_pDCR->m_NumPatches, m_pDCR->m_NumSegments, m_pDCR->m_NumVertices );
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
            geo::VizObject( &m_SRC_Object, m_VizIS );
            geo::VizObject( &m_VIZ_Object, m_VizIS );
            geo::VizObject( &m_EXT_Object, m_VizIS );
            geo::VizObject( &m_SIM_Object, m_VizIS, geo::eODDF_Shape | geo::eSDDF_Interior | geo::eSDDF_DCR );// | geo::eSDDF_Topology );
            Redraw_DCR();
            Redraw_EMB(); //\todo BEFORE Redraw_BVH() because it REFITS IT!!
#ifdef __ENABLE_BVH_TEST
            geo::VizBVH( m_BVH, m_VizIS, geo::eBVHDDF_Default | geo::eBVHDDF_Levels );
#endif
        }

    void Redraw_EMB()
        {
            /*\todo We do not have access to this info, should do it inside geo::VizObject( m_EMB_Object ), but THERE the info is also hidden...
            geo::EditableEmbeddedVerticesInTriangles m_EVIT;
            for( unsigned int it_tevs=0; it_tevs < m_EVIT.m_NumTEVS; it_tevs++ )
            {
                const geo::BakedTriangleEmbeddedVertexSubset& btevs( m_EVIT.m_vecTEVS[it_tevs] );
                for( unsigned int it_ev=0; it_ev < btevs.m_NumEmbeddedVertices; it_ev++ )
                {
                    geo::feature_index_type vid( m_EVIT.m_vecEmbeddedVID[ btevs.m_FirstEmbeddedVertexIndex + it_ev ] );
                    VIZ_SEGMENT2( m_VizIS, btevs.m_TriangleBarycenter, m_EMB_Object.GetShape()->V_Pos_0(vid), 1, Vec4f(1,1,1,1) );
                }
            }
            */

#ifdef __ENABLE_EMB_SubdSIM
            // Distort Deformer, redraw embedded (automatically synced) and restore Deformer afterwards
            if( 0 != m_SIM_Object.GetShape() && 0 != m_EMB_Object.GetShape() )
            {
                mal::SetRandomSeed(666);
                // Create a NEW Subdivided SIM with original nodes in
                // distorted positions, subdivide it, copy ALL node
                // positions to m_EMB_SubdSIM and use this one to
                // update EMB (internally embedded to exactly that
                // object, therefore we cannot rebuild it directly!)
                geo::EditableMeshSolidShape2 tmp_subd_shape;
                tmp_subd_shape.Set( m_SIM_Shape );
                // Distort SIM nodes
                tmp_subd_shape.BeginEdition();
                {
                    for( uint32 it_v=0; it_v < m_SIM_Shape.GetNumV(); it_v++ )
                        tmp_subd_shape.SetVertex( it_v, m_SIM_Shape.V_Pos_0(it_v) + m_Params.m_SRC_NaturalSize * m_Params.m_EMB_DistortionScale * mal::RandomUnitVec<geo::Real,2>() );
                }
                tmp_subd_shape.EndEdition();
                // Subd SIM keeping original nodes fixed (non-smoothed)
                for( uint32 it_subd=0; it_subd<m_Params.m_EMB_Subd_Iter; it_subd++ )
                    tmp_subd_shape.Subdivide( m_SIM_Shape.GetNumV() );
                // Save deformer SDOF and transform
                uint32 num_deformer_sdof = m_EMB_SubdSIM_Object.GetShape()->GetNumSDOF();
                geo::Vec2* vec_saved_deformer_sdof = new geo::Vec2[ num_deformer_sdof ];
                memcpy( vec_saved_deformer_sdof, m_EMB_SubdSIM_Object.GetVecSDOF(), sizeof(geo::Vec2)*num_deformer_sdof );
                geo::Transform2 deformer_tr = m_EMB_SubdSIM_Object.GetTransform();
                // Distort deformer SDOF
                for( uint32 it_sdof=0; it_sdof < num_deformer_sdof; it_sdof++ )
                    m_EMB_SubdSIM_Object.GetVecSDOF()[it_sdof] = tmp_subd_shape.GetVecDefaultSDOF()[it_sdof];
                //TEMP m_EMB_SubdSIM_Object.GetVecSDOF()[it_sdof] += m_Params.m_SRC_NaturalSize * m_Params.m_EMB_DistortionScale * mal::RandomUnitVec<geo::Real,2>();
                m_EMB_SubdSIM_Object.SetTransform( geo::Transform2( geo::Vec2( 3 * m_Params.m_SRC_NaturalSize, 0 ),
                                                                    geo::Mat2x2::Identity() ) );

                //TEMP Testing difference between classified-EMB and partitioned-DCR
                m_EMB_SubdSIM_Shape.AddDCR( &m_VIZ_Shape, geo::Transform2::Identity() );
                //TEMP
                // Draw distorted m_SIM_Object and restore SDOF afterwards
                geo::VizObject( &m_EMB_SubdSIM_Object, m_VizIS, geo::eODDF_Shape | geo::eSDDF_Interior );//| geo::eSDDF_DCR );
                geo::VizObject( &m_EMB_Object, m_VizIS, geo::eODDF_Shape | geo::eSDDF_Boundary | geo::eSDDF_Vertices );

                // Restore undistorted SDOF
                memcpy( m_EMB_SubdSIM_Object.GetVecSDOF(), vec_saved_deformer_sdof, sizeof(geo::Vec2)*num_deformer_sdof );
                m_EMB_SubdSIM_Object.SetTransform( deformer_tr );
                //\note m_EMB_Object transform is OUT OF SYNC here, but should be recomputed when accessed
                //m_EMB_Object.SyncEmbedding();
                // Dealloc
                delete[] vec_saved_deformer_sdof;
            }
#else
            // Distort Deformer, redraw embedded (automatically synced) and restore Deformer afterwards
            if( 0 != m_SIM_Object.GetShape() && 0 != m_EMB_Object.GetShape() )
            {
                // Save deformer SDOF and transform
                uint32 num_deformer_sdof = m_SIM_Object.GetShape()->GetNumSDOF();
                geo::Vec2* vec_saved_deformer_sdof = new geo::Vec2[ num_deformer_sdof ];
                memcpy( vec_saved_deformer_sdof, m_SIM_Object.GetVecSDOF(), sizeof(geo::Vec2)*num_deformer_sdof );
                geo::Transform2 deformer_tr = m_SIM_Object.GetTransform();
                // Distort deformer SDOF
                for( uint32 it_sdof=0; it_sdof < num_deformer_sdof; it_sdof++ )
                    m_SIM_Object.GetVecSDOF()[it_sdof] += m_Params.m_SRC_NaturalSize * m_Params.m_EMB_DistortionScale * mal::RandomUnitVec<geo::Real,2>();
                static float s_Angle(0);
                //TEMP: s_Angle += mal::Deg2Rad(5.0f);
                m_SIM_Object.SetTransform( geo::Transform2( geo::Vec2( 3 * m_Params.m_SRC_NaturalSize, 0 ),
                                                            mal::GRotation2x2_From( s_Angle ) ) );//mal::Mat2x2_Frgeo::Mat2x2::Identity() ) );
                // Draw distorted m_SIM_Object and restore SDOF afterwards
                geo::VizObject( &m_SIM_Object, m_VizIS, Flags32(geo::eODDF_Everything | geo::eSDDF_Everything).Disable(geo::eSDDF_DCR) );
                geo::VizObject( &m_EMB_Object, m_VizIS, geo::eODDF_Shape | geo::eSDDF_Boundary | geo::eSDDF_Vertices );
#  ifdef __ENABLE_BVH_TEST
                switch( m_Params.m_BVH_Method )
                {
                case Params::eBVHM_BV_E: m_BVH.Refit( boost::bind<void>( &geo::GEBV_MeshSolidShape2_E<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                         &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
                case Params::eBVHM_BV_BSLAB: m_BVH.Refit( boost::bind<void>( &geo::GEBV_MeshSolidShape2_BSlab_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                             &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
                case Params::eBVHM_BV_BDOP: m_BVH.Refit( boost::bind<void>( &geo::GEBV_MeshSolidShape2_BDOP_Safe<bvh_type::entry_index_type,bvh_type::bv_type>,
                                                                            &m_SIM_Shape, m_SIM_Object.GetTransform(), m_SIM_Object.GetVecSDOF(), _1, _2) ); break;
                default: break;
                }
#  endif
                // Restore undistorted SDOF
                memcpy( m_SIM_Object.GetVecSDOF(), vec_saved_deformer_sdof, sizeof(geo::Vec2)*num_deformer_sdof );
                m_SIM_Object.SetTransform( deformer_tr );
                //\note m_EMB_Object transform is OUT OF SYNC here, but should be recomputed when accessed
                //m_EMB_Object.SyncEmbedding();
                // Dealloc
                delete[] vec_saved_deformer_sdof;
            }
#endif //__ENABLE_EMB_SubdSIM
        }

    void Redraw_DCR()
        {
            /*TEMP: This is done automatically in VizObject( m_SIM_Object )
            if( m_pDCR == 0 ) return;
            geo::VizDCR( m_pDCR,
                         &m_SIM_Shape,
                         geo::Transform2( geo::Vec2( 2 * m_Params.m_SRC_NaturalSize, 0 ),
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
            // Register embedded shape //\todo This IGNORES EMBEDDING, as it's in the Object, not Shape, but saves the detailed/refined Polygonal shape
            sl.Register( m_EMB_Shape, "EMB" );
            // Save SL
            sl.Save( (m_OutputName + ".txt").c_str(), false );
            sl.Save( (m_OutputName + ".bin").c_str(), true ); //\todo redundant, keep while testing txt/bin formats
            // Save Params
            {
                //\todo InitArchetype may be neceessary if called from command-line, where no SetView() happens
                if( s_Cibulet_ArchetypeLibrary.IsEmpty() )
                {
                    Params::InitArchetype( s_Cibulet_ArchetypeLibrary );
                }
                util::ItemStream params_is;
                params_is.BeginComplex( "AppParams", eType_Property_Group );
                {
                    s_Cibulet_ArchetypeLibrary.ExportInstance( "Archetype_Cibulet_Params", &m_Params, params_is );
                }
                params_is.EndComplex();
                params_is.SaveTxt( (m_OutputName + ".params").c_str() );
            }
        }

    void OnLoadCall( void *p_user_data )
        {
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
        }

    void OnSubdCall( void *p_user_data )
        {
            m_SIM_Shape.Subdivide();
            m_SIM_Object.ResetShape( &m_SIM_Shape );
            m_Params.m_DirtyFlags.Enable( eDF_EMB );
            m_DirtyFlags.Enable( eDF_EMB );
        }

    void OnCropSimExtCall( void *p_user_data )
        {
#ifdef __ENABLE_CROP_SIM_EXT
            //-- Gather vec_pid of SIM_Ext elements and remove them
            std::vector<bool> vecIsCrustElement( m_SIM_Shape.GetNumP(), false );

            // Move both SIM and VIZ to origin, saving old Tr
            geo::Transform2 sim_tr = m_SIM_Object.GetTransform();
            geo::Transform2 viz_tr = m_VIZ_Object.GetTransform();
            m_SIM_Object.SetTransform( geo::Transform2::Identity() );
            m_VIZ_Object.SetTransform( geo::Transform2::Identity() );

            // Iterate over VIZ edges, testing them against SIM elements and flagging them as Crust.
            //\todo CHECK ROBUSTNESS and stability of edge-element tests (inverted segment returns same result?)
            geo::Vec2 a0( m_VIZ_Shape.V_Pos( 0, m_VIZ_Object.GetVecSDOF() ) );
            for( unsigned int it_edge=1; it_edge<m_VIZ_Shape.GetNumVertices(); it_edge++ )
            {
                geo::Vec2 a1( m_VIZ_Shape.V_Pos( it_edge, m_VIZ_Object.GetVecSDOF() ) );
                for( unsigned int it_p=0; it_p<m_SIM_Shape.GetNumP(); it_p++ )
                {
                    if( !vecIsCrustElement[it_p] ) //already-in-crust elements don't need to be tested again
                    {
                        unsigned int it_he( m_SIM_Shape.P_FirstHEID(it_p) );
                        geo::Vec2 b0( m_SIM_Shape.V_Pos( m_SIM_Shape.HE_OriginVID(it_he), m_SIM_Object.GetVecSDOF() ) );
                        do
                        {
                            it_he = m_SIM_Shape.HE_Next(it_he);
                            geo::Vec2 b1( m_SIM_Shape.V_Pos( m_SIM_Shape.HE_OriginVID(it_he), m_SIM_Object.GetVecSDOF() ) );
                            geo::Real lambdaA,lambdaB;
                            vecIsCrustElement[it_p] = geo::np::Intersection_Segment2_Segment2( a0, a1,
                                                                                               b0, b1,
                                                                                               lambdaA, lambdaB );
                            b0 = b1;
                        } while( it_he != m_SIM_Shape.P_FirstHEID(it_p)
                                 && !vecIsCrustElement[it_p] );
                    }
                }
                a0 = a1;
            }

            // Iterate over SIM elements and, if not flagged as Crust,
            // test if their barycenter (or any other point) is inside
            // VIZ. If not, add to array of EXT elements to be removed
            std::vector<geo::feature_index_type> vecExtPID;
            for( unsigned int it_p=0; it_p<m_SIM_Shape.GetNumP(); it_p++ )
                if( !vecIsCrustElement[it_p]
                    && !geo::np::Overlap_Point2_Polygonal2( m_SIM_Shape.P_Barycenter(it_p, m_SIM_Object.GetVecSDOF()),
                                                            &m_VIZ_Shape, m_VIZ_Object.GetTransform(), m_VIZ_Object.GetVecDOF(),
                                                            0 ) )
                    vecExtPID.push_back( it_p );

            // Restore transforms
            m_SIM_Object.SetTransform( sim_tr );
            m_VIZ_Object.SetTransform( viz_tr );
#else
            std::vector<geo::feature_index_type> vecExtPID;
            vecExtPID.push_back(0);
            vecExtPID.push_back(1);
            vecExtPID.push_back(2);
            vecExtPID.push_back(3);
#endif
            m_SIM_Shape.RemovePolygons( vecExtPID.size(), &vecExtPID[0] );
            m_SIM_Object.ResetShape( &m_SIM_Shape );
            m_Params.m_DirtyFlags.Enable( eDF_EMB | eDF_DCR );
            m_DirtyFlags.Enable( eDF_EMB | eDF_DCR );
        }

#ifdef __ENABLE_BVH_TEST
    bool TestBVH( const Vec2f& p ) const
    {
        bvh_type::bv_type bv( p );
        std::vector< bvh_type::entry_index_type > vec_entry;
        bool bHit = m_BVH.Test( bv, vec_entry );
        if( bHit ) APP_LOG("TestBVH hit!");
        return bHit;
    }
#endif

private:
    util::ItemStream m_VizIS;
    util::ItemStream m_TweakIS;
    sfr::IView *m_pAppView;

    // Args
    std::string m_FileNameSVG;
    std::string m_PathElementNameSVG;
    std::string m_OutputName;

    // State
    bool m_bIsActive;
    Flags32 m_DirtyFlags;

    // SRC
    geo::EditablePathShape2 m_SRC_Shape;
    geo::GObjectSS<geo::EditablePathShape2> m_SRC_Object;

    // VIZ( SRC )
    geo::EditablePolygonalShape2 m_VIZ_Shape;
    geo::GObjectSS<geo::EditablePolygonalShape2> m_VIZ_Object;

    // EXT( SRC )
    geo::EditablePolygonalShape2 m_EXT_Shape;
    geo::GObjectSS<geo::EditablePolygonalShape2> m_EXT_Object;

    // SIM( EXT, INT, \todo VIZ => CD )
    geo::EditableMeshSolidShape2 m_SIM_Shape;
    geo::GObjectSS<geo::EditableMeshSolidShape2> m_SIM_Object;

    // EMB( VIZ, SIM )
    geo::EditableMeshSolidShape2 m_EMB_SubdSIM_Shape;
    geo::GObjectSS<geo::EditableMeshSolidShape2> m_EMB_SubdSIM_Object;
    geo::EditablePolygonalShape2 m_EMB_Shape;
    geo::GObjectSS<geo::EditablePolygonalShape2> m_EMB_Object;

    // CD( VIZ, SIM )
    const geo::DCR_MeshSolidShape2* m_pDCR; //TEMP

#ifdef __ENABLE_BVH_TEST
    // BVH( CD )
    bvh_type m_BVH;
#endif

};

#endif //CIBULET_APP_TASK_H
