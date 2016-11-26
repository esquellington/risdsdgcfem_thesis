#include "PropertyTreeControl.h"
#include "SliderControl.h"
#include <stdio.h> //TEMPORAL: For snprintf
#include <vector>
#include <utility> //for std::pair

#define __ENABLE_STATS_MONITOR
#ifdef __ENABLE_STATS_MONITOR
#  include <fstream> //to output plots
#endif

namespace sfr { namespace gui {

class IPropertyControl: public HorizontalLayout
{
public:
    IPropertyControl( PropertyIt property_it )
    : m_PropertyIT(property_it)
        {
            AddChild( &m_labelValue );
            m_labelValue.SetLabel("unknown");
            m_labelValue.SetColor( gfx::Color(0,0,0,0) );
            m_labelValue.SetFontColor( gfx::Color(1,1,1,1) );
            m_labelValue.SetFieldLength( 8 );
            /*\note CANNOT call RecomputeLabel() here, because it
              requires the virtual method ValueToString() to be
              available, and it's NOT, because we're in base class'
              ctor, the derived class is NOT yet constructed and thus,
              ValueToString() specialization is NOT called.*/
        }
    virtual ~IPropertyControl() {}

    bool Update( float dt )
    {
        // Refresh value if property_it is Touched() from outside
        if( m_PropertyIT.IsTouched() )
        {
            RecomputeLabel();
            SyncEditor(); //to sync editing gui to changed value, eg: slider pos
            m_PropertyIT.Untouch();
            //SFR_LOG_WARNING( "Touched %s", m_PropertyIT.GetName() );
        }
        // Update layout
        return HorizontalLayout::Update(dt);
    }

    inline void Touched( void *p_user_data )
        {
            /*if(!mNotifyChangeDelegate.is_empty()) mNotifyChangeDelegate(mVarId);*/
            // Touch and recompute label immediately
            m_PropertyIT.Touch();
            RecomputeLabel();
        }

protected:
    void RecomputeLabel()
        {
            char str[LabelControl::cMaxLength];
            ValueToString( str, LabelControl::cMaxLength );
            m_labelValue.SetLabel( str );
        }
    virtual void ValueToString( char *str, int max_length ) const = 0;//{ snprintf( str, max_length, "unknown" ); }
    virtual void SyncEditor() {}

protected:
    PropertyIt m_PropertyIT;
    LabelControl m_labelValue;
};

class Unknown_PC: public IPropertyControl
{
public:
    Unknown_PC( PropertyIt property_it ) : IPropertyControl(property_it) { RecomputeLabel(); }
    ~Unknown_PC() {}
protected:
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "unknown" ); }
};

class Bool_PC: public IPropertyControl
{
public:
    Bool_PC( PropertyIt property_it )
    : IPropertyControl(property_it)
    , m_buttonToggle("(!)")
        {
            AddChild(&m_buttonToggle);
            m_buttonToggle.SetOnButtonCall( std::bind1st( std::mem_fun(&Bool_PC::Toggle), this) );
            RecomputeLabel();
        }
    ~Bool_PC() {}

    inline void Toggle( void *p_user_data )
        {
            bool bValue( m_PropertyIT.Get<bool>() );
            m_PropertyIT.Set<bool>( !bValue );
            Touched(0);
        }

protected:
    void ValueToString( char *str, int max_length ) const
        {
            snprintf( str, max_length, (m_PropertyIT.Get<bool>()) ? "true " : "false" );
        }

private:
    ButtonControl m_buttonToggle;
};

//---- Simple NIR (\todo DEPRECATED!!)
template <typename T>
class GNIR_PC: public IPropertyControl
{
private:
    typedef GProperty_NumberInRange<T> property_nir_type;
public:
    GNIR_PC( PropertyIt property_it )
    : IPropertyControl(property_it)
    , m_Slider( 100, property_it.Get<property_nir_type>() )
        {
            AddChild( &m_Slider );
            m_Slider.SetOnChangeCall( std::bind1st( std::mem_fun(&IPropertyControl::Touched), this) );
        }
    ~GNIR_PC() {}

    inline T Get() const { return m_PropertyIT.Get<property_nir_type>().m_Value; }
    void SyncEditor() { m_Slider.OnValueChanged(); }

private:
    GNIR_SliderControl<T> m_Slider;
};

class NIR_int32_PC: public GNIR_PC<int32>
{
public:
    NIR_int32_PC( PropertyIt property_it ) : GNIR_PC<int32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7d", Get() ); }
};

class NIR_uint32_PC: public GNIR_PC<uint32>
{
public:
    NIR_uint32_PC( PropertyIt property_it ) : GNIR_PC<uint32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7u", Get() ); }
};

class NIR_float32_PC: public GNIR_PC<float32>
{
public:
    NIR_float32_PC( PropertyIt property_it ) : GNIR_PC<float32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4f", Get() ); }
};

class NIR_float64_PC: public GNIR_PC<float64>
{
public:
    NIR_float64_PC( PropertyIt property_it ) : GNIR_PC<float64>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4lf", Get() ); }
};

//---- Complex NIR
template <typename T>
class GNIR2_PC: public IPropertyControl
{
public:
    GNIR2_PC( PropertyIt property_it )
    : IPropertyControl(property_it)
    , m_Slider( 100,
                property_it.GetSubItem().Find("value").Get<T>(),
                property_it.GetSubItem().Find("min").Get<T>(),
                property_it.GetSubItem().Find("max").Get<T>() )
        {
            m_ValueIt = property_it.GetSubItem().Find("value");
            AddChild( &m_Slider );
            m_Slider.SetOnChangeCall( std::bind1st( std::mem_fun(&IPropertyControl::Touched), this) );
        }
    ~GNIR2_PC() {}

    /*
    inline void Touched( void *p_user_data )
        {
            SFR_LOG_WARNING( "Touched NIR2 %s with value %f", m_PropertyIT.GetName(), m_ValueIt.Get<float32>() );
            IPropertyControl::Touched(p_user_data);
            SFR_ASSERT( m_PropertyIT.IsTouched() );
        }
    */

    inline T Get() const { return m_ValueIt.Get<T>(); }
    void SyncEditor() { m_Slider.OnValueChanged(); }

private:
    GSliderControl<T> m_Slider;
    PropertyIt m_ValueIt;
};

class NIR2_int32_PC: public GNIR2_PC<int32>
{
public:
    NIR2_int32_PC( PropertyIt property_it ) : GNIR2_PC<int32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7d", Get() ); }
};

class NIR2_uint32_PC: public GNIR2_PC<uint32>
{
public:
    NIR2_uint32_PC( PropertyIt property_it ) : GNIR2_PC<uint32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7u", Get() ); }
};

class NIR2_float32_PC: public GNIR2_PC<float32>
{
public:
    NIR2_float32_PC( PropertyIt property_it ) : GNIR2_PC<float32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4f", Get() ); }
};

class NIR2_float64_PC: public GNIR2_PC<float64>
{
public:
    NIR2_float64_PC( PropertyIt property_it ) : GNIR2_PC<float64>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4lf", Get() ); }
};

//---- Basic number types
template <typename T>
class GNumber_PC: public IPropertyControl
{
public:
    GNumber_PC( PropertyIt property_it ) : IPropertyControl(property_it) {}
    ~GNumber_PC() {}

    inline T Get() const { return m_PropertyIT.Get<T>(); }
};

class Number_int32_PC: public GNumber_PC<int32>
{
public:
    Number_int32_PC( PropertyIt property_it ) : GNumber_PC<int32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7d", Get() ); }
};

class Number_uint32_PC: public GNumber_PC<uint32>
{
public:
    Number_uint32_PC( PropertyIt property_it ) : GNumber_PC<uint32>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7u", Get() ); }
};

class Number_float32_PC: public GNumber_PC<float32>
{
public:
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4f", Get() ); }
    Number_float32_PC( PropertyIt property_it ) : GNumber_PC<float32>(property_it) { RecomputeLabel(); }
};

class Number_float64_PC: public GNumber_PC<float64>
{
public:
    Number_float64_PC( PropertyIt property_it ) : GNumber_PC<float64>(property_it) { RecomputeLabel(); }
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4lf", Get() ); }
};

//\todo Pair used as sample <time,value> pair
template <typename T1, typename T2>
class GPair_PC: public IPropertyControl
{
public:
    GPair_PC( PropertyIt property_it ) : IPropertyControl(property_it) {}
    ~GPair_PC() {}

    inline GPair<T1,T2> Get() const { return m_PropertyIT.Get< GPair<T1,T2> >(); }
};

#ifdef __ENABLE_STATS_MONITOR //\todo THIS should be generalized to All properties, at least all Non-Editable properties, used as stats
class Pair_float32_float32_PC: public GPair_PC<float32,float32>
{
public:
    void ValueToString( char *str, int max_length ) const { snprintf( str, max_length, "%7.4f", Get().m_Second ); }

    Pair_float32_float32_PC( PropertyIt property_it )
    : GPair_PC<float32,float32>(property_it)
    , m_buttonToggleRecord("O"), m_buttonPlot("~")
    , m_bIsRecording(false)
        {
            AddChild(&m_buttonToggleRecord);
            AddChild(&m_buttonPlot);
            m_buttonToggleRecord.SetOnButtonCall( std::bind1st( std::mem_fun(&Pair_float32_float32_PC::ToggleRecord), this) );
            m_buttonPlot.SetOnButtonCall( std::bind1st( std::mem_fun(&Pair_float32_float32_PC::Plot), this) );
            RecomputeLabel();
        }
    bool Update( float dt )
        {
            if( m_PropertyIT.IsTouched() )
            {
                RecomputeLabel();
                m_PropertyIT.Untouch();
            }
            if( m_bIsRecording )
                m_vecSamples.push_back( std::make_pair( Get().m_First, Get().m_Second ) );
            // Update layout
            return HorizontalLayout::Update(dt);
        }
    inline void ToggleRecord( void *p_user_data )
        {
            m_bIsRecording = !m_bIsRecording;
            if( m_bIsRecording ) { m_buttonToggleRecord.SetLabel("\""); m_vecSamples.clear(); }
            else m_buttonToggleRecord.SetLabel("O");
        }

    inline void Plot( void *p_user_data )
        {
            if( m_vecSamples.empty() ) return;
            std::string file_name( std::string("plot_testSaphyre_") + m_PropertyIT.GetName() + ".txt" );
            std::ofstream fs;
            fs.open( file_name.c_str() );
            for( std::vector< std::pair<float32,float32> >::const_iterator it_sample = m_vecSamples.begin(); it_sample != m_vecSamples.end(); it_sample++ ) fs << it_sample->first << " " << it_sample->second << std::endl;
            fs.close();
            system( (std::string("gnuplot -e \"plot '") + file_name + "'; \" -p").c_str() );
        }
private:
    ButtonControl m_buttonToggleRecord;
    ButtonControl m_buttonPlot;
    std::vector< std::pair<float32,float32> > m_vecSamples;
    bool m_bIsRecording;
};
#endif

//---- Complex Enum
class Enum_PC: public IPropertyControl
{
public:
    Enum_PC( PropertyIt property_it )
    : IPropertyControl(property_it)
    , m_ValueIndex(0)
        {
            m_ValueIt = property_it.GetSubItem().Find("value");
            AddChild( &m_Layout );
            util::ItemStream::ItemIt vec_value_it( property_it.GetSubItem().Find("vec_value") );
            util::ItemStream::ItemIt vec_name_it( property_it.GetSubItem().Find("vec_name") );
            SFR_ASSERT( vec_value_it.IsValid() && vec_value_it.IsArray() && vec_value_it.GetType() == pla_type_id<uint32>::value
                        && vec_name_it.IsValid() && vec_name_it.IsArray() && vec_name_it.GetType() == pla_type_id<String32>::value
                        && vec_value_it.GetArrayCount() == vec_name_it.GetArrayCount() );
            unsigned int count( vec_value_it.GetArrayCount() );
            const uint32 *vec_value( vec_value_it.GetArrayPtr<uint32>() );
            const String32 *vec_name( vec_name_it.GetArrayPtr<String32>() );
            for( unsigned int i=0; i<count; i++ )
            {
                // Init value index
                if( m_ValueIt.Get<uint32>() == vec_value[i] ) m_ValueIndex = i;
                // Create button with value index as user data
                ButtonControl *pBC = new ButtonControl( vec_name[i],
                                                        reinterpret_cast<void*>( static_cast<machine_uint_type>( i ) ) );
                                                        //reinterpret_cast<void*>( static_cast<machine_uint_type>( vec_value[i] ) ) );
                pBC->SetOnButtonCall( std::bind1st( std::mem_fun(&Enum_PC::SetEnumValueIndex), this) );
                m_Layout.AddChild( pBC );
            }
            RecomputeLabel();
        }
    ~Enum_PC()
        {
            m_Layout.RemoveAndDeleteChildren();
        }

    inline void SetEnumValueIndex( void *p_user_data )
        {
            m_ValueIndex = static_cast<uint32>( reinterpret_cast<machine_uint_type>(p_user_data) );
            m_ValueIt.Set<uint32>( m_PropertyIT.GetSubItem().Find("vec_value").GetArrayPtr<uint32>()[m_ValueIndex] );
            Touched(0);
        }

    /* TEMP: To trace touch calls
    inline void Touched( void *p_user_data )
        {
            SFR_LOG_WARNING( "Touched NIR2 %s with value %f", m_PropertyIT.GetName(), m_ValueIt.Get<float32>() );
            IPropertyControl::Touched(p_user_data);
            SFR_ASSERT( m_PropertyIT.IsTouched() );
        }
    */

    inline uint32 Get() const { return m_ValueIt.Get<uint32>(); }
    void SyncEditor() { }//m_Slider.OnValueChanged(); }
private:

    void ValueToString( char *str, int max_length ) const
        {
            snprintf( str, max_length, "%s",
                      util::ItemStream::ItemIt( m_PropertyIT ).GetSubItem().Find("vec_name").GetArrayPtr<String32>()[m_ValueIndex].GetStr() );
        }

private:
    HorizontalLayout /*VerticalLayout*/ m_Layout;
    PropertyIt m_ValueIt;
    uint32 m_ValueIndex;
};

//---- Complex Flags
class Flags_PC: public IPropertyControl
{
public:
    Flags_PC( PropertyIt property_it )
    : IPropertyControl(property_it)
        {
            m_ValueIt = property_it.GetSubItem().Find("value");
            AddChild( &m_Layout );
            util::ItemStream::ItemIt vec_value_it( property_it.GetSubItem().Find("vec_value") );
            util::ItemStream::ItemIt vec_name_it( property_it.GetSubItem().Find("vec_name") );
            SFR_ASSERT( vec_value_it.IsValid() && vec_value_it.IsArray() && vec_value_it.GetType() == pla_type_id<int32>::value
                        && vec_name_it.IsValid() && vec_name_it.IsArray() && vec_name_it.GetType() == pla_type_id<String32>::value
                        && vec_value_it.GetArrayCount() == vec_name_it.GetArrayCount() );
            unsigned int count( vec_value_it.GetArrayCount() );
            const int32 *vec_value( vec_value_it.GetArrayPtr<int32>() );
            const String32 *vec_name( vec_name_it.GetArrayPtr<String32>() );
            Flags32 value( Get() );
            for( unsigned int i=0; i<count; i++ )
            {
                // Create button with value index as user data
                ButtonControl *pBC = new ButtonControl( vec_name[i],
                                                        reinterpret_cast<void*>( static_cast<machine_uint_type>( vec_value[i] ) ) );
                pBC->SetOnButtonCall( std::bind1st( std::mem_fun(&Flags_PC::ToggleFlag), this) );
                pBC->SetElastic(false); //button behaves like a checkbox
                pBC->SetPressed( value.Test( vec_value[i] ) );
                m_Layout.AddChild( pBC );
            }
            RecomputeLabel();
        }
    ~Flags_PC()
        {
            m_Layout.RemoveAndDeleteChildren();
        }

    inline void ToggleFlag( void *p_user_data )
        {
            int32 flag( static_cast<int32>( reinterpret_cast<machine_uint_type>(p_user_data) ) );
            m_ValueIt.Set<Flags32>( Get().Toggle( flag ) );
            Touched(0);
        }

    /* TEMP: To trace touch calls
    inline void Touched( void *p_user_data )
        {
            SFR_LOG_WARNING( "Touched NIR2 %s with value %f", m_PropertyIT.GetName(), m_ValueIt.Get<float32>() );
            IPropertyControl::Touched(p_user_data);
            SFR_ASSERT( m_PropertyIT.IsTouched() );
        }
    */

    inline Flags32 Get() const { return m_ValueIt.Get<Flags32>(); }
    void SyncEditor() { }
private:

    void ValueToString( char *str, int max_length ) const
        {
            snprintf( str, max_length, "%7x", (uint32)Get() );
        }

private:
    HorizontalLayout /*VerticalLayout*/ m_Layout;
    PropertyIt m_ValueIt;
};

PropertyTreeControl::PropertyTreeControl( PropertyIt property_tree_it, D_SyncCall sync_call )
: m_PropertyTreeIT(property_tree_it)
, m_TreeControl(10,5,0)
, m_SyncCall(sync_call)
, m_PropertyTreeAPI(this)
{
    AddChild( &m_TreeControl );
    Rebuild( property_tree_it );
}

void PropertyTreeControl::Clear()
{
    m_TreeControl.Clear( TreeItemId(TreeControl::cRootItemId) );
}

void PropertyTreeControl::Rebuild( PropertyIt property_tree_it )
{
    Clear();
    //\note Root property is NOT added as a subtree to avoid unnecessary structure
    for( PropertyIt sub_it=property_tree_it.GetSubItem(); sub_it.IsValid(); ++sub_it )
        if( sub_it.IsComplex() ) AddComplexProperty( sub_it, TreeItemId(TreeControl::cRootItemId) );
        else AddSimpleProperty( sub_it, TreeItemId(TreeControl::cRootItemId) );
    // Unroll tree once populated
    m_TreeControl.SetRolled( TreeItemId(TreeControl::cRootItemId), false );

    // Roll 1st level properties
    TreeItemId tid( m_TreeControl.GetFirstChildItem( TreeItemId(TreeControl::cRootItemId) ) );
    while( tid != 0 )
    {
        m_TreeControl.SetRolled( tid, true );
        tid = m_TreeControl.GetNextSiblingItem( tid );
    }
}

void PropertyTreeControl::AddComplexProperty( PropertyIt pit, TreeItemId tid )
{
    SFR_ASSERT( pit.IsComplex() );
    switch( pit.GetType() )
    {
    case eType_Property_NIR:
        AddComplexProperty_NIR( pit, tid );
        break;
    case eType_Property_Enum:
        {
            TreeItemId property_tid = m_TreeControl.Add( pit.GetName(), 0, tid );
            m_TreeControl.AddControl( property_tid, new Enum_PC(pit) );
        }
        break;
    case eType_Property_Flags:
        {
            TreeItemId property_tid = m_TreeControl.Add( pit.GetName(), 0, tid );
            m_TreeControl.AddControl( property_tid, new Flags_PC(pit) );
        }
        break;
    case eType_Property_Object:
    case eType_Property_Group:
    default:
        {
            SFR_ASSERT( eType_Property_Object == pit.GetType()
                        || eType_Property_Group == pit.GetType() );
            TreeItemId group_tid = m_TreeControl.Add( pit.GetName(), 0, tid );
            for( PropertyIt sub_it=pit.GetSubItem(); sub_it.IsValid(); ++sub_it )
                if( sub_it.IsComplex() ) AddComplexProperty( sub_it, group_tid );
                else AddSimpleProperty( sub_it, group_tid );
            // Unroll group once populated
            m_TreeControl.SetRolled( group_tid, false );
        }
        break;
    }
}

void PropertyTreeControl::AddComplexProperty_NIR( PropertyIt pit, TreeItemId tid )
{
    // Create specific property control
    IPropertyControl *pPC(0);
    switch( pit.GetSubItem().Find("value").GetType() )
    {
    case eType_Int32: pPC = new NIR2_int32_PC(pit); break;
    case eType_UInt32: pPC = new NIR2_uint32_PC(pit); break;
    case eType_Float32: pPC = new NIR2_float32_PC(pit); break;
    case eType_Float64: pPC = new NIR2_float64_PC(pit); break;
    default: pPC = new Unknown_PC(pit); break;
    }
    TreeItemId property_tid = m_TreeControl.Add( pit.GetName(), 0, tid );
    m_TreeControl.AddControl( property_tid, pPC );
}

void PropertyTreeControl::AddSimpleProperty( PropertyIt pit, TreeItemId tid )
{
    // Create specific property control
    IPropertyControl *pPC(0);
    switch( pit.GetType() )
    {
    case eType_Bool32: pPC = new Bool_PC(pit); break;

    case eType_Property_NIR_Int32: pPC = new NIR_int32_PC(pit); break;
    case eType_Property_NIR_UInt32: pPC = new NIR_uint32_PC(pit); break;
    case eType_Property_NIR_Float32: pPC = new NIR_float32_PC(pit); break;
    case eType_Property_NIR_Float64: pPC = new NIR_float64_PC(pit); break;

    case eType_Int32: pPC = new Number_int32_PC(pit); break;
    case eType_UInt32: pPC = new Number_uint32_PC(pit); break;
    case eType_Float32: pPC = new Number_float32_PC(pit); break;
    case eType_Float64: pPC = new Number_float64_PC(pit); break;

#ifdef __ENABLE_STATS_MONITOR
    case eType_Pair_Float32_Float32: pPC = new Pair_float32_float32_PC(pit); break;
#endif

    default: pPC = new Unknown_PC(pit); break;
    }
    TreeItemId property_tid = m_TreeControl.Add( pit.GetName(), 0, tid );
    m_TreeControl.AddControl( property_tid, pPC );
}

PropertyTreeControl::~PropertyTreeControl()
{
}

bool PropertyTreeControl::Update( float dt )
{
    // Sync client properties to PTC-touched data
    if( !m_SyncCall.empty() ) m_SyncCall( m_PropertyTreeIT );
    // Sync PTC properties to client-touched data
    return FrameLayout::Update(dt);
}

}} //namespace sfr::gui
