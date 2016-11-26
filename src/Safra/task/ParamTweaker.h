#ifndef SFR_TASK_PARAM_TWEAKER_H
#define SFR_TASK_PARAM_TWEAKER_H

#include <Safra/core/ITask.h>
#include <Safra/core/IView.h>

#include <pla_types.h>
#include <util/ItemStream.h>

#include "../../../../local/AntTweakBar/include/AntTweakBar.h"

namespace sfr
{

class ParamTweaker: public ITask
{
public:
    ParamTweaker( const util::ItemStream *p_is, IView *p_view )
    : m_pIS(p_is)
    , m_pView(p_view)
    , m_pBar(0)
    {
        if( m_pView->GetGUI() )
            m_pBar = static_cast<TwBar*>( m_pView->GetGUI()->CreatePanel("sfr::ParamTweaker",0)->GetATB() );
    }
    ~ParamTweaker() {}

    void Run( unsigned int tick )
    {
        if( !m_pIS || !m_pBar )
            return;        
        
        //Read new entries from ParamItemStream
        char str[256];
        for( util::ItemStream::ItemIt it = m_pIS->Begin(); it.IsValid(); ++it )
        {
            if( it.IsPointer() )
            {
                // Direct internal param pointers with no more info
                switch( it.GetType() )
                {                    
                case eTypeFloat:
                    {
                        const float *p_param = it.GetPtr<float>();
                        TwAddVarRW( m_pBar, it.GetName(), TW_TYPE_FLOAT, (void*)p_param,
                                    " min=0.0 max=10000 step=0.5 keyIncr=z keyDecr=Z help='HILFE!' ");
                    }                    
                    break;
                default:
                    break;
                }
            }            
            else if( it.IsComplex()
                     && it.GetType() == eTypeParamGroup )
            {
                for( util::ItemStream::ItemIt it2 = it.GetSubItem(); it2.IsValid(); ++it2 )
                {                    
                    switch( it2.GetType() )
                    {                    
                    case eTypeParamInt32:
                        {
                            const GSimpleParameter<int32> &param = it2.Get< GSimpleParameter<int32> >();
                            sprintf( str, " min=%f max=%f step=%f", (float)param.m_Min, (float)param.m_Max, (float)param.m_Step );
                            TwAddVarRW( m_pBar, it2.GetName(), TW_TYPE_INT32, param.m_Ptr, str );
                        }                    
                        break;
                    case eTypeParamFloat:
                        {
                            const GSimpleParameter<float> &param = it2.Get< GSimpleParameter<float> >();
                            sprintf( str, " min=%f max=%f step=%f", (float)param.m_Min, (float)param.m_Max, (float)param.m_Step );
                            TwAddVarRW( m_pBar, it2.GetName(), TW_TYPE_FLOAT, param.m_Ptr, str );
                        }
                        break;
                    default:
                        break;
                    }
                }
            }
        }
        // Clear ParamItemStream when tweakers have been created
        //m_pIS->Clear();
        m_pIS = NULL; //!<TEMPORAAAAAAAAAAAAAAAAAAAAAAAAAAL!!
    }
    
private:
    const util::ItemStream *m_pIS;
    IView *m_pView;
    TwBar *m_pBar;
};

} // namespace sfr

#endif // SFR_TASK_PARAM_TWEAKER_H
