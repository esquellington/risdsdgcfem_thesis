#include <Safra/util/Observable.h>
#include <Safra/util/IObserver.h>

namespace sfr { namespace Util
{

void Observable::AddObserver( IObserver *p_observer )
{
    m_listObservers.push_back( p_observer );
}

void Observable::RemoveObserver( IObserver *p_observer )
{
    //todo
}

void Observable::NotifyChanged() const
{
    for( std::list<IObserver*>::const_iterator it = m_listObservers.begin();
         it != m_listObservers.end();
         ++it )
        (*it)->NotifyChanged( this );
}
    
} } // namespace sfr::Util               
