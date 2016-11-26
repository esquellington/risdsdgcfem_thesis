#ifndef SFR_UTIL_IOBSERVER_H
#define SFR_UTIL_IOBSERVER_H

#include <list>

namespace sfr { namespace Util
{

class Observable;

class IObserver
{
public:
    IObserver() {}
    virtual ~IObserver() {}

    //! \name Observables management methods
    //@{
    void AddObservable( Observable *p_observable )
    {
        m_listObservables.push_back( p_observable );
        p_observable->AddObserver( this );
    }

    void RemoveObservable( const Observable *p_observable )
    {
        //todo
    }
    //@}

    
    virtual void NotifyChanged( const Observable *p_observable )
    {
        // Maybe find and confirm observable in list...
    }
    
private:
    std::list<const Observable *> m_listObservables;
};

} } // namespace sfr::Util

#endif //SFR_UTIL_IOBSERVER_H
