#ifndef SFR_UTIL_OBSERVABLE_H
#define SFR_UTIL_OBSERVABLE_H

#include <list>

namespace sfr { namespace Util
{

class IObserver;

class Observable
{
    friend class IObserver;

public:
    Observable() {}
    ~Observable() {}

protected:
    //! \name Observable interface accessible only to subclasses and to IObserver
    //@{
    void AddObserver( IObserver *p_observer );
    void RemoveObserver( IObserver *p_observer );
    void NotifyChanged() const;
    //@}
    
private:
    std::list<IObserver *> m_listObservers;
};

} } // namespace sfr::Util               

#endif //SFR_UTIL_OBSERVABLE_H
