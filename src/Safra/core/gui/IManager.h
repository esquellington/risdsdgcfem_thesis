#ifndef SFR_CORE_GUI_IMANAGER_H
#define SFR_CORE_GUI_IMANAGER_H

namespace sfr { namespace gui
{

class IPanel
{
public:
    IPanel() {}
    virtual ~IPanel() {}

    virtual void *GetATB() = 0; //!< Ugly hack to avoid wrapping the whole AntTweakBar stuff
};

class IManager
{
public:
    IManager() {}
    virtual ~IManager() {}

    virtual IPanel *CreatePanel( const char *name, unsigned int flags ) = 0;
    virtual void Clear() = 0;
    virtual bool Draw() = 0;
};

} } // namespace sfr::gui

#endif // SFR_CORE_GUI_IMANAGER_H
