#ifndef SFR_SKR_H
#define SFR_SKR_H

#include <pla_types.h>
#include <util/ItemStream.h>
#include <util/ItemStreamSerialization.h>
#include <string>
#include <map>

namespace skr
{

//! VarSeq sequence definition
class VarSeq
{
public:
    static const unsigned int cMaxVarSeqDataMem = (1<<10);
    static const unsigned int cMaxVarSeqNameMem = (1<<10);

public:
    VarSeq() : m_pIS(0) {}
    ~VarSeq() { if(m_pIS) delete m_pIS; }

    void Init() { m_pIS = new util::ItemStream(cMaxVarSeqDataMem,cMaxVarSeqNameMem); }
    void Init( uint32 max_data_mem, uint32 max_name_mem )
    { m_pIS = new util::ItemStream(max_data_mem,max_name_mem); }

    //! Concatenable var addition
    template <typename T>
    inline VarSeq &Add( const char *var_name, const T& value )
    {
        m_pIS->Write( var_name, value );
        return *this;
    }
    template <typename T>
    inline VarSeq &AddArray( const char *var_name, const T* data, unsigned int num_values )
    {
        m_pIS->WriteArray( var_name, data, num_values );
        return *this;
    }
    inline VarSeq &AddString( const char *var_name, const char* str )
    {
        m_pIS->WriteString( var_name, str ); //includes '\0' ending mark
        return *this;
    }

    template <typename T>
    inline const T &Get( const char *var_name ) const
    {
        return m_pIS->Find(var_name).Get<T>();
    }
    template <typename T>
    inline T SafeGet( const char *var_name, const T &default_value ) const
    {
        return m_pIS->Find(var_name).SafeGet<T>( default_value );
    }
    template <typename T>
    inline const T *GetArrayPtr( const char *var_name ) const
    {
        return m_pIS->Find(var_name).GetArrayPtr<T>();
    }
    template <typename T>
    inline const T *SafeGetArrayPtr( const char *var_name, const T *default_ptr ) const
    {
        return m_pIS->Find(var_name).SafeGetArrayPtr<T>(default_ptr);
    }
    inline uint32 GetArrayCount( const char *var_name ) const
    {
        return m_pIS->Find(var_name).GetArrayCount();
    }
    inline const char *GetString( const char *var_name ) const
    {
        return m_pIS->Find(var_name).GetString();
    }
    inline const char *SafeGetString( const char *var_name, const char* default_str ) const
    {
        return m_pIS->Find(var_name).SafeGetString(default_str);
    }

    inline void Clear() { m_pIS->Clear(); }

    inline const util::ItemStream &GetIS() const { return *m_pIS; }

    /*! Matches a secondary varseq to a primary varseq and stores in "this"
      - Vars that are in primary and secondary get copied from secondary
      - Vars that are in primary but not in secondary get copied from primary
      - Vars that are in secondary but not in primary are ignored
    */
    void Match( const VarSeq &primary, const VarSeq &secondary );

private:
    util::ItemStream *m_pIS;
};

//! skr-Compatible method prototype
typedef bool(*MethodPtr)( VarSeq &user_context, VarSeq &results, const VarSeq &params );

class LanguageDef
{
public:
    static const unsigned int cMaxMethodNameLength = 32;
    static const unsigned int cMaxParamNameLength = 32;
    static const unsigned int cMaxNumberDigits = 16;
    static const unsigned int cMaxStringLength = 128;
    static const unsigned int cMaxArrayLength = 1024;

public:

    class MethodDef
    {

    public:
        MethodDef() : m_MethodPtr(0) {}
        ~MethodDef() {}

        //! Inits method definition and returns param definition varseq to be filled
        VarSeq &Init( const char *group, const char *name, const char *help, MethodPtr method_ptr )
        {
            m_Group = group;
            m_Name = name;
            m_Help = help;
            m_MethodPtr = method_ptr;
            m_ParamsDef.Init();
            return m_ParamsDef;
        }

        inline const std::string &GetGroup() const { return m_Group; }
        inline const std::string &GetName() const { return m_Name; }
        inline const std::string &GetHelp() const { return m_Help; }
        inline MethodPtr GetMethodPtr() const { return m_MethodPtr; }
        inline const VarSeq &GetParamsDef() const { return m_ParamsDef; }

    private:
        std::string m_Group;
        std::string m_Name;
        std::string m_Help;
        MethodPtr m_MethodPtr;
        VarSeq m_ParamsDef;
    };

public:
    LanguageDef();
    ~LanguageDef() {}

    VarSeq &DefineMethod( const char *group, const char *name, const char *help, MethodPtr method_ptr );
    const MethodDef *FindMethod( const char *full_name ) const;

private:
    static bool Method_Help( VarSeq &user_context, VarSeq &results, const VarSeq &params );
    static bool Method_Run( VarSeq &user_context, VarSeq &results, const VarSeq &params );
    static bool Method_Abort( VarSeq &user_context, VarSeq &results, const VarSeq &params );

private:
    typedef std::map< std::string, MethodDef > MapMethods;
    MapMethods m_mapMethods;
};


/*! Script execution context.
  \todo MAY BECOME RECURSIVE/STACK
*/
class Context
{
public:
    Context();
    ~Context() {}

    void Clear();

    const LanguageDef::MethodDef *m_pMethodDef;
    VarSeq m_Params;
    VarSeq m_MatchedParams;
    VarSeq m_Results;
};

/*! Parses a text buffer and interprets it according to a language
    definition.
*/
class Interpreter
{
public:

public:
    Interpreter( const LanguageDef &language_def )
    : m_rLanguageDef(language_def)
    {
        m_UserContext.Init(); //!\todo maybe allow user to specify max_context_mem
        m_UserContext.Add( "p_interpreter", reinterpret_cast<void_pointer_type>(this) ); //Required for command _run
    }
    ~Interpreter() {}

    //! Runs a script
    bool Execute( const char *script, unsigned int length = 0 );

    //! Complete command or param name
    void Autocomplete( const char *script, unsigned int length );

    //! Get user context to add/read variables
    VarSeq &GetUserContext() { return m_UserContext; }

    const LanguageDef &GetLanguageDef() const { return m_rLanguageDef; }

private:
    bool InterpretScript( const char *script, unsigned int length );

    const char *ParseStatement( const char *str );
    const char *ParseComment( const char *str );
    const char *ParseCommandName( const char *str );
    const char *ParseParamList( const char *str );
    const char *ParseParam( const char *str );
    const char *ParseNumber( const char *str, double &value );
    const char *ParseAndStoreNumber( const char *param_name, const char *str );
    const char *ParseAndStoreArrayNumber( const char *param_name, const char *str );
    const char *ParseAndStoreString( const char *param_name, const char *str );

    bool Execute( const LanguageDef::MethodDef &method,
                  VarSeq &user_context,
                  const VarSeq &params,
                  VarSeq &results );

private:
    //! _run( file="filename" );
    static bool Method_Run( VarSeq &user_context, VarSeq &results, const VarSeq &params );

private:
    const LanguageDef &m_rLanguageDef;
    Context m_Context;
    VarSeq m_UserContext;
};

} // namespace skr

#endif // SFR_SKR_H
