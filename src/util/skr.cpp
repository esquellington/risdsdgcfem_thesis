#include "skr.h"
#include <util/StringUtils.h>
#include <string.h>
#include <stdio.h> //TEMPORAL: For fopen/fread

namespace skr
{

//---- VarSeq Implementation
void VarSeq::Match( const VarSeq &primary, const VarSeq &secondary )
{
    //For each primary var
    for( util::ItemStream::ItemIt it_primary=primary.m_pIS->Begin();
         it_primary.IsValid();
         ++it_primary )
    {
        /*Copy from secondary if exists and matches type...or set to primary otherwise
          Match(a,b) = Same full type
          Arrays must:
          - Match in size (same full type)
          - Be longer (excess values will be silently ignored)
          - Be strings with possibly different lengths
        */
        const char *name = it_primary.GetName();
        util::ItemStream::ItemIt it_secondary = secondary.m_pIS->Find(name);
        if( it_secondary.IsValid()
            &&
            ( it_secondary.GetFullType() == it_primary.GetFullType()
              ||
              ( it_secondary.GetType() == pla_type_id<char>::value
                && it_secondary.GetType() == it_primary.GetType()
                && it_secondary.IsArray() && it_primary.IsArray() )
              ||
              ( it_secondary.GetType() == it_primary.GetType()
                && it_secondary.IsArray() && it_primary.IsArray()
                && it_secondary.GetArrayCount() >= it_primary.GetArrayCount() ) )
            )
            m_pIS->WriteItem( name, it_secondary );
        else
        {
            m_pIS->WriteItem( name, it_primary );
            if( it_secondary.IsValid() )
                std::cout << "[ERROR] Wrong parameter type for "
                          << it_secondary.GetName() << std::endl;
        }
    }

    /* Add "FreeParams" from secondary. Any unmatched param could
      be considered a free param, but better require a prefix "%"
      check prefix to avoid confusion. No type-check needed for
      FreeParams.
    */
    for( util::ItemStream::ItemIt it_secondary=secondary.m_pIS->Begin();
         it_secondary.IsValid();
         ++it_secondary )
    {
        const char *name = it_secondary.GetName();
        if( name[0] == '%' )
        {
            if( primary.m_pIS->Find(name).IsValid() )
                std::cout << "[ERROR] FreeParam " << name << "is also declared as Primary. Avoid '%' prefix in Primary params!" << std::endl;
            else
                m_pIS->WriteItem( name, it_secondary );
        }
    }
}

//---- LanguageDef Implementation
LanguageDef::LanguageDef()
{
    // Define builtin methods
    DefineMethod("_","help","prints help text",&Method_Help)
        .Add("language_def",reinterpret_cast<void_pointer_type>(this))
        .AddString("concept","*");
    // Run script. Requires "interpreter" to be defined in user context
    DefineMethod("_","run","runs skr file",&Method_Run)
        .AddString("file","");
    // Run script. Requires "interpreter" to be defined in user context
    DefineMethod("_","abort","aborts program immediately",&Method_Abort)
        .Add("rv",double(0));
}

VarSeq &LanguageDef::DefineMethod( const char *group,
                                   const char *name,
                                   const char *help,
                                   MethodPtr method_ptr )
{
    MethodDef &md = m_mapMethods.insert(std::make_pair(std::string(group) + name,MethodDef())).first->second;
    return md.Init(group,name,help,method_ptr);
}

const LanguageDef::MethodDef *LanguageDef::FindMethod( const char *full_name ) const
{
    MapMethods::const_iterator it = m_mapMethods.find(std::string(full_name));
    if( it != m_mapMethods.end() )
        return &it->second;
    else
        return 0;
}


bool LanguageDef::Method_Help( VarSeq &user_context, VarSeq &results, const VarSeq &params )
{
    //std::cout << "Method_Help()" << std::endl;
    const LanguageDef *ld = reinterpret_cast<const LanguageDef*>( params.Get<void_pointer_type>("language_def") );
    std::string concept( params.GetString("concept") );
    std::cout << "Help on concept: " << concept << std::endl;
    if( concept == "*" )
    {
        for( MapMethods::const_iterator it_method=ld->m_mapMethods.begin();
             it_method != ld->m_mapMethods.end();
             ++it_method )
        {
            std::cout << "\tMethod: " << it_method->second.GetGroup() << it_method->second.GetName()
                      << "\tDox: " << it_method->second.GetHelp() << std::endl;
        }
    }
    else
    {
        MapMethods::const_iterator it_method( ld->m_mapMethods.find(concept) );
        if( it_method != ld->m_mapMethods.end() )
            std::cout << "\tMethod: " << it_method->second.GetGroup() << it_method->second.GetName()
                      << "\tDox: " << it_method->second.GetHelp() << std::endl;
        else
            std::cout << "\tConcept Unknown!!" << std::endl;
    }
    return true;
}

bool LanguageDef::Method_Run( VarSeq &user_context, VarSeq &results, const VarSeq &params )
{
    //std::cout << "Method_Run()" << std::endl;
    //std::cout << "UC = " << std::endl << user_context.GetIS() << std::endl;
    Interpreter *pInterpreter( reinterpret_cast<Interpreter*>( user_context.Get<void_pointer_type>("p_interpreter") ) );
    std::string file_name( params.GetString("file") );
    if( file_name != "" )
    {
        FILE *pFile = fopen( file_name.c_str(), "r" );
        if( pFile )
        {
            //std::cout << "Running skr file: " << file_name << std::endl;
            // Get filesize
            fseek( pFile , 0 , SEEK_END );
            size_t file_size = ftell( pFile );
            rewind( pFile );
            //std::cout << "file_size: " << file_size << std::endl;
            // Read into memory buffer
            char *file_buffer = new char[file_size+1];
            if( 0 == file_buffer )
            {
                std::cout << "OUT OF MEMORY allocating file buffer of size " << file_size << std::endl;
                fclose( pFile );
                return false;
            }
            fread( file_buffer, 1, file_size, pFile );
            file_buffer[file_size] = '\0';
            // Run skr file
            //\todo we may need to save Context here
            pInterpreter->Execute( file_buffer, file_size );
            //\todo we may need to recover Context here, as the results are "_run" results, NOT last _runned script command results...
            // Cleanup
            delete [] file_buffer;
            fclose( pFile );
            return true;
        }
        else
        {
            std::cout << "Cannot open skr file " << file_name << std::endl;
            return false;
        }
    }
    else
    {
        std::cout << "Empty skr filename!" << std::endl;
        return false;
    }
}

bool LanguageDef::Method_Abort( VarSeq &user_context, VarSeq &results, const VarSeq &params )
{
    int return_value( static_cast<int>( params.Get<double>("rv") ) );
    UTIL_LOG_WARNING( "Aborting with rv = %d", return_value );
    exit( return_value );
}

//---- Context Implementation

Context::Context()
{
    m_pMethodDef = 0;
    m_Params.Init();
    m_MatchedParams.Init();
    m_Results.Init();
}
void Context::Clear()
{
    m_pMethodDef = 0;
    m_Params.Clear();
    m_MatchedParams.Clear();
    m_Results.Clear();
}


//---- Interpreter Implementation

/*! Runs a script
  Performs the following sequence of tasks:
  - Split input text in tokens
    - Deduce token type from first character (lexical analysis)
      - Find method to be executed
      - Build method parameter VarSeq
      - Execute method
    - Stop at statement end ';'

    Grammar:
      - Script => { Statement }*
      - Statement => command_name ParamsList ;
      - ParamList => { \lambda | ( Param {,Param}* ) }
      - Param => _param_name Value
      - Value => { Number | Array<Number> | String | $var_name }
      - Numer => {digit}+{ \lambda | .{digit}+ }
      - Array => [ Number {,Number}* ]
      - String => ".."
      - Comment => /.....\n

      Parser:
      - ParseScript
        - ParseComment
        - ParseStatement
          - ParseCommandName
          - ParseParams
            - ParseParamName
            - ParseValue (switch type)
*/
bool Interpreter::Execute( const char *script, unsigned int length )
{
    unsigned int actual_length = length;
    if( 0 == length ) actual_length = strlen(script);
    //std::cout << "Interpreting [" << actual_length << "] bytes:\n>> " << script << std::endl;
    // Interpret script
    return InterpretScript(script,actual_length);
}

//!< Complete command or param name
void Interpreter::Autocomplete( const char *script, unsigned int length )
{
}

bool Interpreter::InterpretScript( const char *script, unsigned int length )
{
    // Interpret script
    bool bResult(true);
    const char *p_statement = util::Skip(script," ;");
    while( bResult && 0 != p_statement )
    {
        // Init statement context
        m_Context.Clear();

        // Parse Comment (may be empty, may consume all tokens too)
        const char *p_comment_end = ParseComment(p_statement);
        if( 0 == p_comment_end )
            return bResult;
        else if( p_comment_end == p_statement ) //no comment
        {
            // Parse statement
            const char *p_end = ParseStatement(p_comment_end);

            // Execute if successful
            if( 0 != p_end )
            {
                if( 0 != m_Context.m_pMethodDef )
                    bResult = Execute( *m_Context.m_pMethodDef,
                                       m_UserContext,
                                       m_Context.m_Params,
                                       m_Context.m_Results );
                else
                    std::cout << "[ERROR] Unknown command_name" << std::endl;
                //Skip to next statement
                p_statement = util::Skip(p_end," ;");
            }
            else
            {
                std::cout << "[ERROR] Cannot parse Statement '" << p_statement << "'" << std::endl;
                bResult = false;
            }
        }
        else //Skip to next statement
            p_statement = util::Skip(p_comment_end," ;");
    }
    return bResult;
}

const char *Interpreter::ParseStatement( const char *str )
{
    // Parse CommandName (mandatory)
    const char *p_name_end = ParseCommandName(str);
    if( 0 == p_name_end ) return 0;

    // Parse ParamList (optional)
    const char *p_params_end = ParseParamList(p_name_end);
    return p_params_end;
}

const char *Interpreter::ParseComment( const char *str )
{
    // Skip spaces, find '/', then skip whole line if found, and
    // return p_end AFTER the comment-line
    const char *p_begin = util::Skip(str,' ');
    if( *str == '/' )
    {
        const char *p_end = util::Find(p_begin,"\n");
        if( p_end != 0 )
        {
            /* Just ignore comments...
            unsigned int length = p_end - p_begin;
            if( length < 128 )
            {
                char tmp_str[128];
                strncpy( tmp_str, p_begin, length );
                tmp_str[length] = '\0';
                std::cout << "[COMMENT] '" << tmp_str << "'" << std::endl;
            }
            */
            ++p_end;
        }
        return p_end;
    }
    else
        return p_begin;
}

const char *Interpreter::ParseCommandName( const char *str )
{
    const char *p_begin = util::Skip(str,' ');
    if( 0 == p_begin ) return 0;
    const char *p_end = util::Find(p_begin," (;");
    if( 0 == p_end ) return 0;
    unsigned int length = p_end - p_begin;
    if( util::IsLetter(*p_begin) )
    {
        if( length < LanguageDef::cMaxMethodNameLength )
        {
            char command_name[LanguageDef::cMaxMethodNameLength];
            strncpy( command_name, p_begin, length );
            command_name[length] = '\0';
            m_Context.m_pMethodDef = m_rLanguageDef.FindMethod(command_name);
            if( 0 == m_Context.m_pMethodDef )
                std::cout << "[ERROR] Unknown command_name '"
                          << command_name << "'" << std::endl;
            return p_end;
        }
        else
        {
            std::cout << "[ERROR] command_name too long > "
                      << LanguageDef::cMaxMethodNameLength << std::endl;
            return 0;
        }
    }
    else
    {
        std::cout << "[ERROR] Cannot parse command_name '" << p_begin << "'" <<std::endl;
        return 0;
    }
    return p_end;
}

const char *Interpreter::ParseParamList( const char *str )
{
    // Search for beginning of param-list or end of statement
    const char *p_begin = util::Find(str,";(");
    // If '(' not found BEFORE end of statement ';', param-list is void.
    if( 0 == p_begin || *p_begin == ';' ) return str;

    //\todo length = find ')', m_Context.m_Params.Parse( p_begin, length );

    // Otherwise, parse params individually until end of param-list
    p_begin = util::Skip(p_begin," (");
    while( 0 != p_begin && *p_begin != ')' )
    {
        const char *p_end = ParseParam(p_begin);
        if( 0 == p_end ) return 0;
        p_begin = util::Skip(p_end,", ");
    }
    // Jump last ')' if found
    if( 0 != p_begin ) p_begin++;
    return p_begin;
}

const char *Interpreter::ParseParam( const char *str )
{
    // Extract param name
    const char *p_name_begin = util::Skip(str," (");
    const char *p_name_end = util::Find(p_name_begin," =");
    if( 0 == p_name_end ) return 0;

    // Check if its a valid token
    char param_name[LanguageDef::cMaxParamNameLength];
    unsigned int length = p_name_end - p_name_begin;
    if( util::IsLetter(*p_name_begin) )
    {
        if( length < LanguageDef::cMaxParamNameLength )
        {
            strncpy( param_name, p_name_begin, length );
            param_name[length] = '\0';
            //std::cout << "Parsed param name '" << param_name << "'" << std::endl;
        }
        else
        {
            std::cout << "[ERROR] Param name too long > "
                      << LanguageDef::cMaxParamNameLength << std::endl;
            return 0;
        }
    }
    else
        return 0;

    // Extract value and store with proper name and type
    const char *p_value_begin = util::Skip(p_name_end," =");
    const char *p_value_end = 0;
    switch( *p_value_begin )
    {
    case '[': p_value_end = ParseAndStoreArrayNumber(param_name,p_value_begin); break;
    case '"': p_value_end = ParseAndStoreString(param_name,p_value_begin); break;
    case '$': break; //p_value_end = ParseAndStoreVariable(param_name,p_value_begin) break;
    default:
        if( util::IsDigit(*p_value_begin) || util::IsNumericSymbol(*p_value_begin) )
            p_value_end = ParseAndStoreNumber(param_name,p_value_begin);
        break;
    }
    return p_value_end;
}

const char *Interpreter::ParseNumber( const char *str, double &value )
{
    // Extract number substring
    const char *p_begin = util::Skip(str," ");
    const char *p_end = util::Skip(p_begin,"0123456789+-.");
    char number_str[LanguageDef::cMaxNumberDigits];
    unsigned int length = p_end - p_begin;
    if( length < LanguageDef::cMaxNumberDigits )
    {
        strncpy( number_str, p_begin, length );
        number_str[length] = '\0';
        //std::cout << "Parsed number '" << number_str << "'" << std::endl;
        value = atof(number_str);
        return p_end;
    }
    else
    {
        value = 0;
        std::cout << "[ERROR] Number too long > "
                  << LanguageDef::cMaxNumberDigits << std::endl;
        return 0;
    }
}

const char *Interpreter::ParseAndStoreNumber( const char *param_name, const char *str )
{
    double value;
    const char *p_end = ParseNumber(str,value);
    m_Context.m_Params.Add(param_name,value);
    return p_end;
}

const char *Interpreter::ParseAndStoreArrayNumber( const char *param_name, const char *str )
{
    double value_array[LanguageDef::cMaxArrayLength];
    unsigned int num_values = 0;

    //std::cout << "ParseAndStoreArrayNumber: " << str << std::endl;

    // Skip array opening
    const char *p_begin = util::Skip(str," [");
    while( 0 != p_begin && *p_begin != ']' )
    {
        if( num_values < LanguageDef::cMaxArrayLength )
        {
            const char *p_end = ParseNumber(p_begin,value_array[num_values]);
            if( 0 == p_end )
            {
                std::cout << "[ERROR] begin = " << p_begin << std::endl;
                return 0;
            }
            p_begin = util::Skip(p_end,", ");
        }
        num_values++;
    }

    if( num_values >= LanguageDef::cMaxArrayLength )
        std::cout << "[ERROR] Array too long! " << p_begin << std::endl;

    // Store array
    m_Context.m_Params.AddArray(param_name,value_array,num_values);

    // Jump last ']' if found
    if( 0 != p_begin ) p_begin++;
    return p_begin;
}

const char *Interpreter::ParseAndStoreString( const char *param_name, const char *str )
{
    // Extract substring
    const char *p_begin = util::Skip(str,'"');
    if( 0 == p_begin ) return 0;
    const char *p_end = util::Find(p_begin,"\"");
    if( 0 == p_end ) return 0;
    char tmp_str[LanguageDef::cMaxStringLength];
    unsigned int length = p_end - p_begin;
    if( length < LanguageDef::cMaxStringLength )
    {
        strncpy( tmp_str, p_begin, length );
        tmp_str[length] = '\0';
        //std::cout << "Parsed string '" << tmp_str << "'" << std::endl;
        m_Context.m_Params.AddString(param_name,tmp_str);
        return p_end+1; //jump closing '"'
    }
    else
    {
        std::cout << "[ERROR] String too long > " << LanguageDef::cMaxStringLength << std::endl;
        return 0;
    }
}

bool Interpreter::Execute( const LanguageDef::MethodDef &method,
                           VarSeq &user_context,
                           const VarSeq &params,
                           VarSeq &results )
{
    //std::cout << "Executing Method: " << method.GetName() << std::endl;
    //Match params with paramsdef
    m_Context.m_MatchedParams.Match( method.GetParamsDef(), params );
    //std::cout << "Params = " << std::endl << params.GetIS() << std::endl;
    //std::cout << "MatchedParams = " << std::endl << m_Context.m_MatchedParams.GetIS() << std::endl;
    //Execute
    bool bResult = method.GetMethodPtr()( user_context, results, m_Context.m_MatchedParams );
    //std::cout << "Results = " << std::endl << results.GetIS() << std::endl;
    //Return result
    return bResult;
}

#ifdef __ENABLE_TEST
class TestLanguage: public LanguageDef
{
public:
    TestLanguage()
    {
        double value_array[4] = {0.1, 0.2, 0.3, 0.4};
        DefineMethod("testmethod","tests method definition",&Method_Test)
            .Add("arg1",int32(0))
            .Add("arg2",float(1))
            .AddArray("arg3",value_array,4);

        DefineMethod("particle3d","Creates a Particle 3D",&Method_CreateParticle)
            .Add("mass",double(1))
            .Add("radius",double(0.1))
            .AddArray("pos",value_array,3)
            .AddArray("vel",value_array,3);
    }
    ~TestLanguage() {}

private:
    static bool Method_Test( VarSeq &user_context, VarSeq &results, const VarSeq &params )
    {
        results.Add("res",double(42.23));
        return true;
    }

    static bool Method_CreateParticle( VarSeq &user_context, VarSeq &results, const VarSeq &params )
    {
        results.Add("particle_id",int(1));
        return true;
    }
};

void Test_Script()
{
    TestLanguage test;
    Interpreter interp(test);
    char script[] = "help (); testmethod (arg1 = 255, arg3 = [ 1 0.5 -2.33 ] ); particle3d(pos=[5,3,1]);";
    bool bResult = interp.Execute(0,script,strlen(script));
    std::cout << "Interpreter returned: " << ((bResult)?"true":"false") << std::endl;
};
#endif //__ENABLE_TEST

} // namespace skr
