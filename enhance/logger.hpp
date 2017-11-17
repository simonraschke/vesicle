#pragma once

#ifndef PISPECTOR_VERSION
    #define PISPECTOR_VERSION "unknown"
#endif

// #include "program_options.hpp"
#include "singleton.hpp"
#include <mutex>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <cassert>



enum class SEPERATOR : char 
{ 
    WHITESPACE = ' ', 
    COMMA = ',', 
    SEMICOLON = ';', 
    SLASH = '/', 
    BACKSLASH = '\\', 
    NONE = '\0' ,
    NEWLINE = '\n'
};



struct SYMBOL
{
    template<SEPERATOR sep>
    static constexpr inline char get() 
    {
        return static_cast<char>(sep);
    }
};



struct Logger
  : public enhance::Singleton<Logger>
{
    friend struct enhance::Singleton<Logger>;
    
    
    
    template<SEPERATOR sep = SEPERATOR::WHITESPACE, typename ... Args>
    void write_new_line( Args&& ... args ); 
    
    
    template<SEPERATOR sep = SEPERATOR::NONE, typename ... Args>
    void write( Args&& ... args ); 
    
    
    std::_Put_time<char> wallTime();
    
    
private:
    Logger();
    ~Logger();
    
    
    template <SEPERATOR sep, typename... Args>
    constexpr inline typename std::enable_if<sep!=SEPERATOR::NONE>::type write_args(Args&&... args);
    
    
    template <SEPERATOR sep, typename... Args>
    constexpr inline typename std::enable_if<sep==SEPERATOR::NONE>::type write_args(Args&&... args); 
    
    
    Logger(const Logger&) = delete;
    Logger& operator = (const Logger&) = delete;
    
    
    std::mutex mutex { };
    std::ofstream logfile;
    
};




template<SEPERATOR sep, typename ... Args>
void Logger::write_new_line( Args&& ... args ) 
{
    assert(&getInstance());
    assert(logfile.is_open());
    std::lock_guard<std::mutex> lock(mutex);
    logfile << SYMBOL::get<SEPERATOR::NEWLINE>();
    logfile << wallTime();
    write_args<sep>(args...);
    logfile.flush();
}




template<SEPERATOR sep, typename ... Args>
void Logger::write( Args&& ... args ) 
{
    assert(&getInstance());
    assert(logfile.is_open());
    std::lock_guard<std::mutex> lock(mutex);
    write_args<sep>(args...);
    logfile.flush();
}




template <SEPERATOR sep, typename... Args>
constexpr inline typename std::enable_if<sep!=SEPERATOR::NONE>::type Logger::write_args(Args&&... args)
{
    assert(&getInstance());
    assert(logfile.is_open());
    using expander = int[];
    (void) expander {0, (void(logfile << SYMBOL::get<sep>() << std::forward<Args>(args)),0)...};
}




template <SEPERATOR sep, typename... Args>
constexpr inline typename std::enable_if<sep==SEPERATOR::NONE>::type Logger::write_args(Args&&... args) 
{
    assert(&getInstance());
    assert(logfile.is_open());
    using expander = int[];
    (void) expander {0, (void(logfile << std::forward<Args>(args)),0)...};
}
