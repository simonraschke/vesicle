#pragma once

#include <iostream>
#include <sstream>
#include <string>



namespace enhance
{
    template<typename Bindable_type>
    std::string streamBindableToString( Bindable_type& bindable )
    {
        return static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << bindable).str();
    }

    template<typename Bindable_type>
    const char* streamBindableToChar( Bindable_type& bindable )
    {
        return static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << bindable).str().c_str();
    }
}
