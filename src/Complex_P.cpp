/// ------------------------------------------
/// @file Complex_P.cpp
///
/// @brief Source file for polar complex number structures and associated
/// functions
/// ------------------------------------------

#include "../inc/Complex_P.h"

std::ostream& operator<<(std::ostream& os, const Complex_P_t& com)
{
    if (com.m_mag < 0)
    {
        os << "-";
    }
    else
    {
        os << "+";
    }
    os << std::to_string(fabs(com.m_mag)) << " ∠ " << std::to_string(com.m_arg);

    return os;
}