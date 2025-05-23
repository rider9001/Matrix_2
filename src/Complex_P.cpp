/// ------------------------------------------
/// @file Complex_P.cpp
///
/// @brief Source file for polar complex number structures and associated
/// functions
/// ------------------------------------------

#include "../inc/Complex_P.h"

///--------------------------------------------------------
Complex_P_t operator+(const Complex_P_t& lcom, const Complex_P_t& rcom)
{
    Complex_C_t cartL{lcom.real(), lcom.imaginary()};
    Complex_C_t cartR{rcom.real(), rcom.imaginary()};
    Complex_C_t res = cartL + cartR;
    return Complex_P_t{res.absolute(), res.argument()};
}

///--------------------------------------------------------
Complex_P_t operator+(const Complex_P_t& lcom, const double& rreal)
{
    Complex_P_t polR{rreal};
    return lcom + polR;
}

///--------------------------------------------------------
Complex_P_t operator+(const double& lreal, const Complex_P_t& rcom)
{
    // + is commutative so use the other arragement
    return rcom + lreal;
}

///--------------------------------------------------------
void operator+=(Complex_P_t& lcom, const Complex_P_t& rcom)
{
    Complex_C_t cartL{lcom.real(), lcom.imaginary()};
    Complex_C_t cartR{rcom.real(), rcom.imaginary()};
    Complex_C_t res = cartL + cartR;
    lcom = Complex_P_t{res.absolute(), res.argument()};
}

///--------------------------------------------------------
void operator+=(Complex_P_t& lcom, const double& rreal)
{
    Complex_P_t polR{rreal};
    lcom = lcom + polR;
}

///--------------------------------------------------------
Complex_P_t operator-(const Complex_P_t& lcom, const Complex_P_t& rcom)
{
    Complex_C_t cartL{lcom.real(), lcom.imaginary()};
    Complex_C_t cartR{rcom.real(), rcom.imaginary()};
    Complex_C_t res = cartL - cartR;
    return Complex_P_t{res.absolute(), res.argument()};
}

///--------------------------------------------------------
Complex_P_t operator-(const Complex_P_t& lcom, const double& rreal)
{
    Complex_P_t polR{rreal};
    return lcom - polR;
}

///--------------------------------------------------------
Complex_P_t operator-(const double& lreal, const Complex_P_t& rcom)
{
    Complex_P_t polR{lreal};
    return polR - rcom;
}

///--------------------------------------------------------
Complex_P_t operator-(const Complex_P_t& com)
{
    return {-com.m_mag, com.m_arg};
}

///--------------------------------------------------------
void operator-=(Complex_P_t& lcom, const Complex_P_t& rcom)
{
    Complex_C_t cartL{lcom.real(), lcom.imaginary()};
    Complex_C_t cartR{rcom.real(), rcom.imaginary()};
    Complex_C_t res = cartL - cartR;
    lcom = Complex_P_t{res.absolute(), res.argument()};
}

///--------------------------------------------------------
void operator-=(Complex_P_t& lcom, const double& rreal)
{
    Complex_P_t polR{rreal};
    lcom = lcom - polR;
}

///--------------------------------------------------------
Complex_P_t operator*(const Complex_P_t& lcom, const Complex_P_t& rcom)
{
    return Complex_P_t{lcom.m_mag * rcom.m_mag, lcom.m_arg + rcom.m_arg};
}

///--------------------------------------------------------
Complex_P_t operator*(const Complex_P_t& lcom, const double& rreal)
{
    return Complex_P_t{lcom.m_mag * rreal, lcom.m_arg};
}

///--------------------------------------------------------
Complex_P_t operator*(const double& lreal, const Complex_P_t& rcom)
{
    return rcom * lreal;
}

///--------------------------------------------------------
void operator*=(Complex_P_t& lcom, const Complex_P_t& rcom)
{
    lcom = lcom * rcom;
}

///--------------------------------------------------------
void operator*=(Complex_P_t& lcom, const double& rreal)
{
    lcom = lcom * rreal;
}

///--------------------------------------------------------
Complex_P_t operator/(const Complex_P_t& lcom, const Complex_P_t& rcom)
{
    return Complex_P_t{lcom.m_mag / rcom.m_mag, lcom.m_arg - rcom.m_arg};
}

///--------------------------------------------------------
Complex_P_t operator/(const Complex_P_t& lcom, const double& rreal)
{
    return Complex_P_t{lcom.m_mag / rreal, lcom.m_arg};
}

///--------------------------------------------------------
Complex_P_t operator/(const double& lreal, const Complex_P_t& rcom)
{
    Complex_C_t cartR{rcom.real(), rcom.imaginary()};
    Complex_C_t res = lreal / cartR;
    return Complex_P_t{res.absolute(), res.argument()};
}

///--------------------------------------------------------
void operator/=(Complex_P_t& lcom, const Complex_P_t& rcom)
{
    lcom = lcom / rcom;
}

///--------------------------------------------------------
void operator/=(Complex_P_t& lcom, const double& rreal)
{
    lcom = lcom / rreal;
}

///--------------------------------------------------------
bool operator==(const Complex_P_t& lcom, const Complex_P_t& rcom)
{
    return (lcom.m_mag == rcom.m_mag) and (lcom.m_arg == rcom.m_arg);
}

///--------------------------------------------------------
bool operator==(const Complex_P_t& lcom, const double& rreal)
{
    return (lcom.m_mag == rreal) and (lcom.m_arg == 0);
}

///--------------------------------------------------------
bool operator==(const double& lreal, const Complex_P_t& rcom)
{
    return rcom == lreal;
}

///--------------------------------------------------------
bool operator!=(const Complex_P_t& lcom, const Complex_P_t& rcom)
{
    return !(lcom == rcom);
}

///--------------------------------------------------------
bool operator!=(const Complex_P_t& lcom, const double& rreal)
{
    return !(lcom == rreal);
}

///--------------------------------------------------------
bool operator!=(const double& lreal, const Complex_P_t& rcom)
{
    return !(rcom == lreal);
}

///--------------------------------------------------------
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
    os << std::to_string(fabs(com.m_mag)) << "∠ ";

    if (com.m_arg < 0)
    {
        os << "-";
    }
    else
    {
        os << "+";
    }
    os << std::to_string(fabs(com.m_arg) / M_PI) << "π";

    return os;
}

///--------------------------------------------------------
double Complex_P_t::real() const
{
    return m_mag * cos(m_arg);
}

///--------------------------------------------------------
double Complex_P_t::imaginary() const
{
    return m_mag * sin(m_arg);
}

///--------------------------------------------------------
Complex_P_t powReal(const Complex_P_t& base, const double& raise)
{
    return Complex_P_t
                    {
                        pow(base.m_mag, raise),
                        base.m_arg * raise
                    };
}