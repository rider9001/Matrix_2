/// ------------------------------------------
/// @file Complex.cpp
///
/// @brief Header file for Complex_Cart_t number structure
/// ------------------------------------------
#pragma once

#define _USE_MATH_DEFINES
#include "../inc/Complex.h"

///--------------------------------------------------------
Complex_C_t operator+(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{lcom.real + rcom.real, lcom.imagine + rcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator+(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.real + rreal, lcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator+(const double& lreal, const Complex_C_t& rcom)
{
    // + is commutative so use the other arragement
    return rcom + lreal;
}

Complex_C_t operator+=(const Complex_C_t& lcom, Complex_C_t& rcom)
{
    return Complex_C_t{lcom.real + rcom.real, lcom.imagine + rcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator-(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{lcom.real - rcom.real, lcom.imagine - rcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator-(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.real - rreal, lcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator-(const double& lreal, const Complex_C_t& rcom)
{
    return Complex_C_t{lreal - rcom.real, -rcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator-=(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{lcom.real - rcom.real, lcom.imagine - rcom.imagine};
}

///--------------------------------------------------------
Complex_C_t operator*(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{
        lcom.real * rcom.real - lcom.imagine * rcom.imagine,
        lcom.real * rcom.imagine + lcom.imagine * rcom.real
    };
}

///--------------------------------------------------------
Complex_C_t operator*(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.real * rreal, lcom.imagine * rreal};
}

///--------------------------------------------------------
Complex_C_t operator*(const double& lreal, const Complex_C_t& rcom)
{
    // * is commutative so use the other arragement
    return rcom * lreal;
}

///--------------------------------------------------------
Complex_C_t operator/(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    double div = pow(rcom.real, 2) + pow(rcom.imagine, 2);
    return Complex_C_t{
        (lcom.real * rcom.real + lcom.imagine + rcom.imagine) / div,
        (lcom.imagine * rcom.real - lcom.real * rcom.imagine) / div
    };
}

///--------------------------------------------------------
Complex_C_t operator/(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.real / rreal, lcom.imagine / rreal};
}

///--------------------------------------------------------
Complex_C_t operator/(const double& lreal, const Complex_C_t& rcom)
{
    double div = pow(rcom.real, 2) + pow(rcom.imagine, 2);
    return Complex_C_t{
        (lreal * rcom.real) / div,
        (-lreal * rcom.imagine) / div
    };
}

///--------------------------------------------------------
bool operator==(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return (lcom.real == rcom.real) and (lcom.imagine == rcom.imagine);
}

///--------------------------------------------------------
bool operator==(const Complex_C_t& lcom, const double& rreal)
{
    return (lcom.real == rreal) and (lcom.imagine == 0);
}

///--------------------------------------------------------
bool operator==(const double& lreal, const Complex_C_t& rcom)
{
    return rcom == lreal;
}

///--------------------------------------------------------
bool operator!=(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return !(lcom == rcom);
}

///--------------------------------------------------------
bool operator!=(const Complex_C_t& lcom, const double& rreal)
{
    return !(lcom == rreal);
}

///--------------------------------------------------------
bool operator!=(const double& lreal, const Complex_C_t& rcom)
{
    return !(lreal == rcom);
}

///--------------------------------------------------------
Complex_C_t conjugate(const Complex_C_t& com)
{
    return Complex_C_t{com.real, -com.imagine};
}

///--------------------------------------------------------
double absolute(const Complex_C_t& com)
{
    return sqrt(pow(com.real, 2) + pow(com.imagine, 2));
}

///--------------------------------------------------------
double argument(const Complex_C_t& com)
{
    // Implemenation of atan that allows for +/- inf inputs
    // and scales output to [0, 2pi] range
    if (com.real == 0 and com.imagine == 0)
    {
        return 0;
    }

    double preRotation;
    if (com.real > 0 and com.imagine > 0)
    {
        // upper left quad
        preRotation = 0;
    }
    else if (com.real < 0 and com.imagine > 0)
    {
        // upper right quad
        preRotation = M_PI_2;
    }
    else if (com.real < 0 and com.imagine < 0)
    {
        // lower left quad
        preRotation = -M_PI_2;
    }
    else if (com.real > 0 and com.imagine < 0)
    {
        // lower left quad
        preRotation = 0;
    }

    return atan(com.imagine / com.real) + preRotation;
}

///--------------------------------------------------------
Complex_C_t raiseEComplex(const Complex_C_t& com)
{
    // e^(b+ic) = (e^b)(e^(ic)) = (e^b)((cos c) + i(sin c))
    // e^(b+ic) = e^b * cos(c) + i * e^b * sin(c)
    double eb = exp(com.real);

    return Complex_C_t{
        eb * cos(com.imagine),
        eb * sin(com.imagine)
    };
}

///--------------------------------------------------------
Complex_C_t powComplex(const Complex_C_t& base, const Complex_C_t& raise)
{
    /*
    see here for explanation: https://math.stackexchange.com/q/476998
    (a+ib) ^ (c+id) = e^( ln(r)*(c+id) + iθ*(c+id) )
    r = abs(a+ib)
    θ = arg(a+ib)

    real = ln(r)*c - d*θ
    imagine = i( d*ln(r) + cθ )
    eRaiseComplex(real + imagine)
    */

    double logAbs = log(absolute(base));
    double arg = argument(base);

    return raiseEComplex(Complex_C_t{
    logAbs * raise.real - raise.imagine * arg,
    logAbs * raise.imagine + raise.real * arg
    });
}