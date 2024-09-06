/// ------------------------------------------
/// @file Complex.h
///
/// @brief Header file for complex number structures and associated
/// functions
/// ------------------------------------------

#include <cmath>
#include <stdexcept>

/// @brief Complex number structure, cartesian form
typedef struct Complex_C_t
{
    double real = 0;
    double imagine = 0;
};


///--------------------------------------------------------
/// @brief Overload of +, adds two complex numbers together
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator+(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +, adds a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator+(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of +, adds a complex and real number
///
/// @param lcom left hand complex
/// @param lreal right hand real
///
/// @return resulting complex number
Complex_C_t operator+(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator+=(const Complex_C_t& lcom, Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -, subtracts two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator-(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -, subtracts a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator-(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of -, subtracts real number and a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator-(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -=
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator-=(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *, multiplies two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator*(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *, multiplies a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator*(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of *, multiplies a real number and a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator*(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /, divides two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator/(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /, divides a complex by a real
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator/(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of /, divides a real number by a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator/(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of ==, are complex numbers equal?
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator==(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of ==, is a complex and a real equal?
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return equality boolean
bool operator==(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of ==, is a complex and a real equal?
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator==(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of !=, are complex numbers unequal?
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return inequality boolean
bool operator!=(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of !=, is a complex and a real unequal?
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return equality boolean
bool operator!=(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of !=, is a complex and a real unequal?
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator!=(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Find the conjugate of the input com
///
/// @param com complex number input
///
/// @return conjugate of input
Complex_C_t conjugate(const Complex_C_t& com);

/// @brief Find the absolute value of the input com
///
/// @param com complex number input
///
/// @return argument
double absolute(const Complex_C_t& com);

///--------------------------------------------------------
/// @brief Find the argument of the input com
///
/// @param com complex number input
///
/// @return argument in radians, [-pi, pi] range
/// counterclockwise relative to 1+0i
///
/// @note returns 0 if com = 0+0i
double argument(const Complex_C_t& com);

///--------------------------------------------------------
/// @brief Raises eulers number by a complex value
///
/// @param com complex to raise e by
///
/// @return complex result of e^com
Complex_C_t raiseEComplex(const Complex_C_t& com);

///--------------------------------------------------------
/// @brief Raises a complex number by another complex number
///
/// @param base complex number to raise
/// @param raise power to raise the base to
///
/// @return resulting complex number
Complex_C_t powComplex(const Complex_C_t& base, const Complex_C_t& raise);
