/// ------------------------------------------
/// @file Poly.h
///
/// @brief Header file for polynomial related functions
/// ------------------------------------------
#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <utility>
#include <stdint.h>
#include <stddef.h>
#include <iostream>
#include <exception>
#include <algorithm>
#include <limits>

#include "Complex.h"

/// @brief Smallest number to allow in the starting values for either real or imaginary compontent
#define SMALLEST_ALLOWED_START_VAL 1.0E-12

/// @brief maximum number of durand-kerner iterations that may be performed
#define MAX_DK_ITERATIONS 1048576 // 2^20

/// @brief minimum differnce that must be exceeded by at least 1 root between durand-kerner iterations to continue
#define MIN_DIFF_CONV_TEST 1.0E-9

/// @brief Polynomial represented as the list of coefficents
typedef std::vector<Complex_C_t> Poly_Coeff_t;

/// @brief Polynomial factor, i.e: (2x-(4+3i)) -> <2, {4,3}>
typedef std::pair<double, Complex_C_t> Poly_factor_t;

/// ------------------------------------------
/// @brief Operator implementation for adding two polynomial coefficent sets
///
/// @param coeffL left polynomial coeff
/// @param coeffR right polynomial coeff
///
/// @return resulting polynomial as coeff list
Poly_Coeff_t operator+(const Poly_Coeff_t& coeffL, const Poly_Coeff_t& coeffR);

/// ------------------------------------------
/// @brief Operator implementation for subtracting two polynomial coefficent sets
///
/// @param coeffL left polynomial coeff
/// @param coeffR right polynomial coeff
///
/// @return resulting polynomial as coeff list
Poly_Coeff_t operator-(const Poly_Coeff_t& coeffL, const Poly_Coeff_t& coeffR);

/// ------------------------------------------
/// @brief Operator implementation for multiplying two polynomial coefficent sets
///
/// @param coeffL left polynomial coeff
/// @param coeffR right polynomial coeff
///
/// @return resulting polynomial as coeff list
Poly_Coeff_t operator*(const Poly_Coeff_t& coeffL, const Poly_Coeff_t& coeffR);

/// ------------------------------------------
/// @brief Cast to ostream for printing polynomial
///
/// @param os output stream
/// @param poly compressed polynomial to print
///
/// @return output ostream
std::ostream& operator<<(std::ostream& os, const Poly_Coeff_t& poly);

/// ------------------------------------------
/// @brief Cast to ostream for printing factors
///
/// @param os output stream
/// @param factors factors to print
///
/// @return output ostream
std::ostream& operator<<(std::ostream& os, const std::vector<Poly_factor_t>& factors);

/// ------------------------------------------
/// @brief Compresses a list of factors into the minimal form
/// e.g: (x-3)(x+2) -> x^2 - x - 6
/// Returns a list of coefficents, with the index indicating the power
/// E.g for the example above: x^2 - x - 6 -> {-6, -1, 1}
///
/// Supports complex factors: (x+2+3i) -> {1, {2,3}}
///
/// @param factorList list of factors stored as pairs, (2x-3) -> {2, -3}
///
/// @return list of coefficents for the compressed polynomial
Poly_Coeff_t CompressFactors(const std::vector<Poly_factor_t>& factorList);

/// ------------------------------------------
/// @brief Using a compressed polynomial coefficent list, return the output for a value of x
///
/// @param x input value
/// @param compressedPoly compressed polynomial to use as function
///
/// @return output of polynomial function for x
Complex_C_t getValCompressedPoly(const Complex_C_t x, const Poly_Coeff_t& compressedPoly);

/// ------------------------------------------
/// @brief Factorize the given complex polynomial into all roots
/// Will not filter non-unique roots
///
/// WARN: Works for 95% of all polynomials, but is not garenteed to converge
/// Seems to not like high (<3-5 depends on the polynomial) coefficents on the highest rank
///
/// Uses Durand-Kerner method: https://youtu.be/5JcpOj2KtWc
///
/// @param compressedPoly complex polynomial
///
/// @return factor list
std::vector< std::pair<double, Complex_C_t> > FactorizePoly(const Poly_Coeff_t& compressedPoly);