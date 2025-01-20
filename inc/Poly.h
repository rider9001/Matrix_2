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

#include "Complex.h"

/// ------------------------------------------
/// @brief Prints a compressed polynomial
/// Single line with endl
///
/// @param poly compressed polynomial to print
void printCompressedPoly(const std::vector<double>& poly);

/// ------------------------------------------
/// @brief Compresses a list of factors into the minimal form
/// e.g: (x-3)(x+2) -> x^2 - x - 6
/// Returns a list of coefficents, with the index indicating the power
/// E.g for the example above: x^2 - x - 6 -> {-6, -1, 1}
///
/// @param factorList list of factors stored as pairs, (2x-3) -> <2, -3>
///
/// @return list of coefficents for the compressed polynomial
std::vector<double> CompressFactors(const std::vector< std::pair<double, double> >& factorList);

/// @brief Using a compressed polynomial coefficent list, return the output for a value of x
///
/// @param x input value
/// @param compressedPoly compressed polynomial to use as function
///
/// @return output of polynomial function for x
double getValCompressedPoly(const double x, const std::vector<double>& compressedPoly);