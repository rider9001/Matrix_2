/// ------------------------------------------
/// @file Poly.cpp
///
/// @brief Source file for polynomial related functions
/// ------------------------------------------

#include "../inc/Poly.h"

/// ------------------------------------------
Poly_Coeff_t operator+(const Poly_Coeff_t& coeffL, const Poly_Coeff_t& coeffR)
{
    size_t highestRank = std::max(coeffL.size(), coeffR.size());
    size_t shorterRank = std::min(coeffL.size(), coeffR.size());
    bool LHigherRank = coeffL.size() > coeffR.size();

    Poly_Coeff_t outCoeff(highestRank);
    for (size_t i = 0; i < highestRank; i++)
    {
        if (i < shorterRank)
        {
            outCoeff.at(i) = coeffL.at(i) + coeffR.at(i);
        }
        else
        {
            if (LHigherRank)
            {
                outCoeff.at(i) = coeffL.at(i);
            }
            else
            {
                 outCoeff.at(i) = coeffR.at(i);
            }
        }
    }

    return outCoeff;
}

/// ------------------------------------------
Poly_Coeff_t operator-(const Poly_Coeff_t& coeffL, const Poly_Coeff_t& coeffR)
{
    size_t highestRank = std::max(coeffL.size(), coeffR.size());
    size_t shorterRank = std::min(coeffL.size(), coeffR.size());
    bool LHigherRank = coeffL.size() > coeffR.size();

    Poly_Coeff_t outCoeff(highestRank);
    for (size_t i = 0; i < highestRank; i++)
    {
        if (i < shorterRank)
        {
            outCoeff.at(i) = coeffL.at(i) - coeffR.at(i);
        }
        else
        {
            if (LHigherRank)
            {
                outCoeff.at(i) = coeffL.at(i);
            }
            else
            {
                 outCoeff.at(i) = -coeffR.at(i);
            }
        }
    }

    return outCoeff;
}

/// ------------------------------------------
Poly_Coeff_t operator*(const Poly_Coeff_t& coeffL, const Poly_Coeff_t& coeffR)
{
    Poly_Coeff_t outCoeff(coeffL.size() + coeffR.size() - 1, 0);

    for (size_t i = 0; i < coeffL.size(); i++)
    {
        for (size_t j = 0; j < coeffR.size(); j++)
        {
            outCoeff.at(i+j) += coeffL.at(i) * coeffR.at(j);
        }
    }

    return outCoeff;
}

/// ------------------------------------------
std::ostream& operator<<(std::ostream& os, const Poly_Coeff_t& poly)
{
    for (size_t power = 0; power < poly.size(); power++)
    {
        if (poly.at(power) == 0.0)
        {
            continue;
        }

        if (power == 0)
        {
            os << poly.at(power) << " ";
            continue;
        }

        if (power == 1.0)
        {
            os << poly.at(power) << "x ";
            continue;
        }

        os << poly.at(power) << "x^" << power << " ";
    }

    return os;
}

/// ------------------------------------------
std::ostream& operator<<(std::ostream& os, const Poly_factors_t& factors)
{
    for (auto factor : factors)
    {
        os << "(";

        if (factor.first == 1)
        {
            os << "x";
        }
        else
        {
            os << factor.first << "x";
        }

        os << factor.second << ")";
    }

    return os;
}

/// ------------------------------------------
Poly_Coeff_t CompressFactors(const Poly_factors_t& factorList)
{
    // Highest polynomial rank is equal to the number of factors
    Poly_Coeff_t compressedPoly(factorList.size() + 1);

    for (size_t curPower = 0; curPower < compressedPoly.size(); curPower++)
    {
        // First and last powers are expections, being the product of either the x or non-x factor components
        if (curPower == 0)
        {
            Complex_C_t product = 1;
            for (auto factor : factorList)
            {
                product *= factor.second;
            }

            compressedPoly.at(curPower) = product;
            continue;
        }

        // Last power execption
        if (curPower == compressedPoly.size() - 1)
        {
            Complex_C_t product = 1;
            for (auto factor : factorList)
            {
                product *= factor.first;
            }

            compressedPoly.at(curPower) = product;
            continue;
        }

        Complex_C_t sum_product = 0;
        for (size_t curFactor = 0; curFactor < factorList.size(); curFactor++)
        {
            Complex_C_t product = factorList.at(curFactor).first;
            size_t x_factors_left = curPower - 1;
            size_t curIdx = (curFactor + 1 == factorList.size()) ? 0 : curFactor + 1;

            while (curIdx != curFactor)
            {
                if (x_factors_left > 0)
                {
                    product *= factorList.at(curIdx).first;
                    x_factors_left--;
                }
                else
                {
                    product *= factorList.at(curIdx).second;
                }

                curIdx++;
                if (curIdx == factorList.size()) curIdx = 0;
            }

            sum_product += product;
        }

        compressedPoly.at(curPower) = sum_product;
    }

    return compressedPoly;
}

/// ------------------------------------------
Complex_C_t getValCompressedPoly(const Complex_C_t x, const Poly_Coeff_t& compressedPoly)
{
    Complex_C_t sum_factors = 0;
    for (size_t i = 0; i < compressedPoly.size(); i++)
    {
        if (compressedPoly.at(i) != 0.0)
        {
            sum_factors += compressedPoly.at(i) * powReal(x, i);
        }
    }

    return sum_factors;
}

/// ------------------------------------------
Poly_factors_t FactorizePoly(const Poly_Coeff_t& compressedPoly)
{
    size_t maxRank = compressedPoly.size() - 1;

    if (maxRank < 2)
    {
        throw std::invalid_argument("Polynomials below rank 2 have trivial solutions, and also break this algorithm,"
                                    "might implement rank 1 at some point.");
    }

    Complex_C_t first_nonzero_coeff;
    for (auto coeff : compressedPoly)
    {
        if (coeff != 0)
        {
            first_nonzero_coeff = coeff;
            break;
        }
    }

    // Create a distribution circle for initial values
    const double radius = pow( first_nonzero_coeff.absolute() / compressedPoly.at(maxRank).absolute(), (1.0 / maxRank) );
    const double base_angle = (2 * M_PI) / maxRank;
    const double offset = M_PI / (2 * maxRank);

    // std::cout << "radius: " << radius << ", base angle: " << base_angle << ", offset: " << offset << std::endl;

    Poly_Coeff_t nextValues(maxRank);
    for (size_t i = 0; i < nextValues.size(); i++)
    {
        nextValues.at(i) = polarToCart({radius, (i * base_angle) + offset});

        // Fixes an issue with very small values breaking some math functions
        if (fabs(nextValues.at(i).m_real) < SMALLEST_ALLOWED_START_VAL)
        {
            nextValues.at(i).m_real = 0.0;
        }

        if (fabs(nextValues.at(i).m_imagine) < SMALLEST_ALLOWED_START_VAL)
        {
            nextValues.at(i).m_imagine = 0.0;
        }
    }

    // std::cout << "starting values: " << std::endl;
    // for(size_t i = 0; i < nextValues.size(); i++)
    // {
    //     std::cout << i+1 << ": " << nextValues.at(i) << std::endl;
    // }

    Poly_Coeff_t currentValues = nextValues;
    size_t iter_count = 0;
    bool all_converged = false;

    while (iter_count < MAX_DK_ITERATIONS && !all_converged)
    {
        iter_count++;

        for (size_t i = 0; i < currentValues.size(); i++)
        {
            Complex_C_t curVal = currentValues.at(i);

            Complex_C_t sub_product = 1.0;
            for (size_t j = 0; j < currentValues.size(); j++)
            {
                if (j != i)
                {
                    sub_product *= curVal - currentValues.at(j);
                }
            }

            nextValues.at(i) = curVal - (getValCompressedPoly(curVal, compressedPoly) / sub_product);
        }

        for (size_t i = 0; i < currentValues.size(); i++)
        {
            if (std::abs( currentValues.at(i).absolute() - nextValues.at(i).absolute() ) < MIN_DIFF_CONV_TEST)
            {
                all_converged = true;
            }
            else
            {
                // if any roots fail this test, continue iterations
                all_converged = false;
                break;
            }
        }

        currentValues = nextValues;
    }

    //std::cout << "Completed in " << iter_count << " iterations" << std::endl;

    std::vector< std::pair<double, Complex_C_t> > factors(maxRank);
    for(size_t i = 0; i < maxRank; i++)
    {
        factors.at(i) = {1, -nextValues.at(i)};
    }

    return factors;
}