/// ------------------------------------------
/// @file Poly.cpp
///
/// @brief Source file for polynomial related functions
/// ------------------------------------------

#include "../inc/Poly.h"

/// ------------------------------------------
std::ostream& operator<<(std::ostream& os, const std::vector<Complex_C_t>& poly)
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
std::ostream& operator<<(std::ostream& os, const std::vector< std::pair<double, Complex_C_t> >& factors)
{
    for (auto factor : factors)
    {
        os << "(";

        if (factor.first == 1)
        {
            os << "x ";
        }
        else
        {
            os << factor.first << "x ";
        }

        os << factor.second << ")";
    }

    return os;
}

/// ------------------------------------------
std::vector<Complex_C_t> CompressFactors(const std::vector< std::pair<double, Complex_C_t> >& factorList)
{
    // Highest polynomial rank is equal to the number of factors
    std::vector<Complex_C_t> compressedPoly(factorList.size() + 1);

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
Complex_C_t getValCompressedPoly(const Complex_C_t x, const std::vector<Complex_C_t>& compressedPoly)
{
    Complex_C_t sum_factors = 0;
    for (size_t i = 0; i < compressedPoly.size(); i++)
    {
        if (compressedPoly.at(i) != 0.0)
        {
            sum_factors += compressedPoly.at(i) * powComplex(x, i);
        }
    }

    //std::cout << compressedPoly << ", f(" << x << ") = " << sum_factors << std::endl;

    return sum_factors;
}

/// ------------------------------------------
std::vector< std::pair<double, Complex_C_t> > FactorizePoly(const std::vector<Complex_C_t>& compressedPoly, MAX_FACTOR_ITER max_itr_flag)
{
    size_t maxRank = compressedPoly.size() - 1;

    if (maxRank < 2)
    {
        throw std::invalid_argument("Polynomials below rank 2 have trivial solutions, and also break this algorithm,"
                                    "might implement rank 1 at some point.");
    }

    // Create a distribution circle for initial values
    const double radius = pow( compressedPoly.at(0).absolute() / compressedPoly.at(maxRank).absolute(), (1.0 / maxRank) );
    const double base_angle = (2 * M_PI) / maxRank;
    const double offset = M_PI / (2 * maxRank);

    std::cout << "Radius: " << radius << ", base angle: " << base_angle << ", offset: " << offset << std::endl;

    std::vector<Complex_C_t> nextValues(maxRank);
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

    // Prevent conjugate pairs and values on the real line
    // for (size_t i = 0; i < nextValues.size(); i++)
    // {
    //     if (nextValues.at(i).m_imagine == 0)
    //     {
    //         nextValues.at(i).m_imagine += OFFSET_VAL;
    //     }

    //     for (size_t j = 0; j < nextValues.size(); j++)
    //     {
    //         // for any starting value where the conjugate is too close to another starting point
    //         // adjust the imaginary component slightly
    //         if(j != i &&
    //            fabs(nextValues.at(i).conjugate().m_imagine - nextValues.at(j).m_imagine) < CONJUGATE_PROX_LIM &&
    //            nextValues.at(i).m_real == nextValues.at(j).m_real
    //           )
    //         {
    //             nextValues.at(j).m_imagine += OFFSET_VAL;
    //         }
    //     }
    // }

    std::cout << "Starting values: ";
    for (auto val : nextValues)
    {
        std::cout << val << " , ";
    }
    std::cout << std::endl;

    std::vector<Complex_C_t> currentValues = nextValues;
    size_t iter_count = 0;

    while (iter_count < max_itr_flag)
    {
        for (size_t i = 0; i < currentValues.size(); i++)
        {
            Complex_C_t curVal = currentValues.at(i);

            Complex_C_t sub_product = 1;
            for (size_t j = 0; j < currentValues.size(); j++)
            {
                if (j != i)
                {
                    sub_product *= curVal - currentValues.at(j);
                }
            }

            nextValues.at(i) = curVal - (getValCompressedPoly(curVal, compressedPoly) / sub_product);

            //std::cout << "(" << curVal << ") - (" << getValCompressedPoly(curVal, compressedPoly) << " / " <<  sub_product << ") = " << nextValues.at(i) << std::endl;
        }

        iter_count++;

        std::cout << "iteration " << iter_count << " done, values:" << std::endl;
        for (size_t j = 0; j < nextValues.size(); j++)
        {
            std::cout << "Root " << (j+1) << ": " << nextValues.at(j) << std::endl;
        }

        currentValues = nextValues;
    }

    std::cout << "Finished after " << iter_count << " iterations" << std::endl;

    std::vector< std::pair<double, Complex_C_t> > factors(maxRank);
    for(size_t i = 0; i < maxRank; i++)
    {
        factors.at(i) = {1, nextValues.at(i)};
    }

    return factors;
}