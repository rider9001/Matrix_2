/// ------------------------------------------
/// @file Poly.cpp
///
/// @brief Source file for polynomial related functions
/// ------------------------------------------

#include "../inc/Poly.h"

/// ------------------------------------------
void printCompressedPoly(const std::vector<double>& poly)
{
    for (size_t power = 0; power < poly.size(); power++)
    {
        if (power == 0)
        {
            std::cout << poly.at(power) << " + ";
            continue;
        }
        else if (power == 1)
        {
            std::cout << poly.at(power) << "x + ";
            continue;
        }

        std::cout << poly.at(power) << "x^" << power;
        if (power != poly.size() - 1)
        {
            std::cout << " + ";
        }
    }

    std::cout << std::endl;
}

/// ------------------------------------------
std::vector<double> CompressFactors(const std::vector< std::pair<double, double> >& factorList)
{
    // Highest polynomial rank is equal to the number of factors
    std::vector<double> compressedPoly(factorList.size() + 1);

    for (size_t curPower = 0; curPower < compressedPoly.size(); curPower++)
    {
        // First and last powers are expections, being the product of either the x or non-x factor components
        if (curPower == 0)
        {
            double product = 1;
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
            double product = 1;
            for (auto factor : factorList)
            {
                product *= factor.first;
            }

            compressedPoly.at(curPower) = product;
            continue;
        }

        double sum_product = 0;
        for (size_t curFactor = 0; curFactor < factorList.size(); curFactor++)
        {
            double product = factorList.at(curFactor).first;
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
double getValCompressedPoly(const double x, const std::vector<double>& compressedPoly)
{
    double sum_factors = 0;
    for (size_t i = 0; i < compressedPoly.size(); i++)
    {
        sum_factors += compressedPoly.at(i) * pow(x, i);
    }

    return sum_factors;
}