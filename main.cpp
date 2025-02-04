/// ------------------------------------------
/// @file main.cpp
///
/// @brief Main starting point for Matrix testing
/// ------------------------------------------

#include <iostream>

#include "inc/Matrix.h"
#include "inc/Complex.h"
#include "inc/Poly.h"

using std::cout;
using std::endl;

int main()
{
    // auto mat = Matrix<Complex_C_t>({
    //     {{1,1},{1,2}},
    //     {{3,2},{2,1}}
    // });

    // Complex_P_t degTest{2, M_PI};
    // cout << degTest << ", " << polarToCart(degTest) << endl;
    // Complex_P_t test2{6, M_PI / 2};
    // Complex_P_t res = degTest / test2;
    // cout << res << ", " << polarToCart(res) << endl;

    Poly_Coeff_t test_poly = {-16,2,1};

    Poly_Coeff_t test2 = {10,1,1};

    auto testout = test2 * test_poly;

    cout << "(" << test2 << ") * (" << test_poly << ") = " << testout << endl;

    // cout << "Finding roots of: " << test_poly << endl;
    // auto factors = FactorizePoly(test_poly);
    // cout << factors << endl;

    return 0;
}