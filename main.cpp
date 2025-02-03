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

    std::vector<Complex_C_t> test_poly = {0,-16,0,0,0,4};

    cout << "Finding roots of: " << test_poly << endl;

    auto factors = FactorizePoly(test_poly);

    cout << factors << endl;

    return 0;
}