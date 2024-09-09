/// ------------------------------------------
/// @file main.cpp
///
/// @brief Main starting point for Matrix testing
/// ------------------------------------------

#include <iostream>

#include "inc/Matrix.h"
#include "inc/Complex.h"

using std::cout;
using std::endl;

int main()
{
    auto mat = Matrix<Complex_C_t>({
        {{1,1},{1,2}},
        {{3,2},{2,1}}
    });

    cout << mat.transpose() << endl;

    return 0;
}