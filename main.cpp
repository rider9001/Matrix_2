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
    Matrix<double> test = {{3,8,7},{4,8,6},{1,28,1}};

    auto QR_set = test.qr_decompose();
    cout << "-----A------" << endl;
    cout << test << endl;

    cout << "-----Q------" << endl;
    cout << QR_set.first << endl;

    cout << "-----R------" << endl;
    cout << QR_set.second << endl;

    /*
    vec = {2,4,3};
    test.setCol(1, vec);
    col = test.getCol(1);
    cout << "Col: ";
    for (auto num : col)
    {
        cout << num << ", ";
    }
    cout << endl;

    cout << test << endl;

    cout << "----------------------" << endl;

    Vector<double> vect({2,3,4});
    Vector<double> vect2({5,6,7});

    cout << vect << endl;
    cout << vect*vect2 << endl;
    cout << vect.crossR3(vect2) << endl;

    cout << vect.magnitude() << endl;
    cout << vect2.magnitude() << endl;

    cout << vect.cosineAng(vect2) << endl;

    vect2 = vect;
    cout << vect2 << endl;
    */

    return 0;
}