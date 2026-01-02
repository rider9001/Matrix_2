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
    Matrix<double> test = {{1,2,3},{4,3,2},{9,1,1}};

    auto row = test.getRow(1);
    auto col = test.getCol(1);

    cout << "Row: ";
    for (auto num : row)
    {
        cout << num << ", ";
    }
    cout << endl;

    cout << "Col: ";
    for (auto num : col)
    {
        cout << num << ", ";
    }
    cout << endl;

    std::vector<double> vec = {2,1,3};
    test.setRow(1, vec);

    row = test.getRow(1);
    cout << "Set row: ";
    for (auto num : row)
    {
        cout << num << ", ";
    }
    cout << endl;

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

    return 0;
}