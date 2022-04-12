#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "./sources/Matrix.hpp"

using namespace std;
using namespace zich;

int main(){
    cout << "Please select your values to matrix:" << endl;

    Matrix mat({1},1,1);
    
    cin >> mat;

    cout << "See some examples of opartors:" << endl;
    cout << endl <<"The matrix: " << endl << mat << endl;
    cout << endl << "Operaor of Unari (-): " << endl << (mat * (-1)) << endl;
    cout << endl << "Operates multiplication in scalar: " << endl << (3*mat) << endl;
    
    return 0;
}

/*
    Help: 

    Martix for examples:
    [1 1 1 1], [1 1 1 1]
    [1 2 3 4], [5 6 7 8], [1 3 5 7]
    [1], [2], [3]

    Run:
    clang++-9 Main.cpp ./sources/Matrix.hpp ./sources/Matrix.cpp
*/
