#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <map>

unsigned int factorial(unsigned int val) {
    unsigned int result = 1;
    for (unsigned int i = 0; i < val - 1; i++) {
        result *= val - i;
    }
    return result;
}

struct Vector {
    unsigned int nElements;
    std::vector<unsigned int> elements;
};

struct Row {
    unsigned int nRows;
    std::map<std::pair<Vector, Vector>, double>  row;
};

class Matrix {
    private:
        unsigned int dim;
        std::pair<std::vector<Row>, double> rows;

    public:
        Matrix(unsigned int dim) {
            this->dim = dim;
        }

        double getIndexValue(unsigned int i, unsigned int j) {

        }
 
};

#endif