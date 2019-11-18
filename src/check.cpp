#include <iostream>
#include <memory> // unique_ptr
#include "couf.hpp"

using namespace std;


int main()
{
    vector<array<double,5>> v;

    v.push_back(array<double,5>{1,2,3,4,5});
    v.push_back(array<double,5>{1,2,3,4,5});
    v.push_back(array<double,5>{1,2,3,4,5});
    v.push_back(array<double,5>{1,2,3,4,5});

    couf::write_to_file(v, "test.dat");

}

