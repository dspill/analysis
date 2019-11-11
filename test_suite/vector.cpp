#include "trajectory.hpp"
#include <ctime>
#include <algorithm>

#define DEBUG

using namespace std;
/**
 * This program is meant to test the code for inconsistencies.
 */
int main()
{
    vector<Real3D> a = {
        Real3D{1,2,3},
        Real3D{2,3,4}
    };
    vector<Real3D> b = {
        Real3D{2,4,6},
        Real3D{4,6,8}
    };
    couf::my_assert(a*2 == b);

    b += a;
    couf::my_assert(a*3 == b);

    b -= a;
    couf::my_assert(a*2 == b);

    b /= 2;
    couf::my_assert(a == b);

    b *= 6;
    couf::my_assert(a*6 == b);
    return 0;
}

