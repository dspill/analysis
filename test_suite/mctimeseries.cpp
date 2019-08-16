#include "trajectory.hpp"
#include "mctimeseries.hpp"
#include <ctime>
#include <algorithm>
#include <cassert>
#include <boost/timer.hpp>

#define DEBUG

using namespace std;
/**
 * This program is meant to test the code for inconsistencies.
 */
int main()
{
    MCTimeseries<double> ts("dummy.dat");
    assert(ts.consistent());
    assert(ts.timestep() == 2e3);
    assert(ts.initial_time() == 2e3);

    vector<Timeseries<double>> v;
    for(size_t i = 0; i < ts.size(); ++i)
    {
        v.emplace_back(Timeseries<double>("dummy.dat", i));
        assert(v[i] == ts.timeseries(i));
    }
        cout << "===== MCTimeseries clear =====\n";
}
