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
    {
        Trajectory traj("dummy_traj.xyz", 2);

        Timeseries<double> ts1 = 
            traj.timeseries_mean(&Molecule::end_to_end_squared);
        Timeseries<double> ts2 = 
            traj.timeseries_mean(&Molecule::end_to_end_squared);

        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1 == ts2);
        couf::my_assert(!(ts1 != ts2));

        ts2 *= 2;
        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1*2 == ts2);

        for(size_t i = 0; i < ts1.size(); ++i)
        {
            couf::my_assert(2*ts1[i] == ts2[i]);
        }

        ts2 += ts1*4;
        couf::my_assert(ts1*6 == ts2);
        couf::my_assert(!(ts1*5 == ts2));
        couf::my_assert(ts1.size() == ts2.size());

        ts2 -= ts1*2;
        couf::my_assert(ts1*4 == ts2);
        couf::my_assert(ts1.size() == ts2.size());

        ts2 /= 4;
        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1 == ts2);

        ts1.write("/tmp/temp_timeseries.dat");
        ts2.clear();
        ts2.read("/tmp/temp_timeseries.dat");
        couf::my_assert(ts1 == ts2);
    }


    {
        Trajectory traj("dummy_traj.xyz", 2);
        Timeseries<Real3D> ts1 = 
            traj.timeseries_mean(&Molecule::end_to_end);
        Timeseries<Real3D> ts2 = 
            traj.timeseries_mean(&Molecule::end_to_end);

        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1 == ts2);
        couf::my_assert(!(ts1 != ts2));

        ts2 *= 2;
        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1*2 == ts2);

        for(size_t i = 0; i < ts1.size(); ++i)
        {
            couf::my_assert(2*ts1[i] == ts2[i]);
        }

        ts2 += ts1*4;
        couf::my_assert(ts1*6 == ts2);
        couf::my_assert(!(ts1*5 == ts2));
        couf::my_assert(ts1.size() == ts2.size());

        ts2 -= ts1*2;
        couf::my_assert(ts1*4 == ts2);
        couf::my_assert(ts1.size() == ts2.size());

        ts2 /= 4;
        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1 == ts2);
    }

    {
        Trajectory traj("dummy_traj.xyz", 2);
        Timeseries<vector<Real3D>> ts1 = 
            traj.timeseries(&Molecule::end_to_end);
        Timeseries<vector<Real3D>> ts2 = 
            traj.timeseries(&Molecule::end_to_end);

        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1 == ts2);
        couf::my_assert(!(ts1 != ts2));


        ts2 *= 2;
        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1*2 == ts2);

        ts2 += ts1*4;
        couf::my_assert(ts1*6 == ts2);
        couf::my_assert(!(ts1*5 == ts2));
        couf::my_assert(ts1.size() == ts2.size());

        ts2 -= ts1*2;
        couf::my_assert(ts1*4 == ts2);
        couf::my_assert(ts1.size() == ts2.size());

        ts2 /= 4;
        couf::my_assert(ts1.size() == ts2.size());
        couf::my_assert(ts1 == ts2);
    }
    return 0;
}
