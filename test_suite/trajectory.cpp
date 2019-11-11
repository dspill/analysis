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
    Trajectory traj("dummy_traj.xyz");
    couf::my_assert(traj.frames_read() == 1);

    Frame frame0 = *traj;

    couf::my_assert(traj.is_null() == false);
    couf::my_assert(traj.index() == 0);
    couf::my_assert(traj.timestep() == 0.);
    couf::my_assert(traj.number_of_molecules() == 1);

    traj.advance();
    couf::my_assert(traj.frames_read() == 2);
    couf::my_assert(traj.index() == 1);

    size_t current_index = traj.index();
    couf::my_assert(traj.size() == 32); 
    couf::my_assert(traj.index() == current_index);

    traj.reset();
    couf::my_assert(traj.frames_read() == 1);
    couf::my_assert(traj.index() == 0);

    traj.advance();
    couf::my_assert(traj.frames_read() == 2);
    Frame frame1 = *traj;
    traj.advance();
    couf::my_assert(traj.frames_read() == 3);
    Frame frame2 = *traj;

    traj.reset();
    couf::my_assert(traj.frames_read() == 1);
    couf::my_assert(*traj == frame0);
    ++traj;
    couf::my_assert(traj.frames_read() == 2);
    couf::my_assert(*traj == frame1);
    traj += 1;
    couf::my_assert(traj.frames_read() == 3);
    couf::my_assert(*traj == frame2);

    traj.reset();
    couf::my_assert(traj.frames_read() == 1);
    traj.advance(2);
    couf::my_assert(traj.frames_read() == 2);
    couf::my_assert(*traj == frame2);
    traj.move_to(1);
    couf::my_assert(traj.index() == 1);
    couf::my_assert(traj.frames_read() == 3);
    couf::my_assert(*traj == frame1);
    traj.move_to(2);
    couf::my_assert(traj.index() == 2);
    couf::my_assert(traj.frames_read() == 4);
    couf::my_assert(*traj == frame2);

    traj.reset();
    couf::my_assert(traj[0] == frame0);
    couf::my_assert(traj[1] == frame1);
    couf::my_assert(traj[2] == frame2);
    couf::my_assert(*traj == frame0);
    ++traj;
    couf::my_assert(traj[0] == frame0);
    couf::my_assert(traj[1] == frame1);
    couf::my_assert(traj[2] == frame2);
    couf::my_assert(*traj == frame1);

    traj.reset();
    while(!traj.is_null())
    {
        frame1 = *++traj;
    }


    // analysis part TODO

    return 0;
}

