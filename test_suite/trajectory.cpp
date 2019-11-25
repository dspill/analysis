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
    constexpr size_t length = 32;
    constexpr char* infile = "dummy_traj.xyz";
    {
        Trajectory t(infile);
        t.advance(2 * length);
        assert(t.index() == length - 1);
    }

    {
        Trajectory t(infile);
        assert(t.is_good());
        assert(t.frames_read() == 1);
        Frame frame0 = *t;
        assert(t.index() == 0);
        assert(t.timestep() == 0.);
        assert(t.number_of_molecules() == 1);
        assert(t.is_good());
    }

    // advance()
    {
        Trajectory t(infile);
        assert(t.frames_read() == 1);

        t.advance();
        assert(t.frames_read() == 2);
        assert(t.index() == 1);

        t.advance(2);
        assert(t.frames_read() == 3);
        assert(t.index() == 3);

        while(!(t.is_null())) t.advance();
        assert(t.index() == length - 1);

        t.advance();
        assert(t.index() == length - 1);
        assert(t.is_null());

        t.advance(10);
        assert(t.index() == length - 1);
        assert(t.is_null());
    }
    {
        Trajectory t(infile);
        t.advance(2 * length);
        assert(t.index() == length - 1);
        assert(t.frames_read() == 2);
        assert(t.is_null());
    }

    // move_to();
    {
        Trajectory t(infile);
        t.move_to(5);

        assert(t.index() == 5);
        assert(t.is_good());

        t.move_to(4);
        assert(t.index() == 4);
        assert(t.is_good());

        t.move_to(0);
        assert(t.index() == 0);
        assert(t->coordinate(0) == Real3D(1,1,1));
        assert(t.is_good());

        t.move_to(2 * length);
        assert(t.index() == length - 1);
        assert(t.is_null());
    }

    // size();
    {
        Trajectory t(infile);
        t.advance(3);

        const size_t current_index = t.index();
        assert(t.size() == length);
        assert(t.is_good());
        assert(t.index() == current_index);
    }

    // go_to_last_frame();
    {
        Trajectory t(infile);

        t.move_to(length - 1);
        assert(!t.clear_ahead());
        assert(t.is_null());
        const Frame f = *t;
        t.go_to_last_frame();
        assert(t.index() == length - 1);
        assert(t.is_null());
        assert(f == *t);

        t.reset();
        t.go_to_last_frame();
        assert(f == *t);
        assert(t.index() == length - 1);
        assert(t.is_null());
    }


    // reset();
    {
        Trajectory t(infile);
        t.advance(3);
        t.reset();
        assert(t.is_good());
        assert(t.frames_read() == 1);
        assert(t.index() == 0);
    }

    // operator *
    {
        Trajectory t(infile);
        Frame frame0 = *t;
        t.advance();
        Frame frame1 = *t;
        t.advance();
        Frame frame2 = *t;

        t.reset();
        assert(*t == frame0);
        ++t;
        assert(*t == frame1);
        t += 1;
        assert(*t == frame2);

        t.reset();
        t.advance(2);
        assert(*t == frame2);
        t.move_to(1);
        assert(t.index() == 1);
        assert(*t == frame1);
        t.move_to(2);
        assert(t.index() == 2);
        assert(*t == frame2);
    }

    // operator []
    {
        Trajectory t(infile);
        Frame frame0 = *t;
        t.advance();
        Frame frame1 = *t;
        t.advance();
        Frame frame2 = *t;

        assert(t[0] == frame0);
        assert(t[1] == frame1);
        assert(t[2] == frame2);

        assert(*t == frame0);
        ++t;
        assert(t[0] == frame0);
        assert(t[1] == frame1);
        assert(t[2] == frame2);
        assert(*t == frame1);
    }

    // operator ++
    {
        size_t i = 0;
        for(Trajectory t(infile); t.is_good(); ++t)
        {
            assert(t.index() == i++);
        }
        assert(i == length);
    }

    // operator +=
    {
        constexpr size_t step = 7;
        assert(length % step != 0);
        Trajectory t(infile);
        while(t.is_good())
        {
            t += step;
        }
        assert(t.index() == length - 1);
        assert(t.frames_read() == length / step + 2); // very last frame is read too
    }

    // analysis part TODO

    return 0;
}

