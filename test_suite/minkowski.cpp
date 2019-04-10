#include "trajectory.hpp"
#include <ctime>
#include <algorithm>
#include <boost/timer.hpp>

#define DEBUG

using namespace std;
/**
 * This program is meant to test the code for inconsistencies.
 */
int main()
{
    {
        array<double, 6> mfs, temp;

        // N'
        Frame frame(Real3D(2));

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {0., 0., 0., 0., 0., 0.};
        couf::my_assert(mfs == temp);

        // A
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {1., 3., 3., 3., 1., 1.};
        couf::my_assert(mfs == temp);

        // B
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {2., 4., 2., 2., 0., 0.};
        couf::my_assert(mfs == temp);

        // C
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 1., 0.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {2., 6., 6., 2., 2., -2.};
        couf::my_assert(mfs == temp);

        // D
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {2., 6., 6., 6., 2., -6.};
        couf::my_assert(mfs == temp);

        // E
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{1., 1., 0.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {3., 5., 1., 1., -1., -1.};
        couf::my_assert(mfs == temp);

        // F
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {3., 7., 5., 1., 1., -3.};
        couf::my_assert(mfs == temp);

        // G
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{0., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {3., 9., 9., -3., 3., -1.};
        couf::my_assert(mfs == temp);

        // I
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 0.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {4., 4., 0., 0., 0., 0.};
        couf::my_assert(mfs == temp);

        // M
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {4., 6., 0., 0., -2., -2.};
        couf::my_assert(mfs == temp);

        // H
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {4., 6., 0., 0., -2., -2.};
        couf::my_assert(mfs == temp);

        // L
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {4., 8., 4., -4., 0., 0.};
        couf::my_assert(mfs == temp);

        // J
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 1.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {4., 8., 4., -4., 0., 0.};
        couf::my_assert(mfs == temp);

        // K
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {4., 12., 12., -12., 4., 4.};
        couf::my_assert(mfs == temp);

        // G'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {5., 9., 3., -9., -1., 3.};
        couf::my_assert(mfs == temp);

        // F'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {5., 7., -1., -5., -3., 1.};
        couf::my_assert(mfs == temp);

        // E'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {5., 5., -1., -1., -1., -1.};
        couf::my_assert(mfs == temp);

        // D'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {6., 6., -6., -6., -6., 2.};
        couf::my_assert(mfs == temp);

        // C'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {6., 6., -2., -6., -2., 2.};
        couf::my_assert(mfs == temp);

        // B'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{1., 1., 1.});
        frame.add_particle(Real3D{0., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {6., 4., -2., -2., 0., 0.};
        couf::my_assert(mfs == temp);

        // A'
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {7., 3., -3., -3., 1., 1.};
        couf::my_assert(mfs == temp);

        // N
        frame.clear();
        frame.set_box(Real3D(2));

        frame.add_particle(Real3D{0., 0., 0.});
        frame.add_particle(Real3D{1., 0., 0.});
        frame.add_particle(Real3D{0., 1., 0.});
        frame.add_particle(Real3D{1., 1., 0.});
        frame.add_particle(Real3D{0., 0., 1.});
        frame.add_particle(Real3D{1., 0., 1.});
        frame.add_particle(Real3D{0., 1., 1.});
        frame.add_particle(Real3D{1., 1., 1.});

        mfs = minkowski_functionals(frame, 2, 0.5);
        mfs /= 8;

        mfs[0] *= 8;
        mfs[1] *= 24;
        mfs[2] *= 24;
        mfs[3] *= 24;
        mfs[4] *= 8;
        mfs[5] *= 8;

        temp = {8., 0., 0., 0., 0., 0.};
        couf::my_assert(mfs == temp);
    }

    /* check additivity of Minkowski functionals. Both configurations
     * must not overlap.
     */
    {
        Frame frame1, frame2, frame3;

        frame1.set_box(Real3D(20.));
        frame2.set_box(Real3D(20.));
        frame3.set_box(Real3D(20.));

        frame1.add_particle(Real3D{1.234, 2.540, 3.229});
        frame1.add_particle(Real3D{4.295, 3.683, 8.204});
        frame1.add_particle(Real3D{9.239, 2.794, 4.204});

        frame2.add_particle(Real3D{15.293, 18.123, 16.843});
        frame2.add_particle(Real3D{18.290, 12.594, 12.329});
        frame2.add_particle(Real3D{14.943, 12.539, 11.693});

        frame3.add_particle(Real3D{1.234, 2.540, 3.229});
        frame3.add_particle(Real3D{4.295, 3.683, 8.204});
        frame3.add_particle(Real3D{9.239, 2.794, 4.204});
        frame3.add_particle(Real3D{15.293, 18.123, 16.843});
        frame3.add_particle(Real3D{18.290, 12.594, 12.329});
        frame3.add_particle(Real3D{14.943, 12.539, 11.693});

        couf::my_assert(
                minkowski_functionals(frame1, 20, 0.5) +
                minkowski_functionals(frame2, 20, 0.5) ==
                minkowski_functionals(frame3, 20, 0.5)
                );
    }

    // calculate mfs for sphere and cube
    {
        Real3D r;

        const double radius = 40.;
        size_t size = 3*radius;

        // sphere
        double thr = 0.015;
        Frame frame{Real3D((double) size)};
        frame.make_sphere(radius);
        auto mfs = minkowski_functionals(frame, size, .5, 's');

        // cube
        thr = 10e-8;
        frame.clear();
        frame.set_box(Real3D((double) size));
        frame.make_cube(radius);
        mfs = minkowski_functionals(frame, size, .5, 'c');
        couf::my_assert(abs(mfs[0] / pow(radius,3) - 1.) < thr);
        couf::my_assert(abs(mfs[1] / (6*radius*radius) - 1.) < thr);
        couf::my_assert(abs(mfs[2] / (3.*M_PI*radius) - 1.) < thr);
        couf::my_assert(abs(mfs[3] / (3.*M_PI*radius) - 1.) < thr);
        couf::my_assert(abs(mfs[4] - 1.) < thr);
        couf::my_assert(abs(mfs[5] - 1.) < thr);

        // shell output
        if(false)
        {
            frame.clear();
            frame.set_box(Real3D((double) size));
            frame.make_sphere(radius);
            mfs = minkowski_functionals(frame, size, .5, 'c');

            cout << "sphere:\n";
            cout << (abs(mfs[0] / (4./3.*M_PI*pow(radius,3)))) << '\n';
            cout << (abs(mfs[1] / (4*M_PI*radius*radius))) << '\n';
            cout << (abs(mfs[2] / (4.*M_PI*radius))) << '\n';
            cout << (abs(mfs[3] / (4.*M_PI*radius))) << '\n';
            cout << (abs(mfs[4])) << '\n';
            cout << (abs(mfs[5])) << '\n';

            frame.clear();
            frame.set_box(Real3D((double) size));
            frame.make_cube(radius);
            mfs = minkowski_functionals(frame, size, .5, 'c');

            cout << "cube:\n";
            cout << (abs(mfs[0] / pow(radius,3))) << '\n';
            cout << (abs(mfs[1] / (6*radius*radius))) << '\n';
            cout << (abs(mfs[2] / (3.*M_PI*radius))) << '\n';
            cout << (abs(mfs[3] / (3.*M_PI*radius))) << '\n';
            cout << (abs(mfs[4])) << '\n';
            cout << (abs(mfs[5])) << '\n';
        }
    }
    return 0;

}
