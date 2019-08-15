#include "trajectory.hpp"
#include <ctime>
#include <algorithm>
#include <boost/timer.hpp>
#include <cassert>

#define DEBUG

using namespace std;
/**
 * This program is meant to test the code for inconsistencies.
 */
int main()
{
    // TODO frame without velocities
    // empty constructor
    {
        Frame frame;

        assert(frame.box() == Real3D(0.));
        assert(frame.number_of_molecules() == 0);
        assert(frame.particles_per_molecule() == 0);
        assert(frame.has_velocities() == false);
        assert(frame.is_null());
        assert(frame.size() == 0);
    }

    // box constructor
    {
        Frame frame(Real3D{1., 2., 3.});
        assert(frame.box() == Real3D(1., 2., 3.));
        assert(frame.number_of_molecules() == 0);
        assert(frame.particles_per_molecule() == 0);
        assert(frame.has_velocities() == false);
        assert(frame.is_null());
        assert(frame.size() == 0);
    }


    // read from file
    {
        Frame frame0("dummy.xyz", 2);
        Frame frame;
        frame.read_xyz("dummy.xyz", 2);
        assert(frame0 == frame);
        assert(frame.box() == Real3D(97., 98., 99.));
        assert(frame.particles_per_molecule() == 2);
        assert(frame.number_of_molecules() == 2);
        assert(frame.has_velocities());
        assert(!frame.is_null());
        assert(frame.size() == 4);

        assert(frame.coordinate(0) == Real3D(4., 5., 6.));
        assert(frame.velocity(0) == Real3D(7., 8., 9.));

        assert(*frame.c_cbegin() == frame.coordinate(0));
        assert(*frame.v_cbegin() == frame.velocity(0));

        assert(frame.coordinate(3) == Real3D(22., 23., 24.));
        assert(frame.velocity(3) == Real3D(25., 26., 27.));

        assert(*(frame.c_cbegin() + 3) == frame.coordinate(3));
        assert(*(frame.v_cbegin() + 3) == frame.velocity(3));

        frame.set_box(Real3D(98., 99., 100.));
        assert(frame.box() == Real3D(98., 99., 100.));

        frame.set_particles_per_molecule(4);
        assert(frame.number_of_molecules() == 1);
        assert(frame.particles_per_molecule() == 4);

        frame.set_particles_per_molecule(1);
        assert(frame.number_of_molecules() == 4);
        assert(frame.particles_per_molecule() == 1);

        frame.set_particles_per_molecule(2);
        assert(frame.number_of_molecules() == 2);
        assert(frame.particles_per_molecule() == 2);

        frame.set_number_of_molecules(4);
        assert(frame.number_of_molecules() == 4);
        assert(frame.particles_per_molecule() == 1);

        frame.set_number_of_molecules(1);
        assert(frame.number_of_molecules() == 1);
        assert(frame.particles_per_molecule() == 4);

        frame.set_number_of_molecules(2);
        assert(frame.number_of_molecules() == 2);
        assert(frame.particles_per_molecule() == 2);

        // copy constructor
        Frame temp_frame{frame};
        assert(temp_frame == frame);

        // clearing frame
        frame.clear();
        assert(frame.box() == Real3D(0.));
        assert(frame.number_of_molecules() == 0);
        assert(frame.has_velocities() == false);
        assert(frame.is_null());
        assert(frame.size() == 0);
    }

    // adding and removing particles
    {
        Frame frame;

        assert(frame.particles_per_molecule() == 0);
        frame.add_particle(Real3D(28., 29., 30.), Real3D(31., 32., 33.));
        assert(frame.coordinate(0) == Real3D(28., 29., 30.));
        assert(frame.velocity(0) == Real3D(31., 32., 33.));
        assert(frame.particles_per_molecule() == 1);

        frame.add_particle(Real3D(34., 35., 36.));
        assert(frame.particles_per_molecule() == 2);
        assert(frame.coordinate(1) == Real3D(34., 35., 36.));
        assert(frame.velocity(1) == Real3D(0.));

        frame.set_particles_per_molecule(1);
        assert(frame.particles_per_molecule() == 1);
        assert(frame.number_of_molecules() == 2);

        frame.remove_particle(1);
        assert(frame.size() == 1);

        frame.remove_particle(0);
        assert(frame.size() == 0);

        frame.add_particle(Real3D(28., 29., 30.), Real3D(31., 32., 33.));
        frame.add_particle(Real3D(34., 35., 36.));

        // scaling and folding
        frame.add_particle(frame.box());
        frame.fold();
        assert(frame.coordinate(frame.size() - 1) == Real3D(0.));

        frame.set_box(Real3D(100.));
        frame.remove_particle(frame.size() - 1);
        frame.add_particle(frame.box());
        frame.fold();
        assert(frame.coordinate(frame.size() - 1) == Real3D(0.));

        size_t old_size = frame.size();
        frame.remove_particle(frame.size() - 1);
        assert(frame.size() == old_size - 1);

        Frame old_frame = frame;
        frame.scale_box(2.);
        {
            auto cit = frame.c_cbegin();
            auto ocit = old_frame.c_cbegin();
            while(cit != frame.c_cend())
            {
                assert(*cit++ == 2. * *ocit++);
            }
            assert(ocit == old_frame.c_cend());
        }

        // logical equality
        {
            Frame frame1{10.}, frame2{10.};

            frame1.add_particle(Real3D{1., 2., 3.});
            frame2.add_particle(Real3D{1., 2., 3.});

            frame1.add_particle(Real3D{2., 3., 4.});
            frame2.add_particle(Real3D{2., 3., 4.});

            assert(frame1 == frame2);
        }
        {
            frame.write_xyz("/tmp/temporary_frame.xyz");
            Frame frame2;
            frame2.read_xyz("/tmp/temporary_frame.xyz", 
                    frame.particles_per_molecule());

            assert(frame == frame2);
        }
    }

    // adding and removing molecules
    {
        Frame frame0, frame1;
        frame0.add_particle(Real3D(1., 2., 3.));
        frame0.add_particle(Real3D(2., 3., 4.));
        frame0.add_particle(Real3D(3., 4., 5.));
        frame0.add_particle(Real3D(4., 5., 6.));
        frame0.set_number_of_molecules(2);

        frame1.set_particles_per_molecule(frame0.particles_per_molecule());
        frame1.add_molecule(frame0.molecule(0), Real3D(1., 1., 1.));
        frame1.add_molecule(frame0.molecule(1), Real3D(1., 1., 1.));

        assert(frame1.size() == frame0.size());
        assert(frame1.number_of_molecules() == frame0.number_of_molecules());
        assert(frame1.coordinate(0) == Real3D(2., 3., 4.));
        assert(frame1.coordinate(1) == Real3D(3., 4., 5.));
        assert(frame1.coordinate(2) == Real3D(4., 5., 6.));
        assert(frame1.coordinate(3) == Real3D(5., 6., 7.));
    }

    // frame.multiply() TODO
    {
        Frame frame0, frame1;
        frame0.add_particle(Real3D(1., 2., 3.));
        frame0.add_particle(Real3D(2., 3., 4.));
        frame0.add_particle(Real3D(3., 4., 5.));
        frame0.add_particle(Real3D(4., 5., 6.));

        frame1 = frame0.multiply();

        assert(frame1.particles_per_molecule() == frame0.particles_per_molecule());
        assert(frame1.number_of_molecules() == 8*frame0.number_of_molecules());
        assert(frame1.size() == 8*frame0.size());
    }


    // TODO frame.slice()
    // TODO frame.divide_box()
    // TODO frame.crop_box()
    // TODO frame.mean_squared_displacement()

    /* Lattice projection */
    {
        // TODO velocity lattice
        // TODO finegraining

        size_t size = 10;
        size_t n_sites = pow(size, 3);
        double * lattice = new double[n_sites]{0.};
        double * lattice2 = new double[n_sites]{0.};

        Frame frame{(double) size};
        assert(frame.box() == Real3D((double) size));

        Real3D c0 = Real3D(0.5);
        Real3D c1 = Real3D(1.);
        Real3D c2 = Real3D{1., 2., 3.};

        frame.add_particle(c0);
        frame.add_particle(c1);
        frame.add_particle(c2);

        read_lattice(frame, lattice, size);

        assert(lattice[0 + size * (0 + size * 0) == 1./8.]);
        assert(lattice[1 + size * (0 + size * 0) == 1./8.]);
        assert(lattice[0 + size * (1 + size * 0) == 1./8.]);
        assert(lattice[1 + size * (1 + size * 0) == 1./8.]);
        assert(lattice[0 + size * (0 + size * 1) == 1./8.]);
        assert(lattice[1 + size * (0 + size * 1) == 1./8.]);
        assert(lattice[0 + size * (1 + size * 1) == 1./8.]);
        assert(lattice[1 + size * (1 + size * 1) == 1./8.]);

        assert(lattice[(int)c1[2] + size * ((int)c1[1] + size * (int)c1[1]) == 1.]);
        assert(lattice[(int)c2[2] + size * ((int)c2[1] + size * (int)c2[2]) == 1.]);

        // read / write
        couf::write_3d_array_to_file(lattice, "/tmp/assert_lattice.dat", size);
        read_lattice("/tmp/assert_lattice.dat", lattice2, size);

        for(size_t i = 0; i < n_sites; ++i)
            assert(lattice[i] == lattice2[i]); 


        frame.write_lattice("/tmp/assert_lattice.dat", size);
        read_lattice("/tmp/assert_lattice.dat", lattice2, size);

        for(size_t i = 0; i < n_sites; ++i)
        {
            assert(lattice[i] == lattice2[i]); 
        }

        delete [] lattice;
        delete [] lattice2;
    }

    /* rotations */
    {
        Frame frame1("dummy.xyz", 2);
        Frame frame2("dummy.xyz", 2);

        const double thr = 10e-12;
        double angle = 1.234;
        frame1.rotate_x(angle);
        frame2.rotate(Real3D{1., 0., 0.}, angle);

        //cout << "rotated:\n";
        for(size_t i = 0; i < frame1.size(); ++i)
        {
            assert((frame1.coordinate(i) - frame2.coordinate(i)).abs()
                    < thr);
        }

        angle = 4.321;
        frame1.rotate_y(angle);
        frame2.rotate(Real3D{0., 1., 0.}, angle);
        for(size_t i = 0; i < frame1.size(); ++i)
        {
            assert((frame1.coordinate(i) - frame2.coordinate(i)).abs()
                    < thr);
        }

        angle = 2.341;
        frame1.rotate_z(angle);
        frame2.rotate(Real3D{0., 0., 1.}, angle);
        for(size_t i = 0; i < frame1.size(); ++i)
        {
            assert((frame1.coordinate(i) - frame2.coordinate(i)).abs()
                    < thr);
        }
    }

    cout << "===== Frame clear =====\n";
    return 0;
}

