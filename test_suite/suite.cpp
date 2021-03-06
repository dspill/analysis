#include "trajectory.hpp"
#include <ctime>
#include <algorithm>

#define DEBUG

using namespace std;

void test(bool val)
{
    if(!val) throw runtime_error("something went wrong");
}

/**
 * This program is meant to test the code for inconsistencies.
 */
int main()
{
    /* test Molecule class */
    {
        vector<Real3D> coordinates;
        vector<Real3D> velocities;

        Real3D c0 = Real3D(0.);
        Real3D c1 = Real3D(1.);
        Real3D c2 = Real3D{1., 2., 3.};

        Real3D v0 = Real3D(0.);
        Real3D v1 = Real3D(1.);
        Real3D v2 = Real3D{2., 3., 4.};

        coordinates.push_back(c0);
        coordinates.push_back(c1);
        coordinates.push_back(c2);

        velocities.push_back(v0);
        velocities.push_back(v1);
        velocities.push_back(v2);

        // empty molecule
        {
            Molecule mol{coordinates.cbegin(), coordinates.cbegin(), 
                velocities.cbegin(), velocities.cbegin()};
            test(mol.size() == 0);

            test(mol.has_velocities() == false);
            test(mol.is_null() == true);
        }

        Molecule mol{coordinates.cbegin(), coordinates.cend(), 
            velocities.cbegin(), velocities.cend()};

        test(mol.size() == 3);
        test(mol.has_velocities() == true);
        test(mol.is_null() == false);

        test(mol.coordinate(0) == c0);
        test(mol.coordinate(1) == c1);
        test(mol.coordinate(2) == c2);

        test(mol.velocity(0) == v0);
        test(mol.velocity(1) == v1);
        test(mol.velocity(2) == v2);

        test(mol.bond(0) == c1 - c0);
        test(mol.bond(1) == c2 - c1);

        test(mol.bond_length(0) == (c1 - c0).abs());
        test(mol.bond_length(1) == (c2 - c1).abs());

        test(mol.mean_bond_length() 
                == ((c1 - c0).abs() + (c2 - c1).abs())/2.);

        test(mol.max_bond_length() 
                == max((c1 - c0).abs(), (c2 - c1).abs()));

        test(mol.end_to_end() == c2 - c0);
        test(mol.end_to_end_squared() == (c2 - c0).sqr());

        {
            Real3D rcm((c0 + c1 + c2) / 3.);
            test(mol.center_of_mass() == rcm);


            double rg_sq{0.};
            for(auto c : coordinates)
            {
                rg_sq += (c - rcm).sqr();
            }
            rg_sq /= coordinates.size();
            test(mol.radius_of_gyration_squared() == rg_sq);
        }

        vector<Real3D> coordinates2;

        coordinates2.push_back(Real3D(0.));
        coordinates2.push_back(Real3D(0.));
        coordinates2.push_back(Real3D(0.));

        Molecule mol2{coordinates2.cbegin(), coordinates2.cend()};

        test(mol2.size() == 3);
        test(mol2.has_velocities() == false);

        // mean squared displacement
        {
            double msd = 0.;
            auto it = coordinates.cbegin();
            auto it2 = coordinates2.cbegin();

            while(it != coordinates.cend())
            {
                msd += (*it - *it2).sqr();
                ++it;
                ++it2;
            }
            msd /= coordinates.size();

            test(mol2.mean_squared_displacement(mol) == msd);
        }
        cout << "===== Molecule clear =====\n";
    }


    /* test Frame class */
    {
        // TODO frame without velocities

        // empty constructor
        {
            Frame frame;

            test(frame.box() == Real3D(0.));
            test(frame.number_of_molecules() == 0);
            test(frame.particles_per_molecule() == 0);
            test(frame.has_velocities() == false);
            test(frame.is_null());
            test(frame.size() == 0);
        }

        // box constructor
        {
            Frame frame(Real3D{1., 2., 3.});
            test(frame.box() == Real3D(1., 2., 3.));
            test(frame.number_of_molecules() == 0);
            test(frame.particles_per_molecule() == 0);
            test(frame.has_velocities() == false);
            test(frame.is_null());
            test(frame.size() == 0);
        }


        // read from file
        {

            Frame frame;
            frame.read_xyz("dummy.xyz", 2);
            test(frame.box() == Real3D(97., 98., 99.));
            test(frame.particles_per_molecule() == 2);
            test(frame.number_of_molecules() == 2);
            test(frame.has_velocities());
            test(!frame.is_null());
            test(frame.size() == 4);

            test(frame.coordinate(0) == Real3D(4., 5., 6.));
            test(frame.velocity(0) == Real3D(7., 8., 9.));

            test(*frame.c_cbegin() == frame.coordinate(0));
            test(*frame.v_cbegin() == frame.velocity(0));

            test(frame.coordinate(3) == Real3D(22., 23., 24.));
            test(frame.velocity(3) == Real3D(25., 26., 27.));

            test(*(frame.c_cbegin() + 3) == frame.coordinate(3));
            test(*(frame.v_cbegin() + 3) == frame.velocity(3));

            frame.set_box(Real3D(98., 99., 100.));
            test(frame.box() == Real3D(98., 99., 100.));

            frame.set_particles_per_molecule(4);
            test(frame.number_of_molecules() == 1);
            test(frame.particles_per_molecule() == 4);

            frame.set_particles_per_molecule(1);
            test(frame.number_of_molecules() == 4);
            test(frame.particles_per_molecule() == 1);

            frame.set_particles_per_molecule(2);
            test(frame.number_of_molecules() == 2);
            test(frame.particles_per_molecule() == 2);

            frame.set_number_of_molecules(4);
            test(frame.number_of_molecules() == 4);
            test(frame.particles_per_molecule() == 1);

            frame.set_number_of_molecules(1);
            test(frame.number_of_molecules() == 1);
            test(frame.particles_per_molecule() == 4);

            frame.set_number_of_molecules(2);
            test(frame.number_of_molecules() == 2);
            test(frame.particles_per_molecule() == 2);

            // copy constructor
            Frame temp_frame{frame};
            test(temp_frame == frame);

            // clearing frame
            frame.clear();
            test(frame.box() == Real3D(0.));
            test(frame.number_of_molecules() == 0);
            test(frame.has_velocities() == false);
            test(frame.is_null());
            test(frame.size() == 0);
        }

        // adding and removing particles
        {
            Frame frame;

            test(frame.particles_per_molecule() == 0);
            frame.add_particle(Real3D(28., 29., 30.), Real3D(31., 32., 33.));
            test(frame.coordinate(0) == Real3D(28., 29., 30.));
            test(frame.velocity(0) == Real3D(31., 32., 33.));
            test(frame.particles_per_molecule() == 1);

            frame.add_particle(Real3D(34., 35., 36.));
            test(frame.particles_per_molecule() == 2);
            test(frame.coordinate(1) == Real3D(34., 35., 36.));
            test(frame.velocity(1) == Real3D(0.));

            frame.set_particles_per_molecule(1);
            test(frame.particles_per_molecule() == 1);
            test(frame.number_of_molecules() == 2);
            // TODO move from number of molecules to particles_per_molecule

            frame.remove_particle(1);
            test(frame.size() == 1);

            frame.remove_particle(0);
            test(frame.size() == 0);

            frame.add_particle(Real3D(28., 29., 30.), Real3D(31., 32., 33.));
            frame.add_particle(Real3D(34., 35., 36.));

            // scaling and folding
            frame.add_particle(frame.box());
            frame.fold();
            test(frame.coordinate(frame.size() - 1) == Real3D(0.));

            frame.set_box(Real3D(100.));
            frame.remove_particle(frame.size() - 1);
            frame.add_particle(frame.box());
            frame.fold();
            test(frame.coordinate(frame.size() - 1) == Real3D(0.));

            size_t old_size = frame.size();
            frame.remove_particle(frame.size() - 1);
            test(frame.size() == old_size - 1);

            Frame old_frame = frame;
            frame.scale_box(2.);
            {
                auto cit = frame.c_cbegin();
                auto ocit = old_frame.c_cbegin();
                while(cit != frame.c_cend())
                {
                    test(*cit++ == 2. * *ocit++);
                }
                test(ocit == old_frame.c_cend());
            }

            // logical equality
            {
                Frame frame1{10.}, frame2{10.};

                frame1.add_particle(Real3D{1., 2., 3.});
                frame2.add_particle(Real3D{1., 2., 3.});

                frame1.add_particle(Real3D{2., 3., 4.});
                frame2.add_particle(Real3D{2., 3., 4.});

                test(frame1 == frame2);
            }
            {
                frame.write_xyz("/tmp/temporary_frame.xyz");
                Frame frame2;
                frame2.read_xyz("/tmp/temporary_frame.xyz", 
                        frame.particles_per_molecule());

                test(frame == frame2);
            }
        }


        // TODO frame.slice()
        // TODO frame.multiply()
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
            test(frame.box() == Real3D((double) size));

            Real3D c0 = Real3D(0.5);
            Real3D c1 = Real3D(1.);
            Real3D c2 = Real3D{1., 2., 3.};

            frame.add_particle(c0);
            frame.add_particle(c1);
            frame.add_particle(c2);

            read_lattice(frame, lattice, size);

            test(lattice[0 + size * (0 + size * 0) == 1./8.]);
            test(lattice[1 + size * (0 + size * 0) == 1./8.]);
            test(lattice[0 + size * (1 + size * 0) == 1./8.]);
            test(lattice[1 + size * (1 + size * 0) == 1./8.]);
            test(lattice[0 + size * (0 + size * 1) == 1./8.]);
            test(lattice[1 + size * (0 + size * 1) == 1./8.]);
            test(lattice[0 + size * (1 + size * 1) == 1./8.]);
            test(lattice[1 + size * (1 + size * 1) == 1./8.]);

            test(lattice[(int)c1[2] + size * ((int)c1[1] + size * (int)c1[1]) == 1.]);
            test(lattice[(int)c2[2] + size * ((int)c2[1] + size * (int)c2[2]) == 1.]);

            // read / write
            couf::write_3d_array_to_file(lattice, "/tmp/test_lattice.dat", size);
            read_lattice("/tmp/test_lattice.dat", lattice2, size);

            for(size_t i = 0; i < n_sites; ++i)
                test(lattice[i] == lattice2[i]); 


            frame.write_lattice("/tmp/test_lattice.dat", size);
            read_lattice("/tmp/test_lattice.dat", lattice2, size);

            for(size_t i = 0; i < n_sites; ++i)
            {
                test(lattice[i] == lattice2[i]); 
            }

            delete [] lattice;
            delete [] lattice2;
        }

        /* Minkowski functionals */
        {
            {
                vector<double> mfs, temp;

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);

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
                test(mfs == temp);
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

                test(
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


                test(abs(mfs[0] / (4./3.*M_PI*pow(radius,3)) - 1.) < thr);
                test(abs(mfs[1] / (4*M_PI*radius*radius) - 1.) < thr);
                test(abs(mfs[2] / (4.*M_PI*radius) - 1.) < thr);
                test(abs(mfs[3] / (4.*M_PI*radius) - 1.) < thr);
                test(abs(mfs[4] - 1.) < thr);
                test(abs(mfs[5] - 1.) < thr);


                // cube
                thr = 10e-8;
                frame.clear();
                frame.set_box(Real3D((double) size));
                frame.make_cube(radius);
                mfs = minkowski_functionals(frame, size, .5, 'c');

                test(abs(mfs[0] / pow(radius,3) - 1.) < thr);
                test(abs(mfs[1] / (6*radius*radius) - 1.) < thr);
                test(abs(mfs[2] / (3.*M_PI*radius) - 1.) < thr);
                test(abs(mfs[3] / (3.*M_PI*radius) - 1.) < thr);
                test(abs(mfs[4] - 1.) < thr);
                test(abs(mfs[5] - 1.) < thr);

            }
        }

        cout << "===== Frame clear =====\n";
    }


    /* test Trajectory class */
    {
        Trajectory traj("dummy_traj.xyz");
        Frame frame0 = *traj;

        test(traj.is_null() == false);
        test(traj.index() == 0);
        test(traj.timestep() == 0.);
        test(traj.number_of_molecules() == 1);

        traj.advance();
        test(traj.index() == 1);

        size_t current_index = traj.index();
        test(traj.size() == 1323);
        test(traj.index() == current_index);

        traj.reset();
        test(traj.index() == 0);

        traj.advance();
        Frame frame1 = *traj;
        traj.advance();
        Frame frame2 = *traj;

        traj.reset();
        test(*traj == frame0);
        ++traj;
        test(*traj == frame1);
        traj += 1;
        test(*traj == frame2);

        traj.reset();
        traj.advance(2);
        test(*traj == frame2);
        traj.move_to(1);
        test(*traj == frame1);
        traj.move_to(2);
        test(*traj == frame2);

        traj.reset();
        test(traj[0] == frame0);
        test(traj[1] == frame1);
        test(traj[2] == frame2);
        test(*traj == frame0);
        ++traj;
        test(traj[0] == frame0);
        test(traj[1] == frame1);
        test(traj[2] == frame2);
        test(*traj == frame1);

        traj.reset();
        while(!traj.is_null())
        {
            frame1 = *++traj;
        }


        // analysis part TODO

        cout << "===== Trajectory clear =====\n";
    }

    cout << "===== ALL CLEAR =====\n";

    // TODO test timeseries
    return 0;
}

