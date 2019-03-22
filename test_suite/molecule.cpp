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
            couf::my_assert(mol.size() == 0);

            couf::my_assert(mol.has_velocities() == false);
            couf::my_assert(mol.is_null() == true);
        }

        // single atom molecule
        {
            Molecule mol{coordinates.cbegin(), coordinates.cbegin() + 1, 
                velocities.cbegin(), velocities.cbegin() + 1};
            couf::my_assert(mol.size() == 1);

            couf::my_assert(mol.has_velocities() == true);
            couf::my_assert(mol.is_null() == false);

            couf::my_assert(mol.end_to_end() == Real3D(0.));
            couf::my_assert(mol.end_to_end_squared() == 0.);
        }

        Molecule mol{coordinates.cbegin(), coordinates.cend(), 
            velocities.cbegin(), velocities.cend()};

        couf::my_assert(mol.size() == 3);
        couf::my_assert(mol.has_velocities() == true);
        couf::my_assert(mol.is_null() == false);

        couf::my_assert(mol.coordinate(0) == c0);
        couf::my_assert(mol.coordinate(1) == c1);
        couf::my_assert(mol.coordinate(2) == c2);

        couf::my_assert(mol.velocity(0) == v0);
        couf::my_assert(mol.velocity(1) == v1);
        couf::my_assert(mol.velocity(2) == v2);

        couf::my_assert(mol.bond(0) == c1 - c0);
        couf::my_assert(mol.bond(1) == c2 - c1);

        couf::my_assert(mol.bond_length(0) == (c1 - c0).abs());
        couf::my_assert(mol.bond_length(1) == (c2 - c1).abs());

        couf::my_assert(mol.mean_bond_length() 
                == ((c1 - c0).abs() + (c2 - c1).abs())/2.);

        couf::my_assert(mol.max_bond_length() 
                == max((c1 - c0).abs(), (c2 - c1).abs()));

        couf::my_assert(mol.end_to_end() == c2 - c0);
        couf::my_assert(mol.end_to_end_squared() == (c2 - c0).sqr());

        {
            Real3D rcm((c0 + c1 + c2) / 3.);
            couf::my_assert(mol.center_of_mass() == rcm);


            double rg_sq{0.};
            for(auto c : coordinates)
            {
                rg_sq += (c - rcm).sqr();
            }
            rg_sq /= coordinates.size();
            couf::my_assert(mol.radius_of_gyration_squared() == rg_sq);
        }

        vector<Real3D> coordinates2;

        coordinates2.push_back(Real3D(0.));
        coordinates2.push_back(Real3D(0.));
        coordinates2.push_back(Real3D(0.));

        Molecule mol2{coordinates2.cbegin(), coordinates2.cend()};

        couf::my_assert(mol2.size() == 3);
        couf::my_assert(mol2.has_velocities() == false);

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

            couf::my_assert(mol2.mean_squared_displacement(mol) == msd);
        }
        cout << "===== Molecule clear =====\n";
    }
    return 0;
}

