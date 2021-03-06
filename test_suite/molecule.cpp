#include "trajectory.hpp"
#include <ctime>
#include <algorithm>
#include <cassert>

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
            assert(mol.size() == 0);

            assert(mol.has_velocities() == false);
            assert(mol.is_null() == true);
        }

        // single atom molecule
        {
            Molecule mol{coordinates.cbegin(), coordinates.cbegin() + 1, 
                velocities.cbegin(), velocities.cbegin() + 1};
            assert(mol.size() == 1);

            assert(mol.has_velocities() == true);
            assert(mol.is_null() == false);

            assert(mol.end_to_end() == Real3D(0.));
            assert(mol.end_to_end_squared() == 0.);
        }

        Molecule mol{coordinates.cbegin(), coordinates.cend(), 
            velocities.cbegin(), velocities.cend()};

        assert(mol.size() == 3);
        assert(mol.has_velocities() == true);
        assert(mol.is_null() == false);

        assert(mol.coordinate(0) == c0);
        assert(mol.coordinate(1) == c1);
        assert(mol.coordinate(2) == c2);

        assert(mol.velocity(0) == v0);
        assert(mol.velocity(1) == v1);
        assert(mol.velocity(2) == v2);

        assert(mol.bond(0) == c1 - c0);
        assert(mol.bond(1) == c2 - c1);

        assert(mol.bond_length(0) == (c1 - c0).abs());
        assert(mol.bond_length(1) == (c2 - c1).abs());

        assert(mol.mean_bond_length() 
                == ((c1 - c0).abs() + (c2 - c1).abs())/2.);

        assert(mol.max_bond_length() 
                == max((c1 - c0).abs(), (c2 - c1).abs()));



        assert(mol.end_to_end() == c2 - c0);
        assert(mol.end_to_end_squared() == (c2 - c0).sqr());

        {
            Real3D rcm((c0 + c1 + c2) / 3.);
            assert(mol.center_of_mass() == rcm);


            double rg_sq{0.};
            for(auto c : coordinates)
            {
                rg_sq += (c - rcm).sqr();
            }
            rg_sq /= coordinates.size();
            assert(mol.radius_of_gyration_squared() == rg_sq);
        }

        assert(mol.rouse_mode_0() == mol.rouse_mode(0));
        assert((mol.center_of_mass() 
                    - mol.rouse_mode(0)/sqrt(mol.size())).abs() < 10e-12);

        vector<Real3D> coordinates2;

        coordinates2.push_back(Real3D(0.));
        coordinates2.push_back(Real3D(0.));
        coordinates2.push_back(Real3D(0.));

        Molecule mol2{coordinates2.cbegin(), coordinates2.cend()};

        assert(mol2.size() == 3);
        assert(mol2.has_velocities() == false);

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

            assert(mol2.mean_squared_displacement(mol) == msd);
        }
        cout << "===== Molecule clear =====\n";
    }
    return 0;
}

