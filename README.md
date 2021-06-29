This project was developed for the analysis of Molecular Dynamics trajectories produced by the simulation package ESPResSo++.
The main focus is on trajectories in the `.xyz` format which are concatenations of individual `.xyz` files.
As simulations of large systems can easily produce trajectories in the tens or even hundreds of gigabytes, efficiency becomes an important factor in processing these files.
Because of this, C++ was chosen as a programming language to develop the analysis tools which are based on a hierarchy of libraries: 

`trajectory.hpp` defines a trajectory reader which can efficiently iterate large `xyz` trajectories.
It starts reading the file at the first configuration and allows to jump ahead without explicitly reading every line into memory.

From every position a `Frame` object can be read which contains particle types, coordinates and optionally velocities.
This class is defined in the fole `frame.hpp` which also defines various functions to modify and analyze the configurations.
This includes the calculation of the dynamic structure factor by mapping the configuration onto a density lattice and performing a Fast Fourier Transform, as well as the calculation of the Minkowski functionals in three dimensions.

By the `Molecule` class introduced in `molecule.hpp`, a `Frame` can be interpreted as a set of linear chain molecules (polymers).
Because of the molecules' chain-like nature, it is enough to store references to the first and last particles.
Therefore, a `Molecule` object can operate on the data stored in a `Frame` object while needing very little additional memory for itself.
The `Molecule` class provides various functions to analyze the properties of polymer molecules.

Properties calculated from `Frame` objects and `Molecule` objects at successive points in time can be stored and processed with the `Timeseries` class.

The classes introduced above allow to easily create compact tailor-made programs for the analysis and manipulation of large `xyz` trajectories.

In the below example the average squared end-to-end vector of all molecules is calculated at successive configurations in the trajectory and written to an output file:

~~~{.cpp}
#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // set output precision
    constexpr size_t precision = 10;
    // get input file from command line
    const char *infile  = argv[1];
    // get molecule size from command line argument
    const size_t particles_per_molecule = static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm", "0")));
    // open output file
    ofstream outfile("ete_squared.dat");
    outfile.precision(precision);

    // define trajectory reader with molecule size
    Trajectory traj{infile, particles_per_molecule};
    // iterate Frames 
    while(!traj.is_null())
    {
        outfile << scientific << setw(precision + 2);
        // write step to outfile
        outfile << traj.index() << ' ';
        // calculate the average squared end-to-end vector of all molecules in 
        // the current frame and write to outfile
        outfile << traj->mean(&Molecule::end_to_end_squared) << '\n';
        // advance to next frame according to command line arguments
        traj.loop_advance(argc, argv);
    }
    outfile.close();
    exit(0);
}
~~~
