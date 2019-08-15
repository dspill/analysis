#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include <algorithm>  // std::sort

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << endl;
        cout << argv[0] << endl;
        cout << "  <char*> infile" << endl;
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << endl;
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << endl;
        cout << "  --max    <int>    (opt) maximum frame to read" << endl;
        cout << "  --ppm     <int>   particles_per_molecule" << '\n';
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    size_t particles_per_molecule {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm", "0")))};

    if(particles_per_molecule == 0)
    {
        cerr << "Warning: particles_per_molecule not given\n";
    }

    char outfile[256];

    Trajectory traj(infile, particles_per_molecule);
    cout << "particles_per_molecule: " << traj->particles_per_molecule() << '\n';
    cout << "number_of_molecules: " << traj->number_of_molecules() << '\n';

    /* loop through frames */
    while(!traj.is_null())
    {
        /* file output */
        sprintf(outfile, "multiplied_frame_%05zd.xyz", traj.index());
        //traj->multiply().write_xyz(outfile);
        Frame f = *traj;
        f = f.multiply();
        f.consistent();
        f.write_xyz(outfile);

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }
    exit(0);
}

