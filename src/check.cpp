#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // infile
    const char *infile  = argv[1];
    // frames to skip at the beginning

    char outfile[256];
    sprintf(outfile, "test_%s", infile);
    remove(outfile);

    Trajectory traj(infile, particles_per_molecule);

    /* loop through frames */
    while(!traj.is_null())
    {
        printf("\rprocessing frame %zd ", traj.index());
        cout.flush();

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }
    exit(0);
}

