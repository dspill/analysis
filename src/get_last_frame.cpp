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
        cout << "  --ppm    <int>    particles per molecule" << '\n';
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    size_t particles_per_molecule {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm")))};

    Trajectory traj(infile, particles_per_molecule);
    size_t size = traj.size();
    cout << size << " frames in total\n";
    traj.advance(size - 1);
    traj->write_xyz("last_frame.xyz");
    exit(0);
}

