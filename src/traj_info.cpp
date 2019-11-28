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
        cout << "  -v print individual frame infos" << endl;
        cout << "  -n        <int>   number_of_molecules" << '\n';
        exit(1);
    }

    size_t n_molecules {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "-n", "0")))};

    // infile
    const char *infile  = argv[1];
    const int v         = atoi(couf::parse_arguments(argc, argv, "-v"));
    size_t size = 0;

    Trajectory traj(infile, n_molecules);
    if(v) 
    {
        /* loop through frames */
        while(traj.is_good())
        {
            cout << traj;
            if(n_molecules) 
                printf("largest_bond = %8.2e\n", traj->max_bond_length());
            ++traj;
            ++size;
        }
    }
    else size = traj.size();
    printf("n_frames    = %zd\n", size);

    exit(0);
}

