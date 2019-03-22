#include "trajectory.hpp"
#include "couf.hpp"

#define BIN_WIDTH 0.1;
using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage:\n";
        cout << argv[0] << '\n';
        cout << "Mandatory:\n";
        cout << "  <char*> infile\n";
        cout << "  --ls     <int>    lattice size\n";
        cout << "Optional:"<< endl;
        cout << "  --bw     <double> (= 0.1) bin width\n";
        cout << "  --lc     <double> (= 1.)  lattice constant (for coarse-/fine-grained lattices)\n";
        exit(1);
    }

    // infile
    const char *infile        = argv[1];
    double bin_width          = atof(couf::parse_arguments(argc, argv, "--bw"));
    const size_t lattice_size = atof(couf::parse_arguments(argc, argv, "--ls"));
    double lattice_constant   = atof(couf::parse_arguments(argc, argv, "--lc"));
    
    if(bin_width == 0.) bin_width = .1;
    if(lattice_constant == 0.) lattice_constant = 1.;

    /* format output */
    vector<vector<double>> struc_fac =
        structure_factor(infile, lattice_size, lattice_constant, bin_width);

    /* file output */
    couf::write_to_file(struc_fac, "dsf.dat");

    exit(0);
}
