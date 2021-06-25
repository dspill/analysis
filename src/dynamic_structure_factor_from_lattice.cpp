#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "fftw3.h"
#include "limits.h"
#include <algorithm>    // std::sort
#include <sys/stat.h>   // mkdir

#define QMAX 10.

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  --sl     <int>          lattice side length" << '\n';
        cout << "  --of     <char*>  (opt) outfile name" << '\n';
        cout << "  --bw     <double> (opt) bin width" << '\n';
        cout << "  --dim    <int>    (2)   spacial dimension\n";
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // outfile
    char outfile[256];

    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw", "0"));
    int original_side_length = atof(couf::parse_arguments(argc, argv, "--ls", "0"));
    size_t dim          = static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--dim", "2")));

    const double lattice_constant = 1.;
    /* ----------- input done ---------- */

    /* get lattice size */
    const size_t side_length = round(original_side_length / lattice_constant);
    const double norm = 1.; 

    cout << "infile                = " << infile << '\n';
    cout << "bin width             = " << bin_width << '\n';
    cout << "spacial dimension     = " << dim << '\n';
    cout << "original lattice size = " << original_side_length << '\n';
    cout << "scaled lattice size   = " << side_length << '\n';
    cout << "lattice constant      = " << lattice_constant << '\n';
    if (original_side_length < 1)
        throw runtime_error("lattice side length must be given");

    /* compute structure factor */
    vector<array<double, 2>> struc_fac =
        structure_factor(infile, side_length, lattice_constant, bin_width, norm, dim);

    /* file output */
    if(couf::is_given(argc, argv, "--of"))
        strcpy(outfile, couf::parse_arguments(argc, argv, "--of"));
    else
    {
        sprintf(outfile, "dsf_%s", infile);
    }
    couf::write_to_file(struc_fac, outfile);

    exit(0);
}
