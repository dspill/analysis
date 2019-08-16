#include "trajectory.hpp"
#include "couf.hpp"

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  --ppm     <int> particles per molecule" << '\n';
        cout << "  --step    <int> (opt) only evaluate every step-th frame" << '\n';
        cout << "  --offset  <int> (opt) number of lines to skip at the beginning" << '\n';
        exit(1);
    }

    // infile
    const char *infile {argv[1]};
    size_t particles_per_molecule {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm")))};
    size_t step        {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--step", "1")))};
    size_t max         {static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")))};
    const int offset   {atoi(couf::parse_arguments(argc, argv, "--offset", "0"))};

    if(particles_per_molecule == 0) throw runtime_error("You have to give the chain length.");
    if(offset) cout << "Skipping " << offset << " frames." << '\n';
    if(!couf::is_given(argc, argv, "--max"))
        max = std::numeric_limits<int>::max();
    /* ----------- input done ---------- */

    // open infile
    Trajectory traj{infile, particles_per_molecule};
    traj += offset;

    std::vector<Timeseries<std::vector<Real3D>>> vec 
        = traj.timeseries_set(
                {&Molecule::end_to_end, &Molecule::center_of_mass},
                step, max
                );

    // open outfile
    const char outfile_name[256] = "autocorrelation_function.dat";
    remove(outfile_name);
    FILE *outfile;
    outfile = fopen(outfile_name, "w");

    for(size_t span = 0; span < vec[0].size(); ++span)
    {
        fprintf(outfile, "%4zd ", span * step);
        for(Timeseries ts: vec)
        {
            double norm = ts.autocorrelation_function(0);
            fprintf(outfile, "%16.9e", ts.autocorrelation_function(span)/norm);
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);

    exit(0);
}
