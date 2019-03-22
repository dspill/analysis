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
        cout << "  -n        <int>   number_of_molecules" << '\n';
        cout << "  -o        <char*> (opt) output file" << '\n';
        cout << "  --step    <int>   (opt) only evaluate every step-th frame" << '\n';
        cout << "  --dt      <double>(opt) timestep" << '\n';
        cout << "  --cf_step <int>   (opt) take steps of cf_step when evaluating." << '\n';
        cout << "  --offset  <int>   (opt) number of lines to skip at the beginning" << '\n';
        exit(1);
    }

    // infile
    const char * infile {argv[1]};
    const char * outfile = couf::parse_arguments(argc, argv, "-o");
    size_t n_molecules {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "-n")))};
    size_t step        {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--step")))};
    double timestep    {atof(couf::parse_arguments(argc, argv, "--dt"))};
    size_t max         {static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")))};
    const int offset   {atoi(couf::parse_arguments(argc, argv, "--offset"))};

    if(outfile[0] == '0') outfile = "traj_timeseries.dat";
    if(n_molecules == 0) throw runtime_error("You have to give the number of molecules.");
    if(step == 0) step = 1;
    if(timestep == 0.) timestep = 1.;
    if(offset) cout << "Skipping " << offset << " frames." << '\n';
    if(!couf::is_given(argc, argv, "--max")) max = std::numeric_limits<size_t>::max();
    /* ----------- input done ---------- */

    // open infile
    Trajectory traj{infile, n_molecules};
    traj.advance(offset);

    std::vector<std::vector<double>> data;
    while(!traj.is_null() && traj.index() <= max)
    {
        cout << "\rreading frame " << traj.index();
        cout.flush();
        std::vector<double> values;
        
        values.push_back(traj->mean(&Molecule::end_to_end_squared));
        values.push_back(traj->mean(&Molecule::radius_of_gyration_squared));

        data.push_back(values);
        traj += step;
    }
    cout << '\n';
    Timeseries<std::vector<double>> ts{data};
    ts.set_timestep(step*timestep);
    ts.set_comment("r_e^2 r_g^2");
    ts.write(outfile);

    exit(0);
}
