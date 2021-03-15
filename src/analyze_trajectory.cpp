#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "limits.h"
#include "fftw3.h"
#include <algorithm>   // sort
#include <sys/stat.h>  // mkdir

#define QMAX 10.

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  --ppm    <int>    (opt) particles_per_molecule" << '\n';
        cout << "  --thr    <double> (opt) threshold" << '\n';
        cout << "  --cg     <int>    (opt) factor by which the lattice is "
            "coarse grained" << '\n';
        cout << "  --fg     <int>    (opt) factor by which the lattice is "
            "fine grained" << '\n';
        cout << "  --step   <int>    (opt) only evaluate every step-th frame"
            << '\n';
        cout << "  --offset <int>    (opt) number of lines to skip at the "
            "beginning" << '\n';
        cout << "  --max    <int>    (opt) highest number of frame to read" 
            << '\n';
        cout << "  --exp    <double>\n";
        cout << "  --natural_units <bool>\n";
        exit(1);
    }

    constexpr size_t precision = 10;

    // infile
    const char *infile  = argv[1];
    // outfile
    char outfile[256];

    // frames to skip at the beginning
    const size_t particles_per_molecule = static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm", "0")));
    double threshold    = atof(couf::parse_arguments(argc, argv, "--thr", "-1"));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg", "1"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg", "1"));
    bool natural_units  = atoi(couf::parse_arguments(argc, argv, "--natural_units", "1"));
    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw"));

    Trajectory traj{infile, particles_per_molecule};
    double lattice_constant = (double) cg_factor / fg_factor;
    if(particles_per_molecule <= 0)
        cout << "\nWARNING: Degree of polymerization not given\n";

    /* ----------- input done ---------- */

    cout << "infile:    " << infile << '\n';
    cout << "threshold: " << threshold << '\n';
    cout << "cg_factor: " << cg_factor << '\n';
    cout << "fg_factor: " << fg_factor << '\n';
    cout << "bin_width: " << bin_width << '\n';

    /* get lattice size */
    const size_t original_side_length = round(traj->box()[0]);

    for(auto dirname : vector<const char*>{"./dsf", "./mnk"})
    {
        continue; // TODO
        if(couf::is_directory(dirname)) continue;
        else
        {
            if(mkdir(dirname, 0777) == -1)
                throw std::runtime_error(string("could not create directory")
                        + string(dirname));
            else
                cout << "created directory " << dirname << '\n';
        }
    }

    ofstream obs_file("observables.dat");
    obs_file << "# traj.index() <R_e^2> <R_g^2> R_MSD(t,0) <bondlen>\n";
    obs_file.precision(precision);

    const Frame frame_0 = *traj;

    /* loop through frames */
    size_t side_length;
    while(!traj.is_null())
    {
        cout << "reading frame " << traj.index() << '\n';
        /* compute observables */
        obs_file << scientific << setw(precision + 2);
        obs_file << traj.index() << ' ';
        obs_file << traj->mean(&Molecule::end_to_end_squared) << ' ';
        obs_file << traj->mean(&Molecule::radius_of_gyration_squared) << ' ';
        obs_file << traj->mean_squared_displacement(frame_0) << ' ';
        //obs_file << traj->mean(&Molecule::mean_bond_length) << ' ';
        obs_file << endl;

        traj.loop_advance(argc, argv); // TODO
        continue; // TODO
        // skip the rest

        /* compute structure factor */
        {
            cout << "dsf ";
            lattice_constant = (double) cg_factor / fg_factor;
            side_length = round(original_side_length / lattice_constant);
            vector<array<double, 2>> struc_fac =
                structure_factor(*traj, side_length, lattice_constant,
                        bin_width, traj->size());

            sprintf(outfile, "./dsf/dsf_%05zd.dat", traj.index());
            couf::write_to_file(struc_fac, outfile);
        }

        /* compute minkowski_functionals */
        vector<double> thresholds{.6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4};
        //vector<double> thresholds{1.};
        vector<double> cgs {1, 2, 4, 8};
        cout << "mnk ";
        for(const size_t& cg : cgs)
        {
            //continue; // TODO
            cout << cg << ' ';
            lattice_constant = (double) cg;
            side_length = round(original_side_length / lattice_constant);
            const size_t number_of_sites = pow(side_length, 3);
            double* const lattice = new double[number_of_sites]();
            read_lattice(*traj, lattice, side_length);
            for(const double& thr : thresholds)
            {
                double eff_thr = thr;

                sprintf(outfile,
                        "./mnk/minkowski_functionals_fg%d_cg%zd_thr%3.2f.dat",
                        1, cg, eff_thr);
                if(traj.index() == 0) 
                {
                    remove(outfile);
                }
                ofstream mnk_file(outfile, ofstream::app);

                /* compute minkowski functionals */
                array<double, 6> mfs = minkowski_functionals(
                        lattice, side_length, eff_thr, 'c', natural_units);

                /* file output */
                mnk_file << fixed;
                mnk_file << scientific << setw(6);
                mnk_file << traj.index();

                mnk_file.precision(precision);
                for(auto mf : mfs)
                {
                    mnk_file << scientific << setw(precision + 8) << mf;
                }
                mnk_file << endl;
                mnk_file.close();
            }
            delete[] lattice;
        }
        cout << '\n';

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }
    obs_file.close();
    exit(0);
}
