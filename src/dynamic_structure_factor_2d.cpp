#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "fftw3.h"
#include "limits.h"
#include <algorithm>    // std::sort

#define KMAX 10.

using namespace std;

int get_linear_index(const int i, const int j, const int lattice_size, 
        const bool center = true)
{
    int inew, jnew, itemp, jtemp, index;

    if(center)
    {
        /* center the image by exchanging quadrants diagonally */
        itemp = (i + lattice_size/2) % lattice_size;
        jtemp = (j + lattice_size/2) % lattice_size;
    }
    else
    {
        itemp = i;
        jtemp = j;
    }

    if(jtemp > lattice_size/2) 
    {
        jnew = lattice_size - jtemp;

        if(itemp == 0) inew = 0;
        else inew = lattice_size - itemp;
    }
    else
    {
        inew = itemp;
        jnew = jtemp;
    }

    index = jnew + (lattice_size/2 + 1) * inew;
    return index;
}

/* function to write square matrix of fftw_complex to a file */
void write_cmatrix(const fftw_complex *lattice_transformed, const int lattice_size, const char *filename)
{
    /* write transform */
    int index;
    double re, im, abs = 0;
    ofstream outstream(filename);

    /* reversed loop order in order to obtain standard column major format */
    //for(int j = 0; j < lattice_size/2 + 1; ++j)
    for(int j = 0; j < lattice_size; ++j)
    {
        for(int i = 0; i < lattice_size; ++i)
        {
            index = get_linear_index(i, j, lattice_size);
            //index = j + lattice_size * i;

            re = lattice_transformed[index][0];
            /* caution: negative part of im must be taken depending on where we
             * are. Does not matter, if we are only interested in abs. */
            im = lattice_transformed[index][1];
            abs = re*re + im*im;

            outstream << setw(15) << sqrt(abs)/lattice_size << ' ';
        }
        outstream << endl;
    }
    outstream.close();
}

/* iterate points in 2d lattice and write them into linear array in order of
 * increasing distance to origin */
void linearize_lattice(const fftw_complex *lattice_transformed, 
        const int lattice_size, const double lattice_constant,
        vector<vector<double>> &coordinates, const double norm = 1., 
        const double bin_width = 0.01)
{
    if(!coordinates.empty())
    {
        fprintf(stderr, "ERROR: Give me an empty vector!\n");
    }

    int index, inew, jnew;
    double k, S, re, im, abs, dist;
    const double k_0 = 2. * M_PI / (lattice_size * lattice_constant);

    const double max = KMAX;
    //const double max = k_0 * lattice_size / sqrt(2.);
    const int n_bins = ceil(max / bin_width);
    int bin;

    double *histogram = (double*) calloc(n_bins, sizeof(double));
    double *abscissa  = (double*) calloc(n_bins, sizeof(double));
    int *count        = (int*)    calloc(n_bins, sizeof(int));

    /* build the histogram */
    for(int i = 0; i < lattice_size; ++i)
    {
        for(int j = 0; j < lattice_size/2 + 1; ++j)
        {
            /*  skip trivial mode */
            if(i == 0 && j == 0) continue;

            index = j + (lattice_size/2 + 1) * i;
            re = lattice_transformed[index][0];
            im = lattice_transformed[index][1];
            abs = re*re + im*im;

            if(i < lattice_size / 2) inew = i;
            else inew = (lattice_size - 1) - i;

            if(j < lattice_size / 2) jnew = j;
            else jnew = (lattice_size - 1) - j;

            dist = sqrt(inew*inew + jnew*jnew);

            k = k_0 * dist;
            S = abs / norm;

            if(k < max)
            {
                bin = (int) (k / bin_width);
                histogram[bin] += S;
                abscissa[bin]  += k;
                ++count[bin];
            }
        }
    }

    /* build the vector */
    vector<double> tuple(2);
    for(int i = 0; i < n_bins; ++i)
    {
        //tuple[0] = bin_width * i;

        if(count[i] > 0){
            // tuple[0] = abscissa[i] / count[i]; // TODO average abscissa
            tuple[0] = bin_width * (i + 0.5);
            tuple[1] = histogram[i] / count[i];
        }
        else
        {
            tuple[0] = bin_width * (i + 0.5);
            tuple[1] = 0.;
        }

        coordinates.push_back(tuple);
    }
    free(histogram);
    free(count);

    return;
}

int main(int argc, char **argv)
{
    //cout << Frame::minimal_distance(Real3D(1.,0,-1), Real3D(9.,0,4), Real3D(10)) << endl;
    //cout << Frame::fold(Real3D(10.,5.,-2.), Real3D(10)) << endl;
    //exit(1);

    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << endl;
        cout << argv[0] << endl;
        cout << "  <char*> infile" << endl;
        cout << "  --dt     <double> (opt) md timestep" << endl;
        cout << "  --bw     <double> (opt) bin width" << endl;
        cout << "  --cg     <int>    (opt) factor by which the lattice is coarse grained" << endl;
        cout << "  --fg     <int>    (opt) factor by which the lattice is fine grained" << endl;
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << endl;
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << endl;
        cout << "  --max    <int>    (opt) highest number opf frame to read" << endl;
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // frames to skip at the beginning
    const int offset    = atoi(couf::parse_arguments(argc, argv, "--offset"));
    int max       = atoi(couf::parse_arguments(argc, argv, "--max"));
    double timestep     = atof(couf::parse_arguments(argc, argv, "--dt"));
    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw"));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg"));
    int step            = atoi(couf::parse_arguments(argc, argv, "--step"));

    if(timestep == 0.0) timestep = 1;
    if(bin_width == 0.0) bin_width = 0.01;
    if(cg_factor == 0) cg_factor = 1;
    if(fg_factor == 0) fg_factor = 1;
    if(max == 0) max = std::numeric_limits<int>::max();;

    if(step == 0) step = 1;
    if(offset) cout << "Skipping " << offset << " frames." << endl;

    const double lattice_constant = (double) cg_factor / fg_factor;

    const int n_molecules = 1;
    char outfile[256];

    ifstream instream(infile);
    int n_frames = Trajectory::skip_frames(instream, offset);
    Frame frame = Trajectory::read_frame(instream, n_molecules);
    if(frame.is_null()) exit(1);

    /* get lattice size and allocate arrays */
    const int original_lattice_size = round(frame.get_box()[0]);
    const int lattice_size = round(original_lattice_size / lattice_constant);
    const int number_of_sites = lattice_size * lattice_size;
    cout << "original lattice size = " << original_lattice_size << endl;
    cout << "scaled lattice size   = " << lattice_size << endl;
    cout << "lattice constant      = " << lattice_constant << endl;
    if(lattice_size % cg_factor != 0)
    {
        fprintf(stderr, "ERROR: lattice_size must be integer multiple of cg_factor\n");
        exit(1);
    }

    /* allocate arrays */
    fftw_complex *lattice_transformed =
        fftw_alloc_complex(lattice_size * (lattice_size / 2 + 1) * sizeof(fftw_complex));
    double *lattice = fftw_alloc_real(number_of_sites * sizeof(double));
    /* generate fftw plan */
    fftw_plan plan = fftw_plan_dft_r2c_2d(lattice_size, lattice_size, lattice,
            lattice_transformed, FFTW_ESTIMATE);

    /* loop through frames */
    while(!frame.is_null() && n_frames <= max)
    {
        printf("\rReading frame %d ", n_frames);
        cout.flush();

        // reset lattice
        for(int i = 0; i < number_of_sites; ++i) lattice[i] = 0.;
        frame.read_2d_lattice(lattice, lattice_size);
        /* do the transformation */
        fftw_execute(plan);

        /* format output */
        vector<vector<double>> coordinates;
        linearize_lattice(lattice_transformed, lattice_size, lattice_constant,
                coordinates, frame.get_number_of_particles(), bin_width);

        /* file output */
        sprintf(outfile, "dsf_%05d.dat", n_frames);
        couf::write_to_file(coordinates, outfile);

        if(atoi(couf::parse_arguments(argc, argv, "--bf")))
        {
            printf("\nBrute forcing\n");
            //sprintf(outfile, "dsf_%d_brute_force.dat", n_frames);
            //couf::write_to_file(
            //), outfile
            //); 
            cout << frame.structure_factor_2d(0.25, 32);
        }

        /* advance to next frame */
        n_frames += Trajectory::skip_frames(instream, step - 1);
        if(frame.is_null() || n_frames > max) break;

        frame = Trajectory::read_frame(instream, n_molecules);
        ++n_frames;
    }
    instream.close();
    int frames_read = n_frames - offset;
    frames_read = frames_read / step + (frames_read % step != 0);

    printf("\n");
    printf("n_frames    = %d\n", n_frames);
    printf("step        = %d\n", step);
    printf("offset      = %d\n", offset);
    printf("frames_read = %d\n", frames_read);

    fftw_destroy_plan(plan);
    fftw_free(lattice);
    fftw_free(lattice_transformed);

    exit(0);
}
