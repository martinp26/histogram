/*
 * histogram
 *
 * A tool to generate histograms out of sampled data.
 *
 * It can be used like a normal unix filter, that is, it reads from stdin and
 * writes to stdout.
 *
 * There is also raw output available for two-dimensional
 * data, using the --raw8 resp. --raw16 switches.  This is useful for
 * visualization purposes.  E.g., these images can be used as a height-map
 * for povray.
 *
 * Using 'convert', this data can be converted to a
 * image-program--viewable format:
 *  $> convert -flip -depth 8 -size 1000x1000 gray:test.raw test.pgm
 *
 * Adapt as needed (depth, resolution, and file format), flip brings
 * the (0, 0) corner to the lower left corner of the resulting image.
 * For raw output the first row in the input data is used as the
 * x-axis in raw output data.
 *
 * (C) 2004, 2007 by Martin Pohlack.  Released under the GPL.
 */

#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <arpa/inet.h>

#define HELP \
"This program calculates a histogram of a sequence of number tupels \
with dimensionality d \
read from the standard input until an end of file is encountered. \
The results are printed on the standard output as lines of the format:\n\
<bin midpoint d1> <bin midpoint d2> ... <bin midpoint dn> \
<number of sample points in subinterval>\n\
(as gnuplot likes it)\n\
\n\
Usage: %s [-r] [-d <dimensions>] [-l <low bound d1> -h <high bound d1> \
-w <# of bins in d1> ... [ ... for d2 [ ... for d3]]] \n\
\n\
 -r .............. compute relative frequencies rather than absolute ones.\n\
 -d <int> ........ input data has this dimensionality\n\
 -l <double> ..... low bound for histogram\n\
 -h <double> ..... high bound for histogram\n\
 -w <int> ........ use this amount of bins for histogram\n\
 --raw8, --raw16 . Use raw output instead of ascii.  The resulting file is a\n\
                   grey scale raw image with 8 resp. 16 bit per pixel.\n\
                   This option requires -r and scales the output to 2^8 resp.\n\
                   2^16.  Furthermore this option works only with -d 2.\n\
                   The resolution of the image is determined via -w switches\n\
 -o .............. Omit the output of leading and trailing empty bins\n\
                   (Does not work together with raw modi)\n\
 -q .............. Be quiet.\n\
\n\
 You should specify as many <d, l, w> tupels as you specified dimensions.\n\
Example 1: %s -r -d 1 -l -5.0 -h 5.0 -w 10 < in.dat > out.dat\n\
Example 2: %s -d2 -l0 -h2 -w4 -l-1 -h1 -w50 < 2d_in.dat > 2d_out.dat\n"

#define MAX_DIM  100     // maximal dimensions supported
#define MAX_LINE 1024    // maximal line length in input files
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

double low[MAX_DIM];     // specified limits
double hi[MAX_DIM];      // specified limits
double size[MAX_DIM];    // specified limits

double min_ob[MAX_DIM];  // observed values
double max_ob[MAX_DIM];  // observed values

int bins[MAX_DIM];       // stores how many bins per dimension are used
int dims = 1;
int relative = 0;
int omit_outer_zero = 0;
int verbose = 1;

int raw8  = 0;
int raw16 = 0;

int min_bin_seen[MAX_DIM];  // keep track of outer limits for each dimension
int max_bin_seen[MAX_DIM];

/* Parse input options and do some sanity checks
 */
static int parseOptions(int argc, char *argv[])
{
    struct option long_options[] =
    {
        {"raw8", 0, 0, 0},
        {"raw16", 0, 0, 0},
        {0, 0, 0, 0}
    };

    int c, option_index;
    int dim = 0;

    while ((c = getopt_long(argc, argv, "l:h:w:rd:oq", long_options,
                            &option_index)) != EOF)
    {
        switch(c)
        {
        case 0:
        {
            switch (option_index)
            {
            case 0:
                raw8 = 1;
                break;
            case 1:
                raw16 = 1;
                break;
            default:
                fprintf(stderr, "Invalid long option specified!\n");
                return 1;
            }
            break;
        }
        case 'l':
            low[dim] = atof(optarg);
            break;
        case 'h':
            hi[dim] = atof(optarg);
            break;
        case 'w':
            bins[dim] = atoi(optarg);
            dim++;
            break;
        case 'd':
            dims = atoi(optarg);
            break;
        case 'r':
            relative = 1;
            break;
        case 'o':
            omit_outer_zero = 1;
            break;
        case 'q':
            verbose = 0;
            break;
        default:
            fprintf(stderr, "Invalid option specified: '%c'\n", c);
            return 1;
        }
    }

    if (argv[optind] != NULL)
    {
        fprintf(stderr, "Too many arguments!\n");
        return 1;
    }

    if (dims < 1 || dims >= MAX_DIM)
    {
        fprintf(stderr, "Wrong dimensions specified: '%d',"
                "should be between 1 and %d!\n", dims, MAX_DIM);
        return 1;
    }

    if ((raw8 || raw16) && (dims != 2 || relative == 0))
    {
        fprintf(stderr, "-raw* needs other parameters\n");
        return 1;
    }

    if (raw8 && raw16)
    {
        fprintf(stderr, "You cannot have both, raw8 and raw16, pick one!\n");
        return 1;
    }


    for (dim = 0; dim < dims; dim++)
    {
        if ((low[dim] >= hi[dim]) || (bins[dim] < 1))
        {
            fprintf(stderr, "Wrong range arguments!\n");
            return 1;
        }
    }
    return 0;
}


/* Write accessor function to multidimensional array
 */
static void set_pos(int * field, int value, int dims, int * size, int * pos)
{
    long ind, offset, dim;

    // address into a field of variable dimensionality and size
    offset = 1;
    ind = 0;
    for (dim = 0; dim < dims; dim++)
    {
        ind += pos[dim] * offset;
        offset *= size[dim];
    }
    field[ind] = value;

    // update min / max bins
    for (dim = 0; dim < dims; dim++)
    {
        if (pos[dim] < min_bin_seen[dim])
            min_bin_seen[dim] = pos[dim];
        if (pos[dim] > max_bin_seen[dim])
            max_bin_seen[dim] = pos[dim];
    }
}


/* Read accessor function to multidimensional array
 */
static int get_pos(int * field, int dims, int * size, int * pos)
{
    long ind, offset, dim;

    // address into a field of variable dimensionality and size
    offset = 1;
    ind = 0;
    for (dim = 0; dim < dims; dim++)
    {
        ind += pos[dim] * offset;
        offset *= size[dim];
    }
    return field[ind];
}


/* Advance a position, starting with the lowest dimension,
 *   return != 0 on overflow
 */
static int adv_pos(int dims, int * pos)
{
    int dim;
    for (dim = 0; dim < dims; dim++) // increase first row first
    {
        pos[dim]++;
        if (pos[dim] >= bins[dim])   // check for overflow ...
            pos[dim] = 0;            // ... and reset in case
        else
            break;
    }
    if (dim >= dims) // 'ich habe fertig'
        return 1;
    else
        return 0;
}


/* Import data from stdin
 */
static void import(int * field, int * count)
{
    int line_nr = 0, dim, in_range;
    int oor = 0;                  // count samples 'out of range'
    int pos[MAX_DIM];
    double vals[MAX_DIM];

    while (1)
    {
        char line[MAX_LINE];
        char *s, *next_s, *old_s;

        // get a line
        s = fgets(line, MAX_LINE, stdin);
        if (s == NULL)
            break; // EOF??
        line_nr++;
        if (s[0] == '#')
            continue;  // skip comment lines

        // parse the line
        next_s = line;
        old_s = line;
        for (dim = 0; dim < dims; dim++)
        {
            vals[dim] = strtod(old_s, &next_s);
            if (old_s == next_s)
            {
                fprintf(stderr, "Error parsing this line (%d): '%s'\n",
                        line_nr, line);
                fprintf(stderr, "Stopping import here ...\n");
                goto exit_import;  // no break as we leave 2 levels
            }
            old_s = next_s;
        }

        // insert data if in range
        in_range = 1;
        for (dim = 0; dim < dims; dim++)
        {
            min_ob[dim] = MIN(min_ob[dim], vals[dim]);
            max_ob[dim] = MAX(max_ob[dim], vals[dim]);
            in_range = in_range && (vals[dim] >= low[dim]) &&
                                   (vals[dim] < hi[dim]);
        }
        if (in_range)
        {
            int val;
            // convert values to bin addresses ...
            for (dim = 0; dim < dims; dim++)
            {
                pos[dim] = (int)floor(((vals[dim] - low[dim]) / size[dim]));
            }
            // ... and increment the corresponding point in the field
            val = get_pos(field, dims, bins, pos);
            set_pos(field, val + 1, dims, bins, pos);
        }
        else
            oor++;
        (*count)++;
    }
 exit_import:

    if (verbose)
    {
        fprintf(stderr, "Ranges of values read: ");
        for (dim = 0; dim < dims; dim++)
        {
            fprintf(stderr, "[%g, %g], ", min_ob[dim], max_ob[dim]);
        }
        fprintf(stderr, "\n");
    }

    if (oor > 0)
        fprintf(stderr, "Lost '%d' tupels because they were out of the"
                        " specified range'\n", oor);

    if ((*count - oor) <= 0)
    {
        fprintf(stderr, "No input data received, giving up!\n");
        exit(1);
    }
    else
        if (verbose)
        {
            fprintf(stderr, "Read '%d' tupels, '%d' were in the "
                    "specified range\n", *count, *count - oor);
        }
}


/* Output data to stdout.
 */
static void output(int * field, int count)
{
    int dim;
    int pos[MAX_DIM];
    int max = 0;             // maximal value for output scaling

    /***********************
     * prepare output
     ***********************/
    // reset pos pointer
    for (dim = 0; dim < dims; dim++)
        pos[dim] = 0;

    // compute max now
    if (raw8 || raw16)
    {
        do
        {
            max = MAX(max, get_pos(field, dims, bins, pos));
        } while (! adv_pos(dims, pos));
        if (verbose)
            fprintf(stderr, "Maximum value found: '%d'\n", max);
    }

    /***********************
     * output
     ***********************/

    // fixup min / max bins to allow one outer occurrence of 0
    for (dim = 0; dim < dims; dim++)
    {
        if (min_bin_seen[dim] > 0)
            min_bin_seen[dim]--;
        if (max_bin_seen[dim] < bins[dim] - 1)
            max_bin_seen[dim]++;
    }

    // reset pos pointer
    for (dim = 0; dim < dims; dim++)
    {
        if (omit_outer_zero)
            pos[dim] = min_bin_seen[dim];
        else
            pos[dim] = 0;
    }
    while(1)
    {
        if (raw8 || raw16)
        {
            if (raw8)
            {
                uint8_t out;

                out = (((double)get_pos(field, dims, bins, pos)) / max) * 255;
                write(1, &out, 1);
            }
            else if (raw16)
            {
                uint16_t out;

                out = (((double)get_pos(field, dims, bins, pos)) / max) * 65535;
                out = htons(out);  // convert to big endian 
                write(1, &out, 2);
            }

            if (adv_pos(dims, pos))
                break;
        }
        else
        {
            // output one line
            for (dim = 0; dim < dims; dim++)
                printf("%lf\t",
                       low[dim] + ((pos[dim] + 0.5) *
                                   (hi[dim] - low[dim])) / (double)bins[dim]);
            // now the value
            if (relative)
            {
                double bin_sizes = 1.0;
                for (dim = 0; dim < dims; dim++)
                    bin_sizes *= size[dim];
                printf("%le\n", ((double)get_pos(field, dims, bins, pos)) /
                                  (bin_sizes * count));
            }
            else
                printf("%d\n", get_pos(field, dims, bins, pos));

            // update pos and maybe print a spacer line
            for (dim = dims - 1; dim >= 0; dim--)
            {
                pos[dim]++;
                // check for position overflow ...
                if (pos[dim] >= bins[dim] ||
                    (omit_outer_zero && pos[dim] > max_bin_seen[dim]))
                {
                    if (omit_outer_zero)   // ... and reset in case
                        pos[dim] = min_bin_seen[dim];
                    else
                        pos[dim] = 0;
                    if (dim == dims - 1)
                        printf("\n");
                }
                else
                    break;
            }
            if (dim < 0) // 'ich habe fertig'
                break;
        }
    }
}


int main(int argc, char *argv[])
{
    int count = 0, dim = 0, temp;
    int * field;

    // get args
    if (parseOptions(argc, argv))
    {
        fprintf(stderr, "\n");
        fprintf(stderr, HELP, argv[0], argv[0], argv[0]);
        exit(1);
    }

    // compute bin borders etc.
    for (dim = 0; dim < dims; dim++)
    {
        size[dim] = (hi[dim] - low[dim]) / bins[dim];
    }

    if (verbose)
    {
        fprintf(stderr, "Using %d dimensions, relative = '%d' with:\n",
                dims, relative);
        for (dim = 0; dim < dims; dim++)
        {
            fprintf(stderr, "  [%lf, %lf), bin_count = '%d', bin_size ="
                    " '%lf'\n", low[dim], hi[dim], bins[dim], size[dim]);
        }
    }

    /***********************
     * create array
     ***********************/
    // compute array size for histogram buffer ...
    temp = 1;
    for (dim = 0; dim < dims; dim++)
    {
        temp *= bins[dim];
    }
    // ... and allocate it, clear it
    field = calloc(temp, sizeof(*field));

    // init. min and max statistics
    for (dim = 0; dim < dims; dim++)
    {
        min_ob[dim] = HUGE_VAL;
        max_ob[dim] = -HUGE_VAL;
        min_bin_seen[dim] = INT_MAX;
        max_bin_seen[dim] = -1;
    }

    import(field, &count);

    output(field, count);

    return 0;
}
