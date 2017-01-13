// Expects raw sequences (only uppercase 'A', 'C', 'G' or 'T') on a
// line-by-line basis


#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>


 // Global size fitting in small memory model for k-mers upto 12
static size_t const USIZE = 1 << (2 * 12);

// k-mer length
static size_t K = 4;
static size_t SIZE = 1 << (2 * K); // 4 ^ K


// Count table
static size_t count[USIZE][4] = {{0}};
// Probability table
static double prob[USIZE][4];


static size_t count_g[USIZE][4] = {{0}};
static size_t count_e[USIZE][4] = {{0}};


static double prob_g[USIZE][4];
static double prob_e[USIZE][4];



// Symbolic constants used for printing
static char const ALPHA[4] =
{
    'A',
    'C',
    'G',
    'T'
}; // ALPHA


// Map a DNA base to an index (trivial)
static inline size_t to_idx(char const ch)
{
    switch (ch)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
    } // switch
    return 0;
} // to_idx


// Build a count table. Expects strings consisting of 'A', 'C', 'G' or
// 'T' only!
static inline void accumulate(char const* const string,
                              size_t const      length)
{
    size_t key = 0;
    for (size_t i = 0; i < K; ++i)
    {
        key <<= 2;
        key |= to_idx(string[i]);
    } // for

    for (size_t i = K; i < length; ++i)
    {
        count[key % SIZE][to_idx(string[i])] += 1;
        key <<= 2;
        key |= to_idx(string[i]);
    } // for
} // accumulate


// Calculate the log2 probabilities from a count table
static inline void normalize(size_t const count[USIZE][4],
                             double       prob[USIZE][4])
{
    for (size_t i = 0; i < SIZE; ++i)
    {
        double const SUM = static_cast<double>(count[i][0] + 1) +
                           static_cast<double>(count[i][1] + 1) +
                           static_cast<double>(count[i][2] + 1) +
                           static_cast<double>(count[i][3] + 1);
        prob[i][0] = log2(static_cast<double>(count[i][0] + 1) / SUM);
        prob[i][1] = log2(static_cast<double>(count[i][1] + 1) / SUM);
        prob[i][2] = log2(static_cast<double>(count[i][2] + 1) / SUM);
        prob[i][3] = log2(static_cast<double>(count[i][3] + 1) / SUM);
    } // for
} // normalize


// Compare a probability table to a count table
static inline double compare(double const prob[USIZE][4],
                             size_t const count[USIZE][4])
{
    double sum = 0.0;

    for (size_t i = 0; i < SIZE; ++i)
    {
        double const prod = prob[i][0] * static_cast<double>(count[i][0]) +
                            prob[i][1] * static_cast<double>(count[i][1]) +
                            prob[i][2] * static_cast<double>(count[i][2]) +
                            prob[i][3] * static_cast<double>(count[i][3]);
        sum += prod;
    } // for
    return sum;
} // compare


// DEBUGGING: Do *NOT* use the print functions for large K
static inline void print_count(void)
{
#if !defined(NDEBUG)
    for (size_t i = 0; i < SIZE; ++i)
    {
        for (size_t j = 0; j < K; ++j)
        {
            fprintf(stderr, "%c", ALPHA[(i >> ((K - j - 1) * 2)) & 3]);
        } // for
        fprintf(stderr, ": %9ld %9ld %9ld %9ld\n", count[i][0],
                                                   count[i][1],
                                                   count[i][2],
                                                   count[i][3]);
    } // for
#endif
} // print_count


// DEBUGGING: Do *NOT* use the print functions for large K
static inline void print_prob(void)
{
#if !defined(NDEBUG)
    for (size_t i = 0; i < SIZE; ++i)
    {
        for (size_t j = 0; j < K; ++j)
        {
            fprintf(stderr, "%c", ALPHA[(i >> ((K - j - 1) * 2)) & 3]);
        } // for
        fprintf(stderr, ": %9f %9f %9f %9f\n", prob[i][0],
                                               prob[i][1],
                                               prob[i][2],
                                               prob[i][3]);
    } // for
#endif
} // print_prob


// Save a probability table in binary form to stdout
static inline size_t write_model(void)
{
    return fwrite(count[0], sizeof(count[0][0]), SIZE * 4, stdout);
} // write_model


// Load a probability table in binary form from stdin
static inline void read_model(void)
{
    size_t const read = fread(count, sizeof(count[0][0]), SIZE * 4, stdin);
    if (read != SIZE * 4)
    {
        fprintf(stderr, "fread failed\n");
        exit(EXIT_FAILURE);
    } // if
} // read_model


// Check whether a character is a valid base
static inline bool is_base(char const ch)
{
    return ch == 'A' ||
           ch == 'C' ||
           ch == 'G' ||
           ch == 'T';
} // is_base


// Calculate the context key from a given start position
static inline size_t calc_key(char const* const string,
                              size_t const      length,
                              size_t const      start)
{
    if (start + K >= length)
    {
        return -1;
    } // if

    size_t key = 0;
    for (size_t i = start; i < start + K; ++i)
    {
        if (!is_base(string[i]))
        {
            return -1;
        } // if
        key <<= 2;
        key |= to_idx(string[i]);
    } // for
    return key;
} // calc_key



static inline void complement(size_t count[USIZE][4])
{
    for (size_t i = 0; i < SIZE; ++i)
    {
        size_t key = 0;
        for (size_t j = 0; j < 4; ++j)
        {
            key = (~j) & 3;
            for (size_t k = 0; k < K - 1; ++k)
            {
                key <<= 2;
                key |= (~(i >> (2 * k))) & 3;
            } // for
            size_t const idx = (~(i >> ((K - 1) * 2))) & 3;

            if (key >= i)
            {
                count[i][j] += count[key][idx];
                count[key][idx] = count[i][j];
            } // if
        } // for
    } // for
} // complement


// Build the count table with Ns
static inline size_t accumulateN(char const* const string,
                                 size_t const      length)
{
    size_t i = 0;
    size_t cnt = 0;


    while (i < length)
    {
        // ignore non-bases (N)
        while (i < length && !is_base(string[i]))
        {
            i += 1;
        } // while

        // Calculate the context key from the first valid base
        size_t key = calc_key(string, length, i);
        i += K;
        // Skip invalid keys (context between Ns is too small)
        if (key == static_cast<size_t>(-1))
        {
            continue;
        } // if

        // Accumulate the counts
        while (i < length && is_base(string[i]))
        {
            count[key % SIZE][to_idx(string[i])] += 1;
            cnt += 1;
            key <<= 2;
            key |= to_idx(string[i]);
            i += 1;
        } // while
    } // while
    return cnt;
} // accumulateN


// Clear the count table
static inline void clear(void)
{
    for (size_t i = 0; i < SIZE; ++i)
    {
        count[i][0] = 0;
        count[i][1] = 0;
        count[i][2] = 0;
        count[i][3] = 0;
    } // for
} // clear


// Entry point
int main(int argc, char *argv[])
{

/* COMBINE TWO MODELS
    if (argc < 4)
    {
        fprintf(stderr, "usage: %s K\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    K = atoi(argv[1]);
    SIZE = 1 << (2 * K); // 4 ^ K

    if (K > 12)
    {
        fprintf(stderr, "K = %ld is too large\n", K);
        return EXIT_FAILURE;
    } // if

    {
    FILE* file = fopen(argv[2], "rb");
    if (file == 0)
    {
        fprintf(stderr, "fopen failed\n");
        return EXIT_FAILURE;
    } // if

    size_t const read = fread(count_g, sizeof(count_g[0][0]), SIZE * 4, file);
    fclose(file);
    if (read != SIZE * 4)
    {
        fprintf(stderr, "fread failed\n");
        return EXIT_FAILURE;
    } // if
    }

    {
    FILE* file = fopen(argv[3], "rb");
    if (file == 0)
    {
        fprintf(stderr, "fopen failed\n");
        return EXIT_FAILURE;
    } // if

    size_t const read = fread(count_e, sizeof(count_e[0][0]), SIZE * 4, file);
    fclose(file);
    if (read != SIZE * 4)
    {
        fprintf(stderr, "fread failed\n");
        return EXIT_FAILURE;
    } // if
    }

    for (size_t i = 0; i < SIZE; ++i)
    {
        count[i][0] = count_g[i][0] + count_e[i][0];
        count[i][1] = count_g[i][1] + count_e[i][1];
        count[i][2] = count_g[i][2] + count_e[i][2];
        count[i][3] = count_g[i][3] + count_e[i][3];
    } // for

    complement(count);

    print_count();

    write_model();
    return EXIT_SUCCESS;

    // lees twee modellen genome / exome
    // maak reverse complement genome /exome



    ssize_t read = 0;
    while (read != -1)
    {
        char *string = 0;
        size_t length = 0;
        // this is POSIX compliant (MS Windows is not)
        read = getline(&string, &length, stdin);
        if (read > 0)
        {
            clear();
            char *here = string; // strchr(string, '|');
            if (here != 0)
            {
                size_t const len = strlen(here + 1) - 1;
                accumulateN(here, len);
            } // if

            print_count();
            complement();
            print_count();
        } // if
        free(string);
    } // while


    return EXIT_SUCCESS;

*/ // END COMBINE TWO MODELS

// 
    if (argc < 4)
    {
        fprintf(stderr, "usage: %s K\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    K = atoi(argv[1]);
    SIZE = 1 << (2 * K); // 4 ^ K

    if (K > 12)
    {
        fprintf(stderr, "K = %ld is too large\n", K);
        return EXIT_FAILURE;
    } // if

    {
    FILE* file = fopen(argv[2], "rb");
    if (file == 0)
    {
        fprintf(stderr, "fopen failed\n");
        return EXIT_FAILURE;
    } // if

    size_t const read = fread(count_g, sizeof(count_g[0][0]), SIZE * 4, file);
    fclose(file);
    if (read != SIZE * 4)
    {
        fprintf(stderr, "fread failed\n");
        return EXIT_FAILURE;
    } // if
    }

    {
    FILE* file = fopen(argv[3], "rb");
    if (file == 0)
    {
        fprintf(stderr, "fopen failed\n");
        return EXIT_FAILURE;
    } // if

    size_t const read = fread(count_e, sizeof(count_e[0][0]), SIZE * 4, file);
    fclose(file);
    if (read != SIZE * 4)
    {
        fprintf(stderr, "fread failed\n");
        return EXIT_FAILURE;
    } // if
    }

    normalize(count_g, prob_g);
    normalize(count_e, prob_e);


    ssize_t read = 0;
    while (read != -1)
    {
        char *string = 0;
        size_t length = 0;
        // this is POSIX compliant (MS Windows is not)
        read = getline(&string, &length, stdin);
        if (read > 0)
        {
            clear();
            char *here = string; // strchr(string, '|');
            if (here != 0)
            {
                size_t const len = strlen(here + 1) - 1;
                double const norm = accumulateN(here, len);

                if (norm > 0.0)
                {
                    //normalize(count, prob);
                    double const v_g = compare(prob_g, count) / norm;
                    double const v_e = compare(prob_e, count) / norm;
                    //double const v_s = compare(prob, count) / norm;

                    if (v_e > v_g)
                    {
                        fputc('>', stdout);
                        fwrite(string, sizeof(string[0]), here - string, stdout);
                        fputc('\n', stdout);
                        fputs(here + 1, stdout);
                    } // if
                } // if
            } // if
        } // if
        free(string);
    } // while




    return EXIT_SUCCESS;



/*
// DEBUG


    if (argc < 2)
    {
        fprintf(stderr, "usage: %s K\n", argv[0]);
        return EXIT_FAILURE;
    } // if

    K = atoi(argv[1]);
    SIZE = 1 << (2 * K); // 4 ^ K

    if (K > 12)
    {
        fprintf(stderr, "K = %ld is too large\n", K);
        return EXIT_FAILURE;
    } // if

    read_model();
    print_count();

    return -2;




    ssize_t read = 0;
    while (read != -1)
    {
        char *string = 0;
        size_t length = 0;
        // this is POSIX compliant (MS Windows is not)
        read = getline(&string, &length, stdin);
        if (read > 0)
        {
            // last character of the string is \n
            accumulate(string, read - 1);
        } // if
        free(string);
    } // while
    //normalize();

    // Do *NOT* use for large K
    print_count();

    write_model();
    return EXIT_SUCCESS;
*/
} // main

