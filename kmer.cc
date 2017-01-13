// Expects raw sequences (only uppercase 'A', 'C', 'G' or 'T') on a
// line-by-line basis


#include <cmath>
#include <cstdlib>
#include <cstdio>


 // Global size fitting in small memory model for k-mers upto 12
static size_t const USIZE = 1 << (2 * 12);

// k-mer length
static size_t K = 4;
static size_t SIZE = 1 << (2 * K); // 4 ^ K


// Count table
static size_t count[USIZE][4] = {{0}};
// Probability table
static double prob[USIZE][4];


// Symbolic constants used for printing
static char const ALPHA[4] =
{
    'A',
    'C',
    'T',
    'G'
}; // ALPHA


// Map a DNA base to an index (bitwise magic)
static inline size_t to_idx(char const ch)
{
    // A -> 00 (0)
    // C -> 01 (1)
    // G -> 11 (3) NB: different order
    // T -> 10 (2)
    return (ch >> 1) & 3;
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
static inline void normalize(void)
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
static inline double compare(void)
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
} // compare_a


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


// Entry point
int main(int argc, char *argv[])
{
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
    //print_prob();

    write_model();
    return EXIT_SUCCESS;
} // main

