//#include <google/sparse_hash_map>
//#include <google/dense_hash_map>
//#include "MurmurHash3.cpp"
#include <iostream>
#include <climits>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <time.h>
#include "metrohash64.cpp"
#include <stdint.h>
#include <unordered_map>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <string.h>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <stack>
#include <limits.h>
#include <map>
#include <bitset>
#include <ctime>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include <random>
#include <cinttypes>
//#include "dna_test.h"
#include "ntHashIterator.hpp"

#include<chrono>

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"

using spp::sparse_hash_map;
typedef sparse_hash_map<uint64_t, uint32_t> SMap;

using namespace std;

using namespace chrono;

//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)        // The C header file kseq.h is a small library for parsing the FASTA/FASTQ format.
                            // For ordinary file I/O, you can use KSEQ_INIT(gzFile, gzread) to set the type of
                            // file handler and the read() function.
                            // FMI: http://lh3lh3.users.sourceforge.net/parsefastq.shtml


unsigned trailing_zeros(unsigned n) {
    return n ? __builtin_ctz(n) : -1;
}


unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}

void print_help()
{

    cout << "KmerEst [options] -f <fasta/fastq> -k <k-mer length>  -s <sample size> -o <output file>"    << endl
    << "  -h               help"                                   << endl
    << "  -f <file>       Input sequence file "                << endl
    << "  -k <k-mer size >        kmer size (default 31) "        << endl
    << "  -s <sample size>        sample size (default 25m)"        << endl
     << "  -c coverage>       coverage (default 64)"        << endl
    << "  -o         	  Prefix of the Output file " << endl;

    exit(0);
}

void parse_input(int argc, char **argv, int &k, int &maxSampleCount, int &coverage, string &inpFile, string &outpFile)
{
    for (int c = 1; c < argc; c++)  // parse command-line input
    {
        if(!strcmp(argv[c], "-h"))      // help prompt
            print_help();
        else if(!strcmp(argv[c], "-k"))
        {
            k = atoi(argv[c+1]);        // k-mer length
            c++;
        }
        else if(!strcmp(argv[c], "-f"))
        {
            inpFile = argv[c+1];              // input file
            c++;
        }
        else if(!strcmp(argv[c], "-s"))
        {
            maxSampleCount = atoi(argv[c+1]);        // sample size
            c++;
        }
        else if(!strcmp(argv[c], "-c"))
        {
            coverage = atoi(argv[c+1]);      // coverage
            c++;
        }
        else if(!strcmp(argv[c], "-o"))
        {
            outpFile = argv[c+1];           // output file
            c++;
        }
    }
}

void process_sequence(char *s, int &th, uint64_t &no_kmers, int &count, int k, int maxSampleCount, vector<SMap> &MAP)
{
    uint64_t hash = 0;

    ntHashIterator itr(s, 1, k); // ntHash iterator to iterate over the read sequence and provide
                                                // ntHash values for each of its k-mers of length n; initialized with
                                                // the first length-n window on the sequence

    while(itr != itr.end())                            // iterate until the last k-mer window
    {
        hash = (*itr)[0];                               // get the ntHash value
        ++no_kmers;                                     // one more k-mer read
        uint8_t tz = trailing_zeros(hash);              // #trailing_zeroes of this k-mer

        if(tz >= th)                                    // if #trailing_zeroes is greater than or equal to threshold
        {                                               // then sample this k-mer
            if(MAP[tz].find(hash) != MAP[tz].end())     // k-mer already present in hash map
                MAP[tz][hash] += 1;                     // increment k-mer count
            else                                        // k-mer absent in hash map
            {
                MAP[tz].insert(make_pair(hash, 1));     // insert k-mer into hash map
                ++count;                                // increment #samples_present by one

                //cout << "\r" << "count: " << count << flush;// << endl;

                /*
                    I guess there is a potential bug present in this conditional check. The conditional should be
                    (count >= k) instead of (count == k).

                    Assume that you have reached the max sample count 'k', and now have dropped the s'th hash map
                    (MAP[th] here). Hence, the number of samples present 'count' drop down to
                    (count - MAP[th].size()). Suppose that the hash map MAP[th] did not have any entries present
                    prior to deletion. Hence, 'count' retains the same value 'k' as earlier. Now at the next
                    iteration, suppose that you encounter a new k-mer, insert it to one of the hash maps, increment
                    'count' by 1; so it goes to k + 1. Only now the code flow will go to the conditional.

                    You will never drop any hash map anymore in the lifetime of this execution, and space-usage
                    will not be optimized as theorized.

                    However, my proposal of the conditional (count >= k) also should present very unlikely yet
                    potential bug(s) too, in case of a stream of hash values with very long suffixes of trailing
                    zeroes, and the lower order entries in the hash map arrays being all empty.

                    Best solution is to use a loop of the following format:
                        while(count == k):
                            drop hash maps and update count
                */
                if(count == maxSampleCount)                              // max sample count reached
                {                                           // one hash map will be dropped now
                    cout << "Samples count reached " << count << endl;

                    int cnt = MAP[th].size();               // size of the hash map to be dropped

                    cout << "Dropping a hash map of size " << cnt << endl;

                    count = count - cnt;                    // #samples_present after dropping corresponding hash map
                    SMap().swap(MAP[th]);                   // drop hash map from memory
                    //MAP[th].clear(); //MAP[th].resize(0);

                    ++th;                                   // increment threshold (s parameter)
                    cout  << "New samples count: " << count << endl;
                }
            }
        }

        ++itr;      // go over to the next window
    }
}

int main(int argc, char** argv)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    if(argc == 1)
    {
        cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
        exit(0);
    }

    int k = 31;                     // default k-mer length
    int maxSampleCount = 25000000;  // default sample size
    int coverage = 64;                   // default coverage (maximum k-mer frequency we are interested in)
    string inpFile = "", outpFile = "";       // input and output FASTA file names

    parse_input(argc, argv, k, maxSampleCount, coverage, inpFile, outpFile);

    if (inpFile.empty()  || outpFile.empty())  // empty file(s) mentioned
        print_help();

    FILE *inpFilePtr;
    kseq_t *seq;
    /*
        The C header file kseq.h is a small library for parsing the FASTA/FASTQ format.
        Function kseq_read() reads one sequence and fills the kseq_t struct which is:

            typedef struct {
                size_t l, m;
                char *s;
            } kstring_t;

            typedef struct {
                kstring_t name, comment, seq, qual;
                int last_char;
                kstream_t *f;
            } kseq_t;

        // FMI: http://lh3lh3.users.sourceforge.net/parsefastq.shtml
    */

    inpFilePtr = fopen(inpFile.c_str(), "r");     // file pointer for input FASTA file
    if(inpFilePtr == Z_NULL){
      cout << "File: " << inpFile << " does not exist" << endl;
      exit(1);
    }


    /*
        The C header file kseq.h is a small library for parsing the FASTA/FASTQ format.
        Function kseq_init() is used to initialize the parser and kseq_destroy() to destroy it.
        Function kseq_read() reads one sequence and fills the kseq_t struct.
        FMI: http://lh3lh3.users.sourceforge.net/parsefastq.shtml
    */
    seq = kseq_init(fileno(inpFilePtr));    // seq is the FASTA input parser

    /*
        sparse_hash_map is distinguished from other hash-map implementations by its stingy use of memory and by the
        ability to save and restore contents to disk. On the other hand, this hash-map implementation, while still
        efficient, is slower than other hash-map implementations.
        FMI: http://goog-sparsehash.sourceforge.net/doc/sparse_hash_map.html
    */
    vector<SMap> MAP(64);           // array of hash maps for sampled k-mers

    cout << "read the Sequences .. " << endl;

    int th = 0;                     // sample-size adaptation parameter 's';
                                    // (the threshold count of the trailing zeroes for hash values)
    uint64_t total = 0;             // count of sequences read
    uint64_t no_kmers = 0;          // count of k-mers read
    int count = 0;                  // count of samples present in the hash maps currently

    double diskReadTime = 0;
    while(true)   // read a sequence
    {
        high_resolution_clock::time_point readStart = high_resolution_clock::now();

        if(kseq_read(seq) < 0)
            break;

        high_resolution_clock::time_point readEnd = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double>>(readEnd - readStart);

        diskReadTime += time_span.count();

        ++total;    // one more sequence read
        //cout << "\r" << total << " processing ..." << flush;

        process_sequence(seq -> seq.s, th, no_kmers, count, k, maxSampleCount, MAP);
    }


    cout << "th: " << th << endl;                       // final value of the sampling parameter s
    cout << "No. of sequences: " << total << endl;      // total sequences read

    FILE *outpFilePtr = fopen(outpFile.c_str(), "w");                // file pointer for output file
    uint32_t csize = 0; //MAP.size();
    for(int i = th; i < 64; i++)
        csize += MAP[i].size();                         // total number of samples present in the hash maps;
                                                        // isn't it the same as 'count'?

    cout << "Number of distinct k-mers present in the hash maps: " << count << endl;
    cout << "Total size of the hash maps: " << csize << endl;

    unsigned long F0 = csize * pow(2, (th));            // Approximate number of distinct k-mers encountered;
                                                        // note that, csize is the number of distinct samples present
                                                        // in the hash maps, and we have ignored 'th' number of bits
                                                        // from each hash value (taken only the hashes with all
                                                        // trailing s bits being zero);
                                                        // considering a uniform distributions of bits in each of
                                                        // those 'th' bits, there are 2^th equally likely prefixes
                                                        // possible for each sample k-mer present.

                                                        // Another way of interpretation is that, the final sampling
                                                        // rate is 1/2^(th) i.e. we have kept one sample per 2^th samples;
                                                        // hence, scale csize by 2^th to get approximate distinct
                                                        // k-mer count.

    cout << "F0: " << F0 << endl;

    fprintf(outpFilePtr, "F1\t%lu\n", no_kmers);
    fprintf(outpFilePtr, "F0\t%lu\n", F0);

    cout << endl;
    cout << "total sequences: " << total << endl;       // total sequences read
    cout << "no_kmer: " << no_kmers << endl;            // total k-mers read

    //unsigned long freq[65];
   unsigned long *freq = new unsigned long[coverage];        // k-mer frequency distribution table;
                                                        // only interested in the k-mers with frequency <= coverage
   for(int i = 1; i <= coverage; i++)
        freq[i] = 0;

    for(int i = th; i < 64; i++)    // iterate over the hash maps (first 'th' maps have been dropped during sampling)
        for(auto& p : MAP[i])           // for each sample in hash map i
            if(p.second <= coverage)             // if its frequency does not exceed the coverage
                freq[p.second]++;               // add this k-mer's frequency to the distribution


    cout << "final th (s-value): " << th << endl;
    for(int i = 1; i <= coverage; i++)
    {
        unsigned long fff = (freq[i] * pow(2, th)); // approximation of f_i (scaled by 2^th, as the final sampling
                                                    // rate is 1 / 2^th)
        fprintf(outpFilePtr, "f%d\t%lu\n", i, fff);
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    //double elapsedSecs = double(endTime - beginTime) / CLOCKS_PER_SEC;
    double elapsedSecs = time_span.count();

    cout << "\n\nTime taken = " << elapsedSecs << " seconds\n" << endl;
    fprintf(outpFilePtr, "\n\nTime taken = %lf seconds\n", elapsedSecs);

    cout << "Disk read time " << diskReadTime << endl;

    fclose(outpFilePtr);

    // add kseq_t destroyer here

    return 0;
}
