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

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"

using spp::sparse_hash_map;

using namespace std;
//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)        // The C header file kseq.h is a small library for parsing the FASTA/FASTQ format.
                            // For ordinary file I/O, you can use KSEQ_INIT(gzFile, gzread) to set the type of
                            // file handler and the read() function.
                            // FMI: http://lh3lh3.users.sourceforge.net/parsefastq.shtml
std::map<char, char> mapp = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}};


static const int MultiplyDeBruijnBitPosition[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

unsigned trailing_zeros(unsigned n) {
    return n ? __builtin_ctz(n) : -1;
}

static const char basemap[255] =
    {
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };


unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}

void printHelp()
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

int main(int argc, char** argv)
{

    if(argc == 1)
    {
      cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
      exit(0);
    }
    int n=31;                       // k-mer length
    int s = 25000000;               // sample size
    int cov = 64;                   // coverage (?)
    string f = "", outf = "";       // input and output FASTA files
    for (int c = 1; c < argc; c++)  // parse command-line input
        {

            if (!strcmp(argv[c], "-h"))       { printHelp(); }
            else if (!strcmp(argv[c], "-k"))     { n = atoi(argv[c+1]); c++; }  // k-mer length
            else if (!strcmp(argv[c], "-f"))    { f = argv[c+1]; c++; }         // input file
            else if (!strcmp(argv[c], "-s"))    { s = atoi(argv[c+1]); c++; }   // sample size
            else if (!strcmp(argv[c], "-c"))    { cov = atoi(argv[c+1]); c++; } // coverage (?)
            else if (!strcmp(argv[c], "-o")) { outf = argv[c+1]; c++; }         // output file
        }

       if (f.empty()  || outf.empty())  // empty file(s) mentioned
        {
          printHelp();
        }

    FILE *fp;
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

    int l;

    fp = fopen(f.c_str(), "r");     // file pointer for input FASTA file
    if( fp == Z_NULL){
      cout<<"File: "<< f  << " does not exist" <<endl;
      exit(1);
    }


    /*
        The C header file kseq.h is a small library for parsing the FASTA/FASTQ format.
        Function kseq_init() is used to initialize the parser
        and kseq_destroy() to destroy it.
        Function kseq_read() reads one sequence and fills the kseq_t struct
        FMI: http://lh3lh3.users.sourceforge.net/parsefastq.shtml
    */
    seq = kseq_init(fileno(fp));    // seq is the FASTA input parser
    int k = s;                      // sample size

    typedef sparse_hash_map<uint64_t, uint32_t> SMap;
    vector<SMap> MAP(64);           // array of hashmaps for sampled k-mers

    cout << "read the Sequences .. " << endl;

    int th = 0;                     // sample-size adaptation parameter 's';
                                    // the threshold count of the trailing zeroes for hash values
    uint64_t total = 0;             // count of sequences read
    uint64_t no_kmers = 0;          // count of k-mers read
    int count = 0;                  //
    uint64_t hash=0, fhVal=0, rhVal=0;
    while ((l = kseq_read(seq)) >= 0)   // read a sequence
    {
        ++total;    // one sequence read
        //cout << "\r" << total << " processing ..." << flush;

        int len = strlen(seq->seq.s);           // length of the sequence
        ntHashIterator itr(seq->seq.s, 1, n);   // ntHash iterator to iterate over the read sequence and provide
                                                // ntHash values for each k-mer of length n; initialized with the first
                                                // length-n window on the sequence

        while (itr != itr.end())                            // iterate until the last window
        {
            hash = (*itr)[0];                               // get the ntHash value
            ++no_kmers;                                     // one more k-mer encountered
            uint8_t tz = trailing_zeros(hash);              // #trailing_zeroes of this k-mer
            if(tz >= th)                                    // if #trailing_zeroes is greater than or equal to the
                                                            // adaptive sampling parameter 's' (threshold value)
            {
                if(MAP[tz].find(hash) != MAP[tz].end())     // k-mer already present in hash map
                    MAP[tz][hash] += 1;                     // increment k-mer count
                else                                        // k-mer absent in hash map
                {
                  MAP[tz].insert(make_pair(hash, 1));       // insert k-mer into hash map
                  ++count;                                  // increment #samples_present by one
                  //cout << "\r" << "count: " << count << flush;// << endl;

                    /*
                        I guess there is a potential bug present in this conditional check. The conditional should be
                        (count >= k) instead of (count == k).

                        Assume that you have reached the max sample count 'k', and now have dropped the s'th hash map
                        (MAP[th] here). Hence, the number of samples present 'count' drop down to
                        (count - MAP[th].size()). Suppose that the hash map MAP[th] did not have any entries present
                        prior to deletion. Hence, 'count' retains the same value 'k' as earlier. Now at the next
                        iteration, suppose that you encounter a new k-mer, insert it to one of the hash maps accordingly,
                        and increment 'count' by 1. So it goes to k + 1.

                        BOOM! You will never drop any hash map anymore in the lifetime of this execution, and space usage
                        will not be optimized as theorized.

                        However, my proposal of the conditional (count >= k) also should present very unlikely yet
                        potential bugs too, in case of a stream of hash values with very long suffixes of trailing
                        zeroes, and the lower order entries in the hash map arrays being all empty.

                        Best solution is to use a loop of the following format:
                            while(count == k):
                                drop hash maps and update count
                    */
                  if(count == k){                           // max sample count reached;
                                                            // one hash map will be dropped now
                    int cnt = MAP[th].size();               // size of the hash map to be dropped
                    count = count - cnt;                    // #samples_present after dropping hash map
                    SMap().swap(MAP[th]);                   // drop hash map from memory
                    //MAP[th].clear(); //MAP[th].resize(0);
                    ++th;                                   // increment threshold (s parameter)
                    cout  << "count: " << count << endl;
                  }
                }
            }
	    ++itr;      // go over to the next window
        }
    }


    cout << "th: " << th << endl;                       // final value of sampling parameter s
    cout << "No. of sequences: " << total << endl;      // total sequences read
    FILE *fo = fopen(outf.c_str(), "w");                // file pointer for output file
    uint32_t csize = 0; //MAP.size();
    for(int i=th; i<64; i++) csize += MAP[i].size();    // total number of samples present in the hash maps;
                                                        // isn't it the same as 'count'?

    unsigned long F0 = csize * pow(2, (th));            // approximate number of distinct k-mers encountered;
                                                        // note that, csize is the number of distinct samples present
                                                        // in the hash maps, and we have ignored 'th' number of bits
                                                        // from each hash value; considering a uniform distributions of
                                                        // bits in each of those 'th' bits, there are 2^th equally likely
                                                        // prefixes possible for each sample k-mer present.
                                                        // another way to interpret is that, the final sampling rate is
                                                        // 1/2^th i.e. we have kept one sample per 2^th samples; hence,
                                                        // scale csize by 2^th to get approximate distinct k-mer count

    cout << "F0: " << F0 << endl;
    fprintf(fo, "F1\t%lu\n", no_kmers);
    fprintf(fo, "F0\t%lu\n", F0);
    cout << endl;
    cout << "total: " << total << endl;
    cout << "no_kmer: " << no_kmers << endl;            // total k-mers encountered
    //unsigned long freq[65];
   unsigned long *freq = new unsigned long[cov];        // k-mer frequency distribution table;
                                                        // only interested in the k-mers with frequency <= coverage
   for(int i=1; i<=cov; i++) freq[i] = 0;
    //unsigned long tot = 0;
    //int xx = 0;
    for(int i=th; i<64; i++)    // iterate over the hash maps (first 'th' maps have been dropped during sampling)
    {
      for(auto& p: MAP[i])      // for each sample in hash map i
      {
        if(p.second <= cov)         // if its frequency does not exceed the coverage input
            freq[p.second]++;       // add this k-mer's frequency to the distribution
      }
    }


    cout << "th: " << th << endl;
    for(int i=1; i<=cov; i++){
      unsigned long fff = (freq[i]*pow(2, th));     // approximation of f_i (scaled by 2^th, as the final sampling
                                                    // rate is 1 / 2^th)
      fprintf(fo, "f%d\t%lu\n", i, fff);
    }
    fclose(fo);

    return 0;
}
