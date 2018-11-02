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


// New headers added
#include<ctime>
#include<mutex>
#include<thread>

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"

// New macros defined
#define THREAD_COUNT 8
#define MAX_SEQ_LEN 1025

using spp::sparse_hash_map;
typedef sparse_hash_map<uint64_t, uint32_t> SMap;

using namespace std;
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

void parse_input(int argc, char **argv, int &k, int &maxSampleCount, int &coverage, string &inpFile, string &outpFile)
{
    for (int c = 1; c < argc; c++)  // parse command-line input
    {
        if(!strcmp(argv[c], "-h"))      // help prompt
            printHelp();
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


void process_sequence(char *s, int &th, uint64_t &no_kmers, int &count, int k, int maxSampleCount, vector<SMap> &MAP,
                      mutex &thLock, mutex &no_kmersLock, mutex &countLock, vector<mutex> &mapLock, mutex &mapDropLock,
                      bool *threadFree, int threadNum)
{
    uint64_t hash = 0;
    int kmerCount = 0;
    int newSample = 0;

    ntHashIterator itr(s, 1, k);    // ntHash iterator to iterate over the read sequence and provide
                                    // ntHash values for each of its k-mers of length n; initialized with
                                    // the first length-n window on the sequence

    while (itr != itr.end())                            // iterate until the last k-mer window
    {
        //if(count > maxSampleCount)
        //    cout << "\n\n\nCluster fuck\n\n\n";

        hash = (*itr)[0];                               // get the ntHash value
        kmerCount++;                                     // one more k-mer read


        // can add micro-optimizations of bit operations here
        uint8_t tz = trailing_zeros(hash);              // #trailing_zeroes of this k-mer


        if(tz >= th)                             // if #trailing_zeroes is greater than or equal to threshold
        {                                               // then sample this k-mer
            mapLock[tz].lock();
            auto p = MAP[tz].find(hash), e = MAP[tz].end();
            //mapLock[tz].unlock();

            // mapLock[tz].lock();

            if(p != e)     // k-mer already present in hash map
            {
                // mapLock[tz].lock();

                MAP[tz][hash] += 1;                     // increment k-mer count

                mapLock[tz].unlock();
            }
            else                                        // k-mer absent in hash map
            {
                //countLock.lock();

                //mapLock[tz].lock();

                MAP[tz].insert(make_pair(hash, 1));     // insert k-mer into hash map

                mapLock[tz].unlock();

                //cout << "Thread " << threadNum << "took hold of count lock" << endl;

                newSample++;


                //cout << "Thread " << threadNum << "releasing hold of count lock" << endl;

            }
        }

        ++itr;      // go over to the next window
    }

    no_kmersLock.lock();
    no_kmers += kmerCount;
    no_kmersLock.unlock();


    countLock.lock();

    count += newSample;

    if(count >= maxSampleCount)          // max sample count reached
    {                                           // one hash map will be dropped now
        cout << "Samples count reached " << count << endl;

        //mapDropLock.lock();
        mapLock[th].lock();

        int cnt = MAP[th].size();               // size of the hash map to be dropped
        SMap().swap(MAP[th]);                   // drop hash map from memory

        mapLock[th].unlock();

        cout << "Dropping a hash map of size " << cnt << endl;

        //countLock.lock();
        count -= cnt;                    // #samples_present after dropping corresponding hash map
        //countLock.unlock();

        //MAP[th].clear(); //MAP[th].resize(0);
        //thLock.lock();
        ++th;
        //thLock.unlock();

        // mapDropLock.unlock();
                                  // increment threshold (s parameter)
        cout  << "New samples count: " << count << endl;
    }

    countLock.unlock();

    threadFree[threadNum] = true;
    //pthread_join(*this);
}


/*
    Use a stack / queue based method instead, later.
*/
int get_free_thread(bool *threadFree)
{
    int i;

    for(i = 0;; i = (i + 1) % THREAD_COUNT)
        if(threadFree[i])
            break;

    return i;
}

int main(int argc, char** argv)
{
    clock_t beginTime = clock();

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
        printHelp();

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
        sparse_hash_map is distinguished from other hash-map implementations by its stingy use of memory and by the
        ability to save and restore contents to disk. On the other hand, this hash-map implementation, while still
        efficient, is slower than other hash-map implementations.
        FMI: http://goog-sparsehash.sourceforge.net/doc/sparse_hash_map.html
    */
    vector<SMap> MAP(64);           // array of hash maps for sampled k-mers

    mutex thLock, no_kmersLock, countLock, mapDropLock;
    vector<mutex> mapLock(64);

    cout << "read the Sequences .. " << endl;

    int th = 0;                     // sample-size adaptation parameter 's';
                                    // (the threshold count of the trailing zeroes for hash values)
    uint64_t total = 0;             // count of sequences read
    uint64_t no_kmers = 0;          // count of k-mers read
    int count = 0;                  // count of samples present in the hash maps currently

    /*
        The C header file kseq.h is a small library for parsing the FASTA/FASTQ format.
        Function kseq_init() is used to initialize the parser and kseq_destroy() to destroy it.
        Function kseq_read() reads one sequence and fills the kseq_t struct.
        FMI: http://lh3lh3.users.sourceforge.net/parsefastq.shtml
    */
    seq = kseq_init(fileno(inpFilePtr));    // seq is the FASTA input parser

    bool done = false;
    thread T[THREAD_COUNT];

    kseq_t *S[THREAD_COUNT];
    bool threadFree[THREAD_COUNT];
    // stack<int> freeThreads;
    bool threadUsedOnce[THREAD_COUNT];
    char sequences[THREAD_COUNT][MAX_SEQ_LEN + 1];

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        S[i] = kseq_init(fileno(inpFilePtr));
        threadFree[i] = true;
        threadUsedOnce[i] = false;
    }

    double ioTime = 0;

    while(true)
    {
        int t = get_free_thread(threadFree);

        threadFree[t] = false;

        kstream_t *originalStream = seq -> f;
        seq = S[t];
        seq -> f = originalStream;

        clock_t ioStart = clock();

        if(kseq_read(seq) < 0)
        {
            //cout << "\n\n\n\nBreaking Bad. All hells breaking loose at sequence " << total << "\n\n\n\n";
            // cout << "seq : " << seq -> seq.s << endl;
            break;
        }

        clock_t ioEnd = clock();

        ioTime += (double)(ioEnd - ioStart);

        S[t] = seq;

        total++;

        if(threadUsedOnce[t])
            //pthread_join(T[t]);
        {
            T[t].join();
        }
        threadUsedOnce[t] = true;

        T[t] = thread(&process_sequence, S[t] -> seq.s, ref(th), ref(no_kmers), ref(count), k, maxSampleCount,
                      ref(MAP), ref(thLock), ref(no_kmersLock), ref(countLock), ref(mapLock), ref(mapDropLock),
                      threadFree, t);

    }

    for(int i = 0; i < THREAD_COUNT; ++i)
        if(!threadFree[i])
            T[i].join();

    /*
    while(kseq_read(seq) >= 0)   // read a sequence
    {
        ++total;    // one more sequence read
        //cout << "\r" << total << " processing ..." << flush;
        process_sequence(seq -> seq.s, th, no_kmers, count, k, maxSampleCount, MAP,
                         thLock, no_kmersLock, countLock, mapLock);
    }
    */

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

    clock_t endTime = clock();
    double elapsedSecs = double(endTime - beginTime) / CLOCKS_PER_SEC;

    cout << "\n\nTime taken = " << elapsedSecs << " seconds\n" << endl;
    fprintf(outpFilePtr, "\n\nTime taken = %lf seconds\n", elapsedSecs);
    cout << "IO time = " << ioTime / CLOCKS_PER_SEC << endl;

    fclose(outpFilePtr);

    // add kseq_t destroyer here

    return 0;
}
