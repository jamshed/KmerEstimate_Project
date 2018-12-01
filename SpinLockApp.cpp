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

#include<ctime>
#include<thread>
#include<mutex>
#include<condition_variable>
#include<chrono>
#include<atomic>

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"
#define MAX_THREADS 1003
// #define THREAD_COUNT 4
// #define MAX_LEN 1024

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

void parse_input(int argc, char **argv, int &k, int &maxSampleCount, int &threadCount, int &coverage,
                 bool &memUnconstrained, string &inpFile, string &outpFile)
{
    for (int c = 1; c < argc; c++)              // parse command-line input
    {
        if(!strcmp(argv[c], "-h"))              // help prompt
            print_help();
        else if(!strcmp(argv[c], "-k"))
        {
            k = atoi(argv[c+1]);                // k-mer length
            c++;
        }
        else if(!strcmp(argv[c], "-f"))
        {
            inpFile = argv[c+1];                // input file
            c++;
        }
        else if(!strcmp(argv[c], "-s"))
        {
            maxSampleCount = atoi(argv[c+1]);   // sample size
            c++;
        }
        else if(!strcmp(argv[c], "-c"))
        {
            coverage = atoi(argv[c+1]);         // coverage
            c++;
        }
        else if(!strcmp(argv[c], "-o"))
        {
            outpFile = argv[c+1];               // output file
            c++;
        }
        else if(!strcmp(argv[c], "-m"))         // whether unconstraining the memory
            memUnconstrained = true;
        else if(!strcmp(argv[c], "-t"))
        {
            threadCount = atoi(argv[c+1]);      // thread count
            c++;
        }
    }
}


void process_sequence(char *s, int k, int maxSampleCount, int &currSampleCount, int &th, uint64_t &kmerCount,
                      vector<SMap> &MAP)
{
    uint8_t tz;
    int dropSize;
    uint64_t hashVal;

    ntHashIterator itr(s, 1, k);    // ntHash iterator to iterate over the read sequence and provide
                                    // ntHash values for each of its k-mers of length k; initialized with
                                    // the first length-k window on the sequence

                                    // why does parameter 1 equal to 1 instead of 0? (legacy code)

    while(itr != itr.end())                             // iterate until the last k-mer window
    {
        kmerCount++;

        hashVal = (*itr)[0];                            // get the ntHash value
        tz = trailing_zeros(hashVal);                   // #trailing_zeroes of this k-mer

        if(tz >= th)                                    // if #trailing_zeroes is greater than or equal to threshold.
        {                                               // then sample this k-mer
            auto p = MAP[tz].find(hashVal);

            if(p != MAP[tz].end())                      // k-mer already present in hash map
                p -> second++;                          // increment k-mer count
            else                                        // k-mer absent in hash map
            {
                MAP[tz].insert(make_pair(hashVal, 1));  // insert k-mer into hash map
                currSampleCount++;                      // increment #samples_present by one

                while(currSampleCount >= maxSampleCount)    // max sample count reached
                {                                           // some hash maps will be dropped now
                    //cout << "Samples count reached " << currSampleCount << endl;

                    dropSize = MAP[th].size();              // size of the hash map to be dropped
                    SMap().swap(MAP[th]);                   // drop hash map from memory

                    //cout << "Dropping a hash map of size " << dropSize << endl;

                    currSampleCount -= dropSize;            // #samples_present after dropping corresponding hash map
                    //MAP[th].clear(); //MAP[th].resize(0);

                    ++th;                                   // increment threshold (s parameter)

                    //cout  << "New samples count: " << currSampleCount << endl;
                }
            }
        }

        ++itr;      // go over to the next window
    }
}


struct DistributedCount
{
    int currSampleCount = 0;
    int th = 0;
    uint64_t kmerCount = 0;
    vector<SMap> MAP = vector<SMap>(64);

    /*
        sparse_hash_map is distinguished from other hash-map implementations by its stingy use of memory and by the
        ability to save and restore contents to disk. On the other hand, this hash-map implementation, while still
        efficient, is slower than other hash-map implementations.
        FMI: http://goog-sparsehash.sourceforge.net/doc/sparse_hash_map.html
    */
};

volatile bool readFinished;
volatile int seqAvailable[MAX_THREADS];
volatile bool threadFree[MAX_THREADS];
char *S[MAX_THREADS];

#define ABSENT 0
#define PRESENT 1
#define NO_MORE -1

//mutex bufferLock[THREAD_COUNT];
//condition_variable cv[THREAD_COUNT];

void thread_operation(DistributedCount &memory, int threadID, int k, int maxSampleCount)
{
    cout << "Thread " << threadID << " initiated." << endl;

    while(!readFinished || seqAvailable[threadID] == PRESENT)
    {
        //unique_lock<mutex> lck(bufferLock[threadID]);

        while(seqAvailable[threadID] == ABSENT);
        /*{
            cv[threadID].wait(lck);

            if(readFinished)
                break;
        }*/

        if(seqAvailable[threadID] == PRESENT)
        {
            process_sequence(S[threadID], k, maxSampleCount, memory.currSampleCount, memory.th, memory.kmerCount,
                             memory.MAP);

            //seqAvailable[threadID] = false;
            seqAvailable[threadID] = ABSENT;

            threadFree[threadID] = true;
            //cv[threadID].notify_one();
        }
    }
}

void consolidate_outputs(DistributedCount *distCount, DistributedCount &output, int maxSampleCount, int threadCount)
{
    int &currSampleCount = output.currSampleCount;
    int &th = output.th;
    uint64_t &kmerCount = output.kmerCount;
    vector<SMap> &MAP = output.MAP;

    for(int i = 0; i < threadCount; ++i)
    {
        kmerCount += distCount[i].kmerCount;

        cout << "th for thread " << i << " is " << distCount[i].th << endl;

        for(int j = distCount[i].th; j < 64; ++j)
        {
            for(auto &entry : distCount[i].MAP[j])
            {
                auto p = MAP[j].find(entry.first);

                if(p != MAP[j].end())
                    p -> second += entry.second;
                else
                {
                    MAP[j].insert(entry);
                    currSampleCount++;
                }
            }

            SMap().swap(distCount[i].MAP[j]);
        }
    }

    th = distCount[0].th;
}

void round_robin(string &inpFile, DistributedCount &output, int k, int maxSampleCount, int threadCount,
                 bool memUnconstrained)
{
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

    FILE *inpFilePtr = fopen(inpFile.c_str(), "r");     // file pointer for input FASTA file

    if(inpFilePtr == Z_NULL)
    {
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

    cout << "read the Sequences .. " << endl;

    readFinished = false;

    uint64_t total = 0;             // count of sequences read
    thread T[MAX_THREADS];
    kseq_t *seqReads[MAX_THREADS];

    DistributedCount distCount[MAX_THREADS];

    for(int i = 0; i < threadCount; ++i)
    {
        threadFree[i] = true;
        seqAvailable[i] = false;

        seqReads[i] = kseq_init(fileno(inpFilePtr));

        T[i] = thread(&thread_operation, ref(distCount[i]), i, k,
                      memUnconstrained ? maxSampleCount : maxSampleCount / threadCount);
    }

    int threadIdx = 0;  // thread index
    while(true)         // read a sequence
    {
        //unique_lock<mutex> lck(bufferLock[threadIdx]);

        while(!threadFree[threadIdx]);
            //cv[threadIdx].wait(lck);

        kstream_t *originalStream = seq -> f;
        seq = seqReads[threadIdx];
        seq -> f = originalStream;

        if(kseq_read(seq) < 0)
            break;

        seqReads[threadIdx] = seq;

        ++total;        // one more sequence read

        //cout << "sequence " << total << " read" << endl;
        //cout << "assigned to thread " << threadIdx << endl;

        threadFree[threadIdx] = false;

        //strcpy(&S[t][0], seq -> seq.s);
        S[threadIdx] = seqReads[threadIdx] -> seq.s;

        //seqAvailable[threadIdx] = true;
        seqAvailable[threadIdx] = PRESENT;
        //cv[threadIdx].notify_one();

        threadIdx = (threadIdx + 1) % threadCount;        //get to next thread
    }

    readFinished = true;

    cout << "read finished: " << readFinished << endl;

    for(int i = 0; i < threadCount; ++i)
    {
        cout << "joining thread " << i << endl;
        //cv[i].notify_one();
        seqAvailable[i] = NO_MORE;
        T[i].join();
    }

    consolidate_outputs(distCount, output, maxSampleCount, threadCount);

    cout << "No. of sequences: " << total << endl;      // total sequences read
}

void get_output(FILE *outpFilePtr, DistributedCount &output, int coverage)
{
    int &currSampleCount = output.currSampleCount;
    int &th = output.th;
    uint64_t &kmerCount = output.kmerCount;
    vector<SMap> &MAP = output.MAP;

    cout << "th: " << th << endl;                       // final value of the sampling parameter s

    uint32_t csize = 0;
    for(int i = th; i < 64; i++)
        csize += MAP[i].size();                         // total number of samples present in the hash maps;
                                                        // isn't it the same as 'count'?

    cout << "Number of distinct k-mers present in the hash maps: " << currSampleCount << endl;
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

    fprintf(outpFilePtr, "F1\t%lu\n", uint64_t(kmerCount));
    fprintf(outpFilePtr, "F0\t%lu\n", F0);

    cout << endl;
    //cout << "total sequences: " << total << endl;       // total sequences read
    cout << "no_kmer: " << kmerCount << endl;            // total k-mers read

   unsigned long *freq = new unsigned long[coverage];        // k-mer frequency distribution table;
                                                        // only interested in the k-mers with frequency <= coverage
   for(int i = 1; i <= coverage; i++)
        freq[i] = 0;

    for(int i = th; i < 64; i++)        // iterate over the hash maps (first 'th' maps have been dropped during sampling)
        for(auto& p : MAP[i])           // for each sample in hash map i
            if(p.second <= coverage)    // if its frequency does not exceed the coverage
                freq[p.second]++;       // add this k-mer's frequency to the distribution


    cout << "final th (s-value): " << th << endl;
    for(int i = 1; i <= coverage; i++)
    {
        unsigned long fff = (freq[i] * pow(2, th)); // approximation of f_i (scaled by 2^th, as the final sampling
                                                    // rate is 1 / 2^th)
        fprintf(outpFilePtr, "f%d\t%lu\n", i, fff);
    }

    fprintf(outpFilePtr, "\n\nfFinal value of th = %d\n", th);
}

int main(int argc, char** argv)
{
    // clock_t beginTime = clock();
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    if(argc == 1)
    {
        cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
        exit(0);
    }

    int k = 31;                             // default k-mer length
    int maxSampleCount = 25000000;          // default sample size
    int threadCount = 1;                    // default threads count
    int coverage = 64;                      // default coverage (maximum k-mer frequency we are interested in)
    bool memUnconstrained = false;          // whether the maxSampleCount is being distributed over the threads
    string inpFile = "", outpFile = "";     // input and output FASTA file names

    parse_input(argc, argv, k, maxSampleCount, threadCount, coverage, memUnconstrained, inpFile, outpFile);

    if (inpFile.empty()  || outpFile.empty())   // empty file(s) mentioned
        print_help();

    DistributedCount output;

    round_robin(inpFile, output, k, maxSampleCount, threadCount, memUnconstrained);


    FILE *outpFilePtr = fopen(outpFile.c_str(), "w");   // file pointer for output file

    get_output(outpFilePtr, output, coverage);

    //clock_t endTime = clock();
    high_resolution_clock::time_point t2 = high_resolution_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    //double elapsedSecs = double(endTime - beginTime) / CLOCKS_PER_SEC;
    double elapsedSecs = time_span.count();


    cout << "\n\nTime taken = " << elapsedSecs << " seconds\n" << endl;
    fprintf(outpFilePtr, "\n\nTime taken = %lf seconds\n", elapsedSecs);

    fclose(outpFilePtr);

    // add kseq_t destroyer here

    return 0;
}
