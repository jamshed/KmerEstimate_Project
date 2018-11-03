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

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"
#define THREAD_COUNT 1
#define MAX_LEN 1024

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


mutex countLock, no_kmersLock;
vector<mutex> mapLock(65);

int th;
uint64_t no_kmers;
int k;
int maxSampleCount;
vector<SMap> MAP(64);


void process_sequence(char *s, int &count)
{
    uint64_t hash = 0;
    uint64_t kc = 0;

    ntHashIterator itr(s, 1, k); // ntHash iterator to iterate over the read sequence and provide
                                                // ntHash values for each of its k-mers of length n; initialized with
                                                // the first length-n window on the sequence

    while(itr != itr.end())                            // iterate until the last k-mer window
    {
        hash = (*itr)[0];                               // get the ntHash value
        kc++;
        uint8_t tz = trailing_zeros(hash);              // #trailing_zeroes of this k-mer

        if(tz >= th)                                    // if #trailing_zeroes is greater than or equal to threshold
        {                                               // then sample this k-mer
            mapLock[tz].lock();

            if(MAP[tz].find(hash) != MAP[tz].end())     // k-mer already present in hash map
            {
                MAP[tz][hash] += 1;                     // increment k-mer count
                mapLock[tz].unlock();
            }
            else                                        // k-mer absent in hash map
            {
                MAP[tz].insert(make_pair(hash, 1));     // insert k-mer into hash map
                mapLock[tz].unlock();

                countLock.lock();

                ++count;                                // increment #samples_present by one

                while(count >= maxSampleCount)                              // max sample count reached
                {                                           // one hash map will be dropped now
                    cout << "Samples count reached " << count << endl;

                    mapLock[th].lock();
                    int cnt = MAP[th].size();               // size of the hash map to be dropped
                    SMap().swap(MAP[th]);                   // drop hash map from memory
                    mapLock[th].unlock();

                    cout << "Dropping a hash map of size " << cnt << endl;

                    count = count - cnt;                    // #samples_present after dropping corresponding hash map
                    //MAP[th].clear(); //MAP[th].resize(0);

                    ++th;                                   // increment threshold (s parameter)
                    cout  << "New samples count: " << count << endl;
                }

                countLock.unlock();
            }
        }

        ++itr;      // go over to the next window
    }

    //kmerCountLock.lock();
    no_kmersLock.lock();
    no_kmers += kc;
    no_kmersLock.unlock();
    //kmerCountLock.unlock();
}

bool readFinished;
bool seqAvailable[THREAD_COUNT];
bool threadFree[THREAD_COUNT];

mutex bufferLock[THREAD_COUNT];
condition_variable cv[THREAD_COUNT];

void garbage()
{
    for(int i = 0; i < 100000; ++i);
}

void thread_operation(char *threadBuffer, int threadID, int &count)
{
    //cout << "Thread " << threadID << " initiated." << endl;

    while(!readFinished)
    {
        /*
        unique_lock<mutex> lck(bufferLock[threadID]);

        while(!seqAvailable[threadID])
            cv[threadID].wait(lck);

        if(readFinished)
            break;
            */

        if(seqAvailable[threadID])
        {
            //bufferLock[threadID].lock();

            process_sequence(threadBuffer, ref(count));

            seqAvailable[threadID] = false;
            threadFree[threadID] = true;

            //bufferLock[threadID].unlock();
        }
        else
            sleep(0.00001);
    }
}

int get_free_thread(bool *threadFree)
{
    int t = 0;

    //for(; !threadFree[i]; i = (i + 1) % THREAD_COUNT);
    while(true)
    {
        if(threadFree[t])
            break;

        t = (t + 1) % THREAD_COUNT;
        if(!t)
            sleep(0.00001);
    }

    return t;
}

int main(int argc, char** argv)
{
    clock_t beginTime = clock();

    if(argc == 1)
    {
        cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
        exit(0);
    }

    k = 31;                     // default k-mer length
    maxSampleCount = 25000000;  // default sample size
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
    //vector<SMap> MAP(64);           // array of hash maps for sampled k-mers

    cout << "read the Sequences .. " << endl;

    //int th = 0;                     // sample-size adaptation parameter 's';
                                    // (the threshold count of the trailing zeroes for hash values)
    uint64_t total = 0;             // count of sequences read
    int count;
    //uint64_t no_kmers = 0;          // count of k-mers read
    //int count = 0;                  // count of samples present in the hash maps currently

    double diskReadTime = 0;

    //thread thrd;
    //bool readFinished = false;
    //bool threadFree[THREAD_COUNT];
    //bool seqAvailable[THREAD_COUNT];
    thread T[THREAD_COUNT];
    char S[THREAD_COUNT][MAX_LEN + 1];

    mutex no_kmersLock;


    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        threadFree[i] = true;
        seqAvailable[i] = false;

        T[i] = thread(&thread_operation, &S[i][0], i, ref(count));
    }

    while(true)   // read a sequence
    {
        double readStart = clock();

        if(kseq_read(seq) < 0)
            break;

        double readEnd = clock();

        diskReadTime += readEnd - readStart;

        ++total;    // one more sequence read

        //cout << "sequence " << total << " read" << endl;

        int t = get_free_thread(threadFree);

        //cout << "assigned to thread " << t << ". Its sequences availability = " << seqAvailable[t] << endl;

        //unique_lock<mutex> lck(bufferLock[t]);

        threadFree[t] = false;

        //bufferLock[t].lock();
        //lock_buffer(t);
        strcpy(&S[t][0], seq -> seq.s);

        //bufferLock[t].unlock();
        //unlock_buffer(t);

        seqAvailable[t] = true;

        //cv[t].notify_one();
    }

    readFinished = true;

    cout << "read finished: " << readFinished << endl;

    for(int i = 0; i < THREAD_COUNT; ++i)
    {
        cout << "joining thread " << i << endl;
        //cv[i].notify_one();
        T[i].join();
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

    clock_t endTime = clock();
    double elapsedSecs = double(endTime - beginTime) / CLOCKS_PER_SEC;

    cout << "\n\nTime taken = " << elapsedSecs << " seconds\n" << endl;
    fprintf(outpFilePtr, "\n\nTime taken = %lf seconds\n", elapsedSecs);

    cout << "Disk read time " << diskReadTime / CLOCKS_PER_SEC << endl;

    fclose(outpFilePtr);

    // add kseq_t destroyer here

    return 0;
}
