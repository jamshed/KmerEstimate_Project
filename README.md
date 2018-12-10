Compile and Run
------------------------------

Compile original single-threaded algorithm implementation:

		g++ -o kmerEst kmerCountEstimate.cpp -std=c++11 -O3 -march=native
		
Compile parallel algorithm implementation(s):

		g++ -o monitor Monitor.cpp -std=c++11 -O3 -march=native -pthread
		
		g++ -o spinlock Spinlock.cpp -std=c++11 -O3 -march=native -pthread



Run original (serial) algorithm implementation:

		./kmerEst -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>
  
  
  
Run parallel (multi-threaded) implementation (Monitor i.e. lock and conditional-variable based approach):

		./monitor -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -t <thread count> -c <coverage> -o <out.txt>
		
		
Run parallel (multi-threaded) implementation (Spinlock based approach):

		./spinlock -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -t <thread count> -c <coverage> -o <out.txt>
  


Datasets
------------------------------
https://drive.google.com/drive/folders/1Nbp2Q38OY2_Y-CKe7pcxU6b3pLVnYuEL?usp=sharing
