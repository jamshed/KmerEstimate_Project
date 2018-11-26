Compile and Run
------------------------------
Compile:


		g++ -o kmerEst kmerCountEstimate.cpp -std=c++11 -O3 -march=native -pthread

Run original single-threaded implementation:

		./kmerEst -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>
  
  
  
Run multi-threaded implementation (Conditional-variable based approach):

		./distCount -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -t <thread count> -c <coverage> -o <out.txt>
		
		
Run multi-threaded implementation (Spin-lock based approach):

		./spin -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -t <thread count> -c <coverage> -o <out.txt>
  


Datasets
------------------------------
https://drive.google.com/drive/folders/1Nbp2Q38OY2_Y-CKe7pcxU6b3pLVnYuEL?usp=sharing
