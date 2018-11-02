#include <iostream>
#include <thread>

void tmp_thr(double a) { std::cout<<a<<std::endl; }

int main()
{
     int j;

     for(j=0;j<10000000;j++)
     {
         std::thread tEx(tmp_thr,j);
         tEx.join();
     }
}
