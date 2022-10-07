#include <iostream>
#include <algorithm>
#include <climits>
#include <cassert>
#include <omp.h>
#include <fstream>
#include <string>
#include <cstring>
using namespace std;

// int v[10000010];
// int t[10000010];
// int v2[10000010];

// int mergesort_count = 0;
// int merge_count = 0;

int main(int argc, char** argv) {

    // read inputs
    if(argc != 5) {
        cout << "Usage: ./a.out seed length basecase nthreads" << endl;
        return 1;
    }
    srand(atoi(argv[1]));
    int length = atoi(argv[2]);
    int bc = atoi(argv[3]);
    int nthreads = atoi(argv[4]);
    
    

    // fstream myfile;


    // freopen(argv[5], "w", stdout);
    // myfile << "Writing this to a file.\n";
    // printf("test!\n");

    // // allocate memory
    int* v = new int[length]; // array to be sorted
    int* t = new int[length]; // temporary workspace
    int* v2 = new int[length]; // copy for checking

    // initialize input randomly
    for (int i = 0; i < length; i++) {
        v[i] = v2[i] = rand() % length;
    }
    omp_set_num_threads(nthreads);
    double start_sort = omp_get_wtime();
    #pragma omp parallel   // [Enze]
    {
        sort(v2, v2+length);
    }
    double elapsed_sort = omp_get_wtime() - start_sort;
    cout << "nthreads:   \t" << nthreads << endl;
    cout << "time:   \t" << elapsed_sort << endl;

    delete [] v2;
    delete [] t;
    delete [] v;

    return 0;
}