#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "vptree.h"
#include <math.h>
#include <time.h>
#include "validation.h"
#include <omp.h>

int main(int argc,char **argv){
    if(argc!=4){
        printf("Correct Usage: .\\%s n d Nthreads\n",argv[0]);
        printf("Where n : Number of Particles\nd : Dimensions\nNthreads = Number of Threads\n");
        return 1;
    }
    int n=atoi(argv[1]);
    int d=atoi(argv[2]);
    //nthreads =atoi(argv[3]);
    //omp_set_num_threads(nthreads);
    struct timeval startwtime, endwtime;

    double *data = (double *)malloc(n*d*sizeof(double));
    srand((unsigned int)time(NULL));
    for(int i=0;i<n;i++){            
            for(int j=0;j<d;j++){
                data[i*d+j]=((float)rand()/(float)(RAND_MAX)) * 1000.58725;
        }
    }
    vptree *tree;
    //printf("Size %ld\n",sizeof(vptree  ));
    gettimeofday (&startwtime, NULL); 
    tree=buildvp(data,n,d);
    gettimeofday (&endwtime, NULL); 
    double t = (double)((endwtime.tv_usec - startwtime.tv_usec)
				/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    printf("Elapsed: %f seconds\n",t);
    printf("\n\n");
    verify(tree,data,n,d);
    printf("\n");
    deleteTree(&tree,d);
    free(data);

    return 0;

}
