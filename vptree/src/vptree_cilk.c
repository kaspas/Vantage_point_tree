#include <math.h>
#include <stdlib.h> 
#include <time.h> 
#include "vptree.h"
#include <stdio.h>
#include <string.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h>

#include <sys/time.h>

#define LIMIT 20000 

void swap_double(double *a, double *b)
{
    double t = *a;
    *a = *b;
    *b = t;
}
void swap_int(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}
int partition(int *list,double *arr, int l, int r) 
{ 
    double x = arr[r];
    int i = l; 
    for (int j = l; j <= r - 1; j++) { 
        if (arr[j] <= x) { 
            swap_double(&arr[i], &arr[j]); 
            swap_int(&list[i],&list[j]);
            i++; 
        } 
    } 
    swap_double(&arr[i], &arr[r]); 
    swap_int(&list[i],&list[r]);
    return i; 
} 
  
double kthSmallest(int *list,double *arr, int l, int r, int k) 
{ 

    if (k > 0 && k <= r - l + 1) { 
  
        int index = partition(list,arr, l, r); 

        if (index - l == k - 1) 
            return arr[index]; 
  
        if (index - l > k - 1)  
            return kthSmallest(list,arr, l, index - 1, k); 
  
        return kthSmallest(list,arr, index + 1, r,(k - index + l - 1)); 
    } 
  
    return -1; 
} 

void dis(int *list,double *distance,double *X,int n,int d){
    
    if(n*d>LIMIT){
        cilk_for(int i=0;i<n-1;i++){ //let cilk decide what is best with work stealing
            double temp=0;
            double temp2;
            for(int j=0;j<d;j++){
                temp2=X[list[n-1]*d+j]-X[list[i]*d+j];
                temp+=temp2*temp2;
            }
            distance[i]=sqrt(temp);
        }
    }
    else{
        for(int i=0;i<n-1;i++){ //let cilk decide what is best with work stealing
            double temp=0;
            double temp2;
            for(int j=0;j<d;j++){
                temp2=X[list[n-1]*d+j]-X[list[i]*d+j];
                temp+=temp2*temp2;
            }
            distance[i]=sqrt(temp);
        }
    }
}
// type definition of vptree
// ========== LIST OF ACCESSORS
//! Build vantage-point tree given input dataset X
/*!
    \param X    Input data points, stored as [n-by-d] array
    \param n    Number of data points (rows of X)
    \param d    Number of dimensions (columns of X)
    \return     The vantage-point tree
*/

vptree * buildvp(double *X,int n,int d)
{
    int *list = (int *)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)list[i]=i;
    vptree *tree = vpbuild(X,list,n,d);//build tree function 
    free(list);
    return tree;   
    
}
vptree *vpbuild(double *X,int *list,int n,int d){
    vptree *tree=(vptree*)malloc(sizeof(vptree));
    if(n==1){

        tree->index=list[n-1];
        tree->VP=&X[list[n-1]*d];
        tree->median=0;
        tree->inner=NULL;
        tree->outer=NULL;
        return tree;
    }
    else if(n<=0){
        free(tree);
        return NULL;
    }
    else{
        //int num = rand()%n; 
        int num=n-1;
        tree->VP=&X[list[num]*d];
        
        double *distance=(double *)malloc((n-1)*sizeof(double));
        dis(list,distance,X,n,d);
        
        int k=floor((n)/2);
        tree->index=list[num];
        tree->median=kthSmallest(list,distance,0,n-2,k);
        
        free(distance);
        if(n*d>2*LIMIT){
          tree->inner=cilk_spawn vpbuild(X,list,k,d); //  vpbuild can be executed in parallel
        tree->outer=vpbuild(X,&list[k],n-k-1,d); // continue to this recursion   
        }
        else{
            tree->inner=vpbuild(X,list,k,d); //  vpbuild can be executed in parallel
            tree->outer=vpbuild(X,&list[k],n-k-1,d); // continue to this recursion 
        }
        
        cilk_sync; //wait for synchronization
        return tree;   
    }
}
//! Return vantage-point subtree with points inside radius
/*!
    \param  node    A vantage-point tree
    \return         The vantage-point subtree
*/

vptree * getInner(vptree * T){
   return T->inner;
}

//! Return vantage-point subtree with points outside radius
/*!
    \param  node    A vantage-point tree
    \return         The vantage-point subtree
*/

vptree * getOuter(vptree * T){
    return T->outer;
}

//! Return median of distances to vantage point
/*!
    \param  node    A vantage-point tree
    \return         The median distance
*/

double getMD(vptree * T){
    return T->median;
}

//! Return the coordinates of the vantage point
/*!
    \param node     A vantage-point tree
    \return         The coordinates [d-dimensional vector]
*/

double* getVP(vptree * T){
    return T->VP;
}

//! Return the index of the vantage point
/*!
    \param node     A vantage-point tree
    \return         The index to the input vector of data points
 */

int getIDX(vptree * T){
    return T->index;
}

