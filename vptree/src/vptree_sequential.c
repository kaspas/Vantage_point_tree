#include <math.h>
#include <stdlib.h> 
#include <time.h> 
#include "vptree.h"
#include <stdio.h>
#include <string.h>

#include <sys/time.h>


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
/*
    always take right element as pivot
    so that doesnt need swap 
*/
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
  /*
  Quickselect with recursive calls
  k must be valid
  else it returns -1
  */
double kthSmallest(int *list,double *arr, int l, int r, int k) 
{ 

    if (k > 0 && k <= r - l + 1) { 
  
        int index = partition(list,arr, l, r); 

        if (index - l == k - 1) 
            return arr[index]; //it is equal return median value 
  
        if (index - l > k - 1)  
            return kthSmallest(list,arr, l, index - 1, k); //if index greater than k set the boundaries to the left side of the array
  
        return kthSmallest(list,arr, index + 1, r,(k - index + l - 1)); //set the boundaries to the right side and set the k correct
    } 
  
    return -1; 
} 
/*
    calculating L2 norm with VP the last element in X array,
    since X doesnt change the last element index is stored at list 
*/
void dis(int *list,double *distance,double *X,int n,int d){
    for(int i=0;i<n-1;i++){
        double temp=0;
        double temp2;
        for(int j=0;j<d;j++){
            temp2=X[list[n-1]*d+j]-X[list[i]*d+j];
            temp+=temp2*temp2;
        }
        distance[i]=sqrt(temp);
    }
}
//! Build vantage-point tree given input dataset X
/*!

  \param X      Input data points, stored as [n-by-d] array
  \param n      Number of data points (rows of X)
  \param d      Number of dimensions (columns of X)
  \return The vantage-point tree
*/
vptree * buildvp(double *X,int n,int d)
{
    int *list = (int *)malloc(n*sizeof(int)); //initializing index array
    for (int i=0;i<n;i++)list[i]=i; //setting list values to follow from 0 to n-1 
    vptree *tree = vpbuild(X,list,n,d);//build tree function 
    free(list);//free list array , i don't need him any more
    return tree;    
}
//! Build vantage-point tree given input dataset X
/*!

  \param X      Input data points, stored as [n-by-d] array
  \param list   Vector with indexes for the Input data points(initially having values 0 to n-1)
  \param n      Number of data points (rows of X)
  \param d      Number of dimensions (columns of X)
  \return The vantage-point tree
*/

vptree *vpbuild(double *X,int *list,int n,int d){
    vptree *tree=(vptree*)malloc(sizeof(vptree)); // allocating memory for vptree structure
    if(n==1){ //one point. so keep the index and the pointer VP and median distance is zero since it is leaf 
        tree->index=list[n-1];
        tree->VP=&X[list[n-1]*d];
        tree->median=0;
        tree->inner=NULL;//there is no inner points so NULL
        tree->outer=NULL;//there is no outer points so NULL
        return tree;
    }
    else if(n<=0){//it shouldn't be there free tree and return
        free(tree);
        return NULL;
    }
    else{
        //int num = rand()%n; 
        int num=n-1;
        tree->VP=&X[list[num]*d];// take as VP last element 
        double *distance=(double *)malloc((n-1)*sizeof(double));
        dis(list,distance,X,n,d);    //calculate distances  
        int k=floor((n)/2);
        tree->index=list[num]; //save index
        tree->median=kthSmallest(list,distance,0,n-2,k); //split in half the Input data 
        //keep the indexes that changed in list  
        free(distance);//free distance
        /*
            recursive call to the same function till there are no other points left
            inner -> takes X, list and the size if from (previous left value) till kth element
            outer -> set the list to kth element for starting point till n-k-1 elemeth  
        */
        tree->inner=vpbuild(X,list,k,d);
        tree->outer=vpbuild(X,&list[k],n-k-1,d);
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

