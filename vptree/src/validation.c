#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "validation.h"
#include "vptree.h"

void verify(vptree * T,double *X,int n,int d)
{
    printf("Validation Procedure...\n");
    double *pvp=getVP(T);
    int index=getIDX(T);
    for(int i=0;i<d;i++){
        if(X[index*d+i]!=pvp[i]){
            printf("Wrong Index on Root Node.\n");
        }
    }
    if(validate(T,getInner(T),X,d,0) == 0 && validate(T,getInner(T),X,d,1) == 0){
        printf("Validation Completed.\n");
        printf("No index Error Found.\n");
        printf("No element in wrong subtree found.\n");
        printf("No error Found. \n");
    }
    else{
        printf("Validation Completed\n");
        printf("Errors occured.\n");
        printf("Check your implementation\n");
    }
}

double distance(double *X,double *Y,int d){
    double D=0;
    for(int i=0;i<d;i++){
        D+=(X[i]-Y[i])*(X[i]-Y[i]);
    }
    return sqrt(D);
}

int validate(vptree *P, vptree *C,double *X,int d,int isInner){
    double *pvp,*cvp;
    pvp=getVP(P);
    cvp=getVP(C);
    int id = getIDX(C);
    double md = getMD(P);
    for (int i=0;i<d;i++){
        if(X[id*d+i]!=cvp[i]){
            printf("Wrong Index on Node.\n");
            return -1;
        }
    }
    double D = distance(pvp,cvp,d);
        if(md<(D) && isInner==0){
            printf("Node at wrong place it is Inner but it should be Outer.\n");
            return -1;
        }
    if(md>(D) && isInner==1){
        printf("Node at wrong place it is Outer but it should be Inner.\n");
        return -1;
    }
    if(getInner(C)==NULL && getOuter(C)!=NULL){
        return validate(C,getOuter(C),X,d,1);
    }
    else if(getInner(C)!=NULL && getOuter(C)==NULL){
        return validate(C,getInner(C),X,d,0);
    }
    else if(getInner(C)==NULL && getOuter(C)==NULL){
        return 0;
    }
    else{

        return (validate(C,getInner(C),X,d,0)+validate(C,getOuter(C),X,d,1));
    }
}