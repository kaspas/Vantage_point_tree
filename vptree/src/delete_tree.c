#include "vptree.h"
#include <stdio.h>
#include <stdlib.h>

void _deleteTree(vptree* node,int d) 
{ 
    if (node == NULL) return; 
    
    /* first delete both subtrees */
    _deleteTree(getInner(node),d); 
    _deleteTree(getOuter(node),d); 
        
     //then delete the node */
    /*printf("{");
            for(int j=0;j<d;j++){
                if(j<d-1)printf("%.2f  ",node->VP[j]);
                else printf("%.2f}\n",node->VP[j]);
        }
        printf("Index Number : %d\n",getIDX(node));
        printf("Median Distance : %.2f\n",node->median);
    */
   free(node); 
} 
  
/* Deletes a tree and sets the root as NULL */
void deleteTree(vptree** node_ref,int d) 
{ 
    
  _deleteTree(*node_ref,d); 
  
  *node_ref = NULL; 
} 