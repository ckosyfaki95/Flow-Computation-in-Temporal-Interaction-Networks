/* Chrysanthi Kosyfaki, PhD Candidate, University of Ioannina, Greece */

#include "computeflow.h"

int main(int argc, char **argv)
{
	int i,j,k;
	FILE *f; // graph input file
	struct Graph G;
	
	struct DAG *G2;
	struct DAG *retDAG=NULL;
	
	struct Edge** edgearray;
	int numedges;

    clock_t t;
    double time_taken;
    
    double flow;
	
	edgearray = (struct Edge **)malloc(MAXEDGES*sizeof(struct Edge *)); //this array holds all edges that exist in valid paths 


	if (argc != 3)
	{
		printf("filename and source-id expected as arguments. Exiting...\n");
		return -1;
	}
	
	f = fopen(argv[1],"r");
	if (f==NULL)
	{
		printf("ERROR: file %s does not exist. Exiting...\n",argv[1]);
		return -1;
	}	
	
    if (read_graph(&G, f)==-1)
		return -1;
	
	int source = atoi(argv[2]);   
	int sink = source;
	
	int pathlen = 4;
	// write to edgearray distinct edges on paths from src to dest
	int totinter = findPaths2(G,source,sink,pathlen,edgearray,&numedges); 
	
	printf("Extracted paths: numedges=%d, totinter=%d\n",numedges, totinter);
	
	if (numedges==0) {
		printf("No paths found\n");
		return -1;
	} 
	
	
	//convert edgearray to DAG	
	G2 = edgearray2DAG(edgearray, numedges, sink, 0);
	printf("Created DAG: numnodes=%d\n",G2->numnodes);

	printf("processing DAG \n");
	printf("Greedy is running \n");
	t = clock();
	flow=computeFlowGreedy(*G2);
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Computed flow: %f\n", flow);
    printf("Total time of execution: %f seconds\n", time_taken);
    
    if (totinter<10000) {
		printf("LP is running \n");
		t = clock();    
		flow=computeFlowLP(*G2);
		t = clock() - t;
		time_taken = ((double)t)/CLOCKS_PER_SEC;
		printf("Computed flow: %f\n", flow);
		printf("Total time of execution: %f seconds\n", time_taken);
    }
	
    int *order = topoorder(G2);
	
	if (order) {
		printf("Preprocessing is running \n");
		int numdeletedinter=0;
		int numdeletededges=0;
		int numdeletednodes=0;
		t = clock(); 
		retDAG = preprocessDAG(G2, order, &numdeletedinter, &numdeletededges, &numdeletednodes);
	
		t = clock() - t;
		time_taken = ((double)t)/CLOCKS_PER_SEC;
		printf("Preprocessing time: %f seconds\n", time_taken);
		printf("Total deletions: interactions=%d, edges=%d, nodes=%d\n",numdeletedinter,numdeletededges,numdeletednodes);

		if(retDAG==NULL)
		{
			printf("DAG has no flow: sink disconnected\n");
		}
		else if (totinter-numdeletedinter<10000) {
			printf("LP after preprocessing is running \n");
			t = clock();    
			flow=computeFlowLP(*retDAG);
			printf("Computed flow: %f\n", flow);
			t = clock() - t;
			time_taken = ((double)t)/CLOCKS_PER_SEC;
			printf("Total time of execution: %f seconds\n", time_taken);
			
			printf("Graph simplification after preprocessing is running \n");
			t = clock(); 
			flow=compFlow(*retDAG, 0);
			t = clock() - t;
			time_taken = ((double)t)/CLOCKS_PER_SEC;
			printf("Computed flow: %f\n", flow);
			printf("Total time of execution: %f seconds\n", time_taken);
		}
	
		free(order);
	}

    if (G2 == retDAG) {
    	if (retDAG != NULL) 
    		freeDAG(G2); //no need to free retDAG
    }
    else {
        if (G2 != NULL) freeDAG(G2);
	    if (retDAG != NULL) freeDAG(retDAG);
    }
	free(edgearray);
	
    return 0;
}
