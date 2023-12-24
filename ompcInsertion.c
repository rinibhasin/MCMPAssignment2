#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>
#include<stdbool.h>
#include <omp.h>
#include <stdlib.h>


struct TourData {
    int* tour;
    double tourSize;
};


struct TourData cheapestInsertion(double **distanceMatrix, int numOfCoords, int startingNode)
{
    int visitedCount = 0;

    int *tour = (int*)malloc((numOfCoords+1)*sizeof(int));
    bool *visited = (bool*)malloc(numOfCoords*sizeof(bool));
    int m=0;
    for(m=0; m<numOfCoords; m++)
    {
        tour[m] = 0;
        visited[m] =0;
    }

    tour[numOfCoords] = 0;

    // Initialise with the first vertex
    tour[0] = startingNode;
    visited[startingNode] = true;
    visitedCount++;

    // Find the nearest vertex
    double minimumDistance = DBL_MAX;

    int nearestVertex;
    int i = 0;

    int noOfThreads = omp_get_max_threads();

    for(i = 0 ; i <numOfCoords; i++)
    {
        if(i != startingNode && distanceMatrix[startingNode][i]< minimumDistance)
        {
            minimumDistance = distanceMatrix[startingNode][i];
            nearestVertex = i;
        }
    }

    // Add the nearest vertex in the tour
    tour[1]= nearestVertex;
    visited[nearestVertex] = true;
    visitedCount++; // 2
    tour[2] = startingNode;

    double *minimumAdditionalCosts = (double*)malloc(noOfThreads*sizeof(double));
    int *positions = (int*)malloc(noOfThreads*sizeof(int));
    int *nearestVertexes = (int*)malloc(noOfThreads*sizeof(int));

    int y =0;
    while(visitedCount < numOfCoords)
    {
        for(y =0;y< noOfThreads; y++)
        {
            minimumAdditionalCosts[y] = DBL_MAX;
            positions[y] =0;
            nearestVertexes[y]=0;
        }

        double minimumAdditionalCost = DBL_MAX;

        int minN;
        int minUnvisited;
        double additionalCost;
        int j;int threadID;
        // tour = {0,1}

            #pragma omp parallel for collapse(2) private(i,j, additionalCost, threadID) shared(visited, distanceMatrix, minimumAdditionalCosts, positions,nearestVertexes)
            for (i = 0; i < visitedCount; i++) {
                // unvisited nodes
                for (j = 0; j < numOfCoords; j++) {
                    threadID  = omp_get_thread_num();
                    // check for unvisited nodes
                    if (!visited[j]) {
                        // j =2
                        additionalCost = distanceMatrix[j][tour[i]] + distanceMatrix[j][tour[i + 1]] -
                                         distanceMatrix[tour[i]][tour[i + 1]];
                        if (additionalCost < minimumAdditionalCosts[threadID]) {

                            minimumAdditionalCosts[threadID] = additionalCost;
                            positions[threadID] = i;
                            nearestVertexes[threadID] = j;
                        }
                    }
                }
            }

        int x=0;
        for(x =0; x< noOfThreads; x++)
        {
            if(minimumAdditionalCosts[x]< minimumAdditionalCost)
            {
                minimumAdditionalCost = minimumAdditionalCosts[x];
                minN = positions[x];
                minUnvisited = nearestVertexes[x];
            }
        }

        // Make space to add unvisited node to computed index
        for(i = visitedCount; i > minN; i--)
        {
            tour[i+1] = tour[i];
        }


        // add the node to tour
        visited[minUnvisited] = true;
        tour[minN+1] = minUnvisited;


        visitedCount++;

    }

    double totalLength = 0;
    for ( i = 0; i <=numOfCoords; i++) {
        printf("%d ", tour[i]);
        if(i>0) {
            totalLength += distanceMatrix[tour[i]][tour[i - 1]];
        }
    }
    printf("%f", totalLength);

    struct TourData tourData;

    tourData.tour = malloc((numOfCoords+1) * sizeof(int));

    int count =0;
    for(count =0; count< numOfCoords+1; count++)
    {
        tourData.tour[count] = tour[count];
    }

    tourData.tourSize =totalLength;

    return tourData;
}



