#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <float.h>
#include<string.h>
#include<stdbool.h>

struct TourData {
    int* tour;
    double tourSize;
};

struct TourData nearestAddition(double **distanceMatrix, int numOfCoords, int startingNode){

    int visitedCount = 0;

    // Memory allocation for tour and visited Array
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    bool *visited = malloc(numOfCoords * sizeof(bool));

    //Initialising tour and visited array to empty
    int p=0;
    for(p=0; p<numOfCoords; p++)
    {
     tour[p] = 0;
     visited[p] = false;
    }

    tour[numOfCoords] = 0;


    int nearest;

    int noOfThreads = omp_get_max_threads();

    double *minimumAdditionalCosts = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *nearestNodes = (int *) malloc(noOfThreads * sizeof(int));

    // Starting tour with node received from main method
    tour[0] = startingNode;
    // Marked the node as visited and increased count of visited nodes
    visited[startingNode] = true;
    visitedCount++;

    double minimumDistance = DBL_MAX;

    // Finding the nearest vertex
    int i;
    for (i = 0; i < numOfCoords; i++)
    {
        if (!visited[i])
        {
            if (distanceMatrix[startingNode][i] < minimumDistance)
            {
                minimumDistance = distanceMatrix[startingNode][i];
                nearest = i;
            }
        }
    }

    // Adding the nearest vertex in the tour
    tour[1] = nearest;
    visited[nearest] = true;
    visitedCount++;

    // Adding the starting node again to complete the partial tour
    tour[2] = startingNode;

    while (visitedCount < numOfCoords)
    {
        int y=0;
        for (y = 0; y <noOfThreads; y++)
        {
            minimumAdditionalCosts[y] = DBL_MAX;
            positions[y] = 0;
            nearestNodes[y] = 0;
        }

        double minimumAdditionalCost = DBL_MAX;
        int minN;
        int minUnvisited;
        int positionToAdd, position;

        int j=0;
        int threadID;

        #pragma omp parallel for collapse(2) private(i, j, threadID) shared(visited, distanceMatrix, minimumAdditionalCosts, positions, nearestNodes)
        for ( i = 0; i < visitedCount; i++)
        {
            for ( j = 0; j < numOfCoords; j++)
            {
                threadID = omp_get_thread_num();
                if (!visited[j])
                {
                    double additionalCost = distanceMatrix[tour[i]][j];
                    if (additionalCost < minimumAdditionalCosts[threadID])
                    {
                        minimumAdditionalCosts[threadID] = additionalCost;
                        positions[threadID] = i;
                        nearestNodes[threadID] = j;
                    }
                }
            }
        }


        int x = 0;
        for (x = 0; x < noOfThreads; x++)
        {

            if (minimumAdditionalCosts[x] < minimumAdditionalCost)
            {
                minimumAdditionalCost = minimumAdditionalCosts[x];
                minN = positions[x];
                minUnvisited = nearestNodes[x];
            }
        }

        int indexPrevious = minN == 0 ? visitedCount - 1 : minN - 1;
        int indexNext = minN + 1;

        double distanceBefore = distanceMatrix[tour[minN]][minUnvisited] + distanceMatrix[tour[indexPrevious]][minUnvisited] -
                distanceMatrix[tour[minN]][tour[indexPrevious]];

        double distanceAfter = distanceMatrix[tour[minN]][minUnvisited] + distanceMatrix[tour[indexNext]][minUnvisited] -
                distanceMatrix[tour[minN]][tour[indexNext]];

        if (distanceAfter < distanceBefore)
        {
            minN = indexNext;
        } else
        {
            minN = indexPrevious + 1;
        }

        visitedCount++;

        for (i = visitedCount; i > minN; i--)
        {
            tour[i] = tour[i - 1];
        }

        tour[minN] = minUnvisited;
        visited[minUnvisited] = true;
    }

    struct TourData tourData;
    tourData.tour = malloc((numOfCoords + 1) * sizeof(int));

    int count =0;
    for(count =0; count< numOfCoords+1; count++)
    {
        tourData.tour[count] = tour[count];
    }

    double totalLength = 0;
//    int i;
    for (i = 0; i <=numOfCoords; i++)
    {
        if(i>0)
        {
            totalLength += distanceMatrix[tour[i]][tour[i - 1]];
        }
    }

    tourData.tourSize = totalLength;

    return tourData;
}

