#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include<omp.h>
#include <float.h>


struct TourData {
    int* tour;
    double tourSize;
};


struct TourData nearestAddition(double **distances, int numOfCoords, int startingNode) {

    int visitedCount = 0;
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    bool *visited = malloc(numOfCoords * sizeof(bool));

    int p =0;
    for(p =0; p< numOfCoords; p++)
    {
     tour[p] = 0;
     visited[p] = false;
    }

    tour[numOfCoords] = 0;

    int nearest;
    struct TourData result;
    result.tour = malloc((numOfCoords + 1) * sizeof(int));

    int noOfThreads = omp_get_max_threads();

    double *minimumAdditionalCosts = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *nearestNodes = (int *) malloc(noOfThreads * sizeof(int));

    tour[0] = startingNode;
    visited[startingNode] = true;
    visitedCount++;

    double minimumDistance = DBL_MAX;

    for (int i = 0; i < numOfCoords; i++) {
        if (!visited[i]) {
            if (distances[startingNode][i] < minimumDistance) {
                minimumDistance = distances[startingNode][i];
                nearest = i;
            }
        }
    }

    tour[1] = nearest;
    visited[nearest] = true;
    visitedCount++;

    tour[2] = startingNode;

    while (visitedCount < numOfCoords)
    {
        double minAdditionCost = DBL_MAX;
        double minimum_Cost = DBL_MAX;
        int min_position;
        int min_Unvisited_node;
        int positionToAdd, position;
        int y=0;
        double max = DBL_MAX;
        int threadID;

        for (y = 0; y < noOfThreads; y++) {
            minimumAdditionalCosts[y] = DBL_MAX;
            positions[y] = 0;
            nearestNodes[y] = 0;
        }
        int i=0, j=0;
        #pragma omp parallel for collapse(2) private(i, j, threadID) shared(visited, distances, minimumAdditionalCosts, positions, nearestNodes)
        for ( i = 0; i < visitedCount; i++) {
            for ( j = 0; j < numOfCoords; j++) {
                threadID = omp_get_thread_num();
                if (!visited[j]) {

                    double additionalCost = distances[tour[i]][j];
                    if (additionalCost < minimumAdditionalCosts[threadID]) {
                        minimumAdditionalCosts[threadID] = additionalCost;
                        positions[threadID] = i;
                        nearestNodes[threadID] = j;
                    }
                }
            }
        }


        int x = 0;
        for (x = 0; x < noOfThreads; x++) {

            if (minimumAdditionalCosts[x] < minimum_Cost) {
                minimum_Cost = minimumAdditionalCosts[x];
                min_position = positions[x];
                min_Unvisited_node = nearestNodes[x];
            }
        }

        int indexBefore = min_position == 0 ? visitedCount - 1 : min_position - 1;
        int indexAfter = min_position + 1;

        double distanceAfter =
                distances[tour[min_position]][min_Unvisited_node] + distances[tour[indexAfter]][min_Unvisited_node] -
                distances[tour[min_position]][tour[indexAfter]];
        double distanceBefore =
                distances[tour[min_position]][min_Unvisited_node] + distances[tour[indexBefore]][min_Unvisited_node] -
                distances[tour[min_position]][tour[indexBefore]];

        if (distanceAfter < distanceBefore) {
            min_position = indexAfter;
        } else {
            min_position = indexBefore + 1;
        }

        visitedCount++;
        visited[min_Unvisited_node] = true;

        for (i = visitedCount; i > min_position; i--) {
            tour[i] = tour[i - 1];
        }
        tour[min_position] = min_Unvisited_node;


    }

    for (int i = 0; i <= numOfCoords; i++) {
        printf("%d ", tour[i]);
        result.tour[i] = tour[i];
    }

    double totalLength = 0;
    for ( i = 0; i <=numOfCoords; i++) {
        printf("%d ", tour[i]);
        if(i>0) {
            totalLength += distances[tour[i]][tour[i - 1]];
        }
    }

    result.tourSize = totalLength;

//    double cost = 0.0;
//    for (int i = 0; i < numOfCoords - 1; i++) {
//        int fromVertex = tour[i];
//        int toVertex = tour[i + 1];
//
//        cost += distances[fromVertex][toVertex];
//
//    }
//
//    cost += distances[tour[numOfCoords - 1]][tour[0]];

//    result.tourSize = totalLength;

    return result;
}

