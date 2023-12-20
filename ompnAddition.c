#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include<omp.h>
#include <float.h>


int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

struct TourData {
    int* tour;
    double tourSize;
};


double **createDistanceMatrix(double **coords, int numOfCoords){
    int i, j;

    double **dMatrix = (double **)malloc(numOfCoords * sizeof(double));

    for(i = 0; i < numOfCoords; i++){
        dMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

#pragma omp parallel for collapse(2)
    for(i = 0; i < numOfCoords; i++){
        for(j = 0; j < numOfCoords; j++){
            double diffX = coords[i][0] - coords[j][0];
            double diffY = coords[i][1] - coords[j][1];
            dMatrix[i][j] = sqrt((diffX * diffX) + (diffY * diffY));
        }
    }

    return dMatrix;
}

struct TourData nearestAddition(double **distances, int numOfCoords, int startingNode) {

    int visitedCount = 0;

    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    bool *visited = malloc(numOfCoords * sizeof(bool));

    // Initialising the tour and visited array
    int p =0;
    for(p =0; p< numOfCoords; p++)
    {
     tour[p] = 0;
     visited[p] = false;
    }

    tour[numOfCoords] = 0;


    // Creating a struct as we need to return two values tour and tour size
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
    int nearest;

    for (int i = 0; i < numOfCoords; i++) {
        if (!visited[i]) {
            if (distances[startingNode][i] < minimumDistance) {
                minimumDistance = distances[startingNode][i];
                nearest = i;
            }
        }
    }


    // Adding the second node
    tour[1] = nearest;
    visited[nearest] = true;
    visitedCount++;

    // Adding starting node at the end to complete the partial tour
    tour[2] = startingNode;

    while (visitedCount < numOfCoords)
    {
        double minimumAdditionalCost = DBL_MAX;
        int minimumPosition;
        int minimumUnvisited;
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
        for ( i = 0; i < visitedCount; i++)
        {
            for ( j = 0; j < numOfCoords; j++)
            {
                threadID = omp_get_thread_num();
                if (!visited[j])
                {
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

            if (minimumAdditionalCosts[x] < minimumAdditionalCost) {
                minimumAdditionalCost = minimumAdditionalCosts[x];
                minimumPosition = positions[x];
                minimumUnvisited = nearestNodes[x];
            }
        }

        int indexBefore = minimumPosition == 0 ? visitedCount - 1 : minimumPosition - 1;
        int indexAfter = minimumPosition + 1;

        double distanceAfter =
                distances[tour[minimumPosition]][minimumUnvisited] + distances[tour[indexAfter]][minimumUnvisited] -
                distances[tour[minimumPosition]][tour[indexAfter]];
        double distanceBefore =
                distances[tour[minimumPosition]][minimumUnvisited] + distances[tour[indexBefore]][minimumUnvisited] -
                distances[tour[minimumPosition]][tour[indexBefore]];

        if (distanceAfter < distanceBefore) {
            minimumPosition = indexAfter;
        } else {
            minimumPosition = indexBefore + 1;
        }


        // Making space in the array
        for (i = visitedCount; i > minimumPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[minimumPosition] = minimumUnvisited;
        visited[minimumUnvisited] = true;
        visitedCount++;

    }


    // copy tour to the return object
    for (int i = 0; i <= numOfCoords; i++) {
        result.tour[i] = tour[i];
    }

    // Calculate tour cost
    double totalLength = 0;
    int i =0;
    for (i = 0; i <=numOfCoords; i++) {
        printf("%d ", tour[i]);
        if(i>0) {
            totalLength += distances[tour[i]][tour[i - 1]];
        }
    }

    result.tourSize = totalLength;

    return result;
}


int main(int argc, char *argv[]){

    if(argc != 3){
        printf("Program should be called as ./program <coordFile> <outFileName>");
        return 1;
    }


    //Argument setup for file and output
    char filename[500];
    char outFileName[500];

    strcpy(filename, argv[1]);
    strcpy(outFileName, argv[2]);

    //Reading files and setting up the distance matrix
    int numOfCoords = readNumOfCoords(filename);
    double **coords = readCoords(filename, numOfCoords);
    double shortestTour = DBL_MAX;
    int *shortestTourArray =  (int *)malloc((numOfCoords+1) * sizeof(int *));


    double tStart = omp_get_wtime();

    /*Program starts*/

    double **distances = createDistanceMatrix(coords, numOfCoords);

    for(int top = 0; top<numOfCoords; top++) {

        struct TourData tempTour  = nearestAddition(distances, numOfCoords, top);
        int currentTour = tempTour.tourSize;

        if(currentTour < shortestTour)
        {
            printf("Shortest tour: %f\n", shortestTour);
            printf("Current Tour: %f\n", currentTour);
            shortestTour = currentTour;

            printf("Found tour shorter than current tour");
            printf("\n");
            printf("shortest tour now starting with %d", tempTour.tour[0]);
            int copy =0;
            for(copy =0; copy <numOfCoords+1; copy++)
            {

                shortestTourArray[copy] = tempTour.tour[copy];
            }
        }
        else{
            printf("\n");
            printf("Current tour size:");
            printf("%f", tempTour.tourSize);
            printf("\n");

            printf("Shortest tour size:");
            printf("%f", shortestTour);
            printf("\n");
            printf("Found tour longer than current tour");
            printf("\n");
        }

    }

    /*Program ends*/

    double tEnd = omp_get_wtime();

    printf("\nTook %f milliseconds", (tEnd - tStart) * 1000);

    if (writeTourToFile(shortestTourArray, numOfCoords + 1, outFileName) == NULL){
        printf("Error");
    }

    //Free memory
    for(int i = 0; i < numOfCoords; i++){
        free(distances[i]);
    }

    free(distances);
//	free(tour);
}
