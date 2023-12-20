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

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);
struct TourData {
    int *tour;
    double tourSize;
};


struct TourData nearestAddition(double **distanceMatrix, int numOfCoords, char *outputfile, int startingNode) {

    int visitedCount = 0;
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    bool *visited = malloc(numOfCoords * sizeof(bool));

    int nearest;
    struct TourData result;
    result.tour = malloc((numOfCoords + 1) * sizeof(int));

    int noOfThreads = omp_get_max_threads();
   
    double *minimumAdditionalCosts = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *nearestVertexes = (int *) malloc(noOfThreads * sizeof(int));


    tour[0] = startingNode;
    visited[startingNode] = true;
    visitedCount++;
    
    int currentCity = startingNode;
    double minimumDistance = DBL_MAX;
    int i = 0;

    for (i = 0; i < numOfCoords; i++)
    {
        if (!visited[i])
        {
            if (distanceMatrix[currentCity][i] < minimumDistance)
            {
                minimumDistance = distanceMatrix[currentCity][i];
                nearest = i;
            }
        }
    }

    tour[1] = nearest;
    visited[nearest] = true;
    visitedCount++;

    // complete the partial tour
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
            nearestVertexes[y] = 0;
        }
        int i=0, j=0;
        #pragma omp parallel for collapse(2) private(i, j, threadID) shared(visited, distanceMatrix, minimumAdditionalCosts, positions, nearestVertexes)
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
                        nearestVertexes[threadID] = j;
                    }
                }
            }
        }

        int x = 0;
        for (x = 0; x < noOfThreads; x++) {

            if (minimumAdditionalCosts[x] < minimum_Cost) {
                minimum_Cost = minimumAdditionalCosts[x];
                min_position = positions[x];
                min_Unvisited_node = nearestVertexes[x];
            }
        }

        int indexBefore = min_position == 0 ? visitedCount - 1 : min_position - 1;
        int indexAfter = min_position + 1;

        double distanceAfter =
                distanceMatrix[tour[min_position]][min_Unvisited_node] + distanceMatrix[tour[indexAfter]][min_Unvisited_node] -
                distanceMatrix[tour[min_position]][tour[indexAfter]];
        double distanceBefore =
                distanceMatrix[tour[min_position]][min_Unvisited_node] + distanceMatrix[tour[indexBefore]][min_Unvisited_node] -
                distanceMatrix[tour[min_position]][tour[indexBefore]];

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
    printf(" \n \n \n ");
//    }

    double tourCost = 0.0;
// Iterate over the tour to calculate the cost
    for (int i = 0; i < numOfCoords - 1; i++) {
        int fromVertex = tour[i];
        int toVertex = tour[i + 1];

        // Accumulate the distance between consecutive vertices
        cost += distanceMatrix[fromVertex][toVertex];

    }

// Add the distance from the last vertex back to the starting vertex
    cost += distanceMatrix[tour[numOfCoords - 1]][tour[0]];

    tour[numOfCoords + 1] = cost;
    result.tourSize = cost;

// Return the tour
    return result;

}

double calculateDistance(double x1, double y1, double x2, double y2) {
    // Calculate the Euclidean distance between two points.
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}

double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {

    int i =0; int j =0;
    double x1 =0;
    double x2 =0;
    double y1 =0;
    double y2 =0;
    double distance = 0;

#pragma omp parallel for collapse(2) private(i, j, x1, y1, x2, y2, distance) shared(numOfCoords)
    for (i = 0; i < numOfCoords; i++) {
        for (j = 0; j < numOfCoords; j++) {

            x1 = coordinates[i][0];
            y1 = coordinates[i][1];
            x2 = coordinates[j][0];
            y2 = coordinates[j][1];

            distance = calculateDistance(x1, y1, x2, y2);

            distanceMatrix[i][j] = distance;

        }
    }

    return distanceMatrix;
}

int main(int argc, char *argv[]) {

    // Taking default file names if user didn't provide input
    char *fileName = "512_coords.coord";
    char *outputfile = "output2.txt";

    if (argc > 1) {
        fileName = argv[1];
        outputfile = argv[2];
    }

    double start, end;
    double time_taken;

//    start = omp_get_wtime();;

    int numOfCoords = readNumOfCoords(fileName);

    double **coordinates = readCoords(fileName, numOfCoords);

    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));

    int i = 0;
#pragma omp parallel for private(i)
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);

    double shortestTour = DBL_MAX;

    int *shortestTourArray =  (int *)malloc((numOfCoords+1) * sizeof(int *));

    for(i = 0; i< 512; i++)
    {
        struct TourData tempTour = nearestAddition_TSP(distanceMatrix, numOfCoords, outputfile, i);
        int currentTour = tempTour.tourSize;

        if(currentTour < shortestTour)
        {
            shortestTour = currentTour;
            printf("\n");
            printf("Found tour shorter than current tour");

            printf("shortest tour now starting with %d", tempTour.tour[0]);
            //make copy
            int j =0;
            for(j =0; j <numOfCoords+1; j++)
            {
                shortestTourArray[j] = tempTour.tour[j];
            }
        }

        else
        {
            printf("Current tour size:");
            printf("%f", tempTour.tourSize);
            printf("\n");

            printf("Shortes tour size:");
            printf("%f", shortestTour);
            printf("\n");
            printf("Found tour longer than current tour");
            printf("");
        }
    }

    printf("Writing the shortest tour to file\n: ");

    printf("Starting with %d: ", shortestTourArray[0]);

    writeTourToFile(shortestTourArray, numOfCoords+1, outputfile);

//    end = omp_get_wtime();
//    time_taken = (end - start);

    printf("Nearest insertion for %d nearestVertexes completed\n", numOfCoords);
//    printf("The time taken is %fs .\n", time_taken);

    // Free memory
    for (i = 0; i < numOfCoords; i++) {
        free(coordinates[i]);
    }
    free(coordinates);

    for (i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
//    cleanupStruct(&tempTour);
    return 0;
}
