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

struct InsertionResult {
    int *tour;
    double tourcost;
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



// Function to solve the TSP using Nearest addition method
struct InsertionResult nearestAddition_TSP(double **distances, int numOfCoords, int startingNode) {

    int visitedCount = 1;
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
    struct InsertionResult result;
    result.tour = malloc((numOfCoords + 1) * sizeof(int));

    int noOfThreads = omp_get_max_threads();

    double *threads_min_distance = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *nodes = (int *) malloc(noOfThreads * sizeof(int));

    // Step 1 - Start off with a vertex V0
    tour[0] = startingNode;
    visited[startingNode] = true;
    int currentCity = startingNode;
    double minDistance = DBL_MAX;

    // Step 2 - Find a vertex Vi such that dist(V0, Vi) is minimal
    for (int i = 0; i < numOfCoords; i++) {
        if (!visited[i]) {
            if (distances[currentCity][i] < minDistance) {
                minDistance = distances[currentCity][i];
                nearest = i;
            }
        }
    }

    // Add Vi to the tour
    tour[1] = nearest;
    visited[nearest] = true;
    visitedCount++;

    tour[visitedCount] = startingNode;
//    visitedCount++;
// Iterate through the rest of the nodes/coords
    while (visitedCount < numOfCoords) {
        // Initialize variables for finding the next vertex to add
        double minAdditionCost = DBL_MAX;
        double minimum_Cost = DBL_MAX;
        int min_position;
        int min_Unvisited_node;
        int positionToAdd, position;
        int i = 0, y=0;
        double max = DBL_MAX;
        int thread_ID;
//double additionalCost = 0.0;

        for (y = 0; y < noOfThreads; y++) {
            threads_min_distance[y] = DBL_MAX;
            positions[y] = 0;
            nodes[y] = 0;
        }
        int n=0, k=0;
        // Step 3 - For all vertices vn in the partial tour
#pragma omp parallel for collapse(2) private(n, k, thread_ID) shared(visited, distances, threads_min_distance, positions, nodes)
        for ( n = 0; n < visitedCount; n++) {
            for ( k = 0; k < numOfCoords; k++) {
                thread_ID = omp_get_thread_num();
                if (!visited[k]) {

                    double additionalCost = distances[tour[n]][k];
                    if (additionalCost < threads_min_distance[thread_ID]) {
                        threads_min_distance[thread_ID] = additionalCost;
                        positions[thread_ID] = n;
                        nodes[thread_ID] = k;
                    }
                }
            }
        }
        /*
      *  Check through each thread's memory location and update the minimumAdditionalCost, position and
      *  best vertex
      */
        int x = 0;
        for (x = 0; x < noOfThreads; x++) {

            if (threads_min_distance[x] < minimum_Cost) {
                minimum_Cost = threads_min_distance[x];
                min_position = positions[x];
                min_Unvisited_node = nodes[x];
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

// Add the starting vertex at the end to complete the tour

    for (int i = 0; i <= numOfCoords; i++) {
        printf("%d ", tour[i]);
        result.tour[i] = tour[i];
    }
    printf(" \n \n \n ");
//    }

    double cost = 0.0;
// Iterate over the tour to calculate the cost
    for (int i = 0; i < numOfCoords - 1; i++) {
        int fromVertex = tour[i];
        int toVertex = tour[i + 1];

        // Accumulate the distance between consecutive vertices
        cost += distances[fromVertex][toVertex];

    }

// Add the distance from the last vertex back to the starting vertex
    cost += distances[tour[numOfCoords - 1]][tour[0]];

    tour[numOfCoords + 1] = cost;
    result.tourcost = cost;

// Return the tour
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

        struct InsertionResult tempTour  = nearestAddition_TSP(distances, numOfCoords, top);
        int currentTour = tempTour.tourcost;

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
            printf("%f", tempTour.tourcost);
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
