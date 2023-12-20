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

//#pragma omp parallel for collapse(2) private(i, j, x1, y1, x2, y2, distance) shared(numOfCoords)
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



// Function to solve the TSP using Nearest addition method
struct TourData nearestAddition_TSP(double **distances, int numOfCoords, int starting_coord) {

    int currentSize = 1;
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    bool *visited_nodes = malloc(numOfCoords * sizeof(bool));
    int nearest;
    struct TourData result;
    result.tour = malloc((numOfCoords + 1) * sizeof(int));

    /* Get the number of Threads */
    int noOfThreads = omp_get_max_threads();
    /* Set up arrays to reserve memory locations for the threads using malloc
     * Arrays are equal to the noOfThreads allocated for them
     */
    double *threads_min_distance = (double *) malloc(noOfThreads * sizeof(double));
    int *positions = (int *) malloc(noOfThreads * sizeof(int));
    int *nodes = (int *) malloc(noOfThreads * sizeof(int));

    // Step 1 - Start off with a vertex V0
    tour[0] = starting_coord;
    visited_nodes[starting_coord] = true;
    int currentCity = starting_coord;
    double minDistance = DBL_MAX;

    // Step 2 - Find a vertex Vi such that dist(V0, Vi) is minimal
    for (int i = 0; i < numOfCoords; i++) {
        if (!visited_nodes[i]) {
            if (distances[currentCity][i] < minDistance) {
                minDistance = distances[currentCity][i];
                nearest = i;
            }
        }
    }

    // Add Vi to the tour
    tour[1] = nearest;
    visited_nodes[nearest] = true;
    currentSize++;

    tour[currentSize] = starting_coord;
//    currentSize++;
// Iterate through the rest of the nodes/coords
    while (currentSize < numOfCoords) {
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
#pragma omp parallel for collapse(2) private(n, k, thread_ID) shared(visited_nodes, distances, threads_min_distance, positions, nodes)
        for ( n = 0; n < currentSize; n++) {
            for ( k = 0; k < numOfCoords; k++) {
                thread_ID = omp_get_thread_num();
                if (!visited_nodes[k]) {

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

        int indexBefore = min_position == 0 ? currentSize - 1 : min_position - 1;
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

        currentSize++;
        visited_nodes[min_Unvisited_node] = true;

        for (i = currentSize; i > min_position; i--) {
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
    result.tourSize = cost;

// Return the tour
    return result;

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
//#pragma omp parallel for private(i)
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);

    double shortestTour = DBL_MAX;
    int *shortestTourArray =  (int *)malloc((numOfCoords+1) * sizeof(int *));

    for(i = 0; i< numOfCoords; i++)
    {
        struct TourData tempTour = nearestAddition_TSP(distanceMatrix, numOfCoords, i);
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
