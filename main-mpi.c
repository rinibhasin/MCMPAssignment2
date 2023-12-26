#include <float.h>
#include <stdlib.h>
#include<stdbool.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

struct TourData {
    int *tour;
    double tourSize;
};

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

double calculateDistance(double x1, double y1, double x2, double y2);


//int *cheapestInsertion_TSP(double **distances, int numOfCoords, char *outputfile, int starting_coord);
//struct InsertionResult cheapestInsertion(double **distanceMatrix, int numOfCoords,char *outputfile,  int starting_coord);
//struct TourData (double **distances, int numOfCoords, char *outputfile, int starting_coord);

//struct InsertionResult farthestInsertion_TSP(double **distances, int numOfCoords, char *outputfile, int starting_coord);


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

int main(int argc, char *argv[])
{
    int myRank, commSize;
    double start, end;
    double cpu_time_used;
    int* nodes = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    char *coordinatefile = argv[1];
    char *outputfile = argv[2];
    double **distanceMatrix;
    double **coordinates;
    int **finalarrays;
    double *flattenedDistanceMatrix = NULL;
    int  i=0;
    int numOfCoords = 0;


    // Only one the root process will read input file and send data to other MPI processes
    if (myRank == 0) {

        printf("Rank : %d, Reading coordinates", myRank);
        numOfCoords = readNumOfCoords(coordinatefile);
        coordinates = readCoords(coordinatefile, numOfCoords);
        nodes = (int*)malloc(numOfCoords * sizeof(int));
        int i;
        for ( i = 0; i < numOfCoords; i++) {
            nodes[i] = i;
        }
        distanceMatrix = (double **) malloc(numOfCoords * sizeof(double *));

        finalarrays = (int **) malloc((numOfCoords + 1) * sizeof(int *));

        for (i = 0; i < numOfCoords + 1; i++) {
            finalarrays[i] = (int *) malloc((numOfCoords + 1) * sizeof(int *));
        }

        for (i = 0; i < numOfCoords; i++) {
            distanceMatrix[i] = (double *) malloc(numOfCoords * sizeof(double));
        }

        // Calculate the Distance matrix
        distanceMatrix = createDistanceMatrix(coordinates, numOfCoords);
        flattenedDistanceMatrix = (double *)malloc(numOfCoords * numOfCoords * sizeof(double));
        int count = 0, j=0;
        for (i = 0; i < numOfCoords; i++)
        {
            for (j = 0; j < numOfCoords; j++)
            {
                flattenedDistanceMatrix[count++] = distanceMatrix[i][j];
            }
        }
    }


// Allocate memory for flattenedDistanceMatrix on all MPI processes excpet the root node
// flattenedDistanceMatrix = (double *)malloc(numOfCoords * numOfCoords * sizeof(double));

    if (myRank != 0) {
        flattenedDistanceMatrix = (double *)malloc(numOfCoords * numOfCoords * sizeof(double));
    }

    // Broadcast the value of numOfCoords to all processes
    MPI_Bcast(&numOfCoords, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(flattenedDistanceMatrix, numOfCoords * numOfCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // Reshape the distance matrix
    if (myRank != 0) {
        distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
        for (i = 0; i < numOfCoords; i++) {
            distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
        }
        int index = 0,j=0;
        for (i = 0; i < numOfCoords; i++) {
            for (j = 0; j < numOfCoords; j++) {
                distanceMatrix[i][j] = flattenedDistanceMatrix[index++];
                printf("rescaled matrix : %f\n", distanceMatrix[i][j]);
            }
        }
    }

    printf("Rank : %d\n",myRank);
    // Scatter the points to all processes
    int coordsPerProcess = (numOfCoords + commSize - 1) / commSize; // Ceiling division

    // Dividing the coordinates between different MPI processes
    int startingCoord = myRank * coordsPerProcess;
    int endingCoord = (myRank + 1) * coordsPerProcess;
    endingCoord = (endingCoord > numOfCoords) ? numOfCoords : endingCoord;
    float *localNodes = (float*)malloc(coordsPerProcess * sizeof(float));
    printf("comm_size : %d\n", commSize);
    printf("coordsPerProcess= %d\n", coordsPerProcess);
    printf("startingCoord  : %d\n", startingCoord);
    printf("endingCoord : %d\n", endingCoord);
    double max_cost = DBL_MAX;
    double local_min_cost = DBL_MAX;
    int r = 0;
    // Each process independently calculates its part of the tour using nearest addition
    for ( r = startingCoord; r < endingCoord; r++) {
        int counter = 0, j;
        int rowid = 0;
        int *tour = (int *) malloc((numOfCoords + 1) * sizeof(int));
        struct InsertionResult result;
        int i;
        printf("Rank %d: Processing start = %d\n", myRank, r);
        // Check if coordinates and distanceMatrix are not NULL
        if ( distanceMatrix == NULL) {
            printf("Rank: %d, Coordinates or distanceMatrix is NULL\n", myRank);
            //      MPI_Abort(MPI_COMM_WORLD, 1);  // Terminate MPI execution
        }
        result = nearestAddition_TSP(distanceMatrix, numOfCoords, outputfile, r);

//        if(r == 1){
//            writeTourToFile(result.tour, numOfCoords+1, outputfile);
//        }
        double currentTourCost = result.tourcost;
        if(currentTourCost < local_min_cost){
            local_min_cost = currentTourCost;
        }

        printf(" Total cost = %f \n",result.tourcost);
        printf("Rank %d: Finished processing start = %d\n", myRank, r);
    }
// Combine results from different MPI processes and find the global minimum
    double global_min;
    printf(" Total cost before reduce  = %f \n",local_min_cost);
    MPI_Reduce(&local_min_cost, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    // Print the result
    if (myRank== 0) {
        printf("Global minimum tour cost: %f\n", global_min);
        for (i = 0; i < numOfCoords; i++) {
            free(distanceMatrix[i]);
        }
        for (i = 0; i < numOfCoords; i++) {
            free(coordinates[i]);
        }

        for (i = 0; i <= numOfCoords; i++) {
            free(finalarrays[i]);
        }
        free(flattenedDistanceMatrix);
    }
    MPI_Finalize();
}

