#include <float.h>
#include <stdlib.h>
#include<stdbool.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>


struct TourData {
    int* tour;
    double tourSize;
};

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);
double **createDistanceMatrix(double **coords, int numOfCoords);
double sqrt(double arg);
struct TourData farthestInsertion(double **dMatrix, int numOfCoords, int top);
struct TourData cheapestInsertion(double **dMatrix, int numOfCoords, int top);
struct TourData nearestAddition(double **distances, int numOfCoords, int startingNode);


int main(int argc, char *argv[]){

    // MPI variables setup
    int myRank, commSize;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    //Argument setup for file and output
    char *filename;
    char *outFileName1;
    char *outFileName2;
    char *outFileName3;

    filename = argv[1];
    outFileName1 = argv[2];
    outFileName2 = argv[3];
    outFileName3 = argv[4];


    int numOfCoords = 0;

    double **distanceMatrix;

    double *flattenedDistanceMatrix = NULL;

    // Only one the root process will read input file and send data to other MPI processes
    if (myRank == 0) {

        numOfCoords = readNumOfCoords(filename);
        double **coords = readCoords(filename, numOfCoords);
        distanceMatrix = createDistanceMatrix(coords, numOfCoords);
        flattenedDistanceMatrix = (double *)malloc(numOfCoords * numOfCoords * sizeof(double));

        int count = 0, j=0;
        int i;
        for (i = 0; i < numOfCoords; i++)
        {
            for (j = 0; j < numOfCoords; j++)
            {
                flattenedDistanceMatrix[count++] = distanceMatrix[i][j];
            }
        }


    }

    // Allocate memory for flattenedDistanceMatrix on all MPI processes except the root node
    if (myRank != 0) {
        flattenedDistanceMatrix = (double *)malloc(numOfCoords * numOfCoords * sizeof(double));
    }

    // Broadcast the value of numOfCoords to all processes
    MPI_Bcast(&numOfCoords, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(flattenedDistanceMatrix, numOfCoords * numOfCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // Reshape the distance matrix
    if (myRank != 0)
    {
        distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
        int i;
        for (i = 0; i < numOfCoords; i++)
        {
            distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
        }
        int index = 0,j=0;
        for (i = 0; i < numOfCoords; i++)
        {
            for (j = 0; j < numOfCoords; j++)
            {
                distanceMatrix[i][j] = flattenedDistanceMatrix[index++];
                printf("Reshaped matrix : %f\n", distanceMatrix[i][j]);
            }
        }
    }

    int coordsPerProcess = (numOfCoords + commSize - 1) / commSize; // Ceiling division

    // Dividing the coordinates between different MPI processes
    int startingCoord = myRank * coordsPerProcess;
    int endingCoord = (myRank + 1) * coordsPerProcess;
    endingCoord = (endingCoord > numOfCoords) ? numOfCoords : endingCoord;

    printf("Commsize : %d\n", commSize);
    printf("Coordinates per process= %d\n", coordsPerProcess);
    printf("Starting coordinate  : %d\n", startingCoord);
    printf("Ending coordinate : %d\n", endingCoord);

    //Reading files and setting up the distance matrix


    // Variables for storing the shortest tour for each implementations and the corrresponding tour arrays
    double shortestTourFarthest = DBL_MAX;
    int *shortestTourArrayFarthest =  (int *)malloc((numOfCoords+1) * sizeof(int *));


    double tStart = omp_get_wtime();

    int top;

    for(top = startingCoord; top<endingCoord; top++)
    {

        struct TourData tempTourFarthest  = farthestInsertion(distanceMatrix, numOfCoords, top);
        int currentTourFarthest = tempTourFarthest.tourSize;

        if(currentTourFarthest < shortestTourFarthest)
        {

            shortestTourFarthest = currentTourFarthest;
            // Copying the array to keep track to write to output file later
            int copy=0;
            for(copy =0; copy <numOfCoords+1; copy++)
            {

                shortestTourArrayFarthest[copy] = tempTourFarthest.tour[copy];
            }
        }

        printf(" Total cost = %f \n",tempTourFarthest.tourSize);
        printf("Rank %d: Finished processing starting corrdinate = %d\n", myRank, top);


    }

    double global_min;
    printf(" Total cost before reduce  = %f \n", shortestTourFarthest);
    MPI_Reduce(&shortestTourFarthest, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);


    if (myRank == 0)
    {
        printf("Global minimum tour cost: %f\n", global_min);
        for (i = 0; i < numOfCoords; i++) {
            free(distanceMatrix[i]);
        }
//        for (i = 0; i < numOfCoords; i++) {
//            free(coordinates[i]);
//        }


        free(flattenedDistanceMatrix);
    }

    double tEnd = omp_get_wtime();

    printf("\nTook %f milliseconds", (tEnd - tStart) * 1000);
    printf("Writing tour to file farthest %s\n", outFileName1);


    if (writeTourToFile(shortestTourArrayFarthest, numOfCoords + 1, outFileName1) == NULL){
        printf("Error");
    }







//    Free memory

//    for(int i = 0; i < numOfCoords; i++){
//        free(distanceMatrix[i]);
//    }

    free(distances);
//    free(shortestTourArrayCheapest);
//    free(shortestTourArrayNearest);
//    free(shortestTourArrayFarthest);

}

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


