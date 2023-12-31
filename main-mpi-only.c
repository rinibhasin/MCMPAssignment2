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
struct TourData cheapestInsertion(double **distanceMatrix, int numOfCoords, int top);

//struct TourData cheapestInsertion(double **dMatrix, int numOfCoords, int top);
//struct TourData nearestAddition(double **distances, int numOfCoords, int startingNode);


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
    if (myRank == 0)
    {
        numOfCoords = readNumOfCoords(filename);
        double **coords = readCoords(filename, numOfCoords);

        distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
        int i = 0;
        for (i = 0; i < numOfCoords; i++) {
            distanceMatrix[i] = (double *) malloc(numOfCoords * sizeof(double));
        }
        distanceMatrix = createDistanceMatrix(coords, numOfCoords);
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

    MPI_Bcast(&numOfCoords, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate memory for flattenedDistanceMatrix on all MPI processes except the root node
    if (myRank != 0)
    {
        flattenedDistanceMatrix = (double *)malloc(numOfCoords * numOfCoords * sizeof(double));
    }

    // Broadcast the value of numOfCoords to all processes
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


    printf("Going to have a nervous breakdown");

    printf("Comm size : %d\n", commSize);
    printf("Coordinates per process= %d\n", coordsPerProcess);
    printf("Starting coordinate  : %d\n", startingCoord);
    printf("Ending coordinate : %d\n", endingCoord);


    // Variables for storing the shortest tour for each implementations and the corrresponding tour arrays
    double shortestTourFarthest = DBL_MAX;
    double shortestTourCheapest = DBL_MAX;

    int *shortestTourArrayFarthest =  (int *)malloc((numOfCoords+1) * sizeof(int *));
    int *shortestTourArrayCheapest =  (int *)malloc((numOfCoords+1) * sizeof(int *));


    double tStart = omp_get_wtime();

    double start = MPI_Wtime();

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
        printf("Rank %d: Finished processing starting coordinate = %d\n", myRank, top);


        struct TourData tempTourCheapest  = cheapestInsertion(distanceMatrix, numOfCoords, top);
        int currentTourCheapest = tempTourCheapest.tourSize;

        if(currentTourCheapest < shortestTourCheapest)
        {

            shortestTourCheapest = currentTourCheapest;
            // Copying the array to keep track to write to output file later
            int copy=0;
            for(copy =0; copy <numOfCoords+1; copy++)
            {

                shortestTourArrayCheapest[copy] = tempTourCheapest.tour[copy];
            }
        }

        printf(" Total cost = %f \n",tempTourFarthest.tourSize);
        printf("Rank %d: Finished processing starting coordinate = %d\n", myRank, top);
    }


    int *gatheredToursFarthest = NULL;
    double *gatheredTourCostsFarthest = NULL;

    int *gatheredToursCheapest = NULL;
    double *gatheredTourCostsCheapest = NULL;

    if (myRank == 0)
    {
        gatheredToursFarthest = (int *)malloc(commSize * (numOfCoords + 1) * sizeof(int));
        gatheredTourCostsFarthest = (double *)malloc(commSize * sizeof(double));

        gatheredToursCheapest = (int *)malloc(commSize * (numOfCoords + 1) * sizeof(int));
        gatheredTourCostsCheapest = (double *)malloc(commSize * sizeof(double));
    }

    MPI_Gather(shortestTourArrayFarthest, numOfCoords + 1, MPI_INT, gatheredToursFarthest, numOfCoords + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&shortestTourFarthest, 1, MPI_DOUBLE, gatheredTourCostsFarthest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    MPI_Gather(shortestTourArrayCheapest, numOfCoords + 1, MPI_INT, gatheredToursCheapest, numOfCoords + 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&shortestTourCheapest, 1, MPI_DOUBLE, gatheredTourCostsCheapest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        int processId=0, tourIdFarthest=0;
        int tourIdCheapest=0;
        double minimumCostFarthest = DBL_MAX;
        double minimumCostCheapest = DBL_MAX;

        int **finalResultFarthest = (int **)malloc(commSize * sizeof(int *));
        int **finalResultCheapest = (int **)malloc(commSize * sizeof(int *));


        for (processId = 0; processId < commSize; processId++)
        {
            finalResultFarthest[processId] = (int *)malloc((numOfCoords + 1) * sizeof(int));
            finalResultCheapest[processId] = (int *)malloc((numOfCoords + 1) * sizeof(int));

            printf("Tour from process %d: ", processId);
            int i;
            for (i = 0; i <= numOfCoords; i++)
            {
                finalResultFarthest[processId][i] = gatheredToursFarthest[processId * (numOfCoords + 1) + i];
                finalResultCheapest[processId][i] = gatheredToursCheapest[processId * (numOfCoords + 1) + i];

                printf("%d ", finalResultFarthest[processId][i]);
            }

            if (gatheredTourCostsFarthest[processId] < minimumCostFarthest ||
                (gatheredTourCostsFarthest[processId] == minimumCostFarthest &&
                finalResultFarthest[processId][0] < finalResultFarthest[tourIdFarthest][0]))
            {
                minimumCostFarthest = gatheredTourCostsFarthest[processId];
                tourIdFarthest = processId;
            }

            if (gatheredToursCheapest[processId] < minimumCostCheapest ||
                (gatheredTourCostsCheapest[processId] == minimumCostCheapest &&
                 finalResultCheapest[processId][0] < finalResultCheapest[tourIdCheapest][0]))
            {
                minimumCostCheapest = gatheredTourCostsCheapest[processId];
                tourIdCheapest = processId;
            }
        }

        printf("Tour id farthest : %d\n", tourIdFarthest);
        printf("Cost: %f\n", minimumCostFarthest);
        int i;
        for (i = 0; i <= numOfCoords; i++)
        {
            printf("%d ", finalResultFarthest[tourIdFarthest][i]);
            printf("%d ", finalResultCheapest[tourIdCheapest][i]);
        }

        writeTourToFile(finalResultFarthest[tourIdFarthest], numOfCoords+1 , outFileName1);
        writeTourToFile(finalResultCheapest[tourIdCheapest], numOfCoords+1 , outFileName2);

        for (processId = 0; processId < commSize; processId++)
        {
            free(finalResultFarthest[processId]);
            free(finalResultCheapest[processId]);
        }

        free(finalResultFarthest);
        free(finalResultCheapest);

    }


    double tEnd = omp_get_wtime();

    double end = MPI_Wtime();

    printf("\nTook %f milliseconds", (tEnd - tStart) * 1000);
    printf("\nTook %f seconds MPI time", (end - start));
    printf("Writing tour to file farthest %s\n", outFileName1);


//    if (writeTourToFile(shortestTourArrayFarthest, numOfCoords + 1, outFileName1) == NULL){
//        printf("Error");
//    }


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


