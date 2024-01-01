#include <float.h>
#include <stdlib.h>
#include<stdbool.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>


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


    //Argument setup for file and output
    char *filename;
    char *outFileName1;
    char *outFileName2;
    char *outFileName3;

    filename = argv[1];
    outFileName1 = argv[2];
    outFileName2 = argv[3];
    outFileName3 = argv[4];


    //Reading files and setting up the distance matrix
    int numOfCoords = readNumOfCoords(filename);
    double **coords = readCoords(filename, numOfCoords);


    // Variables for storing the shortest tour for each implementations and the corrresponding tour arrays
    double shortestTourFarthest = DBL_MAX;
    int *shortestTourArrayFarthest =  (int *)malloc((numOfCoords+1) * sizeof(int *));

    double shortestTourCheapest = DBL_MAX;
    int *shortestTourArrayCheapest =  (int *)malloc((numOfCoords+1) * sizeof(int *));

    double shortestTourNearest = DBL_MAX;
    int *shortestTourArrayNearest =  (int *)malloc((numOfCoords+1) * sizeof(int *));

    double tStart = omp_get_wtime();

    double **distances = createDistanceMatrix(coords, numOfCoords);

    for(int top = 0; top<numOfCoords; top++)
    {
        struct TourData tempTourFarthest  = farthestInsertion(distances, numOfCoords, top);
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


        struct TourData tempTourCheapest  = cheapestInsertion(distances, numOfCoords, top);
        int currentTourCheapest = tempTourCheapest.tourSize;

        if(currentTourCheapest < shortestTourCheapest)
        {

            shortestTourCheapest = currentTourCheapest;
            printf("Found tour shorter than current tour");
            printf("\n");
            printf("Shortest tour now starting with %d", tempTourCheapest.tour[0]);
            // Copying the array to keep track to write to output file later
            int copy=0;
            for(copy =0; copy <numOfCoords+1; copy++)
            {
                shortestTourArrayCheapest[copy] = tempTourCheapest.tour[copy];
            }
        }

        struct TourData tempTourNearest  = nearestAddition(distances, numOfCoords, top);
        int currentTourNearest = tempTourNearest.tourSize;

        if(currentTourNearest < shortestTourNearest)
        {
//            printf("Nearest");
//            printf("Shortest tour: %f\n", shortestTourNearest);
//            printf("Current Tour: %f\n", currentTourNearest);
            shortestTourNearest = currentTourNearest;
//            printf("Found tour shorter than current tour");
//            printf("\n");
//            printf("Shortest tour now starting with %d", tempTourCheapest.tour[0]);
            // Copying the array to keep track to write to output file later
            int copy=0;
            for(copy =0; copy <numOfCoords+1; copy++)
            {
                shortestTourArrayNearest[copy] = tempTourNearest.tour[copy];
            }
        }
    }

    double tEnd = omp_get_wtime();

    printf("\nTook %f milliseconds", (tEnd - tStart) * 1000);
    printf("Writing tour to file farthest %s\n", outFileName1);

    printf("\n");
    printf("Shortest tour Farthest %f", shortestTourArrayFarthest);

    printf("\n");
    printf("Shortest tour Cheapest %f", shortestTourArrayCheapest);

    printf("\n");
    printf("Shortest tour Cheapest %f", shortestTourArrayNearest);



    printf("Writing tour to file cheapest%s\n", outFileName2);


    if (writeTourToFile(shortestTourArrayCheapest, numOfCoords + 1, outFileName1) == NULL){
        printf("Error");
    }


    if (writeTourToFile(shortestTourArrayFarthest, numOfCoords + 1, outFileName2) == NULL){
        printf("Error");
    }

    printf("Writing tour to file nearest %s\n", outFileName3);


    if (writeTourToFile(shortestTourArrayNearest, numOfCoords + 1, outFileName3) == NULL){
        printf("Error");
    }

//    Free memory
    for(int i = 0; i < numOfCoords; i++){
        free(distances[i]);
    }

    free(distances);
    free(shortestTourArrayCheapest);
//    free(shortestTourArrayNearest);
    free(shortestTourArrayFarthest);

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


