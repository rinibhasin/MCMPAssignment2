#include <float.h>
#include <stdlib.h>
#include<stdbool.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

void farthestInsertion(double **distanceMatrix, int totalCoordinates, char *outputfile) {

    int visitedCount = 0;

    int *tour = malloc((totalCoordinates+1) * sizeof(int));
    bool *visited = malloc(totalCoordinates * sizeof(bool));

    // Initialise with the first vertex
    tour[0] = 0;
    visited[0] = true;
    visitedCount++;

    // Find the farthest vertex
    double maximumDistance = DBL_MIN;

    int farthestVertex;
    int i = 0;
    for(i = 1 ; i <totalCoordinates; i++)
    {
        if(distanceMatrix[0][i]> maximumDistance)
        {
            maximumDistance = distanceMatrix[0][i];
            farthestVertex = i;
        }
    }

    tour[1] = farthestVertex;
    visited[farthestVertex] = true;
    visitedCount++; //2
    tour[2] = 0;

    int maxThreads = omp_get_max_threads();

    double *farthestDistances = (double*)malloc(maxThreads*sizeof(double));
    int *positions = (int*)malloc(maxThreads*sizeof(int));
    int *farthestNodes = (int*)malloc(maxThreads*sizeof(int));

    double minimumDistance = DBL_MAX;

    while (visitedCount < totalCoordinates) {
        double farthestDistance = 0; // unvisited max
        int farthestNode = 0;
        int j, x;
        int threadID;
        for(x = 0; x < maxThreads; x++)
        {
            farthestDistances[x] = 0.0;
            positions[x] =0;
            farthestNodes[x]=0;
        }

#pragma omp parallel for collapse(2) private(i, j, threadID) shared(visited, distanceMatrix, farthestDistances, positions, farthestNodes)
        for (i = 0; i < visitedCount; i++)
        {
            for (j = 0; j < totalCoordinates; j++)
            {
                threadID = omp_get_thread_num();
                if (!visited[j])
                {
                    double currentDistance = distanceMatrix[j][tour[i]];
                    if (currentDistance >  farthestDistances[threadID])
                    {
                        farthestDistances[threadID] = distanceMatrix[tour[i]][j];
                        farthestNodes[threadID] = j;
                    }
                }
            }
        }

        int m=0;
        for(m =0; m< totalCoordinates; m++)
        {
            if(farthestDistances[m] > farthestDistance)
            {
                farthestDistance = farthestDistances[m];
                farthestNode = farthestNodes[m];
            }
        }

        double minimumAdditionalCost = DBL_MAX;
        int minN;

        for (i = 0; i < visitedCount; i++)
        {
            double additionalCost = distanceMatrix[farthestNode][tour[i]]+ distanceMatrix[farthestNode][tour[i+1]] - distanceMatrix[tour[i]][tour[i+1]];

            if (additionalCost < minimumAdditionalCost)
            {
                minimumAdditionalCost = additionalCost;
                minN = i; // where to insert
                additionalCost = additionalCost;
            }
        }

        // Make space to add unvisited node to computed index
        for(i = visitedCount; i > minN; i--)
        {
            tour[i+1] = tour[i];
        }

        // Add the node to the tour
        tour[minN+1] = farthestNode;
        visited[farthestNode] = true;
        visitedCount++;
    }

    int tourLength = visitedCount+1;
    writeTourToFile(tour, tourLength, outputfile);
}

double calculateDistance(double x1, double y1, double x2, double y2) {
    return sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}


double **calculateDistanceMatrix(double **coordinates, int numOfCoords, double **distanceMatrix) {
    int i, j;
    double x1, x2, y1, y2;
    double distance;

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

    // taking default file names if user didn't provide input
    char *fileName = "9_coords.coord";
    char *outputfile = "output.txt";

    if (argc > 1) {
        fileName = argv[1];
        outputfile = argv[2];

    }

    double start, end;
    double time_taken;

    int numOfCoords = readNumOfCoords(fileName);
    double **coordinates = readCoords(fileName, numOfCoords);


    double **distanceMatrix = (double **) malloc(numOfCoords * sizeof(double *));
    int i;
#pragma omp parallel for private(i)
    for (i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *) malloc(numOfCoords * sizeof(double));
    }

    start = omp_get_wtime();;
    distanceMatrix = calculateDistanceMatrix(coordinates, numOfCoords, distanceMatrix);
    farthestInsertion(distanceMatrix, numOfCoords, outputfile);
    end = omp_get_wtime();

    time_taken = (end - start);

    printf("The time taken is %fs\n", time_taken);
}



