#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <omp.h>
#include <pthread.h>

#include "shared.h"
#include "alert.h"
#include "utils.h"
#include "balloon.h"
#include "node.h"

/*
Alert thread to receive all alerts from ocean nodes and log alert into the log file.
*/
void *alertThread(void *args)
{
    // printf("[Alert Thread successfully created]\n");

    AlertThread *data = (AlertThread *)args;
    MPI_Comm comm = data->comm;
    int *iteration = data->iteration;
    double threshold_diff_magnitude = data->threshold_diff_magnitude;
    double threshold_diff_distance = data->threshold_diff_distance;

    // Find total buffer size for MPI_Pack
    int sizeFloat,
        sizeTime, sizeCoord, totalSize;
    MPI_Pack_size(4, MPI_FLOAT, comm, &sizeFloat);
    MPI_Pack_size(40, MPI_CHAR, comm, &sizeTime);
    MPI_Pack_size(2, MPI_INT, comm, &sizeCoord);
    totalSize = sizeFloat + sizeTime + sizeCoord;
    int size = 2 * sizeof(int) + 5 * totalSize;
    char packBufferRecv[size];

    // Continue only if base station has not send termination signal
    while (*(data->sentinel) != SENTINEL)
    {
        MPI_Status probe_status;
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 1, comm, &flag, &probe_status);

        // Continue only if alert actually receive message from an ocean node
        if (flag == 1)
        {
            MPI_Recv(packBufferRecv, size, MPI_PACKED, MPI_ANY_SOURCE, 1, comm, NULL);

            // Unpack message received - format [number of neighbors] [number of matching neighbors] [source node reading] [adjacent node readings]
            int pos = 0;
            int matchings, totalNeighbor;
            int sourceCoord[2];
            int neighborCoord[8];

            // Total neighbor, matchings
            MPI_Unpack(packBufferRecv, size, &pos, &totalNeighbor, 1, MPI_INT, comm);
            MPI_Unpack(packBufferRecv, size, &pos, &matchings, 1, MPI_INT, comm);
            // printf("BASE - received alert from node with %d matchings out of %d neighbors\n", matchings, totalNeighbor);

            SensorReading nodeReading;
            SensorReading neighborReading[totalNeighbor];

            // Source node
            MPI_Unpack(packBufferRecv, size, &pos, &nodeReading.magnitude, 1, MPI_FLOAT, comm);
            MPI_Unpack(packBufferRecv, size, &pos, &nodeReading.latitude, 1, MPI_FLOAT, comm);
            MPI_Unpack(packBufferRecv, size, &pos, &nodeReading.longitude, 1, MPI_FLOAT, comm);
            MPI_Unpack(packBufferRecv, size, &pos, &nodeReading.depth, 1, MPI_FLOAT, comm);
            MPI_Unpack(packBufferRecv, size, &pos, &nodeReading.datetime, 40, MPI_CHAR, comm);
            MPI_Unpack(packBufferRecv, size, &pos, sourceCoord, 2, MPI_INT, comm);
            // printf("BASE - received from source node, mag: %f, lat: %f, long: %f, depth: %f\n", nodeReading.magnitude, nodeReading.latitude, nodeReading.longitude, nodeReading.depth);

            // Adjacent nodes
            for (int j = 0; j < totalNeighbor; j++)
            {
                char bufferSingle[totalSize];
                MPI_Unpack(packBufferRecv, size, &pos, bufferSingle, totalSize, MPI_CHAR, comm);

                int pos1 = 0;

                // printf("BASE - trying to read adj node %d data\n", j);
                MPI_Unpack(bufferSingle, size, &pos1, &neighborReading[j].magnitude, 1, MPI_FLOAT, comm);
                MPI_Unpack(bufferSingle, size, &pos1, &neighborReading[j].latitude, 1, MPI_FLOAT, comm);
                MPI_Unpack(bufferSingle, size, &pos1, &neighborReading[j].longitude, 1, MPI_FLOAT, comm);
                MPI_Unpack(bufferSingle, size, &pos1, &neighborReading[j].depth, 1, MPI_FLOAT, comm);
                MPI_Unpack(bufferSingle, size, &pos1, &neighborReading[j].datetime, 40, MPI_CHAR, comm);
                MPI_Unpack(bufferSingle, size, &pos1, &neighborCoord[2 * j], 2, MPI_INT, comm);
                // printf("BASE - received from adj node %d , mag: %f, lat: %f, long: %f, depth: %f\n", j, neighborReading[j].magnitude, neighborReading[j].latitude, neighborReading[j].longitude, neighborReading[j].depth);
            }
            // printf("BASE - read %d neighbours, onto compare balloon to reporting\n", totalNeighbor);
            compareBalloonToReportingNodeReadings(iteration, matchings, totalNeighbor, &nodeReading, neighborReading, sourceCoord, neighborCoord, threshold_diff_magnitude, threshold_diff_distance);
        }
    }
    // printf("[Alert Thread successfully exited]\n");
    pthread_exit(NULL);
}

/*
Remove existing log file and write header for the log file.
*/
void writeAlertLogHeader(double threshold_diff_magnitude, double threshold_diff_distance, int dims[2])
{
    remove("logfile.txt");
    FILE *fptr = fopen("logfile.txt", "w");

    fprintf(fptr, "FIT3143 Assignment 2\n\n");
    fprintf(fptr, "MA Lab 6 Team 1\n");
    fprintf(fptr, "Editha Karina Hermawan 32179081\n");
    fprintf(fptr, "Lionie Annabella Wijaya 31316115\n\n");

    fprintf(fptr, "Ocean nodes dimension: %dx%d\n", dims[0], dims[1]);
    fprintf(fptr, "Coordinate different threshold (km): %f\n", threshold_diff_distance);
    fprintf(fptr, "Magnitude different threshold: %f\n", threshold_diff_magnitude);
    fprintf(fptr, "Earthquake magnitude threshold: %f\n", EARTHQUAKE_MAGNITUDE_THRESHOLD);

    fprintf(fptr, "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    fclose(fptr);
}

/*
Write a single log into log file
*/
void writeAlertLog(int *iter, int conclusive, int matching, int totalNeighbor, SensorReading *nodeReading, SensorReading adjacentReadings[], SensorReading *balloonReading, int sourceCoord[2], int neighborCoord[8], float threshold_diff_magnitude, float threshold_diff_distance)
{
    FILE *fptr = fopen("logfile.txt", "a");
    time_t currtime = time(NULL);
    char currtimestr[40];
    timeToString(&currtime, currtimestr);

    fprintf(fptr, "Iteration: %d\n", (*iter));
    fprintf(fptr, "Logged time: %s\n", currtimestr);
    fprintf(fptr, "Alert type: %s\n", conclusive == 1 ? "Conclusive" : "Inconclusive");
    fprintf(fptr, "\n");

    fprintf(fptr, "Reporting Node:\n");
    fprintf(fptr, "Coord: (%d,%d) | Date and time: %s| Seismic coord: (%f,%f)| Magnitude: %f| Depth: %f\n", sourceCoord[0], sourceCoord[1], nodeReading->datetime, nodeReading->longitude, nodeReading->latitude, nodeReading->magnitude, nodeReading->depth);
    fprintf(fptr, "\n");

    fprintf(fptr, "%d out of %d Adjacent Nodes matching to Reporting Node:\n", matching, totalNeighbor);
    for (int i = 0; i < totalNeighbor; i++)
    {
        ReadingComparison *comparison = compareTwoReadings(nodeReading, &adjacentReadings[i], threshold_diff_magnitude, threshold_diff_distance);
        fprintf(fptr, "Coord: (%d,%d) | Date and time: %s| Seismic coord: (%f,%f), Diff(coord,km): %f | Magnitude: %f, Diff(mag) %f| Depth: %f\n", neighborCoord[i * 2], neighborCoord[i * 2 + 1], adjacentReadings[i].datetime, adjacentReadings[i].longitude, adjacentReadings[i].latitude, comparison->diff_dist, adjacentReadings[i].magnitude, comparison->diff_mag, adjacentReadings[i].depth);
    }
    fprintf(fptr, "\n");

    ReadingComparison *comparison = compareTwoReadings(nodeReading, balloonReading, threshold_diff_magnitude, threshold_diff_distance);
    fprintf(fptr, "Balloon seismic reporting time: %s\n", balloonReading->datetime);
    fprintf(fptr, "Balloon seismic reporting coord: (%f,%f)\n", balloonReading->longitude, balloonReading->latitude);
    fprintf(fptr, "Balloon seismic reporting coord diff. with reporting node: %f\n", comparison->diff_dist);
    fprintf(fptr, "Balloon seismic reporting magnitude: %f\n", balloonReading->magnitude);
    fprintf(fptr, "Balloon seismic reporting magnitude diff. with reporting node: %f\n", comparison->diff_mag);
    fprintf(fptr, "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    fclose(fptr);
}