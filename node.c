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
Node generate thread to generate reading for node per time interval given, ask information from adjacent nodes if magnitude exceeds
threshold, and alerts the base station if 2 or more adjacent nodes matches current reading.
*/
void *nodeGenerateThread(void *args)
{
    // printf("[Generate Thread successfully created]\n");

    GenerateThread *data = (GenerateThread *)args;

    MPI_Comm mainComm = data->mainComm;
    MPI_Comm currCartcomm = data->cartComm;
    int nbrs[4], dims[2], coords[2];
    memcpy(nbrs, data->nbrs, 4 * sizeof(int));
    memcpy(dims, data->dims, 2 * sizeof(int));
    memcpy(coords, data->coords, 2 * sizeof(int));
    int size = data->totalSize;
    float threshold_diff_magnitude = data->threshold_diff_magnitude;
    float threshold_diff_distance = data->threshold_diff_distance;
    int currRank = data->newRank;
    int baseRank = data->baseRank;
    char packBufferRecv[size];

    pthread_mutex_t lock;
    struct timespec start, end;
    double time_taken;

    // Continue only if base station has not send termination signal
    while (*(data->sentinel) != SENTINEL)
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
        // printf("G - Node generate sentinel (current: %d)\n", *(data->sentinel));

        pthread_mutex_lock(&lock);
        generateNodeReading(data->reading, dims, coords);
        pthread_mutex_unlock(&lock);

        // printf("G - Sensor x: %d, y: %d read = {lat: %f, long: %f, depth:  %f, mag: %f}\n", coords[0], coords[1], data->reading->latitude, data->reading->longitude, data->reading->depth, data->reading->magnitude);

        // Send request to adjacent nodes if reading exceed threshold
        if (data->reading->magnitude > EARTHQUAKE_MAGNITUDE_THRESHOLD)
        {
            // printf("G - Request from adjacent nodes\n");
            int pos = 0;
            MPI_Request request[4];

            // Send source rank to adjacent nodes
            for (int k = 0; k < 4; k++)
            {
                if (nbrs[k] >= 0)
                {
                    MPI_Send(&currRank, 1, MPI_INT, nbrs[k], 0, currCartcomm);
                    // printf("G - Current rank %d sent signal to rank %d\n", currRank, nbrs[k]);
                }
            }

            MPI_Status status[4];
            SensorReading neighborReading[4];
            ReadingComparison *compareResult[4];
            char *adjacentNodeBuffer[4];
            int matchings = 0;
            int totalNeighbour = 0;

            // Receive information from adjacent nodes
            for (int k = 0; k < 4; k++)
            {
                if (nbrs[k] >= 0)
                {
                    pos = 0;
                    // printf("G - Node source rank %d waiting for receive\n", currRank);

                    MPI_Status probe_status;
                    int flag = 0;
                    MPI_Iprobe(MPI_ANY_SOURCE, 1, currCartcomm, &flag, &probe_status);

                    // Continue only if generate node actually receive message from receive node
                    if (flag == 1)
                    {
                        MPI_Recv(packBufferRecv, size, MPI_PACKED, nbrs[k], 1, currCartcomm, &status[k]);

                        // Unpack message received - format [magnitude] [latitude] [longitude] [depth] [datetime] [coords]
                        MPI_Unpack(packBufferRecv, size, &pos, &neighborReading[k].magnitude, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(packBufferRecv, size, &pos, &neighborReading[k].latitude, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(packBufferRecv, size, &pos, &neighborReading[k].longitude, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(packBufferRecv, size, &pos, &neighborReading[k].depth, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(packBufferRecv, size, &pos, &neighborReading[k].datetime, 40, MPI_CHAR, currCartcomm);

                        // Compare and calculate matching adjacent nodes
                        compareResult[k] = compareTwoReadings(data->reading, &neighborReading[k], threshold_diff_magnitude, threshold_diff_distance);
                        if (compareResult[k]->matching == 1)
                        {
                            matchings += 1;
                        }

                        char *bufferSingle = (char *)malloc(sizeof(char) * size);
                        memcpy(bufferSingle, packBufferRecv, size);
                        adjacentNodeBuffer[totalNeighbour] = bufferSingle;

                        /*
                        pos = 0;
                        MPI_Unpack(adjacentNodeBuffer[totalNeighbour], size, &pos, &neighborReading[k].magnitude, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(adjacentNodeBuffer[totalNeighbour], size, &pos, &neighborReading[k].latitude, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(adjacentNodeBuffer[totalNeighbour], size, &pos, &neighborReading[k].longitude, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(adjacentNodeBuffer[totalNeighbour], size, &pos, &neighborReading[k].depth, 1, MPI_FLOAT, currCartcomm);
                        MPI_Unpack(adjacentNodeBuffer[totalNeighbour], size, &pos, &neighborReading[k].datetime, 40, MPI_CHAR, currCartcomm);
                        */

                        // printf("G - Receive at source, Current rank %d received from rank %d, match: %d, diff mag: %f, diff dist: %f\n", currRank, nbrs[k], compareResult[k]->matching, compareResult[k]->diff_mag, compareResult[k]->diff_dist);

                        totalNeighbour += 1;
                    }
                }
            }

            // Send alert to base station
            if (matchings >= 2)
            {

                // Pack message to send - format [totalSize] [totalNeighbor] [sourceNode] [adjacentNodes]
                int totalSize = sizeof(int) * 2 + size * (totalNeighbour + 1);
                char packBufferSent[totalSize];

                pos = 0;

                // Total neighbor, matchings
                MPI_Pack(&totalNeighbour, 1, MPI_INT, packBufferSent, totalSize, &pos, mainComm);
                MPI_Pack(&matchings, 1, MPI_INT, packBufferSent, totalSize, &pos, mainComm);

                // Source node
                MPI_Pack(&data->reading->magnitude, 1, MPI_FLOAT, packBufferSent, size, &pos, mainComm);
                MPI_Pack(&data->reading->latitude, 1, MPI_FLOAT, packBufferSent, totalSize, &pos, mainComm);
                MPI_Pack(&data->reading->longitude, 1, MPI_FLOAT, packBufferSent, totalSize, &pos, mainComm);
                MPI_Pack(&data->reading->depth, 1, MPI_FLOAT, packBufferSent, totalSize, &pos, mainComm);
                MPI_Pack(data->reading->datetime, 40, MPI_CHAR, packBufferSent, totalSize, &pos, mainComm);
                MPI_Pack(coords, 2, MPI_INT, packBufferSent, totalSize, &pos, mainComm);

                // Adjacent nodes
                for (int i = 0; i < totalNeighbour; i++)
                {
                    MPI_Pack(adjacentNodeBuffer[i], size, MPI_CHAR, packBufferSent, totalSize, &pos, mainComm);
                }

                MPI_Send(packBufferSent, totalSize, MPI_PACKED, baseRank, 1, mainComm);

                for (int i = 0; i < totalNeighbour; i++)
                {
                    free(adjacentNodeBuffer[i]);
                }
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &end);

        time_taken = (end.tv_sec - start.tv_sec) * 1e9;
        time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
        if (time_taken < CYCLE_NODE)
        {
            usleep((CYCLE_NODE - time_taken) * 1e6);
            // printf("node finished sleeping, start %f end %f timetaken: %lf rank: %d in cartcomm. Time left: %f\n", start, end, time_taken, currRank, (CYCLE_NODE - time_taken));
        }
    }
    // printf("[Generate Thread successfully exited]\n");
    pthread_exit(NULL);
}

/*
Node receive thread to receive request from an adjacent source node and return information asked.
*/
void *nodeRecvThread(void *args)
{
    // printf("[Receive Thread successfully created]\n");

    RecvThread *data = (RecvThread *)args;
    MPI_Comm currCartcomm = data->cartComm;
    SensorReading *reading = data->reading;
    int coords[2];
    memcpy(coords, data->coords, 2 * sizeof(int));
    int totalSize = data->totalSize;
    int sourceRank;
    int my_rank = data->rank;
    char packBuffer[totalSize];

    MPI_Status status;
    MPI_Request request;

    // Continue only if base station has not send termination signal
    while (*(data->sentinel) != SENTINEL)
    {
        MPI_Status probe_status;
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, 0, currCartcomm, &flag, &probe_status);

        // Continue only if receive node actually receive message from generate node
        if (flag == 1)
        {
            // printf("R - Adjacent node received sentinel (current: %d)\n", *(data->sentinel));
            MPI_Recv(&sourceRank, 1, MPI_INT, MPI_ANY_SOURCE, 0, currCartcomm, &status);
            // printf("R - Adjacent node with rank %d received signal from source rank %d\n", my_rank, sourceRank);

            int pos = 0;

            // Pack message to send - format [magnitude] [latitude] [longitude] [depth] [datetime] [coords]
            MPI_Pack(&reading->magnitude, 1, MPI_FLOAT, packBuffer, totalSize, &pos, currCartcomm);
            MPI_Pack(&reading->latitude, 1, MPI_FLOAT, packBuffer, totalSize, &pos, currCartcomm);
            MPI_Pack(&reading->longitude, 1, MPI_FLOAT, packBuffer, totalSize, &pos, currCartcomm);
            MPI_Pack(&reading->depth, 1, MPI_FLOAT, packBuffer, totalSize, &pos, currCartcomm);
            MPI_Pack(reading->datetime, 40, MPI_CHAR, packBuffer, totalSize, &pos, currCartcomm);
            MPI_Pack(coords, 2, MPI_INT, packBuffer, totalSize, &pos, currCartcomm);

            MPI_Send(packBuffer, totalSize, MPI_PACKED, sourceRank, 1, currCartcomm);
        }
    }
    // printf("[Receive Thread successfully exited]\n");
    pthread_exit(NULL);
}