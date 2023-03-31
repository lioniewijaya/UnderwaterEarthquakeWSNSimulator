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

BalloonReadings balloonReadings =
    {.capacity = 10,
     .size = 0,
     .front = 0,
     .rear = 0,
     .readings = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

/*
Balloon thread to generate reading into static array per time interval given.
*/
void *balloonThread(void *args)
{
    // printf("[Balloon Thread is successfully created]\n");

    BalloonThread *data = (BalloonThread *)args;

    pthread_mutex_t lock;
    struct timespec start, end;
    double time_taken;

    // Continue only if base station has not send termination signal
    while (*(data->sentinel) != SENTINEL)
    {
        pthread_mutex_lock(&lock);
        clock_gettime(CLOCK_MONOTONIC, &start);

        // Balloon array is full, overwrite oldest reading
        if (balloonReadings.size == balloonReadings.capacity)
        {
            free(balloonReadings.readings[balloonReadings.front]);
            balloonReadings.readings[balloonReadings.front] = NULL;
            balloonReadings.front = (balloonReadings.front + 1) % balloonReadings.capacity;
            balloonReadings.size -= 1;
        }

        // Generate a new reading to add to balloon array
        SensorReading *reading = (SensorReading *)malloc(sizeof(SensorReading *) + sizeof(char) * 40);
        generateReading(reading, EARTHQUAKE_MAGNITUDE_THRESHOLD, START_LAT, START_LONG, END_LAT, END_LONG);
        balloonReadings.readings[balloonReadings.rear] = reading;
        balloonReadings.rear = (balloonReadings.rear + 1) % balloonReadings.capacity;
        balloonReadings.size += 1;

        pthread_mutex_unlock(&lock);
        clock_gettime(CLOCK_MONOTONIC, &end);

        time_taken = (end.tv_sec - start.tv_sec) * 1e9;
        time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;
        if (time_taken < CYCLE_NODE)
        {
            usleep((CYCLE_NODE - time_taken) * 1e6);
        }

        // printf("Balloon Thread added a reading, current properties - rear: %d, front: %d, size:%d\n", balloonReadings.rear, balloonReadings.front, balloonReadings.size);
        //  for (int i = 0; i < balloonReadings.size; i++)
        //  {
        //    printf("Balloon array element %d: [%s] mag %f dep %f lat %f long %f\n", i, balloonReadings.readings[i]->datetime, balloonReadings.readings[i]->magnitude, balloonReadings.readings[i]->depth, balloonReadings.readings[i]->latitude, balloonReadings.readings[i]->longitude);
        //  }
    }
    // printf("[Balloon Thread successfully exited]\n");
    pthread_exit(NULL);
}

/*
Compare a node reading to the three latest balloon reading before writing into log file.
*/
void compareBalloonToReportingNodeReadings(int *iter, int matching, int totalNeighbor, SensorReading *nodeReading, SensorReading adjacentReadings[], int sourceCoord[2], int neighborCoord[8], double threshold_diff_magnitude, double threshold_diff_distance)
{
    SensorReading *balloonReading;
    ReadingComparison *comparison;
    int conclusive = 0;

    int front = balloonReadings.front;
    int rear = balloonReadings.rear;

    int loops = balloonReadings.size >= 3 ? 3 : balloonReadings.size;
    int i = 0, curr;
    if (balloonReadings.size == 10)
    {
        curr = rear - 1 >= 0 ? rear - 1 : 9;
    }
    else
    {
        curr = rear - 1 >= 0 ? rear - 1 : 0;
    }
    // printf("Balloon comparing to source node - front %d rear %d size %d, calculated curr %d\n", balloonReadings.front, balloonReadings.rear, balloonReadings.size, curr);

    // Compare the latest 3 readings in balloon array only
    while (i < loops)
    {
        if (balloonReadings.readings[curr] != NULL)
        {
            comparison = compareTwoReadings(balloonReadings.readings[curr], nodeReading, threshold_diff_magnitude, threshold_diff_distance);
            balloonReading = balloonReadings.readings[curr];
            if (comparison->matching)
            {
                conclusive = 1;
                break;
            }
        }

        // printf("Balloon comparing to source node - capacity %d, old rear: %d, new rear: %d", balloonReadings.size, curr, (curr - 1) % balloonReadings.capacity);

        if (balloonReadings.size == 10)
        {
            curr = curr - 1 >= 0 ? curr - 1 : 9;
        }
        else
        {
            curr = curr - 1 >= 0 ? curr - 1 : 0;
        }

        i += 1;
    }

    // Write to alert log
    if (balloonReading != NULL)
    {
        writeAlertLog(iter, conclusive, matching, totalNeighbor, nodeReading, adjacentReadings, balloonReading, sourceCoord, neighborCoord, threshold_diff_magnitude, threshold_diff_distance);
    }
}