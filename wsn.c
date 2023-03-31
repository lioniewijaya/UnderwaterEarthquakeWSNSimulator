/*
    FIT3143 Parallel Computing
    Assignment 2

    Student: Editha Karina Hermawan 32179081
    Student: Lionie Annabella Wijaya 31316115

    To compile:
      make all
      or
      mpicc wsn.c utils.c alert.c balloon.c node.c -lm -lpthread -o wsn
    To run via command line argument in VM:
      version 1 - fixed iterations
      mpirun --oversubscribe -np [number of process] wsn [grid row] [grid column] [difference in magnitude threshold] [difference in distance threshold] [number of iterations]
      version 2 - programs stops during runtime after user set a sentinel value of 2 in sentinel.txt
      mpirun --oversubscribe -np [number of process] wsn [grid row] [grid column] [difference in magnitude threshold] [difference in distance threshold] 
    To run via CAAS:
      srun wsn [grid row] [grid column] [difference in magnitude threshold] [difference in distance threshold] [number of iterations]
*/

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

int main(int argc, char **argv)
{
  int my_rank, color, newRank, tmpRank;
  int p, i, iterations, provided = 0, sentinel = 0, currentIteration = 1;
  time_t start, end;
  double time_taken;
  pthread_t nodeGenerateThr, nodeRecvThr, balloonThr, alertThr;
  int sizeFloat, sizeTime, sizeCoord, totalSize;
  int dims[2] = {-1, -1};
  float threshold_diff_magnitude, threshold_diff_distance;

  // MPI Setup
  MPI_Comm cartcomm = NULL;
  MPI_Comm tmpcartcomm = NULL;
  MPI_Group all_grp, node_grp;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_group(MPI_COMM_WORLD, &all_grp);

  srand(time(NULL) + my_rank);

  // Get user input
  if (my_rank == 0)
  {
    dims[0] = atoi(argv[1]);
    dims[1] = atoi(argv[2]);
    threshold_diff_magnitude = atof(argv[3]);
    threshold_diff_distance = atof(argv[4]);
    printf("Dimension %dx%d | Difference in magnitude threshold %f | Difference in distance threshold %f | ", dims[0], dims[1], threshold_diff_magnitude, threshold_diff_distance);

    if (argc == 6)
    {
      iterations = atoi(argv[5]);
      printf("Iterations %d\n", iterations);
    }
    else
    {
      printf("Iterations is based on sentinel value set during runtime, termination only after user save a sentinel value of 2 in sentinel.txt file\n");
    }
  }

  // Share user input values to all ranks
  MPI_Bcast(&dims[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dims[1], 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&threshold_diff_magnitude, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&threshold_diff_distance, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Ensure correct number of processes is used
  if (dims[0] * dims[1] != p - 1)
  {
    if (my_rank == 0)
    {
      printf("Ensure (mxn)+1 is equal to number of processes!\n");
    }
    MPI_Finalize();
    return 0;
  }

  if (my_rank == p - 1)
  {
    color = 1;
  }
  else
  {
    color = 2;
  }
  MPI_Comm_split(MPI_COMM_WORLD, color, my_rank, &tmpcartcomm);

  int coords[2], nbrs[4], periods[2] = {0, 0};

  if (my_rank != p - 1)
  {
    // Create new communicator for cartesian grid
    MPI_Cart_create(tmpcartcomm, 2, dims, periods, 0, &cartcomm);
    // Assign n-1 processes into cartesian to be assigned as ocean nodes
    MPI_Comm_rank(cartcomm, &newRank);
    MPI_Cart_coords(cartcomm, newRank, 2, coords);
    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
    i = 0;

    // printf("rank= %d coords= %d %d  neighbors(u,d,l,r)= %d %d %d %d\n", newRank, coords[0], coords[1], nbrs[UP], nbrs[DOWN], nbrs[LEFT], nbrs[RIGHT]);
    // printf("new rank: %d, my_rank:%d\n", newRank, my_rank);

    // Find total buffer size for MPI_Pack
    MPI_Pack_size(4, MPI_FLOAT, cartcomm, &sizeFloat);
    MPI_Pack_size(40, MPI_CHAR, cartcomm, &sizeTime);
    MPI_Pack_size(2, MPI_INT, cartcomm, &sizeCoord);
    totalSize = sizeCoord + sizeFloat + sizeTime;

    char *packBuffer = (char *)malloc(totalSize * sizeof(char));

    int flag;
    MPI_Request request[4];
    MPI_Request recvReq;
  }
  // Ensure balloon and base station are created first before ocean nodes
  else
  {
    // Write log header
    writeAlertLogHeader(threshold_diff_distance, threshold_diff_magnitude, dims);

    /* 2.0) Simulating the balloon seismic sensor */
    BalloonThread *argsBalloon = (BalloonThread *)malloc(sizeof(BalloonThread));
    argsBalloon->sentinel = &sentinel;
    pthread_create(&balloonThr, NULL, (void *)balloonThread, (void *)argsBalloon);

    /* 3.0) Simulating the base station */
    AlertThread *argsAlert = (AlertThread *)malloc(sizeof(AlertThread));
    argsAlert->comm = MPI_COMM_WORLD;
    argsAlert->iteration = &currentIteration;
    argsAlert->threshold_diff_magnitude = threshold_diff_magnitude;
    argsAlert->threshold_diff_distance = threshold_diff_distance;
    argsAlert->sentinel = &sentinel;
    pthread_create(&alertThr, NULL, (void *)alertThread, (void *)argsAlert);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /* 1.0) Simulating the seafloor seismic sensor node */
  if (my_rank != p - 1)
  {
    // Each ocean node stores one latest reading, new generated reading will overwrite previous reading - to start off, a random reading is generated here to avoid null comparison
    SensorReading *reading = (SensorReading *)malloc(sizeof(SensorReading *) + sizeof(char) * 40);
    generateNodeReading(reading, dims, coords);
    MPI_Barrier(cartcomm);

    // Generate thread to generate and receive info from adjacent nodes
    GenerateThread *argsGen = (GenerateThread *)malloc(sizeof(GenerateThread) + sizeof(int) * 8);
    argsGen->cartComm = cartcomm;
    argsGen->mainComm = MPI_COMM_WORLD;
    argsGen->reading = reading;
    memcpy(argsGen->nbrs, nbrs, 4 * sizeof(int));
    memcpy(argsGen->dims, dims, 2 * sizeof(int));
    memcpy(argsGen->coords, coords, 2 * sizeof(int));
    argsGen->totalSize = totalSize;
    argsGen->threshold_diff_magnitude = threshold_diff_magnitude;
    argsGen->threshold_diff_distance = threshold_diff_distance;
    argsGen->newRank = newRank;
    argsGen->baseRank = p - 1;
    argsGen->sentinel = &sentinel;
    pthread_create(&nodeGenerateThr, NULL, (void *)nodeGenerateThread, (void *)argsGen);

    // Receive thread to receive signal indicating request to share info)
    RecvThread *argsRecv = (RecvThread *)malloc(sizeof(RecvThread) + sizeof(int) * 2);
    argsRecv->cartComm = cartcomm;
    argsRecv->reading = reading;
    memcpy(argsRecv->coords, coords, 2 * sizeof(int));
    argsRecv->totalSize = totalSize;
    argsRecv->rank = my_rank;
    argsRecv->sentinel = &sentinel;
    pthread_create(&nodeRecvThr, NULL, (void *)nodeRecvThread, (void *)argsRecv);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Ensure balloon, base station, and ocean nodes are created first before simulation starts
  if (my_rank != p - 1)
  {
    // Received sentinel from base station to terminate threads
    // printf("****** Waiting to receive sentinel (current: %d) from base to this rank %d\n", sentinel, my_rank);
    MPI_Recv(&sentinel, 1, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, NULL);
    // printf("****** Received sentinel (current: %d) from base to this rank %d\n", sentinel, my_rank);

    // Clean up threads
    pthread_join(nodeGenerateThr, NULL);
    pthread_join(nodeRecvThr, NULL);
  }
  else
  {
    // Start simulation
    printf("--------- Starting Earthquake Detection Simulation\n");

    if (argc == 6)
    {
      // Simulate n iterations
      while (currentIteration <= iterations)
      {
        printf("********* Currently in iteration %d ... \n", currentIteration);
        sleep(CYCLE_NODE);
        currentIteration += 1;
      }
      sentinel = SENTINEL;
    }
    else
    {
      while (1)
      {
        FILE *sentinelTxt = fopen("sentinel.txt", "r");
        fscanf(sentinelTxt, "%d", &sentinel);
        fclose(sentinelTxt);

        if (sentinel == SENTINEL)
        {
          printf("User has set a sentinel value to terminate program in sentinel.txt\n");
          break;
        }

        printf("********* Currently in iteration %d ... sentinel.txt has no sentinel value of %d set yet\n", currentIteration, SENTINEL);
        sleep(CYCLE_NODE);
        currentIteration += 1;
      }
    }
    currentIteration -= 1;

    // Simulation is completed, set sentinel and send to ocean nodes to terminate threads
    for (int i = 0; i < my_rank; i++)
    {
      MPI_Send(&sentinel, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
      // printf("****** Sent sentinel (current: %d) from base %d to other ranks\n", my_rank, sentinel);
    }

    // Clean up threads
    pthread_join(alertThr, NULL);
    pthread_join(balloonThr, NULL);

    printf("--------- Earthquake Detection Simulation Completed\n");
  }

  MPI_Finalize();
  return 0;
}