#ifndef BALLOON_H
#define BALLOON_H

void *balloonThread(void *args);
void compareBalloonToReportingNodeReadings(int *iter, int matching, int totalNeighbor, SensorReading *nodeReading, SensorReading adjacentReadings[], int sourceCoord[2], int neighborCoord[8], double threshold_diff_magnitude, double threshold_diff_distance);

extern BalloonReadings balloonReadings;

#endif