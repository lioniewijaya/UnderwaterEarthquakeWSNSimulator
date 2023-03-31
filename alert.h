#ifndef ALERT_H
#define ALERT_H

void *alertThread(void *args);
void writeAlertLogHeader(double threshold_diff_magnitude, double threshold_diff_distance, int dims[2]);
void writeAlertLog(int *iter, int conclusive, int matching, int totalNeighbor, SensorReading *nodeReading, SensorReading adjacentReadings[], SensorReading *balloonReading, int sourceCoord[2], int neighborCoord[8], float threshold_diff_magnitude, float threshold_diff_distance);

#endif