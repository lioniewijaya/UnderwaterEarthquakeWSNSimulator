#ifndef UTILS_H
#define UTILS_H

double deg2rad(double);
double rad2deg(double);
float random_float(float min, float max);
double diff_distance(double lat1, double lon1, double lat2, double lon2, char unit);
double diff_magnitude(double mag1, double mag2);
void timeToString(time_t *currTime, char *currTimeStr);
time_t stringToTime(char *theString);

void generateReading(SensorReading *reading, double minmag, double lat1, double lon1, double lat2, double lon2);
void generateNodeReading(SensorReading *reading, int dims[], int coords[]);
ReadingComparison *compareTwoReadings(SensorReading *readingA, SensorReading *readingB, float threshold_diff_magnitude, float threshold_diff_distance);

#endif