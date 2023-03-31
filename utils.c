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
Compute a random float between range given.
Reference: https://stackoverflow.com/questions/13408990/how-to-generate-random-float-number-in-c 
*/
float random_float(float min, float max)
{
    float scale = rand() / (float)RAND_MAX;
    return min + scale * (max - min);
}

/*
Convert degree to radian.
Reference: https://www.geodatasource.com/developers/c
*/
double deg2rad(double deg)
{
    return (deg * pi / 180);
}

/*
Convert radian to degree.
Reference: https://www.geodatasource.com/developers/c
*/
double rad2deg(double rad)
{
    return (rad * 180 / pi);
}

/*
Compute difference in distance between two points based on longitude and latitude.
Reference: https://www.geodatasource.com/developers/c
*/
double diff_distance(double lat1, double lon1, double lat2, double lon2, char unit)
{
    double theta, dist;
    if ((lat1 == lat2) && (lon1 == lon2))
    {
        return 0;
    }
    else
    {
        theta = lon1 - lon2;
        dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
        dist = acos(dist);
        dist = rad2deg(dist);
        dist = dist * 60 * 1.1515;
        switch (unit)
        {
        case 'M':
            break;
        case 'K':
            dist = dist * 1.609344;
            break;
        case 'N':
            dist = dist * 0.8684;
            break;
        }

        // printf("Diff dist calculation pos 1 %f %f and pos 2 %f %f results with %f\n", lat1, lon1, lat2, lon2, dist);
        return dist;
    }
}

/*
Compute difference in magnitude.
*/
double diff_magnitude(double mag1, double mag2)
{
    // printf("Diff mag calculation %f and %f results with %f\n", mag1, mag2, fabs(mag1 - mag2));
    return fabs(mag1 - mag2);
}

/* Convert time format into */
void timeToString(time_t *currTime, char *currTimeStr)
{
    struct tm *info;
    info = localtime(currTime);
    strftime(currTimeStr, 40, "%Y/%m/%d %H:%M:%S", info);

    // printf("time to string: %s", currTimeStr);
}

/*
Convert string format into time format.
*/
time_t stringToTime(char *currTimeStr)
{
    struct tm newTime = {0};
    struct tm *info;
    char newStr[40];
    int year = 0, month = 0, day = 0, hour = 0, min = 0, sec = 0;
    sscanf(currTimeStr, "%4d/%2d/%2d %2d:%2d:%2d", &year, &month, &day, &hour, &min, &sec);

    newTime.tm_year = year - 1900;
    newTime.tm_mon = month - 1;
    newTime.tm_mday = day;
    newTime.tm_hour = hour;
    newTime.tm_min = min;
    newTime.tm_sec = sec;

    time_t t = mktime(&newTime);
    info = localtime(&t);
    strftime(newStr, 40, "%Y/%m/%d %H:%M:%S", info);

    // printf("string to time: Formatted date & time : |%s|\n", newStr);
    return t;
}

/*
Generate a random reading.
*/
void generateReading(SensorReading *reading, double minmag, double lat1, double lon1, double lat2, double lon2)
{
    reading->magnitude = random_float(minmag, MAX_MAG);
    reading->depth = random_float(MIN_DEPTH, MAX_DEPTH);
    reading->latitude = random_float(lat1, lat2);
    reading->longitude = random_float(lon1, lon2);

    time_t currtime = time(NULL);
    timeToString(&currtime, reading->datetime);

    // printf("generated [%s] mag %f dep %f lat %f long %f\n",reading->datetime,reading->magnitude,reading->depth,reading->latitude,reading->longitude);
}

/*
Generate a random reading for ocean node.
*/
void generateNodeReading(SensorReading *reading, int dims[], int coords[])
{
    float node_width = (END_LONG - START_LONG) / dims[1];
    float node_height = (END_LAT - START_LAT) / dims[0];
    float min_lat = START_LAT + (node_height * coords[1]);
    float min_lon = START_LONG + (node_width * coords[0]);
    float overlap_coord_lon = 0.2 * node_width;
    float overlap_coord_lat = 0.2 * node_height;
    float max_lat = coords[1] == dims[0] - 1 ? min_lat + node_height : min_lat + node_height + overlap_coord_lat;
    float max_lon = coords[0] == dims[1] - 1 ? min_lon + node_width : min_lon + node_width + overlap_coord_lon;

    return generateReading(reading, MIN_MAG, min_lat, min_lon, max_lat, max_lon);
}

/*
Compare two readings.
*/
ReadingComparison *compareTwoReadings(SensorReading *readingA, SensorReading *readingB, float threshold_diff_magnitude, float threshold_diff_distance)
{
    ReadingComparison *comparison = (ReadingComparison *)malloc(sizeof(ReadingComparison *));
    comparison->diff_mag = diff_magnitude(readingA->magnitude, readingB->magnitude);
    comparison->diff_dist = diff_distance(readingA->latitude, readingA->longitude, readingB->latitude, readingB->longitude, 'K');
    comparison->longitude = readingB->longitude;
    comparison->latitude = readingB->latitude;
    comparison->matching = comparison->diff_mag <= threshold_diff_magnitude && comparison->diff_dist <= threshold_diff_distance ? 1 : 0;
    // printf("compare two readings:\n latA: %f, longA: %f, magA: %f\n latB:%f, longB:%f, magB:%f\n\n", readingA->latitude, readingA->longitude, readingA->magnitude, readingB->latitude, readingB->longitude, readingB->magnitude);

    return comparison;
}