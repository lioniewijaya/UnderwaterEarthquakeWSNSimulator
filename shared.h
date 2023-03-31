/* Structures */

typedef struct
{
  char datetime[40];
  float magnitude;
  float depth;
  float longitude;
  float latitude;
} SensorReading;

typedef struct
{
  int capacity;
  int size;
  int front;
  int rear;
  SensorReading *readings[10];
} BalloonReadings;

typedef struct
{
  MPI_Comm comm;
  int *iteration;
  int *sentinel;
  double threshold_diff_magnitude;
  double threshold_diff_distance;
} AlertThread;

typedef struct
{
  int *sentinel;
} BalloonThread;

typedef struct
{
  MPI_Comm cartComm;
  MPI_Comm mainComm;
  int nbrs[4];
  int dims[2];
  int coords[2];
  int newRank;
  int baseRank;
  int totalSize;
  float threshold_diff_magnitude;
  float threshold_diff_distance;
  int *sentinel;
  SensorReading *reading;
} GenerateThread;

typedef struct
{
  MPI_Comm cartComm;
  SensorReading *reading;
  int coords[2];
  int totalSize;
  int *sentinel;
  int rank;
} RecvThread;

typedef struct
{
  float diff_mag;
  float diff_dist;
  float longitude;
  float latitude;
  int matching;
} ReadingComparison;

/* Constants */

#define pi 3.14159265358979323846
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3
#define START_LAT 20.0
#define START_LONG 20.0
#define END_LAT 40.0
#define END_LONG 40.0
#define MIN_DEPTH 20
#define MAX_DEPTH 170
#define MIN_MAG 3.0
#define MAX_MAG 7.0
#define CYCLE_NODE 2.0
#define EARTHQUAKE_MAGNITUDE_THRESHOLD 2.5
// #define DIFF_MAG_THRESHOLD 1.0
// #define DIFF_DISTANCE_THRESHOLD 1500.0
#define SENTINEL 2