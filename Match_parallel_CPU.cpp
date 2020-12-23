/*This program:
* -Reads data from the LiDAR OS1-16
* -Transforms the data from polar to Cartesian coordinates
* -Rotates and translates the point cloud
* -Computes the correspondence between the 2 clouds by finding the closest point (NNS) -> Match
* PARALLEL PROGRAM USING PTHREADS LIBRARY
By: Carlos Huapaya*/
#define HAVE_STRUCT_TIMESPEC

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "mkl.h"//only for measuring time

#define NUM_POINTS 16384//number of points in the cloud
#define NUMTH 8//number of procesors to use

//struct to handle the point cloud
typedef struct
{
	float* P;
	float* Q;
	int* idx;//correspondence
	int len;
}POINTCLOUD;

POINTCLOUD cloudstr;//define the struct
pthread_t thread[NUMTH];//threads
pthread_attr_t tattr;//Attribute variable

void* Match(void* arg)
{
	int offset = (int)(arg);
	
	//variables for each thread
	int len = cloudstr.len;
	int start = len * offset;
	int end = start + len;
	float* P = cloudstr.P;
	float* Q = cloudstr.Q;
	int* idx = cloudstr.idx;
	float dist;
	float max_dist = 1000000000.0f;//maximum distance

	for (int i = start; i < end; i++)
	{
		/*each thread is in charge of a group of p points
		so MUTEX variables are not used here*/
		max_dist = 1000000000.0f;
		for (int j = 0; j < NUM_POINTS; j++)
		{
			dist = (P[0 + i * 3] - Q[0 + j * 3]) * (P[0 + i * 3] - Q[0 + j * 3]) +
				(P[1 + i * 3] - Q[1 + j * 3]) * (P[1 + i * 3] - Q[1 + j * 3]) +
				(P[2 + i * 3] - Q[2 + j * 3]) * (P[2 + i * 3] - Q[2 + j * 3]);
			if (dist < max_dist)
			{
				max_dist = dist;
				idx[i] = j;
			}
		}
	}

	pthread_exit(NULL);
	return NULL;
}

int main(void)
{
	/////////////////////Block 1: Open and read the file with 64 LiDAR data packets////////////////
	int i = 0;
	const int N_LINE = 128;
	char line[N_LINE];

	FILE* document;
	fopen_s(&document, "Donut_1024x16.csv", "r");//open the file
	if (!document) {
		perror("File opening failed");
		return (-1);
	}
	int j = 1;//number of line
	int channel = 2;
	int azimuth_block = 0;
	int lidar_packet = 0;
	int word;
	int offset = 0;
	float r[NUM_POINTS] = {};//ranges
	unsigned long int encoder_count = 0;//to read the initial encoder counter value (then it grows by 88 ticks)
	unsigned long int aux = 0;
	while (fgets(line, N_LINE, document) != NULL)
	{
		//get the intial encoder counter value
		if (j == 13) encoder_count = atoi(line);
		if (j == 14) encoder_count = atoi(line) << 8 | encoder_count;

		//read ranges
		word = 17 + 12 * channel + 788 * azimuth_block + 12608 * lidar_packet;
		if (j == word) aux = (unsigned long int) atoi(line);
		if (j == word + 1) aux = (unsigned long int) (atoi(line) << 8) | aux;
		if (j == word + 2) aux = (unsigned long int) ((atoi(line) & 0x0000000F) << 16) | aux;

		if (j > (word + 2))
		{
			r[offset] = (float)aux;
			offset++;
			channel += 4;
		}
		if (channel >= 64)
		{
			channel = 2;
			azimuth_block++;
		}
		if (azimuth_block >= 16)
		{
			azimuth_block = 0;
			lidar_packet++;
		}
		if (lidar_packet >= 64) break;
		j++;
	}
	//printf("%d\n", encoder_count);
	//for (i = 0; i < 5000; i++) printf("%.3f\n", r[i]);
	fclose(document);

	//read the beam_intrinsics file
	fopen_s(&document, "beam_intrinsics.csv", "r");//open the file
	if (!document) {
		perror("File opening failed");
		return (-1);
	}

	float altitude[16] = {};
	float azimuth[16] = {};
	j = 1;
	while (fgets(line, N_LINE, document) != NULL)
	{
		//read altitute angles
		if (j == 2) offset = 0;
		if (j >= 2 && j <= 65)
		{
			if (j % 4 == 0)
			{
				altitude[offset] = (float)atof(line);
				offset++;
			}
		}

		//read azimuth angles
		if (j == 68) offset = 0;
		if (j >= 68 && j <= 131)
		{
			if ((j - 66) % 4 == 0)
			{
				azimuth[offset] = (float)atof(line);
				offset++;
			}
		}
		j++;
	}
	fclose(document);
	/*printf("Altitute angles:\n");
	for (i = 0; i < 16; i++) printf("%.3f\n", altitude[i]);
	printf("\nAzimuth angles:\n");
	for (i = 0; i < 16; i++) printf("%.3f\n", azimuth[i]);*/

	///////////////////////////////////////////End of Block 1/////////////////////////////////

	///////////////////////Block 2: Conversion to Cartesian coordinates///////////////////////

	float* P = (float*)malloc((size_t)3 * (size_t)NUM_POINTS * sizeof(float));
	if (P == NULL)
	{
		printf("Error: Unable to allocate memory for the point cloud\n");
		return (-1);
	}
	float theta = 0.0;//azimuth angle
	float phi = 0.0;//altitude angle
	//data format: x1y1z1 x2y2z2 x3y3z3 ....
	unsigned long int counter;
	for (i = 0; i < NUM_POINTS; i++)
	{
		azimuth_block = i / 16;
		counter = (encoder_count + azimuth_block * 88) % 90112;
		channel = i % 16;
		theta = (float)(2 * M_PI * (counter / 90112.0 + azimuth[channel] / 360.0));
		phi = (float)(2 * M_PI * altitude[channel] / 360.0);
		P[0 + 3 * i] = (float)(r[i] * cos(theta) * cos(phi));//x
		P[1 + 3 * i] = (float)(-r[i] * sin(theta) * cos(phi));//y
		P[2 + 3 * i] = (float)(r[i] * sin(phi));//z
	}
	//////////////////////////////End of Block 2//////////////////////////////

	//////////////////////////////Block 3: Compute the rotation and translation//////////////////////////////

	//converted cloud
	float* Q = (float*)malloc((size_t)3 * (size_t)NUM_POINTS * sizeof(float));

	//rotation matrix and translation vector
	float* T = (float*)malloc(3 * sizeof(float));
	float* R = (float*)malloc(9 * sizeof(float));
	if (T == NULL || R == NULL)
	{
		printf("Error: Unable to allocate memory for T or R\n");
		return (-1);
	}

	//Translation values
	T[0] = 0.8f;//x
	T[1] = -0.3f;//y
	T[2] = 0.2f;//z

	//Rotation values (rad)
	float rx = 0.2f;//axis x
	float ry = -0.2f;//axis y
	float rz = 0.05f;//axis z

	float cx = (float)cos(rx); float cy = (float)cos(ry); float cz = (float)cos(rz);
	float sx = (float)sin(rx); float sy = (float)sin(ry); float sz = (float)sin(rz);
	R[0] = cy * cz; R[1] = (cz * sx * sy) + (cx * sz); R[2] = -(cx * cz * sy) + (sx * sz);
	R[3] = -cy * sz; R[4] = (cx * cz) - (sx * sy * sz); R[5] = (cx * sy * sz) + (cz * sx);
	R[6] = sy; R[7] = -cy * sx; R[8] = cx * cy;

	for (int i = 0; i < NUM_POINTS; i++)
	{
		Q[0 + i * 3] = R[0 + 0 * 3] * P[0 + i * 3] +
			R[0 + 1 * 3] * P[1 + i * 3] +
			R[0 + 2 * 3] * P[2 + i * 3] + T[0];
		Q[1 + i * 3] = R[1 + 0 * 3] * P[0 + i * 3] +
			R[1 + 1 * 3] * P[1 + i * 3] +
			R[1 + 2 * 3] * P[2 + i * 3] + T[1];
		Q[2 + i * 3] = R[2 + 0 * 3] * P[0 + i * 3] +
			R[2 + 1 * 3] * P[1 + i * 3] +
			R[2 + 2 * 3] * P[2 + i * 3] + T[2];
	}
	//////////////////////////////////////End of Block 3//////////////////////////////////////

	//////////////////////////////Block 4: Compute the match//////////////////////////////////////

	//correspondence for P points
	int* idx = (int*)malloc(NUM_POINTS * sizeof(int));

	//Initialize the struct
	cloudstr.P = P;
	cloudstr.Q = Q;
	cloudstr.len = NUM_POINTS / NUMTH;
	cloudstr.idx = idx;

	//Initialize attribute as joinable
	pthread_attr_init(&tattr);
	pthread_attr_setdetachstate(&tattr, PTHREAD_CREATE_JOINABLE);

	int rc;
	void* status;
	float start, end;

	start = dsecnd();
	//Create threads
	for (int i = 0; i < NUMTH; i++)
	{
		rc = pthread_create(&thread[i], &tattr, Match, (void*)i);
		if (rc)
		{
			printf("Error when creating thread #%d: %d\n", i + 1, rc);
			return (-1);
		}
	}
	pthread_attr_destroy(&tattr);//free attribute

	//Wait until all threads are done (join)
	for (int i = 0; i < NUMTH; i++)
	{
		rc = pthread_join(thread[i], &status);
		if (rc)
		{
			printf("Error when joining thread #%d: %d\n", i + 1, rc);
			return (-1);
		}
	}
	end = dsecnd();
	idx = cloudstr.idx;
	float seconds = end - start;
	fopen_s(&document, "Correspondence_parallel.csv", "w");//open the file
	if (!document) {
		perror("File opening failed");
		return (-1);
	}
	printf("Correspondence:\n");
	for (int i = 0; i < NUM_POINTS; i++)
		printf("%d: %d\n", i + 1, idx[i]);//1-based index (for checking with MATLAB)
	for (int i = 0; i < NUM_POINTS; i++)
		fprintf(document, "%d: %d\n", i + 1, idx[i]);//1-based index (for checking with MATLAB)
	printf("\nElapsed time: %f ms", 1000.0 * seconds);

	//////////////////////////////////////End of Block 4//////////////////////////////////////
	
	free(P), free(Q), free(R), free(T);
	free(idx);
	pthread_exit(NULL);
	return 0;
}