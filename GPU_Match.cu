/*Este programa recibe un archivo CSV con 64 LiDAR data packets y
devuelve un vector de 16384 valores en double con informacion de radios de los puntos escaneados*/
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define NUM_POINTS 16384

__global__
void Conversion(float* r, unsigned long int* encoder_count, float* altitude, float* azimuth, float* point_cloud)
{
	int azimuth_block, channel;
	unsigned long int counter;
	float theta, phi;

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	azimuth_block = i / 16;
	counter = (encoder_count[0] + azimuth_block * 88) % 90112;
	channel = i % 16;
	theta = (float)(2 * M_PI * (counter / 90112.0 + azimuth[channel] / 360.0));
	phi = (float)(2 * M_PI * altitude[channel] / 360.0);
	point_cloud[0 + 3 * i] = (float)(r[i] * cos(theta) * cos(phi));//x
	point_cloud[1 + 3 * i] = (float)(-r[i] * sin(theta) * cos(phi));//y
	point_cloud[2 + 3 * i] = (float)(r[i] * sin(phi));//z
}

__global__
void RyT(float* R, float* T, float* P, float* Q)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	Q[0 + i * 3] = R[0 + 0 * 3] * P[0 + i * 3] + R[0 + 1 * 3] * P[1 + i * 3] + R[0 + 2 * 3] * P[2 + i * 3] + T[0];
	Q[1 + i * 3] = R[1 + 0 * 3] * P[0 + i * 3] + R[1 + 1 * 3] * P[1 + i * 3] + R[1 + 2 * 3] * P[2 + i * 3] + T[1];
	Q[2 + i * 3] = R[2 + 0 * 3] * P[0 + i * 3] + R[2 + 1 * 3] * P[1 + i * 3] + R[2 + 2 * 3] * P[2 + i * 3] + T[2];
}

__global__
void Match(float* P, float* Q, int q_points, int* idx)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	float min = 100000;
	float d;

	float xp = P[0 + i * 3];
	float yp = P[1 + i * 3];
	float zp = P[2 + i * 3];

	float xq, yq, zq;
	int j;
	for (j = 0; j < q_points / 2; j++)
	{
		xq = Q[0 + j * 3];
		yq = Q[1 + j * 3];
		zq = Q[2 + j * 3];
		d = (xp - xq) * (xp - xq) + (yp - yq) * (yp - yq) + (zp - zq) * (zp - zq);
		if (d < min)
		{
			min = d;
			idx[i] = j;
		}
	}

	for (j = j; j < q_points; j++)
	{
		xq = Q[0 + j * 3];
		yq = Q[1 + j * 3];
		zq = Q[2 + j * 3];
		d = (xp - xq) * (xp - xq) + (yp - yq) * (yp - yq) + (zp - zq) * (zp - zq);
		if (d < min)
		{
			min = d;
			idx[i] = j;
		}
	}
}

int main(void)
{
	///////Block 1: Open and read the file with 64 LiDAR data packets//////
	//int i = 0;
	const int N_LINE = 128;
	char line[N_LINE];

	FILE* document;
	document = fopen("Donut_1024x16.csv", "r");
	if (!document) {
		perror("File opening failed");
		return (-1);
	}

	float* h_r = NULL;//radios
	size_t bytes_r = NUM_POINTS * sizeof(float);
	h_r = (float*)malloc(bytes_r);
	unsigned long int h_encoder_count = 0;//initial encoder counter (then grows with 88 ticks)

	int offset = 0;
	unsigned long int word = 0;

	int channel = 2;
	int azimuth_block = 0;
	int lidar_packet = 0;
	int idx_line;//indice palabra a leer
	int j = 1;//numero de linea
	while (fgets(line, N_LINE, document) != NULL)
	{
		//get the first values of the encoder counter
		if (j == 13) h_encoder_count = atoi(line);
		if (j == 14) h_encoder_count = atoi(line) << 8 | h_encoder_count;

		//read the ranges
		idx_line = 17 + 12 * channel + 788 * azimuth_block + 12608 * lidar_packet;
		if (j == idx_line) word = (unsigned long int) atoi(line);
		if (j == idx_line + 1) word = (unsigned long int) (atoi(line) << 8) | word;
		if (j == idx_line + 2) word = (unsigned long int) ((atoi(line) & 0x0000000F) << 16) | word;

		if (j > (idx_line + 2))//go to next channel
		{
			h_r[offset] = (float)word;
			offset++;
			channel += 4;
		}
		if (channel >= 64)//go to next azimuth block
		{
			channel = 2;
			azimuth_block++;
		}
		if (azimuth_block >= 16)//go to next lidar packet
		{
			azimuth_block = 0;
			lidar_packet++;
		}
		if (lidar_packet >= 64) break;//done
		j++;
	}
	fclose(document);

	//printf("%ld\n",h_encoder_count);
	//for(i=0;i<100;i++) printf("%.3f\n",h_r[i]);

	document = fopen("beam_intrinsics.csv", "r");
	if (!document) {
		perror("File opening failed");
		return (-1);
	}

	float* h_altitude = NULL;
	float* h_azimuth = NULL;
	size_t bytes_angles = 16 * sizeof(float);//16 channels
	h_altitude = (float*)malloc(bytes_angles);
	h_azimuth = (float*)malloc(bytes_angles);

	j = 1;
	while (fgets(line, N_LINE, document) != NULL)
	{
		//leer altitute angles
		if (j == 2) offset = 0;
		if (j >= 2 && j <= 65)
		{
			if (j % 4 == 0)
			{
				h_altitude[offset] = (float)atof(line);
				offset++;
			}
		}

		//leer azimuth angles
		if (j == 68) offset = 0;
		if (j >= 68 && j <= 131)
		{
			if ((j - 66) % 4 == 0)
			{
				h_azimuth[offset] = (float)atof(line);
				offset++;
			}
		}
		j++;
	}
	fclose(document);

	///////End of Block 1///////

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float milliseconds1 = 0;//for "Conversion"
	float milliseconds2 = 0;//for "RyT"

	//Optimize occupancy
	int GridSize = 64;
	int BlockSize = 256;

	cudaError_t err = cudaSuccess;//for checking errors in kernels

	///////Block 2: Conversion to Cartesian coordinates///////

	//allocate memory for storing the cloud (format: x1y1z1 x2y2z2 x3y3z3 ...)
	float* h_P = (float*)malloc(3 * bytes_r);
	float* d_P = NULL;

	float* d_r = NULL;
	float* d_azimuth = NULL;
	float* d_altitude = NULL;
	unsigned long int* d_encoder_count;
	cudaMalloc(&d_P, 3 * bytes_r);
	cudaMalloc(&d_r, bytes_r);
	cudaMalloc(&d_azimuth, bytes_angles);
	cudaMalloc(&d_altitude, bytes_angles);
	cudaMalloc(&d_encoder_count, sizeof(unsigned long int));

	//move data to GPU
	cudaMemcpy(d_r, h_r, bytes_r, cudaMemcpyHostToDevice);
	cudaMemcpy(d_azimuth, h_azimuth, bytes_angles, cudaMemcpyHostToDevice);
	cudaMemcpy(d_altitude, h_altitude, bytes_angles, cudaMemcpyHostToDevice);
	cudaMemcpy(d_encoder_count, &h_encoder_count, sizeof(unsigned long int), cudaMemcpyHostToDevice);

	//Launch "Conversion" kernel
	cudaEventRecord(start);

	Conversion << <GridSize, BlockSize >> > (d_r, d_encoder_count, d_altitude, d_azimuth, d_P);
	err = cudaGetLastError();
	if (err != cudaSuccess) printf("Error in Conversion kernel: %s\n", cudaGetErrorString(err));
	cudaDeviceSynchronize();

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&milliseconds1, start, stop);
	printf("Conversion kernel's elapsed time: %.3f ms\n", milliseconds1);

	///////End of Block 2///////

	///////Block 3: Compute the rotation and translation///////

	//converted cloud
	float* h_Q = (float*)malloc(3 * bytes_r);
	float* d_Q = NULL;
	cudaMalloc(&d_Q, 3 * bytes_r);

	//rotation matrix and translation vector
	float* h_R = (float*)malloc(9 * sizeof(float));
	float* h_T = (float*)malloc(3 * sizeof(float));
	float* d_R, * d_T;
	cudaMalloc(&d_R, 9 * sizeof(float));
	cudaMalloc(&d_T, 3 * sizeof(float));

	//Translation values
	h_T[0] = 0.8f;//x
	h_T[1] = -0.3f;//y
	h_T[2] = 0.2f;//z

	//Rotation values (rad)
	float rx = 0.2f;//axis x
	float ry = -0.2f;//axis y
	float rz = 0.05f;//axis z

	float cx = (float)cos(rx); float cy = (float)cos(ry); float cz = (float)cos(rz);
	float sx = (float)sin(rx); float sy = (float)sin(ry); float sz = (float)sin(rz);
	h_R[0] = cy * cz; h_R[1] = (cz * sx * sy) + (cx * sz); h_R[2] = -(cx * cz * sy) + (sx * sz);
	h_R[3] = -cy * sz; h_R[4] = (cx * cz) - (sx * sy * sz); h_R[5] = (cx * sy * sz) + (cz * sx);
	h_R[6] = sy; h_R[7] = -cy * sx; h_R[8] = cx * cy;

	//move data to GPU
	cudaMemcpy(d_R, h_R, 9 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_T, h_T, 3 * sizeof(float), cudaMemcpyHostToDevice);

	cudaEventRecord(start);

	RyT << <GridSize, BlockSize >> > (d_R, d_T, d_P, d_Q);
	err = cudaGetLastError();
	if (err != cudaSuccess) printf("Error in RyT kernel: %s\n", cudaGetErrorString(err));
	cudaDeviceSynchronize();

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&milliseconds2, start, stop);
	printf("RyT kernel's elapsed time: %.3f ms\n", milliseconds2);
	///////End of Block 3///////

	///////Block 4: Match///////
	int* d_idx;
	cudaMalloc(&d_idx, NUM_POINTS*sizeof(int));
	
	cudaEventRecord(start);
	
	Match << <GridSize, BlockSize >> > (d_P, d_Q, NUM_POINTS, d_idx);
	
	err = cudaGetLastError();
	if (err != cudaSuccess) printf("Error in Match kernel: %s\n", cudaGetErrorString(err));
	cudaDeviceSynchronize();

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milliseconds2, start, stop);
	printf("Match kernel's elapsed time: %.3f ms\n", milliseconds2);
	
	///////End of Block 4///////

	printf("Success!\n");

	//Free variables
	free(h_P), cudaFree(d_P);
	free(h_Q), cudaFree(d_Q);
	free(h_T), free(h_R);
	cudaFree(d_T), cudaFree(d_R);
	free(h_r), free(h_altitude), free(h_azimuth);
	cudaFree(d_r), cudaFree(d_altitude), cudaFree(d_azimuth), cudaFree(d_encoder_count);
	cudaFree(d_idx);
	
	return 0;
}
