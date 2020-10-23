/*Este programa recibe un archivo CSV con 64 LiDAR data packets y 
devuelve un vector de 16384 valores en double con informacion de radios de los puntos escaneados*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
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

int main(void)
{
	///////Bloque 1: Abrir y leer los archivos Donut y beam_intrinsics y  //////
	int i = 0;
	const int N_LINE = 128;//numero maximo de caracteres a leer en cada linea
	char line[N_LINE];

	FILE* document;
	document = fopen("Donut_1024x16.csv", "r");//abrir el archivo
	if (!document) {//revisar si fue correctamente abierto
		perror("File opening failed");
		return 0;
	}

	float* h_r = NULL;//radios
	size_t bytes_r = NUM_POINTS*sizeof(float);
	h_r = (float*)malloc(bytes_r);
	unsigned long int h_encoder_count = 0;//contador inicial del encoder (luego crece en 88 ticks)
	
	int offset = 0;
	unsigned long int word = 0;

	int channel = 2;
	int azimuth_block = 0;
	int lidar_packet = 0;
	int idx_line;//indice palabra a leer
	int j = 1;//numero de linea
	while (fgets(line, N_LINE, document) != NULL)
	{
		//obtener el primer valor de encoder_count
		if (j == 13) h_encoder_count = atoi(line);
		if (j == 14) h_encoder_count = atoi(line) << 8 | h_encoder_count;

		//leer radios del archivo Donut
		idx_line = 17 + 12 * channel + 788 * azimuth_block + 12608 * lidar_packet;
		if (j == idx_line) word = (unsigned long int) atoi(line);
		if (j == idx_line + 1) word = (unsigned long int) (atoi(line) << 8) | word;
		if (j == idx_line + 2) word = (unsigned long int) ((atoi(line) & 0x0000000F)<<16) | word;

		if (j > (idx_line + 2))//si se leyo el radio, pasar al sgte channel
		{
			h_r[offset] = (float)word;
			offset++;
			channel += 4;
		}
		if (channel >= 64)//si se terminaron los channels del bloque, pasar al sgte azimuth block
		{
			channel = 2;
			azimuth_block++;
		}
		if (azimuth_block >= 16)//si se terminaron los azimuth blocks, pasar al sgte lidar packet
		{
			azimuth_block = 0;
			lidar_packet++;
		}
		if (lidar_packet >= 64) break;//si se terminaron los lidar packets, salir
		j++;
	}
	fclose(document);

	//printf("%ld\n",h_encoder_count);
	//for(i=0;i<100;i++) printf("%.3f\n",h_r[i]);

	//lectura del archivo beam_intrinsics
	document = fopen("beam_intrinsics.csv", "r");//abrir el archivo
	if (!document) {//revisar si fue correctamente abierto
		perror("File opening failed");
		return 0;
	}

	float *h_altitude = NULL;
	float *h_azimuth = NULL;
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

	//for(i=0;i<16;i++) printf("%.3f\n",h_altitude[i]);
	//for(i=0;i<16;i++) printf("%.3f\n",h_azimuth[i]);

	///////Fin del Bloque 1///////

	///////Bloque 2: Conversion a coordenadas cartesianas///////

	//reservar memoria para el puntero de salida
	float *h_point_cloud = NULL;
	h_point_cloud = (float*)malloc(3 * bytes_r);

	//declaracion de variables y reserva de memoria en el GPU
	float *d_point_cloud = NULL;//arreglo con los puntos en coordenadas cartesianas 
	//formato: x1y1z1 x2y2z2 x3y3z3 ...
	float *d_r = NULL;
	float *d_azimuth = NULL;
	float *d_altitude = NULL;
	unsigned long int* d_encoder_count;
	cudaMalloc(&d_point_cloud, 3*bytes_r);
	cudaMalloc(&d_r, bytes_r);
	cudaMalloc(&d_azimuth, bytes_angles);
	cudaMalloc(&d_altitude, bytes_angles);
	cudaMalloc(&d_encoder_count, sizeof(unsigned long int));

	//mover data a GPU
	cudaMemcpy(d_r,h_r,bytes_r,cudaMemcpyHostToDevice);
	cudaMemcpy(d_azimuth,h_azimuth,bytes_angles,cudaMemcpyHostToDevice);
	cudaMemcpy(d_altitude,h_altitude,bytes_angles,cudaMemcpyHostToDevice);
	cudaMemcpy(d_encoder_count,&h_encoder_count,sizeof(unsigned long int),cudaMemcpyHostToDevice);

	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	int BlockSize = NUM_POINTS/16;
	int GridSize = 16;
	//lanzar el kernel 
	cudaEventRecord(start);

	Conversion<<<GridSize,BlockSize>>>(d_r, d_encoder_count, d_altitude, d_azimuth, d_point_cloud);
	cudaDeviceSynchronize();

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	printf("Kernel's elapsed time: %.3f ms\n",milliseconds);

	//mover data de salida al CPU
	cudaMemcpy(h_point_cloud, d_point_cloud, 3 * bytes_r, cudaMemcpyDeviceToHost);
	///////Fin del Bloque 2///////

	///////Bloque 3: Escribir los puntos en un documento de salida (Output_file.csv)///////
	//abrir el documento a llenar
	document = fopen("Output_file.csv", "w");
	if (!document) {
		perror("File opening failed");
		return 0;
	}

	//llenar el documento con datos
	for (i = 0; i < NUM_POINTS; i++)
	{
		for (j = 0; j < 2; j++) fprintf(document, "%.4f, ", h_point_cloud[j + i * 3]);
		fprintf(document, "%.4f\n ", h_point_cloud[j + i * 3]);
	}

	fclose(document);

	printf("Success!\n");
	///////Fin del Bloque 3///////

	//liberar memoria
	free(document);
	free(h_r), free(h_altitude), free(h_azimuth), free(h_point_cloud);
	cudaFree(d_r), cudaFree(d_altitude), cudaFree(d_azimuth), cudaFree(d_point_cloud), cudaFree(d_encoder_count);

	return 1;
}
