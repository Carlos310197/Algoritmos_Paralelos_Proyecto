/*Este programa recibe un archivo CSV con 64 LiDAR data packets y 
devuelve un vector de 16384 valores en double con informacion de radios de los puntos escaneados*/

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define NUM_POINTS 16384

int main(void)
{
	///////Bloque 1: Abrir y Leer el archivo con 64 LiDAR data packets//////
	int i = 0;
	const int N_LINE = 128;//numero maximo de caracteres a leer en cada linea
	char line[N_LINE];

	FILE* document;
	fopen_s(&document, "Donut_1024x16.csv", "r");//abrir el archivo
	if (!document) {//revisar si fue correctamente abierto
		perror("File opening failed");
		return 0;
	}
	int j = 1;//numero de linea
	int channel = 2;
	int azimuth_block = 0;
	int lidar_packet = 0;
	int word;//indice de la palabra a leer
	int offset = 0;
	float r[NUM_POINTS] = {};//radios
	unsigned long int encoder_count = 0;//para leer el contador inicial del encoder (luego crece en 88 ticks)
	unsigned long int aux = 0;
	while (fgets(line, N_LINE, document) != NULL)
	{
		//obtener el primer valor de encoder_count
		if (j == 13) encoder_count = atoi(line);
		if (j == 14) encoder_count = atoi(line) << 8 | encoder_count;

		//leer radios
		word = 17 + 12 * channel + 788 * azimuth_block + 12608 * lidar_packet;
		if (j == word) aux = (unsigned long int) atoi(line);
		if (j == word + 1) aux = (unsigned long int) (atoi(line) << 8) | aux;
		if (j == word + 2) aux = (unsigned long int) ((atoi(line) & 0x0000000F)<<16) | aux;

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

	//lectura del archivo beam_intrinsics
	fopen_s(&document, "beam_intrinsics.csv", "r");//abrir el archivo
	if (!document) {//revisar si fue correctamente abierto
		perror("File opening failed");
		return 0;
	}

	float altitude[16] = {};
	float azimuth[16] = {};
	j = 1;
	while (fgets(line, N_LINE, document) != NULL)
	{
		//leer altitute angles
		if (j == 2) offset = 0;
		if (j >= 2 && j <= 65)
		{
			if (j % 4 == 0)
			{
				altitude[offset] = (float)atof(line);
				offset++;
			}
		}

		//leer azimuth angles
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

	///////Fin del Bloque 1///////

	///////Bloque 2: Conversion a coordenadas cartesianas///////

	float point_cloud[3 * NUM_POINTS] = {};//xyz coordinates
	float theta = 0.0;//azimuth angle
	float phi = 0.0;//altitude angle
	//data format: x1y1z1 x2y2z2 x3y3z3 ....
	/*channel = 0;
	for (i = 0; i < NUM_POINTS; i++)
	{
		
		theta = (float)(2 * M_PI * ((encoder_count / 90112.0) + (azimuth[channel] / 360.0)));
		phi = (float)(2 * M_PI * altitude[channel] / 360.0);
		point_cloud[0 + 3 * i] = (float)(r[i] * cos(theta) * cos(phi));//x
		point_cloud[1 + 3 * i] = (float)(-r[i] * sin(theta) * cos(phi));//y
		point_cloud[2 + 3 * i] = (float)(r[i] * sin(phi));//z

		channel++;
		if (channel >= 16)//siguiente azimuth block
		{
			channel = 0;
			encoder_count += 88;
		}
		if (encoder_count >= 90112) encoder_count = 0;
	}*/
	unsigned long int counter;
	for (i = 0; i < NUM_POINTS; i++)
	{
		azimuth_block = i / 16;
		counter = (encoder_count + azimuth_block * 88) % 90112;
		channel = i % 16;
		theta = (float)(2 * M_PI * (counter / 90112.0 + azimuth[channel] / 360.0));
		phi = (float)(2 * M_PI * altitude[channel] / 360.0);
		point_cloud[0 + 3 * i] = (float)(r[i] * cos(theta) * cos(phi));//x
		point_cloud[1 + 3 * i] = (float)(-r[i] * sin(theta) * cos(phi));//y
		point_cloud[2 + 3 * i] = (float)(r[i] * sin(phi));//z
	}
	///////Fin del Bloque 2///////

	///////Bloque 3: Escribir los puntos en un documento///////
	fopen_s(&document, "Output_file.csv", "w");
	if (!document) {
		perror("File opening failed");
		return 0;
	}

	for (i = 0; i < NUM_POINTS; i++)
	{
		for (j = 0; j < 2; j++) fprintf(document, "%.3f, ", point_cloud[j + i * 3]);
		fprintf(document, "%.3f\n ", point_cloud[j + i * 3]);
	}
	fprintf(document, "\n");

	fclose(document);

	printf("Success!\n");
	///////Fin del Bloque 3///////
	return 1;
}