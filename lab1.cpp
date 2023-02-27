#include <stdio.h>
#include <math.h>
#pragma warning(disable : 4996)


#define ROW			15	// n2 = ROW - i
#define COL			15	// n1 = j + 1
#define WINDOW_SIZE 3

void PrintMatrix(FILE* file, int* matrix[], const char* name)
{
	fprintf(file, "======= Matrix %s =======\n", name);
	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			fprintf(file, "%d\t", matrix[i * COL + j]);
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n\n");
}

void main(void)
{
	//          		  	 [n1][n2]     [n1][n2]
	// in signal without noise: [4][6] = 305, [12][8] = 483
	const int signal_noise[15][15] =
	{   {12,30,54,42,48,68,103,117,68,62,55,62,56,27,12},
		{6,28,47,76,76,94,135,162,177,145,63,96,104,67,38},
		{12,41,59,94,122,184,206,234,258,208,180,162,116,95,54},
		{28,74,131,128,182,252,304,370,369,286,251,241,201,120,73},
		{29,84,122,153,218,299,393,447,413,380,345,272,206,153,71},
		{26,88,129,178,243,354,880,505,514,449,375,343,304,188,94},
		{17,90,178,256,277,380,507,605,617,542,462,406,298,175,98},
		{50,127,206,267,372,515,634,721,682,601,544,456,360,224,98},
		{23,69,150,219,308,387,535,631,659,542,479,392,276,195,148},
		{22,85,144,231,275,383,474,536,558,517,800,331,277,139,60},
		{35,90,130,196,256,359,403,518,486,387,346,307,232,127,87},
		{28,64,90,132,151,249,300,351,372,331,285,231,162,75,57},
		{18,36,52,103,143,184,245,287,279,228,198,119,118,75,45},
		{0,36,74,102,106,153,191,194,238,145,159,100,97,55,30},
		{0,0,15,16,12,4,1,0,15,25,15,5,0,0,0} 
	};


	int filter[3][3] = { 
				  {1, 2, 1},
				  {2, 5, 2},
				  {1, 2, 1} 
	};

	const int amount = 17;
	const double m = 1.25;

	double nu = 0.0;

	int signal[15][15];
	int delta[15][15];

	FILE* file = fopen("rash3.xls","w");

	// anisotropic filter
	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			signal[i][j] = 0;
			delta[i][j] = 0;
		}
	}

	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			if (i == 0 || i == ROW - 1 || j == 0 || j == COL - 1)
			{
				signal[i][j] = signal_noise[i][j];
			}
			else
			{
				int sg = 0;
				for (int k1 = 0; k1 < WINDOW_SIZE; k1++)
				{
					for (int k2 = 0; k2 < WINDOW_SIZE; k2++)
					{
						sg += filter[k1][k2] * signal_noise[i - 1 + k1][j - 1 + k2];
					}
				}
				signal[i][j] = sg / amount;
			}
		}
	}

	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			delta[i][j] = signal_noise[i][j] - signal[i][j];
		}
	}

	fprintf(file, "Anisotropic filter\n\n");
	PrintMatrix(file, (int**)signal_noise, "Signal-noise");
	PrintMatrix(file, (int**)signal, "Signal");
	PrintMatrix(file, (int**)delta, "Delta (Signal-noise - Signal)");
	fprintf(file, "\n\n");

	// statistic filter (nu = m * sigma)
	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			signal[i][j] = 0;
			delta[i][j] = 0;
		}
	}

	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			int    sum1 = 0;
			double sum2 = 0.0;

			for (int k1 = 0; k1 < WINDOW_SIZE; k1++)
			{
				for (int k2 = 0; k2 < WINDOW_SIZE; k2++)
				{
					sum1 += signal_noise[i - 1 + k1][j - 1 + k2];
				}
			}

			double G = (double)sum1 / pow((double)WINDOW_SIZE, 2);

			for (int k1 = 0; k1 < WINDOW_SIZE; k1++)
			{
				for (int k2 = 0; k2 < WINDOW_SIZE; k2++)
				{
					sum2 += pow((double)signal_noise[i - 1 + k1][j - 1 + k2] - G, 2);
				}
			}

			double D = sum2 / (pow((double)WINDOW_SIZE, 2) - 1);

			nu = m * sqrt(D);

			signal[i][j] = (signal_noise[i][j] - G) < nu ? signal_noise[i][j] : (int)G;
		}
	}

	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			delta[i][j] = signal_noise[i][j] - signal[i][j];
		}
	}

	fprintf(file, "Statistic filter\n\n");
	PrintMatrix(file, (int**)signal_noise, "Signal-noise");
	PrintMatrix(file, (int**)signal, "Signal");
	PrintMatrix(file, (int**)delta, "Delta (Signal-noise - Signal)");
	fprintf(file, "\n\n");

	fclose(file);
}
