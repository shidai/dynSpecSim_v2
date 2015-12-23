// Tools to study simulated dynamic spectrum  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "dynSpecTool.h"

int main (int argc, char* argv[])
{
	controlStruct control;

	char fname[1024];   // dynamic spectrum
	char pname[1024];   // read in parameter file
	char oname[1024];   // output file
	int i;
	int index, n;

	float *data;

	// read options
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-p") != 0 )
			{
				n++;
			}
		}
		else if (strcmp(argv[i],"-p")==0)
		{
			strcpy(pname,argv[++i]);
			printf ("Parameters are in %s\n", pname);
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
			printf ("Output are in %s\n", oname);
		}
	}

		
	initialiseControl(&control);
	readParams (pname,oname,&control);

	data = (float*)malloc(sizeof(float)*n*control.nsub*control.nchan);
	
	for (i = index; i < index + n; i++)
	{
		// get the data file name
		strcpy(fname,argv[i]);
		//printf ("File %s\n", fname);
		readData(fname,data,i-index);
	}

	makePlot (data, n, &control);

	free(data);

	return 0;
}

