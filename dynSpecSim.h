#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
////#include <gsl/gsl_rng.h>
////#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "T2toolkit.h"
#include "cpgplot.h"

typedef struct acfStruct {
	double cFlux;  // pulsar flux density
	double whiteLevel;  // white noise level

	double phaseGradient;
	double cFreq; // observing central frequency
	double bw; // observing bandwidth
	double f0; // scintillation bandwidth
	double tint; // integration time
	double t0; // scintillation time-scale
	int nchn;
	int nsubint;

	int ns; // sampling number of spatial scale 
	int nf; // sampling number of frequency scale
	double size[2]; // sampling boundary
	double steps;
	double stepf;
	double *s; // spatial scale
	double *f; // bw scale
	double *acf2d;  // ACF 
	double *psrt;  // power spactrum
	fftw_complex *eField; // complex electric field
	fftw_complex *intensity;  // intensity 
	double **dynSpec; // dynamic spectrum 
	double **dynSpecWindow; // dynamic spectrum window, nchn*nsubint dimension
	float *dynPlot; // dynamic spectrum for pgplot
} acfStruct;

typedef struct controlStruct {
	char oname[1024]; // output file name

	double cFreq;  // observing central frequency
	double tsub;  // subintegration time
	int nsub;  // number of subintegrations
	double chanBW;  // subchannel bandwidth
	int nchan; // number of subchannels

	double scint_freqbw;   // scintillation bandwidth
	double scint_ts;       // scintillation timescale

	double whiteLevel;   // white noise level, mJy
	double cFlux;        // flux density of pulsars, mJy
}controlStruct;

int idft2d (acfStruct *acfStructure);
int dft2d (acfStruct *acfStructure, fftw_complex *out);
int calACF (acfStruct *acfStructure);
int power (acfStruct *acfStructure);
void deallocateMemory (acfStruct *acfStructure);
void allocateMemory (acfStruct *acfStructure);
int simDynSpec (acfStruct *acfStructure);
int calculateScintScale (acfStruct *acfStructure, controlStruct *control);
void preAllocateMemory (acfStruct *acfStructure);
float find_peak_value (int n, float *s);
int calSize (acfStruct *acfStructure, double *size, double *ratio);
int windowSize (acfStruct *acfStructure, double *size);
int readParams(char *fname, char *oname, controlStruct *control);
void initialiseControl(controlStruct *control);

void heatMap (acfStruct *acfStructure);
void palett(int TYPE, float CONTRA, float BRIGHT);
int plotDynSpec (char *pname);

int calculateScintScale (acfStruct *acfStructure, controlStruct *control)
{
	FILE *fin;
	int i, j;

	printf ("Starting simulating dynamic spectrum\n");
	acfStructure->cFlux = control->cFlux; // mJy
	acfStructure->whiteLevel = control->whiteLevel; // mJy
	acfStructure->cFreq = control->cFreq; // MHz
	acfStructure->bw = fabs(control->chanBW*control->nchan); // MHz
	acfStructure->f0 = control->scint_freqbw;  // MHz
	acfStructure->tint = control->nsub*control->tsub;  // s
	acfStructure->t0 = control->scint_ts; // MHz
	printf ("Scintillation bandwidth: %lf (MHz)\n", acfStructure->f0);
	printf ("Scintillation time-scale: %lf (s)\n", acfStructure->t0);

	acfStructure->nchn = control->nchan;
	acfStructure->nsubint = control->nsub;

	preAllocateMemory (acfStructure);
	allocateMemory (acfStructure);

	calACF (acfStructure);
	power (acfStructure);
	simDynSpec (acfStructure);

	if ((fin=fopen(control->oname,"w"))==NULL)
	{
		printf ("Can't open output file!\n");
		exit(1);
	}

	fprintf(fin,"INFO nsub nchn bandwidth tint cFreq\n");
	fprintf(fin,"START %d %d %lf %lf %lf\n",acfStructure->nsubint,acfStructure->nchn,acfStructure->bw,acfStructure->tint,acfStructure->cFreq);

	for (i=0;i<acfStructure->nchn;i++)
	{
		for (j=0;j<acfStructure->nsubint;j++)
		{
			fprintf(fin,"%d %d %lf\n", i, j, acfStructure->dynPlot[i*acfStructure->nsubint+j]);
		}
	}

	if (fclose(fin))
	{
		printf ("Can't close output file!\n");
		exit(1);
	}

	return 0;
}

int calACF (acfStruct *acfStructure)
{
	printf ("Calculating ACF\n");
	int i,j;
	int ns = acfStructure->ns;
	int nf = acfStructure->nf;
	double steps = acfStructure->steps;
	double stepf = acfStructure->stepf;
	double *acf;
	acf = (double *)malloc(sizeof(double)*ns*nf);

	for (i = 0; i < ns; i++)
	{
		acfStructure->s[i] = -acfStructure->size[1]+i*steps;
	}

      	for (i = 0; i < nf; i++)
	{
		acfStructure->f[i] = -acfStructure->size[0]+i*stepf;
	}

	double rand;
	rand = acfStructure->phaseGradient;
	//printf ("%lf\n",rand);
	
	int n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			//acf[n] = exp(-pow((pow(acfStructure->s[j],2.5)+pow(acfStructure->f[i],1.5)),2.0/3.0));
			//acf[n] = exp(-pow((pow(fabs(acfStructure->s[j]),2.5)+pow(fabs(acfStructure->f[i]),1.5)),2.0/3.0));
			acf[n] = exp(-pow((pow(fabs(acfStructure->s[j]+2.0*rand*0.4*acfStructure->f[i]),2.5)+pow(fabs(acfStructure->f[i]),1.5)),2.0/3.0));
			n++;
		}
	}

	//////////////////////////////////////////////////////////////////

	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			if (i >= (int)ceil(nf/2.0) && j >= (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i-(int)ceil(nf/2.0))+(j-(int)ceil(ns/2.0))];
			}
			else if  (i >= (int)ceil(nf/2.0) && j < (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i-(int)ceil(nf/2.0))+(j+(int)floor(ns/2.0))];
			}
			else if  (i < (int)ceil(nf/2.0) && j < (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i+(int)floor(nf/2.0))+(j+(int)floor(ns/2.0))];
			}
			else if  (i < (int)ceil(nf/2.0) && j >= (int)ceil(ns/2.0))
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i+(int)floor(nf/2.0))+(j-(int)ceil(ns/2.0))];
			}
		}
	}

	free(acf);

	return 0;
}


int dft2d (acfStruct *acfStructure, fftw_complex *out)
{
	int i;
	int n0 = acfStructure->nf;
	int n1 = acfStructure->ns;
	double *in;
	in = (double *)malloc(sizeof(double)*n0*n1);

	fftw_plan p;
	
	p = fftw_plan_dft_r2c_2d (n0, n1, in, out, FFTW_ESTIMATE);
	//p = fftw_plan_dft_r2c_2d (n0, n1, in, out, FFTW_MEASURE);

	for (i = 0; i < n0*n1; i++)
	{
		in[i] = acfStructure->acf2d[i];
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	free (in);
  
	return 0;
}

int idft2d (acfStruct *acfStructure)
{
	int i;
	int n0 = acfStructure->nf;
	int n1 = acfStructure->ns;

	fftw_complex *in;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n0*n1);

	fftw_plan p;
	
	p = fftw_plan_dft_2d (n0, n1, in, acfStructure->intensity, FFTW_BACKWARD, FFTW_ESTIMATE);
	//p = fftw_plan_dft_2d (n0, n1, in, acfStructure->intensity, FFTW_BACKWARD, FFTW_MEASURE);

	for (i = 0; i < n0*n1; i++)
	{
		in[i][0] = acfStructure->eField[i][0];
		in[i][1] = acfStructure->eField[i][1];
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	fftw_free (in);
  
	return 0;
}

int power (acfStruct *acfStructure)
{
	int i;
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;
	/////////////////////////////////////////////////////////////////////////////////

  fftw_complex *out;
	
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*((int)(ns/2)+1));
	
	dft2d (acfStructure, out);

	int n, j;
	n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			if (j < (int)(ns/2)+1)
			{
				acfStructure->psrt[n] = sqrt(sqrt(pow(out[i*((int)(ns/2)+1)+j][0],2.0)+pow(out[i*((int)(ns/2)+1)+j][1],2.0)));
			}
			else
			{
				if (i == 0)
				{
					acfStructure->psrt[n] = sqrt(sqrt(pow(out[i*((int)(ns/2)+1)+ns-j][0],2.0)+pow(out[i*((int)(ns/2)+1)+ns-j][1],2.0)));
				}
				else
				{
					acfStructure->psrt[n] = sqrt(sqrt(pow(out[(nf-i)*((int)(ns/2)+1)+ns-j][0],2.0)+pow(out[(nf-i)*((int)(ns/2)+1)+ns-j][1],2.0)));
				}
			}
			//printf ("%lf ", acfStructure->psrt[n]);
			n++;
		}
		//printf ("\n");
	}

	fftw_free(out); 

	return 0;
}

void allocateMemory (acfStruct *acfStructure)
{
	int ns, nf;
	int nchn, nsubint;

	ns = acfStructure->ns;
	nf = acfStructure->nf;
	nchn = acfStructure->nchn;
	nsubint = acfStructure->nsubint;
	
	acfStructure->s = (double *)malloc(sizeof(double)*ns);
	acfStructure->f = (double *)malloc(sizeof(double)*nf);
	acfStructure->acf2d = (double *)malloc(sizeof(double)*ns*nf);
	acfStructure->psrt = (double *)malloc(sizeof(double)*ns*nf);
	acfStructure->eField = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*ns);
	acfStructure->intensity = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*ns);
	
	acfStructure->dynSpec = (double **)fftw_malloc(sizeof(double *)*nf);
	acfStructure->dynSpecWindow = (double **)fftw_malloc(sizeof(double *)*nchn);

	int i;
	for (i = 0; i < nf; i++)
	{
		acfStructure->dynSpec[i] = (double *)fftw_malloc(sizeof(double)*ns);
	}

	for (i = 0; i < nchn; i++)
	{
		acfStructure->dynSpecWindow[i] = (double *)fftw_malloc(sizeof(double)*nsubint);
	}

	acfStructure->dynPlot = (float *)malloc(sizeof(float)*nsubint*nchn);
}

void deallocateMemory (acfStruct *acfStructure)
{
	int nf = acfStructure->nf;
	int nchn = acfStructure->nchn;
	
	fftw_free(acfStructure->s); 
	fftw_free(acfStructure->f); 
	fftw_free(acfStructure->eField); 
	fftw_free(acfStructure->intensity); 
	free(acfStructure->acf2d);
	free(acfStructure->psrt);
	free(acfStructure->dynPlot);

	int i;
	for (i = 0; i < nf; i++)
	{
		free(acfStructure->dynSpec[i]);
	}

	for (i = 0; i < nchn; i++)
	{
		free(acfStructure->dynSpecWindow[i]);
	}

	free(acfStructure->dynSpec);
}

int simDynSpec (acfStruct *acfStructure)
{
	printf ("Simulating dynamic spectrum\n");
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;
	int nchn = acfStructure->nchn;
	int nsubint = acfStructure->nsubint;

	long seed;
	
	int i;
	int j;
	seed = TKsetSeed();
	//printf ("seed %ld\n",seed);

	for (i = 0; i < nf*ns; i++)
	{
		//acfStructure->eField[i][0] = acfStructure->psrt[i];
		//acfStructure->eField[i][1] = acfStructure->psrt[i];
		acfStructure->eField[i][0] = acfStructure->psrt[i]*TKgaussDev(&seed);
		acfStructure->eField[i][1] = acfStructure->psrt[i]*TKgaussDev(&seed);
		//printf ("%lf\n",TKgaussDev(&seed));
	}

	/////////////////////////////////////////////////////////////////////////////////

	// ifft
	idft2d (acfStructure);

	// form the matrix and normalize
	int n = 0;
	double sum = 0.0;

	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			acfStructure->dynSpec[i][j] = pow(acfStructure->intensity[n][1]/(nf*ns),2.0)+pow(acfStructure->intensity[n][0]/(nf*ns),2.0);
			//fprintf (fp, "%lf  ", acfStructure->dynSpec[i][j]);
			sum += pow(acfStructure->intensity[n][1]/(nf*ns),2.0)+pow(acfStructure->intensity[n][0]/(nf*ns),2.0);
			n++;
		}
		//fprintf (fp, "\n");
	}

	sum = sum/n;
	printf ("Normalization %.10lf\n",sum);

	// choose a subwindow
	double rand, rand2;
	//printf ("seed %ld\n",seed);
	rand = TKgaussDev(&seed);
	rand2 = rand - floor(rand);
	//printf ("rand %lf\n",rand);

	int nf0 = (int)(rand2*(nf-nchn));
	int ns0 = (int)(rand2*(ns-nsubint));
	
	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < nsubint; j++)
		{
			acfStructure->dynSpecWindow[i][j] = acfStructure->dynSpec[i+nf0][j+ns0]/sum;
			acfStructure->dynPlot[i*nsubint+j] = (float)(acfStructure->dynSpecWindow[i][j]*acfStructure->cFlux+acfStructure->whiteLevel*TKgaussDev(&seed));   // add in noise here
			//printf ("noise rand %lf\n",TKgaussDev(&seed));
			//acfStructure->dynPlot[i*nsubint+j] = (float)(acfStructure->dynSpecWindow[i][j]);
			//fprintf (fp, "%.10lf  ", acfStructure->dynSpec[i][j]/sum);
		}
	}

	return 0;
}

int windowSize (acfStruct *acfStructure, double *size)
{
	double bw, f0, tint, t0;
	bw = acfStructure->bw;
	f0 = acfStructure->f0;
	tint = acfStructure->tint;
	t0 = acfStructure->t0;

	if ( (bw/f0) > 6 )
	{
		size[0] = bw/f0;
	}
	else
	{
		size[0] = 6.0;
	}
			  
	if ( (tint/t0) > 6 )
	{
		size[1] = tint/t0;
	}
	else
	{
		size[1] = 6.0;
	}

	double ratio[2];
	calSize (acfStructure, size, ratio);
	//printf ("f0 ratio: %lf\n", ratio[0]);
	//printf ("s0 ratio: %lf\n", ratio[1]);

	while (ratio[0] >= 1e-6 || ratio[1] >= 1e-6)
	{
		size[0] = 1.05*size[0];
		size[1] = 1.05*size[1];
		calSize (acfStructure, size, ratio);
		//printf ("f0 ratio: %lf\n", ratio[0]);
		//printf ("s0 ratio: %lf\n", ratio[1]);
	}

	acfStructure->size[0] = size[0];
	acfStructure->size[1] = size[1];

	return 0;
}

int calSize (acfStruct *acfStructure, double *size, double *ratio)
{
	int i;
	double rand = acfStructure->phaseGradient;
	double steps = acfStructure->steps;
	double stepf = acfStructure->stepf;

	int nf = (int)(size[0]*2/stepf)+1;
	int ns = (int)(size[1]*2/steps)+1;

	float  s[ns], acfs[ns], smax;
	float  f[nf], acff[nf], fmax;
	//double c; // value at the center

	for (i = 0; i < ns; i++)
	{
		s[i] = -size[1]+i*steps;
	}

  for (i = 0; i < nf; i++)
	{
		f[i] = -size[0]+i*stepf;
	}

	//c = exp(-pow((pow(fabs(s[(int)(ns/2)]+2.0*rand*0.4*f[(int)(nf/2)]),2.5)+pow(fabs(f[(int)(nf/2)]),1.5)),2.0/3.0));

	for (i = 0; i < nf; i++)
	{
		acff[i] = exp(-pow((pow(fabs(s[0]+2.0*rand*0.4*f[i]),2.5)+pow(fabs(f[i]),1.5)),2.0/3.0));
	}

	for (i = 0; i < ns; i++)
	{
		acfs[i] = exp(-pow((pow(fabs(s[i]+2.0*rand*0.4*f[0]),2.5)+pow(fabs(f[0]),1.5)),2.0/3.0));
	}

	smax = find_peak_value (ns, acfs);
	fmax = find_peak_value (nf, acff);

	ratio[0] = fmax;
	ratio[1] = smax;

	return 0;
}

float find_peak_value (int n, float *s)
{
	int i;
	float temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	
	float a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
															
		c = (a >= b ? a : b);
		
		temp[i+1] = c;
	
	}
	
	return temp[n-1];
}
															
void preAllocateMemory (acfStruct *acfStructure)
{
	long seed;
	//seed = -1420603014;
	seed = TKsetSeed();
	printf ("Seed %ld\n",seed);
	acfStructure->phaseGradient = TKgaussDev(&seed);
	//acfStructure->phaseGradient = 0.0;
	printf ("Phase gradient: %lf\n", acfStructure->phaseGradient);

	double bw, f0, tint, t0;
	int nchn, nsubint;

	bw = acfStructure->bw;
	f0 = acfStructure->f0;
	tint = acfStructure->tint;
	t0 = acfStructure->t0;

	nchn = acfStructure->nchn;
	nsubint = acfStructure->nsubint;

	double steps = (tint/t0)/nsubint;
	double stepf = (bw/f0)/nchn;
	acfStructure->steps = steps;
	acfStructure->stepf = stepf;
	//printf ("%lf\n", steps);
	//printf ("%lf\n", stepf);
	
	double size[2]; // sampling boundary
	windowSize (acfStructure, size);
	printf ("f0 size: %lf\n", size[0]);
	printf ("s0 size: %lf\n", size[1]);

	int nf = (int)(size[0]*2/stepf)+1;
	int ns = (int)(size[1]*2/steps)+1;
	//printf ("%d\n", ns);
	//printf ("%d\n", nf);

	acfStructure->ns = ns;
	acfStructure->nf = nf;
}

int readParams(char *fname, char *oname, controlStruct *control)
{
	FILE *fin;
	char param[1024];
	int endit=-1;
	int finished=0;

	// define the output file name
	strcpy(control->oname,oname);

	///////////////////////////////////////
	if ((fin=fopen(fname,"r"))==NULL)
	{
		printf ("Can't open file!\n");
		exit(1);
	}

	printf("Reading parameters...\n");

      	// Find the start observation
	while (!feof(fin))
	{
		if (fscanf(fin,"%s",param)==1)
		{
			if (strcasecmp(param,"START_OBS")==0)
			{
				endit=0;
				break;
			}
		}
		else 
			return 1;
	}

	printf("Got to this bit with %d\n",endit);
	if (endit==-1)
		return 1;

	do
	{
		fscanf(fin,"%s",param);
		if (strcasecmp(param,"END_OBS")==0)
			endit=1;
		else
		{
			//if (strcasecmp(param,"PHEAD")==0)
			//	fscanf(fin,"%s",control->primaryHeaderParams);
			//else if (strcasecmp(param,"SRC")==0)
			//     	fscanf(fin,"%s",control->src);
			//else if (strcasecmp(param,"EXACT_EPHEMERIS")==0)
			//	fscanf(fin,"%s",control->exact_ephemeris);
			//else if (strcasecmp(param,"TEMPLATE")==0)
			//	fscanf(fin,"%s",control->template);
			if (strcasecmp(param,"SCINT_TS")==0)
				fscanf(fin,"%lf",&(control->scint_ts));
			else if (strcasecmp(param,"SCINT_FREQBW")==0)
				fscanf(fin,"%lf",&(control->scint_freqbw));	  
			//else if (strcasecmp(param,"FILE")==0)
			//	fscanf(fin,"%s",control->fname);
			//else if (strcasecmp(param,"TYPE")==0)
			//	fscanf(fin,"%s",control->type);
			//else if (strcasecmp(param,"STT_IMJD")==0)
			//	fscanf(fin,"%d",&(control->stt_imjd));
			//else if (strcasecmp(param,"STT_SMJD")==0)
			//	fscanf(fin,"%lf",&(control->stt_smjd));
			//else if (strcasecmp(param,"STT_OFFS")==0)
      			//	fscanf(fin,"%lf",&(control->stt_offs));
			else if (strcasecmp(param,"TSUB")==0)
      				fscanf(fin,"%lf",&(control->tsub));
			else if (strcasecmp(param,"CFREQ")==0)
      				fscanf(fin,"%lf",&(control->cFreq));
			else if (strcasecmp(param,"CHAN_BW")==0)
      				fscanf(fin,"%lf",&(control->chanBW));
			else if (strcasecmp(param,"NCHAN")==0)
			      	fscanf(fin,"%d",&(control->nchan));
			//else if (strcasecmp(param,"NBIN")==0)
      			//	fscanf(fin,"%d",&(control->nbin));
			//else if (strcasecmp(param,"NPOL")==0)
      			//	fscanf(fin,"%d",&(control->npol));
			else if (strcasecmp(param,"NSUB")==0)
      				fscanf(fin,"%d",&(control->nsub));
			//else if (strcasecmp(param,"SEGLENGTH")==0)
      			//	fscanf(fin,"%lf",&(control->segLength));
			//else if (strcasecmp(param,"NFREQ_COEFF")==0)
      			//	fscanf(fin,"%d",&(control->nfreqcoeff));
			//else if (strcasecmp(param,"NTIME_COEFF")==0)
			//      	fscanf(fin,"%d",&(control->ntimecoeff));
			else if (strcasecmp(param,"WHITE_LEVEL")==0)
      				fscanf(fin,"%lf",&(control->whiteLevel));
			//else if (strcasecmp(param,"TSYS")==0)
      			//	fscanf(fin,"%lf",&(control->tsys));
			//else if (strcasecmp(param,"TSKY")==0)
      			//	fscanf(fin,"%lf",&(control->tsky));
			//else if (strcasecmp(param,"GAIN")==0)
      			//	fscanf(fin,"%lf",&(control->gain));
			else if (strcasecmp(param,"CFLUX")==0)
			      	fscanf(fin,"%lf",&(control->cFlux));
			//else if (strcasecmp(param,"SI")==0)
			//	fscanf(fin,"%lf",&(control->si));
			//else if (strcasecmp(param,"phaseResolvedSI")==0)
			//{
			//	fscanf(fin,"%s",control->phaseResolvedSI);
			//	control->simProf = 1;
			//}
		}
	} while (endit==0);

	//////////////////////////////////////////////////////////////////////
	//if (control->tsys != 0.0 && control->tsky != 0.0 && control->gain != 0.0 && control->whiteLevel == 0)
	//{
	//	control->radioNoise = (control->tsys+control->tsky)/(control->gain)/sqrt(2.0*(control->tsubRequested/control->nbin)*(fabs(control->obsBW)/control->nchan));
	//}
	//else if (control->tsys == 0.0 && control->tsky == 0.0 && control->gain == 0.0 && control->whiteLevel != 0)
	//{
	//	control->radioNoise = control->whiteLevel;
	//}
	//else 
	//{
	//	printf ("Double definiation of radio-meter noise!\n");
	//	exit (1);
	//}
	//printf ("Nchan: %d; Tsys: %lf; Tsky: %lf; Gain: %lf; Radio-meter noise: %lf mJy\n", control->nchan, control->tsys, control->tsky, control->gain, control->radioNoise);

	if (fclose(fin))
	{
		printf ("Can't close file.\n");
		exit(1);
	}

	return finished;
}

void initialiseControl(controlStruct *control)
{
	//strcpy(control->primaryHeaderParams,"UNKNOWN");
	//strcpy(control->exact_ephemeris,"UNKNOWN");
	//strcpy(control->src,"UNKNOWN");
	strcpy(control->oname,"UNKNOWN");
	control->tsub = 0;
	
	// Standard defaults
	//strcpy(control->type,"PSR");
	//control->nbin = 128;
	control->nchan = 32;
	//control->npol = 1;
	control->nsub = 8;
	control->cFreq = 1400.0;
	control->chanBW = -1;   // MHz
	//control->segLength = 48000;
	//control->nfreqcoeff = 16;
	//control->ntimecoeff = 16;
	//control->stt_imjd = 55000;
	//control->stt_smjd = 5234.0;
	//control->stt_offs = 0.1234;
	control->tsub = 120;  // second
	control->whiteLevel = 0;   // mJy
	control->scint_ts  = 0.0;  // second
	control->scint_freqbw = 0.0;   // MHz
	//control->tsys = 0.0;
	//control->tsky = 0.0;
	//control->gain = 0.0;
	control->cFlux = 0.0;   // mJy
	//control->si = 0.0;
	//control->radioNoise = 0.0;
	//control->bat = 0;
	//control->simProf = 0; // default: do not simulate profile with phase-resolved SI
  // Note that the DM comes from the ephemeris
}

void heatMap (acfStruct *acfStructure)
{
	//int i,j;                     
	//int dimx = acfStructure.ns;
	//int dimy = acfStructure.nf; // dimensions 
	//float tab[dimx*dimy];       // value
	char caption[1024];
	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", acfStructure->cFreq, "MHz", "BW:", acfStructure->bw, "MHz", "Length:", acfStructure->tint, "s");
  
	float zmin,zmax;            /* min et max des valeurs de la fonction */
	float tr[6];                /* matrice utilisee par pgimag */

	int dimx = acfStructure->nsubint;
	int dimy = acfStructure->nchn;
	double bw = acfStructure->bw;
  

	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};

	double f1 = acfStructure->cFreq-bw/2.0-2.0*bw/dimy; // MHz
	double f2 = acfStructure->cFreq+bw/2.0-2.0*bw/dimy; // MHz
	//printf ("f1 f2: %lf %lf\n", f1, f2);

	zmin=0; 
	zmax=find_peak_value (dimx*dimy, acfStructure->dynPlot);
	//double f1 = 1241; // MHz
	//double f2 = 1497; // MHz
	/*The transformation matrix TR is used to calculate the world
	coordinates of the center of the "cell" that represents each
	array element. The world coordinates of the center of the cell
	corresponding to array element A(I,J) are given by:
	X = TR(1) + TR(2)*I + TR(3)*J
	Y = TR(4) + TR(5)*I + TR(6)*J
	Usually TR(3) and TR(5) are zero -- unless the coordinate
	transformation involves a rotation or shear.  The corners of the
	quadrilateral region that is shaded by PGIMAG are given by
	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
  
	//tr[0]=0;
	//tr[1]=(float)(dimy)/dimx;
	//tr[2]=0;
	//tr[3]=0;
	//tr[4]=0;
	//tr[5]=1;
  
	tr[0]=-0.5;
       	tr[1]=1;
       	tr[2]=0;
      	tr[3]=f2+0.5;
      	tr[4]=0;
      	tr[5]=-bw/dimy;

	// plot 
	//cpgbeg(0,"?",1,1);
	cpgbeg(0,"2/xs",1,1);
      	cpgsch(1.2); // set character height
      	cpgscf(2); // set character font
	//cpgswin(0,dimx,164,132); // set window
	cpgswin(0,dimx,f2,f1); // set window
	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
	//cpgbox("BCTSIN",10,5,"BCTSIN",50,5);
      	cpglab("Subintegration","Frequency (MHz)",caption);
      	//cpglab("Subintegration","Frequency (MHz)","Freq: 150.0 MHz BW: -32.000 MHz Length: 960.0 s");
	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
	//palett(3, -0.4, 0.3);
	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	cpgimag(acfStructure->dynPlot,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
	//cpggray(acfStructure->dynPlot,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);

	cpgend();
} 

void palett(int TYPE, float CONTRA, float BRIGHT)
{
//-----------------------------------------------------------------------
// Set a "palette" of colors in the range of color indices used by
// PGIMAG.
//-----------------------------------------------------------------------
	float GL[] = {0.0, 1.0};
	float GR[] = {0.0, 1.0};
	float GG[] = {0.0, 1.0};
	float GB[] = {0.0, 1.0};
	float RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
	float RR[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
	float RG[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
	float RB[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
	float HL[] = {0.0, 0.2, 0.4, 0.6, 1.0};
	float HR[] = {0.0, 0.5, 1.0, 1.0, 1.0};
	float HG[] = {0.0, 0.0, 0.5, 1.0, 1.0};
	float HB[] = {0.0, 0.0, 0.0, 0.3, 1.0};
	float WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
	float WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
	float WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
	float WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
	float AL[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
	float AR[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	float AG[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
	float AB[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      
	if (TYPE == 1)
	{   
		//-- gray scale
		cpgctab(GL, GR, GG, GB, 2, CONTRA, BRIGHT);
	}
	else if (TYPE == 2) 
	{
		//-- rainbow
		cpgctab(RL, RR, RG, RB, 9, CONTRA, BRIGHT);
	}
	else if (TYPE == 3) 
	{
		//-- heat
		cpgctab(HL, HR, HG, HB, 5, CONTRA, BRIGHT);
	}
	else if (TYPE == 4) 
	{
		//-- weird IRAF
		cpgctab(WL, WR, WG, WB, 10, CONTRA, BRIGHT);
	}
	else if (TYPE == 5) 
	{
		//-- AIPS
		cpgctab(AL, AR, AG, AB, 20, CONTRA, BRIGHT);
	}
}

int plotDynSpec (char *pname)
{
	FILE *fin;
	char start[128];
	double tint,bw,cFreq;
	int nsub,nchn;
	float val;
	int n1, n2;
	float *dynSpec;

	int i;

	int dimx;
	int dimy;
	float zmin,zmax;            /* min et max des valeurs de la fonction */
	float tr[6];                /* matrice utilisee par pgimag */

	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
	char caption[1024];

	if ((fin=fopen(pname,"r"))==NULL)
	{
		printf ("Can't open dynamic spectrum!\n");
		exit(1);
	}

      	// Find the start of dynamic spectrum, which contains basic info
	while (!feof(fin))
	{
		if (fscanf(fin,"%s %d %d %lf %lf %lf",start,&nsub,&nchn,&bw,&tint,&cFreq)==6)
		{
			if (strcasecmp(start,"START")==0)
			{
				break;
			}
		}
		//else 
		//	return 1;
	}

	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", cFreq, "MHz", "BW:", bw, "MHz", "Length:", tint, "s");

	dynSpec = (float*)malloc(sizeof(float)*nsub*nchn);

	i = 0;
	while (fscanf(fin,"%d %d %f", &n1, &n2, &val)==3)
	{
		dynSpec[i] = val;
		i++;
	}

	if (fclose(fin))
	{
		printf ("Can't close dynamic spectrum!\n");
		exit(1);
	}

	dimx = nsub;
	dimy = nchn;
  

	zmin=0; 
	zmax=find_peak_value (dimx*dimy, dynSpec);

	double f1 = cFreq-bw/2.0-2.0*bw/dimy; // MHz
	double f2 = cFreq+bw/2.0-2.0*bw/dimy; // MHz
	//double f1 = 1241; // MHz
	//double f2 = 1497; // MHz
	/*The transformation matrix TR is used to calculate the world
	coordinates of the center of the "cell" that represents each
	array element. The world coordinates of the center of the cell
	corresponding to array element A(I,J) are given by:
	X = TR(1) + TR(2)*I + TR(3)*J
	Y = TR(4) + TR(5)*I + TR(6)*J
	Usually TR(3) and TR(5) are zero -- unless the coordinate
	transformation involves a rotation or shear.  The corners of the
	quadrilateral region that is shaded by PGIMAG are given by
	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
  
	//tr[0]=0;
	//tr[1]=(float)(dimy)/dimx;
	//tr[2]=0;
	//tr[3]=0;
	//tr[4]=0;
	//tr[5]=1;
  
	tr[0]=-0.5;
       	tr[1]=1;
       	tr[2]=0;
      	tr[3]=f2+0.5;
      	tr[4]=0;
      	tr[5]=-bw/dimy;
 
	// plot 
	//cpgbeg(0,"?",1,1);
	cpgbeg(0,"2/xs",1,1);
      	cpgsch(1.2); // set character height
      	cpgscf(2); // set character font
	cpgswin(0,dimx,f2,f1); // set window
	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
      	cpglab("Subintegration","Frequency (MHz)",caption);
	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
	//palett(3, -0.4, 0.3);
	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
	//cpggray(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	cpgimag(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);

	cpgend();

	free(dynSpec);

	return 0;
} 

