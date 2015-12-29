#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
////#include <gsl/gsl_rng.h>
////#include <gsl/gsl_randist.h>
//#include "T2toolkit.h"
#include "cpgplot.h"
#include "T2toolkit.h"

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

int readParams(char *fname, controlStruct *control);
void initialiseControl(controlStruct *control);

void palett(int TYPE, float CONTRA, float BRIGHT);
int readData (char *fname, float *data, int n);
int makePlot (float *data, int n, controlStruct *control, char *dname);
int histogram (float *data, int n, float *x, float *val, float low, float up, int step);

int readParams(char *fname, controlStruct *control)
{
	FILE *fin;
	char param[1024];
	int endit=-1;
	int finished=0;

	// define the output file name
	//strcpy(control->oname,oname);

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

//int plotDynSpec (char *pname)
//{
//	FILE *fin;
//	char start[128];
//	double tint,bw,cFreq;
//	int nsub,nchn;
//	float val;
//	int n1, n2;
//	float *dynSpec;
//
//	int i;
//
//	int dimx;
//	int dimy;
//	float zmin,zmax;            /* min et max des valeurs de la fonction */
//	float tr[6];                /* matrice utilisee par pgimag */
//
//	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
//	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
//	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
//	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
//	char caption[1024];
//
//	if ((fin=fopen(pname,"r"))==NULL)
//	{
//		printf ("Can't open dynamic spectrum!\n");
//		exit(1);
//	}
//
//      	// Find the start of dynamic spectrum, which contains basic info
//	while (!feof(fin))
//	{
//		if (fscanf(fin,"%s %d %d %lf %lf %lf",start,&nsub,&nchn,&bw,&tint,&cFreq)==6)
//		{
//			if (strcasecmp(start,"START")==0)
//			{
//				break;
//			}
//		}
//		//else 
//		//	return 1;
//	}
//
//	sprintf (caption, "%s %.2f %s %s %.2f %s %s %.2f %s", "Freq:", cFreq, "MHz", "BW:", bw, "MHz", "Length:", tint, "s");
//
//	dynSpec = (float*)malloc(sizeof(float)*nsub*nchn);
//
//	i = 0;
//	while (fscanf(fin,"%d %d %f", &n1, &n2, &val)==3)
//	{
//		dynSpec[i] = val;
//		i++;
//	}
//
//	if (fclose(fin))
//	{
//		printf ("Can't close dynamic spectrum!\n");
//		exit(1);
//	}
//
//	dimx = nsub;
//	dimy = nchn;
//  
//
//	zmin=0; 
//	zmax=find_peak_value (dimx*dimy, dynSpec);
//
//	double f1 = cFreq-bw/2.0-2.0*bw/dimy; // MHz
//	double f2 = cFreq+bw/2.0-2.0*bw/dimy; // MHz
//	//double f1 = 1241; // MHz
//	//double f2 = 1497; // MHz
//	/*The transformation matrix TR is used to calculate the world
//	coordinates of the center of the "cell" that represents each
//	array element. The world coordinates of the center of the cell
//	corresponding to array element A(I,J) are given by:
//	X = TR(1) + TR(2)*I + TR(3)*J
//	Y = TR(4) + TR(5)*I + TR(6)*J
//	Usually TR(3) and TR(5) are zero -- unless the coordinate
//	transformation involves a rotation or shear.  The corners of the
//	quadrilateral region that is shaded by PGIMAG are given by
//	applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).*/
//  
//	//tr[0]=0;
//	//tr[1]=(float)(dimy)/dimx;
//	//tr[2]=0;
//	//tr[3]=0;
//	//tr[4]=0;
//	//tr[5]=1;
//  
//	tr[0]=-0.5;
//       	tr[1]=1;
//       	tr[2]=0;
//      	tr[3]=f2+0.5;
//      	tr[4]=0;
//      	tr[5]=-bw/dimy;
// 
//	// plot 
//	//cpgbeg(0,"?",1,1);
//	cpgbeg(0,"2/xs",1,1);
//      	cpgsch(1.2); // set character height
//      	cpgscf(2); // set character font
//	cpgswin(0,dimx,f2,f1); // set window
//	//cpgsvp(0.1,0.9,0.1,0.9); // set viewport
//      	//cpgenv(1,dimx,f1,f2,0,0); // set window and viewport and draw labeled frame
//	cpgbox("BCTSIN",4,4,"BCTSIN",16,8);
//      	cpglab("Subintegration","Frequency (MHz)",caption);
//	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
//	//palett(3, -0.4, 0.3);
//	//cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	//cpggray(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
//	cpgimag(dynSpec,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
//
//	cpgend();
//
//	free(dynSpec);
//
//	return 0;
//} 

int readData (char *fname, float *data, int n)
{
	FILE *fin;
	char start[128];
	double tint,bw,cFreq;
	int nsub,nchn;
	float val;
	int n1, n2;

	int i;

	if ((fin=fopen(fname,"r"))==NULL)
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

	i = 0;
	while (fscanf(fin,"%d %d %f", &n1, &n2, &val)==3)
	{
		data[n*nsub*nchn+i] = val;
		i++;
	}

	if (fclose(fin))
	{
		printf ("Can't close dynamic spectrum!\n");
		exit(1);
	}

	return 0;
} 

int makePlot (float *data, int n, controlStruct *control, char *dname)
{
	int i, j;
	int nsub = control->nsub;
	int nchan = control->nchan;
	float sum;
	//double bw = control->chanBW*control->nchan;
	//double cFreq = control->cFreq;

	float *flux; // flux densities
	float *x;    // x axis of the flux variation
	float *xHis; // x axis of the histogram
	float *val;  // data value of the histogram
	int step = 500; // steps in the histogram

	long seed;
	float *noise; // white noise
	float *xHisN; // x axis of the histogram
	float *valN;  // noise value of the histogram

	flux = (float*)malloc(sizeof(float)*n);
	x = (float*)malloc(sizeof(float)*n);
	noise = (float*)malloc(sizeof(float)*n*nsub*nchan);

	seed = TKsetSeed();
	for (i=0;i<n;i++)
	{
		sum = 0;
		for (j=0;j<nsub*nchan;j++)
		{
			sum += data[i*nsub*nchan+j];
			noise[i*nsub*nchan+j]=TKgaussDev(&seed);
			//printf ("%d %f %f\n", i*nsub*nchan+j, data[i*nsub*nchan+j], noise[i*nsub*nchan+j]);
		}
		flux[i] = sum/(nsub*nchan);
		x[i] = (float)(i);
		//printf ("%f %f\n", x[i], flux[i]);
	}

	// make plots

	char caption[1024];

	sprintf (caption, "%s", "Flux density variation");

	// plot 
	//cpgbeg(0,"/xs",1,1);
	cpgbeg(0,dname,1,2);
      	cpgsch(1.2); // set character height
      	cpgscf(2); // set character font
      	cpgenv(0,n,0,3,0,1); // set window and viewport and draw labeled frame
	//cpgsvp(0.1,0.5,0.5,0.9); // set viewport
	//cpgswin(0.0,10.0,0.0,2.0); // set window, must be float number
	//cpgbox("BCTSN",0.0,0,"BCTSN",0.0,0);
      	cpglab("Time","Flux (mJy)",caption);
	cpgpt(n,x,flux,9);

	////////////////////////////
	// make histogram
	xHis = (float*)malloc(sizeof(float)*step);
	val = (float*)malloc(sizeof(float)*step);
	xHisN = (float*)malloc(sizeof(float)*step);
	valN = (float*)malloc(sizeof(float)*step);
	histogram (noise, nsub*nchan*n, xHisN, valN, -5.0, 5.0, step);
	histogram (data, nsub*nchan*n, xHis, val, -5.0, 5.0, step);

	sprintf (caption, "%s", "Flux density histogram");
      	cpgenv(-5,5,0,4500,0,1); // set window and viewport and draw labeled frame
      	cpglab("Flux (mJy)","Number",caption);
	cpgbin(step,xHis,val,0);
	cpgsci(2);
	cpgbin(step,xHisN,valN,0);
	cpgend();
	///////////////////////////////////////////////////////
	
	free(flux);
	free(x);
	free(xHis);
	free(val);
	free(noise);
	free(xHisN);
	free(valN);

	return 0;
}

int histogram (float *data, int n, float *x, float *val, float low, float up, int step)
{
	int i,j,count;
	float width;
	float *temp;

	temp = (float*)malloc(sizeof(float)*(step+1));

	width = (up-low)/step;
	for (i=0; i<step; i++)
	{
		x[i] = low + i*width + width/2.0;
	}

	for (i=0; i<=step; i++)
	{
		temp[i] = low + i*width;
	}

	for (i=0; i<step; i++)
	{
		count = 0;
		for (j=0; j<n; j++)
		{
			if (data[j]>=temp[i] && data[j]<temp[i+1])
			{
				count += 1;
			}
		}
		val [i] = count;
	}

	free(temp);
	return 0;
}
