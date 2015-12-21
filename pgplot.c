#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "dynSpecSim.h"

//#define MAX(a,b) ( (a)>(b) ? (a) : (b))
//#define MIN(a,b) ( (a)<(b) ? (a) : (b))


void heatMap (float *tab, acfStruct *acfStructure)
{
  //int i,j;                     
  //int dimx = acfStructure.ns;
	//int dimy = acfStructure.nf; // dimensions 
  //float tab[dimx*dimy];       // value
  float zmin,zmax;            /* min et max des valeurs de la fonction */
  float tr[6];                /* matrice utilisee par pgimag */

	int dimx = acfStructure->nsubint;
	int dimy = acfStructure->nchn;
	double bw = acfStructure->bw;
  zmin=0; zmax=1;

	double f1 = 1241; // MHz
	double f2 = 1497; // MHz
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
  tr[3]=f2;
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
	cpgbox("BCTSIN",10,5,"BCTSIN",50,5);
  cpglab("Subintegration","Frequency (MHz)","Freq: 1369.0 MHz BW: -256.000 Length: 3840.0 S/N: 1000.0");
	//cpgtick(0,0,1024,0,1.0/64,0.1,0.2,0,0,"1");
	//palett(3, -0.4, 0.3);
  //cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);
  cpggray(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);

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

