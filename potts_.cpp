/******************************************************************************

  ANOTHER CELLULAR POTTS MODEL:

    This program simulates a number of cells on a quasi-2D grid.

******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
// HEAD: INCLUDES AND DEFINES

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  Note that if you IMPORT something, you best be
//  sure that N, numCells, and OFFSETS are correct
//  as defined here:

#define  IMPORT           0	// 0 or 1
#define  COLLAGEN_OFFSET  1	// 0 or 1
#define  PERIMETER_OFFSET 1	// 0 or 1

///////////////////////////////////////////////////////////////////////////////
// HEAD: GLOBAL PARAMETERS 

  // Random Number Generation

  int seed = 			0;

  // Simulation Length ( = Loops * Flips ) 

  const int numLoops =		5000;
  const int numFlips =		pow(2,5);
  const int numPrint = 		5;

  // Energies 

  const double beta =		1.0;

  const double J_air =		-0.75;
  const double J_cel =		2.0;

  const double J_col =		0.0;

  const double L_vol = 		0.04;
  const double L_per =		0.0;
  const double L_blb =		0.0;

  // System 

  const int N =			60;

  const int numCells =		1;
  const int numCollagen =	0;

  const double cellSpawn = 	8.0;
  const double cellRadius =	8.0;

  const int collagenWidth =	1;

///////////////////////////////////////////////////////////////////////////////
// HEAD: GLOBAL VARIABLES

  int lattice[N][N][2] = {0};
  int cellVolume[numCells+1];
  int cellPerimeter[numCells+1];
  int cellVolumeList[numCells+1][9000][2];
  int cellPerimeterList[numCells+1][9000][2];

  double totalEnergy;
  double avg[3],dev[3];

  const double targetVolume = 3.141593*cellRadius*cellRadius;

///////////////////////////////////////////////////////////////////////////////
// HEAD: FUNCTIONS

  #include "potts_print_.h"
  #include "potts_spawn_.h"
  #include "potts_energy_.h"
  #include "potts_flip_.h"
  #include "potts_analysis_.h"

///////////////////////////////////////////////////////////////////////////////
// MAIN 

int main(int argc, char *argv[])
{

  char   dname[100];
  char   fname[100];

  /* Seed random number generator */

  if(seed==0)
    seed=time(0);
  srand(seed);

  /* Create output directory */

  if(argv[1]==NULL)
    strcpy(dname,"output");
  else
    strcpy(dname,argv[1]);
  strcpy(fname,"mkdir ");
  strcat(fname,dname);
  strcat(fname," 2>/dev/null");
  system(fname);

  /* Print log file */

  strcpy(fname,dname);
  strcat(fname,"/aaa_log_.dat");
  printLog(fname);

  /* Let there be life */

  #if IMPORT
  printf("\n  Warning: importing cancer.\n\n");
  readCells();
  readCollagen();
  #else
  printf("\n  Warning: creating cancer.\n\n");
  putCells();
  putCollagen();
  #endif

  strcpy(fname,dname);
  strcat(fname,"/lattice_0_.dat");
  printLattice(fname);
  strcpy(fname,dname);
  strcat(fname,"/collagen_.dat");
  printCollagen(fname);

  /* Calculate the initial energy */

  measureCells();

  totalEnergy=Hamiltonian();

  strcpy(fname,dname);
  strcat(fname,"/aaa_energy_.dat");
  FILE* efile=fopen(fname,"w");
  fprintf(efile,"%5d ",0);
  fprintf(efile,"%8.3lf %8.3lf ",totalEnergy,0.0);
  fprintf(efile,"%8.3lf %8.3lf ",avg[0],dev[0]);
  fprintf(efile,"%8.3lf %8.3lf ",avg[1],dev[1]);
  fprintf(efile,"%8.3lf %8.3lf ",avg[2],dev[2]);
  fprintf(efile,"\n");
  fflush(efile);

  /* Perform flips and conditionally accept the change via the Metropolis algorithm */

  int accepted;

  printf("  Mutating: %d x 2^%d spin flips...\n\n",numLoops,(int)log2((double)numFlips));

  for(int outerCount=1; outerCount<=numLoops; outerCount++){

    if(outerCount%15==0){
      printf("\r    count = %d",outerCount);
      fflush(stdout);
    }

    accepted=0;
    for(int count=0; count<numFlips; count++)
      accepted+=flip();

    measureCells();

    fprintf(efile,"%5d ",outerCount);
    fprintf(efile,"%8.3lf %8.3lf ",totalEnergy,(double)accepted/(double)numFlips);
    fprintf(efile,"%8.3lf %8.3lf ",avg[0],dev[0]);
    fprintf(efile,"%8.3lf %8.3lf ",avg[1],dev[1]);
    fprintf(efile,"%8.3lf %8.3lf ",avg[2],dev[2]);
    fprintf(efile,"\n");
    fflush(efile);

    if(outerCount%numPrint==0){
      sprintf(fname,"%s/lattice_%d_.dat",dname,outerCount/numPrint);
      printLattice(fname);
    }

  }

  printf("\r    finished!              \n\n");

  fclose(efile);

  /*** Show some stats ***/

  printf("  Final statistics...\n\n");
  printf("    acceptance ratio : %8.3lf\n\n",(double)accepted/(double)(numFlips));
  printf("    cell volume      : %8.3lf +/- %7.3lf\n",avg[0],dev[0]);
  printf("    cell perimeter   : %8.3lf +/- %7.3lf\n",avg[1],dev[1]);
  printf("    cell anisotropy  : %8.3lf +/- %7.3lf\n\n",avg[2],dev[2]);
  printf("  Done.\n\n");

  return 0;

}

















































