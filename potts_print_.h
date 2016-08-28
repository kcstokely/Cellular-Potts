/*******************************************************************************/
/*** PRINTING FUNCTIONS ***/

  void  printLog(char*);
  void  printCells(char*);
  void  printLattice(char*);
  void  printCollagen(char*);

/*******************************************************************************/
/*** Prints a log file ***/

void printLog(char *fname)
{
  FILE *pfile;
  pfile=fopen(fname,"w");
  fprintf(pfile,"SEED\t\t%d\n\n",seed);
  fprintf(pfile,"LOOPS\t\t%d\n",numLoops);
  fprintf(pfile,"FLIPS\t\t%d\n",numFlips);
  fprintf(pfile,"PRINT\t\t%d\n\n",numPrint);
  fprintf(pfile,"BETA\t\t%lf\n\n",beta);
  fprintf(pfile,"AIR\t\t%lf\n",J_air);
  fprintf(pfile,"CELL\t\t%lf\n\n",J_cel);
  fprintf(pfile,"COLLAGEN\t%lf\n\n",J_col);
  fprintf(pfile,"VOLUME\t\t%lf\n",L_vol);
  fprintf(pfile,"PERIMETER\t%lf\n",L_per);
  fprintf(pfile,"BLOBULAR\t%lf\n\n",L_blb);
  fprintf(pfile,"LATTICE EDGE\t%d\n",N);
  fprintf(pfile,"NUM CELLS\t%d\n",numCells);
  fprintf(pfile,"NUM COLLAGEN\t%d\n",numCollagen);
  fprintf(pfile,"CELL SPAWN\t%lf\n",cellSpawn);
  fprintf(pfile,"CELL RADIUS\t%lf\n",cellRadius);
  fprintf(pfile,"COLLAGEN WIDTH\t%d\n\n",collagenWidth);
  fclose(pfile);
}

/*******************************************************************************/
/*** Prints the cells to a text file ***/

void printCells(char *fname)
{
  FILE *printTo;
  printTo=fopen(fname,"w");
  for(int cell=1;cell<=numCells;cell++){
    fprintf(printTo,"%d\n",cellVolume[cell]);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
        if(lattice[i][j][0]==cell)
          fprintf(printTo,"%9d%9d%9d\n",i,j,2);
      }
    }
  }
  fclose(printTo);
}

/*******************************************************************************/
/*** Prints the lattice to a text file ***/

void printLattice(char *fname)
{
  FILE *pfile;
  pfile=fopen(fname,"w");
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      //air
      if(lattice[i][j][0]==0){
        if(lattice[i][j][1]==0)
          fprintf(pfile,"%d ",0);
        else
          fprintf(pfile,"%d ",COLLAGEN_OFFSET*numCells+COLLAGEN_OFFSET);
      }
      //cell
      else{
        for(int p=0;p<cellPerimeter[lattice[i][j][0]];p++){
          if(cellPerimeterList[lattice[i][j][0]][p][0]==i && cellPerimeterList[lattice[i][j][0]][p][1]==j){
            fprintf(pfile,"%d ",lattice[i][j][0]+PERIMETER_OFFSET*COLLAGEN_OFFSET+PERIMETER_OFFSET*numCells);
            goto yep;
          }
        }
        fprintf(pfile,"%d ",lattice[i][j][0]);
        yep:;
      }
    }
    fprintf(pfile,"\n");
  }
  fclose(pfile);
}

/*******************************************************************************/
/*** Prints the lower lattice only to a text file ***/

void printCollagen(char *fname)
{
  FILE *pfile;
  pfile=fopen(fname,"w");
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++)
      fprintf(pfile,"%d ",lattice[i][j][1]);
    fprintf(pfile,"\n");
  }
  fclose(pfile);
}

/*******************************************************************************/




























