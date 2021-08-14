#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#include "utils.h"
#include "fatorLU.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

int main(int argc, char *argv[])
 {
  int i,j,n,code,count;
  double t1,t2;

  t1 = timestamp();

  if ( Init(argc,argv) != 0 ) return -1; // interpretação dos argumentos e inicialização das varíaveis

  LIKWID_MARKER_INIT;
  
  LIKWID_MARKER_START("T1");

  AlocaMemoria(&mat);
  lerValoresTabeladas(&mat);
  decompose(&mat);

  for ( count = 0; count < mat.m ; count++ )
   {
    lerFuncaoTabelada(&mat);
    solve(&mat);
    resultado(&mat,0);
    MinimosQuadrados(&mat);
    decompose(&mat);
    solve(&mat);
    resultado(&mat,1);
   }

  LIKWID_MARKER_STOP("T1");

  Finalizar(); // liberar a memória e fechar o arquivo de saída
  LIKWID_MARKER_CLOSE;

  t2 = timestamp() - t1;
  printf("\ntempo: %lf\n",t2);
  return 0;
 }


