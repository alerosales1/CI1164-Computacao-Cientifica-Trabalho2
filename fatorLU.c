#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include "utils.h"
#include "fatorLU.h"
//#include <likwid.h>

// Liberaçao da memória alocada
void liberaMatriz (MATRIZ *mat)
 {
  int i,n;

  n = mat->n;

  if ( mat->A1 != NULL ) free(mat->A1);
  if ( mat->A2 != NULL ) free(mat->A2);

  if ( mat->LU1 != NULL ) free(mat->LU1);
  if ( mat->LU2 != NULL ) free(mat->LU2);

  if ( mat->b1 != NULL ) free(mat->b1);
  if ( mat->b2 != NULL ) free(mat->b2);

  if ( mat->y != NULL )  free(mat->y);
  if ( mat->x != NULL )  free(mat->x);
  //if ( mat->xi != NULL )  free(mat->xi);
  //if ( mat->P != NULL )  free(mat->P);
  //if ( mat->v != NULL )  free(mat->v);
  //if ( mat->a != NULL )  free(mat->a);
  //if ( mat->f != NULL )  free(mat->f);
 }

// finalizar o programa
void Finalizar()
 {
  liberaMatriz(&mat); // liberação da memória alocada

  fflush(fout); // esvaziar o buffer
  if (fout != stdout) fclose(fout); // fechar o arquivo de saída
 }

// sair do programa com erro
// str : mensagem de erro
void sair (char *str)
 {
  printf("%s\n",str); // affichar a mensagem de erro
  Finalizar();        // liberar a memória e fechar o arquivo de saída
  exit (-1);
 }

long double power(long double x,int n)
 {
  long double buf;
  int i;

  buf = 1.0;
  for ( i = 0 ; i < n ; i++ ) buf *= x;
  return buf;
 }

int lerFuncaoTabelada (MATRIZ *mat)
 {
  int i,n = mat->n;
  long double buf;

  for ( i = 0 ; i < n ; i++ )
   {
    if ( scanf("%Lf",&buf) != 1 ) sair("Failed to read double"); // leitura de um valor da função tabelada: Y(Xi)
    //buf *= 1000;
    mat->f[i] = buf;
    mat->b1[i] = buf;
   }

  mat->A  = mat->A1;
  mat->LU = mat->LU1;
  mat->b  = mat->b1;
 }

int AlocaMemoria(MATRIZ *mat)
 {
  int test,i,j,m,n;
  long double buf;

  n = m = 0;

  if ( scanf("%d",&n) != 1 ) sair("Failed to read integer"); // leitura do número de valores tabelados

  if ( n == 0 ) return 0; // finalizar o programa
  if ( n < 2 ) sair("quantia de pontos incorreta"); 

  if ( n != mat->n ) // se o número lido e anterior ao número anterior
   {
    mat->n = n;     // novo valor de n
    test = 1;       // realocar a memória
   }
  else test = 0; // n e m não mudaram, não há necessidade de realocar a memória

  if ( scanf("%d",&m) != 1 ) sair("Failed to read integer"); // leitura do número de funções f(x) tabeladas

  if ( m < 1 )  sair("quantia de funções incorreta"); 
  mat->m = m;     // novo valor de m

  // aloca mémoria para a matriz A
  if ( mat->A1 == NULL ) mat->A1 = (long double *) calloc (n*n,sizeof(long double));
  else if ( test ) mat->A1 = (long double *) realloc (mat->A1,n*n*sizeof(long double));

  if ( mat->A1 == NULL ) sair("erro de alocação memória");

  // aloca mémoria para a matriz A
  if ( mat->A2 == NULL ) mat->A2 = (long double *) calloc (n*n,sizeof(long double));
  else if ( test ) mat->A2 = (long double *) realloc (mat->A2,n*n*sizeof(long double));

  if ( mat->A2 == NULL ) sair("erro de alocação memória");

  // aloca mémoria para a matriz LU
  if ( mat->LU1 == NULL ) mat->LU1 = (long double *) calloc (n*n,sizeof(long double));
  else if ( test ) mat->LU1 = (long double *) realloc (mat->LU1,n*n*sizeof(long double));

  if ( mat->LU1 == NULL ) sair("erro de alocação memória");

  // aloca mémoria para a matriz LU
  if ( mat->LU2 == NULL ) mat->LU2 = (long double *) calloc (n*n,sizeof(long double));
  else if ( test ) mat->LU2 = (long double *) realloc (mat->LU2,n*n*sizeof(long double));

  if ( mat->LU2 == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor dos valores tabelados
  if ( mat->x == NULL ) mat->x = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->x = (long double *) realloc (mat->x,n*sizeof(long double));

  if ( mat->x == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->y == NULL ) mat->y = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->y = (long double *) realloc (mat->y,n*sizeof(long double));

  if ( mat->y == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->b1 == NULL ) mat->b1 = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->b1 = (long double *) realloc (mat->b1,n*sizeof(long double));

  if ( mat->b1 == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->b2 == NULL ) mat->b2 = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->b2 = (long double *) realloc (mat->b2,n*sizeof(long double));

  if ( mat->b2 == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->v == NULL ) mat->v = (long double *) calloc (2*n-1,sizeof(long double));
  else if ( test ) mat->v = (long double *) realloc (mat->v,(2*n-1)*sizeof(long double));

  if ( mat->v == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->a == NULL ) mat->a = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->a = (long double *) realloc (mat->a,n*sizeof(long double));

  if ( mat->a == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->f == NULL ) mat->f = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->f = (long double *) realloc (mat->f,n*sizeof(long double));

  if ( mat->f == NULL ) sair("erro de alocação memória");

  // aloca mémoria para o vetor das funcoes tabeladas
  if ( mat->xi == NULL ) mat->xi = (long double *) calloc (n,sizeof(long double));
  else if ( test ) mat->xi = (long double *) realloc (mat->xi,n*sizeof(long double));

  if ( mat->xi == NULL ) sair("erro de alocação memória");

  if ( mat->P == NULL ) mat->P = (int *) calloc (n+1,sizeof(int));
  else if ( test ) mat->P = (int *) realloc(mat->P,(n+1)*sizeof(int));

  if ( mat->P == NULL ) sair("erro de alocação de memória - vetor permutação");
 }

int lerValoresTabeladas(MATRIZ *mat)
 {
  int i,j,n;
  long double buf;

  n = mat->n;

  // preenchimento das linhas da matriz
  for ( i = 0 ; i < n ; i++ ) 
   {
    if ( scanf("%Lf",&buf) != 1 ) sair("Failed to read double"); // leitura de um valor tabelado: Xi
    //buf *= 1000;
    mat->xi[i] = buf;  // salvar o valor i da tabela

    mat->A1[i*n+n-1] = 1.0; // mat->A[i][n-1] = 1.0;

    // calcular os outros coeficientes da matriz: A[i][j] = Xi^(n-j-1)
    for ( j = n-2 ; j >= 0 ; j-- ) mat->A1[i*n+j] = buf*mat->A1[i*n+j+1]; // mat->A[i][j] = buf*A[i][j+1];
 
    mat->A  = mat->A1;
    mat->b  = mat->b1;
    mat->LU = mat->LU1;
   } 

  return 1;      
 }

void MinimosQuadrados (MATRIZ *mat)
 {
  int i,j,k,n;
  long double buf;

  n = mat->n;

  for ( k = 0 ; k <= 2*n-2 ; k++ ) 
   {
    buf = 0.0;
    for ( i = 0; i < n ; i++ ) buf += power(mat->xi[i],k);
    mat->v[k] = buf;
   }
 
  for ( k = 0 ; k < n ; k++ )
    for ( i = 0; i <= k; i++ ) mat->A2[i*n+k-i] = mat->v[k];

  for ( k = 1; k < n ; k++ )
    for ( i = k; i < n ; i++ )
     {
      j = k+n-1-i;
      mat->A2[i*n+j] = mat->v[i+j];
     }
  
  for ( k = 0 ; k < n ; k++ )
   {
    buf = 0.0;
    for ( i = 0; i < n ; i++ ) buf += mat->f[i]*power(mat->xi[i],k);
    mat->a[k] = buf;
   }

  for ( k = 0 ; k < n ; k++ ) mat->b2[k] = mat->a[k];

  mat->A  = mat->A2;
  mat->b  = mat->b2;
  mat->LU = mat->LU2;
 }

int LUPdecompose (MATRIZ *mat) // decomposição LU com pivoteamento
 {
  double tempo = timestamp();
  int i,j,k,n,lpivo,troca;
  long double pivo,buf,m,det;

  fprintf(fout,"com pivoteamento\n");
  n = mat->n;

  for ( i = 0 ; i < n ; i++ )  // para não modificar a matriz A    
    for ( j = 0 ; j < n ; j++ ) mat->LU[(i+1)*n+j+1] = mat->A[i*n+j];//mat->LU[i+1][j+1] = mat->A[i][j]; // indices de 1..n

  for ( i = 1 ; i <= n ; i++ ) mat->P[i] = i;

  for ( k = 1 ; k < n ; k++ )
   {
    pivo = mat->LU[k*n+k]; // mat->LU[k][k];
    lpivo = k;

    // determinando o pivo
    for ( i = k+1 ; i <= n ; i++ )
      if ( fabs(mat->LU[i*n+k]) > fabs(pivo) ) // if ( fabs(mat->LU[i][k]) > fabs(pivo) )
       {
        pivo = mat->LU[i*n+k]; // mat->LU[i][k];
        lpivo = i;
       }

    if ( pivo == 0 ) return -1; // matriz singular

    if ( lpivo != k )
     {
      //pivotando P
      troca = mat->P[k];
      mat->P[k] = lpivo;
      lpivo = troca;

      //pivotando as linhas da matriz
      for ( j = 1 ; j <= n ; j++ )
       {
        buf = mat->LU[k*n+j]; // mat->LU[k][j];
        mat->LU[k*n+j] = mat->LU[lpivo*n+j]; // mat->LU[k][j] = mat->LU[lpivo][j];
        mat->LU[lpivo*n+j] = buf; // mat->LU[lpivo][j] = buf;
       }
     }

    for ( i = k+1 ; i <= n ; i++ )
     {
      if ( mat->LU[k*n+k] == 0 ) return -1; // divisão por zero // if ( mat->LU[k][k] == 0 )
      m = mat->LU[i*n+k] / mat->LU[k*n+k]; //       m = mat->LU[i][k] / mat->LU[k][k];
      mat->LU[i*n+k] = m; // mat->LU[i][k] = m;

      for ( j = k+1 ; j <= n ; j++ ) mat->LU[i*n+j] -= m * mat->LU[k*n+j]; // mat->LU[i][j] -= m * mat->LU[k][j];
     }
   }

  for ( i = 0 ; i < n ; i++ )
    for ( j = 0 ; j < n ; j++ ) mat->LU[i*n+j] = mat->LU[(i+1)*n+j+1]; // restaurar os indices de 0 .. n-1
                                                                   // mat->LU[i][j] = mat->LU[i+1][j+1];
  mat->t_LU = timestamp() - tempo; // tempo para a triangularização

  // det|A| = det|L.U| = det|L|.det|U| = 1 x det|U| = det|U|
  // O determinante de uma matriz triangular é o produto dos elementos da diagonal principal.
  // É impossível inverter a matriz se o determinante = 0
  det = 1.0;
  for ( i = 0 ; i < n ; i++ ) det *= mat->LU[i*n+i]; // mat->LU[i][i];
  if ( det == 0.0 ) return -2;

  return 0;
 }

int LUdecompose (MATRIZ *mat) // decomposição LU sem otimização
 {
  double tempo = timestamp();
  int i,j,k,n;
  long double det,buf;
  //printf("sem otimização\n");

  n = mat->n;

  for ( i = 0 ; i < n ; i++ )  // para não modificar a matriz A    
    for ( j = 0 ; j < n ; j++ ) mat->LU[i*n+j] = mat->A[i*n+j];//mat->LU[i][j] = mat->A[i][j]; 

  for ( k = 0 ; k < n ; k++ )
   {   // elementos da matriz U                         
    for ( j = k ; j < n ; j++ ) 
     {
      buf = 0.0;

      for ( i = 0 ; i < k ; i++ )  buf += mat->LU[k*n+i]*mat->LU[i*n+j];  //buf += mat->LU[k][i]*mat->LU[i][j];  

      mat->LU[k*n+j] -= buf;//mat->LU[k][j] -= buf;
     }

        // elementos da matriz L
    for ( i = k+1 ; i < n ; i++ ) 
     {
      buf = 0.0;

      for ( j = 0 ; j < k ; j++ )  buf += mat->LU[i*n+j]*mat->LU[j*n+k];  //buf += mat->LU[i][j]*mat->LU[j][k];  

      if ( mat->LU[k*n+k] == 0.0 ) return -1;//if ( mat->LU[k][k] == 0.0 ) return -1;
      mat->LU[i*n+k] -= buf;//mat->LU[i][k] -= buf;
      mat->LU[i*n+k] /= mat->LU[k*n+k];//mat->LU[i][k] /= mat->LU[k][k];
      //printf("LU[%d][%d] = %f\n",k,k,(float) mat->LU[k*n+k]);
     }
   }

  mat->t_LU = timestamp() - tempo; // tempo para a triangularização

  // det|A| = det|L.U| = det|L|.det|U| = 1 x det|U| = det|U|
  // O determinante de uma matriz triangular é o produto dos elementos da diagonal principal.
  // É impossível inverter a matriz se o determinante = 0
  det = 1.0;
  for ( i = 0 ; i < n ; i++ ) det *= mat->LU[i*n+i];//det *= mat->LU[i][i];
  if ( det == 0.0 ) return -2;

  return 0;
 }

int LUdecompose_otimizado (MATRIZ *mat) 
 {
  double tempo = timestamp();
  int i,j,k,n,end;
  long double det,buf,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10;
  //fprintf(fout,"em cache\n");
  n = mat->n;

  for ( k = 0 ; k < n ; k++ )
   {                             
    for ( j = k ; j < n ; j++ ) // elementos da matriz U
     {
      end = k - (k % 10);   
      sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = 0.0;     
      buf = mat->A[k*n+j]; // mat->A[k][j];

      for ( i = 0 ; i < end ; i += 10 )  
       {
        sum1  += mat->LU[k*n+i]*mat->LU[i*n+j];     // mat->LU[k][i]*mat->LU[i][j];
        sum2  += mat->LU[k*n+i+1]*mat->LU[(i+1)*n+j];  
        sum3  += mat->LU[k*n+i+2]*mat->LU[(i+2)*n+j];  
        sum4  += mat->LU[k*n+i+3]*mat->LU[(i+3)*n+j];  
        sum5  += mat->LU[k*n+i+4]*mat->LU[(i+4)*n+j];
        sum6  += mat->LU[k*n+i+5]*mat->LU[(i+5)*n+j];  
        sum7  += mat->LU[k*n+i+6]*mat->LU[(i+6)*n+j];  
        sum8  += mat->LU[k*n+i+7]*mat->LU[(i+7)*n+j];
        sum9  += mat->LU[k*n+i+8]*mat->LU[(i+8)*n+j];  
        sum10 += mat->LU[k*n+i+9]*mat->LU[(i+9)*n+j];    
       }

      buf -= ( sum10 + sum9 + sum8 + sum7 + sum6 + sum5 + sum4 + sum3 + sum2 + sum1 );

      for ( i = end ; i < k ; i++ ) buf -= mat->LU[k*n+i]*mat->LU[i*n+j];  // mat->LU[k][i]*mat->LU[i][j];

      mat->LU[k*n+j] = buf; // mat->LU[k][j] = buf;
     }

    //printf("\nelementos matriz L\n");
    for ( i = k+1 ; i < n ; i++ ) // elementos da matriz L
     {
      end = k - (k % 10);   
      sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = 0.0;     
      buf = mat->A[i*n+k]; // mat->A[i][k];

      for ( j = 0 ; j < end ; j += 10 )
       {
        sum1  += mat->LU[i*n+j]*mat->LU[j*n+k];      // mat->LU[i][j]*mat->LU[j][k];
        sum2  += mat->LU[i*n+j+1]*mat->LU[(j+1)*n+k];  
        sum3  += mat->LU[i*n+j+2]*mat->LU[(j+2)*n+k];  
        sum4  += mat->LU[i*n+j+3]*mat->LU[(j+3)*n+k];  
        sum5  += mat->LU[i*n+j+4]*mat->LU[(j+4)*n+j];
        sum6  += mat->LU[i*n+j+5]*mat->LU[(j+5)*n+k];  
        sum7  += mat->LU[i*n+j+6]*mat->LU[(j+6)*n+k];  
        sum8  += mat->LU[i*n+j+7]*mat->LU[(j+7)*n+k];  
        sum9  += mat->LU[i*n+j+8]*mat->LU[(j+8)*n+k];  
        sum10 += mat->LU[i*n+j+9]*mat->LU[(j+9)*n+k];    
       }

      buf -= ( sum10 + sum9 + sum8 + sum7 + sum6 + sum5 + sum4 + sum3 + sum2 + sum1 );

      for ( j = end ; j < k ; j++ ) buf -= mat->LU[i*n+j]*mat->LU[j*n+k];  // mat->LU[i][j]*mat->LU[j][k];

      if ( mat->LU[k*n+k] == 0.0 ) return -1; // if ( mat->LU[k][k] == 0.0 ) return -1;
      mat->LU[i*n+k] = buf/mat->LU[k*n+k]; // mat->LU[i][k] = buf/mat->LU[k][k];
     }
   }

  mat->t_LU = timestamp() - tempo; // tempo para a triangularização

  // det|A| = det|L.U| = det|L|.det|U| = 1 x det|U| = det|U|
  // O determinante de uma matriz triangular é o produto dos elementos da diagonal principal.
  // É impossível inverter a matriz se o determinante = 0
  det = 1.0;
  for ( i = 0 ; i < n ; i++ ) det *= mat->LU[i*n+i]; // mat->LU[i][i];
  if ( det == 0.0 ) return -2;

  return 0;
 }

// inversão da matriz A pelo metodo da fatoração LU
void LUSolve (MATRIZ *mat)
 {
  int i,j,k,n;
  double t1,t2,t3;
  long double max,norma;

  n = mat->n;
  mat->tm_y = 0.0;
  mat->tm_x = 0.0;

  max = 0.0;

  t1 = timestamp();

  // resolver Ly = b
  for ( i = 0; i < n; i++) 
   {
    mat->y[i] = mat->b[i]; 
    for ( k = 0; k < i; k++) mat->y[i] -= mat->y[k] * mat->LU[i*n+k];//mat->LU[i][k];
   }

  t2 = timestamp();
  mat->tm_y += t2 - t1;

  // resolver Ux = y
  for ( i = n - 1; i >= 0; i-- )  
   {
    mat->x[i] = mat->y[i];
    for ( k = i + 1; k < n; k++ ) mat->x[i] -= mat->x[k] * mat->LU[i*n+k];//mat->LU[i][k];

    mat->x[i] /= mat->LU[i*n+i];//mat->LU[i][i];
   }

  t3 = timestamp();
  mat->tm_x += t3 - t2;
 }

// inversão da matriz A pelo metodo da fatoração LU
void LUSolve_otimizado (MATRIZ *mat)
 {
  int i,j,k,n,end;
  double t1,t2,t3;
  long double max,norma;
  long double buf,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10;

  n = mat->n;
  mat->tm_y = 0.0;
  mat->tm_x = 0.0;

  max = 0.0;

  t1 = timestamp();

  // resolver Ly = b
  for ( i = 0; i < n; i++) 
   {
    sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = 0.0;     
    mat->y[i] = mat->b[i];
    end = i - (i % 10);   
    buf = mat->y[i];

    for ( k = 0; k < end; k+= 10 )
     {
      //mat->y[i] -= mat->y[k] * mat->LU[i*n+k];//mat->LU[i][k];
      sum1  += mat->y[k]*mat->LU[i*n+k];     
      sum2  += mat->y[k+1]*mat->LU[i*n+k+1];  
      sum3  += mat->y[k+2]*mat->LU[i*n+k+2];  
      sum4  += mat->y[k+3]*mat->LU[i*n+k+3];  
      sum5  += mat->y[k+4]*mat->LU[i*n+k+4];
      sum6  += mat->y[k+5]*mat->LU[i*n+k+5];  
      sum7  += mat->y[k+6]*mat->LU[i*n+k+6];  
      sum8  += mat->y[k+7]*mat->LU[i*n+k+7];  
      sum9  += mat->y[k+8]*mat->LU[i*n+k+8];  
      sum10 += mat->y[k+9]*mat->LU[i*n+k+9];    
     }

    buf -= ( sum10 + sum9 + sum8 + sum7 + sum6 + sum5 + sum4 + sum3 + sum2 + sum1 );
    for ( k = end; k < i; k++) buf -= mat->y[k] * mat->LU[i*n+k];//mat->LU[i][k];
    mat->y[i] = buf;
   }

  t2 = timestamp();
  mat->tm_y += t2 - t1;

  // resolver Ux = y
  for ( i = n - 1; i >= 0; i-- )  
   {
    sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = 0.0;     
    mat->x[i] = mat->y[i];
    //end = n-i - ((n-i) % 10);   
    buf = mat->x[i];
    end = i+1;
/*
    for ( k = i + 1; k < end; k += 10 )
     {
      //mat->x[i] -= mat->x[k] * mat->LU[i*n+k];//mat->LU[i][k];
      sum1  += mat->x[k]*mat->LU[i*n+k];     
      sum2  += mat->x[k+1]*mat->LU[i*n+k+1];  
      sum3  += mat->x[k+2]*mat->LU[i*n+k+2];  
      sum4  += mat->x[k+3]*mat->LU[i*n+k+3];  
      sum5  += mat->x[k+4]*mat->LU[i*n+k+4];
      sum6  += mat->x[k+5]*mat->LU[i*n+k+5];  
      sum7  += mat->x[k+6]*mat->LU[i*n+k+6];  
      sum8  += mat->x[k+7]*mat->LU[i*n+k+7];  
      sum9  += mat->x[k+8]*mat->LU[i*n+k+8];  
      sum10 += mat->x[k+9]*mat->LU[i*n+k+9];    
     }

    buf -= ( sum10 + sum9 + sum8 + sum7 + sum6 + sum5 + sum4 + sum3 + sum2 + sum1 );
*/
    for ( k = end; k < n; k++ ) buf -= mat->x[k] * mat->LU[i*n+k];//mat->LU[i][k];
    mat->x[i] = buf;

    mat->x[i] /= mat->LU[i*n+i];//mat->LU[i][i];
   }

  t3 = timestamp();
  mat->tm_x += t3 - t2;
 }

void resultado(MATRIZ *mat,int code)
 {
  int i,j,n,m;

  n = mat->n;

  if ( code == 0 )
   {
    for ( j = n-1 ; j >=0 ; j-- ) printf("%.19f ",(float) mat->x[j]);
    mat->A  = mat->A2;
    mat->b  = mat->b2;
    mat->LU = mat->LU2;
   }
  else
   {
    for ( j = 0 ; j < n ; j++ ) printf("%.19f ",(float) mat->x[j]);
    mat->A  = mat->A1;
    mat->b  = mat->b1;
    mat->LU = mat->LU1;
   }

  printf("\n");
 }

int Init(int argc, char *argv[])
 {
  int i;
  fout = NULL;
  decompose = LUdecompose; // decomposição LU sem otimização
  solve = LUSolve;
  mat.n = 0; // inicialização da dimensão da matriz

  if ( argc > 1 )
    for ( i = 1 ; i < argc ; i++ )
     {
      if ( argv[i][0] == '-' )
       {
        switch ( argv[i][1] ) 
         {
          case 'o' : decompose = LUdecompose_otimizado; // decomposição LU otimizada
                     solve = LUSolve_otimizado; break; 
          case 'f' : if ( i+1 < argc) { i++; fout = fopen(argv[i],"w"); } break; // abrir o aquivo de saída
          default  : // mensagem de erro
                     printf("argumento desconhecido: %s\n%s -p -o <arquivo_saida>\n",argv[i],argv[0]);
                     printf("-f <arquivo_saida> : arquivo de saída opcional\n");
                     printf("-o : otimização\n"); // decomposição LU com otimização
                     return -1;
         }
       }
      else
       {      // mensagem de erro
        printf("argumento desconhecido: %s\n%s -p -o <arquivo_saida>\n",argv[i],argv[0]);
        printf("-f <arquivo_saida> : arquivo de saída opcional\n");
        printf("-o : otimização\n"); // decomposição LU com otimização
        return -1;
       }
     }
 
  if (fout == NULL) fout = stdout; // não foi definido arquivo de saída

  return 0;
 }

