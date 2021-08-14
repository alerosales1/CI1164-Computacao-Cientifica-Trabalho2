
typedef struct
 {
  int n;          // número de valores tabelados
  int m;          // número de funções f(x) tabeladas
  long double *A;      // ponteiro sobre a matriz A1 ou A2
  long double *A1;     // matriz de tamanho n*n para a interpolação polinomial
  long double *A2;     // matriz de tamanho n*n para o metódo dos mínimos quadrados
  long double *LU;     // ponteiro sobre a matriz LU1 ou LU2
  long double *LU1;    // LU = L - ID + U : matrizes L e U para a interpolação polinomial
  long double *LU2;    // LU = L - ID + U : matrizes L e U para o metódo dos mínimos quadrados
  long double *a;      // matriz dos coeficientes dos polinômios gerados pela Interpolação
  long double *b;      // ponteiro sobre o vetor b1 ou b2
  long double *b1;     // vetor dos coeficientes para a interpolação polinomial
  long double *b2;     // vetor dos coeficientes para o metódo dos mínimos quadrados
  long double *y;      // Ly = b
  long double *x;      // Ux = y
  long double *xi;     // valores tabeladas
  long double *v;      // vetor dos xi^k -> ( 1 , xi , xi² , xi³ , .... , xi^(n-1) )
  long double *f;      // vetor solução
  int *P;         // vetor de permutação
  double t_LU;    // tempo para a triangularização
  double tm_y;    // tempo médio do cálculo de y em Ly = b, em milisegundos
  double tm_x;    // tempo médio do cálculo de x em Ux = y, em milisegundos
 } MATRIZ;

MATRIZ mat;
FILE *fout; // ponteiro sobre o arquivo de saída. Na ausência de um arquivo de saída, fout = stdout
int (*decompose)(MATRIZ *mat); // ponteiro sobre a função de decomposição com ou sem otimização
void (*solve)(MATRIZ *mat);    // ponteiro sobre a função de resolução por LU

// Liberaçao da memória alocada
void liberaMatriz (MATRIZ *mat);
 
// finalizar o programa
void Finalizar();

// sair do programa com erro
// str : mensagem de erro
void sair (char *str);

//alocação de memória
int AlocaMemoria(MATRIZ *mat);

// cálculo de x^n com loong double
long double power(long double x,int n);

// leitura dos dados com stdin e alocaçao de memória 
int lerValoresTabeladas(MATRIZ *mat);

// leitura de uma função tabelada
int lerFuncaoTabelada (MATRIZ *mat);

// decomposição LU sem pivotamento
// return 0 se a decomposição foi efetuada
int LUdecompose (MATRIZ *mat);

// decomposição LU com pivoteamento
int LUPdecompose (MATRIZ *mat); 

// resolve o sistema usando a decomposição LU sem otimização
void LUSolve (MATRIZ *mat);

// resolve o sistema usando a decomposição LU com otimização
void LUSolve_otimizado (MATRIZ *mat);

// cálculo a matriz dos mínimos quadrados
void MinimosQuadrados (MATRIZ *mat);

// imprimir a matriz A e a Matriz IA inversa da matriz A
void resultado (MATRIZ *mat,int code);

// incializar o programa
int Init(int argc, char *argv[]);

