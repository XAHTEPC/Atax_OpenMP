#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>

/* Размерность проблемы. */
const int M = 5000;  // Примерное значение для M
const int N = 5000;  // Примерное значение для N

/* Инициализация массивов. */
void init_array(int m, int n, double A[m][n], double x[n]) {
    int i, j;
    double fn = (double)n;

    for (i = 0; i < n; i++) {
        x[i] = 1 + (i / fn);
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = (double)((i + j) % n) / (5 * m);
        }
    }
}

/* Вывод массива для проверки. */
void print_array(int n, double y[n]) {
    int i;
    for (i = 0; i < n; i++) {
        if (i % 20 == 0) printf("\n");
        printf("%lf ", y[i]);
    }
    printf("\n");
}

/* Основное вычислительное ядро. */
void kernel_atax(int m, int n, double A[m][n], double x[n], double y[n], double tmp[m]) {
    int i, j;

    /* Инициализация вектора y. */
    for (i = 0; i < n; i++) {
        y[i] = 0.0;
    }

    /* Основной цикл вычислений. */
    for (i = 0; i < m; i++) {
        tmp[i] = 0.0;
        for (j = 0; j < n; j++) {
            tmp[i] += A[i][j] * x[j];
        }
        for (j = 0; j < n; j++) {
            y[j] += A[i][j] * tmp[i];
        }
    }
}

void kernel_atax_paral(int m, int n, double A[m][n], double x[n], double y[n], double tmp[m]) {

    #pragma scop
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        
    }
   #pragma omp parallel for
    for(int i = 0; i < m; i++) {
    	tmp[i] = 0.0;
        for (int j = 0; j < n; j++) {
            tmp[i] += A[i][j] * x[j];
        }
    }

   
    for (int i = 0; i < m; i++) {   
         #pragma omp parallel for
        for (int j = 0; j < n; j++) {
            y[j] += A[i][j] * tmp[i];
        }
    }
    #pragma endscop
    
    
    /*#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
    }

    #pragma omp parallel 
    {
        double tmp_local[m]; 
        double y_local[n]; 

        for (j = 0; j < n; j++) {
            y_local[j] = 0.0;
        }

        #pragma omp for 
        for (i = 0; i < m; i++) {
            tmp_local[i] = 0.0;
            for (j = 0; j < n; j++) {
                tmp_local[i] += A[i][j] * x[j];
            }

            for (j = 0; j < n; j++) {
                y_local[j] += A[i][j] * tmp_local[i];
            }
        }

        
        #pragma omp single
        for (j = 0; j < n; j++) {
            y[j] += y_local[j];
        }
      }*/
    
}

int main() {
    /* Задание размеров. */
    int m = M;
    int n = N;
    //clock_t start, stop;

    /* Выделение памяти под массивы. */
    double (*A)[N] = (double(*)[N]) malloc(m * N * sizeof(double));
    double *x = (double*) malloc(n * sizeof(double));
    double *y = (double*) malloc(n * sizeof(double));
    double *tmp = (double*) malloc(m * sizeof(double));

    /* Инициализация массивов. */
    init_array(m, n, A, x);

    /* Запуск основного вычисления. */
    struct timeval start, stop;
    gettimeofday(&start, NULL);
    //start = clock();
    kernel_atax(m, n, A, x, y, tmp);

    /* Вывод результата. */
    //print_array(n, y);
    //stop = clock();
    gettimeofday(&stop, NULL); 
    //double cpu_time_used_seq = ((double)(stop - start)) / CLOCKS_PER_SEC;
    //printf("Time without parallel: %f \n", cpu_time_used_seq);
    double elapsed_time = (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0;
    printf("Time without parallel: %f seconds\n", elapsed_time);


    double (*A2)[N] = (double(*)[N]) malloc(m * N * sizeof(double));
    double *x2 = (double*) malloc(n * sizeof(double));
    double *y2 = (double*) malloc(n * sizeof(double));
    double *tmp2 = (double*) malloc(m * sizeof(double));

    init_array(m, n, A2, x2);
    gettimeofday(&start, NULL); 
    //start = clock();
    kernel_atax_paral(m, n, A2, x2, y2, tmp2);
    gettimeofday(&stop, NULL); 
    //stop = clock();
    //cpu_time_used_seq = ((double)(stop - start)) / CLOCKS_PER_SEC;
    //printf("Time with parallel: %f \n", cpu_time_used_seq);
    elapsed_time = (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0;
    printf("Time without parallel: %f seconds\n", elapsed_time);

    /* Освобождение памяти. */
    free(A);
    free(x);
    free(y);
    free(tmp);
    free(A2);
    free(x2);
    free(y2);
    free(tmp2);

    return 0;
}
