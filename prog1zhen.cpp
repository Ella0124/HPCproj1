#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


#define MASTER 0


double randu() {
  return (double)rand() / RAND_MAX;
}

int dboard(int N) {
  int p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  int n = N / p;

  int m = 0;
  double sq2 = sqrt(2);


  for (int i = 0; i < n; i++) {
    double r = sqrt(randu());
    double theta = randu() * 2 * M_PI;

    double x = r * cos(theta);
    double y = r * sin(theta);

    if (fabs(x) * sq2 <= 1 && fabs(y) * sq2 <= 1) {
      m++;
    }
  }

  return m;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int N, R;
  if (rank == MASTER) {
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &R);
  }
  MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&R, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

  srand(rank);

  double pi, pi_sum = 0;

  double start, end;
  double total_time;
  double max_total_time_sum = 0;

  for (int i = 0; i < R; i++) {

    int m;
    int M = 0;

    start = MPI_Wtime();

    // timed //
    m = dboard(N);
    MPI_Reduce(&m, &M, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    // timed - end //

    end = MPI_Wtime();

    total_time = end - start;
    // printf("%lf\n", total_time);

    double max_total_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);


    if (rank == MASTER) {
      // printf("Max total time: %lf\n", max_total_time);

      max_total_time_sum += max_total_time;

      pi = 2 * (double)N / M;
      // printf("%lf\n", pi);
      pi_sum += pi;
    }
  }

  double pi_avg;
  if (rank == MASTER) {
    pi_avg = pi_sum / R;

    printf("N=%d, R=%d, P=%d, PI=%lf\n", N, R, p, pi_avg);
    printf("Time=%lf\n", max_total_time_sum);
  }


  MPI_Finalize();
  return 0;
}

