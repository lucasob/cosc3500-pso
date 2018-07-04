#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <omp.h>
#include <mpi.h>

using namespace std;

// MAIN FUNCTIONS AND VARIABLES
#define number_of_particles 10000000
#define time_steps 10000000
#define global_avg (number_of_particles/4)

double particle_xvalues[number_of_particles];
double particle_yvalues[global_avg];
double particle_fitni[global_avg];
double particle_maxi[4];

double global_best_fitness[4];
double global_best_x[1];
double global_best_y[1];

// THE FUNCTION TO FIND GLOBAL MAXIMUM OF
double opt_func(double x, double y) {
  //-(x^2) - (y^2)
  return -1*pow(x-2, 2) - pow(y-2, 2) + 6;
}
// CREATE THE ARRAY OF PARTICLES
void makeXParticles(void) {
  srand(time(NULL));
  for (int i = 0; i < number_of_particles; i++) {
    double tempx = rand() % 100;
    if (rand() % 2 > 0) tempx = -tempx;
    particle_xvalues[i] = tempx;
  }
}

int main(int argc, char *argv[]) {
  // MPI NEEDED VARIABLES
  int rank, size, avg, i;
  int start, end;

  // PSO OCCURS
  MPI_Status status;
  MPI_Request srqust, rrqust;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  avg = number_of_particles/size;

  // Begins at Rank 0, and splits the particle values into SIZE chunks
  // and sends away for work to be done
  if (rank == 0) {
    for (i = 1; i < size; i++) {
      end = (i+1)*avg - 1;
      start = i*avg;

      MPI_Isend(&particle_xvalues[start], avg, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, &srqust);
      MPI_Wait(&srqust, &status);
    }
  } else {
    MPI_Irecv(particle_xvalues, avg, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD, &rrqust);
    MPI_Wait(&rrqust, &status);
  }

  // CREATE y positions for particles
  srand(time(NULL));
  for (int i = 0; i < number_of_particles; i++) {
    double tempy = rand() % 100;
    if (rand() % 2 > 0) tempy = -tempy;
    particle_yvalues[i] = i*tempy;
  }

  // MAIN PSO COMPONENT
  double z, mymax, myBX, myBY;
  int k;
  double personal_best_x[avg];
  double personal_best_y[avg];
  double personal_best_fitness[avg];

  for (int i = 0; i < avg; i++) {
    personal_best_fitness[i] = opt_func(particle_xvalues[i], particle_yvalues[i]);
  }

  while (k < time_steps) {
    srand((unsigned)(time(NULL)));
    #pragma omp parallel
    {
      #pragma omp for
      for (int i = 0; i < avg; i++) {
        // EVALUATE THE FUNCTION FOR THESE VALUES
        z = opt_func(particle_xvalues[i], particle_yvalues[i]);
        // IS Z > THIS LOCAL BEST?
        if (z > mymax) {
          mymax = z;
          myBX = particle_xvalues[i];
          myBY = particle_yvalues[i];
        }
        // IS Z MY OWN BEST?
        if (z > personal_best_fitness[i]) {
          personal_best_fitness[i] = z;
          personal_best_x[i] = particle_xvalues[i];
          personal_best_y[i] = particle_yvalues[i];
        }

        // FUNCTION COEFFICIENTS FOR UPDATING VELOCITY
        double a = (double)rand()/(double)RAND_MAX;
        double b = (double)rand()/(double)RAND_MAX;

        // UPDATING X AND Y VELOCITIES
        double xupd = particle_xvalues[i] + a*1.5*(personal_best_x[i] - particle_xvalues[i]) + b*2*(mymax - personal_best_fitness[i]);
        double yupd = particle_yvalues[i] + a*1.5*(personal_best_y[i] - particle_yvalues[i]) + b*2*(mymax - personal_best_fitness[i]);

        particle_xvalues[i] = xupd;
        particle_yvalues[i] = yupd;
      }
    }
    // WRITE THE CURRENT BEST X, Y AND FITNESS OF THIS NEIGHBOURHOOD
    //particle_progression << myBX << "," << myBY << "," << mymax << endl;
    // NEXT TIME STEP
    k++;
  }

  // RECEIVE FINAL VALUES
  if (rank == 0) {
    global_best_fitness[0] = mymax;
    global_best_x[0] = myBX;
    global_best_y[0] = myBY;

    for (int i = 1; i < size; i++) {
      MPI_Irecv(&mymax, 1, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, &rrqust);
      global_best_fitness[i] = mymax;
      MPI_Wait(&rrqust, &status);
    }
  } else {
    MPI_Isend(&mymax, 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD, &srqust);
    MPI_Wait(&srqust, &status);
  }

  if (rank == 0) {
    double avg_best = global_best_fitness[0];
    cout << "=================" << endl;
    printf("Global max at: x = %f, y = %f from rank 0.\n", global_best_x[0], global_best_y[0]);
    for (int i = 1; i < size; i++) {
      printf("Rank %d max value: %f\n", i, global_best_fitness[i]);
      avg_best += global_best_fitness[i];
    }
    printf("Average best position across each of the ranks: %f", (avg_best)/size);
  }

  MPI_Finalize();

  // LORD HAVE MERCY ON MY POOR SOUL
  return 0;
}
