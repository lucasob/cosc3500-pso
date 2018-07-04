#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_PARTICLES 1000

//function degree (highest degree of the function)
int function_degree = 4;
// provide a function to actually one-off input a function to optimise
double optimisable(int);

// result variables
double global_best_position;
double global_best_fitness;

// number of iterations
int time_steps = 1000;

// particles: each row is a particle
// each particle has current best health, current velocity, current x and y pos
// BEST POSITION[0]; CURRENT VELOCITY[1]; CURRENT X[2]; BEST FITNESS[3];
double particles[4][NUM_PARTICLES];

void makeParticles(void) {
  srand(time(NULL));

  // assign random x and y starting positions
  for (int i = 0; i < NUM_PARTICLES; i++) {
    double temp = rand()%10;
    if (rand() % 2 > 0) particles[2][i] = temp;
    else particles[2][i] = -temp;

    // for the first run the best fitness and location will be the initial point
    particles[0][i] = particles[2][i];
    particles[3][i] = optimisable(i);
  }
}

void PSO(void) {
  // trying to find global maximum
  // stop case will be $time_steps increments of "time"
  // OPTIMISATION: factored out multiple loops
  int k = 0;
  while (k < time_steps) {
    srand((unsigned)time(NULL));
    double fit;
    // evaluate fitness and evaluate
    for (int i = 0; i < NUM_PARTICLES; i++) {
      fit = optimisable(i);
      // check for personal best
      if (fit > particles[3][i]) {
        particles[0][i] = particles[2][i];
        particles[3][i] = fit;
      }
      // check for global best fitness
      if (fit > global_best_fitness) {
        global_best_fitness = fit;
        global_best_position = particles[0][i];
      }
      //update velocity
      double a = (double)rand()/(double)RAND_MAX;
      double b = (double)rand()/(double)RAND_MAX;
      particles[1][i] += a*(particles[0][i]-particles[2][i]) + b*(global_best_position-particles[2][i]);

      //update position according to new velocity
      particles[2][i] += particles[1][i];
    }

    k++;
  }
}


int main() {
  // create particles to move through the problem space
  makeParticles();

  // set initial best health and best position of particle 1 generated above
  global_best_position = particles[0][0];
  global_best_fitness = optimisable(3);

  PSO();

  printf("----------------------------------\n");

  printf("Optimal Point occurs at x = %f\n", global_best_position);

  return 0;
}

double optimisable(int x) {
  // -(x-2)^(4) + 3(x)^(2)
  return -1 * pow((particles[2][x] - 2), function_degree) + 3 * pow(particles[2][x], 2);
}
