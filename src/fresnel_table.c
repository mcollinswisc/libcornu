#include <pthread.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "progress.h"
#include "fresnel.h"

/**
 * The number of entries to be claimed by a single thread.
 */
const int thread_block = 16;

/**
 * Error bound to which to calculate table entries.
 */
double epsilon = 1e-8;

/**
 * Frequency to update progress meter
 */
int progress_print_interval = 32;


/**
 * The tables.
 */
double *S, *C;

/**
 * The next unclaimed table entry.
 */
int next;

/**
 * The last table entry being considered.
 */
int end;

/**
 * The next power of 2 greater than end
 */
int endpow;

pthread_mutex_t queue_lock;

// Progress tracking
time_t start_time;
int progress, total;
pthread_mutex_t progress_lock;

/**
 * Calculate floor(log2(n)) using integer arithmetic
 */
static int ilog2(int n)
{
  int m, l;

  for(m = 1, l = 0; m < n; m <<= 1, l++);

  if(m == n) {
    return l;
  }
  else {
    return l-1;
  }
}

/**
 * Get the table entry to be calculated for job i
 *
 * @note N must be a power of 2
 */
static int entry(int i, int N)
{
  int j, level;

  if(i == 0) {
    return 0;
  }

  level = ilog2(i);
  j = i % (1 << level);
  level = 1 << (ilog2(N) - level - 1);

  return level * (2*j + 1);
}

static void *thread_calculate(void *p)
{
  int i, j, task;
  pthread_t self;

  self = pthread_self();

  while(1) {
    // Get the next task
    pthread_mutex_lock(&queue_lock);

    // Quit when nothing left
    if(next >= endpow) {
      pthread_mutex_unlock(&queue_lock);
      break;
    }

    // Claim next task
    task = next;
    next += thread_block;
    pthread_mutex_unlock(&queue_lock);

    // Perform task
    for(i = task; i < task + thread_block; i++) {
      j = entry(i, endpow);
      if(j >= end) {
        continue;
      }

      S[j] = fresnel_sine_series_exact(M_PI * j, epsilon);
      C[j] = fresnel_cosine_series_exact((2*j + 1) * M_PI/2, epsilon);

      // Update progress
      pthread_mutex_lock(&progress_lock);
      progress++;

      if(progress % progress_print_interval == 0) {
        print_progress(start_time, progress, total, 2);
      }
      pthread_mutex_unlock(&progress_lock);
    }
  }

  return NULL;
}

int main(int argc, char **argv)
{
  pthread_t *threads;
  int i, nthreads, N;

  const char *out_path;
  FILE *tout;

  assert(argc == 4);

  assert(sscanf(argv[1], "%d", &N) == 1);
  assert(sscanf(argv[2], "%d", &nthreads) == 1);
  out_path = argv[3];

  // Print out arguments
  printf("Precomputing table of Fresnel integrals\n");
  printf("N = %d\n", N);
  printf("nthreads = %d\n", nthreads);
  printf("Smax = %lg\n", M_PI * (N+1));
  printf("Cmax = %lg\n", (2 * (N+1) + 1) * M_PI/2);
  printf("outpath = \"%s\"\n", out_path);


  pthread_mutex_init(&queue_lock, NULL);
  pthread_mutex_init(&progress_lock, NULL);

  // Allocate space
  S = (double *)malloc(sizeof(double) * N);
  C = (double *)malloc(sizeof(double) * N);
  for(i = 0; i < N; i++) {
    S[i] = C[i] = -1;
  }

  // Initialize queue
  next = 0;
  end = N;

  endpow = (1 << ilog2(end));
  if(endpow != end) {
    endpow *= 2;
  }

  // Initialize progress
  start_time = time(NULL);
  progress = 0;
  total = N;

  progress_print_interval = N / 1024;
  if(progress_print_interval < 4) {
    progress_print_interval = 4;
  }

  // Start threads
  threads = (pthread_t *)malloc(sizeof(pthread_t) * nthreads);
  for(i = 0; i < nthreads; i++) {
    pthread_create(&threads[i], NULL, thread_calculate, NULL);
  }

  // Wait on threads
  for(i = 0; i < nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  // Output table
  tout = fopen(out_path, "w+");

  fprintf(tout, "// fresnel_table.h\n");
  fprintf(tout, "// Table auto-generated with command:\n");
  fprintf(tout, "//   ");
  for(i = 0; i < argc; i++) {
    fprintf(tout, "%s ", argv[i]);
  }
  fprintf(tout, "\n");

  fprintf(tout, "\n#define FRESNEL_TABLE_ENTRIES (%d)\n\n", N);

  fprintf(tout, "static const double fresnel_table_S[FRESNEL_TABLE_ENTRIES] = {\n");
  for(i = 0; i < N; i++) {
    assert(S[i] >= 0);
    fprintf(tout, "  %.10e,\n", S[i]);
  }
  fprintf(tout, "};\n\n");

  fprintf(tout, "static const double fresnel_table_C[FRESNEL_TABLE_ENTRIES] = {\n");
  for(i = 0; i < N; i++) {
    assert(C[i] >= 0);
    fprintf(tout, "  %.10e,\n", C[i]);
  }
  fprintf(tout, "};\n\n");

  fclose(tout);

  return 0;
}
