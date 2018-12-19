#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <pthread.h>

long long number_circle_count;
long long n;
long thread_count;

pthread_mutex_t mutex;

void* Thread_toss(void* rank){
  long my_rank = (long)rank;
  long long my_n = n / thread_count;
  long long my_number_circle_count = 0;

  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_real_distribution<float> unif(-1.0, 1.0);

  for (long long i = 0; i < my_n; i++){
    float x = unif(gen);
    float y = unif(gen);
    if ((x*x + y*y) < 1.0) my_number_circle_count++;
  }

  pthread_mutex_lock(&mutex);
  number_circle_count += my_number_circle_count;
  pthread_mutex_unlock(&mutex);

  return NULL;
}

int main(int argc, char **argv){
  long thread;
  pthread_t* thread_handles;
  
  thread_count = strtol(argv[1], NULL, 10);
  n = strtoll(argv[2], NULL, 10);
  thread_handles = (pthread_t*)malloc(thread_count * sizeof(pthread_t));
  pthread_mutex_init(&mutex, NULL);
  number_circle_count = 0;
  
  for (thread = 0; thread < thread_count; thread++)
    pthread_create(&thread_handles[thread], NULL, Thread_toss, (void*)thread);

  for (thread = 0; thread < thread_count; thread++)
    pthread_join(thread_handles[thread], NULL);

  double pi_estimate = 4 * number_circle_count / ((double) n);
  std::cout << pi_estimate << std::endl;
  
  pthread_mutex_destroy(&mutex);
  free(thread_handles);

  return 0;
}