#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define USAGE "./random_gauss.x n_points"
#define pi 3.14159265359

float *reserva(int n_puntos);
void print_array(float * array, int n_puntos);
void copy(float *origen, float *destino, int n_puntos);
float *crearArreglo(int n_points);
void *asignar(float *x, float *y);
float *my_model(float *t,float a, float d, float b, float c, int lineas);
void *model(float *y_obs,float *a_walk, float *d_walk, float *b_walk, float *c_walk, float *l_walk, int n_points, int n_burn);
float likelihood(float *y_obs, float *y_model);
float randomGauss();
void escribirArchivo(int n_burn, int n_points, float *a_walk, float *d_walk, float *b_walk, float *c_walk, float *l_walk);


int main(int argc, char **argv){
  float *year, *spots, *a_walk, *d_walk, *b_walk, *c_walk, *l_walk, *y_init;
  float *b;
  int lineas,j,n_points,n_burn;

  //parametros
  n_points = atoi(argv[1]);
  n_burn = atoi(argv[2]);
  // lineas numero de lineas del archivo
  lineas = 196;
  year = crearArreglo(lineas);
  spots = crearArreglo(lineas);
  a_walk = crearArreglo(n_points);
  d_walk = crearArreglo(n_points);
  b_walk = crearArreglo(n_points);
  c_walk = crearArreglo(n_points);
  l_walk = crearArreglo(n_points);
  y_init = crearArreglo(lineas);
  //la funcion asignar asigna los valores de año y manchas de los datos a los arreglos respectivos de year y spots
  asignar(year, spots);
  //primer valor aleatorio para cada walk
  srand48(time(NULL));
  a_walk[0] = drand48();
  d_walk[0] = drand48();
  b_walk[0] = drand48();
  c_walk[0] = drand48();
  //y_init  
  y_init = my_model(year,a_walk[0],d_walk[0],b_walk[0],c_walk[0],lineas);
  
  //primer valor de l_walk
  l_walk[0] = likelihood(spots, y_init);
  
  model(spots, a_walk, d_walk, b_walk, c_walk, l_walk, n_points, n_burn);
  
  return 0;
}


float *reserva(int n_puntos){
  float *array;
  int i;
  if(!(array = malloc(n_puntos * sizeof(float)))){
    printf("Problema en reserva\n");
    exit(1);
  }
  for(i=0;i<n_puntos;i++){
    array[i] = 0.0;
  }
  return array;
}

void print_array(float * array, int n_puntos){
  int i;
  for(i=0;i<n_puntos;i++){
    printf("%f\n", array[i]);
  }
}


void copy(float *origen, float *destino, int n_puntos){
  int i;
  for(i=0;i<n_puntos;i++){
    destino[i] = origen[i];
  }
}

float *crearArreglo(int n_points){
  int i;
  float *arregloRespuesta;
  if(!(arregloRespuesta = malloc(n_points * sizeof(float)))){
    printf("Problema en reserva\n");
    exit(1);
  }
  //Inicializar en 0
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=0.0;
    }
  return arregloRespuesta;
}

void *asignar(float *x, float *y){
  FILE *in;
  int i;
  float anno, man;
  char filename[100] = "data_new.dat";
  char linea[3000];
  in = fopen(filename,"r");
  while(fgets(linea,sizeof(linea),in)!=0){
    sscanf(linea,"%f %f",&anno,&man);
    x[i] = anno;
    y[i] = man;
    i++;
  }
}

float *my_model(float *t,float a, float d, float b, float c, int lineas){
  float *model;
  int i;
  model = crearArreglo(lineas);
  for(i=0;i<lineas;i++){
    model[i] = a*cos((2*pi/d)*t[i]+b)+c;
  }
  return model;
}

float likelihood(float *y_obs, float *y_model){
  float chi_squared, suma;
  int i;
  suma = 0.0;
  for(i=0;i<196;i++){
    suma = suma + (y_obs[i]-y_model[i])*(y_obs[i]-y_model[i]);
  }
  chi_squared = (1.0/2.0)*suma;
  return chi_squared;
}

float randomGauss(){
  float random;

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  random = gsl_ran_gaussian(r, 0.1);
  gsl_rng_free (r);
  return random;
}

void *model(float *y_obs, float *a_walk, float *d_walk, float *b_walk, float *c_walk, float *l_walk, int n_points, int n_burn){
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  int i, contador;
  float a_prime,d_prime,b_prime,c_prime,l_prime,l_init,gamma,alpha,beta;
  float *nuevo,*y_init;
  y_init=crearArreglo(196);
  nuevo = crearArreglo(196);
  for(i=0;i<n_points-1;i++){
    
    a_prime = gsl_ran_gaussian(r, 0.1) + a_walk[i];

    y_init = my_model(y_init,a_walk[i],d_walk[i],b_walk[i],c_walk[i],196);
    nuevo = my_model(y_init,a_prime,d_walk[i],b_walk[i],c_walk[i],196);

    l_prime = likelihood(y_obs, nuevo);
    l_init = likelihood(y_obs,y_init);
    gamma = l_prime-l_init;

    if(gamma >= 0.0){
      a_walk[i+1] = a_prime;
      
    }else{
      beta = drand48();
      alpha = exp(gamma);
      if(beta<=alpha){
	a_walk[i+1] = a_prime;
	
      }else{
	a_walk[i+1] = a_walk[i];
      }
    }

    d_prime = gsl_ran_gaussian(r, 0.1) + d_walk[i];

    y_init = my_model(y_init,a_walk[i],d_walk[i],b_walk[i],c_walk[i],196);
    nuevo = my_model(y_init,a_walk[i],d_prime,b_walk[i],c_walk[i],196);

    l_prime = likelihood(y_obs, nuevo);
    l_init = likelihood(y_obs,y_init);
    gamma = l_prime-l_init;

    if(gamma >= 0.0){
      d_walk[i+1] = d_prime;
    }else{
      beta = drand48();
      alpha = exp(gamma);
      if(beta<=alpha){
	d_walk[i+1] = d_prime;
      }else{
	d_walk[i+1] = d_walk[i];
      }
    }

    b_prime = gsl_ran_gaussian(r, 0.1) + b_walk[i];
        
    y_init = my_model(y_init,a_walk[i],d_walk[i],b_walk[i],c_walk[i],196);
    nuevo = my_model(y_init,a_walk[i],d_walk[i],b_prime,c_walk[i],196);

    l_prime = likelihood(y_obs, nuevo);
    l_init = likelihood(y_obs,y_init);
    gamma = l_prime-l_init;

    if(gamma >= 0.0){
      b_walk[i+1] = b_prime;
    }else{
      beta = drand48();
      alpha = exp(gamma);
      if(beta<=alpha){
	b_walk[i+1] = b_prime;
      }else{
	b_walk[i+1] = b_walk[i];
      }
    }

    c_prime = gsl_ran_gaussian(r, 0.1) + c_walk[i];

    y_init = my_model(y_init,a_walk[i],d_walk[i],b_walk[i],c_walk[i],196);
    nuevo = my_model(y_init,a_walk[i],d_walk[i],b_walk[i],c_prime,196);

    l_prime = likelihood(y_obs, nuevo);
    l_init = likelihood(y_obs,y_init);
    gamma = l_prime-l_init;

    if(gamma >= 0.0){
      c_walk[i+1] = c_prime;
    }else{
      beta = drand48();
      alpha = exp(gamma);
      if(beta<=alpha){
	c_walk[i+1] = c_prime;
      }else{
	c_walk[i+1] = c_walk[i];
      }
    }

    nuevo = my_model(y_init,a_walk[i+1],d_walk[i+1],b_walk[i+1],c_walk[i+1],196);
    l_walk[i+1] = likelihood(y_obs,y_init);
   
    if(i>n_burn){
      printf("%f %f %f %f %f\n", a_walk[i+1], d_walk[i+1],b_walk[i+1],c_walk[i+1],l_walk[i+1]);
    }
  }
  //escribirArchivo(n_burn,n_points,a_walk,d_walk, b_walk, c_walk, l_walk);
}


void escribirArchivo(int n_burn, int n_points, float *a_walk, float *d_walk, float *b_walk, float *c_walk, float *l_walk){
  FILE *in;
  int i, var, test;
  char filename[100]="data.dat";

  in = fopen(filename, "w");
  for(i=n_burn;i<n_points;i++){
    fprintf(in, "%f %f %f %f %f\n", a_walk[i+1], d_walk[i+1],b_walk[i+1],c_walk[i+1],l_walk[i+1]);
  }
  fclose(in);
}
