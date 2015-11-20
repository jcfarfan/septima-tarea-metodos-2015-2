#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define USAGE "./random_gauss.x n_points"

float *pasoRungerKutta4(float *ks,float posicionesX, float posicionesY, float h,float alpha, float beta, float gamma, float epsilon);
void *loopRungeKutta4(float *ks, float *x, float *y, float *t,float h,float alpha, float beta, float gamma, float epsilon, int n);
float *crearArregloCero(int n_points);
float xprime(float a, float b, float alpha, float beta);
float yprime(float a, float b, float gamma, float epsilon);
//void *model(float h, float *ks, float *datost, float  *datosx, float *datosy,float *alpha_walk, float *beta_walk, float *gamm_walk, float *epsilon_walk, float *l_walk, int n_points, int n_burn, float n);
void *asignar(float *datost, float *datosx, float *datosy);
float randomGauss();


int main(int argc, char **argv){
  float *x,*y,*t,*ks,*alpha_walk, *beta_walk, *gamma_walk, *epsilon_walk, *datost, *datosx, *datosy;
  float h,alpha,beta,gamma,epsilon;
  int n_points,i,n_puntos,n_burn;

  //96 lineas tienen los datos
  n_points = 96;
  //parametros entrados por consola
  n_puntos = atoi(argv[1]);
  n_burn = atoi(argv[2]);
  //arreglos para guardar los datos del tiempo, presas y depredadores
  datost = crearArregloCero(96);
  datosx = crearArregloCero(96);
  datosy = crearArregloCero(96);
  //asignar los datos del archivo a los arreglos datost, datosx, datosy
  asignar(datost,datosx,datosy);
  //arreglos para guardar los mejores parametros estimados
  alpha_walk = crearArregloCero(n_puntos);
  beta_walk = crearArregloCero(n_puntos);
  gamma_walk = crearArregloCero(n_puntos);
  epsilon_walk = crearArregloCero(n_puntos);
  //valores aleatorios para el modelo
  srand48(time(NULL));
  alpha_walk[0] = drand48()*0.2-0.1;
  beta_walk[0] = drand48()*0.2-0.1;
  gamma_walk[0] = drand48()*0.2-0.1;
  epsilon_walk[0] = drand48()*0.2-0.1;
  //inicializar los arreglos 
  alpha = alpha_walk[0];
  beta = beta_walk[0];
  gamma = gamma_walk[0];
  epsilon = epsilon_walk[0];
  //Arreglos para el RungeKutta
  ks = crearArregloCero(2);
  x = crearArregloCero(n_points);
  x[0] = 15.0;
  y = crearArregloCero(n_points);
  y[0] = 13.0;
  t = crearArregloCero(n_points);
  t[0] = 0.006;
  t[95] = 0.799;
  h = (t[95]-t[0])/n_points ;
  //RungeKutta
  loopRungeKutta4(ks,x,y,t,h,alpha,beta,gamma,epsilon,n_points);
  //Rungekutta y markov deben ir en un loop

  for(i=0;i<n_points;i++){
    printf("%f %f %f\n", t[i],x[i],y[i]);
  }

  return 0;
}


float *crearArregloCero(int n_points){
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

float *pasoRungerKutta4(float *ks,float posicionesX, float posicionesY, float h,float alpha, float beta, float gamma, float epsilon){
  float k1,k2,k3,k4,kF;
  float l1,l2,l3,l4,lF;
  float inter1x,inter2x,inter3x, inter1y,inter2y,inter3y;

  k1 = xprime(posicionesX, posicionesY,alpha,beta);
  l1 = yprime(posicionesX, posicionesY,gamma,epsilon);
  
  inter1x = posicionesX+0.5*h*k1;
  inter1y = posicionesY+0.5*h*l1;
  k2 = xprime(inter1x,inter1y,alpha,beta);
  l2 = yprime(inter1x,inter1y,gamma,epsilon);
  
  inter2x = posicionesX+0.5*h*k2;
  inter2y = posicionesY+0.5*h*l2;
  k3 = xprime(inter2x,inter2y,alpha,beta);
  l3 = yprime(inter2x,inter2y,gamma,epsilon);

  inter3x = posicionesX+h*k3;
  inter3y = posicionesY+h*l3;
  k4 = xprime(inter3x,inter3y,alpha,beta);
  l4 = yprime(inter3x,inter3y,gamma,epsilon);

  ks[0] = (1.0/6.0)*(k1+2*k2+2*k3+k4);
  ks[1] = (1.0/6.0)*(l1+2*l2+2*l3+l4);   

  return ks;
}

void *loopRungeKutta4(float *ks,float *x, float *y, float *t,float h,float alpha, float beta, float gamma, float epsilon, int n){
  int i;
  float kF,lF; 
  for(i=1; i<n; i++){

    ks = pasoRungerKutta4(ks,x[i-1],y[i-1],h,alpha,beta,gamma,epsilon);

    x[i] = x[i-1]+ ks[0];
    y[i] = y[i-1]+ ks[1];
    t[i] = t[i-1]+h;
  }
}

float xprime(float a, float b, float alpha, float beta){
  float xp;
  xp = a*(alpha-beta*b);
  return xp;
}

float yprime(float a, float b, float gamma, float epsilon){
  float yp;
  yp = -b*(gamma-epsilon*a);
  return yp;
}
/*
void *model(float h,float *ks, float *datost, float  *datosx, float *datosy,float *alpha_walk, float *beta_walk, float *gamma_walk, float *epsilon_walk, float *l_walk, int n_points, int n_burn, float n){
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  int i, contador;
  float alpha_prime,beta_prime,gamma_prime,epsilon_prime,l_prime,l_init,ga1,ga2,al,be;
  float *nuevo_x, *nuevo_y,*x_init, *y_init;
  x_init = crearArregloCero(96);
  y_init = crearArregloCero(96);
  nuevo_x = crearArregloCero(96);
  nuevo_y = crearArregloCero(96);

  for(i=0;i<n_points-1;i++){
    
    alpha_prime = gsl_ran_gaussian(r, 0.1) + alpha_walk[i];
    beta_prime = gsl_ran_gaussian(r, 0.1) + beta_walk[i];
    gamma_prime = gsl_ran_gaussian(r, 0.1) + gamma_walk[i];
    epsilon_prime = gsl_ran_gaussian(r, 0.1) + epsilon_walk[i];

    loopRungeKutta4(ks,x_init,y_init,h, alpha_prime,  beta_prime, gamma_prime, epsilon_prime, n);
    loopRungeKutta4(ks,nuevo_x,nuevo_y,h, alpha_walk[i],  beta_walk[i], gamma_walk[i], epsilon_walk[i], n);
    
    l_prime1 = likelihood(datosx, x_init);
    l_prime2 = likelihood(datosy, y_init);

    l_init1 = likelihood(nuevo_x,x_init);
    l_init2 = likelihood(nuevo_y,y_init);

    ga1 = l_prime1/l_init1;
    ga2 = l_prime2/l_init2;

    if(ga1 >= 1.0){
      alpha_walk[i+1] = alpha_prime;
      beta_walk[i+1] = beta_prime;
      gamma_walk[i+1] = gamma_prime;
      epsilon_walk[i+1] = epsilon_prime;
      l_walk[i+1] = l_prime;
    }else{
      be1 = drand48();
      al1 = exp(gamma);
      if(be1<=al1){
	alpha_walk[i+1] = alpha_prime;
	beta_walk[i+1] = beta_prime;
	gamma_walk[i+1] = gamma_prime;
	epsilon_walk[i+1] = epsilon_prime;
	l_walk[i+1] = l_prime;
      }else{
	alpha_walk[i+1] = alpha_walk[i];
	beta_walk[i+1] = beta_walk[i];
	gamma_walk[i+1] = gamma_walk[i];
	epsilon_walk[i+1] = epsilon_walk[i];
	l_walk[i+1] = l_walk[i];
      }
    }
    
    if(i>n_burn){
      printf("%f %f %f %f %f\n", alpha_walk[i], beta_walk[i],gamma_walk[i],epsilon_walk[i],l_walk[i]);
    }
  }
}
*/
void *asignar(float *datost, float *datosx, float *datosy){
  FILE *in;
  int i;
  float time, prey, predator;
  char filename[100] = "datos.dat";
  char linea[3000];
  in = fopen(filename,"r");
  while(fgets(linea,sizeof(linea),in)!=0){
    sscanf(linea,"%f %f %f",&time,&prey,&predator);
    datost[i] = time;
    datosx[i] = prey;
    datosy[i] = predator;
    i++;
  }
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
