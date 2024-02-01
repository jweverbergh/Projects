#ifndef _POISSON_H_
#define _POISSON_H_

#include <petsc.h>
#include <petscsys.h>

//structure storing petsc vectors
typedef struct {

	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} Poisson_data;

//structure storing problem parameters
typedef struct param {
   double L;
   double H;
   double DT;
   double U;
   double nu;
   double beta;
   double qw;
   double k;
   double g;
   double alpha;
   double l0;
   double Gr;
   double Ec;
   double Re;
   double Pr;
   double CFL;
   double r;
   double T_inf;
   int nx;
   int ny;
   double h;
   double dt;
   int nt;
   double a;
   double D;
   double d;
   double ws;
   double time;
   double dtau;
   int adding_blades;
   double Theta_s;
} param;


PetscErrorCode initialize_poisson_solver(Poisson_data* data, param *p, double **u_star, double **v_star, double **phi);
void poisson_solver(Poisson_data *data, param *p, double **u_star, double **v_star, double **phi);
void free_poisson_solver(Poisson_data* data);

#endif

