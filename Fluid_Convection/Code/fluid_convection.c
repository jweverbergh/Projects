#include "poisson.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

#define M_PI 3.14159265358979323846

//return 1 if the point at (x,y) is in the mixer, 0 otherwise
int xi_omega(param *p, double x, double y){

    x = x - p->L/2;
    y = y - p->L/2;
    double theta = atan2(y,x) - p->ws * p->time;
    double norm = sqrt(x*x + y*y);
    double r_theta = p->a*cos(3*theta);

    return p->adding_blades * ((norm <= r_theta) | (norm <= p->d/2));
}

//return 1 if the point at (x,y) is in the cylinder, 0 otherwise
int xi_omegaCYL(param *p, double x, double y){
    x = x - p->L/2;
    y = y - p->L/2;
    double norm = sqrt(x*x + y*y);

    return p->adding_blades * (norm <= p->D/2);
}

//fill the u_s and v_s array only once for all points in the grid as the angular speed is constant
void init_uv_mixing(param* p, double** u_s, double** v_s){

    double speed;
    double norm;
    double theta;
    double x;
    double y;
    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx; i++){
            x = i*p->h - p->L/2;
            y = (j-1)*p->h+p->h/2 - p->L/2;

            norm = sqrt(x*x + y*y);
            theta = atan2(y,x);
            speed = p->ws * norm;

            u_s[j][i] = - speed * sin(theta);
        }
    }

    for (int j = 1; j < p->ny; j++){
        for (int i = 1; i < p->nx + 1; i++){
            x = (i-1)*p->h+p->h/2  - p->L/2;
            y = j*p->h  - p->L/2;

            norm = sqrt(x*x + y*y);
            theta = atan2(y,x);
            speed = p->ws * norm;

            v_s[j][i] = speed * cos(theta);
        }
    }
}

//compute the advective form of the convection term in the conservation of momentum, and store them in Hx and Hy
void convection_advform(param *p, double **Hx, double **Hy, double **u, double **v){

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx; i++){
            Hx[j][i] = (u[j][i+1] + u[j][i]) * (u[j][i+1] - u[j][i]) / (4*p->h) 
                      +(u[j][i] + u[j][i-1]) * (u[j][i] - u[j][i-1]) / (4*p->h) 
                      +(u[j+1][i] - u[j][i]) * (v[j][i+1] + v[j][i]) / (4*p->h) 
                      +(u[j][i] - u[j-1][i]) * (v[j-1][i+1] + v[j-1][i]) / (4*p->h) ;
        }
    }

    for (int j = 1; j < p->ny; j++){
        for (int i = 1; i < p->nx + 1; i++){
            Hy[j][i] = (v[j+1][i] + v[j][i]) * (v[j+1][i] - v[j][i]) / (4*p->h) 
                      +(v[j][i] + v[j-1][i]) * (v[j][i] - v[j-1][i]) / (4*p->h) 
                      +(v[j][i+1] - v[j][i]) * (u[j+1][i] + u[j][i]) / (4*p->h) 
                      +(v[j][i] - v[j][i-1]) * (u[j+1][i-1] + u[j][i-1]) / (4*p->h) ;
        }
    }
}

//set the ghost point values for u and v
void ghost_point_uv(param *p, double **u, double **v){

    for (int i = 1; i < p->nx; i++){
        u[0][i] = -0.2 * (15*u[1][i] - 5*u[2][i] + u[3][i]); //lower surface
        u[p->ny+1][i] = u[p->ny][i]; // upper surface
    }

    for (int j = 1; j < p->ny; j++){
        v[j][0] = -0.2 * (15*v[j][1] - 5*v[j][2] + v[j][3]); //left surface
        v[j][p->nx+1] = -0.2 * (15*v[j][p->nx] - 5*v[j][p->nx-1] + v[j][p->nx-2]); //right surface
    }
}

//set the ghost point values for Theta
void ghost_point_T(param *p, double **Theta){

    for (int i = 1; i < p->nx + 1; i++){
        Theta[0][i] = p->DT / p->H * p->h + Theta[1][i]; //lower surface
        Theta[p->ny+1][i] = (p->l0 / p->h - 0.5) / (p->l0 / p->h + 0.5) * Theta[p->ny][i]; // upper surface
    }

    for (int j = 1; j < p->ny + 1; j++){
        Theta[j][0] = Theta[j][1]; //left surface
        Theta[j][p->nx+1] = Theta[j][p->nx]; //right surface
    }
}

//set the boundary conditions
void boundary_condition(param *p, double **u, double **v){

    for (int i = 1; i < p->nx+1; i++){
        v[0][i] = 0; // lower surface        
        v[p->ny][i] = 0; // upper surface
    }

    for (int j = 1; j < p->ny+1; j++){
        u[j][0] = 0; // left surface
        u[j][p->nx] = 0; // right surface
    }
}

//for the case 2 with mixer
//compute the first part of the discretization scheme (i.e. the update of u_star and v_star)
// factor = 1.5 gives Adams Bashfort
// factor = -1 gives Explicit euler
void step_1(double factor, param *p, double **Hx, double **Hy, double **u, double **u_star, double **v, double **v_star, double **P, double **Theta, double** u_s, double** v_s){
    int value;
    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx; i++){

            //check whether the point is in the mixer or not (no mixing if value=0)
            value = xi_omega(p, i*p->h, (j-1)*p->h+p->h/2);

            if (value==0){
                u_star[j][i] = u[j][i] + p->dt * (  0.5 * Hx[j][i] 
                                                - (P[j][i+1] - P[j][i]) / p->h
                                                + p->nu * (u[j][i+1] + u[j+1][i] - 4*u[j][i] + u[j-1][i] + u[j][i-1]) / (p->h) / (p->h)
                                                );
            }else{
                u_star[j][i] = (    p->dtau * u[j][i]
                                    + p->dt * u_s[j][i]
                                    + p->dtau * p->dt * (  0.5 * Hx[j][i] 
                                                - (P[j][i+1] - P[j][i]) / p->h
                                                + p->nu * (u[j][i+1] + u[j+1][i] - 4*u[j][i] + u[j-1][i] + u[j][i-1]) / (p->h) / (p->h))
                                ) / (p->dtau + p->dt);
            }
        }
    }

    for (int j = 1; j < p->ny; j++){
        for (int i = 1; i < p->nx + 1; i++){

            //check whether the point is in the mixer or not (no mixing if value=0)
            value = xi_omega(p, (i-1)*p->h+p->h/2, j*p->h);

            if (value==0){
                v_star[j][i] = v[j][i] + p->dt * (  0.5 * Hy[j][i] 
                                               - (P[j+1][i] - P[j][i]) / p->h
                                               + p->nu * (v[j+1][i] + v[j][i+1] - 4*v[j][i] + v[j][i-1] + v[j-1][i]) / (p->h) / (p->h)
                                               + p->beta * 0.5 * (Theta[j+1][i] + Theta[j][i]) * p->g
                                             );
            }else{
                v_star[j][i] = (    p->dtau * v[j][i]
                                    + p->dt * v_s[j][i]
                                    + p->dtau * p->dt * (  0.5 * Hy[j][i] 
                                               - (P[j+1][i] - P[j][i]) / p->h
                                               + p->nu * (v[j+1][i] + v[j][i+1] - 4*v[j][i] + v[j][i-1] + v[j-1][i]) / (p->h) / (p->h)
                                               + p->beta * 0.5 * (Theta[j+1][i] + Theta[j][i]) * p->g)
                                ) / (p->dtau + p->dt);
            }
        }
    }

    convection_advform(p, Hx, Hy, u, v);

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx; i++){
            //check whether the point is in the mixer or not (no mixing if value=0)
            value = xi_omega(p, i*p->h, (j-1)*p->h+p->h/2);

            if (value==0){
                u_star[j][i] -= factor * p->dt * Hx[j][i];
            }else{
                u_star[j][i] -= factor * p->dtau * p->dt * Hx[j][i] / (p->dtau + p->dt);
            }
        }
    }

    for (int j = 1; j < p->ny; j++){
        for (int i = 1; i < p->nx + 1; i++){
            //check whether the point is in the mixer or not (no mixing if value=0)
            value =  xi_omega(p, (i-1)*p->h+p->h/2, j*p->h);

            if (value==0){
                v_star[j][i] -= factor * p->dt * Hy[j][i];
            }else{
                v_star[j][i] -= factor * p->dtau * p->dt * Hy[j][i] / (p->dtau + p->dt);
            }
        }
    }
}

//for the case 2 with mixer
//compute the first part of the discretization scheme (i.e. the update of u_star and v_star)
void step_2(double factor, param *p, double **u, double **u_star, double **v, double **v_star, double **P, double **phi, 
double **Theta, double **lapl_Theta, double **HT_n_1, double **HT_n){

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++){
            //laplacian of Theta
            lapl_Theta[j][i] = (Theta[j][i+1] + Theta[j+1][i] - 4*Theta[j][i] + Theta[j-1][i] + Theta[j][i-1]) / (p->h) / (p->h);
            
            //advective form for Theta
            HT_n[j][i] = (u[j][i]*(Theta[j][i+1]-Theta[j][i]) + u[j][i-1]*(Theta[j][i]-Theta[j][i-1])) / (2*p->h)
                      + (v[j][i]*(Theta[j+1][i]-Theta[j][i]) + v[j-1][i]*(Theta[j][i]-Theta[j-1][i])) / (2*p->h);
        }
    }

    //Theta_(n+1)
    int value;
    double sum_Theta = 0;
    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++){
            //check whether the point is in the mixer or not (no mixing if value=0)
            value = xi_omega(p, (i-1)*p->h+p->h/2, (j-1)*p->h+p->h/2);

            if (value==0){
                Theta[j][i] += p->dt * (
                                        - factor * HT_n[j][i] + 0.5 * HT_n_1[j][i]  
                                        + p->alpha * lapl_Theta[j][i]);
            }else{
                Theta[j][i] = (    p->dtau * Theta[j][i]
                                    + p->dt * p->Theta_s
                                    + p->dtau * p->dt * (
                                                - factor * HT_n[j][i] + 0.5 * HT_n_1[j][i]  
                                                + p->alpha * lapl_Theta[j][i] )
                                ) / (p->dtau + p->dt);
            }

            //compute the average of theta on the cylinder
            value = xi_omegaCYL(p, (i-1)*p->h+p->h/2, (j-1)*p->h+p->h/2);
            if (value==1){
                sum_Theta += Theta[j][i];
            }
        } 
    }

    //update theta_s for the next iteration
    p->Theta_s = sum_Theta*p->h*p->h / (M_PI * p->a * p->a);

    //u_(n+1)
    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx; i++){
            u[j][i] = u_star[j][i] - p->dt * (phi[j][i+1] - phi[j][i]) / p->h;
        }
    }

    //v_(n+1)
    for (int j = 1; j < p->ny; j++){
        for (int i = 1; i < p->nx + 1; i++){
            v[j][i] = v_star[j][i] - p->dt * (phi[j+1][i] - phi[j][i]) / p->h;
        }
    }

    //P_(n+1)
    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++){
            P[j][i] += phi[j][i];
        }
    }
}

//read data file to rerun simulation from another timestep than 0
void readfile(param *p, double **Theta, double **u, double **v, double **P, char *prefix, int iter){
    char filename[80] = {};
    sprintf(filename, "./save_txt/%s%d.txt", prefix, iter);
    FILE *file = fopen(filename, "r");
    int err = 0;
    err += 1;

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++)
            err = fscanf(file, "%lf", Theta[j]+i);
    }

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 0; i < p->nx + 1 ; i++)
            err = fscanf(file, "%lf", u[j]+i);
    }

    for (int j = 0; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1 ; i++)
            err = fscanf(file, "%lf", v[j]+i);
    }

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++)
            err = fscanf(file, "%lf", P[j]+i);
    }
    fclose(file);
}

//write results of simulations in .txt files
void write(param *p, double **Theta, double **u, double **v, double **P, char *prefix, int iter){
    char filename[80] = {};
    sprintf(filename, "./save_txt/%s%d.txt", prefix, iter);
    FILE *file = fopen(filename, "w");

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++)
            fprintf(file, "%.5le ", Theta[j][i]);
        fprintf(file, "\n");
    }

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 0; i < p->nx + 1 ; i++)
            fprintf(file, "%.5lf ", u[j][i]);
        fprintf(file, "\n");
    }

    for (int j = 0; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1 ; i++)
            fprintf(file, "%.5le ", v[j][i]);
        fprintf(file, "\n");
    }

    for (int j = 1; j < p->ny + 1; j++){
        for (int i = 1; i < p->nx + 1; i++)
            fprintf(file, "%.5le ", P[j][i]);
        fprintf(file, "\n");
    }
    fclose(file);
}



int main(int argc, char *argv[]){


    //// PARAMETERS

    //save and write parameters
    int save_results = atoi(argv[1]); //1 to save results in txt files, 0 otherwise
    int verbose = atoi(argv[2]); //1 to print parameters and some intermediate results, 0 otherwise

    //spatial parameters
    param *p = (param *)malloc(sizeof(param));
    p->L = 1.;
    p->H = 1.5 * p->L;
    p->h = atof(argv[3]); //0.005; 
    p->nx = (int) (p->L / p->h);
    p->ny = (int) (p->H / p->h);

    //heat parameters
    p->DT = 1.;
    p->k = 1.;
    p->qw = p->k * p->DT / p->H;

    //other parameters
    p->g = 9.81;
    p->Gr = 1e10;
    p->Pr = 4.;
    p->U = 1.;
    p->nu = p->U * p->H / sqrt(p->Gr);
    p->alpha = p->nu / p->Pr;
    p->beta = pow(p->U, 2) / (p->g * p->DT * p->H);
    p->l0 = 1e-3 * p->H;
    
    //time parameters
    double CFL = atof(argv[6]); //0.5
    p->dt = CFL * p->h / (0.06 + 0.06); 
    double t_end = atof(argv[4]); // adimensional end time
    p->nt = (int) (p->H / p->U / p->dt * t_end); // simulate until adimensional time = t_end

    //parameters for the save files
    double frame_save = 2.; // 2 save per adimensional time
    int n_save =(int) (t_end * frame_save);
    int save_freq = (int) (p->nt / n_save);
    int save_idx = 0;
    int local_idx = save_idx * save_freq;

    //mixer parameters
    p->adding_blades = atoi(argv[5]); //1 if mixer is present, 0 otherwise
    p->a = 3*p->L/10;
    p->D = 2*p->a;
    p->d = p->D/5;
    p->ws = 0.1*p->U/p->H;
    p->dtau = p->dt/1000;
    p->Theta_s = 0;

    //stability
    double r = p->nu*p->dt/p->h/p->h;

    //set the name for mixer or no mixer
    char* prefix = "nomixer";
    if (p->adding_blades == 1){
        prefix = "mixer";
    }

    //print some parameters
    if (verbose==1){
        printf("nx : %d, ny : %d\n", p->nx, p->ny);
        printf("p->beta = %lf\n", p->beta);
        printf("p->U = %lf\n", p->U);
        printf("dt = %lf, nt = %d\n", p->dt, p->nt);
        printf("U/H dt = %lf\n", p->U / p->H * p->dt);
        printf("Mixer : a : %lf, D : %lf, d : %lf, ws : %lf, dtau : %lf\n", p->a, p->D, p->d, p->ws, p->dtau);
        printf("Stability : r = %.16lf\n", r);
        printf("n_save = %d, save_frq = %d \n", n_save, save_freq);
    }

    //save parameters
    if (save_results==1){
        char filename[80] = {};
        sprintf(filename, "./save_txt/%s_param.txt", prefix); 
        FILE *file = fopen(filename, "w"); 
        fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %d %lf %d\n", 
        p->L, p->H, p->DT, p->U, p->nu, p->beta, p->qw, p->k, p->g, p->alpha, p->l0, p->Gr, p->Pr, p->nx, p->ny, p->h, p->dt, p->nt, frame_save, n_save);
        fclose(file);
    }


    //// INITIALIZE ARRAYS


    double **u = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) u[i] = (double *)calloc((p->nx + 1), sizeof(double));

    double **u_star = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) u_star[i] = (double *)calloc((p->nx + 1), sizeof(double));

    double **u_s = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) u_s[i] = (double *)calloc((p->nx + 1), sizeof(double));

    double **Hx = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) Hx[i] = (double *)calloc((p->nx + 1), sizeof(double));

    double **v = (double **)malloc((p->ny + 1) * sizeof(double *));
    for (int i = 0; i < p->ny + 1; i++) v[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **v_star = (double **)malloc((p->ny + 1) * sizeof(double *));
    for (int i = 0; i < p->ny + 1; i++) v_star[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **v_s = (double **)malloc((p->ny + 1) * sizeof(double *));
    for (int i = 0; i < p->ny + 1; i++) v_s[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **Hy = (double **)malloc((p->ny + 1) * sizeof(double *));
    for (int i = 0; i < p->ny + 1; i++) Hy[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **P = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) P[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **phi = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) phi[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **Theta = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) Theta[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double **lapl_Theta = (double **)malloc((p->ny + 2) * sizeof(double *));
    for (int i = 0; i < p->ny + 2; i++) lapl_Theta[i] = (double *)calloc((p->nx + 2), sizeof(double));

    double ***HT = (double ***)malloc(2*sizeof(double **)); 
    for (int j = 0; j < 2; j++) {
        HT[j] = (double **)malloc((p->ny + 2) * sizeof(double *));
        for (int i = 0; i < p->ny + 2; i++) HT[j][i] = (double *)calloc((p->nx + 2), sizeof(double));
    }


    //// SOLVE EQUATIONS


    PetscInitialize(&argc, &argv, 0, 0);

    //initialize data to solve Poisson equation
    Poisson_data* data = (Poisson_data *)malloc(sizeof(Poisson_data));
    initialize_poisson_solver(data, p, u_star, v_star, phi);

    //compute u_s and v_s (velocity of the mixer)
    init_uv_mixing(p, u_s, v_s);

    //set the boundary conditions
    boundary_condition(p, u, v);
    boundary_condition(p, u_star, v_star);

    //read results from another file to begin from this timestep and not from 0, otherwise save initial arrays
    if (save_idx != 0){
        readfile(p, Theta, u, v, P, prefix, save_idx);
    }
    else{
        if (save_results==1){
            write(p, Theta, u, v, P, prefix, 0);
        }
    }

    //initialize the first step of the discretization scheme (with Adams Bashfort)
    p->time = (local_idx + 1) * p->dt; //real time of the simulation; adding 1 to obtain v_s at timestep n+1 in the first discretized equation (to update u_star and v_star)
    ghost_point_uv(p, u, v);
    step_1(-1., p, Hx, Hy, u, u_star, v, v_star, P, Theta, u_s, v_s);
    ghost_point_uv(p, u_star, v_star);
    poisson_solver(data, p, u_star, v_star, phi);
    ghost_point_T(p, Theta);
    step_2(-1., p, u, u_star, v, v_star, P, phi, Theta, lapl_Theta , HT[local_idx%2], HT[(local_idx+1)%2]);

    //record the time of the simulation
    struct timeval t1, t0;
    gettimeofday(&t0, NULL);

    //do the time iterations of the scheme
    int idx_write = save_idx+1;
    for (int i = local_idx + 2; i <= p->nt; i++){
        p->time += p->dt;
        ghost_point_uv(p, u, v);
        step_1(1.5, p, Hx, Hy, u, u_star, v, v_star, P, Theta, u_s, v_s);
        poisson_solver(data, p, u_star, v_star, phi);
        ghost_point_T(p, Theta);
        step_2(1.5, p, u, u_star, v, v_star, P, phi, Theta, lapl_Theta, HT[(i+1)%2], HT[i%2]);

        if ((save_results==1) & (i%save_freq == 0)){
            write(p, Theta, u, v, P, prefix, idx_write);
            idx_write += 1;
        }

        if (verbose==1){
            if (i%50 == 0) printf("[%d / %d] T=%lf, ", i, p->nt, Theta[1][1]);
        }
    }

    //save the time of the simulation
    gettimeofday(&t1, NULL);
    if (verbose==1){
        printf("\nThe simulation took %lu ms\n", (t1.tv_sec - t0.tv_sec) * 1000 + (t1.tv_usec - t0.tv_usec)/1000); 
        printf("\nTheta[1][nx/2] = %lf\n", Theta[1][p->nx/2]); // check to see if blew up or not
        printf("u[1][nx/2] = %lf\n", u[1][p->nx/2]); // check to see if blew up or not
    }

    PetscFinalize();

    
    //// FREE


    for (int i = 0; i < p->ny + 2; i++) free(u[i]);
    free(u);
    for (int i = 0; i < p->ny + 2; i++) free(u_star[i]);
    free(u_star);
    for (int i = 0; i < p->ny + 2; i++) free(u_s[i]);
    free(u_s);
    for (int i = 0; i < p->ny + 2; i++) free(Hx[i]);
    free(Hx);
    for (int i = 0; i < p->ny + 1; i++) free(v[i]);
    free(v);
    for (int i = 0; i < p->ny + 1; i++) free(v_star[i]);
    free(v_star);
    for (int i = 0; i < p->ny + 1; i++) free(v_s[i]);
    free(v_s);
    for (int i = 0; i < p->ny + 1; i++) free(Hy[i]);
    free(Hy);
    for (int i = 0; i < p->ny + 2; i++) free(P[i]);
    free(P);
    for (int i = 0; i < p->ny + 2; i++) free(phi[i]);
    free(phi);
    for (int i = 0; i < p->ny + 2; i++) free(Theta[i]);
    free(Theta);
    for (int i = 0; i < p->ny + 2; i++) free(lapl_Theta[i]);
    free(lapl_Theta);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < p->ny + 2; i++) free(HT[j][i]);
        free(HT[j]);
    }
    free(HT);
    free(p);
}
