#include <mpi.h>
#include "poisson.h"

//Fill in rhs of Poisson equation
/*Called by poisson_solver at each time step*/
void computeRHS(double *rhs, PetscInt rowStart, PetscInt rowEnd, param *p, double **u_star, double **v_star){
    int r = 0;
    for (int j = 1; j < p->ny+1; j++){
        for (int i = 1; i < p->nx+1; i++){
            rhs[r] = p->h / p->dt * (u_star[j][i] - u_star[j][i-1] + v_star[j][i] - v_star[j-1][i]);
            r += 1;
        }
    }
    rhs[p->nx-1] = 0;
}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
void poisson_solver(Poisson_data *data, param *p, double **u_star, double **v_star, double **phi){

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = data->sles;
    Vec b = data->b;
    Vec x = data->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(rhs, rowStart, rowEnd, p, u_star, v_star);
    VecRestoreArray(b, &rhs);

    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);

    VecGetArray(x, &sol);

    int r = 0;
    for (int j = 1; j < p->ny+1; j++){
        for (int i = 1; i < p->nx+1; i++){
            phi[j][i] = sol[r];
            r += 1;
        }
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
void computeLaplacianMatrix(Mat A, int rowStart, int rowEnd, param *p){
    double val;
    int index;
    for (int j = 0; j < p->ny; j++){
        for (int i = 0; i < p->nx; i++){
            index = i + j * p->nx;
            val = 0;
            if (i != 0){
                MatSetValue(A, index, index-1 , 1.0, INSERT_VALUES);
                val += 1.;
            }
            if (i != p->nx-1){
                MatSetValue(A, index, index+1 , 1.0, INSERT_VALUES);
                val += 1.;
            }
            if (j != p->ny-1){
                MatSetValue(A, index, index+p->nx , 1.0, INSERT_VALUES);
                val += 1.;
            }            
            if (j != 0){
                MatSetValue(A, index, index-p->nx , 1.0, INSERT_VALUES);
                val += 1.;
            }
            MatSetValue(A, index, index , -val, INSERT_VALUES);
        }
    }
    
    // fix pressure at lower left corner
    // MatSetValue(A, 0, 0 , 1., INSERT_VALUES);
    // MatSetValue(A, 0, 1 , 0., INSERT_VALUES);
    // MatSetValue(A, 1, 0 , 0., INSERT_VALUES);
    // MatSetValue(A, 0, p->nx , 0., INSERT_VALUES);
    // MatSetValue(A, p->nx, 0 , 0., INSERT_VALUES);

    // lower right corner
    MatSetValue(A, p->nx-1, p->nx-1 , 1., INSERT_VALUES);
    MatSetValue(A, p->nx-1, p->nx-2 , 0., INSERT_VALUES);
    MatSetValue(A, p->nx-2, p->nx-1 , 0., INSERT_VALUES);
    MatSetValue(A, p->nx-1, 2*p->nx-1 , 0., INSERT_VALUES);
    MatSetValue(A, 2*p->nx-1, p->nx-1 , 0., INSERT_VALUES);

}

/*To call during the initialization of your solver, before the begin of the time loop*/
PetscErrorCode initialize_poisson_solver(Poisson_data* data, param *p, double **u_star, double **v_star, double **phi){
    PetscInt rowStart = 0; /*rowStart = 0*/
    PetscInt rowEnd = p->nx * p->ny; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

	int nphi = p->nx * p->ny;

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A,5, NULL); 
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(data->A, rowStart, rowEnd, p);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec,PCLU);
    //KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    //PetscPrintf(PETSC_COMM_WORLD, "Assembly of Matrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
void free_poisson_solver(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
