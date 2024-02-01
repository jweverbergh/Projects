#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

FILE* fichier_out;

//BLAS
double dnrm2_(int* N, double* X, int* incx);
void daxpy_(int* N, double* alpha, double* X, int* incx, double* Y, int* incy);
void dscal_(int* N, double* alpha, double* A, int* incx);
void dcopy_(int* N, double* X, int* incx, double* Y, int* incy);
double ddot_(int* N, double* X, int* incx, double* Y, int* incy);

//Data strutures
typedef struct ijk {  	// Contains one triplet (index, index, value)
	size_t i,j;			// Indices
	double k;			// Value
} ijk ;

typedef struct COO {	// Compressed matrix
	size_t nnz; 		// Number of non zeros in the matrix
	size_t alloc; 		// Number of entries that have been allocated
	int compressed; 	// Status : sorted/compressed or not
	ijk *data;			// Content of the matrix
} COO;

//some utilities for the data structures given by the teachers
int compareijk(const void *_a, const void *_b){
	ijk *a = (ijk*)_a;
	ijk *b = (ijk*)_b;
	if (a->i < b->i) return -1;
	if (a->i > b->i) return 1;
	if (a->j < b->j) return -1;
	if (a->j > b->j) return 1;
	return 0;
}

void allocateCOO(COO** coo, size_t initial){
	*coo = (COO*) malloc (sizeof(COO));
	(*coo)-> compressed = 0;
	(*coo)-> nnz = 0;
	(*coo)-> alloc = initial;
	(*coo)-> data = (ijk*) malloc ((*coo)->alloc*sizeof(ijk));
}

int allocateVector(double** x, int size){
	*x = (double*) malloc(size*sizeof(double));
	if ((*x)==NULL){
		return -1;
	}return 0;
}

void freeCOO(COO** coo){
	if (*coo == NULL) return;
	free ((*coo)->data);
	free (*coo);
	*coo = NULL;;
}

void addToCoo(COO *coo, size_t I, size_t J, double K){
	// if max size is reached: we double the size of the memory allocated
	if (K==0){
		return;
	}
	if (coo-> nnz == coo->alloc){
		ijk *temp = (ijk*) malloc (2*coo->alloc*sizeof(ijk));
		for (size_t i=0; i<coo->alloc;i++){
			temp[i].i = coo->data[i].i;
			temp[i].j = coo->data[i].j;
			temp[i].k = coo->data[i].k;
		}
		free (coo->data);
		coo->data = temp;
		coo->alloc *= 2;
	}
	coo-> compressed = 0;
	coo->data[coo->nnz].i =  I;
	coo->data[coo->nnz].j =  J;
	coo->data[coo->nnz++].k =  K;
}

//utilities for computing a line of A, given the index of line, where the last arguments are all the coefficients
void four_terms(COO *coo, size_t index, int size, double a, double b, double c, double d, double e){
	addToCoo(coo, index, index, a);
	addToCoo(coo, index, index-1, b);
	addToCoo(coo, index, index+1, c);
	addToCoo(coo, index, index-size, d);
	addToCoo(coo, index, index+size, e);
}

void on_the_pipe(COO *coo, size_t index, int size, double coeff, double a, double b, double c, double d, double dt){
	addToCoo(coo, index, index, 1/dt);
	addToCoo(coo, index, index, coeff*(a+b));
	addToCoo(coo, index, index, coeff*(c+d));
	addToCoo(coo, index, index, coeff*(d+b));
	addToCoo(coo, index, index, coeff*(c+a));
	addToCoo(coo, index, index+1, -coeff*(a+b));
	addToCoo(coo, index, index-1, -coeff*(c+d));
	addToCoo(coo, index, index+size, -coeff*(d+b));
	addToCoo(coo, index, index-size, -coeff*(c+a));
}

//computing A, the points on the surface are not considered as unknowns
void compute_A(COO *coo, int M, int N, double L, double P, double D_air, double D_gnd){
	//init
	double dt = P/N;
	double dx = L/(2*M);
	double coeff = 1/(2*dx*dx);
	double coeff_air = D_air / (dx*dx);
	double coeff_gnd = D_gnd / (dx*dx);
	size_t index;
	size_t size = 2*M + 1; //here and only here, size is the number of points in a line of the domain

	double coeff_c;

	//complete A for the Neumann conditions
	coeff_c = 1/dt + 3*coeff_gnd;
	//on the left side
	for (int j=1; j<M-1; j++){
		index = j*size;
		four_terms(coo, index, size, coeff_c, 0, -coeff_gnd, -coeff_gnd, -coeff_gnd);
	}
	//on the right side
	for (int j=1; j<M-1; j++){
		index = 2*M + j*size;
		four_terms(coo, index, size, coeff_c, -coeff_gnd, 0, -coeff_gnd, -coeff_gnd);
	}
	//on the lower side
	for (int i=1; i<2*M; i++){
		index = i;
		four_terms(coo, index, size, coeff_c, -coeff_gnd, -coeff_gnd, 0, -coeff_gnd);
	}
	//on the lower left corner
	coeff_c = 1/dt + 2*coeff_gnd;
	index = 0;
	four_terms(coo, index, size, coeff_c, 0, -coeff_gnd, 0, -coeff_gnd);
	//on the lower right corner
	index = 2*M;
	four_terms(coo, index, size, coeff_c, -coeff_gnd, 0, 0, -coeff_gnd);

	//complete A for points just below the surface (Dirichlet condition)
	coeff_c = 1/dt + 4*coeff_gnd;
	for (int i=1; i<2*M; i++){
		index = i + (M-1)*size;
		four_terms(coo, index, size, coeff_c, -coeff_gnd, -coeff_gnd, -coeff_gnd, 0);
	}
	//on the upper left corner
	coeff_c = 1/dt + 3*coeff_gnd;
	index = size*(M-1);
	four_terms(coo, index, size, coeff_c, 0, -coeff_gnd, -coeff_gnd, 0);
	//on the upper right corner
	index = size*M - 1;
	four_terms(coo, index, size, coeff_c, -coeff_gnd, 0, -coeff_gnd, 0);

	//complete A for the points inside the domain, with the pipe
	double x, y;
	for (int i=1; i<2*M; i++){
		for (int j=1; j<M-1; j++){
			index = i + j*size;
			x = i*dx;
			y = j*dx;

			//inside the pipe
			if (x>4.75 && x<5.25 && y>1.75 && y<2.25){
				four_terms(coo, index, size, 1/dt+4*coeff_air, -coeff_air, -coeff_air, -coeff_air, -coeff_air);
			}
			//on the edges of the pipe, without the corners
			else if (x>4.75 && x<5.25 && y==1.75){
				on_the_pipe(coo, index, size, coeff, D_gnd, D_air, D_gnd, D_air, dt);
			}
			else if (x>4.75 && x<5.25 && y==2.25){
				on_the_pipe(coo, index, size, coeff, D_air, D_gnd, D_air, D_gnd, dt);
			}
			else if (y>1.75 && y<2.25  && x==4.75){
				on_the_pipe(coo, index, size, coeff, D_air, D_air, D_gnd, D_gnd, dt);
			}
			else if (y>1.75 && y<2.25 && x==5.25){
				on_the_pipe(coo, index, size, coeff, D_gnd, D_gnd, D_air, D_air, dt);
			}
			//on the corners of the pipe
			else if (x==4.75 && y==1.75){
				on_the_pipe(coo, index, size, coeff, D_gnd, D_air, D_gnd, D_gnd, dt);
			}
			else if (x==5.25 && y==1.75){
				on_the_pipe(coo, index, size, coeff, D_gnd, D_gnd, D_gnd, D_air, dt);
			}
			else if (x==4.75 && y==2.25){
				on_the_pipe(coo, index, size, coeff, D_air, D_gnd, D_gnd, D_gnd, dt);
			}
			else if (x==5.25 && y==2.25){
				on_the_pipe(coo, index, size, coeff, D_gnd, D_gnd, D_air, D_gnd, dt);
			}
			//around the pipe
			else{
				four_terms(coo, index, size, 1/dt+4*coeff_gnd, -coeff_gnd, -coeff_gnd, -coeff_gnd, -coeff_gnd);
			}
		}
	}
}

//sort lexicographically then compress
void sortAndCompress (COO *coo) {
	if (coo-> compressed == 1) return;

	//sort the data
	qsort(coo-> data, coo-> nnz, sizeof(ijk), compareijk);

	//compress by computing the sum of the entries at the same place
	double sum = 0;
	size_t I = coo-> data[0].i;
	size_t J = coo-> data[0].j;
	int index = 0;
	size_t nnz_new = 0;
	ijk *temp = (ijk*) malloc(coo->alloc*sizeof(ijk));
	for (int ind=0;ind<coo->nnz;ind++){  //loop over the data in coo
		if (I!=coo-> data[ind].i || J!=coo-> data[ind].j){
			temp[index].i = I;
			temp[index].j = J;
			temp[index].k = sum;
			sum = 0;
			I = coo-> data[ind].i;
			J = coo-> data[ind].j;
			index+=1;
			nnz_new += 1;
		}
		sum += coo-> data[ind].k;
	}temp[index].i = I;
	temp[index].j = J;
	temp[index].k = sum;
	nnz_new += 1;
	free (coo->data);
	coo->data = temp;
	coo->nnz = nnz_new;

	//update status
	coo-> compressed = 1;
}

//compute y = A * x (where A=coo)
void matrixVectorProduct(COO *coo, double *x, double *y) {
	sortAndCompress(coo);

	double sum=0;
	size_t I = coo-> data[0].i;
	size_t J = coo-> data[0].j;
	int index = 0;
	for (int ind=0;ind<coo->nnz;ind++){ //loop over the data in coo
		if (I!=coo-> data[ind].i){
			y[index] = sum;
			sum = 0;
			I = coo-> data[ind].i;
			index+=1;
		}
		J = coo-> data[ind].j;
		sum += (coo-> data[ind].k) * x[J];
	}y[index] = sum;
}

//binary search to find the element of A at (i,j)
//'stop' is the index in the coo that is certainly after what we are looking for, to stop earlier
//return the index in A if (i,j) has a non-null element
//return -1 otherwise
int binarySearch(COO *coo, int i, int j, int stop){
	int first = 0;
	int last = stop;
	int middle;
	while (first<=last){
		middle = (first+last)/2;
		if (i < (coo-> data[middle].i)){
			last = middle-1;
		}else if (i > (coo-> data[middle].i)){
			first = middle+1;
		}else{
			if (j < (coo-> data[middle].j)){
				last = middle-1;
			}else if (j > (coo-> data[middle].j)){
				first = middle+1;
			}else{
				return middle;
			}
		}
	}return -1;
}

//compute ILU
//coo is A, and ilu will store the ILU0 factorization of A (algo 10.4 from the notes on Moodle)
void compute_ILU(COO *coo, COO *ilu) {
	sortAndCompress(coo);

	//copy of A into ILU
	for (int i=0; i<(coo-> nnz); i++){
		addToCoo(ilu, coo-> data[i].i, coo-> data[i].j, coo-> data[i].k);
	}

	//algo
	int a_kk, a_kj, j;
	size_t ind_line, ind_col; //indices i and k in the algo
	for (int i=0; i< (ilu-> nnz); i++){ //loop over the data in coo
		ind_line = ilu-> data[i].i;
		ind_col = ilu-> data[i].j;

		//skip the first line of the matrix
		if(ind_line==0){
			continue;
		}else if (ind_col<ind_line){
			a_kk = binarySearch(ilu, (int) ind_col, (int) ind_col, i);

			//compute the new entry by dividing by a_kk
			ilu-> data[i].k = (ilu-> data[i].k)/(ilu-> data[a_kk].k);
			j = i+1;

			//for each element that satisfies the previous condition, loop onto the elements
			//that are in the same line in the matrix, just after this element
			while (j != (ilu-> nnz) && (ilu-> data[j].i) == (ilu-> data[i].i)){
				a_kj = binarySearch(ilu, (int) ind_col, (int) (ilu->data[j].j), i);

				//if a_kj is found, then compute the new entry
				if (a_kj != -1){
					ilu->data[j].k = (ilu-> data[j].k) - (ilu-> data[i].k) * (ilu-> data[a_kj].k);
				}j++;
			}
		}
	}
}

//compute IC (algo 1)
//coo is A, and ic will store the incomplete Cholesky factorization of A (algo from wikipedia, consulted on the 12th of December)
void compute_IC(COO *coo, COO *ic){
	sortAndCompress(coo);

	//algo
	int L_ki, L_kj; //same notation as in the algorithm from wikipedia, except that here, the indices
	//are reversed to construct ic in a sorted way (to use binarySearch)
	double sum, L_ii_value;
	size_t ind_line, ind_col;
	for (int i=0;i<coo-> nnz;i++){  //loop over the data in the coo sorted
		ind_line = coo-> data[i].i;
		ind_col = coo-> data[i].j;

		//compute L_ii
		if (ind_line == ind_col){
			sum = 0;
			for (int k=0;k<ind_line;k++){
				L_ki = binarySearch(ic, k, (int) ind_line, (int) ic-> nnz-1);
				if (L_ki != -1){
					sum += (ic-> data[L_ki].k) * (ic-> data[L_ki].k);
				}
			}
			//store L_ii for the next iterations, to compute L_ij later
			if (coo-> data[i].k - sum < 0){
				return;
			}

			L_ii_value = sqrt(coo-> data[i].k - sum);

			//add L_ii to ic
			addToCoo(ic, ind_line, ind_col, L_ii_value);
		}

		//compute L_ij
		else if (ind_col > ind_line){
			sum = 0;
			for (int k=0;k<ind_line;k++){
				L_ki = binarySearch(ic, k, (int) ind_line, (int) ic-> nnz-1);
				L_kj = binarySearch(ic, k, (int) ind_col, (int) ic-> nnz-1);
				if (L_ki != -1 && L_kj != -1){
					sum += (ic-> data[L_ki].k) * (ic-> data[L_kj].k);
				}
			}

			//add L_ij to ic
			addToCoo(ic, ind_line, ind_col, (1/L_ii_value) * (coo-> data[i].k - sum));
		}
	}

	//copy the transpose in ic, and sort ic (we need to copy for solving, otherwise we must have added a new sorting for the coo (by columns and not by rows)
	size_t nnz_ic = ic->nnz;
	for (int i=0;i<nnz_ic;i++){
		ind_line = ic-> data[i].i;
		ind_col = ic-> data[i].j;
		if (ind_line != ind_col){
			addToCoo(ic, ind_col, ind_line, ic-> data[i].k);
		}
	}
	sortAndCompress(ic);
}


//compute IC (this second algo has the same structure as the one for the ILU0, with some adequate changes)
//coo is A, and ic will store the incomplete Cholesky factorization of A
void compute_IC_2(COO *coo, COO *ic) {
	sortAndCompress(coo);

	//copy of A into ILU
	for (int i=0; i<(coo-> nnz); i++){
		addToCoo(ic, coo-> data[i].i, coo-> data[i].j, coo-> data[i].k);
	}

	//algo
	int a_kk, a_jk, j;
	size_t ind_line, ind_col; //indices i and k in the algo
	for (int i=0; i< (ic-> nnz); i++){ //loop over the data in coo
		ind_line = ic-> data[i].i;
		ind_col = ic-> data[i].j;

		//skip the first line of the matrix
		if(ind_line==0){
			continue;
		}else if (ind_col<ind_line){
			a_kk = binarySearch(ic, (int) ind_col, (int) ind_col, i);

			//compute the new entry by dividing by a_kk
			ic-> data[i].k = (ic-> data[i].k)/sqrt((ic-> data[a_kk].k));
			j = i+1;

			//for each element that satisfies the previous condition, loop onto the elements
			//that are in the same line in the matrix, just after this element
			while (j != (ic-> nnz) && (ic-> data[j].i) == (ic-> data[i].i)){
				a_jk = binarySearch(ic, (int) (ic->data[j].j), (int) ind_col, i);

				//if a_kj is found, then compute the new entry
				if (a_jk != -1){
					ic->data[j].k = (ic-> data[j].k) - (ic-> data[i].k) * (ic-> data[a_jk].k);
				}j++;
			}
		}
	}

	int elem;
	for (int i=0;i<(ic->nnz);i++){
		if ((ic-> data[i].i)==(ic-> data[i].j)){
			ic-> data[i].k = sqrt(ic-> data[i].k);
		}else if((ic-> data[i].i)<(ic-> data[i].j)){
			elem = binarySearch(ic, (ic-> data[i].j), (ic-> data[i].i), (int) ic-> nnz-1);
			if (elem!=1){
				ic-> data[i].k = ic-> data[elem].k;
			}
		}
	}
}

//solve LU x = b
void systemSolve(COO *lu, double *b, double *x, int size, int param) {
	//init
	double sum = 0;
	size_t ind_line, ind_col;
	double *inter;
	int okInter = allocateVector(&inter, size);

	if (okInter!=0){
		return;
	}

	//forward substitution
	for (int i= 0; i<=(lu-> nnz-1);i++){ //loop over the data in coo
		ind_line = lu-> data[i].i;
		ind_col = lu-> data[i].j;
		if (ind_line == ind_col){
			if (param==1){
				inter[ind_col] = (b[ind_col] - sum);
			}
			else if (param==2){
				inter[ind_col] = (b[ind_col] - sum)/(lu-> data[i].k);
			}
			sum = 0;
		}else if (ind_col < ind_line){
			sum += inter[ind_col] * (lu-> data[i].k);
		}
	}

	//backward substitution
	sum = 0;
	for (int i= (lu-> nnz-1);i>=0;i--){ //loop over the data in coo
		ind_line = lu-> data[i].i;
		ind_col = lu-> data[i].j;
		if (ind_line == ind_col){
			x[ind_col] = (inter[ind_col] - sum)/(lu-> data[i].k);
			sum = 0;
		}else if (ind_col > ind_line){
			sum += x[ind_col] * (lu-> data[i].k);
		}
	}free(inter);
}

// conjugate gradients for solving Ax = b, algo from the manual (38.1)
// return 0 if OK
// return 1 if max iterations attained
// return 2 if division by 0
int CG(COO *A, double *b, double *x, size_t itermax, double precision, int size) {
	//we suppose that x=0 when given in the arguments
	//declaration and initialisation
	double *p, *prod, alpha, beta, alpha_num, alpha_den, beta_num, beta_den, nrm_0, nrm_n;
	int inc = 1;
	int nb_it = 0;
	int okP = allocateVector(&p, size);
	int okProd = allocateVector(&prod, size);

	if (okP!=0 || okProd!=0){
		return 0;
	}

	//algo
	dcopy_(&size, b, &inc, p, &inc);
	beta_num = ddot_(&size, b, &inc, b, &inc);

	nrm_0 = dnrm2_(&size, b, &inc); //stores the norm of the initial residual
	nrm_n = 10*nrm_0; //stores the norm of the residual at the nth iteration
	if (nrm_0==0){
		return 2;
	}
	while (nrm_n/nrm_0>precision){

		//we could delete the alpha_num variable but it is clearer like this
		alpha_num = beta_num;
		matrixVectorProduct(A, p, prod);
		alpha_den = ddot_(&size, p, &inc, prod, &inc);

		//check division by zero
		if (alpha_num == 0 || alpha_den == 0){
			return 2;
		}
		alpha = alpha_num/alpha_den;

		daxpy_(&size, &alpha, p, &inc, x, &inc);

		alpha = -alpha;
		daxpy_(&size, &alpha, prod, &inc, b, &inc);

		beta_num = ddot_(&size, b, &inc, b, &inc);
		beta_den = alpha_num;
		beta = beta_num/beta_den;

		dscal_(&size, &beta, p, &inc);
		alpha = 1;
		daxpy_(&size, &alpha, b, &inc, p, &inc);

		nrm_n = dnrm2_(&size, b, &inc);

		//check the number of iterations
		if (nb_it == itermax){
			return 1;
		}
		nb_it+=1;
	}
	free(p);
	free(prod);
	return 0;
}

// preconditioned conjugate gradients for solving Ax = b, PREC ~= A^{-1}, algo from the course notes
// return 0 if OK
// return 1 if max iterations attained
// return 2 if division by 0
int PCG(COO *A, COO *prec, double *b, double *x, size_t itermax, double precision, int size, int param) {
	//we suppose that x=0 when given in the arguments
	//declaration and initialisation
	double *p, *r_tilde, *prod, alpha, beta, alpha_num, alpha_den, beta_num, beta_den, nrm_0, nrm_n;
	int inc = 1;
	int nb_it = 0;
	int okP = allocateVector(&p, size);
	int okProd = allocateVector(&prod, size);
	int okRTilde = allocateVector(&r_tilde, size);

	if (okP!=0 || okProd!=0 || okRTilde!=0){
		return 0;
	}

	//algo
	systemSolve(prec, b, r_tilde, size, param);
	dcopy_(&size, r_tilde, &inc, p, &inc);
	beta_num = ddot_(&size, b, &inc, r_tilde, &inc);

	nrm_0 = dnrm2_(&size, b, &inc); //stores the norm of the initial residual
	nrm_n = 10*nrm_0; //stores the norm of the residual at the nth iteration
	if (nrm_0==0){
		return 2;
	}
	while (nrm_n/nrm_0>precision){

		alpha_num = beta_num;
		matrixVectorProduct(A, p, prod);
		alpha_den = ddot_(&size, p, &inc, prod, &inc);

		//check division by zero
		if (alpha_num == 0 || alpha_den == 0){
			return 2;
		}
		alpha = alpha_num/alpha_den;

		daxpy_(&size, &alpha, p, &inc, x, &inc);

		alpha = -alpha;
		daxpy_(&size, &alpha, prod, &inc, b, &inc);

		systemSolve(prec, b, r_tilde, size, param);

		beta_num = ddot_(&size, b, &inc, r_tilde, &inc);
		beta_den = alpha_num;
		beta = beta_num/beta_den;

		dscal_(&size, &beta, p, &inc);
		alpha = 1;
		daxpy_(&size, &alpha, r_tilde, &inc, p, &inc);

		nrm_n = dnrm2_(&size, b, &inc);

		//check the number of iterations
		if (nb_it == itermax){
			return 1;
		}
		nb_it+=1;
	}
	free(prod);
	free(r_tilde);
	free(p);
	return 0;
}

//compute the function f of the statement
double f(double t, double T_max, double T_min, double P){
	double first = (T_max + T_min)/2;
	double second = (T_max - T_min) * (sin(2*M_PI*t/P))/2;
	return first + second;
}

//reset x to 0 (could be avoided if we do the first iteration of GC with the hand, before the loop
void reset(double *x, int size){
	for (int i=0;i<size;i++){
		x[i] = 0.;
	}
}

//solve the heat equation with iter_time iterations and beginning with T_sol, stores the solution in T_sol
//param = 0 if no conditionning
//param = 1 if conditionning by ILU0
//param = 2 if conditionning by IC
int solveHeatEq(COO *coo, int iter_time, double *T_sol, int param, int size, double precision, size_t itermax, int M, double P, int N, double T_max, double T_min, double coeff_gnd, int center, double *T_center){
	//init
	double dt = P/N;
	int check;
	int inc = 1;
	double cst = 1/dt;
	double eta;
	double mini = 50;
	double maxi = 0;

	COO *ilu, *ic;
	double *T_init;
	int ok = allocateVector(&T_init, size);
	if (ok!=0){
		return 0;
	}
	if (param == 1){
		allocateCOO(&ilu, coo-> alloc);
		compute_ILU(coo, ilu);
	}
	if (param == 2){
		allocateCOO(&ic, coo-> alloc);
		compute_IC(coo, ic);  //or compute_IC_2(coo, ic) to use the other IC implementation, both give the same results
	}

	//solve for iter_time iterations
	for (int i=1;i<=iter_time;i++){
		//copy T_sol into T_init
		dcopy_(&size, T_sol, &inc, T_init, &inc);

		//scan T_init in the center of the pipe, in the file temp.out
		fprintf(fichier_out, "%lf %.16lf\n", (i-1)*dt, T_init[center]);

		//updates the min and the max of T in the center of the pipe
		if (i > iter_time-N && T_init[center] > maxi){ //the first condition imposes that we are on the last periode
			maxi = T_init[center];
		}
		if (i > iter_time-N && T_init[center] < mini){
			mini = T_init[center];
		}

		//reset T_sol to 0
		reset(T_sol, size);

		//divide T_init by dt
		dscal_(&size, &cst, T_init, &inc);

		//complete the part of T_init for the Dirichlet condition
		for (int j=(2*M+1)*(M-1); j<=size-1; j++){
			T_init[j] += coeff_gnd*f((i)*dt, T_max, T_min, P);
		}

		//choose the conditionning and solve
		if (param == 0){
			check = CG(coo, T_init, T_sol, itermax, precision, size);
		}else if (param == 1){
			check = PCG(coo, ilu, T_init, T_sol, itermax, precision, size, param);
		}else if (param == 2){
			check = PCG(coo, ic, T_init, T_sol, itermax, precision, size, param);
		}

		//check if no division by zero and itermax not reached
		if (check == 1){
			return 1;
		}if (check == 2){
			return 2;
		}
	}

	//compute eta and save the value in T_center
	eta = (maxi-mini)/(T_max-T_min);
	T_center[0] = eta;

	//free
	free(T_init);
	if (param == 1){
		freeCOO(&ilu);
	}if (param == 2){
		freeCOO(&ic);
	}return 0;
}

int main(int argc, char *argv[]) {
	//check the arguments
	if (argc!=7){
		return 0;
	}

	//open the file
	if ((fichier_out = fopen(argv[6], "w")) == NULL){
		return 0;
	}

	//parameters
	int DorY = atoi(argv[1]);
	int nb = atoi(argv[2]);

	double L = 10;
	double H = 5;
	double P;
	double T_max;
	double T_min;
	//day-night
	if (DorY == 1){
		P = 24*60*60;
		T_min = 8;
		T_max = 20;
	}
	//summer-winter
	if (DorY == 0){
		P = 24*60*60*365;
		T_min = -5;
		T_max = 25;
	}

	int M = atoi(argv[3]);
	int N = atoi(argv[4]);
	double D_air = 20*pow(10,-6);
	double D_gnd = 0.25*pow(10,-6);
	int size = (2*M+1)*(M);
	double dt = P/N;
	double dx = L/(2*M);
	double coeff_gnd = D_gnd/(dx*dx);
	int center = (5/dx) + (2/dx)*(2*M+1);

	double precision = 1e-12;
	int param = atoi(argv[5]);
	int iter_time = nb*N + 1;
	size_t itermax = 200;
	int check;

	COO *coo;
	double *T_sol, *T_center;

	allocateCOO(&coo, 5);
	int ok = allocateVector(&T_sol, size);
	int ok_center = allocateVector(&T_center, 1);

	if (ok!=0 || ok_center!=0){
		return 0;
	}

	//compute A
	compute_A(coo, M, N, L, P, D_air, D_gnd);

	//compute the initial solution
	for (int i=0;i<size;i++){
		T_sol[i] = f(0, T_max, T_min, P);
	}

	//write a line to keep this line for the last step (with false values but the good length)
	fprintf(fichier_out, "%d %.16lf\n", iter_time, T_center[0]);

	//solve the heat equation
	check = solveHeatEq(coo, iter_time, T_sol, param, size, precision, itermax, M, P, N, T_max, T_min, coeff_gnd, center, T_center);

	//check division by zero and number of iterations reached
	if (check == 1 || check == 2){
		return 0;
	}

	//last step : print the iter_time and eta at the begining of the file
	rewind(fichier_out);
	fprintf(fichier_out, "%d %.16lf\n", iter_time, T_center[0]);

	//close the file and free
	fclose(fichier_out);
	freeCOO(&coo);
	free(T_sol);
	free(T_center);
	return 0;
}
