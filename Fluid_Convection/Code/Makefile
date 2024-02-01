
INC_DIR := -I/home/user/.../Fluid_Convection/Code/lib_petsc/include
LIB_DIR := -L/home/user/.../Fluid_Convection/Code/lib_petsc/lib -Wl,-rpath=/home/user/.../Fluid_Convection/Code/lib_petsc/lib

LIB := -lpetsc

CXX_FLAGS := -O0 -Wall -Werror #-g

#Compilation
all :
	gcc -o project fluid_convection.c poisson.c -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

#Delete of executable files
clean :
	rm project

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
