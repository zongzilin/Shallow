#include<iostream>
#include<tuple>
#include<mpi.h>
#include<cblas.h>
#include<boost/program_options.hpp>
#include<fstream>
#include<cstdlib>

#include"cw_blas.h"
#include"cw_mpi.h"

using namespace std;
ofstream outdata;

int main(int argc, char* argv[]){

    ShallowWater_blas sw;
    ShallowWater_mpi mpi;
    ofstream outdata;

    double dt, T;
    int Nx, Ny, ic;

    // MPI PARAM
    int nprocs;
    int pid;  
    int Left, Right;
    int llx, lly;
    int lhx, lhy;
    int lelg, lelgy;

    sw.promptInput(argc, argv, dt, T, Nx, Ny,ic);

    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid); 

    MPI_Request req[nprocs];
    MPI_Status stats[nprocs];
    MPI_Comm cart_comm;

    tie(Left, Right, cart_comm) = mpi.neighbour_PID(nprocs, pid);

    tie(llx, lhx, lly, lhy, lelg, lelgy) = mpi.domain_decomposition(nprocs, Ny, pid, llx, lhx, lly, lhy); 

    double* lx = new double[lelg];
    double* ly = new double[Ny];
    int offset = 3*lelgy;

    double dx = 1;
    double dy = 1;

    sw.SetParameters(llx, lhx, Ny, dx, dy, lx, ly);

    sw.castSpec(Nx, Ny, dt, T, ic, pid);

    int totlelg = lelg*lelgy;

    double* u = new double[totlelg + 6*lelgy];
    double* v = new double[totlelg + 6*lelgy];   
    double* h = new double[totlelg + 6*lelgy];
    double* h0 = new double[totlelg];

    sw.SetInitialConditions(lelg, lelgy, lx, ly, ic, h0);
          
    cblas_dcopy(lelgy*lelg, h0, 1, h + 3*lelgy, 1);

    mpi.swap_boundary(h, Left, Right, lelg, lelgy, totlelg, cart_comm);

    int klx = 3*Nx ;
    int kux = 3*Nx ;
    int ldax = klx + kux + 1;
    double* ddx = new double[ldax*(lelgy*6 + totlelg)];
    double* ddy = new double[7*(lelgy*6 + totlelg)];

    int kly = lelgy - 1;
    int kuy = lelgy - 1;
    int lday = kly + kuy + 1;

    ddx = sw.gen_deri_x_coeff(dx, Nx, Ny, lelg, lelgy, totlelg);
    ddy = sw.gen_deri_y_coeff(dy, lelg, lelgy, totlelg);

/*     double* dhdx = new double[totlelg + 6*lelgy];
    double* dhdy = new double[totlelg + 6*lelgy];
    sw.deri_x_blas(h, Nx, Ny, lelg, lelgy, totlelg, ddx, dhdx);
    sw.deri_y_blas(h, Nx, Ny, lelg, lelgy, totlelg, ddy, dy, dhdy); */

    int n_t = T/dt;

    for (int i_t = 0; i_t < n_t + 1; ++i_t){

        sw.TimeIntegrate(dx, dy, dt, lelg, lelgy, Nx, Ny, u, v, h, ddx, ddy, Left, Right, cart_comm);

        if (pid == 0){
            cout << "t = " << i_t*dt << endl;
        }        

    }

    sw.final_message(pid);

    sw.write_to_file(Nx, Ny, lelg, lelgy, lx, ly, u, v, h, pid, nprocs);






    MPI_Finalize();
    
    return 0;

}