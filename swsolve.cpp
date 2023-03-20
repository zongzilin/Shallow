#include<iostream>
#include<tuple>
#include<mpi.h>
#include<cblas.h>
#include<boost/program_options.hpp>
#include<fstream>
#include<cstdlib>

#include"cw.h"
#include"cw_mpi.h"

using namespace std;
using std::ofstream;


int main(int argc, char* argv[]){

    ShallowWater sw;
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

    const double L = 100.0;
    double dx = L/Nx;
    double dy = L/Ny;

    sw.SetParameters(llx, lhx, Ny, dx, dy, lx, ly);

    sw.castSpec(Nx, Ny, dt, T, ic, pid);

    int totlelg = lelg*lelgy;
    double* h0 = new double[totlelg];
    double* h = new double[totlelg + 6*lelgy];

    sw.SetInitialConditions(lelg, lelgy, lx, ly, ic, h0);
          
    cblas_dcopy(lelgy*lelg, h0, 1, h + 3*lelgy, 1);

    double* u = new double[totlelg + 6*lelgy];
    double* v = new double[totlelg + 6*lelgy];   

    double* u_rk = new double[totlelg + 6*lelgy];
    double* v_rk = new double[totlelg + 6*lelgy];    
    double* h_rk = new double[totlelg + 6*lelgy];

    double* dudx = new double[totlelg + 6*lelgy];
    double* dudy = new double[totlelg + 6*lelgy];
    double* dvdx = new double[totlelg + 6*lelgy];
    double* dvdy = new double[totlelg + 6*lelgy];

    double* dhdx = new double[totlelg + 6*lelgy]; 
    double* dhdy = new double[totlelg + 6*lelgy];

    double* dudt = new double[totlelg + 6*lelgy];
    double* dvdt = new double[totlelg + 6*lelgy];
    double* dhdt = new double[totlelg + 6*lelgy];

    int n_t = T/dt;
    double rk_co[4] = {1.0/6, 2.0/6, 2.0/6, 1.0/6};
    double rk_tco[4] = {1.0, 1.0/2, 1.0/2, 1.0};

    for (int i_t = 0; i_t < n_t + 1; ++i_t){

        cblas_dcopy(totlelg + offset, u, 1, u_rk, 1);
        cblas_dcopy(totlelg + offset, v, 1, v_rk, 1);
        cblas_dcopy(totlelg + offset, h, 1, h_rk, 1);

        for (int r_k = 0; r_k < 4; ++r_k){

            MPI_Isend(v_rk + offset, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[0]);
            MPI_Irecv(v_rk, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[1]);
            MPI_Isend(v_rk + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[2]);
            MPI_Irecv(v_rk + offset + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[3]);

            MPI_Waitall(4, req, stats);

            MPI_Isend(u_rk + offset, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[0]);
            MPI_Irecv(u_rk, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[1]);
            MPI_Isend(u_rk + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[2]);
            MPI_Irecv(u_rk + offset + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[3]);

            MPI_Waitall(4, req, stats);                

            MPI_Isend(h_rk + offset, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[0]);
            MPI_Irecv(h_rk, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[1]);
            MPI_Isend(h_rk + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[2]);
            MPI_Irecv(h_rk + offset + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[3]);

            MPI_Waitall(4, req, stats);

            sw.deri_x(u_rk, dx, lelg, lelgy, dudx);
            sw.deri_y(u_rk, dy, lelg, lelgy, dudy);

            sw.deri_x(v_rk, dx, lelg, lelgy, dvdx);
            sw.deri_y(v_rk, dy, lelg, lelgy, dvdy);

            sw.deri_x(h_rk, dx, lelg, lelgy, dhdx);
            sw.deri_y(h_rk, dy, lelg, lelgy, dhdy);

            for (int i = 0; i < lelg; ++i){
                for (int j = 0; j < lelgy; ++j){
                    int ind = i*lelgy + j + offset;

                    dudt[ind] = -u[ind]*dudx[ind] - v[ind]*dudy[ind] - 9.810*dhdx[ind];
                    dvdt[ind] = -u[ind]*dvdx[ind] - v[ind]*dudy[ind] - 9.810*dhdy[ind];
                    dhdt[ind] = -h[ind]*dudx[ind] - u[ind]*dhdx[ind] - h[ind]*dvdy[ind] - v[ind]*dhdy[ind];

                }
            }

            for(int i = 0; i < lelg; ++i){
                for(int j = 0; j < lelgy; ++j){
                    int ind = i*lelgy + j + offset;
                    
                    u_rk[ind] = u[ind] + dt*rk_co[r_k]*rk_tco[r_k]*dudt[ind];
                    v_rk[ind] = v[ind] + dt*rk_co[r_k]*rk_tco[r_k]*dvdt[ind];
                    h_rk[ind] = h[ind] + dt*rk_co[r_k]*rk_tco[r_k]*dhdt[ind];
                }
            }
        }

        cblas_dcopy(totlelg + offset, u_rk, 1, u, 1);
        cblas_dcopy(totlelg + offset, v_rk, 1, v, 1);
        cblas_dcopy(totlelg + offset, h_rk, 1, h, 1);

        if (pid == 0){
            cout << i_t*dt << endl;
        }
        
    }

















    double* hout = new double[totlelg];
    double* FFF = new double[Nx*Ny];

    cblas_dcopy(lelgy*lelg, h + offset, 1, hout, 1);
    mpi.gather_to_root(hout, nprocs, lelg, lelgy, totlelg, lx, ly, 0, FFF);

    if (pid == 0){
        outdata.open("u.dat");
        for(int i = 0; i < Nx*Ny; ++i){
            outdata << FFF[i] << endl;
        }
        outdata.close();
    }



    MPI_Finalize();
    
    return 0;

}