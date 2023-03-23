#ifndef cw 
#define cw 

#include<cblas.h>
#include<mpi.h>
#include<cmath>
#include<iostream>
#include<boost/program_options.hpp>
#include<cstdlib>
#include<fstream>

#include"cw_mpi.h"

using namespace std;
using std::ofstream;
namespace po = boost::program_options;


class ShallowWater{

    ofstream outdata;

public:

    void castSpec(int Nx, int Ny, double dt, double T, double ic, int pid){

        if (pid == 0){
            cout << "----------- SPECIFICATIONS -----------" << endl;
            cout << " Nx: " << Nx << "                      " << " Ny: " << Ny << endl;
            cout << " dt: " << dt << "                      " << " T: " << T << endl;
            cout << " Initial Conditions: " << ic << endl;
            cout << "--------------------------------------" << endl;
            cout << "Timestepping ...." << endl;
        }   

    }

    void promptInput(int argc, char* argv[], \
                    double &dt,  double &T, \
                    int &Nx,  int &Ny, \
                    int &ic){
                        
        po::options_description opts(
            "Reads in time stepping, grid and initial conditions info for AERO70011 coursework"
        );

        opts.add_options()
            ("dt", po::value<double>() -> default_value(0.1),
            "time-step size.")
            ("T", po::value<double>() -> default_value(100.0),
            "Total integration time.")
            ("Nx", po::value<int>() -> default_value(100),
            "Number of grid points in x")
            ("Ny", po::value<int>() -> default_value(100),
            "Number of grid points in y")
            ("ic", po::value<int>() -> default_value(1),
            "Index of initial condition to use (1-4)")
            ("help", "Print Help Message");

            po::variables_map vm;
            po::store(po::parse_command_line(argc, argv, opts),vm);

            po::notify(vm);

            if(vm.count("help")){
                cout << opts << endl;
            }

            dt = vm["dt"].as<double>();
            T = vm["T"].as<double>();
            Nx = vm["Nx"].as<int>();
            Ny = vm["Ny"].as<int>();
            ic = vm["ic"].as<int>();
    }

    void SetParameters(int ll, int lh, int ny, double dx, double dy, double* x, double* y){

        for (int i = ll; i < lh; ++i){
            x[i - ll] = i*dx;
        }
        for (int j = 0; j < ny; ++j){
            y[j] = j*dy;
        }
        
    }

    void SetInitialConditions(int lelg, int lelgy, 
                              double* x, double* y,
                              int scen, double* G){

        // X AND Y ARE DECOMPOSED DOMAIN COORDINATES
        // scen IS INITIAL CONDITIONS SCENARIOS
        
        switch (scen){
            case 1:
                
                for (int i = 0; i < lelg; ++i){
                    for (int j = 0; j < lelgy; ++j){
                        G[i*lelgy + j] = 10 + exp(-(x[i] - 50)*(x[i]-50)/25);
                    }
                }
            break;

            case 2:

                for (int i = 0; i < lelg; ++i){
                    for (int j = 0; j < lelgy; ++j){
                        G[i*lelgy + j] = 10 + exp(-(y[j] - 50)*(y[j] - 50)/25);
                    }
                }
            break;

            case 3:

                for (int i = 0; i < lelg; ++i){
                    for (int j = 0; j < lelgy; ++j){
                        G[i*lelgy + j] = 10 + exp(-((x[i] - 50)*(x[i] - 50) + \
                                         (y[j] - 50)*(y[j] - 50))/25);
                    }
                }
            break;

            case 4:

                for (int i = 0; i < lelg; ++i){
                    for (int j = 0; j < lelgy; ++j){
                        G[i*lelgy + j] = 10 + exp(-((x[i] - 25)*(x[i] - 25) + 
                                         (y[j] - 25)*(y[j] - 25))/25) + 
                                         exp(-((x[i] - 75)*(x[i] - 75) + 
                                         (y[j] - 75)*(y[j] - 75))/25);
                    }
                }
            break;
        }

    }

    void deri_y(double* phi, double dy, int lelg, int lelgy, double* dphidy){

        // SIZE OF PHI MUST BE: 3*lelgy + lelgy*lelg + 3*lelgy
        // FIRST 3*lelgy needs to be recieved from left pid
        // LAST 3*lelgy needs to be passed from right pid
        // NEEDS TO PASS A TOTAL OF 6*lelgy from lelgy*lelg to left and right pid
        
        double invdy = 1.0/dy;

        // NEED TO INFLATE phi
/*         double* phi_ghost = new double[lelgy + 6];
        double* dphidy_ghost = new double[lelgy + 6];

        for (int i = 0; i < lelg; ++i){

            int shift = i*lelgy;
            
            cblas_dcopy(3, phi + shift, 1, phi_ghost + lelgy + 3, 1);
            cblas_dcopy(lelgy, phi + shift, 1, phi_ghost + 3, 1);
            cblas_dcopy(3, phi + shift + lelgy - 3, 1, phi_ghost, 1);
            
            for (int j = 3; j < lelgy + 3; ++j){
                
                dphidy_ghost[j] = invdy*((-1.0/60)*phi_ghost[j - 3] 
                                            + (3.0/20)*phi_ghost[j - 2] 
                                            - (3.0/4)*phi_ghost[j - 1] 
                                            + (3.0/4)*phi_ghost[j + 1]  
                                            - (3.0/20)*phi_ghost[j + 2] 
                                            + (1.0/60)*phi_ghost[j + 3]); 

            }

            cblas_dcopy(lelgy, dphidy_ghost + 3, 1, dphidy + shift, 1);
        }
 */

        for(int i = 0; i < lelg; ++i){ // no. of col
            int j = 0; // first row

            int index = i*lelgy + 3*lelgy;

            dphidy[index] = invdy*((-1.0/60)*phi[index +(lelgy - 3)] 
                                        + (3.0/20)*phi[index +(lelgy - 2)] 
                                        - (3.0/4)*phi[index +(lelgy - 1)] 
                                        + (3.0/4)*phi[index + 1]  
                                        - (3.0/20)*phi[index + 2] 
                                        + (1.0/60)*phi[index + 3]);

            j = 1; // second row
            index = i*lelgy + 3*lelgy + j;

            dphidy[index] = invdy*((-1.0/60)*phi[index + (lelgy - 3)] 
                                        + (3.0/20)*phi[index + (lelgy - 2)] 
                                        - (3.0/4)*phi[index - 1] 
                                        + (3.0/4)*phi[index + 1]  
                                        - (3.0/20)*phi[index + 2] 
                                        + (1.0/60)*phi[index + 3]);
            j = 2;
            index = index + 1;

            dphidy[index] = invdy*((-1.0/60)*phi[index + (lelgy - 3)] 
                                        + (3.0/20)*phi[index - 2] 
                                        - (3.0/4)*phi[index - 1] 
                                        + (3.0/4)*phi[index + 1]  
                                        - (3.0/20)*phi[index + 2] 
                                        + (1.0/60)*phi[index + 3]);

        }

        for(int i = 0; i < lelg; ++i){
            for (int j = 3; j < lelgy - 3; ++j){

                int index = i*lelgy + j + 3*lelgy;

                dphidy[index] = invdy*((-1.0/60)*phi[index - 3] 
                                            + (3.0/20)*phi[index - 2] 
                                            - (3.0/4)*phi[index - 1] 
                                            + (3.0/4)*phi[index + 1]  
                                            - (3.0/20)*phi[index + 2] 
                                            + (1.0/60)*phi[index + 3]);
            
            }
        } 

        for(int i = 0; i < lelg; ++i){
            
            int j = lelgy - 3; // third last row
            int index = i*lelgy + 3*lelgy + j;

            dphidy[index] =  invdy*((-1.0/60)*phi[index - 3] 
                                        + (3.0/20)*phi[index - 2] 
                                        - (3.0/4)*phi[index - 1] 
                                        + (3.0/4)*phi[index + 1]  
                                        - (3.0/20)*phi[index + 2] 
                                        + (1.0/60)*phi[index - (lelgy - 3)]);
            
            j = lelgy - 2;
            index = i*lelgy + 3*lelgy + j;
            dphidy[index] =  invdy*((-1.0/60)*phi[index - 3] 
                                    + (3.0/20)*phi[index - 2] 
                                    - (3.0/4)*phi[index - 1] 
                                    + (3.0/4)*phi[index + 1]  
                                    - (3.0/20)*phi[index  - (lelgy - 2)] 
                                    + (1.0/60)*phi[index - (lelgy - 3)]);   

            j = lelgy - 1;    
            index = i*lelgy + 3*lelgy + j; 
            dphidy[index] =  invdy*((-1.0/60)*phi[index - 3] 
                                    + (3.0/20)*phi[index - 2] 
                                    - (3.0/4)*phi[index - 1] 
                                    + (3.0/4)*phi[index - (lelgy - 1)]  
                                    - (3.0/20)*phi[index - (lelgy - 2)] 
                                    + (1.0/60)*phi[index - (lelgy - 3)]);  

        }

    }

    void deri_x(double* phi, double dx, int lelg, int lelgy, double* dphidx){

        double invdx = 1.0/dx;

        for(int i = 0; i < lelg; ++i){
            for (int j = 0; j < lelgy; ++j){

                int index = i*lelgy + j + 3*lelgy;

                dphidx[index] = invdx*((-1.0/60)*phi[index - 3*lelgy]
                                             + (3.0/20)*phi[index - 2*lelgy]
                                             - (3.0/4.0)*phi[index - lelgy]
                                             + (3.0/4.0)*phi[index + lelgy]
                                             - (3.0/20.0)*phi[index + 2*lelgy]
                                             + (1.0/60.0)*phi[index + 3*lelgy]);
            }
        }

    }

    void TimeIntegrate(double dx, double dy, double dt, int lelg, int lelgy,
                       double* u, double* v, double* h,
                       int Left, int Right, 
                       MPI_Comm cart_comm){ 
                     
        
        int totlelg = lelgy*lelg;
        int offset = 3*lelgy;

        double rk_co[4] = {1.0/6, 2.0/6, 2.0/6, 1.0/6};
        double rk_tco[4] = {1.0, 1.0/2, 1.0/2, 1.0};  

        double* dudx = new double[totlelg + 6*lelgy];
        double* dudy = new double[totlelg + 6*lelgy];
        double* dvdx = new double[totlelg + 6*lelgy];
        double* dvdy = new double[totlelg + 6*lelgy];      

        double* dhdx = new double[totlelg + 6*lelgy];
        double* dhdy = new double[totlelg + 6*lelgy];    

        double* u_rk = new double[totlelg + 6*lelgy];
        double* v_rk = new double[totlelg + 6*lelgy];               
        double* h_rk = new double[totlelg + 6*lelgy]; 

        double* dudt = new double[totlelg + 6*lelgy];
        double* dvdt = new double[totlelg + 6*lelgy];               
        double* dhdt = new double[totlelg + 6*lelgy];  

        cblas_dcopy(totlelg + offset*2, u, 1, u_rk, 1);
        cblas_dcopy(totlelg + offset*2, v, 1, v_rk, 1);
        cblas_dcopy(totlelg + offset*2, h, 1, h_rk, 1);         

        for (int r_k = 0; r_k < 4; ++r_k){

            swap_boundary(u_rk,Left, Right, lelg, lelgy, totlelg, cart_comm);
            swap_boundary(v_rk,Left, Right, lelg, lelgy, totlelg, cart_comm);
            swap_boundary(h_rk,Left, Right, lelg, lelgy, totlelg, cart_comm);

            deri_x(u_rk, dx, lelg, lelgy, dudx);
            deri_y(u_rk, dy, lelg, lelgy, dudy);

            deri_x(v_rk, dx, lelg, lelgy, dvdx);
            deri_y(v_rk, dy, lelg, lelgy, dvdy);

            deri_x(h_rk, dx, lelg, lelgy, dhdx);
            deri_y(h_rk, dy, lelg, lelgy, dhdy);

            for (int i = 0; i < lelg; ++i){
                for (int j = 0; j < lelgy; ++j){
                    int ind = i*lelgy + j + offset;

                    dudt[ind] = -u[ind]*dudx[ind] - v[ind]*dudy[ind] - 9.810*dhdx[ind];
                    dvdt[ind] = -u[ind]*dvdx[ind] - v[ind]*dudy[ind] - 9.810*dhdy[ind];
                    dhdt[ind] = -h[ind]*dudx[ind] - u[ind]*dhdx[ind] - h[ind]*dvdy[ind] - v[ind]*dhdy[ind];

                    u_rk[ind] = u[ind] + dt*rk_co[r_k]*rk_tco[r_k]*dudt[ind];
                    v_rk[ind] = v[ind] + dt*rk_co[r_k]*rk_tco[r_k]*dvdt[ind];
                    h_rk[ind] = h[ind] + dt*rk_co[r_k]*rk_tco[r_k]*dhdt[ind];

                }
            }
        }

        cblas_dcopy(totlelg + offset*2, u_rk, 1, u, 1);
        cblas_dcopy(totlelg + offset*2, v_rk, 1, v, 1);
        cblas_dcopy(totlelg + offset*2, h_rk, 1, h, 1);
      
    }

    void write_to_file(int Nx, int Ny, int lelg, int lelgy, double* lx, double* ly,
                        double*u, double* v, double* h, int pid, int nprocs){

    int totlelg = lelg*lelgy;
    int offset = 3*lelgy;

    double* x = new double[Ny];
    double* y = new double[Ny];
    double* X = new double[Nx*Ny];
    double* Y = new double[Nx*Ny];


    if (pid == 0){

        for (int i = 0; i < Ny; ++i){
            x[i] = i;
            for (int j = 0; j < Nx; ++j){
                Y[i*Ny + j] = x[i];
            }
        }


        for (int i = 0; i < lelg; ++i){
            cblas_dcopy(Ny, x, 1, X + i*lelgy, 1);
        }
    }


    double* uout = new double[totlelg];
    double* vout = new double[totlelg];
    double* hout = new double[totlelg];

    double* FFF = new double[Nx*Ny];
    double* GGG = new double[Nx*Ny];
    double* HHH = new double[Nx*Ny];

    cblas_dcopy(lelgy*lelg, h + offset, 1, hout, 1);
    cblas_dcopy(lelgy*lelg, v + offset, 1, vout, 1);
    cblas_dcopy(lelgy*lelg, u + offset, 1, uout, 1);

    gather_to_root(uout, nprocs, lelg, lelgy, totlelg, lx, ly, 0, FFF);
    gather_to_root(vout, nprocs, lelg, lelgy, totlelg, lx, ly, 0, GGG);
    gather_to_root(hout, nprocs, lelg, lelgy, totlelg, lx, ly, 0, HHH);

    if (pid == 0){
        outdata.open("result.dat");
        for(int i = 0; i < Nx*Ny; ++i){
            outdata << X[i] << " " << Y[i] << " " << FFF[i] << " " << GGG[i] << " " << HHH[i] << endl;
        }
        outdata.close();
    }

    }

    void swap_boundary(double* swapbuf, int Left, int Right, int lelg, int lelgy, int totlelg, MPI_Comm cart_comm){

        int offset = 3*lelgy;

        MPI_Request req[4];
        MPI_Status stats[4];

        MPI_Isend(swapbuf + offset, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[0]);
        MPI_Irecv(swapbuf, offset, MPI_DOUBLE, Left, 1, cart_comm, &req[1]);
        MPI_Isend(swapbuf + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[2]);
        MPI_Irecv(swapbuf + offset + totlelg, offset, MPI_DOUBLE, Right, 1, cart_comm, &req[3]);

        MPI_Waitall(3, req, stats);             

    }

    void gather_to_root(double* sendbuf, int nprocs, int lelg, int lelgy, 
                        int totlelg, double* lx, double* ly,
                        int root, double* revcbuf){

        int* revc = new int[nprocs];
        int xstart = lx[0]*lelgy;
        int* disp = new int[nprocs];

        // GATHER DISPLACEMENT AND REVCIEVE COUNT
        MPI_Gather(&totlelg, 1, MPI_INT, revc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&xstart, 1, MPI_INT, disp, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Gatherv(sendbuf, lelg*lelgy, MPI_DOUBLE, revcbuf, revc, disp, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }    

    void final_message(int pid, double tt){
        if (pid == 0){
            cout << "  TOTAL ELASPED TIME:    " << tt << endl;
            cout << "  RUN FINISHED. THIS IS SOLVED USING LOOP-BASED ITERATIONS. " << endl;
        }
    }

};

#endif 