#ifndef cw 
#define cw 

#include<cblas.h>
#include<cmath>
#include<mpi.h>
#include<iostream>
#include<boost/program_options.hpp>
#include<cstdlib>

using namespace std;
namespace po = boost::program_options;

class ShallowWater{

public:

    void castSpec(int Nx, int Ny, double dt, double T, double ic, int pid){

        if (pid == 0){
            cout << "------------ SOLVER START ------------" << endl;
            cout << "                                      " << endl;
            cout << "----------- SPECIFICATIONS -----------" << endl;
            cout << " Nx: " << Nx << "                      " << " Ny: " << Ny << endl;
            cout << " dt: " << dt << "                      " << " T: " << T << endl;
            cout << " Initial Conditions: " << ic << endl;
            cout << "--------------------------------------" << endl;
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

            dphidy[index] = invdy*((-1.0/60)*phi[index + (lelgy - 4)] 
                                        + (3.0/20)*phi[index + (lelgy - 3)] 
                                        - (3.0/4)*phi[index + (lelgy - 2)] 
                                        + (3.0/4)*phi[index + 1]  
                                        - (3.0/20)*phi[index + 2] 
                                        + (1.0/60)*phi[index + 3]);
            j = 2;
            index = index + 1;

            dphidy[index] = invdy*((-1.0/60)*phi[index + (lelgy - 5)] 
                                        + (3.0/20)*phi[index + (lelgy - 4)] 
                                        - (3.0/4)*phi[index + (lelgy - 3)] 
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
            
            index = index + 1;
            dphidy[index] =  invdy*((-1.0/60)*phi[index - 3] 
                                    + (3.0/20)*phi[index - 2] 
                                    - (3.0/4)*phi[index - 1] 
                                    + (3.0/4)*phi[index + 1]  
                                    - (3.0/20)*phi[index - (lelgy - 3)] 
                                    + (1.0/60)*phi[index - (lelgy - 4)]);   

            index = index + 1;     
            dphidy[index] =  invdy*((-1.0/60)*phi[index - 3] 
                                    + (3.0/20)*phi[index - 2] 
                                    - (3.0/4)*phi[index - 1] 
                                    + (3.0/4)*phi[index - (lelgy - 3)]  
                                    - (3.0/20)*phi[index - (lelgy - 4)] 
                                    + (1.0/60)*phi[index - (lelgy - 5)]);  

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
                       int Left, int Right, int nprocs){
        
        const double g = 9.81;

        MPI_Request req[nprocs];
        MPI_Status stats[nprocs];        

        double rk_co[4] = {1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0};
        double rk_tco[4] = {1.0, 1.0/2.0, 1.0/2.0, 1.0};

      
    }

};

#endif 