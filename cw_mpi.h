#ifndef cw_mpi
#define cw_mpi

#include<cstdlib>
#include<mpi.h>
#include<tuple>

using namespace std;

class ShallowWater_mpi{

    public:

    tuple<int, int, MPI_Comm> neighbour_PID(int nprocs, int pid){

        int dims[2] = {0,0};

        MPI_Dims_create(nprocs, 1, dims);

        int periods[2] = {true, true};
        int reorder = true;
        MPI_Comm cart_comm;
        MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, reorder, &cart_comm);
        
        int* neighbour_rank = new int[4];

        MPI_Cart_shift(cart_comm, 1, 1, &neighbour_rank[3], &neighbour_rank[4]);

        int d = 0;
        int u = 1;
        int L = neighbour_rank[3];
        int R = neighbour_rank[4];

        return make_tuple(L, R, cart_comm);

    }

    tuple<int, int, int, int, int, int> domain_decomposition(int nprocs, int n, int pid, 
                            int& start, int& end, int lly, int lhy){
    int r = n % nprocs;
    int k = (n - r)/nprocs;

    if (pid < (n % nprocs)){
        k++;
        start = k*pid;
        end = k*(pid + 1);
    }
    else{
        start = (k + 1)*r + k*(pid - r);
        end = (k + 1)*r + k*(pid - r + 1);
    }

    int lelg = end - start;

    lly = 0;
    lhy = n;

    int lelgy = lhy - lly;

    return make_tuple(start, end, lly, lhy, lelg, lelgy);

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



};


#endif