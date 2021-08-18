#include <iostream>
#include <fstream>
#include <math.h>
#include "mpc_common.hpp"
#include "rti_step.hpp"
#include "casadi_wrapper.hpp"

using namespace std;
int main()
{  

    // define problem size
    model_size size;
    size.nx = 29;
    size.nu = 16;
    size.ny = 45;
    size.nyN= 29;
    size.np = 78;
    size.nbu = 16;
    size.nbx = 29;
    size.nbg = 0;
    size.nbgN = 0;
    size.N = 40;
    size.nbx_idx = new int[size.nbx];   
    int i;         
    for(i=0;i<size.nbx;i++)
        size.nbx_idx[i]=i;
    
    // create a RTI controller
    rti_step_workspace rti_work(size);

    // initial condition and parameters
    rti_work.x0(1) = M_PI;

    for(i=0;i<size.N+1;i++){
        rti_work.QP.in.x(1,i) = M_PI;
    }

    // ref
    for(i=0;i<size.N;i++){
        rti_work.QP.in.y(0,i) = 0.3;
    } 

    // parameters
    for(i=0;i<size.N+1;i++){
        rti_work.QP.in.p(0,i) = 0.3;
    } 

    // weight matrix
    for(i=0;i<6;i++){
        rti_work.QP.in.W(i,i) = 10;
        rti_work.QP.in.W(i+6,i+6) = 1e-2;

        rti_work.QP.in.WN(i,i) = 10;
        rti_work.QP.in.WN(i+6,i+6) = 1e-2;
    }
    rti_work.QP.in.W(12,12) = 0;
    rti_work.QP.in.WN(12,12) = 0;
    for(i=0;i<16;i++){
        rti_work.QP.in.W(i+13,i+13) = 1e-4;

        rti_work.QP.in.WN(i+13,i+13) = 1e-4;
    }
    for(i=0;i<16;i++)
        rti_work.QP.in.W(i+29,i+29) = 1e-4;

    for(i=0;i<16;i++){
        rti_work.QP.in.lbu(i) = -10000; 
        rti_work.QP.in.ubu(i) = 10000;
    }

    rti_work.QP.in.lbx(0) = -2; 
    rti_work.QP.in.lbx(1) = -2;
    rti_work.QP.in.lbx(2) = -2;

    rti_work.QP.in.lbx(3) = -M_PI;
    rti_work.QP.in.lbx(4) = -M_PI/2;
    rti_work.QP.in.lbx(5) = -M_PI;

    rti_work.QP.in.lbx(6) = -2; 
    rti_work.QP.in.lbx(7) = -2;
    rti_work.QP.in.lbx(8) = -2;

    rti_work.QP.in.lbx(9) = -M_PI;
    rti_work.QP.in.lbx(10) = -M_PI;
    rti_work.QP.in.lbx(11) = -M_PI;

    rti_work.QP.in.lbx(12) = 9.8;

    for(i=0;i<16;i++)
        rti_work.QP.in.lbx(i+13) = 0;

    rti_work.QP.in.ubx(0) = 2; 
    rti_work.QP.in.ubx(1) = 2;
    rti_work.QP.in.ubx(2) = 2;

    rti_work.QP.in.ubx(3) = M_PI;
    rti_work.QP.in.ubx(4) = M_PI/2;
    rti_work.QP.in.ubx(5) = M_PI;

    rti_work.QP.in.ubx(6) = 2; 
    rti_work.QP.in.ubx(7) = 2;
    rti_work.QP.in.ubx(8) = 2;

    rti_work.QP.in.ubx(9) = M_PI;
    rti_work.QP.in.ubx(10) = M_PI;
    rti_work.QP.in.ubx(11) = M_PI;

    rti_work.QP.in.ubx(12) = 9.8;

    for(i=0;i<16;i++)
        rti_work.QP.in.ubx(i+13) = 20;  

    rti_work.QP.in.reg = 1E-8;

    // prepare the closed-loop simulation
    rti_work.sample = 0;
    double Tf=2.5, Ts=0.025,t=0;

    double *simu_in[3];
    double *simu_out[1];
    simu_in[2] = rti_work.QP.in.p.col(0).data();

    ofstream myfile;
    myfile.open ("data.txt");

    // start the simulation
    while(t<Tf){
        
        // call RTI solving routine
        rti_work.step();

        // feedback
        simu_in[0] = rti_work.QP.in.x.col(0).data();
        simu_in[1] = rti_work.QP.in.u.col(0).data();
        simu_out[0] = rti_work.x0.data();
        F_Fun(simu_in, simu_out);

        // update sampling and time
        rti_work.sample++;
        t += Ts;

        // store the closed-loop results
        // myfile <<"Sample " << rti_work.sample <<": " << rti_work.x0.transpose() << " | " << rti_work.QP.in.u.col(0) <<" |CPT=: " << rti_work.CPT << "ms" <<endl;
        myfile << rti_work.x0.transpose() << endl;

        // shifting(optional)
        for(i=0;i<size.N-1;i++){
            rti_work.QP.in.x.col(i) = rti_work.QP.in.x.col(i+1);
            rti_work.QP.in.u.col(i) = rti_work.QP.in.u.col(i+1);
        }
        rti_work.QP.in.x.col(size.N-1) = rti_work.QP.in.x.col(size.N);   

    }

    // free memory
    myfile.close();
      
    delete [] size.nbx_idx;
    size.nbx_idx = NULL;
    
    rti_work.free();

    return 0;
}