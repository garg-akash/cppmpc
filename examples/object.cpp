#include <iostream>
#include <fstream>
#include <math.h>
#include "mpc_common.hpp"
#include "rti_step.hpp"
#include "casadi_wrapper.hpp"

using namespace std;
Eigen::MatrixXd skewMat(const Eigen::VectorXd &v){
    Eigen::Matrix3d s;
    s << 0, -v[2], v[1],
        v[2], 0, -v[0],
        -v[1], v[0], 0;
    return s;
}

Eigen::Matrix3d zyx2R(const Eigen::VectorXd &v){ //euler to rotation matrix
    float c_x = cos(v[2]);
    float s_x = sin(v[2]);
    float c_y = cos(v[1]);
    float s_y = sin(v[1]);
    float c_z = cos(v[0]);
    float s_z = sin(v[0]);
  
    Eigen::Matrix3d rot_z, rot_y, rot_x, R;
    rot_z << c_z, -s_z, 0,
            s_z, c_z, 0,
            0, 0, 1;

    rot_y << c_y, 0, s_y,
            0, 1, 0,
            -s_y, 0, c_y;

    rot_x << 1, 0, 0,
            0, c_x, -s_x,
            0, s_x, c_x;

    R << rot_z*rot_y*rot_x;
    return R;       
}

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
    int i, j;         
    for(i=0;i<size.nbx;i++)
        size.nbx_idx[i]=i;
    
    // create a RTI controller
    rti_step_workspace rti_work(size);

    // initial condition and parameters
    // read ref from txt file
    int nrows = 251; int ncols = 45;
    MatrixXd ref = MatrixXd::Zero(nrows,ncols);
    ifstream fin ("./ref.txt");

    if (fin.is_open())
    {
        for (int row = 0; row < nrows; row++)
            for (int col = 0; col < ncols; col++)
            {
                float item = 0.0;
                fin >> item;
                ref(row, col) = item;
            }
        fin.close();
    }
    else
        cout << "ref file not opened" << endl;

    /*ofstream myfile2;
    myfile2.open ("ref2.txt");
    myfile2 << ref.transpose() << endl;*/
    for(i=0;i<size.nx;i++)
        rti_work.x0(i) = ref(0,i);
    for(i=0;i<size.nx;i++){
        for(j=0;j<size.N+1;j++){
            rti_work.QP.in.x(i,j) = ref(0,i);
        }
    }

    // ref
    for(i=0;i<size.ny;i++){
        for(j=0;j<size.N;j++){
            rti_work.QP.in.y(i,j) = ref(j,i);
        }
    } 
    for(i=0;i<size.nyN;i++){
        rti_work.QP.in.yN(i) = ref(size.N,i);
    } 

    // parameters
    MatrixXd M_ = MatrixXd::Zero(6,6);
    M_.block(0,0,3,3) = 0.5*MatrixXd::Identity(3,3);
    M_.block(3,3,3,3) = 1e-4*MatrixXd::Identity(3,3);
    MatrixXd C_ = MatrixXd::Zero(6,6);
    C_.block(0,0,3,3) = skewMat(rti_work.x0.segment(9,3))*M_(0,0);
    C_.block(3,3,3,3) = skewMat(rti_work.x0.segment(9,3))*M_.block(3,3,3,3);
    VectorXd N_(6), gVec(3);
    gVec << 0,0,9.8;
    N_ << M_(0,0)*zyx2R(rti_work.x0.segment(3,3))*gVec, 0.0, 0.0, 0.0;
    Eigen::VectorXd pVEC(size.np);
    pVEC << (Eigen::Map<Eigen::VectorXd>(M_.data(), M_.cols()*M_.rows())),
            (Eigen::Map<Eigen::VectorXd>(C_.data(), C_.cols()*C_.rows())),
            N_;
    for(i=0;i<size.np;i++){
        for(j=0;j<size.N;j++){
            rti_work.QP.in.p(i,j) = pVEC(i);
        }
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
    int n_iter = 0;
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

        // update parameters
        C_.block(0,0,3,3) = skewMat(rti_work.x0.segment(9,3))*M_(0,0);
        C_.block(3,3,3,3) = skewMat(rti_work.x0.segment(9,3))*M_.block(3,3,3,3);
        N_ << M_(0,0)*zyx2R(rti_work.x0.segment(3,3))*gVec, 0.0, 0.0, 0.0;
        pVEC << (Eigen::Map<Eigen::VectorXd>(M_.data(), M_.cols()*M_.rows())),
                (Eigen::Map<Eigen::VectorXd>(C_.data(), C_.cols()*C_.rows())),
                N_;
        for(i=0;i<size.np;i++){
            for(j=0;j<size.N;j++){
                rti_work.QP.in.p(i,j) = pVEC(i);
            }
        }

        // update ref
        for(i=0;i<size.ny;i++){
            for(j=0;j<size.N;j++){
                if(j+n_iter+1 > 251)
                    rti_work.QP.in.y(i,j) = ref(251,i);
                else
                    rti_work.QP.in.y(i,j) = ref(j+n_iter+1,i);
            }
        } 
        for(i=0;i<size.nyN;i++){
            if(n_iter+1+size.N > 251)
                rti_work.QP.in.yN(i) = ref(251,i);
            else
                rti_work.QP.in.yN(i) = ref(n_iter+size.N+1,i);
        }
        n_iter++;
    }

    // free memory
    myfile.close();
      
    delete [] size.nbx_idx;
    size.nbx_idx = NULL;
    
    rti_work.free();

    return 0;
}