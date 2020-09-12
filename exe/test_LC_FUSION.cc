//
// Created by cc on 8/2/20.
//

#include "ReadFiles.h"
#include "DecodeRaw.h"
#include "Solver.h"
INITIALIZE_EASYLOGGINGPP

using namespace PPPLib;

int main(int argc, char** argv){
    string logini_path = SetLogConfPath("");
    int log_level = SetLogLevel(32);
    InitLog(argc,argv,logini_path, log_level);
    tPPPLibConf C;

#if 1
    // GNSS-GSOF AND INS LOOSELY COUPLED
    C.fileC.gsof="/home/cc/dataset/data_carvig_lc/m390311.gsof";
    C.fileC.imu="/home/cc/dataset/data_carvig_lc/m390311.imu";
    C.fileC.sol="/home/cc/dataset/data_carvig_lc/gsof_lc.pos";
    C.mode=MODE_IGLC;
    C.mode_opt=MODE_OPT_GSOF;
    C.insC.imu_type=IMU_M39;
    C.insC.coord_type=IMU_COORD_RFU;
    C.insC.data_format=IMU_FORMAT_RATE;
    C.insC.gyro_val_format=GYRO_FORMAT_DEG;
    C.insC.sample_rate=200.0;
    C.insC.psd_acce=2.33611111111E-07;
    C.insC.psd_gyro=5.72003802085E-09;
    C.insC.psd_ba=1E-07;
    C.insC.psd_bg=2E-12;
    C.insC.lever[0]=-0.818;
    C.insC.lever[1]=-0.010;
    C.insC.lever[2]=-0.010;
#endif

#if 0
    // GNSS-RTK AND INS LOOSELY COUPLED
    C.fileC.rover="../data_cpt0/cpt00870.19o";
    C.fileC.base="../data_cpt0/cpt00870_base.19o";
    C.fileC.brd="../data_cpt0/brdm0870.19p";
    C.fileC.cbias="../data_cpt0/CAS0MGXRAP_20190870000_01D_01D_DCB.BSX";
    C.fileC.sol="../data_cpt0/cpt00870_ppk_lc.pos";
    C.mode=MODE_IGLC;
    C.mode_opt=MODE_OPT_PPK;
    C.insC.imu_type=IMU_NOVTEL_CPT;
    C.insC.coord_type=IMU_COORD_RFU;
    C.insC.data_format=IMU_FORMAT_RATE;
    C.insC.sample_rate=100.0;
#endif

    C.insC.init_pos_unc=30.0;
    C.insC.init_vel_unc=30.0;
    C.insC.init_att_unc=0.174532922;
    C.insC.init_ba_unc=9.80665E-3;
    C.insC.init_bg_unc=4.8481367284E-05;
    C.insC.err_model= false;
    C.insC.correction_time_ba=360.0;
    C.insC.correction_time_bg=360.0;
    C.solC.out_head=true;
    C.solC.out_vel=true;
    C.solC.out_att=true;
    C.solC.out_ba=true;
    C.solC.out_bg=true;
    C.solC.out_ins_mech_frq=100.0;
    C.solC.sol_coord=COORD_XYZ;

    cSolver *solver= nullptr;
    solver=new cFusionSolver(C);

    long t1=clock();

    solver->SolverProcess(C);

    long t2=clock();
    double t=(double)(t2-t1)/CLOCKS_PER_SEC;
    cout<<t<<endl;
}

