//
// Created by cc on 7/20/20.
//

#include "ReadFiles.h"
#include "Solver.h"
#include "PlotFunc.h"
INITIALIZE_EASYLOGGINGPP

using namespace PPPLib;

#define USE_CPT 1
#define USE_M39 0
#define USE_CU  0

int main(int argc, char** argv){
    string logini_path = SetLogConfPath("");
    int log_level = SetLogLevel(1);
    InitLog(argc,argv,logini_path, log_level);
    tPPPLibConf C;

#if USE_CPT
    C.fileC.rover="/home/cc/dataset/data_cpt0/cpt00870.19o";
    C.fileC.base="/home/cc/dataset/data_cpt0/cpt00870_base.19o";
    C.fileC.brd="/home/cc/dataset/data_cpt0/brdm0870.19p";
    C.fileC.cbias="/home/cc/dataset/data_cpt0/CAS0MGXRAP_20190870000_01D_01D_DCB.BSX";
    C.fileC.ref="/home/cc/dataset/data_cpt0/cpt00870.ref";
    C.fileC.sol="/home/cc/dataset/data_cpt0/ppplib_ppk1.pos";
    C.fileC.sol_stat="/home/cc/dataset/data_cpt0/ppplib_ppk1.pos.stat";
    C.gnssC.adj_obs=false;
    C.gnssC.check_dual_phase=false;
    C.gnssC.sample_rate=1.0;
    C.gnssC.gnss_frq[SYS_INDEX_BDS][0]=BDS_B1I;C.gnssC.gnss_frq[SYS_INDEX_BDS][1]=BDS_B2I;C.gnssC.gnss_frq[SYS_INDEX_BDS][2]=BDS_B3I;
    C.gnssC.gnss_frq[NSYS][0]=BDS_B1I;C.gnssC.gnss_frq[NSYS][1]=BDS_B3I;C.gnssC.gnss_frq[NSYS][2]=BDS_B1C;
#endif

#if USE_M39
    C.fileC.rover="/home/cc/dataset/data_carvig_tc/m390311.17o";
    C.fileC.base="/home/cc/dataset/data_carvig_tc/m390311_base.17o";
    C.fileC.brd="/home/cc/dataset/data_carvig_tc/brdm0311.17p";
    C.fileC.cbias="/home/cc/dataset/data_carvig_tc/CAS0MGXRAP_20173110000_01D_01D_DCB.BSX";
    C.fileC.sol="/home/cc/dataset/data_carvig_tc/ppplib_ppk.pos";
    C.fileC.sol_stat="/home/cc/dataset/data_carvig_tc/ppplib_ppk.pos.stat";
    C.gnssC.sample_rate=1.0;
    C.gnssC.adj_obs=true;
    C.gnssC.check_dual_phase=true;
    C.gnssC.gnss_frq[SYS_INDEX_BDS][0]=BDS_B1I;C.gnssC.gnss_frq[SYS_INDEX_BDS][1]=BDS_B2I;C.gnssC.gnss_frq[SYS_INDEX_BDS][2]=BDS_B3I;
    C.gnssC.gnss_frq[NSYS][0]=BDS_B1I;C.gnssC.gnss_frq[NSYS][1]=BDS_B3I;C.gnssC.gnss_frq[NSYS][2]=BDS_B1C;
#endif

#if USE_CU
    C.fileC.rover="/home/cc/dataset/data_cu/cubb0660.20o";
    C.fileC.base="/home/cc/dataset/data_cu/cuaa0660.20o";
    C.fileC.brd="/home/cc/dataset/data_cu/brdm0660.20p";
    C.fileC.cbias="/home/cc/dataset/data_cu/CAS0MGXRAP_20200660000_01D_01D_DCB.BSX";
    C.fileC.sol="/home/cc/dataset/data_cu/ppplib_ppk.pos";
    C.fileC.sol_stat="/home/cc/dataset/data_cu/ppplib_ppk.pos.stat";
    C.gnssC.sample_rate=30.0;
    C.gnssC.adj_obs=false;
    C.gnssC.gnss_frq[SYS_INDEX_BDS][0]=BDS_B1I;C.gnssC.gnss_frq[SYS_INDEX_BDS][1]=BDS_B2I;C.gnssC.gnss_frq[SYS_INDEX_BDS][2]=BDS_B3I;
    C.gnssC.gnss_frq[NSYS][0]=BDS_B1I;C.gnssC.gnss_frq[NSYS][1]=BDS_B3I;C.gnssC.gnss_frq[NSYS][2]=BDS_B1C;

    C.gnssC.check_dual_phase= false;

#endif

    C.mode=MODE_PPK;
    C.mode_opt=MODE_OPT_KINEMATIC;
    C.gnssC.nav_sys=SYS_GPS;
    C.gnssC.frq_opt=FRQ_DUAL;
    C.gnssC.ion_opt=ION_OFF;
    C.gnssC.trp_opt=TRP_SAAS;
    C.gnssC.use_doppler=false;
    C.gnssC.ele_min=15.0;
    C.gnssC.max_pdop=30.0;
    C.gnssC.max_out=20;
    C.gnssC.max_inno=30.0;
    C.gnssC.cs_thres[0]=5.0;
    C.gnssC.cs_thres[1]=0.05;
    C.gnssC.ait_psd[0]=0.0;
    C.gnssC.ar_mode=AR_CONT;
    C.gnssC.ar_thres[0]=3.0;
    C.gnssC.ar_thres[1]=0.25;
    C.gnssC.ar_thres[2]=0.0;
    C.gnssC.ar_thres[3]=1E-9;
    C.gnssC.ar_thres[4]=1E-5;
    C.gnssC.ar_el_mask=10.0;
    C.gnssC.min_sat_num2fix=3;
    C.gnssC.min_sat_num2drop=8;
    C.gnssC.min_lock2fix=0;
    C.gnssC.partial_ar=true;
    C.gnssC.gnss_frq[SYS_INDEX_GPS][0]=GPS_L1;C.gnssC.gnss_frq[SYS_INDEX_GPS][1]=GPS_L2;C.gnssC.gnss_frq[SYS_INDEX_GPS][2]=GPS_L5;
    C.gnssC.gnss_frq[SYS_INDEX_GAL][0]=GAL_E1;C.gnssC.gnss_frq[SYS_INDEX_GAL][1]=GAL_E5a;C.gnssC.gnss_frq[SYS_INDEX_GAL][2]=GAL_E5b;
    C.gnssC.gnss_frq[SYS_INDEX_GLO][0]=GLO_G1;C.gnssC.gnss_frq[SYS_INDEX_GLO][1]=GLO_G2;
    C.gnssC.gnss_frq[SYS_INDEX_QZS][0]=QZS_L1;C.gnssC.gnss_frq[SYS_INDEX_QZS][1]=QZS_L2;C.gnssC.gnss_frq[SYS_INDEX_QZS][2]=QZS_L5;
    C.gnssC.res_qc=1;
    C.solC.out_head=true;
    C.solC.out_err= false;
    C.solC.out_vel=false;
    C.solC.out_stat=true;
    C.solC.sol_coord=COORD_XYZ;

    cSolver *solver= nullptr;
    solver=new cPpkSolver(C);

    long t1=clock();

    if(!C.fileC.rover.empty()){
        cReadGnssObs rover_reader(C.fileC.rover,solver->nav_,solver->rover_obs_,REC_ROVER);
        rover_reader.SetGnssSysMask(C.gnssC.nav_sys);
        rover_reader.Reading();
    }

    if(!C.fileC.brd.empty()){
        cReadGnssBrdEph brd_reader(C.fileC.brd, solver->nav_);
        brd_reader.SetGnssSysMask(C.gnssC.nav_sys);
        brd_reader.Reading();
    }

    if(!C.fileC.cbias.empty()){
        cReadGnssCodeBias cbias_reader(C.fileC.cbias, solver->nav_);
        cbias_reader.Reading();
    }

    cReadRefSol *ref_reader= nullptr;
    if(C.mode_opt==MODE_OPT_KINEMATIC||C.mode>=MODE_IGLC){
        if(!C.fileC.ref.empty()&&C.solC.out_err){

            ref_reader=new cReadRefSol(C.fileC.ref,solver->ref_sols_);
            ref_reader->Reading(1);
        }
    }
    else if(C.mode_opt==MODE_OPT_KINE_SIM){
        tSolInfoUnit ref;
        ref.pos<<5.08465761725042e+06,2.67032540511120e+06,-2.76848090433746e+06;
        ref.vel<<0,0,0;
        solver->ref_sols_.push_back(ref);
    }

    solver->SolverProcess(C);

    long t2=clock();
    double t=(double)(t2-t1)/CLOCKS_PER_SEC;
    cout<<t<<endl;
}


