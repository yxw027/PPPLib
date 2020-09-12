//
// Created by cc on 7/20/20.
//

#include "ReadFiles.h"
#include "Solver.h"
#include "OutSol.h"

INITIALIZE_EASYLOGGINGPP

using namespace PPPLib;

int main(int argc, char** argv){
    string logini_path = SetLogConfPath("");
    int log_level = SetLogLevel(1);
    InitLog(argc,argv,logini_path, log_level);
    tPPPLibConf C;

#if 0
    C.fileC.rover="../data_mgex/faa13350.19o";
    C.fileC.brd="../data_mgex/brdm3350.19p";
    C.fileC.cbias="../data_mgex/CAS0MGXRAP_20193350000_01D_01D_DCB.BSX";
    C.fileC.sol="../data_mgex/faa100870_spp.pos";
    C.gnssC.sample_rate=30.0;
    C.mode_opt=MODE_OPT_KINE_SIM;
#endif

#if 1
    C.fileC.rover="../data_cpt0/cpt00870.19o";
    C.fileC.brd="../data_cpt0/brdm0870.19p";
    C.fileC.cbias="../data_cpt0/CAS0MGXRAP_20190870000_01D_01D_DCB.BSX";
    C.fileC.ref="../data_cpt0/cpt00870.ref";
    C.fileC.sol="../data_cpt0/cpt00870_spp.pos";
    C.gnssC.sample_rate=1.0;
    C.mode_opt=MODE_OPT_KINEMATIC;
#endif

    C.mode=MODE_SPP;
    C.gnssC.nav_sys=SYS_GPS;
    C.gnssC.frq_opt=FRQ_SINGLE;
    C.gnssC.ion_opt=ION_KLB;
    C.gnssC.ele_min=10.0;
    C.gnssC.use_doppler=false;
    C.gnssC.max_pdop=30.0;
    C.gnssC.gnss_frq[SYS_INDEX_GPS][0]=GPS_L1;C.gnssC.gnss_frq[SYS_INDEX_GPS][1]=GPS_L2;C.gnssC.gnss_frq[SYS_INDEX_GPS][2]=GPS_L5;
    C.gnssC.gnss_frq[SYS_INDEX_BDS][0]=BDS_B1I;C.gnssC.gnss_frq[SYS_INDEX_BDS][1]=BDS_B2I;C.gnssC.gnss_frq[SYS_INDEX_BDS][2]=BDS_B3I;
    C.gnssC.gnss_frq[SYS_INDEX_GAL][0]=GAL_E1;C.gnssC.gnss_frq[SYS_INDEX_GAL][1]=GAL_E5a;C.gnssC.gnss_frq[SYS_INDEX_GAL][2]=GAL_E5b;
    C.gnssC.gnss_frq[SYS_INDEX_GLO][0]=GLO_G1;C.gnssC.gnss_frq[SYS_INDEX_GLO][1]=GLO_G2;
    C.gnssC.gnss_frq[SYS_INDEX_QZS][0]=QZS_L1;C.gnssC.gnss_frq[SYS_INDEX_QZS][1]=QZS_L2;C.gnssC.gnss_frq[SYS_INDEX_QZS][2]=QZS_L5;
    C.gnssC.gnss_frq[NSYS][0]=BDS_B1I;C.gnssC.gnss_frq[NSYS][1]=BDS_B3I;C.gnssC.gnss_frq[NSYS][2]=BDS_B1C;
    C.solC.out_head=true;
    C.solC.out_vel=true;
    C.solC.out_err=false;
    C.solC.sol_coord=COORD_XYZ;

    cSolver *solver= nullptr;
    solver=new cSppSolver(C);

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
//        ref.pos<<5.08465761725042e+06,2.67032540511120e+06,-2.76848090433746e+06;
        ref.pos<<-5.24739363972353e+06,-3.07686641833132e+06,-1.91152102536219e+06;
        ref.vel<<0,0,0;
        solver->ref_sols_.push_back(ref);
    }

    solver->SolverProcess(C);

    long t2=clock();
    double t=(double)(t2-t1)/CLOCKS_PER_SEC;
    cout<<t<<endl;
}
