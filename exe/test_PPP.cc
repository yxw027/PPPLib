//
// Created by cc on 8/2/20.
//

#include "ReadFiles.h"
#include "Solver.h"
#include "PlotFunc.h"
INITIALIZE_EASYLOGGINGPP

using namespace PPPLib;

int main(int argc, char** argv){
    string logini_path = SetLogConfPath("");
    int log_level = SetLogLevel(1);
    InitLog(argc,argv,logini_path, log_level);

    tPPPLibConf C;

    // PPP with mgex data
#if 0
    C.mode_opt=MODE_OPT_KINE_SIM;
    C.fileC.rover="/home/cc/dataset/data_mgex/faa13350.19o";
    C.fileC.brd="/home/cc/dataset/data_mgex/brdm3350.19p";
    C.fileC.cbias="/home/cc/dataset/data_mgex/CAS0MGXRAP_20193350000_01D_01D_DCB.BSX";
    C.fileC.sp3[1]="/home/cc/dataset/data_mgex/wum20820.sp3";
    C.fileC.sp3[2]="/home/cc/dataset/data_mgex/wum20821.sp3";
    C.fileC.clk="/home/cc/dataset/data_mgex/wum20820.clk";
    C.fileC.atx="/home/cc/dataset/data_mgex/igs14_2097.atx";
    C.fileC.blq="/home/cc/dataset/data_mgex/ocnload.blq";
    C.fileC.erp="/home/cc/dataset/data_mgex/wum20820.erp";
    C.fileC.sol="/home/cc/dataset/data_mgex/ppplib_ppp.pos";
    C.gnssC.sample_rate=30.0;
    C.gnssC.tid_opt=TID_SOLID;
    C.gnssC.ait_psd[2]=10E-4;
#endif

    // PPP-AR with grg products
#if 0
    C.mode_opt=MODE_OPT_KINE_SIM;
    C.fileC.rover="/home/cc/dataset/data_ppp_ar/harb3350.19o";
    C.fileC.brd="/home/cc/dataset/data_ppp_ar/brdm3350.19p";
    C.fileC.cbias="/home/cc/dataset/data_ppp_ar/CAS0MGXRAP_20193350000_01D_01D_DCB.BSX";
    C.fileC.sp3[1]="/home/cc/dataset/data_ppp_ar/grg20820.sp3";
    C.fileC.sp3[2]="/home/cc/dataset/data_ppp_ar/grg20821.sp3";
    C.fileC.clk="/home/cc/dataset/data_ppp_ar/grg20820.clk";
    C.fileC.atx="/home/cc/dataset/data_ppp_ar/igs14_2097.atx";
    C.fileC.blq="/home/cc/dataset/data_ppp_ar/ocnload.blq";
    C.fileC.erp="/home/cc/dataset/data_ppp_ar/igs19P2082.erp";
    C.fileC.sol="/home/cc/dataset/data_ppp_ar/ppplib_ppp_ar.pos";
    C.gnssC.sample_rate=30.0;
    C.gnssC.tid_opt=TID_SOLID;
    C.gnssC.ait_psd[2]=10E-4;
#endif

    //PPP with real kinematic data
#if 1
    C.mode_opt=MODE_OPT_KINEMATIC;
    C.fileC.rover="/home/cc/dataset/data_cpt0/cpt00870.19o";
    C.fileC.brd="/home/cc/dataset/data_cpt0/brdm0870.19p";
    C.fileC.cbias="/home/cc/dataset/data_cpt0/CAS0MGXRAP_20190870000_01D_01D_DCB.BSX";
    C.fileC.ref="/home/cc/dataset/data_cpt0/cpt00870.ref";
    C.fileC.sp3[0]="/home/cc/dataset/data_cpt0/wum20463.sp3";
    C.fileC.sp3[1]="/home/cc/dataset/data_cpt0/wum20464.sp3";
    C.fileC.sp3[2]="/home/cc/dataset/data_cpt0/wum20465.sp3";
    C.fileC.clk="/home/cc/dataset/data_cpt0/wum20464.clk";
    C.fileC.erp="/home/cc/dataset/data_cpt0/wum20464.erp";
    C.fileC.atx="/home/cc/dataset/data_cpt0/igs14_2097.atx";
    C.fileC.blq="/home/cc/dataset/data_cpt0/ocnload.blq";
    C.fileC.sol="/home/cc/dataset/data_cpt0/ppplib_ppp_igg.pos";
    C.gnssC.sample_rate=1.0;
    C.gnssC.tid_opt=TID_OFF;
    C.gnssC.ait_psd[2]=10E-4;
#endif

    C.mode=MODE_PPP;
    C.dynamic=false;
    C.gnssC.nav_sys=SYS_GPS;
    C.gnssC.frq_opt=FRQ_DUAL;
    C.gnssC.ion_opt=ION_IF;
    C.gnssC.trp_opt=TRP_EST_WET;
    C.gnssC.eph_opt=EPH_PRE;
    C.gnssC.glo_ifcb_opt=GLO_IFCB_QUAD;
    C.gnssC.ele_min=10.0;
    C.gnssC.code_phase_ratio=300.0;
    C.gnssC.meas_err_factor[0]=100.0;
    C.gnssC.meas_err_factor[1]=0.003;
    C.gnssC.meas_err_factor[2]=0.003;
    C.gnssC.meas_err_factor[2]=0.0;
    C.gnssC.max_pdop=30.0;
    C.gnssC.max_out=20;
    C.gnssC.max_prior=30.0;
    C.gnssC.cs_thres[0]=5.0;
    C.gnssC.cs_thres[1]=0.05;
    C.gnssC.ar_mode=AR_PPP_AR_ILS;
    C.gnssC.ar_thres[0]=3.0;
    C.gnssC.ar_thres[1]=0.25;
    C.gnssC.ar_thres[2]=0.0;
    C.gnssC.ar_thres[3]=1E-9;
    C.gnssC.ar_thres[4]=1E-5;
    C.gnssC.gnss_frq[SYS_INDEX_GPS][0]=GPS_L1;C.gnssC.gnss_frq[SYS_INDEX_GPS][1]=GPS_L2;C.gnssC.gnss_frq[SYS_INDEX_GPS][2]=GPS_L5;
    C.gnssC.gnss_frq[SYS_INDEX_BDS][0]=BDS_B1I;C.gnssC.gnss_frq[SYS_INDEX_BDS][1]=BDS_B2I;C.gnssC.gnss_frq[SYS_INDEX_BDS][2]=BDS_B3I;
    C.gnssC.gnss_frq[SYS_INDEX_GAL][0]=GAL_E1;C.gnssC.gnss_frq[SYS_INDEX_GAL][1]=GAL_E5a;C.gnssC.gnss_frq[SYS_INDEX_GAL][2]=GAL_E5b;
    C.gnssC.gnss_frq[SYS_INDEX_GLO][0]=GLO_G1;C.gnssC.gnss_frq[SYS_INDEX_GLO][1]=GLO_G2;
    C.gnssC.gnss_frq[SYS_INDEX_QZS][0]=QZS_L1;C.gnssC.gnss_frq[SYS_INDEX_QZS][1]=QZS_L2;C.gnssC.gnss_frq[SYS_INDEX_QZS][2]=QZS_L5;
    C.gnssC.gnss_frq[NSYS][0]=BDS_B1I;C.gnssC.gnss_frq[NSYS][1]=BDS_B3I;C.gnssC.gnss_frq[NSYS][2]=BDS_B1C;
    C.gnssC.res_qc=0;
    C.solC.out_head=true;
    C.solC.out_vel=true;
    C.solC.out_err=false;
    C.solC.sol_coord=COORD_XYZ;

    cSolver *solver= nullptr;
    solver=new cPppSolver(C);

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
            ref_reader->Reading();
        }
    }
    else if(C.mode_opt==MODE_OPT_KINE_SIM||C.mode_opt==MODE_OPT_STATIC){
        tSolInfoUnit ref;
//        ref.pos<<5.08465761725042e+06,2.67032540511120e+06,-2.76848090433746e+06; //HARB
//        ref.pos<<-2.27982908063864e+06,5.00470646117120e+06,3.21977738845034e+06; //JFNG
        ref.pos<<-5.24739363972353e+06,-3.07686641833132e+06,-1.91152102536219e+06;//FAA1
        Vector3d blh=Xyz2Blh(ref.pos);
        fprintf(stderr,"%10.7f %10.7f %10.7f \n\n",blh[0]*R2D,blh[1]*R2D,blh[2]);
        ref.vel<<0,0,0;
        solver->ref_sols_.push_back(ref);
    }

    solver->SolverProcess(C);

    long t2=clock();
    double t=(double)(t2-t1)/CLOCKS_PER_SEC;
    cout<<t<<endl;
}
