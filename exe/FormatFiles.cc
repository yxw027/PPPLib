//
// Created by cc on 9/5/20.
//

/* The Source code formats imu/ref_sol file to PPPLib standard format */

#include "ReadFiles.h"
#include "DecodeRaw.h"
#include "OutSol.h"


INITIALIZE_EASYLOGGINGPP
using namespace PPPLib;

cImuData imus;

/* turn your reference solution to PPPLib(rtklib) standard format */
void FormatRefSol(tPPPLibConf C,string ref_file,int type){
    vector<tSolInfoUnit> ref_sols;

    cReadRefSol ref_reader(ref_file,ref_sols);

    cOutSol out;
    string out_path=ref_file+"f";
    C.solC.out_stat=false;
    C.solC.out_head=true;
    C.solC.out_ba=true;
    C.solC.out_bg=true;
    C.solC.out_att=true;
    C.solC.out_vel=true;

    if(type){
        ref_reader.Reading(1);
        out.InitOutSol(C,out_path);
        out.WriteHead();
        for(int i=0;i<ref_sols.size();i++){
            if(i%100==0){
                out.WriteSol(ref_sols[i],0);
            }
        }
    }
    else{
        ref_reader.SetDataIdx(1,3,4,5,6,7,8,9,10,11,',',12);
        ref_reader.Reading(0);
        out.InitOutSol(C,out_path);
        out.WriteHead();
        for(int i=0;i<ref_sols.size();i++){
            out.WriteSol(ref_sols[i],0);
        }
    }

    ref_sols.clear();
}

/* turn your imu measurement to PPPLib standard format ------------------------------------------------
 * GPS_week GPS_Wos Acce_x(m/s^2) Acce_y(m/s^2) Acce_z(m/s^2) Gyro_x(deg/s) Gyro_y(deg/s) Gyro_z(deg/s)
 * Note: local navigation frame is NED and body frame is FRD
 * ----------------------------------------------------------------------------------------------------*/
void FormatImu(tPPPLibConf C,string imu_file){
    if(C.insC.imu_type==IMU_M39){
        /* m39 t_tag is second of week, need to add GPS week according gnss data tag */
        cDecodeImuM39 m39_decoder;
        C.insC.coord_type=IMU_COORD_RFU;
        C.insC.data_format=IMU_FORMAT_RATE;
        C.insC.gyro_val_format=GYRO_FORMAT_DEG;
        m39_decoder.DecodeM39(imu_file,C.insC,imus.data_);
    }
    else if(C.insC.imu_type==IMU_NOVTEL_CPT){
        C.insC.coord_type=IMU_COORD_RFU;
        C.insC.data_format=IMU_FORMAT_INCR;
        C.insC.gyro_val_format=GYRO_FORMAT_RAD;
        cReadImu imu_reader(imu_file);
        imu_reader.SetImu(C);
        imu_reader.Reading();
        imus=*imu_reader.GetImus();
    }
    else if(C.insC.imu_type==IMU_MTI_CSV){
        cReadImu imu_reader(imu_file);
        imu_reader.SetImu(C);
        imu_reader.SetDataIdx(1,3,4,5,6,7,8,',',9);
        imu_reader.Reading();
        imus=*imu_reader.GetImus();
#if 1
        C.insC.coord_type=IMU_COORD_RFU;
        C.insC.gyro_val_format=GYRO_FORMAT_RAD;
        C.insC.data_format=IMU_FORMAT_RATE;
        for(int i=0;i<imus.data_.size();i++){
            AdjustImuData(imus.data_[i],C.insC.coord_type,C.insC.data_format,C.insC.gyro_val_format,1.0/C.insC.sample_rate);
        }
#endif
    }

    cOutSol out;
    string out_path=imu_file+"f";
    C.solC.out_stat=false;
    C.solC.out_head=true;
    out.InitOutSol(C,out_path);
    out.WriteImuHead();
    out.imus=&imus;
    out.WriteImuObs();
}


int main(int argc,char** argv){
    string logini_path = SetLogConfPath("");
    int log_level = SetLogLevel(32);
    InitLog(argc,argv,logini_path, log_level);
    tPPPLibConf C;

    Vector3d rr(-2364333.7761,4870287.3700,-3360809.2451);
    Vector3d blh=Xyz2Blh(rr);
    fprintf(stdout,"%9.6f %9.6f %9.6f\n",blh[0]*R2D,blh[1]*R2D,blh[2]);

#if 0
    string imu_file="/home/cc/dataset/data_carvig_lc/m390311.imu";
    C.insC.imu_type=IMU_M39;
    C.insC.sample_rate=200.0;
#endif

#if 1
    string imu_file="/home/cc/dataset/data_cpt0/cpt00870.imu";
    C.insC.imu_type=IMU_NOVTEL_CPT;
    C.insC.sample_rate=100.0;
#endif

#if 0
    string imu_file="/home/cc/dataset/data_urban_hk/mti01180.imu";
    C.insC.imu_type=IMU_MTI_CSV;
    C.insC.sample_rate=100.0;
#endif

    FormatImu(C,imu_file);

    C.mode=MODE_IGLC;
    string ref_file="/home/cc/dataset/data_cpt0/cpt00870.ref";
    FormatRefSol(C,ref_file,1);

    ref_file="/home/cc/dataset/data_urban_hk/cpt01180.ref";
    FormatRefSol(C,ref_file,0);
}
