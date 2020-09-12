//
// Created by cc on 8/3/20.
//

#include "OutSol.h"

namespace PPPLib{

    cOutSol::cOutSol() {}

    cOutSol::cOutSol(PPPLib::tPPPLibConf C) {
        C_=C;
    }

    cOutSol::cOutSol(PPPLib::tPPPLibConf C,vector<tSolInfoUnit>& ref_sols) {C_=C;ref_sols_=ref_sols;}

    cOutSol::~cOutSol() {}

    int cOutSol::MatchRefSol(PPPLib::cTime sol_time) {
        double sow1,sow2;
        int i=0;
        bool stat=false;

        for(i=ref_index_-100<0?0:ref_index_-10;i<ref_sols_.size();i++){
            sow1=ref_sols_.at(i).t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            sow2=sol_time.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<DTTOL){
                ref_index_=i;stat=true;break;
            }
            else if((sow1-sow2)>2.0*DTTOL){
                stat=false;break;
            }
        }
        return stat;
    }

    tSolInfoUnit cOutSol::CompareSol(PPPLib::tSolInfoUnit &sol, PPPLib::tSolInfoUnit &ref_sol) {
        tSolInfoUnit dsol;
        Vector3d ref_blh=Xyz2Blh(ref_sol.pos);
        Vector3d dxyz=sol.pos-ref_sol.pos;
        if(C_.solC.sol_coord==COORD_ENU){
            dsol.pos=Xyz2Enu(ref_blh,dxyz);
        }
        dsol.vel=sol.vel-ref_sol.vel;

        return dsol;
    }

    static double SqrtCovVar(double cov_var){
        return cov_var<0.0?-sqrt(-cov_var):sqrt(cov_var);
    }

    int cOutSol::OutSolStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos, char *buff) {
        if(sol->stat<=SOL_NONE) return 0;

        if(C_.mode==MODE_PPP||C_.mode_opt==MODE_OPT_PPP){
            return 0;
        }

        char *p=buff;

        int week;
        double tow;
        tow=sol->t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   sol->stat,sol->pos[0],sol->pos[1],sol->pos[2],
                   0.0,0.0,0.0);

        if(C_.gnssC.ion_opt>ION_IF_DUAL){

        }

        if(C_.gnssC.trp_opt>TRP_SAAS){

        }


        return (int)(p-buff);
    }

    void cOutSol::WriteSatStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos) {
        char buff[8191+1];

        int n=OutSolStat(sol,sat_infos,buff);
        buff[n]='\0';

        fputs(buff,fout_stat_);
        if(sol->stat<=SOL_NONE) return;

        tSatInfoUnit *sat_info= nullptr;
        int nf=2;
        int i,j,week;
        double tow,amb;
        tow=sol->t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
        for(i=0;i<MAX_SAT_NUM;i++){
            sat_info=&sat_infos[i];
            if(!sat_info->vsat[0]) continue;
            for(j=0;j<nf;j++){
                for (j=0;j<nf;j++) {
                    fprintf(fout_stat_,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f\n",
                               week,tow,sat_info->sat.sat_.id.c_str(),j+1,sat_info->el_az[1]*R2D,sat_info->el_az[0]*R2D,
                               sat_info->post_res_P[j],sat_info->post_res_L[j],sat_info->vsat[j],sat_info->raw_S[j]*0.25,
                               sat_info->fix[j],sat_info->slip[j]&3,sat_info->lock[j],sat_info->outc[j],
                               0,sat_info->rejc[j],sat_info->float_amb[j],sat_info->fix_amb[j],sat_info->raw_mw[0],sat_info->sm_mw[0]);
                }
            }
        }
    }

    int cOutSol::OutEcef(unsigned char *buff, const char *s, tSolInfoUnit &sol) {
        const char *sep=" ";
        char *p=(char *)buff;
        Vector3d blh=Xyz2Blh(sol.pos);

        p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f%s%6.2f",
                   s,sep,sol.pos[0],sep,sol.pos[1],sep,sol.pos[2],sep,sol.stat,sep,
                   sol.valid_sat_num,sep,SQRT(sol.q_pos[0]),sep,SQRT(sol.q_pos[1]),sep,SQRT(sol.q_pos[2]),
                   sep,SqrtCovVar(sol.q_pos[3]),sep,SqrtCovVar(sol.q_pos[4]),sep,SQRT(sol.q_pos[5]),
                   sep,sol.age,sep,sol.ratio,sep,sol.sigma);
        if(C_.mode>=MODE_INS){
            p+=sprintf(p,"%s%4d",sep,sol.ins_stat);
        }

        if(C_.solC.out_vel){
            p+=sprintf(p,"%s%10.5f%s%10.5f%s%10.5f%s%9.5f%s%8.5f%s%8.5f%s%8.5f%s%8.5f%s%8.5f",
                       sep,sol.vel[0],sep,sol.vel[1],sep,sol.vel[2],sep,
                       SQRT(sol.q_vel[0]),sep,SQRT(sol.q_vel[1]),sep,SQRT(sol.q_vel[2]),
                       sep,SqrtCovVar(sol.q_vel[3]),sep,SqrtCovVar(sol.q_vel[4]),sep,
                       SqrtCovVar(sol.q_vel[5]));
        }
        if(C_.solC.out_att&&C_.mode>=MODE_INS){
            p+=sprintf(p,"%s%10.6lf%s%10.6lf%s%10.6lf%s%10.6lf%s%10.6lf%s%10.6lf",
                       sep,sol.att[0]*R2D,sep,sol.att[1]*R2D,sep,sol.att[2]*R2D<0.0?sol.att[2]*R2D+360:sol.att[2]*R2D,
                       sep,SQRT(sol.q_att[0])*R2D,sep,SQRT(sol.q_att[1])*R2D,sep,SQRT(sol.q_att[2])*R2D);
        }
        if(C_.solC.out_ba&&C_.mode>=MODE_IGLC){
            p+=sprintf(p,"%s%10.6lf%s%10.6lf%s%10.6lf",sep,sol.accl_bias[0],sep,sol.accl_bias[1],sep,sol.accl_bias[2]);

        }
        if(C_.solC.out_bg&&C_.mode>=MODE_IGLC){
            p+=sprintf(p,"%s%10.6lf%s%10.6lf%s%10.6lf",sep,sol.gyro_bias[0],sep,sol.gyro_bias[1],sep,sol.gyro_bias[2]);
        }

        p+=sprintf(p,"\n");
        return p-(char*)buff;
    }

    bool cOutSol::InitOutSol(tPPPLibConf C, string file) {
        C_=C;
        fout_=fopen(file.c_str(),"w");
        if(C_.solC.out_stat&&!C_.fileC.sol_stat.empty()){
            fout_stat_=fout_stat_=fopen(C_.fileC.sol_stat.c_str(),"w");
        }
    }

    void cOutSol::WriteHead() {
        if(!C_.solC.out_head) return;
        const char *sep=" ";
        fprintf(fout_,"%s  %-*s%s%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s%s%6s",
                   COMMENTH,20,"GPST",sep,"x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"Q",sep,"ns",sep,
                   "sdx(m)",sep,"sdy(m)",sep,"sdz(m)",sep,"sdxy(m)",sep,
                   "sdyz(m)",sep,"sdzx(m)",sep,"age(s)",sep,"ratio",sep,"sigma");
        if(C_.mode>=MODE_INS){
            fprintf(fout_,"%s%4s",sep,"QINS");
        }
        if(C_.solC.out_vel){
            fprintf(fout_,"%s%10s%s%10s%s%10s%s%9s%s%8s%s%8s%s%8s%s%8s%s%8s",
                       sep,"vn(m/s)",sep,"ve(m/s)",sep,"vu(m/s)",sep,"sdvn",sep,
                       "sdve",sep,"sdvu",sep,"sdvne",sep,"sdveu",sep,"sdvun");
        }
        if(C_.solC.out_att&&C_.mode>=MODE_INS){
            fprintf(fout_,"%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s",
                       sep,"roll(deg)",sep,"pitch(deg)",sep,"yaw(deg)",
                       sep,"sdroll",sep,"sdpitch",sep,"sdyaw");
        }
        if(C_.solC.out_ba&&C_.mode>=MODE_IGLC){
            fprintf(fout_,"%s%10s%s%10s%s%10s",sep,"bax(m/s^2)",sep,"bay(m/s^2)",sep,"baz(m/s^2)");
        }
        if(C_.solC.out_bg&&C_.mode>=MODE_IGLC){
            fprintf(fout_,"%s%10s%s%10s%s%10s",sep,"bgx(rad/s)",sep,"bgy(rad/s)",sep,"bgz(rad/s)");
        }
        fprintf(fout_,"\n");
    }

    void cOutSol::WriteSol(tSolInfoUnit sol,int epoch) {
        tSolInfoUnit dsol;
        dsol=sol;

        unsigned char p[8191+1];
        char s[256];

        int n=OutEcef(p, sol.t_tag.GetTimeStr(3).c_str(),dsol);
        fwrite(p,n,1,fout_);
    }

    void cOutSol::WriteImuHead() {
        if(!C_.solC.out_head) return;
        const char *sep=" ";
        fprintf(fout_,"%s %4s%s%14s%s%14s%s%14s%s%14s%s%14s%s%14s%s%14s",
                COMMENTH,"GPSW",sep,"SOW",sep,"x-acce(m/s^2)",sep,"y-acce(m/s^2)",sep,"z-acce(m/s^2)",sep,"x-gyro(rad/s)",sep,"y-gyro(rad/s)",sep,"z-gyro(rad/s)");
        fprintf(fout_,"\n");

    }

    void cOutSol::WriteImuObs() {
        int n=imus->data_.size();
        cTime t_tag;
        tImuDataUnit data;
        const char *sep=" ";
        int week;
        double sow;

        for(int i=0;i<n;i++){
            data=imus->data_[i];
            sow=data.t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
            fprintf(fout_,"%6d%s%14.4f%s%14.5f%s%14.5f%s%14.5f%s%14.5f%s%14.5f%s%14.5f\n",
                    week,sep,sow,sep,data.acce[0],sep,data.acce[1],sep,data.acce[2],sep,data.gyro[0],sep,data.gyro[1],sep,data.gyro[2]);
        }
        fclose(fout_);
    }

}