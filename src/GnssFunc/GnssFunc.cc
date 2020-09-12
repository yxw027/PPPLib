//
// Created by cc on 7/16/20.
//

#include "GnssFunc.h"

namespace PPPLib{
    cSat::cSat() {}

    cSat::cSat(int sat_no) {
        sat_.no=sat_no;
        SatNo2Prn();
    }

    cSat::cSat(string sat_id) {sat_.id=sat_id;}

    cSat::cSat(int sat_sys,int sat_prn){sat_.sys=sat_sys;sat_.prn=sat_prn;}

    cSat::~cSat() {}

    void cSat::SatPrn2No() {
        if(sat_.prn<=0){
            LOG(ERROR)<<"satellite prn error";
            return;
        }

        switch(sat_.sys){
            case SYS_GPS:
                if(sat_.prn<GPS_MIN_PRN||GPS_MAX_PRN<sat_.prn){sat_.no=0;return;}
                sat_.no=sat_.prn-GPS_MIN_PRN+1;sat_.sys_idx=SYS_INDEX_GPS;break;
            case SYS_BDS:
                if(sat_.prn<BDS_MIN_PRN||BDS_MAX_PRN<sat_.prn){sat_.no=0;return;}
                sat_.no=NUM_GPS_SAT+sat_.prn-BDS_MIN_PRN+1;sat_.sys_idx=SYS_INDEX_BDS;break;
            case SYS_GAL:
                if(sat_.prn<GAL_MIN_PRN||GAL_MAX_PRN<sat_.prn){sat_.no=0;return;}
                sat_.no=NUM_GPS_SAT+NUM_BDS_SAT+sat_.prn-GAL_MIN_PRN+1;sat_.sys_idx=SYS_INDEX_GAL;break;
            case SYS_GLO:
                if(sat_.prn<GLO_MIN_PRN||GLO_MAX_PRN<sat_.prn){sat_.no=0;return;}
                sat_.no=NUM_GPS_SAT+NUM_BDS_SAT+NUM_GAL_SAT+sat_.prn-GLO_MIN_PRN+1;sat_.sys_idx=SYS_INDEX_GLO;break;
            case SYS_QZS:
                if(sat_.prn<QZS_MIN_PRN||QZS_MAX_PRN<sat_.prn){sat_.no=0;return;}
                sat_.no=NUM_GPS_SAT+NUM_BDS_SAT+NUM_GAL_SAT+NUM_GLO_SAT+sat_.prn-QZS_MIN_PRN+1;sat_.sys_idx=SYS_INDEX_QZS;break;
            case SYS_IRN:
                sat_.no=0;break;
            default:
                sat_.no=0;break;
        }
    }

    void cSat::SatNo2Prn() {
        int prn=sat_.no;
        sat_.sys=SYS_NONE;
        if(sat_.no<=0||MAX_SAT_NUM<sat_.no){
//            LOG(ERROR)<<"satellite no. error";
            return;
        }
        else if(prn<=NUM_GPS_SAT){
            sat_.sys=SYS_GPS;prn+=GPS_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_GPS;
        }
        else if((prn-=NUM_GPS_SAT)<=NUM_BDS_SAT){
            sat_.sys=SYS_BDS;prn+=BDS_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_BDS;
        }
        else if((prn-=NUM_BDS_SAT)<=NUM_GAL_SAT){
            sat_.sys=SYS_GAL;prn+=GAL_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_GAL;
        }
        else if((prn-=NUM_GAL_SAT)<=NUM_GLO_SAT){
            sat_.sys=SYS_GLO;prn+=GLO_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_GLO;
        }
        else if((prn-=NUM_GLO_SAT)<=NUM_QZS_SAT){
            sat_.sys=SYS_QZS;prn+=QZS_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_QZS;
        }
        else if((prn-=NUM_QZS_SAT)<=NUM_IRN_SAT){
            sat_.sys=SYS_IRN;prn=0;
        }
        else prn=0;
        sat_.prn=prn;
    }

    void cSat::SatNo2Id() {
        SatNo2Prn();
        string buff;
        switch (sat_.sys){
            case SYS_GPS: sat_.id="G"+Int2Str(2,"0",sat_.prn-GPS_MIN_PRN+1,buff);sat_.sys_idx=SYS_INDEX_GPS; break;
            case SYS_BDS: sat_.id="C"+Int2Str(2,"0",sat_.prn-BDS_MIN_PRN+1,buff);sat_.sys_idx=SYS_INDEX_BDS; break;
            case SYS_GAL: sat_.id="E"+Int2Str(2,"0",sat_.prn-GAL_MIN_PRN+1,buff);sat_.sys_idx=SYS_INDEX_GAL; break;
            case SYS_GLO: sat_.id="R"+Int2Str(2,"0",sat_.prn-GLO_MIN_PRN+1,buff);sat_.sys_idx=SYS_INDEX_GLO; break;
            case SYS_QZS: sat_.id="J"+Int2Str(2,"0",sat_.prn-QZS_MIN_PRN+1,buff);sat_.sys_idx=SYS_INDEX_QZS; break;
            default:
                sat_.id="";break;
        }
    }

    void cSat::SatId2No() {
        sat_.sys=SYS_NONE;

        if(Str2Int((sat_.id.substr(1,2)),sat_.prn)==0) return;
        char code= sat_.id[0];

        switch(code){
            case 'G': sat_.sys=SYS_GPS;sat_.prn+=GPS_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_GPS;break;
            case 'C': sat_.sys=SYS_BDS;sat_.prn+=BDS_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_BDS;break;
            case 'E': sat_.sys=SYS_GAL;sat_.prn+=GAL_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_GAL;break;
            case 'R': sat_.sys=SYS_GLO;sat_.prn+=GLO_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_GLO;break;
            case 'J': sat_.sys=SYS_QZS;sat_.prn+=QZS_MIN_PRN-1;sat_.sys_idx=SYS_INDEX_QZS;break;
            case 'I': sat_.sys=SYS_IRN;break;
            case 'S': sat_.sys=SYS_SBS;break;
            default:
                sat_.sys=SYS_NONE;sat_.prn=0;break;
        }
        SatPrn2No();
    }

    cGnssObs::cGnssObs() {rcv_idx_=REC_ROVER;}

    cGnssObs::~cGnssObs() {}

    void cGnssObs::SetTimeSpan(cTime *ts, cTime *te) {
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cTime* cGnssObs::GetStartTime() { return &ts_;}

    cTime * cGnssObs::GetEndTime() {return &te_;}

    void cGnssObs::SetRcvIdx(RECEIVER_INDEX rcv) {rcv_idx_=rcv;}

    vector<tEpochSatUnit>& cGnssObs::GetGnssObs() {return obs_;}

    tStaInfoUnit* cGnssObs::GetStation() {return &station_;}

    void cGnssObs::CarrierSmoothCode(int smooth_length) {
        tSatObsUnit* sat_obs= nullptr;
        int n[MAX_SAT_NUM][MAX_GNSS_FRQ_NUM]={{0}},sat_no;
        int sys_idx;
        double last_csc_P[MAX_SAT_NUM][MAX_GNSS_FRQ_NUM]={{0}},csc_P[MAX_SAT_NUM][MAX_GNSS_FRQ_NUM]={{0}};
        double tdcp=0.0,lam=0.0;
        for(int i=0;i<epoch_num;i++){
            for(int j=0;j<obs_.at(i).sat_num;j++){
                sat_obs=&obs_.at(i).epoch_data.at(j);
                sat_no=sat_obs->sat.sat_.no;
                sys_idx=sat_obs->sat.sat_.sys_idx;
                for(int k=0;k<MAX_GNSS_FRQ_NUM;k++){
                    if(sat_obs->P[k]==0.0||sat_obs->L[k]==0.0) continue;
                    //周跳
                    if(sat_obs->LLI[k]) n[sat_no-1][k]=0;
                    if(n[sat_no-1][k]==0) csc_P[sat_no-1][k]=sat_obs->P[k];
                    else{
                        lam=CLIGHT/(kGnssFreqs[sys_idx][k]);
                        tdcp=lam*(sat_obs->L[k]-last_csc_P[sat_no-1][k]);
                        csc_P[sat_no-1][k]=sat_obs->P[k]/smooth_length+(csc_P[sat_no-1][k]+tdcp)*(smooth_length-1)/smooth_length;
                    }
                    if(++n[sat_no-1][k]<smooth_length) n[sat_no-1][k]=0;
                    else sat_obs->CSC_P[k]=csc_P[sat_no-1][k];
                    last_csc_P[sat_no-1][k]=sat_obs->L[k];
                }
            }
        }
    }

    cGnssObsOperator::cGnssObsOperator() {}

    cGnssObsOperator::~cGnssObsOperator() {}

    void cGnssObsOperator::ReAlignObs(tPPPLibConf C, tSatInfoUnit &sat_info, tSatObsUnit sat_obs, int f, int raw_idx,int frq_idx,int *glo_frq) {
        sat_info.P_code[f]=sat_obs.code[frq_idx];
        sat_info.raw_P[f]=C.gnssC.csc?sat_obs.CSC_P[frq_idx]:sat_obs.P[frq_idx];
        sat_info.raw_L[f]=sat_obs.L[frq_idx];
        sat_info.raw_D[f]=sat_obs.D[frq_idx];
        sat_info.raw_S[f]=sat_obs.SNR[frq_idx];
        sat_info.LLI[f]=sat_obs.LLI[frq_idx];

        if(sat_info.sat.sat_.sys==SYS_GLO&&glo_frq){
            sat_info.frq[f]=kGnssFreqs[sat_obs.sat.sat_.sys_idx][raw_idx]+glo_frq[sat_info.sat.sat_.prn-1]*kGnssFreqs[sat_obs.sat.sat_.sys_idx][raw_idx+2];
            sat_info.lam[f]=CLIGHT/sat_info.frq[f];
        }
        else{
            sat_info.frq[f]=kGnssFreqs[sat_obs.sat.sat_.sys_idx][raw_idx];
            sat_info.lam[f]=CLIGHT/kGnssFreqs[sat_obs.sat.sat_.sys_idx][raw_idx];
        }

    }

    double cGnssObsOperator::GnssObsIfComb(double obs1, double obs2, double f1, double f2) {
        double alpha=SQR(f1)/(SQR(f1)-SQR(f2)),beta=-SQR(f2)/(SQR(f1)-SQR(f2));
        double if_obs=0.0;
        if(obs1==0.0||obs2==0.0) return 0;

        if_obs=obs1*alpha+beta*obs2;
        return if_obs;
    }

    double cGnssObsOperator::GnssObsMwComb(double obs_P1, double obs_P2,double obs_L1,double obs_L2, double lam1, double lam2) {
        double obs_mw=0.0;
        if(obs_P1==0.0||obs_P2==0.0||obs_L1==0.0||obs_L2==0.0) return 0.0;
        obs_mw=(obs_L1-obs_L2)-(lam2-lam1)/(lam1+lam2)*(obs_P1/lam1+obs_P2/lam2);
        return obs_mw;
    }

    double cGnssObsOperator::GnssObsGfComb(double obs1, double obs2,double lam1,double lam2) {
        if(obs1==0.0||obs2==0.0||lam1==0.0||lam2==0.0) return 0.0;
        return obs1*lam1-obs2*lam2;
    }

    double cGnssObsOperator::GnssObsCmcComb(double obs_P, double obs_L1, double obs_L2, double f1, double f2) {
        double cmc=0.0;
        cmc=obs_P-obs_L1-2*SQR(f1)/(SQR(f1)+SQR(f2))*obs_L2;
        return cmc;
    }

    double cGnssObsOperator::GnssObsTdComb(double cur_obs, double pre_obs) {
        double obs_td=0.0;
        obs_td=cur_obs-pre_obs;
        return obs_td;
    }

    void cGnssObsOperator::MakeGnssObsComb(tPPPLibConf C, GNSS_OBS_COMB type, tSatInfoUnit* sat_info,const tSatInfoUnit previous_sat_info) {
        if(type==COMB_IF){
            double alpha=0.0,beta=0.0;
            int sys=sat_info->sat.sat_.sys;
            if(C.gnssC.frq_opt==FRQ_SINGLE){
                //单频只有一个P、L，且保证P、L不为0
                if(sat_info->cor_P[0]==0.0||sat_info->cor_L[0]==0.0) sat_info->cor_if_P[0]=0.0;
                else sat_info->cor_if_P[0]=0.5*sat_info->cor_P[0]+0.5*sat_info->cor_L[0];
            }
            else if(C.gnssC.frq_opt==FRQ_DUAL){ //双频只有两个P、L,且保证所有P、L不为0
                sat_info->cor_if_P[0]=GnssObsIfComb(sat_info->cor_P[0],sat_info->cor_P[1],sat_info->frq[0],sat_info->frq[1]);
                sat_info->cor_if_L[0]=GnssObsIfComb(sat_info->cor_L[0],sat_info->cor_L[1],sat_info->frq[0],sat_info->frq[1]);
            }
            else if(C.gnssC.frq_opt==FRQ_TRIPLE){
                if(C.gnssC.ion_opt==ION_IF){
                    double f1=sat_info->frq[0],f2=sat_info->frq[1],f3=sat_info->frq[2];
                    if(f1*f2*f3!=0.0){
                        //使用一个3频组合,L1_L2_L5/L1_L2/L1_L5
                        double gam1=(SQR(f1)/SQR(f1)),gam2=SQR(f1)/SQR(f2),gam3=SQR(f1)/SQR(f3);
                        double e=2*(SQR(gam2)+SQR(gam3)-gam2*gam3-gam2-gam3+1.0);
                        double e1=(SQR(gam2)+SQR(gam3)-gam2-gam3)/e;
                        double e2=(SQR(gam3)-gam2*gam3-gam2+1.0)/e;
                        double e3=(SQR(gam2)-gam2*gam3-gam3+1.0)/e;
                        double P1=sat_info->cor_P[0],P2=sat_info->cor_P[1],P3=sat_info->cor_P[2];
                        double L1=sat_info->cor_L[0],L2=sat_info->cor_L[1],L3=sat_info->cor_L[2];
                        if(P1==0.0||P2==0.0||P3==0.0) sat_info->cor_if_P[0]=0.0;
                        else sat_info->cor_if_P[0]=e1*P1+e2*P2+e3*P3;
                        if(L1==0.0||L2==0.0||L3==0.0) sat_info->cor_if_L[0]=0.0;
                        else sat_info->cor_if_L[0]=e1*L1+e2*L2+e3*L3;

                        // 三频观测缺失，使用双频替补
                        if(sat_info->cor_if_P[0]==0.0){
                            if(P1&&P2){
                                sat_info->cor_if_P[1]=GnssObsIfComb(P1,P2,sat_info->frq[0],sat_info->frq[1]);
                            }
                            else if(P1&&P3){
                                sat_info->cor_if_P[2]=GnssObsIfComb(P1,P2,sat_info->frq[0],sat_info->frq[2]);
                            }
                            else sat_info->stat=SAT_NO_USE;
                        }

                        if(sat_info->cor_if_L[0]==0.0){
                            if(P1&&P2){
                                sat_info->cor_if_L[1]=GnssObsIfComb(L1,L2,sat_info->frq[0],sat_info->frq[1]);
                            }
                            else if(P1&&P3){
                                sat_info->cor_if_L[2]=GnssObsIfComb(L1,L2,sat_info->frq[0],sat_info->frq[2]);
                            }
                            else sat_info->stat=SAT_NO_USE;
                        }


                    }
                }
                else if(C.gnssC.ion_opt==ION_IF_DUAL){
                    // IF_L1L2
                    sat_info->cor_if_P[0]=GnssObsIfComb(sat_info->cor_P[0],sat_info->cor_P[1],sat_info->frq[0],sat_info->frq[1]);
                    sat_info->cor_if_L[0]=GnssObsIfComb(sat_info->cor_L[0],sat_info->cor_L[1],sat_info->frq[0],sat_info->frq[1]);
                    // IF_L1L5
                    sat_info->cor_if_P[1]=GnssObsIfComb(sat_info->cor_P[0],sat_info->cor_P[2],sat_info->frq[0],sat_info->frq[2]);
                    sat_info->cor_if_L[1]=GnssObsIfComb(sat_info->cor_L[0],sat_info->cor_L[2],sat_info->frq[0],sat_info->frq[2]);
                }
            }
        }
        else if(type==COMB_GF){
            if(C.gnssC.frq_opt==FRQ_SINGLE) return;
            else if(C.gnssC.frq_opt==FRQ_DUAL){
                sat_info->gf[0]=GnssObsGfComb(sat_info->raw_L[0],sat_info->raw_L[1],sat_info->lam[0],sat_info->lam[1]);
            }
            else if(C.gnssC.frq_opt==FRQ_TRIPLE){

            }
        }
        else if(type==COMB_MW){
            if(C.gnssC.frq_opt==FRQ_DUAL) return;
            else if(C.gnssC.frq_opt==FRQ_DUAL){
                sat_info->raw_mw[0]=GnssObsMwComb(sat_info->cor_P[0],sat_info->cor_P[1],sat_info->raw_L[0],sat_info->raw_L[1],sat_info->lam[0],sat_info->lam[1]);
            }
            else if(C.gnssC.frq_opt==FRQ_TRIPLE){

            }
        }
        else if(type==COMB_MP){
            if(C.gnssC.frq_opt==FRQ_DUAL){
                sat_info->cmc_P=GnssObsCmcComb(sat_info->cor_P[0],sat_info->cor_L[0],sat_info->cor_L[1],sat_info->frq[0],sat_info->frq[1]);
            }
        }
        else if(type==COMB_TDCP){
            //time-difference carrier phase
            if(previous_sat_info.cor_L[0]!=0.0){
                sat_info->tdcp=GnssObsTdComb(sat_info->cor_L[0],previous_sat_info.cor_L[0]);
            }
        }
    }

    double cGnssObsOperator::GnssSdObs(const tSatInfoUnit &sat1, const tSatInfoUnit &sat2, int f, GNSS_OBS type) {
        double obs1=0.0,obs2=0.0;
        if(type==GNSS_OBS_CODE){
            obs1=sat1.raw_P[f];
            obs2=sat2.raw_P[f];
        }
        else if(type==GNSS_OBS_PHASE){
            obs1=sat1.raw_L[f];
            obs2=sat2.raw_L[f];
        }
        return obs1==0.0||obs2==0.0?0.0:obs1-obs2;
    }

    double cGnssObsOperator::GnssObsLinearComb(tPPPLibConf  C,int i, int j, int k, tSatInfoUnit &sat_info, GNSS_OBS type,double* var) {
        double obs1,obs2,obs3;
        double f1,f2,f3;
        if(type==GNSS_OBS_CODE){
            obs1=sat_info.raw_P[0];
            obs2=sat_info.raw_P[1];
            obs3=sat_info.raw_P[2];
        }
        else if(type==GNSS_OBS_PHASE){
            obs1=sat_info.raw_L[0]*sat_info.lam[0];
            obs2=sat_info.raw_L[1]*sat_info.lam[1];
            obs3=sat_info.raw_L[2]*sat_info.lam[2];
        }

        f1=sat_info.frq[0];
        f2=sat_info.frq[1];
        f3=sat_info.frq[2];

        if((i&&!obs1)||(j&&!obs2)||(k&&!obs3)){
            return 0.0;
        }
        double a=C.gnssC.meas_err_factor[1];
        double b=C.gnssC.meas_err_factor[2];
        double sig=0.0;
        if(var){
            sig=C.gnssC.code_phase_ratio*(sqrt(SQR(a)+SQR(b))/sin(sat_info.el_az[0]));
            *var=(SQR(i*f1)+SQR(j*f2)+SQR(k*f3))/SQR(i*f1+j*f2+k*f3)*SQR(sig);
        }

        return (i*f1*obs1+j*f2*obs2+k*f3*obs2)/(i*f1+j*f2+k*f3);
    }

    double cGnssObsOperator::LinearCombLam(int i, int j, int k, tSatInfoUnit &sat_info) {
        double f1=sat_info.frq[0],f2=sat_info.frq[1],f3=sat_info.frq[2];
        return CLIGHT/(i*f1+j*f2+k*f3);
    }

    void cGnssObsOperator::MwCycleSlip(tPPPLibConf C,double sample_dt,double dt,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat,tTime last_time) {
        if(C.gnssC.frq_opt==FRQ_SINGLE) return;
        int del_ep=1;
        double fact=1.0;
        if(sample_dt>0.0) del_ep=Round(dt/sample_dt);
        if(sample_dt<=1.5){
            if(fabs(dt)<=10.0) del_ep=1;
            else if(fabs(dt)<=15.0) del_ep=2;
            else if(fabs(dt)<=22.0) del_ep=3;
        }
        if(sample_dt>=29.5){
            if(del_ep<=2) fact=1.0;
            else if(del_ep<=4) fact=1.25;
            else if(del_ep<=6) fact=1.5;
            else fact=2.0;
        }

        double el=sat_info->el_az[0]*R2D;
        double P1,P2,L1,L2;
        double lam1,lam2;
        int num_loop=C.gnssC.frq_opt==FRQ_TRIPLE?2:1;
        for(int f=0;f<num_loop;f++){
            double w1=0.0,w0=0.0;
            if(sat_info->frq[f]==0.0) continue;
            double a=sat_info->t_tag.TimeDiff(last_time);
            if(a>sample_dt){
                sat_info->sm_mw[f]=sat_info->mw_idx[f]=0;
            }
            if(f==0){
                if(base_sat){
                    P1=GnssSdObs(*sat_info,*base_sat,0,GNSS_OBS_CODE);
                    P2=GnssSdObs(*sat_info,*base_sat,1,GNSS_OBS_CODE);
                    L1=GnssSdObs(*sat_info,*base_sat,0,GNSS_OBS_PHASE);
                    L2=GnssSdObs(*sat_info,*base_sat,1,GNSS_OBS_PHASE);
                    lam1=sat_info->lam[0];lam2=sat_info->lam[1];
                }
                else{
                    P1=sat_info->raw_P[0];P2=sat_info->raw_P[1];
                    L1=sat_info->raw_L[0];L2=sat_info->raw_L[1];
                    lam1=sat_info->lam[0];lam2=sat_info->lam[1];
                }
                w1=GnssObsMwComb(P1,P2,L1,L2,lam1,lam2);
                w0=sat_info->sm_mw[0];
            }
            else if(f==1){
//                if(sat_info->cor_P[2]==0.0||sat_info->raw_L[2]==0.0) continue;
//                w1=GnssObsMwComb(sat_info->raw_P[0],sat_info->raw_P[2],sat_info->raw_L[0],sat_info->raw_L[2],sat_info->lam[0],sat_info->lam[2]);
//                w0=sat_info->sm_mw[1];
            }

            if(w1==0.0||w0==0.0){
                sat_info->sm_mw[f]=sat_info->mw_idx[f]=0;
                continue;
            }

            double thres=0.0;
            if(el>20.0) thres=C.gnssC.cs_thres[0];
            else thres=-C.gnssC.cs_thres[0]*0.1*el+3.0*C.gnssC.cs_thres[0];
            if(el<C.gnssC.ele_min) continue;
            double dmw=w1-w0;
            if(fabs(dmw)>MIN(thres*fact,6.0)){
                LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<"MW DETECT CYCLE SLIP dmw=mw1-mw0 "<<dmw<<"="<<w1<<"-"<<w0<<" THRESHOLD="<<MIN(thres*fact,6.0);
                sat_info->slip[f]|=1;
            }
        }
    }

    void cGnssObsOperator::GfCycleSlip(tPPPLibConf C,double sample_dt,double dt, tSatInfoUnit *sat_info, tSatInfoUnit* base_sat) {
        if(C.gnssC.frq_opt==FRQ_SINGLE) return;
        int del_ep;
        del_ep=Round(fabs(dt/sample_dt));

        double el=sat_info->el_az[0]*R2D;
        int num_loop=C.gnssC.frq_opt==FRQ_TRIPLE?2:1;
        double L1,L2;
        double lam1,lam2;
        for(int f=0;f<num_loop;f++){
            double g1=0.0,g0=0.0;
            if(f==0){
                if(base_sat){
                    L1=GnssSdObs(*sat_info,*base_sat,0,GNSS_OBS_PHASE);
                    L2=GnssSdObs(*sat_info,*base_sat,1,GNSS_OBS_PHASE);
                    lam1=sat_info->lam[0];
                    lam2=sat_info->lam[1];
                }
                else{
                    L1=sat_info->raw_L[0];
                    L2=sat_info->raw_L[1];
                    lam1=sat_info->lam[0];
                    lam2=sat_info->lam[1];
                }
                g1=GnssObsGfComb(L1,L2,lam1,lam2);
                g0=sat_info->gf[0];
            }
            else if(f==1){
//                if(sat_info->cor_L[2]==0.0) continue;
//                g1=GnssObsGfComb(sat_info->raw_L[0],sat_info->raw_L[2],sat_info->lam[0],sat_info->lam[1]);
//                g0=sat_info->gf[1];
            }

            if(g1==0.0||g0==0.0) continue;
            double thres;
            if(el<C.gnssC.ele_min) continue;
            if(el>15.0) thres=C.gnssC.cs_thres[1];
            else thres=-C.gnssC.cs_thres[1]/15.0*el+2.0*C.gnssC.cs_thres[1];
//            thres=C.gnssC.cs_thres[1];

            double dgf=g1-g0;
            if(fabs(dgf)>thres){
                LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<"GF DETECT CYCLE SLIP dgf=gf1-gf0 "<<dgf<<"="<<g1<<"-"<<g0<<" THRESHOLD="<<thres;
                sat_info->slip[0]|=1;
                sat_info->slip[1]|=1;
            }
        }
    }

    static unsigned int GetBitu(const unsigned char *buff,int pos,int len){
        unsigned int bits=0;
        int i;
        for (i=pos;i<pos+len;i++) bits=(bits<<1)+((buff[i/8]>>(7-i%8))&1u);
        return bits;
    }

    static int GetBits(const unsigned char *buff, int pos,int len){
        unsigned int bits=GetBitu(buff,pos,len);
        if (len<=0||32<=len||!(bits&(1u<<(len-1)))) return (int)bits;
        return (int)(bits|(~0u<<len)); /* extend sign */
    }

    static void SetBitu(unsigned char *buff, int pos, int len, unsigned int data)
    {
        unsigned int mask=1u<<(len-1);
        int i;
        if (len<=0||32<len) return;
        for (i=pos;i<pos+len;i++,mask>>=1) {
            if (data&mask) buff[i/8]|=1u<<(7-i%8); else buff[i/8]&=~(1u<<(7-i%8));
        }
    }
    static void SetBits(unsigned char *buff, int pos, int len, int data)
    {
        if (data<0) data|=1<<(len-1); else data&=~(1<<(len-1)); /* set sign bit */
        SetBitu(buff,pos,len,(unsigned int)data);
    }

    bool cGnssObsOperator::LliCycleSlip(tPPPLibConf C, tSatInfoUnit& sat_info,int nf,double tt,RECEIVER_INDEX rcv) {
        int f;
        unsigned int slip=0,LLI;
        int slip_flag=false;

        for(f=0;f<nf;f++){
            if((sat_info.raw_L[f]==0.0||sat_info.LLI[f]==0)){
                continue;
            }
            if(rcv==REC_ROVER) LLI=GetBitu(&sat_info.slip[f],0,2);
            else LLI=GetBitu(&sat_info.LLI[f],2,2);

            if(tt>=0.0){
                if(sat_info.LLI[f]&1){
                    LOG(WARNING)<<sat_info.t_tag.GetTimeStr(1)<<" "<<sat_info.sat.sat_.id<<" L"<<f+1<<" SLIP DETECTED FORWARD";
                }
                slip=sat_info.LLI[f];
            }
            else{
                if(LLI&1){
                    LOG(WARNING)<<sat_info.t_tag.GetTimeStr(1)<<" "<<sat_info.sat.sat_.id<<" L"<<f+1<<" SLIP DETECTED BACKWARD";
                }
                slip=LLI;
            }

            if(rcv==REC_ROVER) SetBitu(&sat_info.slip[f],0,2,sat_info.LLI[f]);
            else SetBitu(&sat_info.slip[f],2,2,sat_info.LLI[f]);

            sat_info.slip[f]|=(unsigned char)slip;
            if(slip) slip_flag=true;
        }
        return slip_flag;
    }

    void cGnssObsOperator::SmoothMw(tPPPLibConf C,PPPLib::tSatInfoUnit *sat_info,tSatInfoUnit *base_sat) {

        sat_info->raw_mw[0]=sat_info->raw_mw[1];

        int num_loop=C.gnssC.frq_opt==FRQ_TRIPLE?2:1;
        double P1,P2,L1,L2,lam1,lam2;
        for(int f=0;f<num_loop;f++){
            double gf=0.0,mw=0.0,w0=0.0;
            if(f==0){
                if(base_sat){
                    P1=GnssSdObs(*sat_info,*base_sat,0,GNSS_OBS_CODE);
                    P2=GnssSdObs(*sat_info,*base_sat,1,GNSS_OBS_CODE);
                    L1=GnssSdObs(*sat_info,*base_sat,0,GNSS_OBS_PHASE);
                    L2=GnssSdObs(*sat_info,*base_sat,1,GNSS_OBS_PHASE);
                    lam1=sat_info->lam[0];lam2=sat_info->lam[1];
                }
                else{
                    P1=sat_info->raw_P[0];P2=sat_info->raw_P[1];
                    L1=sat_info->raw_L[0];L2=sat_info->raw_L[1];
                    lam1=sat_info->lam[0];lam2=sat_info->lam[1];
                }
                gf=GnssObsGfComb(L1,L2,lam1,lam2);
                if(gf!=0.0) sat_info->gf[0]=gf;
                mw=GnssObsMwComb(P1,P2,L1,L2,lam1,lam2);
                if(mw==0.0) continue;
                w0=sat_info->sm_mw[0];
                sat_info->raw_mw[0]=mw;
            }
            else if(f==1){
//                if((gf=GnssObsGfComb(sat_info->raw_L[0],sat_info->raw_L[2],sat_info->lam[0],sat_info->lam[1]))!=0.0) sat_info->gf[1]=gf;
//                if((mw=GnssObsMwComb(sat_info->raw_P[0],sat_info->raw_P[2],sat_info->raw_L[0],sat_info->raw_L[2],sat_info->lam[0],sat_info->lam[2]))==0.0) continue;
//                w0=sat_info->sm_mw[1];
//                sat_info->raw_mw[1]=mw;
            }

            if(sat_info->mw_idx[f]>0){
                int j=sat_info->mw_idx[f];
                double var0=sat_info->var_mw[f];
                double var1=SQR(mw-w0)-var0;
                var1=var0+var1/j;

                sat_info->sm_mw[f]=(w0*j+mw)/(j+1);
                sat_info->mw_idx[f]++;
                sat_info->var_mw[f]=var1;
            }
            else{
                sat_info->sm_mw[f]=mw;
                sat_info->mw_idx[f]++;
                sat_info->var_mw[f]=0.25;
            }
        }
    }

    double GeoDist(Vector3d sat_pos,Vector3d rec_pos,Vector3d& sig_vec){
        double r;
        Vector3d vec=sat_pos-rec_pos;

        r=vec.norm();
        sig_vec=vec/r;
        return r;
    }

    double SatElAz(Vector3d rec_blh,Vector3d sig_vec,Vector2d& el_az){
        double el=PI/2.0,az=0.0;

        if(rec_blh[2]>-WGS84_EARTH_LONG_RADIUS){
            Vector3d enu=Xyz2Enu(rec_blh,sig_vec);
            Vector2d e(enu[0],enu[1]);
            az=e.dot(e)<1E-12?0.0:atan2(enu[0],enu[1]);
            if(az<0.0) az+=2*PI;
            el=asin(enu[2]);
        }
        el_az[0]=el;el_az[1]=az;
        return el;
    }

    double GnssMeasVar(tPPPLibConf C,GNSS_OBS obs_type,tSatInfoUnit sat_info){
        //GPS:BD2_GEO:BD2_IGSO_MEO:BD3:GAL:GLO:QZS=16:1:4:16:16:16
        double fact=obs_type==GNSS_OBS_CODE?C.gnssC.code_phase_ratio:1.0;
        int sys=sat_info.sat.sat_.sys;
        int prn=sat_info.sat.sat_.prn;

        if(sys==SYS_BDS){
            if(std::binary_search(kBD2_GEO,kBD2_GEO+NUM_BD2_GEO,prn)){
                fact*=4.0;
            }
            else if(std::binary_search(kBD2_IGSO,kBD2_IGSO+NUM_BD2_IGSO,prn)||std::binary_search(kBD2_MEO,kBD2_MEO+NUM_BD2_MEO,prn)){
                fact*=2.0;
            }
        }

#if 0
        double el=sat_info.el_az[0]*R2D;
        double var=1.0;
        if(el>=30){
            var=SQR(0.003);
        }
        else{
            var=SQR(0.003)/(4*SQR(sin(el*D2R)));
        }
#else
        double var=0.0;
        double a=C.gnssC.meas_err_factor[1],b=C.gnssC.meas_err_factor[2];
        double sin_el=sin(sat_info.el_az[0]);
        if(obs_type==GNSS_OBS_DOPPLER) a=b=10;
        var=(SQR(a)+SQR(b/sin_el));
#endif

        if(C.gnssC.ion_opt==ION_IF){
            if(C.gnssC.frq_opt==FRQ_SINGLE) fact*=0.5;
            else fact*=3.0;
        }

        if(C.mode==MODE_PPK||C.mode_opt==MODE_OPT_PPK) fact*=sqrt(2.0);

        return SQR(fact)*var;
    }

    static void AstArgs(double t, double *f) {
        static const double fc[][5]={ /* coefficients for iau 1980 nutation */
                { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470 },
                { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149 },
                { 93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417 },
                { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169 },
                { 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939 }
        };
        double tt[4];
        int i,j;

        for (tt[0]=t,i=1; i<4; i++) tt[i]=tt[i-1]*t;
        for (i=0; i<5; i++) {
            f[i]=fc[i][0]*3600.0;
            for (j=0; j<4; j++) f[i]+=fc[i][j+1]*tt[j];
            f[i]=fmod(f[i]*AS2R,2.0*PI);
        }
    }

    static void nut_iau1980(double t,const double *f,double *dpsi,double *deps){
        static const double nut[106][10]={
                { 0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9 },
                { 0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1 },
                { 0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5 },
                { 0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5 },
                { 0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1 },
                { 1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0 },
                { 0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6 },
                { 0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0 },
                { 1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1 },
                { 0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3 },
                { -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0 },
                { 0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0 },
                { -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0 },
                { 1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0 },
                { 0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0 },
                { -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0 },
                { -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0 },
                { 1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0 },
                { -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0 },
                { -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0 },
                { 0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0 },
                { 2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0 },
                { 2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0 },
                { 1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0 },
                { 0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0 },
                { 0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0 },
                { -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0 },
                { 0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0 },
                { 0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0 },
                { -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0 },
                { 0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0 },
                { 1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0 },
                { 0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0 },
                { 2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0 },
                { -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0 },
                { 1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0 },
                { 0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0 },
                { 0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0 },
                { 1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0 },
                { 0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0 },
                { -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0 },
                { 0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0 },
                { 2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0 },
                { 1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0 },
                { 1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0 },
                { 0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0 },
                { 0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0 },
                { 2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0 },
                { 1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0 },
                { 1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0 },
                { 0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0 },
                { 0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0 },
                { 1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0 },
                { 2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0 },
                { 0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0 },
                { 1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0 },
                { 1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0 },
                { -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0 },
                { 0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0 },
                { 1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0 },
                { 3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0 },
                { -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0 },
                { 1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0 },
                { -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0 },
                { 1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0 },
                { -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0 },
                { 0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0 },
                { -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0 },
                { 2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0 },
                { 3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0 },
                { 1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0 },
                { 0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0 },
                { 1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0 },
                { 1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0 },
                { 1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0 },
                { 0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0 },
                { 0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0 },
                { 0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0 },
                { 1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0 },
                { 1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0 },
                { 1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0 },
                { 1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0 },
                { 2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0 },
                { 0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0 },
                { 0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0 },
                { -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0 },
                { 2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0 },
                { 0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0 },
                { 0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0 },
                { 0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0 },
                { 0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0 },
                { 1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0 },
                { 3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0 },
                { -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0 },
                { -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0 },
                { 0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0 },
                { 0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0 },
                { -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0 },
                { 2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0 },
                { 2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0 },
                { 2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0 },
                { 2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0 },
                { 1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0 },
                { -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0 },
                { -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0 },
                { 0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0 }
        };
        double ang;
        int i,j;

        *dpsi=*deps=0.0;

        for (i=0; i<106; i++) {
            ang=0.0;
            for (j=0; j<5; j++) ang+=nut[i][j]*f[j];
            *dpsi+=(nut[i][6]+nut[i][7]*t)*sin(ang);
            *deps+=(nut[i][8]+nut[i][9]*t)*cos(ang);
        }
        *dpsi*=1E-4*AS2R; /* 0.1 mas -> rad */
        *deps*=1E-4*AS2R;
    }

    #define Rx(t,X) do { \
        (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
        (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
    } while (0)

    #define Ry(t,X) do { \
        (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
        (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
    } while (0)

    #define Rz(t,X) do { \
        (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
        (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
    } while (0)

    static void Eci2Ecef(cTime tutc,const double *erpv,Matrix3d& U,double* gmst){
        const double ep2000[]={ 2000,1,1,12,0,0 };
        static cTime tutc_;
        static double gmst_;
        static Matrix3d U_;
        cTime tgps,baseTime;
        double eps,ze,th,z,t,t2,t3,dpsi,deps,gast,f[5];
        double R1[9],R2[9],R3[9];
        Matrix3d R,P,N,W,NP;
        int i;


        if (fabs(tutc.TimeDiff(tutc_.t_))<0.01) { /* read cache */
            U=U_;
            if (gmst) *gmst=gmst_;
            return;
        }
        tutc_=tutc;

        /* terrestrial time */
        tgps=tutc_;
        tgps=tgps.Utc2Gpst();
        baseTime.Epoch2Time(ep2000);
        t=(tgps.TimeDiff(baseTime.t_)+19.0+32.184)/86400.0/36525.0;
        t2=t*t; t3=t2*t;

        /* astronomical arguments */
        AstArgs(t,f);

        /* iau 1976 precession */
        ze=(2306.2181*t+0.30188*t2+0.017998*t3)*AS2R;
        th=(2004.3109*t-0.42665*t2-0.041833*t3)*AS2R;
        z =(2306.2181*t+1.09468*t2+0.018203*t3)*AS2R;
        eps=(84381.448-46.8150*t-0.00059*t2+0.001813*t3)*AS2R;
        Rz(-z,R1); Ry(th,R2); Rz(-ze,R3);
        Matrix3d R19,R29,R39;
        R19=Map<Matrix<double,3,3>>(R1);
        R29=Map<Matrix<double,3,3>>(R2);
        R39=Map<Matrix<double,3,3>>(R3);
        R=R19*R29;
        P=R*R39;
//        MatMulPnt("NN",3,3,3,1.0,R1,R2,0.0,R);
//        MatMulPnt("NN",3,3,3,1.0,R,R3,0.0,P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */

        /* iau 1980 nutation */
        nut_iau1980(t,f,&dpsi,&deps);
        Rx(-eps-deps,R1); Rz(-dpsi,R2); Rx(eps,R3);
        R19=Map<Matrix<double,3,3>>(R1);
        R29=Map<Matrix<double,3,3>>(R2);
        R39=Map<Matrix<double,3,3>>(R3);
        R=R19*R29;
        N=R*R39;
//        MatMulPnt("NN",3,3,3,1.0,R1,R2,0.0,R);
//        MatMulPnt("NN",3,3,3,1.0,R,R3,0.0,N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */

        /* greenwich aparent sidereal time (rad) */
        gmst_=tutc_.Utc2Gmst(erpv[2]);
        gast=gmst_+dpsi*cos(eps);
        gast+=(0.00264*sin(f[4])+0.000063*sin(2.0*f[4]))*AS2R;

        /* eci to ecef transformation matrix */
        Ry(-erpv[0],R1); Rx(-erpv[1],R2); Rz(gast,R3);
        R19=Map<Matrix<double,3,3>>(R1);
        R29=Map<Matrix<double,3,3>>(R2);
        R39=Map<Matrix<double,3,3>>(R3);
        W=R19*R29;
        R=W*R39;
        NP=N*P;
        U_=R*NP;
//        MatMulPnt("NN",3,3,3,1.0,R1,R2,0.0,W);
//        MatMulPnt("NN",3,3,3,1.0,W,R3,0.0,R); /* W=Ry(-xp)*Rx(-yp) */
//        MatMulPnt("NN",3,3,3,1.0,N,P,0.0,NP);
//        MatMulPnt("NN",3,3,3,1.0,R,NP,0.0,U_); /* U=W*Rz(gast)*N*P */

        U=U_;
        if (gmst) *gmst=gmst_;
    }

    static int SunMoonPosEci(cTime tut,Vector3d& rsun,Vector3d& rmoon){
        const double ep2000[]={ 2000,1,1,12,0,0 };
        double t,f[5],eps,Ms,ls,rs,lm,pm,rm,sine,cose,sinp,cosp,sinl,cosl;
        cTime baseTime;

        baseTime.Epoch2Time(ep2000);
        t=tut.TimeDiff(baseTime.t_)/86400.0/36525.0;

        /* astronomical arguments */
        AstArgs(t,f);

        /* obliquity of the ecliptic */
        eps=23.439291-0.0130042*t;
        sine=sin(eps*D2R); cose=cos(eps*D2R);

        /* sun position in eci */
        Ms=357.5277233+35999.05034*t;
        ls=280.460+36000.770*t+1.914666471*sin(Ms*D2R)+0.019994643*sin(2.0*Ms*D2R);
        rs=AU*(1.000140612-0.016708617*cos(Ms*D2R)-0.000139589*cos(2.0*Ms*D2R));
        sinl=sin(ls*D2R); cosl=cos(ls*D2R);
        rsun[0]=rs*cosl;
        rsun[1]=rs*cose*sinl;
        rsun[2]=rs*sine*sinl;
        /* moon position in eci */
        lm=218.32+481267.883*t+6.29*sin(f[0])-1.27*sin(f[0]-2.0*f[3])+
           0.66*sin(2.0*f[3])+0.21*sin(2.0*f[0])-0.19*sin(f[1])-0.11*sin(2.0*f[2]);
        pm=5.13*sin(f[2])+0.28*sin(f[0]+f[2])-0.28*sin(f[2]-f[0])-
           0.17*sin(f[2]-2.0*f[3]);
        rm=WGS84_EARTH_LONG_RADIUS/sin((0.9508+0.0518*cos(f[0])+0.0095*cos(f[0]-2.0*f[3])+
                                        0.0078*cos(2.0*f[3])+0.0028*cos(2.0*f[0]))*D2R);
        sinl=sin(lm*D2R); cosl=cos(lm*D2R);
        sinp=sin(pm*D2R); cosp=cos(pm*D2R);
        rmoon[0]=rm*cosp*cosl;
        rmoon[1]=rm*(cose*cosp*sinl-sine*sinp);
        rmoon[2]=rm*(sine*cosp*sinl+cose*sinp);
    }

    double SunMoonPos(cTime ut1t,const double *erp_val,Vector3d& sun_pos, Vector3d& moon_pos){
        Vector3d rs,rm;
        double gmst0;
        Matrix3d U;
        cTime t=ut1t+erp_val[2];

        SunMoonPosEci(t,rs,rm);
        Eci2Ecef(t,erp_val,U,&gmst0);

        sun_pos=U*rs; moon_pos=U*rm;
        return gmst0;
    }


}