//
// Created by cc on 7/23/20.
//

#include "Solver.h"
#include "ReadFiles.h"
#include "DecodeRaw.h"

namespace PPPLib{
    cSolver::cSolver() {
        epoch_idx_=0;
    }

    cSolver::~cSolver() {}



    int cSolver::GnssObsRes(int post, tPPPLibConf C, double *x) {

    }

    bool cSolver::DetectCodeOutlier(tPPPLibConf C, int post, vector<double> omcs, vector<double> &R) {
        int m,i;
        int sat1,sat2;
        vector<int>idx;
        vector<double>vc;
        cSat sat;
        double s0,k0=ReNorm(0.95),k1=ReNorm(0.99),fact;
        bool flag=false;

        for(m=0;m<NSYS;m++){
            if(C.mode==MODE_PPK||C.mode_opt==MODE_OPT_PPK){
                for(i=0,s0=0.0;i<omcs.size();i++){
                    if(((vflag_[i]>>4)&0xF)==0) continue;

                    sat1=(vflag_[i]>>16)&0xFF;
                    sat2=(vflag_[i]>> 8)&0xFF;

                    if(previous_sat_info_[sat1-1].sat.sat_.sys_idx!=m) continue;
                    if(previous_sat_info_[sat2-1].sat.sat_.sys_idx!=m) continue;
                    s0+=omcs[i];
                    vc.push_back(omcs[i]);
                    idx.push_back(i);
                }
            }
            else if(C.mode==MODE_PPP||C.mode_opt==MODE_OPT_PPP){
                for(i=0,s0=0.0;i<omcs.size();i++){
                    if(((vflag_[i]>>4)&0xF)!=GNSS_OBS_CODE) continue;

                    sat1=(vflag_[i]>> 8)&0xFF;

                    if(previous_sat_info_[sat1-1].sat.sat_.sys_idx!=m) continue;
                    s0+=omcs[i];
                    vc.push_back(omcs[i]);
                    idx.push_back(i);
                }
            }

            if(idx.size()>2){
#if 0
                s0/=idx.size();
                for(i=0;i<idx.size();i++) vc[i]-=s0;
#endif
                MatMul("NT",1,1,idx.size(),1.0/(idx.size()-1),vc.data(),vc.data(),0.0,&s0);

                for(i=0;i<idx.size();i++){
                    sat1=vflag_[idx[i]]>>8&0xFF;
                    sat=cSat(sat1);
                    sat.SatNo2Id();
                    double a=fabs(vc[i])/SQRT(s0);
                    if(fabs(vc[i])/SQRT(s0)>=k1){
                        // zero-weight segment
                        LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<sat.sat_.id<<" PSEUDORANGE RESIDUAL OUTLIER DETECTED, outlier="<<omcs[idx[i]]<<" ZERO-WEIGHT SEGMENT, factor="<<100.0;
                        flag=true;
                    }
                    else if(fabs(vc[i])/SQRT(s0)>=k0){
                        // weight-reduced segment
                        fact=k0/fabs(vc[i]/SQRT(s0))*SQR((k1-fabs(vc[i]/SQRT(s0)))/(k1-k0));
//                        previous_sat_info_[sat1-1].var_factor=1.0/fact;
                        LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<sat.sat_.id<<" PSEUDORANGE RESIDUAL OUTLIER DETECTED, outlier="<<omcs[idx[i]]<<" WEIGHT-REDUCED SEGMENT, factor="<<1.0/fact;
                        flag=true;
                    }
                    else{
                        // weight-reserved segment
//                        previous_sat_info_[sat1-1].var_factor=1.0;
                    }
                }
            }
            idx.clear();vc.clear();
        }
        return flag;
    }

    bool cSolver::IGG3(tPPPLibConf C, int post, vector<double> omcs,vector<double> R) {
//        double s0;
//        MatMul("NT",1,1,omcs.size(),1.0/(omcs.size()-1),omcs.data(),omcs.data(),0.0,&s0);
//
//        int i, sat1;
//        cSat sat;
//        double k0=ReNorm(0.95),k1=ReNorm(0.99),fact;
//        bool flag=false;
//        double obs_type;
//        int ppp=0;
//
//        for(i=0;i<omcs.size();i++){
//            sat1=vflag_[i]>>8&0xFF;
//            sat=cSat(sat1);
//            sat.SatNo2Id();
//            double a=fabs(omcs[i])/SQRT(R[i]);
//
//            if(C.mode==MODE_PPP||C.mode_opt==MODE_OPT_PPP){
//                ppp=1;
//                obs_type=vflag_[i]>>4&0xF;
//            }
//            else{
//
//            }
//
//            if(a>=k1){
//                // zero-weight segment
//                if(obs_type==GNSS_OBS_CODE) previous_sat_info_[sat1-1].c_var_factor[]=1000;
//                else if(obs_type==GNSS_OBS_PHASE) previous_sat_info_[sat1-1].p_var_factor=1000;
//                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<sat.sat_.id<<" "<<(ppp?(obs_type==0?"P":"L"):"RTK")<<" RESIDUAL OUTLIER DETECTED, outlier="<<omcs[i]<<" ZERO-WEIGHT SEGMENT, factor="<<100.0;
//                flag=true;
//            }
//            else if(a>=k0){
//                // weight-reduced segment
//                fact=(k0/fabs(a))*SQR((k1-a)/(k1-k0));
//
//                if(obs_type==GNSS_OBS_CODE) previous_sat_info_[sat1-1].c_var_factor*=1.0/fact;
//                else if(obs_type==GNSS_OBS_PHASE) previous_sat_info_[sat1-1].p_var_factor*=1.0/fact;
//
//                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<sat.sat_.id<<" "<<(ppp?(obs_type==0?"P":"L"):"RTK")<<" RESIDUAL OUTLIER DETECTED, outlier="<<omcs[i]<<" WEIGHT-REDUCED SEGMENT, factor="<<1.0/fact;
//                flag=true;
//            }
//            else{
//                // weight-reserved segment
//                if(obs_type==GNSS_OBS_CODE) previous_sat_info_[sat1-1].c_var_factor=1.0;
//                else if(obs_type==GNSS_OBS_PHASE) previous_sat_info_[sat1-1].p_var_factor=1.0;
//            }
//        }
//
//        return flag;
    }


    int cSolver::GetSingalInd(tSatObsUnit& sat_obs,int *f){
        const tGnssSignal *ps= nullptr;
        int n,p;

        if(sat_obs.L[*f]==0.0||sat_obs.P[*f]==0.0){
            ps=rover_obs_.signal_[sat_obs.sat.sat_.sys_idx];
            for(n=0,p=-1;n<ps->n;n++){
                if(ps->frq[n]==*f+1&&ps->pri[n]&&(p<0||ps->pri[n]>ps->pri[p])&&sat_obs.L[ps->pos[n]]&&sat_obs.P[ps->pos[n]]) p=n;
            }
            if(p<0) return 0;
            *f=ps->pos[p];
            return 1;
        }
        else return 0;
    }

    int cSolver::AdjustObsInd(tSatObsUnit& sat_obs,int *i, int *j, int *k) {
        int info=0;
        info|=i&&GetSingalInd(sat_obs,i);
        info|=j&&GetSingalInd(sat_obs,j);
        info|=k&&GetSingalInd(sat_obs,k);
        return info;
    }

    void cSolver::InitEpochSatInfo(vector<tSatInfoUnit> &sat_infos) {
        for(int j=0;j<MAX_SAT_NUM;j++){
            for(int k=0;k<MAX_GNSS_USED_FRQ_NUM;k++){
                if(epoch_idx_==1) {
                    previous_sat_info_[j].stat=SAT_NO_USE;
                    previous_sat_info_[j].rejc[k]=0;
                }
                previous_sat_info_[j].vsat[k]=0;
                previous_sat_info_[j].outc[k]++;
                previous_sat_info_[j].p_var_factor[k]=previous_sat_info_[j].c_var_factor[k]=1.0;
            }
        }

        for(int j=0;j<sat_infos.size();j++){
            for(int k=0;k<MAX_GNSS_USED_FRQ_NUM;k++){
                sat_infos.at(j).outc[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].outc[k];
                sat_infos.at(j).lock[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].lock[k];
                sat_infos.at(j).rejc[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].rejc[k];
            }
            for(int k=0;k<2;k++){
                sat_infos.at(j).sm_mw[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].sm_mw[k];
                sat_infos.at(j).mw_idx[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].mw_idx[k];
                sat_infos.at(j).var_mw[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].var_mw[k];
                sat_infos.at(j).gf[k]=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].gf[k];
                sat_infos.at(j).phase_wp=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].phase_wp;
            }
            sat_infos.at(j).lc_amb=previous_sat_info_[sat_infos.at(j).sat.sat_.no-1].lc_amb;
        }
    }

    void cSolver::UpdateSatInfo(vector<tSatInfoUnit> &sat_infos) {
        for(int i=0;i<sat_infos.size();i++){
            previous_sat_info_[sat_infos.at(i).sat.sat_.no-1]=sat_infos.at(i);
        }
    }

    void cSolver::UpdateGnssObs(tPPPLibConf C,tEpochSatUnit& epoch_sat_obs,RECEIVER_INDEX rec) {
        int sys,i,j,k;
        int *frqs;
        int f1,f2,f3;

        for(i=0;i<epoch_sat_obs.sat_num;i++){

            tSatInfoUnit sat_info={0};
            sys=epoch_sat_obs.epoch_data.at(i).sat.sat_.sys;
            sat_info.t_tag=epoch_sat_obs.obs_time;
            sat_info.sat=epoch_sat_obs.epoch_data.at(i).sat;
            if(sat_info.sat.sat_.no==4) continue;
            if(sat_info.sat.sat_.no==50) continue;
            sat_info.stat=SAT_USED;
            switch(sys){
                case SYS_BDS:
                    if(sat_info.sat.sat_.prn>18){
                        frqs=C.gnssC.gnss_frq[NSYS];break;
                    }
                    frqs=C.gnssC.gnss_frq[SYS_INDEX_BDS];break;
                case SYS_GAL: frqs=C.gnssC.gnss_frq[SYS_INDEX_GAL];break;
                case SYS_GLO: frqs=C.gnssC.gnss_frq[SYS_INDEX_GLO];break;
                case SYS_QZS: frqs=C.gnssC.gnss_frq[SYS_INDEX_QZS];break;
                default:
                    frqs=C.gnssC.gnss_frq[SYS_INDEX_GPS];break;
            }

            //first frequency
            f1=frqs[0];
            if(C.gnssC.adj_obs) AdjustObsInd(epoch_sat_obs.epoch_data.at(i),&f1,&f2,&f3);
            if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),0,frqs[0],f1,nav_.glo_frq_num);
            else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),0,frqs[0],f1,nullptr);
            sat_info.c_var_factor[0]=1.0;

            //second frequency
            if(C.gnssC.frq_opt==FRQ_DUAL){
                f2=frqs[1];
                if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nav_.glo_frq_num);
                else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nullptr);
                sat_info.c_var_factor[1]=1.0;
            }
            // third frequency
            else if(C.gnssC.frq_opt==FRQ_TRIPLE){
                f2=frqs[1];
                if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nav_.glo_frq_num);
                else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),1,frqs[1],f2,nullptr);
                sat_info.c_var_factor[1]=1.0;

                f3=frqs[2];
                if(sat_info.sat.sat_.sys==SYS_GLO) gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),2,frqs[2],f3,nav_.glo_frq_num);
                else gnss_obs_operator_.ReAlignObs(C,sat_info,epoch_sat_obs.epoch_data.at(i),2,frqs[2],f3,nullptr);
                sat_info.c_var_factor[2]=1.0;
            }

            rec==REC_ROVER?epoch_sat_info_collect_.push_back(sat_info):base_sat_info_collect_.push_back(sat_info);
        }
    }

    void cSolver::CorrGnssObs(tPPPLibConf C,Vector3d& rr) {


        tSatInfoUnit *sat_info= nullptr;
        double pcv_dants[MAX_GNSS_USED_FRQ_NUM]={0},dantr[MAX_GNSS_USED_FRQ_NUM]={0};
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);

            gnss_err_corr_.cbia_model_.InitSatInfo(sat_info, nullptr);
            gnss_err_corr_.cbia_model_.GetCodeBias();
            gnss_err_corr_.cbia_model_.UpdateSatInfo();

            if(C.mode==MODE_PPP){
                if(sat_info->el_az[0]*R2D<C.gnssC.ele_min) continue;
                if(C.gnssC.sat_pcv) gnss_err_corr_.ant_model_.SatPcvCorr(sat_info,rr, pcv_dants);
                if(C.gnssC.rec_ant) gnss_err_corr_.ant_model_.RecAntCorr(sat_info, dantr,REC_ROVER);
                gnss_err_corr_.PhaseWindUp(*sat_info,rr);
            }

            for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++){
                if(sat_info->raw_P[j]==0.0) continue;
                sat_info->cor_P[j]=sat_info->raw_P[j]-sat_info->code_bias[j]-dantr[j]-pcv_dants[j];
                if(sat_info->raw_L[j]==0.0) continue;
                sat_info->cor_L[j]=sat_info->raw_L[j]*sat_info->lam[j]-dantr[j]-pcv_dants[j]-sat_info->phase_wp*sat_info->lam[j];
            }

            if(sat_info->sat.sat_.sys==SYS_BDS&&sat_info->sat.sat_.prn<=18&&(C.mode==MODE_PPP||C.mode_opt==MODE_OPT_PPP)){
                gnss_err_corr_.BD2MultipathModel(sat_info);
                for(int j=0;j<3;j++){
                    if(sat_info->raw_P[j]==0.0) continue;
                    sat_info->cor_P[j]-=sat_info->bd2_mp[j];
                }
            }
            if(C.gnssC.ion_opt==ION_IF||C.gnssC.ion_opt==ION_IF_DUAL)
                gnss_obs_operator_.MakeGnssObsComb(C,COMB_IF,sat_info,previous_sat_info_[sat_info->sat.sat_.no-1]);
        }
    }

    bool cSolver::InitReader(tPPPLibConf C) {
        if(!C.fileC.rover.empty()){
            cReadGnssObs rover_reader(C.fileC.rover,nav_,rover_obs_,REC_ROVER);
            rover_reader.SetGnssSysMask(C.gnssC.nav_sys);
            rover_reader.Reading();
        }

        if(!C.fileC.brd.empty()){
            cReadGnssBrdEph brd_reader(C.fileC.brd, nav_);
            brd_reader.SetGnssSysMask(C.gnssC.nav_sys);
            brd_reader.Reading();
        }

        if(!C.fileC.cbias.empty()){
            cReadGnssCodeBias cbias_reader(C.fileC.cbias, nav_);
            cbias_reader.Reading();
        }
        return true;
    }


    void cSolver::InitSolver(tPPPLibConf C) {
        para_=cParSetting(C);
    }

    bool cSolver::SolverProcess(tPPPLibConf C,int idx) {}

    bool cSolver::SolverEpoch() {}

    bool cSolver::Estimator(tPPPLibConf C) {}

    bool cSolver::SolutionUpdate() {}

    static void InitP(double unc,double unc0,int row_s,int col_s,MatrixXd& P){
        double q=unc==0.0?SQR(unc0):SQR(unc);
        Vector3d vec(q,q,q);
        P.block<3,3>(row_s,col_s)=vec.asDiagonal();
    }

    void cSolver::InitFullPx(tPPPLibConf C) {
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
        Vector3d vec(0,0,0);
        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        int ibg=para_.IndexBg();
        if(para_.NumPos()>0){
            InitP(C.insC.init_pos_unc,UNC_POS,ip,ip,full_Px_);
        }
        if(para_.NumVel()>0){
            InitP(C.insC.init_vel_unc,UNC_VEL,iv,iv,full_Px_);
        }
        if(para_.NumAtt()>0){
            InitP(C.insC.init_att_unc,UNC_ATT,ia,ia,full_Px_);
        }
        if(para_.NumBa()>0){
            InitP(C.insC.init_ba_unc,UNC_BA,iba,iba,full_Px_);
        }
        if(para_.NumBg()>0){
            InitP(C.insC.init_bg_unc,UNC_BG,ibg,ibg,full_Px_);
        }
    }

    void cSolver::InitX(double xi, double var, int idx,double *x,double *p) {
        x[idx]=xi;
        for(int j=0;j<num_full_x_;j++){
//            full_Px_(idx,j)=full_Px_(j,idx)=idx==j?var:0.0;
            p[idx+j*num_full_x_]=p[j+idx*num_full_x_]=idx==j?var:0.0;
        }
    }

    static void SetPsd(double psd,double dt,int row_s,int col_s,MatrixXd& Q){
        Vector3d vec(psd*dt,psd*dt,psd*dt);
        Q.block<3,3>(row_s,col_s)=vec.asDiagonal();
    }

    Eigen::MatrixXd cSolver::InitQ(tPPPLibConf C,double dt) {
        Eigen::MatrixXd Q;
        Q=MatrixXd::Zero(num_full_x_,num_full_x_);

        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        int ibg=para_.IndexBg();
        SetPsd(C.insC.psd_acce,dt,iv,iv,Q);
        SetPsd(C.insC.psd_gyro,dt,ia,ia,Q);
        SetPsd(C.insC.psd_bg,dt,ibg,ibg,Q);
        SetPsd(C.insC.psd_ba,dt,iba,iba,Q);

        return Q;
    }

    cSppSolver::cSppSolver() {
        spp_conf_.mode=MODE_SPP;
        spp_conf_.gnssC.ele_min=10.0;
        spp_conf_.gnssC.frq_opt=FRQ_SINGLE;
        spp_conf_.gnssC.trp_opt=TRP_SAAS;
        spp_conf_.gnssC.ion_opt=ION_KLB;
        spp_conf_.estimator=SOLVE_LSQ;

        para_=cParSetting();
        num_full_x_=para_.GetPPPLibPar(spp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
    }

    cSppSolver::cSppSolver(tPPPLibConf conf) {
        spp_conf_=conf;
        spp_conf_.mode=MODE_SPP;
        spp_conf_.mode_opt=MODE_OPT_KINEMATIC;
        para_=cParSetting(spp_conf_);
        num_full_x_=para_.GetPPPLibPar(spp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
    }

    cSppSolver::~cSppSolver() {}

    void cSppSolver::CorrDoppler(tSatInfoUnit &sat_info, Vector3d &rover_xyz,int f) {
        tSatInfoUnit sat_info0;
        Vector3d rover_blh=Xyz2Blh(rover_xyz);
        Vector3d sat_pos=sat_info.brd_pos,sat_vel=sat_info.brd_vel;
        sat_info0.brd_pos=sat_pos-sat_vel*1.0;
        Vector3d sig_vec0;
        GeoDist(sat_info0.brd_pos,rover_xyz,sig_vec0);
        SatElAz(rover_blh,sig_vec0,sat_info0.el_az);
        gnss_err_corr_.trp_model_.InitSatInfo(&sat_info0,&rover_blh);
        gnss_err_corr_.trp_model_.GetSaasTrp(0.7, nullptr, nullptr);
        gnss_err_corr_.trp_model_.UpdateSatInfo();
        gnss_err_corr_.ion_model_.InitSatInfo(&sat_info0,&rover_blh);
        gnss_err_corr_.ion_model_.UpdateSatInfo();
        double trp_delta=((sat_info.trp_dry_delay[0]+sat_info.trp_wet_delay[0])-(sat_info0.trp_dry_delay[0]+sat_info0.trp_wet_delay[0])/1.0);
        double ion_delta=((sat_info.ion_delay[0]-sat_info0.ion_delay[0])/1.0);
        sat_info.cor_D[f]=sat_info.raw_D[f];
    }

    Vector3d cSppSolver::SatVelSagnacCorr(const Vector3d &sat_vel, const double tau) {
        Vector3d vel_corr(0,0,0);
        double a=OMGE_GPS*tau;

        vel_corr[0]=sat_vel[0]+a*sat_vel[1];
        vel_corr[1]=sat_vel[1]-a*sat_vel[0];
        vel_corr[2]=sat_vel[2];
        return vel_corr;
    }

    int cSppSolver::DopplerRes(tPPPLibConf C,MatrixXd& H_mat, MatrixXd& R_mat,VectorXd& L,VectorXd& x,Vector3d rover_xyz) {
        vector<double>H,omcs,meas_var_vec;
        int sys,sys_mask[NSYS]={0},num_doppler=0;
        double omc,r,meas_var;
        tSatInfoUnit* sat_info= nullptr;
        Vector3d rover_blh=Xyz2Blh(rover_xyz),sig_vec;
        Matrix3d Cen=CalcCen(rover_blh,COORD_ENU);

        int num_used_frq=para_.GetGnssUsedFrqs();

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->sat.sat_.no==4) continue;
            if(sat_info->sat.sat_.no==55) continue;
            sys=sat_info->sat.sat_.sys;
//            if(sys==SYS_BDS&&(1<=sat_info->sat.sat_.prn&&sat_info->sat.sat_.prn<=5||sat_info->sat.sat_.prn>18)) continue;
//            if(sys!=SYS_GPS) continue; // only use GPS dopppler for velocity estimate
            if(sat_info->stat!=SAT_USED) continue;

            LOG_IF(i==0,DEBUG)<<"DOPPLER RESIDUAL"<<"("<<epoch_idx_<<")"<<": "<<sat_info->t_tag.GetTimeStr(1)
                              <<setprecision(12)<<" "<<"VX="<<x[0]<<" "<<"VY="<<x[1]<<" "<<"VZ="<<x[2];

            for(int f=0;f<num_used_frq;f++){
                double doppler=sat_info->raw_D[f];
                double lam=sat_info->lam[f];
                if(doppler==0.0||lam==0.0) continue;
                CorrDoppler(*sat_info,rover_xyz,f);
                Vector3d sat_pos=sat_info->brd_pos;
                Vector3d sat_vel=sat_info->brd_vel;
                Vector3d s_r=(sat_pos-rover_xyz);
                double tau=s_r.norm()/CLIGHT;
                Vector3d sat_vel_corr=SatVelSagnacCorr(sat_vel,tau);
                Vector3d rover_vel(x[0],x[1],x[2]);
                Vector3d e=s_r/s_r.norm();
                Vector3d relative_vel=sat_vel_corr-rover_vel;
                double rate=e.dot(relative_vel);

#if 0
                double cos_el=cos(sat_info->el_az[0]);
                Vector3d a(0.0,0.0,0.0);
                a[0]=sin(sat_info->el_az[1])*cos_el;
                a[1]=cos(sat_info->el_az[1])*cos_el;
                a[2]=sin(sat_info->el_az[0]);
                Vector3d e=Cen.transpose()*a;

                int idx_vel=para_.IndexPos();
                Vector3d vs(0,0,0);
                for(int j=idx_vel;j<idx_vel+3;j++) vs[j-idx_vel]=sat_info->brd_vel[j-idx_vel]-x[j];
                double rate=vs.dot(e)+OMGE_GPS/CLIGHT*(sat_info->brd_vel[1]*rover_xyz[0]+sat_info->brd_pos[1]*x[idx_vel]
                                                       -sat_info->brd_vel[0]*rover_xyz[1]-sat_info->brd_pos[0]*x[idx_vel+1]);
#endif

                int idx_clk_drift=para_.IndexClk(SYS_INDEX_GPS);
                double a=CLIGHT*sat_info->brd_clk[1];
                omc=lam*sat_info->cor_D[f]-(rate+x[idx_clk_drift]-CLIGHT*sat_info->brd_clk[1]);
                meas_var=GnssMeasVar(C,GNSS_OBS_DOPPLER,*sat_info);

                omcs.push_back(omc);
                meas_var_vec.push_back(meas_var);

                for(int j=0;j<num_full_x_;j++) H.push_back(j<3?-e[j]:(j==idx_clk_drift)?1.0:0.0);
                int idx;
                if(sys==SYS_BDS) {idx=idx_clk_drift+SYS_INDEX_BDS;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_BDS]++;}
                else if(sys==SYS_GAL) {idx=idx_clk_drift+SYS_INDEX_GAL;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_GAL]++;}
                else if(sys==SYS_GLO) {idx=idx_clk_drift+SYS_INDEX_GLO;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_GLO]++;}
                else if(sys==SYS_QZS) {idx=idx_clk_drift+SYS_INDEX_QZS;omc-=x[idx];omcs[num_doppler]-=x[idx];H[idx+num_doppler*num_full_x_]=1.0;sys_mask[SYS_INDEX_QZS]++;}
                else sys_mask[SYS_INDEX_GPS]++;

                num_doppler++;

                char buff[MAX_BUFF]={'\0'};
                sprintf(buff,"%s omc=%12.4f, var=%7.3f el=%3.1f clk_drift=%14.3f sat_clk_drift=%14.3f doppler=%10.3f rate=%10.3f",
                        sat_info->sat.sat_.id.c_str(),omc,meas_var,sat_info->el_az[0]*R2D,sys==SYS_GPS?x[idx_clk_drift]:x[idx],sat_info->brd_clk[1]*CLIGHT,-lam*doppler,rate);
                LOG(DEBUG)<<buff;
            }
        }

        int idx_clk=para_.IndexClk(SYS_INDEX_GPS);
        for(int i=0;i<NSYS;i++){
            if(sys_mask[i]) continue;
            omcs.push_back(0.0);
            for(int j=0;j<num_full_x_;j++) H.push_back(j==i+idx_clk?1.0:0.0);
            meas_var_vec.push_back(0.01);
            num_doppler++;
        }

        L=Map<VectorXd>(omcs.data(),num_doppler);
        H_mat=Map<MatrixXd>(H.data(),num_full_x_,num_doppler);
        Map<VectorXd> var_vec(meas_var_vec.data(),num_doppler);
        R_mat=var_vec.asDiagonal();
        H.clear();omcs.clear();meas_var_vec.clear();

        return num_doppler;
    }

    void cSppSolver::EstDopVel(Vector3d rover_xyz) {
        VectorXd L,x(num_full_x_);
        MatrixXd H,R,Q;
        int nv,stat=0;

        for(int i=0;i<num_full_x_;i++) x[i]=0.0;
        for(int i=0;i<iter_;i++){
            if((nv=DopplerRes(spp_conf_,H,R,L,x,rover_xyz))<4){
                break;
            }

//            cout<<H.transpose()<<endl;
//            cout<<R<<endl;
//            cout<<L<<endl;
            stat=lsq_.Adjustment(L,H,R,x,Q,nv,num_full_x_);

            if(stat){
                for(int j=0;j<3;j++) ppplib_sol_.vel[j]=x[j];
                break;
            }
        }
    }

    double cSppSolver::Dops() {
        tSatInfoUnit *sat_info=nullptr;
        MatrixXd H_mat(4,num_valid_sat_),Q,P;
        vector<double>H(4*num_valid_sat_,0.0);
        int n=0;

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat!=SAT_USED) continue;
            double cos_el=cos(sat_info->el_az[0]);
            double sin_el=sin(sat_info->el_az[0]);
            H[4*n]=cos_el*sin(sat_info->el_az[1]);
            H[1+4*n]=cos_el*cos(sat_info->el_az[1]);
            H[2+4*n]=sin_el;
            H[3+4*n++]=1.0;
        }
        if(n<4) return 0.0;
        H_mat=Map<MatrixXd>(H.data(),4,n);
        Q=H_mat*H_mat.transpose();
        P=Q.inverse();
        ppplib_sol_.dops[0]=SQRT(P.trace());
        ppplib_sol_.dops[1]=SQRT(P(0,0)+P(1,1)+P(2,2));
        ppplib_sol_.dops[2]=SQRT(P(0,0)+P(1,1));
        ppplib_sol_.dops[2]=SQRT(P(2,2));

        return ppplib_sol_.dops[1];
    }

    bool cSppSolver::ValidateSol(tPPPLibConf C) {
        tSatInfoUnit* sat_info= nullptr;

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat==SAT_USED&&sat_info->el_az[0]<C.gnssC.ele_min*D2R){
                sat_info->stat=SAT_LOW_EL;
                continue;
            }
        }

        string buff;
        int num_no_used=0;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat!=SAT_USED){
                num_no_used++;
                if(num_no_used==1) buff=epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)+" SATELLITE NO USED: "+sat_info->sat.sat_.id+"("+kGnssSatStatStr[sat_info->stat+1]+") ";
                else{
                    buff+=sat_info->sat.sat_.id+"("+kGnssSatStatStr[sat_info->stat+1]+") ";
                }
            }
        }
        LOG_IF(num_no_used>0,DEBUG)<<buff;

        double vv=lsq_.unit_weight_STD_*((num_L_>num_full_x_)?(num_L_-num_full_x_):num_L_);
        if(num_L_>num_full_x_&&vv>kChiSqr[num_L_-num_full_x_-1]){
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" SPP SOLVE FAILED, CHI-SQR OVERRUN vv: "<<vv<<" THRESHOLD: "<<kChiSqr[num_L_-num_full_x_-1];
            return false;
        }

        double pdop=Dops();

        if(pdop>C.gnssC.max_pdop){
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" SPP SOLVE FAILED, PDOP OVERRUN pdop: "<<pdop<<" THRESHOLD: "<<C.gnssC.max_pdop;
            return false;
        }

        return true;
    }

    bool cSppSolver::PostResidualQc(vector<double>omcs,vector<double>R) {
        // step_1: where bad epoch solution, unit weight std>0.15
        // step_2: who is bad satellite
        // step_3:
        bool flag=false;
        int i,sat_idx,f;
        vector<double>v,avg_v,fabs_v,norm_v;
        double s0=0.0,el,k0=ReNorm(0.95),k1=ReNorm(0.99);

        for(i=0;i<vflag_.size();i++){
            v.push_back(omcs[i]);
            fabs_v.push_back(fabs(omcs[i]));
            s0+=v[i];
        }

        if(lsq_.unit_weight_STD_>0.12){
            auto max_v=max_element(begin(fabs_v),end(fabs_v));
            int idx_max_v=distance(begin(fabs_v),max_v);
            sat_idx=vflag_[idx_max_v]>>4&0xF;
            epoch_sat_info_collect_[sat_idx].stat=SAT_NO_USE;
            flag=true;
        }
        else{
            for(i=0;i<v.size();i++) avg_v.push_back(v[i]-s0/v.size());
            MatMul("NT",1,1,v.size(),1.0/(v.size()-1),avg_v.data(),avg_v.data(),0.0,&s0);
            for(i=0;i<v.size();i++) {
                sat_idx=vflag_[i]>>4&0xF;
                el=epoch_sat_info_collect_[sat_idx].el_az[0];
                if(epoch_sat_info_collect_[sat_idx].c_var_factor[0]!=1.0) continue;
                norm_v.push_back(fabs_v[i]/sin(el));
            }

            auto max_norm_v=max_element(begin(norm_v),end(norm_v));
            int idx_max_norm_v=distance(begin(norm_v),max_norm_v);
            sat_idx=vflag_[idx_max_norm_v]>>4&0xF;
            if(*max_norm_v>k1){
                epoch_sat_info_collect_[sat_idx].stat=SAT_NO_USE;
                flag=true;
            }
            else if(*max_norm_v>k0){
                double fact=(k0/(*max_norm_v))*SQR((k1-(*max_norm_v))/(k1-k0));
                epoch_sat_info_collect_[sat_idx].c_var_factor[0]=1.0/fact;
                flag=true;
            }
        }

        v.clear();avg_v.clear();fabs_v.clear(),norm_v.clear();
        return flag;
    }

    int cSppSolver::GnssObsRes(int post, tPPPLibConf C,double* x) {
        num_L_=0;
        num_valid_sat_=0;
        if(vflag_.size()) vflag_.clear();
        for(int i=0;i<NSYS;i++) sys_mask_[i]=0;

        vector<double>H,omcs,meas_var_vec;
        int sys;
        double omc,r,meas_var;
        bool have_large_res=false;
        tSatInfoUnit* sat_info= nullptr;
        Vector3d rover_xyz(x[0],x[1],x[2]);
        Vector3d rover_blh=Xyz2Blh(rover_xyz),sig_vec;

        int num_used_frq=para_.GetGnssUsedFrqs();
        if(num_used_frq<=0) return 0;

        vector<int>sat_idx;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);

            LOG_IF(i==0,DEBUG)<<"SINGLE POINT POSITION RESIDUAL"<<"("<<epoch_idx_<<"-"<<post<<")"<<": "<<sat_info->t_tag.GetTimeStr(1)
                                       <<setprecision(12)<<" "<<"X="<<rover_xyz[0]<<" "<<"Y="<<rover_xyz[1]<<" "<<"Z="<<rover_xyz[2];

            sys=sat_info->sat.sat_.sys;
            if(sat_info->stat!=SAT_USED) continue;

            for(int f=0;f<num_used_frq;f++){
                double cor_P=0.0;
                if(spp_conf_.gnssC.ion_opt==ION_IF){
                    cor_P=sat_info->cor_if_P[0];
                }
                else{
                    cor_P=sat_info->cor_P[f];
                }

                if(cor_P==0.0){
                    sat_info->stat=SAT_NO_CODE;
                    continue;
                }

                r=GeoDist(sat_info->brd_pos,rover_xyz,sig_vec);
                if(SatElAz(rover_blh,sig_vec,sat_info->el_az)<=0.0){
                    continue;
                }
                double el=sat_info->el_az[0]*R2D;
                if(sat_info->el_az[0]*R2D<C.gnssC.ele_min){
                    continue;
                }

#if 0
                if(sat_info->sat.sat_.sys==SYS_BDS&&sat_info->sat.sat_.prn<19){
                    gnss_err_corr_.BD2MultipathModel(sat_info);
                }
#endif
                gnss_err_corr_.ion_model_.InitSatInfo(sat_info,&rover_blh);
                gnss_err_corr_.ion_model_.GetIonError();
                gnss_err_corr_.ion_model_.UpdateSatInfo();
                gnss_err_corr_.trp_model_.InitSatInfo(sat_info,&rover_blh);
                gnss_err_corr_.trp_model_.GetSaasTrp(0.7, nullptr, nullptr);
                gnss_err_corr_.trp_model_.UpdateSatInfo();

                double sag_err=gnss_err_corr_.SagnacCorr(sat_info->brd_pos,rover_xyz);
                double rec_clk=x[para_.IndexClk(SYS_INDEX_GPS)];
                double sat_clk=sat_info->brd_clk[0]*CLIGHT;
                double trp_err=sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0];
                double ion_err=sat_info->ion_delay[0];
                double cod_bia=sat_info->code_bias[f];

                omc=cor_P-(r+sag_err+(rec_clk-CLIGHT*sat_info->brd_clk[0])+sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0]+sat_info->ion_delay[0]);
                meas_var=GnssMeasVar(C,GNSS_OBS_CODE,*sat_info)+sat_info->brd_eph_var+sat_info->trp_var+sat_info->ion_var;
                meas_var*=sat_info->c_var_factor[f];
                omcs.push_back(omc);
                meas_var_vec.push_back(meas_var);

                int idx_clk=para_.IndexClk(SYS_INDEX_GPS),idx=0;
                for(int j=0;j<num_full_x_;j++) H.push_back(j<3?-sig_vec[j]:(j==idx_clk)?1.0:0.0);
                if(sys==SYS_BDS) {idx=idx_clk+SYS_INDEX_BDS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_BDS]++;}
                else if(sys==SYS_GAL) {idx=idx_clk+SYS_INDEX_GAL;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_GAL]++;}
                else if(sys==SYS_GLO) {idx=idx_clk+SYS_INDEX_GLO;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_GLO]++;}
                else if(sys==SYS_QZS) {idx=idx_clk+SYS_INDEX_QZS;omc-=x[idx];omcs[num_L_]-=x[idx];H[idx+num_L_*num_full_x_]=1.0;sys_mask_[SYS_INDEX_QZS]++;}
                else sys_mask_[SYS_INDEX_GPS]++;
                num_L_++;

                if(post==0) sat_info->prior_res_P[f]=omc;

                char buff[MAX_BUFF]={'\0'};
                sprintf(buff,"%s omc=%12.4f, var=%7.3f el=%3.1f dtr=%14.3f dts=%14.3f trp=%6.3f ion=%6.3f cbias=%6.3f obs=%14.3f range=%14.3f",
                        sat_info->sat.sat_.id.c_str(),omc,meas_var,sat_info->el_az[0]*R2D,sys==SYS_GPS?rec_clk:x[idx],sat_clk,trp_err,ion_err,cod_bia,cor_P,r+sag_err);
                LOG(DEBUG)<<buff;

                if(f==0){
                    sat_idx.push_back(i);
                    num_valid_sat_++;
                }

                vflag_.push_back((sat_info->sat.sat_.no<<8)|(i<<4)|(f));
            }
        }

        if(epoch_idx_>10&&post&&spp_conf_.gnssC.res_qc){
            have_large_res=PostResidualQc(omcs,meas_var_vec);
        }

        int idx_clk=para_.IndexClk(SYS_INDEX_GPS);
        for(int i=0;i<NSYS;i++){
            if(sys_mask_[i]) continue;
            omcs.push_back(0.0);
            for(int j=0;j<num_full_x_;j++) H.push_back(j==i+idx_clk?1.0:0.0);
            meas_var_vec.push_back(0.01);
            num_L_++;
        }

        omc_L_=Map<VectorXd>(omcs.data(),num_L_);
        H_=Map<MatrixXd>(H.data(),num_full_x_,num_L_);
        Map<VectorXd> var_vec(meas_var_vec.data(),num_L_);
        R_=var_vec.asDiagonal();
        H.clear();omcs.clear();meas_var_vec.clear();

        return post==1?have_large_res:num_valid_sat_<=0;
    }

    bool cSppSolver::Estimator(tPPPLibConf C) {
        bool stat;
        bool valid_flag=false;
        ppplib_sol_.stat=SOL_NONE;
        tSatInfoUnit* sat_info= nullptr;

        VectorXd x=full_x_;
        MatrixXd Px=full_Px_;
        for(int i=0;i<iter_;i++){

            if(GnssObsRes(0,spp_conf_,x.data())) continue;
            if(num_valid_sat_<4){
                LOG(WARNING)<<"SPP NO ENOUGH VALID SATELLITE "<<num_valid_sat_;
                return false;
            }
            ppplib_sol_.valid_sat_num=num_valid_sat_;
            stat=lsq_.Adjustment(omc_L_,H_,R_,x,Px,num_L_,num_full_x_);

            if(GnssObsRes(i+1,spp_conf_,x.data())){
                x=full_x_;
                continue;
            }

            if(stat){
                full_x_=x;
                full_Px_=Px;
                valid_flag = ValidateSol(C);
                if(valid_flag){
                    ppplib_sol_.stat=SOL_SPP;
                    ppplib_sol_.sigma=lsq_.unit_weight_STD_;
                }
                else{
                    ppplib_sol_.stat=SOL_NONE;
                }
                break;
            }
            if(i>=iter_){
                LOG(WARNING)<<"SPP ITER OVERRUN";
            }

        }

        if(valid_flag&&C.gnssC.use_doppler){
            Vector3d rover_xyz(x[0],x[1],x[2]);
            EstDopVel(rover_xyz);
        }

        return valid_flag;
    }

    void cSppSolver::InitSolver(tPPPLibConf C) {
        spp_conf_=C;
        para_=cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
        out_=new cOutSol(spp_conf_);
        out_->InitOutSol(spp_conf_,spp_conf_.fileC.sol);
        out_->WriteHead();
        out_->ref_sols_=ref_sols_;
    }

    bool cSppSolver::SolverProcess(tPPPLibConf C,int idx) {
        char buff[MAX_BUFF]={'\0'};
        if(idx!=-1) InitSolver(C);

        int i=0,num_epochs=rover_obs_.epoch_num;

        if(idx!=-1){
            num_epochs=idx+1;
        }
        for(i=idx==-1?0:idx;i<num_epochs;i++){

            epoch_idx_+=1;
            epoch_sat_obs_=rover_obs_.GetGnssObs().at(i);
            LOG(DEBUG)<<"START SOLVING: "<<i+1<<"th EPOCH, ROVER SATELLITE NUMBER: "<<epoch_sat_obs_.sat_num;

            UpdateGnssObs(spp_conf_,epoch_sat_obs_,REC_ROVER);
            if(SolverEpoch()){
                epoch_ok_++;
                SolutionUpdate();
                sprintf(buff,"%s SPP SOLVE SUCCESS POS: %14.3f %14.3f %14.3f VEL: %6.3f %6.3f %6.3f PDOP: %3.1f TOTAL SAT: %3d USED SAT: %3d",
                        ppplib_sol_.t_tag.GetTimeStr(1).c_str(),ppplib_sol_.pos[0],ppplib_sol_.pos[1],ppplib_sol_.pos[2],ppplib_sol_.vel[0],ppplib_sol_.vel[1],ppplib_sol_.vel[2],
                        ppplib_sol_.dops[1],epoch_sat_obs_.sat_num,num_valid_sat_);
                LOG(DEBUG)<<buff;

                out_->WriteSol(ppplib_sol_,epoch_idx_);
                epoch_sat_info_collect_.clear();
            }
            else{
                epoch_fail_++;
                epoch_sat_info_collect_.clear();
            }

        }
        LOG(INFO)<<"TOTAL EPOCH: "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
    }

    bool cSppSolver::SolverEpoch() {
        gnss_err_corr_.eph_model_.EphCorr(epoch_sat_info_collect_);
        Vector3d rr(full_x_[0],full_x_[1],full_x_[2]);
        CorrGnssObs(spp_conf_,rr);

        if(Estimator(spp_conf_)){
            ppplib_sol_.stat=SOL_SPP;
            ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            return true;
        }
        else{
            ppplib_sol_.stat=SOL_NONE;
            if(epoch_sat_info_collect_.size()) ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            return false;
        }
    }

    bool cSppSolver::SolutionUpdate() {
        int m=0;
        int num_pos=para_.NumPos();
        for(int i=0;i<num_pos;i++) ppplib_sol_.pos[i]=full_x_[i];

        int ic=para_.IndexClk(SYS_INDEX_GPS);
        if(sys_mask_[0]==0){
            for(m=1;m<NSYS;m++){
                if(sys_mask_[m]!=0) break;
            }
            for(int i=1;i<NSYS;i++){
                if(i==m) ppplib_sol_.clk_error[i]=full_x_[ic+i];
                else{
                    if(sys_mask_[i]==0){
                        ppplib_sol_.clk_error[i]=0.0;
                        continue;
                    }
                    ppplib_sol_.clk_error[i]=full_x_[ic+i]-full_x_[ic+m]; //ISB
                }
            }
        }
        else{
            for(int i=0;i<NSYS;i++){
                if(sys_mask_[i]){
                    ppplib_sol_.clk_error[i]=full_x_[ic+i];
                }
                else ppplib_sol_.clk_error[i]=0.0;
            }
        }
    }

    cPppSolver::cPppSolver() {}

    cPppSolver::cPppSolver(tPPPLibConf C) {
        ppp_conf_=C;
        ppp_conf_.mode=MODE_PPP;
        ppp_conf_.mode_opt=MODE_OPT_KINEMATIC;
        para_=cParSetting(ppp_conf_);
        num_full_x_=para_.GetPPPLibPar(ppp_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);
    }

    cPppSolver::~cPppSolver() {}

    void cPppSolver::InitSolver(tPPPLibConf C) {
        ppp_conf_=C;
        cReadGnssPreEph clk_reader(C.fileC.clk,nav_);
        clk_reader.Reading(1);

        cReadGnssPreEph orb_reader(C.fileC.sp3[1], nav_);
        for(int i=0;i<3;i++){
            if(C.fileC.sp3[i].empty()) continue;
            orb_reader.file_=C.fileC.sp3[i];
            orb_reader.Reading(0);
        }

        if(C.gnssC.tid_opt){
            cReadGnssErp erp_reader(C.fileC.erp,nav_);
            erp_reader.Reading();

            cReadGnssOcean blq_reader(C.fileC.blq, nav_,C.site_name,REC_ROVER);
            blq_reader.Reading();
        }

        cReadGnssAntex atx_reader(C.fileC.atx,nav_);
        atx_reader.Reading();
        atx_reader.AlignAntPar2Sat(C,*rover_obs_.GetStartTime(),nav_.sta_paras,nav_.sat_ant,nav_.rec_ant);

        para_=cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);

        out_=new cOutSol(ppp_conf_);
        out_->InitOutSol(ppp_conf_,ppp_conf_.fileC.sol);
        out_->WriteHead();
        out_->ref_sols_=ref_sols_;

        tPPPLibConf spp_conf=C;
        spp_conf.mode=MODE_SPP;
        spp_conf.mode_opt=MODE_OPT_KINEMATIC;
        spp_conf.gnssC.ion_opt=ION_KLB;
        spp_conf.gnssC.trp_opt=TRP_SAAS;
        spp_conf.gnssC.frq_opt=FRQ_SINGLE;
        spp_conf.gnssC.eph_opt=EPH_BRD;
        spp_solver_=new cSppSolver(spp_conf);
        spp_solver_->spp_conf_=spp_conf;
        spp_solver_->nav_=nav_;
        spp_solver_->ref_sols_=ref_sols_;
        spp_solver_->InitSolver(spp_conf);
        spp_solver_->spp_conf_.gnssC.res_qc=false;
    }

    bool cPppSolver::SolverProcess(tPPPLibConf C,int idx) {
        double rate=0.0;
        char buff[MAX_BUFF]={'\0'};
        if(idx!=-1) InitSolver(C);

        int i=0,num_epochs=rover_obs_.epoch_num;
        if(idx!=-1){
            num_epochs=idx+1;
        }

        for(i=idx==-1?0:idx;i<num_epochs;i++){

            epoch_sat_obs_=rover_obs_.GetGnssObs().at(i);
            LOG(DEBUG)<<"START PPP SOLVING "<<i+1<<"th EPOCH, ROVER SATELLITE NUMBER "<<epoch_sat_obs_.sat_num;

            epoch_idx_+=1;
            UpdateGnssObs(C,epoch_sat_obs_,REC_ROVER);
            InitEpochSatInfo(epoch_sat_info_collect_);

            if(SolverEpoch()){
                SolutionUpdate();
                sprintf(buff,"%s PPP SOLVE SUCCESS POS: %14.3f %14.3f %14.3f VEL: %6.3f %6.3f %6.3f PDOP: %3.1f TOTAL SAT: %3d USED SAT: %3d",
                        ppplib_sol_.t_tag.GetTimeStr(1).c_str(),ppplib_sol_.pos[0],ppplib_sol_.pos[1],ppplib_sol_.pos[2],ppplib_sol_.vel[0],ppplib_sol_.vel[1],ppplib_sol_.vel[2],
                        ppplib_sol_.dops[1],epoch_sat_obs_.sat_num,num_valid_sat_);
                LOG(DEBUG)<<buff;
            }

            epoch_sat_info_collect_.clear();
        }

        LOG(INFO)<<" TOTAL EPOCH: "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
    }

    bool cPppSolver::SolverEpoch() {
        InitSppSolver();
        if(spp_solver_->SolverEpoch()){

            Spp2Ppp();
            LOG(DEBUG)<<"PPP-SPP SOLVE SUCCESS ROVER, SPP POS "<<spp_solver_->ppplib_sol_.pos.transpose()<<" DOPPLER  VEL "<<spp_solver_->ppplib_sol_.vel.transpose();

            if(ppp_conf_.gnssC.eph_opt==EPH_PRE){
                gnss_err_corr_.eph_model_.EphCorr(epoch_sat_info_collect_);
            }

            CorrGnssObs(ppp_conf_,spp_solver_->ppplib_sol_.pos);
            PppCycSlip(ppp_conf_);
            StateTimeUpdate(ppp_conf_);

            if(Estimator(ppp_conf_)){
                epoch_ok_++;
                LOG(DEBUG)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPP SOLVE SUCCESS ";
                return true;
            }
            else{
                epoch_fail_++;
                ppplib_sol_=spp_solver_->ppplib_sol_;
                LOG(WARNING)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPP SOLVE FAILED ";
                return false;
            }
        }
        else{
            epoch_fail_++;
            ppplib_sol_.stat=SOL_NONE;
            if(epoch_sat_info_collect_.size()) ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            LOG(DEBUG)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPP-SPP SOLVE FAILED ";
            spp_solver_->epoch_sat_info_collect_.clear();
            return false;
        }
    }

    bool cPppSolver::Estimator(tPPPLibConf C) {

        // residuals
        VectorXd x;
        MatrixXd Px;
        int iter;
        for(iter=0;iter<max_iter_;iter++){
            x=full_x_;
            Px=full_Px_;

            if(!GnssObsRes(0,C,x.data())){
                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" MAKE PPP PRIOR RESIDUAL ERROR";
                return false;
            }

            kf_.Adjustment(omc_L_,H_,R_,x,Px,num_L_,num_full_x_);

            if(GnssObsRes(iter+1,C,x.data())) continue;
            else{
                full_x_=x;
                full_Px_=Px;
                ppplib_sol_.stat=SOL_PPP;
                ppplib_sol_.sigma=kf_.unit_weight_STD_;
                UpdateSatInfo(epoch_sat_info_collect_);
                break;
            }

            if(iter>max_iter_){
                LOG(DEBUG)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" PPP ITERATION OVERRUN";
                return false;
            }
        }

        if(ppplib_sol_.stat==SOL_PPP&&C.gnssC.ar_mode>=AR_PPP_AR){

        }

        int sat_no;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            for(int f=0;f<MAX_GNSS_USED_FRQ_NUM;f++){
                if(!epoch_sat_info_collect_[i].vsat[f]) continue;
                sat_no=epoch_sat_info_collect_[i].sat.sat_.no;
                previous_sat_info_[sat_no-1].outc[f]=0;
                previous_sat_info_[sat_no-1].lock[f]++;
//                if(previous_sat_info_[sat_no-1].lock[f]<0||previous_sat_info_[sat_no-1].fix[f]>=2){
//                    previous_sat_info_[sat_no-1].lock[f]++;
//                }
            }
        }

        vflag_.clear();
        return ppplib_sol_.stat;
    }

    bool cPppSolver::SolutionUpdate() {
        ppplib_sol_.stat=SOL_PPP;
        ppplib_sol_.valid_sat_num=num_valid_sat_;
        ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
        int num_pos=para_.NumPos();
        for(int i=0;i<num_pos;i++) ppplib_sol_.pos[i]=full_x_[i];
        int num_clk=para_.NumClk();
        for(int i=para_.IndexClk(SYS_INDEX_GPS);i<para_.IndexClk(SYS_INDEX_GPS)+num_clk;i++) ppplib_sol_.clk_error[i-para_.IndexClk(SYS_INDEX_GPS)]=full_x_[i];
        int ip=para_.IndexPos();
        for(int i=0;i<3;i++) ppplib_sol_.q_pos[ip+i]=full_Px_(ip+i,ip+i);

        out_->WriteSol(ppplib_sol_,epoch_idx_);
//        sol_collect_.push_back(ppplib_sol_);
    }

    void cPppSolver::AmbUpdate(tPPPLibConf C,double tt) {
        if(para_.NumAmb()<=0) return;

        int ia;
        double dt=0.0;
        double bias[64]={0};
        int slip[64]={0};
        int idx_flag[64]={0};
        tSatInfoUnit* sat_info= nullptr;
        int num_use_frqs=para_.GetGnssUsedFrqs();
        if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF) num_use_frqs=1;

        for(int f=0;f<num_use_frqs;f++){
            for(int i=0;i<MAX_SAT_NUM;i++){
                if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF){
                    if(previous_sat_info_[i].outc[0]>5){
                        ia=para_.IndexAmb(0,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                        ia=para_.IndexAmb(1,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                        ia=para_.IndexAmb(2,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                    }
                }
                else{
                    if(previous_sat_info_[i].outc[f]>5){
                        ia=para_.IndexAmb(f,i+1);
                        InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                    }
                };

            }

            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                sat_info=&epoch_sat_info_collect_.at(i);
                ia=para_.IndexAmb(f,sat_info->sat.sat_.no);
                bias[i]=0.0;
                idx_flag[i]=0.0;

                if(sat_info->stat!=SAT_USED) continue;

                // ionosphere-free combination
                if(C.gnssC.ion_opt==ION_IF){
                    if(C.gnssC.frq_opt==FRQ_SINGLE){
                        // single frequency
                        if(sat_info->raw_P[f]!=0.0||sat_info->raw_L[f]!=0.0){
                            bias[i]=sat_info->raw_L[f]*sat_info->lam[f]-sat_info->raw_P[f];
                            slip[i]=sat_info->slip[f];
                            idx_flag[i]=0;
                        }
                    }
                    else if(C.gnssC.frq_opt==FRQ_DUAL){
                        // dual frequency
                        if(sat_info->cor_if_L[f]==0.0||sat_info->cor_if_P[f]==0.0) continue;
                        bias[i]=sat_info->cor_if_L[f]-sat_info->cor_if_P[f];
                        slip[i]=sat_info->slip[f];
                        idx_flag[i]=0;
                    }
                    else if(C.gnssC.frq_opt==FRQ_TRIPLE){
                        // triple frequency with one IF-combination, while missing triple observables,use dual-frequency(L1_L2/L1_L5) instead
                        if(sat_info->cor_if_P[0]!=0.0&&sat_info->cor_if_L[0]!=0.0){
                            // L1_L2_L5
                            bias[i]=sat_info->cor_if_L[0]-sat_info->cor_if_P[0];
                            slip[i]=sat_info->slip[0];
                            idx_flag[i]=0;
                            sat_info->tf_if_idx[0]=1;
                        }
                        else if(sat_info->cor_if_P[1]!=0.0&&sat_info->cor_if_L[1]!=0.0){
                            // L1_L2
                            bias[i]=sat_info->cor_if_L[1]-sat_info->cor_if_P[1];
                            slip[i]=sat_info->slip[0];
                            idx_flag[i]=1;
                            sat_info->tf_if_idx[1]=1;
                        }
                        else if(sat_info->cor_if_P[2]!=0.0&&sat_info->cor_if_L[2]!=0.0){
                            // L1_L5
                            bias[i]=sat_info->cor_if_L[2]-sat_info->cor_if_P[2];
                            slip[i]=sat_info->slip[0];
                            idx_flag[i]=2;
                            sat_info->tf_if_idx[2]=1;
                        }
                        else{
                            sat_info->stat=SAT_NO_USE;
                            continue;
                        }
                    }
                }
                else if(C.gnssC.ion_opt==ION_IF_DUAL&&C.gnssC.frq_opt==FRQ_TRIPLE){
                    // L1_L2
                    if(sat_info->cor_if_P[f]!=0.0&&sat_info->cor_if_L[f]!=0.0){
                        bias[i]=sat_info->cor_if_L[f]-sat_info->cor_if_P[f];
                        slip[i]=sat_info->slip[f];
                    }
                    else{
//                        sat_info->stat=SAT_NO_USE;
                        continue;
                    }
                }
                else if(sat_info->raw_L[f]!=0.0&&sat_info->raw_P[f]!=0.0){
                    slip[i]=sat_info->slip[f];
                    bias[i]=sat_info->raw_L[f]*sat_info->lam[f]-sat_info->raw_P[f];
                }

                if(full_x_[ia]==0.0||slip[i]||bias[i]==0.0) continue;
            }

            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                sat_info=&epoch_sat_info_collect_.at(i);
                dt=spp_solver_->ppplib_sol_.t_tag.TimeDiff(previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);

                // for UC: f==0 L1, f==1 L2, f==2 L3
                // for DF-IF: f==0 L1_L2/L1_L5
                // for TF-IF: f==0 L1_L2_L5, f==1 L1_L2, f==2 L1_L5
                // for TF-IF-DUAL: f==0 L1_L2, f==1 L1_L5
                if(idx_flag[i]==0){
                    ia=para_.IndexAmb(f,sat_info->sat.sat_.no);
                }
                else{
                    ia=para_.IndexAmb(idx_flag[i],sat_info->sat.sat_.no);
                }

                                // UC-L5                                                                    // TF-IF       L1_L2_L5       L1_L5
                if((C.gnssC.frq_opt==FRQ_TRIPLE&&f==FRQ_TRIPLE)||(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF&&(idx_flag[i]==0||idx_flag[i]==2))){
                    //  L1_L2_L5/L1_L5/L5
                    full_Px_(ia,ia)+=3E-5*fabs(dt);
                }
                else {
                    // L1/L2/L1_L2 psd=0.0
                    full_Px_(ia,ia)+=C.gnssC.ait_psd[0]*fabs(dt);
                }

                if(bias[i]==0.0||(full_x_[ia]!=0.0&&!slip[i])) continue;

                InitX(bias[i],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                string s="UC-L";
                if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF){
                    if(sat_info->tf_if_idx[0]==1){
                        s="TF-L_123";
                    }
                    else if(sat_info->tf_if_idx[1]==1){
                        s="TF-L_12";
                    }
                    else if(sat_info->tf_if_idx[2]==1){
                        s="TF-L_15";
                    }
                }
                else if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF_DUAL){
                    if(f==0){
                        s="TF-IF_DUAL-L12";
                    }
                    else if(f==1){
                        s="TF-IF_DUAL-L13";
                    }
                }
                else if(C.gnssC.frq_opt==FRQ_DUAL&&C.gnssC.ion_opt==ION_IF){
                    s="DF-L_12";
                }
                else if(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF){
                    s="SF-L";
                }
                else{
                    s+=to_string(f+1);
                }
                LOG(DEBUG)<<epoch_sat_info_collect_.at(i).t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_.at(i).sat.sat_.id<<" "<<s<<" AMBIGUITY INITIALIZED "<<bias[i];
            }
        }
    }

    void cPppSolver::IonUpdate(tPPPLibConf C,double tt) {
        if(para_.NumIon()<=0) return;

        int ii;
        tSatInfoUnit* sat_info= nullptr;
        double ion;
        for(int i=0;i<MAX_SAT_NUM;i++){
            ii=para_.IndexIon(i+1);
            if(full_x_[ii]!=0.0&&previous_sat_info_[i].outc[0]>120){
                full_x_[ii]=0.0;
            }
        }

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            if(sat_info->stat!=SAT_USED) continue;
            ii=para_.IndexIon(sat_info->sat.sat_.no);
            if(full_x_[ii]==0.0){
                ii=para_.IndexIon(sat_info->sat.sat_.no);
                if(sat_info->cor_P[0]==0.0||sat_info->cor_P[1]==0.0){
                    sat_info->stat=SAT_NO_USE;
                    continue;
                }
                ion=(sat_info->cor_P[0]-sat_info->cor_P[1])/(1.0-SQR(sat_info->lam[1]-sat_info->lam[0]));
                InitX(ion,SQR(60.0),ii,full_x_.data(),full_Px_.data());
            }
            else{
                tt=spp_solver_->ppplib_sol_.t_tag.TimeDiff(previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);
                double sinel=sin(MAX(sat_info->el_az[0],10.0*D2R));
                full_Px_(ii,ii)+=C.gnssC.ait_psd[1]/sinel*fabs(tt);
            }
        }
    }

    void cPppSolver::TrpUpdate(tPPPLibConf C,double tt) {
        if(para_.NumTrp()<=0) return;

        int it=para_.IndexTrp();
        if(full_x_[it]==0.0){
            full_x_[it]=0.175;
            full_Px_(it,it)=SQR(0.3);
            if(C.gnssC.trp_opt==TRP_EST_GRAD){
                full_x_[it+1]=full_x_[it+2]=1E-6;
                full_Px_(it+1,it+1)=full_Px_(it+2,it+2)=SQR(0.01);
            }
            return;
        }
        else{
            full_Px_(it,it)+=C.gnssC.ait_psd[2]*tt;
            if(C.gnssC.trp_opt==TRP_EST_GRAD){
                full_Px_(it+1,it+1)=full_Px_(it+1,it+1)+=C.gnssC.ait_psd[2]*0.1*tt;
            }
        }
    }

    void cPppSolver::GloIfcbUpdate(tPPPLibConf C,double tt) {
        if(para_.NumGloIfcb()<=0) return;

        int ii,j;
        if(C.gnssC.glo_ifcb_opt==GLO_IFCB_LNF){
            ii=para_.IndexGloIfcb(1);
            if(full_x_[ii]==0.0) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
        }
        else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_QUAD){
            ii=para_.IndexGloIfcb(1);
            if(full_x_[ii]==0.0) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
            ii=para_.IndexGloIfcb(2);
            if(full_x_[ii]==0.0) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
        }
        else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT){
            double dtr=spp_solver_->ppplib_sol_.clk_error[SYS_INDEX_GLO];
            int ic=para_.IndexClk(SYS_INDEX_GLO);
            InitX(0.0,0.0,ic,full_x_.data(),full_Px_.data());
            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                if(epoch_sat_info_collect_[i].stat!=SAT_USED) continue;
                ii=para_.IndexGloIfcb(epoch_sat_info_collect_[i].sat.sat_.prn);
                if(full_x_[ii]==0.0){
                    InitX(dtr,SQR(60.0),ii,full_x_.data(),full_Px_.data());
                }
                else full_x_(ii,ii)+=SQR(0.001)*tt;
            }
        }
        else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1FRQ){
            for(int i=0;i<epoch_sat_info_collect_.size();i++){
                if(epoch_sat_info_collect_[i].stat!=SAT_USED) continue;
                for(j=1;j<14;j++){
                    if(GLO_FRQ_NUM[j-1]==nav_.brd_glo_eph[epoch_sat_info_collect_[i].brd_eph_index].frq) break;
                }
                ii=para_.IndexGloIfcb(j);
                if(full_x_[ii]==0.0) InitX(0.1,SQR(60.0),ii,full_x_.data(),full_Px_.data());
            }
        }
    }

    void cPppSolver::ClkUpdate(tPPPLibConf C,double tt) {
        if(para_.NumClk()<=0) return;
        int ic=para_.IndexClk(SYS_INDEX_GPS);
        bool init=false;
        double psd=0.01;

        if(SQR(full_x_[ic])+SQR(full_x_[ic+1])+SQR(full_x_[ic+2])+SQR(full_x_[ic+3])+SQR(full_x_[ic+4])<=0.0){
            init=true;
        }

        if(init){
            for(int i=0;i<NSYS;i++){
                // GPS_clock, GC_ISB, GE_ISB, GR_ISB, GJ_ISB
                if(spp_solver_->ppplib_sol_.clk_error[i]==0.0) continue;
                InitX(spp_solver_->ppplib_sol_.clk_error[i],SQR(60.0),ic+i,full_x_.data(),full_Px_.data());
            }
        }
        else{
            if(spp_solver_->ppplib_sol_.clk_error[SYS_INDEX_GPS]!=0.0){
                // GPS clock error reinitialize
                InitX(spp_solver_->ppplib_sol_.clk_error[SYS_INDEX_GPS],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // BDS,GAL,GLO,QZS ISB modeling
                for(int i=SYS_INDEX_BDS;i<NSYS;i++){
                    if(spp_solver_->ppplib_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    if(full_x_[ic+i]==0.0){
                        InitX(spp_solver_->ppplib_sol_.clk_error[i-SYS_INDEX_GPS],SQR(60.0),ic+i,full_x_.data(),full_Px_.data());
                    }
                    else full_Px_(ic+i,ic+i)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(ppplib_sol_.clk_error[SYS_INDEX_BDS]!=0.0){
                // BDS clock error reinitialize
                ic+=SYS_INDEX_BDS;
                InitX(spp_solver_->ppplib_sol_.clk_error[SYS_INDEX_BDS],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // GAL,GLO,QZS ISB modeling
                for(int i=SYS_INDEX_GAL;i<NSYS;i++){
                    if(spp_solver_->ppplib_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    if(full_x_[ic+i-SYS_INDEX_BDS]==0.0){
                        InitX(spp_solver_->ppplib_sol_.clk_error[i-SYS_INDEX_GPS],SQR(60.0),ic+i-SYS_INDEX_BDS,full_x_.data(),full_Px_.data());
                    }
                    else full_Px_(ic+i,ic+i)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(ppplib_sol_.clk_error[SYS_INDEX_GAL]!=0.0){
                // GAL clock error reinitialize
                ic+=SYS_INDEX_GAL;
                InitX(spp_solver_->ppplib_sol_.clk_error[SYS_INDEX_GAL],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // GLO,QZS ISB modeling
                for(int i=SYS_INDEX_GLO;i<NSYS;i++){
                    if(spp_solver_->ppplib_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    full_Px_(ic+i,ic+i)+=SQR(psd)*tt;
                }
                return ;
            }
            else if(ppplib_sol_.clk_error[SYS_INDEX_GLO]!=0.0){
                // GLO clock error reinitialize
                ic+=SYS_INDEX_GLO;
                InitX(spp_solver_->ppplib_sol_.clk_error[SYS_INDEX_GLO],SQR(60.0),ic,full_x_.data(),full_Px_.data());

                // QZS ISB modeling
                for(int i=SYS_INDEX_QZS;i<NSYS;i++){
                    if(spp_solver_->ppplib_sol_.clk_error[i-SYS_INDEX_GPS]==0.0) continue;
                    full_Px_(ic+i,ic+i)+=SQR(psd)*tt;
                }
                return ;
            }
        }
    }

    void cPppSolver::IfbUpdate(tPPPLibConf C, double tt) {
        if(para_.NumIfb()<=0) return;

        int i,ifcb=para_.IndexIfb(SYS_INDEX_GPS);

        for(i=0;i<NSYS;i++){
            ifcb+=i;
            if(sys_mask_[i]){
                if(full_x_[ifcb]==0.0){
                    InitX(0.01,SQR(60),ifcb,full_x_.data(),full_Px_.data());
                }
                else{
//                    full_Px_(ifcb,ifcb)+=10E-6*tt;
                }
            }
        }
    }

    void cPppSolver::PosUpdate(tPPPLibConf C) {
        int ip=para_.IndexPos();
        Eigen::Vector3d var(SQR(60.0),SQR(60.0),SQR(60.0));

        if((SQR(full_x_[0])+SQR(full_x_[1])+SQR(full_x_[2]))==0.0){
            for(int i=0;i<3;i++) full_x_[i]=spp_solver_->ppplib_sol_.pos[i];
            full_Px_.block<3,3>(0,0)=var.asDiagonal();
            return;
        }

        if(C.mode_opt==MODE_OPT_STATIC){
            // nothing
        }
        else if(C.mode_opt==MODE_OPT_KINE_SIM||C.mode_opt==MODE_OPT_KINEMATIC||C.mode_opt==MODE_OPT_PPP){
            for(int i=0;i<3;i++) full_x_[i]=spp_solver_->ppplib_sol_.pos[i];
            full_Px_.block<3,3>(0,0)=var.asDiagonal();
        }
    }


    void cPppSolver::StateTimeUpdate(tPPPLibConf C) {
        double tt=spp_solver_->ppplib_sol_.t_tag.TimeDiff(ppplib_sol_.t_tag.t_);

        //position
        PosUpdate(C);

        //clock
        ClkUpdate(C,tt);

        //dcb

        //ifb
        IfbUpdate(C,tt);

        //glo ifcb
        GloIfcbUpdate(C,tt);

        //trp
        TrpUpdate(C,tt);

        //ion
        IonUpdate(C,tt);

        //amb
        AmbUpdate(C,tt);

    }

    void cPppSolver::InitSppSolver() {
        spp_solver_->epoch_idx_=epoch_idx_;
        spp_solver_->epoch_sat_info_collect_=epoch_sat_info_collect_;
    }

    void cPppSolver::Spp2Ppp() {
        spp_solver_->SolutionUpdate();
        epoch_sat_info_collect_.clear();
        epoch_sat_info_collect_=spp_solver_->epoch_sat_info_collect_;
        spp_solver_->epoch_sat_info_collect_.clear();
        for(int i=0;i<NSYS;i++) sys_mask_[i]=spp_solver_->sys_mask_[i];
    }

    void cPppSolver::PppCycSlip(tPPPLibConf C) {
        tSatInfoUnit* sat_info= nullptr;
        cTime t=ppplib_sol_.t_tag;
        double dt=C.gnssC.sample_rate;

        if(t.t_.long_time!=0.0) dt=epoch_sat_info_collect_[0].t_tag.TimeDiff(t.t_);

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);

            gnss_obs_operator_.LliCycleSlip(C,*sat_info,2,dt,REC_ROVER);
            if(sat_info->sat.sat_.sys!=SYS_GLO) gnss_obs_operator_.MwCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, nullptr,previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);
            gnss_obs_operator_.GfCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, nullptr);
            if(sat_info->sat.sat_.sys!=SYS_GLO) gnss_obs_operator_.SmoothMw(C,sat_info, nullptr);
        }
    }

    int cPppSolver::GnssObsRes(int post, tPPPLibConf C, double *x) {
        num_L_=0;
        num_valid_sat_=0;
        if(vflag_.size()) vflag_.clear();
        for(int i=0;i<NSYS;i++) sys_mask_[i]=0;

        bool have_larger_res=false;
        char buff[MAX_BUFF]={'\0'};
        vector<double>H,omcs,meas_var_vec,frqs,idxs,types,larger_omcs;
        double omc,r,meas,meas_var,ion=0.0,amb=0.0,fact;
        tSatInfoUnit* sat_info= nullptr;
        Vector3d rover_xyz(x[0],x[1],x[2]);
        cTime t=epoch_sat_info_collect_[0].t_tag;
        Vector3d dr;
        LOG(DEBUG)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<(post==0?" PRIOR(":(" POST("))<<post<<")"<<" RESIDUAL, PRIOR ROVER COORDINATE "<<rover_xyz.transpose();

        if(C.gnssC.tid_opt>TID_OFF){
            gnss_err_corr_.tid_model_.TidCorr(t.Gpst2Utc(),rover_xyz,dr);
            rover_xyz+=dr;
        }

        Vector3d rover_blh=Xyz2Blh(rover_xyz),sig_vec;

        int num_used_frq=para_.GetGnssUsedFrqs();
        int num_used_obs_type=para_.GetNumObsType();
        if(num_used_frq<=0) return 0;
        if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF) num_used_frq=1;

        double sat_clk,cbias,trp_del,alpha,beta;
        int m,idx_isb,idx_gps_clk=para_.IndexClk(SYS_INDEX_GPS),idx_ifcb,ii;
        double base_rec_clk=0.0,isb=0.0,rec_clk=0.0,glo_ifcb=0.0;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){

            sat_info=&epoch_sat_info_collect_.at(i);

            if(sat_info->stat!=SAT_USED) continue;

            if((r=GeoDist(sat_info->pre_pos,rover_xyz,sig_vec))<=0.0){
                sat_info->stat=SAT_NO_USE;
                continue;
            }

            if(SatElAz(rover_blh,sig_vec,sat_info->el_az)<C.gnssC.ele_min*D2R){
                sat_info->stat=SAT_LOW_EL;
                continue;
            }

            sat_info->sagnac=gnss_err_corr_.SagnacCorr(sat_info->pre_pos,rover_xyz);
            sat_info->shapiro=gnss_err_corr_.ShapiroCorr(sat_info->sat.sat_.sys,sat_info->pre_pos,rover_xyz);
            gnss_err_corr_.trp_model_.InitSatInfo(sat_info,&rover_blh);
            gnss_err_corr_.trp_model_.GetTrpError(0.0,x,para_.IndexTrp());
            gnss_err_corr_.trp_model_.UpdateSatInfo();
            gnss_err_corr_.ion_model_.InitSatInfo(sat_info,&rover_blh);
            gnss_err_corr_.ion_model_.GetIonError();
            gnss_err_corr_.ion_model_.UpdateSatInfo();

            int obs_type,frq;
            for(int f=0;f<num_used_frq*num_used_obs_type;f++){
                frq=f/num_used_obs_type;
                obs_type=f%num_used_obs_type;

                //ionosphere-free
                if(C.gnssC.ion_opt==ION_IF){
                    if(C.gnssC.frq_opt==FRQ_TRIPLE){
                        if(sat_info->tf_if_idx[0]==1){
                            // L1_L2_L5
                            meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[0]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[0]:sat_info->cor_D[0]);
                        }
                        else if(sat_info->tf_if_idx[1]==1){
                            // L1_L2
                            meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[1]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[1]:sat_info->cor_D[1]);
                        }
                        else if(sat_info->tf_if_idx[2]==2){
                            // L1_L5
                            meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[2]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[2]:sat_info->cor_D[2]);
                        }
                    }
                    else meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[frq]:sat_info->cor_D[frq]);
                }
                else if(C.gnssC.ion_opt==ION_IF_DUAL&&C.gnssC.frq_opt==FRQ_TRIPLE){
                    meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[frq]:sat_info->cor_D[frq]);
                }
                else{
                    meas=obs_type==GNSS_OBS_CODE?sat_info->cor_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_L[frq]:sat_info->cor_D[frq]);
                }

                if(meas==0.0){
#if 0
                    sat_info->stat=GNSS_OBS_CODE?sat_info->stat=SAT_NO_CODE:sat_info->stat=SAT_NO_PHASE;
                    LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" MISSING  OBS"<<(obs_type==GNSS_OBS_CODE?"P":"L")<<frq+1;
#endif
                    continue;
                }
                ion=0.0;
                amb=0.0;

                //clock and inter-system bias
                if(!post) for(int j=0;j<num_full_x_;j++) H.push_back(j<3?-sig_vec[j]:0.0);
                for(m=0;m<NSYS;m++){
                    if(x[m+idx_gps_clk]!=0.0) break;
                }
                if(sat_info->sat.sat_.sys_idx==m){
                    if(sat_info->sat.sat_.sys==SYS_GLO&&C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) break;
                    if(!post) H[m+idx_gps_clk+num_L_*num_full_x_]=1.0;
                    base_rec_clk=x[m+idx_gps_clk];
                }
                else{
                    for(idx_isb=0;idx_isb<NSYS;idx_isb++){
                        if(idx_isb==m){
                            if(!post) H[m+idx_gps_clk+num_L_*num_full_x_]=1.0;
                        }
                        else{
                            if(idx_isb!=sat_info->sat.sat_.sys_idx) continue;
                            if(x[idx_gps_clk+idx_isb]==0.0) continue;
                            if(sat_info->sat.sat_.sys==SYS_GLO&&C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) continue;
                            if(!post) H[idx_isb+idx_gps_clk+num_L_*num_full_x_]=1.0;
                            isb=x[idx_isb+idx_gps_clk];
                        }
                    }
                }

                rec_clk=base_rec_clk+isb;
#if 1
                // inter-frequency clock bias
                     // TF-UC-L5                                                                       //  TF-IF-DUAL L1_L5
                if((C.gnssC.frq_opt>=FRQ_TRIPLE&&frq==FRQ_TRIPLE&&obs_type==GNSS_OBS_CODE)||(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF_DUAL&&frq==1&&obs_type==GNSS_OBS_CODE)||(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF&&(sat_info->tf_if_idx[1]||sat_info->tf_if_idx[2])&&obs_type==GNSS_OBS_CODE)){
                    int ifcb=para_.IndexIfb(sat_info->sat.sat_.sys_idx);
                    rec_clk+=x[ifcb];
                    if(!post){
                        H[ifcb+num_L_*num_full_x_]=1.0;
                    }
                }
#endif
                // glonass inter-frequency code bias
                if(sat_info->sat.sat_.sys==SYS_GLO&&obs_type==GNSS_OBS_CODE){
                    int jj;
                    double frq=nav_.brd_glo_eph[sat_info->brd_eph_index].frq;
                    if(C.gnssC.glo_ifcb_opt==GLO_IFCB_LNF){
                        idx_ifcb=para_.IndexGloIfcb(1);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=frq;
                        glo_ifcb=frq*full_x_[idx_ifcb];
                        rec_clk+=glo_ifcb;
                    }
                    else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_QUAD){
                        idx_ifcb=para_.IndexGloIfcb(1);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=frq;
                        glo_ifcb=frq*full_x_[idx_ifcb];
                        idx_ifcb=para_.IndexGloIfcb(2);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=SQR(frq);
                        glo_ifcb+=SQR(frq)*full_x_[idx_ifcb];
                        rec_clk+=glo_ifcb;
                    }
                    else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT){
                        idx_ifcb=para_.IndexGloIfcb(sat_info->sat.sat_.prn);
                        rec_clk=full_x_[idx_ifcb];
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=1.0;
                    }
                    else if(C.gnssC.glo_ifcb_opt==GLO_IFCB_1FRQ){
                        for(jj=1;jj<14;jj++){
                            if(frq==GLO_FRQ_NUM[jj]) break;
                        }
                        idx_ifcb=para_.IndexGloIfcb(jj);
                        if(!post) H[idx_ifcb+num_L_*num_full_x_]=1.0;
                    }
                }

                // tropospheric delay
                int it=para_.IndexTrp();
                if(!post) H[it+num_L_*num_full_x_]=sat_info->trp_wet_delay[1];
                if(C.gnssC.trp_opt==TRP_EST_GRAD&&!post){
                    H[it+1+num_L_*num_full_x_]=sat_info->trp_wet_delay[2];
                    H[it+2+num_L_*num_full_x_]=sat_info->trp_wet_delay[3];
                }

                // ionospheric delay
                if(C.gnssC.ion_opt>ION_IF_DUAL){
                    ii=para_.IndexIon(sat_info->sat.sat_.no);
                    fact=SQR(sat_info->frq[0]/sat_info->frq[frq])*(obs_type==GNSS_OBS_PHASE?(-1.0):1.0);
                    ion=x[ii]*fact;

                    if(x[ii]==0.0){
                        LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<"NO VALUE FOR IONOSPHERIC PARAMETER";
                        continue;
                    }
                    if(!post) H[ii+num_L_*num_full_x_]=fact;
                }

                //Ambiguity
                if(obs_type==GNSS_OBS_PHASE||(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF&&obs_type==GNSS_OBS_CODE)){
                    int ia;
                    if(C.gnssC.frq_opt==FRQ_TRIPLE&&C.gnssC.ion_opt==ION_IF){
                        if(sat_info->tf_if_idx[0]==1){
                            ia=para_.IndexAmb(frq,sat_info->sat.sat_.no);
                        }
                        else if(sat_info->tf_if_idx[0]==0&&sat_info->tf_if_idx[1]!=0){
                            ia=para_.IndexAmb(1,sat_info->sat.sat_.no);
                        }
                        else if(sat_info->tf_if_idx[0]==0&&sat_info->tf_if_idx[2]!=0){
                            ia=para_.IndexAmb(2,sat_info->sat.sat_.no);
                        }
                    }
                    else{
                        ia=para_.IndexAmb(frq,sat_info->sat.sat_.no);
                    }

                    amb=x[ia];
                    if(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF&&obs_type==GNSS_OBS_CODE) amb*=0.5;
                    if(amb==0.0) {
                        LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<"NO VALUE FOR AMBIGUITY PARAMETER";
                        continue;
                    }
                    if(!post) H[ia+num_L_*num_full_x_]=1.0;
                }

                sat_clk=sat_info->pre_clk[0]*CLIGHT;
                trp_del=sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0];
                cbias=0.0;
                alpha=SQR(sat_info->frq[0])/(SQR(sat_info->frq[0])-SQR(sat_info->frq[1]));
                beta=-SQR(sat_info->frq[1])/(SQR(sat_info->frq[0])-SQR(sat_info->frq[1]));

                if(obs_type==GNSS_OBS_CODE){
                    if(C.gnssC.ion_opt==ION_IF){
                        cbias=alpha*sat_info->code_bias[0]+beta*sat_info->code_bias[1];
                    }
                    else cbias=sat_info->code_bias[f];
                }

                omc=meas-(r+sat_info->sagnac+rec_clk-sat_clk+trp_del+ion+amb-sat_info->shapiro);
                meas_var=GnssMeasVar(C,obs_type==0?GNSS_OBS_CODE:(obs_type==1?GNSS_OBS_PHASE:GNSS_OBS_DOPPLER),*sat_info)+sat_info->trp_var;
                if(obs_type==GNSS_OBS_CODE) meas_var*=previous_sat_info_[sat_info->sat.sat_.no-1].c_var_factor[frq];
                else if(obs_type==GNSS_OBS_PHASE) meas_var*=previous_sat_info_[sat_info->sat.sat_.no-1].p_var_factor[frq];

                if(obs_type==GNSS_OBS_CODE){
                    sprintf(buff,"%4d el=%3.1f var=%9.5f omc=%8.4f cor_obs=%14.3f range=%14.3f rec_clk=%12.3f sat_clk=%12.3f trp=%7.3f ion=%7.3f shaprio=%7.3f cbias=%7.3f",
                            epoch_idx_,sat_info->el_az[0]*R2D,meas_var,omc,meas,r+sat_info->sagnac,rec_clk,sat_clk,trp_del,ion,sat_info->shapiro,cbias);
                }
                else if(obs_type==GNSS_OBS_PHASE){
                    sprintf(buff,"%4d el=%3.1f var=%9.5f omc=%8.4f cor_obs=%14.3f range=%14.3f rec_clk=%12.3f sat_clk=%12.3f trp=%7.3f ion=%7.3f shaprio=%7.3f amb  =%7.3f phw=%7.3f",
                            epoch_idx_,sat_info->el_az[0]*R2D,meas_var,omc,meas,r+sat_info->sagnac,rec_clk,sat_clk,trp_del,ion,sat_info->shapiro,amb,sat_info->phase_wp);
                }

                LOG(DEBUG)<<sat_info->sat.sat_.id<<" "<<(C.gnssC.ion_opt==ION_IF?((obs_type==GNSS_OBS_CODE?"IF_P":"IF_L")):((obs_type==GNSS_OBS_CODE?"P":"L")))<<frq+1<<" "<<buff;

                if(!post&&omc>C.gnssC.max_prior){
                    if(obs_type==GNSS_OBS_CODE) sat_info->stat=SAT_PRI_RES_C;
                    else sat_info->stat=SAT_PRI_RES_P;
                    LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" DETECTING OUTLIER IN PRIOR RESIDUALS "<<(obs_type==GNSS_OBS_CODE?"P":"L")<<frq+1<<"= "<<omc;
                    continue;
                }

                if(obs_type==GNSS_OBS_PHASE&&!post&&sat_info->sat.sat_.sys==SYS_BDS){
//                    if(fabs(omc-omcs.back())>7.0){
//                        sat_info->stat=SAT_LAR_CP_DIFF;
//                        LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<"("<<epoch_idx_<<")"<<" "<<sat_info->sat.sat_.id<<" CODE AND PHASE PRIOR RESIDUAL SIGNIFICANT DIFFERENCE code_res="<<omc<<" phase_res="<<omcs.back();
//                        continue;
//                    }
                }

                omcs.push_back(omc);
                meas_var_vec.push_back(meas_var);

                if(post&&C.gnssC.res_qc&&fabs(omc)>sqrt(meas_var)*3.0){
//                    have_larger_res=true;
//                    larger_omcs.push_back(omc);
//                    idxs.push_back(i);
//                    frqs.push_back(frq);
//                    types.push_back(obs_type);
//                    LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<"("<<epoch_idx_<<") "<<sat_info->sat.sat_.id<<" "<<(obs_type==GNSS_OBS_CODE?"P":"L")
//                                <<frq+1<<" LARGER POST RESIDUAL res="<<omc<<" thres="<<sqrt(meas_var)*5.0;
                }

                if(!post){
                    if(obs_type==GNSS_OBS_CODE)  sat_info->prior_res_P[frq]=omc;
                    if(obs_type==GNSS_OBS_PHASE) sat_info->prior_res_L[frq]=omc;
                }
                else{
                    if(obs_type==GNSS_OBS_CODE)  sat_info->post_res_P[frq]=omc;
                    if(obs_type==GNSS_OBS_PHASE) sat_info->post_res_L[frq]=omc;
                }

                sat_info->vsat[frq]=1;
                previous_sat_info_[sat_info->sat.sat_.no-1].vsat[frq]=1;

                vflag_.push_back((i<<8)|(obs_type<<4)|(frq));

                num_L_++;
                if(f==0) num_valid_sat_++;
            }
        }

        if(post&&ppp_conf_.gnssC.res_qc&&!(C.gnssC.frq_opt==FRQ_SINGLE&&C.gnssC.ion_opt==ION_IF)){
            if(PppResidualQc(omcs,meas_var_vec)){
                have_larger_res=true;
            }
        }

        if(!post){
            omc_L_=Map<VectorXd>(omcs.data(),num_L_);
            H_=Map<MatrixXd>(H.data(),num_full_x_,num_L_);
            Map<VectorXd> var_vec(meas_var_vec.data(),num_L_);
            R_=var_vec.asDiagonal();
        }

        if(post&&C.gnssC.res_qc&&have_larger_res&&larger_omcs.size()>0){
            for(int i=0;i<larger_omcs.size();i++) larger_omcs[i]=fabs(larger_omcs[i]);
            auto max_omc=max_element(larger_omcs.begin(),larger_omcs.end());
            int idx_max_omc=distance(begin(larger_omcs),max_omc);
            epoch_sat_info_collect_[idxs[idx_max_omc]].stat=types[idx_max_omc]==GNSS_OBS_CODE?SAT_POS_RES_C:SAT_POS_RES_P;
            LOG(WARNING)<<epoch_sat_info_collect_[idxs[idx_max_omc]].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[idxs[idx_max_omc]].sat.sat_.id
                      <<" "<<(types[idx_max_omc]==GNSS_OBS_CODE?"P":"L")<<frqs[idx_max_omc]+1<<" EXCLUDED BY POST RESIDUAL res="<<larger_omcs[idx_max_omc];
        }

        H.clear();omcs.clear();meas_var_vec.clear();types.clear();frqs.clear();idxs.clear();larger_omcs.clear();

        return post?have_larger_res:num_valid_sat_;
    }

    bool cPppSolver::PppResidualQc(vector<double>omcs,vector<double>R) {

        bool flag=false;
        int i,obs_type,f,idx_sat,sat_no,ia;
        vector<double>code_v,phase_v,norm_code_v,norm_phase_v;
        vector<int>idx_code_v,idx_phase_v;
        double el;
        for(i=0;i<vflag_.size();i++){
           obs_type=vflag_[i]>>4&0xF;
           if(obs_type==GNSS_OBS_CODE){
               idx_code_v.push_back(i);
               code_v.push_back(fabs(omcs[i]));
               norm_code_v.push_back(fabs(omcs[i])/SQRT(R[i]));
           }
           else{
               idx_phase_v.push_back(i);
               phase_v.push_back(fabs(omcs[i]));
               norm_phase_v.push_back(fabs(omcs[i])/SQRT(R[i]));
           }
        }


        bool code_flag=false;
        auto max_code_v=max_element(begin(code_v),end(code_v));
        int idx_max_code_v=distance(begin(code_v),max_code_v);
        idx_sat=vflag_[idx_code_v[idx_max_code_v]]>>8&0xFF;
        el=epoch_sat_info_collect_[idx_sat].el_az[0];

        if(*max_code_v>3.0/sin(el)){
            code_flag=true;
        }

        bool norm_code_flag=false;
        auto max_norm_code_v=max_element(begin(norm_code_v),end(norm_code_v));
        int idx_max_norm_code_v=distance(begin(norm_code_v),max_norm_code_v);
        if(*max_norm_code_v>2.0){
            norm_code_flag=true;
        }

        if(code_flag){
            flag=true;
            idx_sat=vflag_[idx_code_v[idx_max_code_v]]>>8&0xFF;
            epoch_sat_info_collect_[idx_sat].stat=SAT_NO_USE;
            f=vflag_[idx_code_v[idx_max_code_v]]&0xF;
            sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
            ia=para_.IndexAmb(f,sat_no);
            LOG(WARNING)<<" ";
            return flag;
        }
        if(norm_code_flag){
            flag=true;
            idx_sat=vflag_[idx_code_v[idx_max_norm_code_v]]>>8&0xFF;
            epoch_sat_info_collect_[idx_sat].stat=SAT_NO_USE;
            f=vflag_[idx_code_v[idx_max_norm_code_v]]&0xF;
            sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
            ia=para_.IndexAmb(f,sat_no);
            InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
            return flag;
        }

#if 1
        if(phase_v.size()){
            auto max_phase_v=max_element(begin(phase_v),end(phase_v));
            int idx_max_phase_v=distance(begin(phase_v),max_phase_v);
            idx_sat=vflag_[idx_phase_v[idx_max_phase_v]]>>8&0xFF;
            f=vflag_[idx_phase_v[idx_max_phase_v]]&0xF;
            el=epoch_sat_info_collect_[idx_sat].el_az[0];
            if(*max_phase_v>0.03/sin(el)){
                int sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
                int ia=para_.IndexAmb(f,sat_no);
                InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                return true;
            }

            auto max_norm_phase_v=max_element(begin(norm_phase_v),end(norm_phase_v));
            int idx_max_norm_phase_v=distance(begin(norm_phase_v),max_norm_phase_v);
            idx_sat=vflag_[idx_phase_v[idx_max_norm_phase_v]]>>8&0xFF;
            f=vflag_[idx_phase_v[idx_max_norm_phase_v]]&0xF;
            if(*max_norm_phase_v>2.0){
                int sat_no=epoch_sat_info_collect_[idx_sat].sat.sat_.no;
                int ia=para_.IndexAmb(f,sat_no);
                InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                return true;
            }
        }

#endif

        return flag;
    }

    bool cPppSolver::PppResIGG3Control(int iter, vector<double> omcs, vector<double> R) {
        double s0;
        MatMul("NT",1,1,omcs.size(),1.0/(omcs.size()-1),omcs.data(),omcs.data(),0.0,&s0);

        int i, sat1,f;
        cSat sat;
        double k0=1.5,k1=3.0,fact;
        bool flag=false;
        double obs_type;

        vector<double>norm_omcs;
        for(i=0;i<omcs.size();i++){
            norm_omcs.push_back(fabs(omcs[i])/SQRT(R[i]));
        }

        auto max_v_norm=max_element(begin(norm_omcs),end(norm_omcs));
        auto idx_max_v_norm=distance(begin(norm_omcs),max_v_norm);
        sat1=vflag_[idx_max_v_norm]>>8&0xFF;
        sat=cSat(sat1);
        sat.SatNo2Id();
        obs_type=vflag_[idx_max_v_norm]>>4&0xF;
        f=vflag_[idx_max_v_norm]&0xF;

        if(*max_v_norm>=k1){
            // zero-weight segment
            if(obs_type==GNSS_OBS_CODE) previous_sat_info_[sat1-1].c_var_factor[f]=10000;
            else if(obs_type==GNSS_OBS_PHASE) previous_sat_info_[sat1-1].p_var_factor[f]=10000;
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<sat.sat_.id<<" "<<(obs_type?"L":"P")<<f+1<<" RESIDUAL OUTLIER DETECTED, v="<<omcs[idx_max_v_norm]<<" norm_v="<<*max_v_norm<<" ZERO-WEIGHT SEGMENT, factor="<<10000.0;
            flag=true;
        }
        else if(*max_v_norm>=k0){
            // weight-reduced segment
            fact=(k0/fabs(*max_v_norm))*SQR((k1-*max_v_norm)/(k1-k0));

            if(obs_type==GNSS_OBS_CODE) previous_sat_info_[sat1-1].c_var_factor[f]=1.0/fact;
            else if(obs_type==GNSS_OBS_PHASE) previous_sat_info_[sat1-1].p_var_factor[f]=1.0/fact;

            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<sat.sat_.id<<" "<<(obs_type?"L":"P")<<f+1<<" RESIDUAL OUTLIER DETECTED, v="<<omcs[idx_max_v_norm]<<" norm_v="<<*max_v_norm<<" WEIGHT-REDUCED SEGMENT, factor="<<1.0/fact;
            flag=true;
        }
        else{
            // weight-reserved segment
            if(obs_type==GNSS_OBS_CODE) previous_sat_info_[sat1-1].c_var_factor[f]=1.0;
            else if(obs_type==GNSS_OBS_PHASE) previous_sat_info_[sat1-1].p_var_factor[f]=1.0;
        }

        if(iter>3) flag=false;
        return flag;
    }

    void cPppSolver::AverageLcAmb() {
        tSatInfoUnit* sat_info= nullptr;
        tLcAmb* amb= nullptr;
        double LC1,LC2,LC3,var1,var2,var3;

        int i,j;
        for(i=0;i<epoch_sat_info_collect_.size();i++){
            sat_info=&epoch_sat_info_collect_.at(i);
            amb=&sat_info->lc_amb;

            if(sat_info->el_az[0]*R2D<ppp_conf_.gnssC.ele_min) continue;
            LC1=gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,1,-1,0,*sat_info,GNSS_OBS_PHASE,nullptr)-gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,1,1,0,*sat_info,GNSS_OBS_CODE,&var1);
            LC2=gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,0,1,-1,*sat_info,GNSS_OBS_PHASE,nullptr)-gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,0,1,1,*sat_info,GNSS_OBS_CODE,&var2);
            LC3=gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,1,-6,5,*sat_info,GNSS_OBS_PHASE,nullptr)-gnss_obs_operator_.GnssObsLinearComb(ppp_conf_,1,1,0,*sat_info,GNSS_OBS_CODE,&var3);

            if(sat_info->slip[0]||sat_info->slip[1]||sat_info->slip[2]||amb->n[0]==0.0||
                fabs(amb->t_tag[0].TimeDiff(epoch_sat_info_collect_[0].t_tag.t_))>300){
                for(j=0;j<3;j++){
                    amb->n[j]=0.0;
                    amb->lc_amb[j]=0.0;
                    amb->var_amb[j]=0.0;
                }
                amb->fix_cnt=0;
                for(j=0;j<MAX_SAT_NUM;j++) amb->flags[j]=0;

                if(LC1){
                    amb->n[0]+=1;
                    amb->lc_amb[0]+=(LC1-amb->lc_amb[0])/amb->n[0];
                    amb->var_amb[0]+=(var1-amb->var_amb[0])/amb->n[0];
                    amb->t_tag[0]=epoch_sat_info_collect_[0].t_tag;
                }
                if(LC2){
                    amb->n[1]+=1;
                    amb->lc_amb[1]+=(LC1-amb->lc_amb[1])/amb->n[1];
                    amb->var_amb[1]+=(var1-amb->var_amb[1])/amb->n[1];
                    amb->t_tag[1]=epoch_sat_info_collect_[1].t_tag;
                }
                if(LC3){
                    amb->n[2]+=1;
                    amb->lc_amb[2]+=(LC1-amb->lc_amb[2])/amb->n[2];
                    amb->var_amb[2]+=(var1-amb->var_amb[2])/amb->n[2];
                    amb->t_tag[2]=epoch_sat_info_collect_[2].t_tag;
                }
            }
        }
    }

    bool cPppSolver::ResolvePppAmb(int nf, double *xa) {
        double elmask;
        int i,j,m=0;
        if(num_valid_sat_<=0||ppp_conf_.gnssC.ion_opt!=ION_IF||nf<2) return false;
        elmask=ppp_conf_.gnssC.ar_el_mask>0.0?ppp_conf_.gnssC.ar_el_mask:ppp_conf_.gnssC.ele_min;
        vector<int>sat_no1,sat_no2,fix_wls;

        AverageLcAmb();

        tSatInfoUnit *sat_info_i=nullptr,*sat_info_j= nullptr;
        for(i=0;i<epoch_sat_info_collect_.size()-1;i++){
            sat_info_i=&epoch_sat_info_collect_.at(i);
            for(j=i+1;j<epoch_sat_info_collect_.size();j++){
                sat_info_j=&epoch_sat_info_collect_.at(j);
                if(!sat_info_i->vsat[0]||!sat_info_j->vsat[0]||sat_info_i->el_az[0]*R2D<elmask||sat_info_j->el_az[0]*R2D<elmask) continue;

                sat_no1.push_back(sat_info_i->sat.sat_.no);
                sat_no2.push_back(sat_info_j->sat.sat_.no);
                fix_wls.push_back(0);
                if(FixWlAmb(sat_no1.back(),sat_no2.back(),&fix_wls.back())) m++;
            }
        }

        if(ppp_conf_.gnssC.ar_mode==AR_PPP_AR){

        }
        else if(ppp_conf_.gnssC.ar_mode==AR_PPP_AR_ILS){

        }

    }

    bool cPppSolver::FixWlAmb(int sat1,int sat2,int *sd_fix_wl) {
        double sd_wl_amb,var_wl,lam_wl=gnss_obs_operator_.LinearCombLam(1,-1,0,previous_sat_info_[sat1-1]);
        tLcAmb *lc_amb1=&previous_sat_info_[sat1-1].lc_amb;
        tLcAmb *lc_amb2=&previous_sat_info_[sat1-1].lc_amb;

        if(!lc_amb1->n[0]||!lc_amb2->n[0]) return false;

        sd_wl_amb=(lc_amb1->lc_amb[0]-lc_amb2->lc_amb[0])/lam_wl+nav_.wide_line_bias[sat1-1]-nav_.wide_line_bias[sat2-1];

        *sd_fix_wl=Round(sd_wl_amb);
        var_wl=(lc_amb1->var_amb[0]/lc_amb1->n[0]+lc_amb2->var_amb[0]/lc_amb2->n[0])/SQR(lam_wl);

        return fabs(*sd_fix_wl-sd_wl_amb)<=ppp_conf_.gnssC.ar_thres[2]&&
        IntAmbConfFunc(*sd_fix_wl,sd_wl_amb,sqrt(var_wl))>=ppp_conf_.gnssC.ar_thres[1];
    }

    bool cPppSolver::FixNlAmbRound(int *sat1, int *sat2, int *fix_wls, int n) {

    }

    bool cPppSolver::FixNlAmbILS(int *sat1, int *sat2, int *fix_wls, int n) {
        double lam_nl=gnss_obs_operator_.LinearCombLam(1,1,0,previous_sat_info_[sat1[0]-1]);
        double lam1,lam2,alpha,beta;
        int i,j,k,m=0,flag[MAX_SAT_NUM]={0},max_flag=0,info;

        VectorXd B1(n),NC(n),s;
        MatrixXd N1(n,2),D(num_full_x_,n),E(n,num_full_x_),Q(n,n);

        for(i=0;i<n;i++){
            if(!AmbLinearDependCheck(sat1[i],sat2[i],flag,&max_flag)) continue;
            lam1=previous_sat_info_[sat1[i]-1].lam[0];
            lam2=previous_sat_info_[sat1[i]-1].lam[1];
            alpha=SQR(lam2)/(SQR(lam2)-SQR(lam1));
            beta=-SQR(lam1)/(SQR(lam2)-SQR(lam1));
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);

            B1[m]=(full_x_[j]-full_x_[k]+beta*lam2*fix_wls[i])/lam_nl;
            N1(0,m)=Round(B1[m]);

            if(fabs(N1(0,m)-B1[m])>ppp_conf_.gnssC.ar_thres[2]) continue;

            D.data()[j+m*num_full_x_]=1.0/lam_nl;
            D.data()[k+m*num_full_x_]=-1.0/lam_nl;
            sat1[m]=sat1[i];
            sat2[m]=sat2[i];
            fix_wls[m++]=fix_wls[i];
        }
        if(m<3) return false;

        MatMul("TN",m,num_full_x_,num_full_x_,1.0,D.data(),full_Px_.data(),0.0,E.data());
        MatMul("NN",m,m,num_full_x_,1.0,E.data(),D.data(),0.0,Q.data());

        if((info=lambda_.IntegerAmb(B1,Q,N1,m,2,s))){
            return false;
        }
        if(s[0]<=0.0) return false;

        ppplib_sol_.ratio=MIN(s[1]/s[0],999.9);

        // ratio test
        if(ppp_conf_.gnssC.ar_thres[0]>0.0&&ppp_conf_.gnssC.ar_thres[0]>ppplib_sol_.ratio){
            return false;
        }

        // nl ro iono-free ambiguity
        for(i=0;i<m;i++){
            NC[i]=alpha*lam1*N1(i,0)+beta*lam2*(N1(i,0)-fix_wls[i]);
        }

        return FixPppSol(sat1,sat2,NC.data(),m);
    }

    bool cPppSolver::AmbLinearDependCheck(int sat1, int sat2, int *flag, int *max_flg) {
        int i;

        if(flag[sat1-1]==0&&flag[sat2-1]==0){
            flag[sat1-1]=flag[sat2-1]=++(*max_flg);
        }
        else if(flag[sat1-1]==0&&flag[sat2-1]!=0){
            flag[sat1-1]=flag[sat2-1];
        }
        else if(flag[sat1-1]!=0&&flag[sat2-1]==0){
            for(i=0;i<MAX_SAT_NUM;i++) if(flag[i]==flag[sat2-1]) flag[i]=flag[sat1-1];
        }
        else if(flag[sat1-1]>flag[sat2-1]){
            for(i=0;i<MAX_SAT_NUM;i++) if(flag[i]==flag[sat1-1]) flag[i]=flag[sat2-1];
        }
        else if(flag[sat1-1]<flag[sat2-1]){

        }
        else return false;

        return true;
    }

    bool cPppSolver::FixPppSol(int *sat1, int *sat2, double *NC, int n) {
        VectorXd v(n);
        MatrixXd H(num_full_x_,n),R(n,n);
        int i,j,k;
        if(n<=0) return false;

        // constraints to fixed ambiguities
        for(i=0;i<n;i++){
            j=para_.IndexAmb(0,sat1[i]);
            k=para_.IndexAmb(0,sat2[i]);
            v[i]=NC[i]-(full_x_[j]-full_x_[k]);
            H.data()[j+i*num_full_x_]=1.0;
            H.data()[k+i*num_full_x_]=-1.0;
            R.data()[i+i*n]=SQR(0.0001);
        }

        kf_.Adjustment(v,H,R,full_x_,full_Px_,n,num_full_x_);

        for(i=0;i<num_real_x_fix_;i++){
            real_x_fix_[i]=full_x_[i];
            for(j=0;j<num_real_x_fix_;j++){
                real_Px_fix_.data()[i+j*num_real_x_fix_]=real_Px_fix_.data()[j+i*num_real_x_fix_]=full_Px_.data()[i+j*num_full_x_];
            }
        }

        for(i=0;i<n;i++){
            previous_sat_info_[sat1[i]-1].lc_amb.flags[sat2[i]-1]=1;
            previous_sat_info_[sat1[i]-1].lc_amb.flags[sat1[i]-1]=1;
        }
        return true;
    }

    cPpkSolver::cPpkSolver() {}

    cPpkSolver::cPpkSolver(tPPPLibConf C) {
        ppk_conf_=C;
        ppk_conf_.mode=MODE_PPK;
        para_=cParSetting(ppk_conf_);
        num_full_x_=para_.GetPPPLibPar(ppk_conf_);
        full_x_=VectorXd::Zero(num_full_x_);
        full_Px_=MatrixXd::Zero(num_full_x_,num_full_x_);

        num_real_x_fix_=para_.GetRealFixParNum(ppk_conf_);
        real_x_fix_=VectorXd::Zero(num_real_x_fix_);
        real_Px_fix_=MatrixXd::Zero(num_real_x_fix_,num_real_x_fix_);
    }

    cPpkSolver::~cPpkSolver() {}

    void cPpkSolver::InitSolver(tPPPLibConf C) {
        ppk_conf_=C;
        cReadGnssObs base_reader(C.fileC.base,nav_,base_obs_,REC_BASE);
        base_reader.SetGnssSysMask(C.gnssC.nav_sys);
        base_reader.Reading();

        para_=cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
        out_=new cOutSol(ppk_conf_);
        out_->InitOutSol(ppk_conf_,ppk_conf_.fileC.sol);
        out_->WriteHead();

        tPPPLibConf spp_conf=C;
        spp_conf.mode=MODE_SPP;
        spp_conf.mode_opt=MODE_OPT_KINEMATIC;
        spp_conf.gnssC.ion_opt=ION_KLB;
        spp_conf.gnssC.trp_opt=TRP_SAAS;
        spp_conf.gnssC.frq_opt=FRQ_SINGLE;
        spp_conf.gnssC.res_qc=0;
        spp_solver_=new cSppSolver(spp_conf);
        spp_solver_->spp_conf_=spp_conf;
        spp_solver_->nav_=nav_;
        spp_solver_->InitSolver(spp_conf);
    }

    bool cPpkSolver::SolverProcess(tPPPLibConf C,int idx) {

        if(idx==-1) InitSolver(C);


#if 0
        //cpt0
        Vector3d blh(34.220254297*D2R,117.143996963*D2R,36.051);
        base_xyz_=Blh2Xyz(blh);
#endif

#if 0
        // m39
        base_xyz_<<-2267796.9640,5009421.6975,3220952.5435;
#endif

#if 0
      base_xyz_<<-2364335.6607,4870281.4902,-3360816.4321;
      Vector3d blh=Xyz2Blh(base_xyz_);
      cout<<blh.transpose()<<endl;
#endif
        int i=idx,num_epochs=rover_obs_.epoch_num;
        if(idx!=-1){
            num_epochs=idx+1;
        }

        for(i=idx==-1?0:idx;i<num_epochs;i++){

            epoch_sat_obs_=rover_obs_.GetGnssObs().at(i);
            LOG(DEBUG)<<"START PPK SOLVING "<<i+1<<"th EPOCH, ROVER SATELLITE NUMBER "<<epoch_sat_obs_.sat_num;

            if(MatchBaseObs(epoch_sat_obs_.obs_time)){

                LOG(DEBUG)<<"MATCH BASE STATION OBSERVATIONS, BASE SATELLITE NUMBER "<<base_epoch_sat_obs_.sat_num;

                if(SolverEpoch()){
                    SolutionUpdate();
                    if(idx!=-1){
                        epoch_sat_info_collect_.clear();
                        base_sat_info_collect_.clear();
                        rover_res.clear();base_res.clear();
                        return true;
                    }
                }

                if(idx!=-1){
                    epoch_sat_info_collect_.clear();
                    base_sat_info_collect_.clear();
                    rover_res.clear();base_res.clear();
                    return false;
                }
            }
            else{
                ppplib_sol_.stat=SOL_NONE;
                ppplib_sol_.t_tag=epoch_sat_obs_.obs_time;
                epoch_fail_++;
                LOG(WARNING)<<"MATCH BASE STATION OBSERVATIONS FAILED";
            }

            epoch_sat_info_collect_.clear();
            base_sat_info_collect_.clear();
            rover_res.clear();base_res.clear();
        }

        LOG(INFO)<<"TOTAL EPOCH: "<<rover_obs_.epoch_num<<", SOLVE SUCCESS EPOCH: "<<epoch_ok_<<", SOLVE FAILED EPOCH: "<<epoch_fail_;
    }

    bool cPpkSolver::SolverEpoch() {
        epoch_idx_+=1;
        base_xyz_=ppk_conf_.gnssC.rb;

        UpdateGnssObs(ppk_conf_,epoch_sat_obs_,REC_ROVER);
        UpdateGnssObs(ppk_conf_,base_epoch_sat_obs_,REC_BASE);
        InitEpochSatInfo(epoch_sat_info_collect_);

        InitSppSolver();
        if(spp_solver_->SolverEpoch()){
            Spp2Ppk();

            LOG(DEBUG)<<"PPK-SPP SOLVE SUCCESS, SPP POS "<<spp_solver_->ppplib_sol_.pos.transpose()<<" DOPPLER  VEL "<<spp_solver_->ppplib_sol_.vel.transpose();

            gnss_err_corr_.eph_model_.EphCorr(base_sat_info_collect_);
            if(ppk_conf_.gnssC.eph_opt==EPH_PRE){
                gnss_err_corr_.eph_model_.EphCorr(epoch_sat_info_collect_);
            }

            if(Estimator(ppk_conf_)){
                epoch_ok_++;
                LOG(DEBUG)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPK SOLVE SUCCESS";
                return true;
            }else{
                epoch_fail_++;
                ppplib_sol_=spp_solver_->ppplib_sol_;
                LOG(WARNING)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPK SOLVE FAILED";
                return false;
            }
        }
        else{
            epoch_fail_++;
            ppplib_sol_.stat=SOL_NONE;
            ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            spp_solver_->epoch_sat_info_collect_.clear();
            LOG(WARNING)<<epoch_sat_obs_.obs_time.GetTimeStr(1)<<" PPK-SPP SOLVE FAILED";
            return false;
        }

    }

    bool cPpkSolver::Estimator(tPPPLibConf C) {

        // select common satellites
        vector<int>ir,ib,cmn_sat_no;

        SelectCmnSat(C,ir,ib,cmn_sat_no);
        // cycle slip
        PpkCycleSlip(C,ir,ib,cmn_sat_no);
        // state time_update
        StateTimeUpdate(C,ir,ib,cmn_sat_no);

        // zero residual for base and rover
        if(!GnssZeroRes(C,REC_BASE,ib,full_x_.data())){
            return false;
        }


        Vector3d rover_xyz(full_x_[0],full_x_[1],full_x_[2]);
        static int ref_sat[NSYS][2*MAX_GNSS_USED_FRQ_NUM]={0};
        VectorXd x;
        for(int iter=0;iter<3;iter++){
            x=full_x_;
            MatrixXd Px=full_Px_;
            if(!GnssZeroRes(C,REC_ROVER,ir,x.data())){
                break;
            }

            if(!GnssDdRes(0,C,ir,ib,cmn_sat_no,x.data(),ref_sat)){
                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"MAKE DOUBLE DIFFERENCE RESIDUAL FAILED";
                return false;
            }

            if(!kf_.Adjustment(omc_L_,H_,R_,x,Px,num_L_,num_full_x_)){
                ppplib_sol_.stat=SOL_SPP;
                return false;
            }

            if(GnssZeroRes(C,REC_ROVER,ir,x.data())){
                if(GnssDdRes(1,C,ir,ib,cmn_sat_no,x.data(),ref_sat)){
                    full_x_=x;
                    full_Px_=Px;
                    ppplib_sol_.stat=SOL_FLOAT;
                    UpdateSatInfo(epoch_sat_info_collect_);
                    break;
                }

            }
            else{
                ppplib_sol_.stat=SOL_SPP;
                return false;
            }
        }


        //PPK-AR
#if 1
        VectorXd xa=x;
        if(ppplib_sol_.stat==SOL_FLOAT){
            if(ResolvePpkAmb(cmn_sat_no,para_.GetGnssUsedFrqs(),xa.data())){

                if(GnssZeroRes(C,REC_ROVER,ir,xa.data())){
                    GnssDdRes(5,C,ir,ib,cmn_sat_no,xa.data(),ref_sat);
                    ppplib_sol_.stat=SOL_FIX;
                }
            }
            else{
                int a=1;
            }
        }

#endif
        int sat_no,ia;
        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            for(int f=0;f<MAX_GNSS_USED_FRQ_NUM;f++){
                if(!epoch_sat_info_collect_[i].vsat[f]) continue;
                sat_no=epoch_sat_info_collect_[i].sat.sat_.no;
                ia=para_.IndexAmb(f,sat_no);
                previous_sat_info_[sat_no-1].outc[f]=0;
                if(previous_sat_info_[sat_no-1].lock[f]<0||previous_sat_info_[sat_no-1].fix[f]>=2){
                    previous_sat_info_[sat_no-1].lock[f]++;
                }
            }
        }

        ir.clear();ib.clear();cmn_sat_no.clear();ddambs_.clear();vflag_.clear();
        return ppplib_sol_.stat;
    }

    bool cPpkSolver::SolutionUpdate() {
        if(ppplib_sol_.stat==SOL_SPP) return false;
        else if(ppplib_sol_.stat==SOL_FLOAT){
            ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            int num_pos=para_.NumPos();
            for(int i=0;i<num_pos;i++) ppplib_sol_.pos[i]=full_x_[i];
            for(int i=0;i<num_pos;i++) ppplib_sol_.q_pos[i]=full_Px_(i,i);
        }
        else if(ppplib_sol_.stat==SOL_FIX){
            ppplib_sol_.t_tag=epoch_sat_info_collect_[0].t_tag;
            int num_pos=para_.NumPos();
            for(int i=0;i<num_pos;i++) ppplib_sol_.pos[i]=real_x_fix_[i];
            for(int i=0;i<num_pos;i++) ppplib_sol_.q_pos[i]=full_Px_(i,i);
        }
        for(int i=0;i<3;i++) spp_solver_->full_x_[i]=full_x_[i];
        ppplib_sol_.valid_sat_num=num_valid_sat_;

        if(ppk_conf_.solC.out_sol) out_->WriteSol(ppplib_sol_,epoch_idx_);
        if(ppk_conf_.solC.out_stat) out_->WriteSatStat(&ppplib_sol_,previous_sat_info_);
    }

    void cPpkSolver::InitSppSolver() {
        spp_solver_->epoch_idx_=epoch_idx_;
        spp_solver_->epoch_sat_info_collect_=epoch_sat_info_collect_;
    }

    void cPpkSolver::Spp2Ppk() {
        spp_solver_->SolutionUpdate();
        epoch_sat_info_collect_.clear();
        epoch_sat_info_collect_=spp_solver_->epoch_sat_info_collect_;
        spp_solver_->epoch_sat_info_collect_.clear();
    }

    bool cPpkSolver::GnssZeroRes(tPPPLibConf C, RECEIVER_INDEX rec,vector<int>sat_idx,double* x) {

        vector<tSatInfoUnit> *sat_collect;
        vector<double> *res;

        sat_collect=rec==REC_ROVER?&epoch_sat_info_collect_:&base_sat_info_collect_;
        res=rec==REC_ROVER?&rover_res:&base_res;
        double omc,r,meas;
        tSatInfoUnit* sat_info= nullptr;

        Vector3d rec_xyz;
        if(rec==REC_ROVER){
            rec_xyz<<x[0],x[1],x[2];
        }
        else{
            rec_xyz=base_xyz_;
        }

        Vector3d rover_blh=Xyz2Blh(rec_xyz);

        int num_used_frq=para_.GetGnssUsedFrqs();
        int num_used_obs_type=para_.GetNumObsType();
        if(num_used_frq<=0) return false;
        int nf=num_used_frq*num_used_obs_type;

        rover_res.resize(sat_idx.size()*num_used_frq*num_used_obs_type,0.0);
        base_res.resize(sat_idx.size()*num_used_frq*num_used_obs_type,0.0);

        for(int i=0;i<sat_idx.size();i++){

            sat_info=&sat_collect->at(sat_idx[i]);
            if(sat_info->stat!=SAT_USED){
                LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" NO USED("<<kGnssSatStatStr[sat_info->stat+1]<<")";
                continue;
            }

            if((r=GeoDist(sat_info->brd_pos,rec_xyz,sat_info->sig_vec))<=0.0){
                LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" NO USED";
                continue;
            }

            if(SatElAz(rover_blh,sat_info->sig_vec,sat_info->el_az)<C.gnssC.ele_min*D2R){
                LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" LOW ELE";
                continue;
            }

            r+=gnss_err_corr_.SagnacCorr(sat_info->brd_pos,rec_xyz);
//            gnss_err_corr_.trp_model_.InitSatInfo(sat_info,&rover_blh);
//            gnss_err_corr_.trp_model_.GetTrpError(0.0,x,para_.IndexTrp());
//            gnss_err_corr_.trp_model_.UpdateSatInfo();
//            gnss_err_corr_.ion_model_.InitSatInfo(sat_info,&rover_blh);
//            gnss_err_corr_.ion_model_.GetIonError();
//            gnss_err_corr_.ion_model_.UpdateSatInfo();
//            gnss_err_corr_.ant_model_.SatPcvCorr(sat_info,rover_blh, nullptr);
//            gnss_err_corr_.ant_model_.RecAntCorr(sat_info, nullptr,REC_ROVER);

            int obs_code,frq;
            for(int f=0;f<num_used_frq*num_used_obs_type;f++){
                //f%num_used_obs_type==0, L, ==1 P, ==2 D,
                //f/num_used_obs_type==1, f1,==2 f2,==3 f3
                frq=f>=num_used_frq?f-num_used_frq:f;
                obs_code=f<num_used_frq?0:1;

                //ionosphere-free
                if(C.gnssC.ion_opt==ION_IF){
//                    meas=obs_type==GNSS_OBS_CODE?sat_info->cor_if_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->cor_if_L[frq]:sat_info->cor_D[frq]);
                    meas=obs_code?sat_info->cor_if_P[frq]:sat_info->cor_if_L[frq];
                }
                    //uncombined and undifference
                else{
                    meas=obs_code?sat_info->raw_P[frq]:sat_info->raw_L[frq]*sat_info->lam[frq];
//                    meas=obs_type==GNSS_OBS_CODE?sat_info->raw_P[frq]:(obs_type==GNSS_OBS_PHASE?sat_info->raw_L[frq]*sat_info->lam[frq]:sat_info->cor_D[frq]);
                }

                if(meas==0.0) {
                    LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<" "<<(rec==REC_BASE?"BASE":"ROVER")<<" "<<sat_info->sat.sat_.id<<" MISSING "<<(obs_code?"P":"L")<<frq+1;
                    continue;
                }

                double sat_clk=C.gnssC.eph_opt==EPH_BRD?sat_info->brd_clk[0]*CLIGHT:sat_info->pre_clk[0]*CLIGHT;
                double trp_del=sat_info->trp_dry_delay[0]+sat_info->trp_wet_delay[0];
                double ion_del=sat_info->ion_delay[0];
                omc=meas-(r-sat_clk);
                res->at(i*num_used_frq*num_used_obs_type+f)=omc;
            }
        }

        res= nullptr;
        return true;
    }

    int cPpkSolver::GnssDdRes(int post, tPPPLibConf C,vector<int>ir,vector<int>ib,vector<int>cmn_sat_no, double *x,int refsat[NSYS][2*MAX_GNSS_USED_FRQ_NUM]) {
        num_L_=0;
        num_valid_sat_=0;
        if(vflag_.size()) vflag_.clear();
        char buff[MAX_BUFF]={'\0'};

        int num_used_frq=para_.GetGnssUsedFrqs();
        int num_used_obs_type=para_.GetNumObsType();
        if(num_used_frq<=0) return false;

        int nf=num_used_frq*num_used_obs_type,ref_sat_idx,j,ia,ref_ia,need_iter=1;
        int frq,obs_code,nobs[NSYS][nf],iobs=0;
        double omc,ambi,ambj,threshadj;
        vector<double> omcs,H,Ri,Rj,R;
        tSatInfoUnit* sat_info= nullptr, *ref_sat_info= nullptr;

        for(int i=0;i<NSYS;i++){
            for(int k=0;k<nf;k++){
                nobs[i][k]=0;
            }
        }

        for(int isys=0;isys<NSYS;isys++){

            for(int ifrq=0;ifrq<nf;ifrq++){
                ambi=0.0;ambj=0.0;
                frq=ifrq>=num_used_frq?ifrq-num_used_frq:ifrq;
                obs_code=ifrq<num_used_frq?0:1;

                for(j=0,ref_sat_idx=-1;j<ir.size();j++){
                    int sys=epoch_sat_info_collect_.at(ir[j]).sat.sat_.sys_idx;
                    if(sys!=isys) continue;
                    if(!ValidObs(j,num_used_frq,ifrq)) continue;
                    if(ref_sat_idx<0||epoch_sat_info_collect_.at(ir[ref_sat_idx]).el_az[0]<epoch_sat_info_collect_.at(ir[j]).el_az[0]) ref_sat_idx=j;
                }
                if(ref_sat_idx<0) continue;
                LOG(DEBUG)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<kGnssSysStr[isys+1]<<(post?" POST":" PRIOR")<<" DOUBLE DIFFERENCE, "<<" REFERENCE SATELLITE "<<epoch_sat_info_collect_[ir[ref_sat_idx]].sat.sat_.id<<" "<<epoch_sat_info_collect_[ir[ref_sat_idx]].el_az[0]*R2D;

                if(refsat){
                    refsat[isys][ifrq]=cmn_sat_no[ref_sat_idx];
                }

                for(int isat=0;isat<cmn_sat_no.size();isat++){
                    sat_info=&epoch_sat_info_collect_[ir[isat]];
                    ref_sat_info=&epoch_sat_info_collect_[ir[ref_sat_idx]];

                    if(isat==ref_sat_idx) {
                        if(post==1&&ifrq<num_used_frq){
                            ref_sat_info->outc[frq]=0;
                            ref_sat_info->lock[frq]++;
                            ref_sat_info->vsat[frq]=1;
                        }
                        continue;
                    }

                    if(sat_info->stat!=SAT_USED) continue;
                    int sys=epoch_sat_info_collect_.at(ir[isat]).sat.sat_.sys_idx;
                    if(sys!=isys) continue;
                    if(!ValidObs(isat,num_used_frq,ifrq)){
                        continue;
                    }

                    if(C.gnssC.check_dual_phase&&C.gnssC.frq_opt>=FRQ_DUAL&&!obs_code){
                        if(!CheckDualFrq(*sat_info)){
                            LOG(DEBUG)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" NO DUAL PHASE OBSERVATIONS";
                            continue;
                        }
                    }


                    omc=(rover_res[ref_sat_idx*nf+ifrq]-base_res[ref_sat_idx*nf+ifrq])-(rover_res[isat*nf+ifrq]-base_res[isat*nf+ifrq]);

                    if(!obs_code){
                        ia=para_.IndexAmb(frq,sat_info->sat.sat_.no);
                        ref_ia=para_.IndexAmb(frq,ref_sat_info->sat.sat_.no);
                        if(C.gnssC.ion_opt!=ION_IF){
                            ambi=ref_sat_info->lam[frq]*x[ref_ia];ambj=sat_info->lam[frq]*x[ia];
                        }
                        else{

                        }
                        omc-=(ambi-ambj);
                    }

                    threshadj=obs_code||(full_Px_(ia,ia)>SQR(30/2))||(full_Px_(ref_ia,ref_ia)>SQR(30/2))?300:1.0;
                    if(C.gnssC.max_inno>0.0&&fabs(omc)>C.gnssC.max_inno*threshadj){
                        sat_info->rejc[frq]++;
                        LOG(WARNING)<<sat_info->t_tag.GetTimeStr(1)<<" "<<sat_info->sat.sat_.id<<" "<<(obs_code?"P":"L")<<frq+1<<" INNOVATION OVERRUN inno="<<omc<<" threshold="<<threshadj;
                        continue;
                    }

                    omcs.push_back(omc);
                    Ri.push_back(GnssMeasVar(C,(obs_code?GNSS_OBS_CODE:GNSS_OBS_PHASE),*ref_sat_info));
                    Rj.push_back(GnssMeasVar(C,(obs_code?GNSS_OBS_CODE:GNSS_OBS_PHASE),*sat_info));
                    if(!post){
                        //position
                        for(int i=0;i<num_full_x_;i++) H.push_back(i<3?(-ref_sat_info->sig_vec[i]+sat_info->sig_vec[i]):0.0);
                        //tropospheric delay
                        if(C.gnssC.ion_opt>=ION_EST){

                        }
                        //ionospheric delay
                        if(C.gnssC.trp_opt>=TRP_EST_WET){

                        }
                        //ambiguity
                        if(!obs_code){
                            if(C.gnssC.ion_opt!=ION_IF){
                                H[ref_ia+num_full_x_*num_L_]=ref_sat_info->lam[frq];
                                H[ia+num_full_x_*num_L_]=-sat_info->lam[frq];
                            }
                            else{
                                H[ref_ia+num_full_x_*num_L_]=1.0;
                                H[ia+num_full_x_*num_L_]=-1.0;
                            }
                        }
                        if(obs_code){
                            sat_info->prior_res_P[frq]=omcs.back();
                        }
                        else if(!obs_code){
                            sat_info->prior_res_L[frq]=omcs.back();
                        }
                    }
                    else{
                        if(obs_code){
                            sat_info->post_res_P[frq]=omcs.back();
                        }
                        else if(!obs_code&&post!=5){
                            sat_info->post_res_L[frq]=omcs.back();
                            if(C.gnssC.ion_opt!=ION_IF){
                                sat_info->float_amb[frq]=ambj/sat_info->lam[frq];
                                ref_sat_info->float_amb[frq]=ambi/ref_sat_info->lam[frq]; //cycle;
                            }
                            else{
                                sat_info->float_amb[frq]=ambj;    //m
                                ref_sat_info->float_amb[frq]=ambi;
                            }
                        }
                    }

                    sat_info->vsat[frq]=1;
                    sat_info->res_idx[frq]=num_L_;
                    previous_sat_info_[sat_info->sat.sat_.no-1].vsat[frq]=1;

                    sprintf(buff,"(%9.5f - %9.5f)-(%9.5f - %9.5f)-(%9.5f - %9.5f)=%9.5f el=%3.1f Ri=%10.6f Rj=%10.6f",rover_res[ref_sat_idx*nf+ifrq],base_res[ref_sat_idx*nf+ifrq],rover_res[isat*nf+ifrq],base_res[isat*nf+ifrq],ambi,ambj,omcs.back(),sat_info->el_az[0]*R2D,Ri.back(),Rj.back());
                    LOG(DEBUG)<<(obs_code?"P":"L")<<frq+1<<": (ROVER_"<<ref_sat_info->sat.sat_.id<<"-BASE_"<<ref_sat_info->sat.sat_.id<<")-(ROVER_"<<sat_info->sat.sat_.id<<"-BASE_"<<sat_info->sat.sat_.id<<")-amb"
                              <<" = "<<buff;

                    vflag_.push_back((ref_sat_info->sat.sat_.no<<16|(sat_info->sat.sat_.no<<8)|(obs_code<<4)|(frq)));
                    num_L_++;
                    if(ifrq==0) num_valid_sat_++;
                    nobs[isys][iobs]++;
                }
                iobs++;
            }
        }

        R.resize(num_L_*num_L_,0.0);
        for(int isys=0,nt=0;isys<NSYS;isys++){
            for(int f=0;f<nf;f++){
                for(int j=0;j<nobs[isys][f];j++){
                    for(int k=0;k<nobs[isys][f];k++){
                        R[nt+k+(nt+j)*num_L_]=j==k?Ri[nt+j]+Rj[nt+j]:Ri[nt+j];
                    }
                }
                nt+=nobs[isys][f];
            }
        }

        if(!post){
            H_=Map<MatrixXd>(H.data(),num_full_x_,num_L_);
            R_=Map<MatrixXd>(R.data(),num_L_,num_L_);
            omc_L_=Map<VectorXd>(omcs.data(),num_L_);
        }

        if(post&&post!=5&&C.gnssC.res_qc==RES_QC_STEP){
//           if(PpkResStepControl(post,ir,cmn_sat_no,omcs,kf_.v_,kf_.Qvv_)){
//               need_iter=0;
//           }
        }

        H.clear();Ri.clear();Rj.clear();R.clear();omcs.clear();
        return post?need_iter:num_L_;
    }

    bool cPpkSolver::ValidObs(int i, int nf, int f) {
        return (rover_res[i*nf*2+f]!=0.0&&base_res[i*nf*2+f]!=0.0&&(f<nf||rover_res[f-nf+i*nf*2]!=0.0&&base_res[f-nf+i*nf*2]!=0.0));
    }

    bool cPpkSolver::MatchBaseObs(cTime t) {
        double sow1,sow2;
        int i,week=0,wod=0,info=false;

        for(i=base_idx_-100<0?0:base_idx_-10;base_idx_<base_obs_.GetGnssObs().size();i++){
            sow1=base_obs_.GetGnssObs().at(i).obs_time.Time2Gpst(&week,&wod,SYS_GPS);
            sow2=t.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<DTTOL){
                base_idx_=i;info=true;
                base_epoch_sat_obs_=base_obs_.GetGnssObs().at(base_idx_);
                break;
            }
            else if((sow1-sow2)>2.0*DTTOL){
                info=false;break;
            }
        }

        return info;
    }

    int cPpkSolver::SelectCmnSat(tPPPLibConf C, vector<int> &ir, vector<int> &ib, vector<int> &cmn_sat_no) {
        string buff;
        tSatInfoUnit rover,base;
        for(int i=0,j=0;i<epoch_sat_info_collect_.size()&&j<base_sat_info_collect_.size();i++,j++){
            rover=epoch_sat_info_collect_.at(i);
            base=base_sat_info_collect_.at(j);
            if(rover.sat.sat_.no<base.sat.sat_.no) j--;
            else if(rover.sat.sat_.no>base.sat.sat_.no) i--;
            else if(rover.el_az[0]>=C.gnssC.ele_min*D2R){
                if(rover.sat.sat_.no!=base.sat.sat_.no) continue;
                ir.push_back(i);
                ib.push_back(j);
                cmn_sat_no.push_back(rover.sat.sat_.no);
            }
        }

        for(int i=0;i<epoch_sat_info_collect_.size();i++){
            buff+=epoch_sat_info_collect_.at(i).sat.sat_.id+" ";
        }
        LOG(DEBUG)<<"PPK ROVER STATION OBSERVED SATELLITES: "<<epoch_sat_info_collect_.size()<<" "<<buff;
        buff.clear();

        for(int i=0;i<base_sat_info_collect_.size();i++){
            buff+=base_sat_info_collect_.at(i).sat.sat_.id+" ";
        }
        LOG(DEBUG)<<"PPK BASE STATION OBSERVED SATELLITES : "<<base_sat_info_collect_.size()<<" "<<buff;
        buff.clear();

        for(int i=0;i<ir.size();i++){
            buff+=epoch_sat_info_collect_.at(ir[i]).sat.sat_.id+" ";
        }

        LOG(DEBUG)<<"PPK ROVER AND BASE COMMON SATELLITE  : "<<ir.size()<<" "<<buff;

        return ir.size();
    }

    void cPpkSolver::PpkCycleSlip(tPPPLibConf C,vector<int>& iu,vector<int>& ib,vector<int>& cmn_sat_no) {

        tSatInfoUnit* sat_info= nullptr;
        tSatInfoUnit* base_sat= nullptr;
        cTime t=ppplib_sol_.t_tag;
        double dt=C.gnssC.sample_rate;
        int f;

        if(t.t_.long_time!=0.0) dt=epoch_sat_info_collect_[0].t_tag.TimeDiff(t.t_);

        for(int i=0;i<cmn_sat_no.size();i++){
            sat_info=&epoch_sat_info_collect_.at(iu[i]);
            base_sat=&base_sat_info_collect_.at(ib[i]);

            gnss_obs_operator_.LliCycleSlip(C,*sat_info,2,dt,REC_ROVER);

            if(gnss_obs_operator_.LliCycleSlip(C,*base_sat,2,dt,REC_BASE)){
                for(int j=0;j<MAX_GNSS_USED_FRQ_NUM;j++) sat_info->slip[j]|=base_sat->slip[j];
            }
            gnss_obs_operator_.MwCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, base_sat,previous_sat_info_[sat_info->sat.sat_.no-1].t_tag.t_);
            gnss_obs_operator_.GfCycleSlip(C,C.gnssC.sample_rate,dt,sat_info, base_sat);

            gnss_obs_operator_.SmoothMw(C,sat_info, base_sat);
        }


    }

    void cPpkSolver::StateTimeUpdate(tPPPLibConf C,vector<int>& ir,vector<int>& ib,vector<int>& cmn_sat_no) {
        double tt=spp_solver_->ppplib_sol_.t_tag.TimeDiff(ppplib_sol_.t_tag.t_);
        //position
        PosUpdate(C);
        //tropospheric delay
        TrpUpdate(C,tt);
        //ionospheric delay
        IonUpdate(C,tt);
        //ambiguity
        AmbUpdate(C,tt,ir,ib,cmn_sat_no);
    }

    void cPpkSolver::PosUpdate(tPPPLibConf  C) {
        if(para_.NumPos()<=0) return;

        Vector3d q(SQR(30),SQR(30),SQR(30));
        int ip=para_.IndexPos();
        if((SQR(full_x_[ip])+SQR(full_x_[ip+1])+SQR(full_x_[ip+1]))==0.0){
            for(int i=0;i<3;i++) InitX(spp_solver_->ppplib_sol_.pos[i],SQR(30.0),ip+i,full_x_.data(),full_Px_.data());
        }else{
            if(C.mode_opt==MODE_OPT_STATIC) return;
            if(C.mode_opt==MODE_OPT_KINE_SIM||C.mode_opt==MODE_OPT_KINEMATIC||C.mode_opt==MODE_OPT_PPK){
                for(int i=0;i<3;i++) InitX(spp_solver_->ppplib_sol_.pos[i],SQR(30.0),ip+i,full_x_.data(),full_Px_.data());
            }
        }
    }

    void cPpkSolver::TrpUpdate(tPPPLibConf C, double tt) {
        if(para_.NumTrp()<=0) return;
    }

    void cPpkSolver::IonUpdate(tPPPLibConf C, double tt) {
        if(para_.NumIon()<=0) return;
    }

    void cPpkSolver::AmbUpdate(tPPPLibConf C, double tt,vector<int>& ir,vector<int>& ib,vector<int>& cmn_sat_no) {
        if(para_.NumAmb()<=0) return;

        int i,j,reset,ia,slip,rejc;
        double offset;
        for(int f=0;f<para_.GetGnssUsedFrqs();f++){

            for(i=0;i<MAX_SAT_NUM;i++){
                ia=para_.IndexAmb(f,i+1);
                reset=previous_sat_info_[i].outc[f]>C.gnssC.max_out;
                if(C.gnssC.ar_mode==AR_INST&&full_x_[ia]!=0.0){
                    InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                }
                else if(reset&&full_x_[ia]!=0.0){
                    InitX(0.0,0.0,ia,full_x_.data(),full_Px_.data());
                    previous_sat_info_[i].outc[f]=0;
                }
                if(C.gnssC.ar_mode!=AR_INST&&reset){

                }
            }

            for(i=0;i<cmn_sat_no.size();i++){
                ia=para_.IndexAmb(f,cmn_sat_no[i]);
                full_Px_(ia,ia)+=SQR(C.gnssC.ait_psd[0])*fabs(tt);
                slip=epoch_sat_info_collect_[ir[i]].slip[f];
                rejc=epoch_sat_info_collect_[ir[i]].rejc[f];
                if(!(slip&1)&&rejc<2) continue;
                full_x_[ia]=0.0;
                epoch_sat_info_collect_[ir[i]].rejc[f]=0;
                epoch_sat_info_collect_[ir[i]].lock[f]=-C.gnssC.min_lock2fix;
            }

            double cp,pr,amb,com_offset;
            vector<double>bias;
            for(i=0,j=0,offset=0.0;i<cmn_sat_no.size();i++){
                if(C.gnssC.ion_opt!=ION_IF){
                    if(C.gnssC.frq_opt==FRQ_DUAL&&C.gnssC.check_dual_phase){
                        if(!CheckDualFrq(epoch_sat_info_collect_[ir[i]])){
                            epoch_sat_info_collect_[ir[i]].stat=SAT_NO_USE;
                            bias.push_back(0.0);
                            continue;
                        }
                    }
                    cp=gnss_obs_operator_.GnssSdObs(epoch_sat_info_collect_[ir[i]],base_sat_info_collect_[ib[i]],f,GNSS_OBS_PHASE);
                    pr=gnss_obs_operator_.GnssSdObs(epoch_sat_info_collect_[ir[i]],base_sat_info_collect_[ib[i]],f,GNSS_OBS_CODE);
                    if(cp==0.0||pr==0.0){
                        if(cp==0.0){
                            LOG(DEBUG)<<epoch_sat_info_collect_[ir[i]].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[i]].sat.sat_.id<<" MISSING L"<<f+1<<" FOR INITIALIZING AMBIGUITY";
                        }
                        if(pr==0.0){
                            LOG(DEBUG)<<epoch_sat_info_collect_[ir[i]].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[i]].sat.sat_.id<<" MISSING P"<<f+1<<" FOR INITIALIZING AMBIGUITY";
//                            epoch_sat_info_collect_[ir[i]].stat=SAT_NO_PR;
                        }
                        bias.push_back(0.0);
                        continue;
                    }
                    amb=cp-pr/epoch_sat_info_collect_[ir[i]].lam[f];
                    bias.push_back(amb);

                }
                else{

                }

                ia=para_.IndexAmb(f,cmn_sat_no[i]);
                if(full_x_[ia]!=0.0){
                    offset+=amb-full_x_[ia]*epoch_sat_info_collect_[ir[i]].lam[f];
                    j++;
                }
            }

            com_offset=j>0?offset/j:0.0;

            for(i=0;i<cmn_sat_no.size();i++){
                ia=para_.IndexAmb(f,cmn_sat_no[i]);
                if(bias[i]==0.0||full_x_[ia]!=0.0) continue;
                InitX((bias[i]-com_offset),SQR(5.0),ia,full_x_.data(),full_Px_.data());
                epoch_sat_info_collect_[ir[i]].lock[f]=-C.gnssC.min_lock2fix;
                LOG(DEBUG)<<epoch_sat_info_collect_.at(ir[i]).t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_.at(ir[i]).sat.sat_.id<<" L"<<f+1<<" AMBIGUITY INITIALIZED "<<(bias[i]-com_offset);
            }
            bias.clear();
        }
    }

    // can improved AR rate
    bool cPpkSolver::CheckDualFrq(tSatInfoUnit &sat_info) {
        for(int i=0;i<2;i++){
            if(sat_info.raw_L[i]==0.0) return false;
        }
        return true;
    }

    // function from carvig
    void cPpkSolver::DetectPhaseOutlier(int post,vector<int>cmn_sat,vector<double> omcs,int refsat[NSYS][2*MAX_GNSS_USED_FRQ_NUM],vector<double>& R) {
        int m,i,j,k,n,f,flag=1,sat1,sat2;
        vector<double>dv;
        vector<int>idx;
        double vv,avg_vv,v0;
        double r0=ReNorm(0.95);

        for(m=0,flag=1;m<cmn_sat.size();m++,flag=1){
#if 1
            for(j=0;j<cmn_sat.size();j++){
                if(!previous_sat_info_[cmn_sat[i]-1].vsat[0]) continue;

                if(refsat[previous_sat_info_[cmn_sat[i]-1].sat.sat_.sys_idx][0]==cmn_sat[i]) continue;
                if(refsat[previous_sat_info_[cmn_sat[i]-1].sat.sat_.sys_idx][1]==cmn_sat[i]) continue;

                k=previous_sat_info_[cmn_sat[i]-1].res_idx[0];
                n=previous_sat_info_[cmn_sat[i]-1].res_idx[1];
                if(((vflag_[k]>>4)&0xF)==1) continue;

                vv=omcs[k]-omcs[n]; //f1-f2
                dv.push_back(vv);
                idx.push_back(i);
            }
#else
#endif
            if(idx.size()<3) break;

            for(i=0;i<idx.size();i++) avg_vv+=dv[i];avg_vv/=idx.size();
            for(i=0;i<idx.size();i++) dv[i]-=avg_vv;

            MatMul("NT",1,1,idx.size(),1.0/(idx.size()-1),dv.data(),dv.data(),0.0,&v0);

            for(i=0;i<idx.size();i++){
                if(fabs(v0)>=0.03&&fabs(dv[i])/SQRT(v0)>=r0){

                    LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<previous_sat_info_[cmn_sat[idx[i]]-1].sat.sat_.id<<" "<<" L1_L2 DD "<<(post?"POST":"PRIOR")<<" RESIDUAL DETECT OUTLIER dv="<<omcs[idx[i]];
                    k=previous_sat_info_[cmn_sat[idx[i]]-1].res_idx[0];
                    n=previous_sat_info_[cmn_sat[idx[i]]-1].res_idx[1];
                    R[k]=R[n]=SQR(100.0);

                    for(f=0;f<MAX_GNSS_USED_FRQ_NUM;f++){
                        previous_sat_info_[cmn_sat[idx[i]]-1].vsat[f]=0;
                    }
                    flag=0;
                }
            }
            if(flag) break;
        }
        idx.clear();dv.clear();
    }

    bool cPpkSolver::PpkResStepControl(int post,vector<int>ir, vector<int> sat_no,vector<double>& omcs, VectorXd& v,MatrixXd& Qvv) {
        vector<double>v_code,v_phase,norm_v_code,norm_v_phase;
        vector<int> code_idx,phase_idx;
        int i,sat1,sat2,f,f1,f2,obs_type,sat_idx;
        bool need_iter=false;

        if(!Qvv.data()||!v.data()) return false;
        for(i=0;i<omcs.size();i++){
            obs_type=(vflag_[i]>>4)&0xF;

            if(obs_type==1){
                code_idx.push_back(i);
                v_code.push_back(fabs(v[i]));
                norm_v_code.push_back(fabs(v[i])/SQRT(fabs(Qvv(i,i))));
            }
            else if(obs_type==0){
                phase_idx.push_back(i);
                v_phase.push_back(fabs(v[i]));
                norm_v_phase.push_back(fabs(v[i])/SQRT(fabs(Qvv(i,i))));
            }
        }

        double thres_code=3.0,thres_norm_code=2.0,el;
        // step-1: test pseudorange residual
        bool code_flag=false;
        auto max_v_code=max_element(begin(v_code),end(v_code));
        int idx_max_v_code=distance(begin(v_code),max_v_code);

        sat1=(vflag_[code_idx[idx_max_v_code]]>>16)&0xFF;
        sat2=(vflag_[code_idx[idx_max_v_code]]>>8)&0xFF;
        f=(vflag_[code_idx[idx_max_v_code]])&0xF;

        int j1;
        for(j1=0;j1<sat_no.size();j1++){
            if(sat_no[j1]==sat2) break;
        }
        el=epoch_sat_info_collect_[ir[j1]].el_az[0];

        if(*max_v_code>thres_code/sin(el)){
            code_flag=true;
        }

        // step-2: test standardized pseudorange residual
        bool norm_code_flag=false;
        auto max_norm_v_code=max_element(begin(norm_v_code),end(norm_v_code));
        int idx_max_norm_v_code=distance(begin(norm_v_code),max_norm_v_code);
        sat1=(vflag_[code_idx[idx_max_norm_v_code]]>>16)&0xFF;
        sat2=(vflag_[code_idx[idx_max_norm_v_code]]>>8)&0xFF;
        f=(vflag_[code_idx[idx_max_norm_v_code]])&0xF;

        if(*max_norm_v_code>thres_norm_code){
            norm_code_flag=true;
        }

        int j2;
        if(norm_code_flag){
            for(j2=0;j2<sat_no.size();j2++){
                if(sat_no[j2]==sat2) break;
            }
        }

        //case 1:
        if(norm_code_flag&&code_flag){
            epoch_sat_info_collect_[j1].stat=SAT_NO_USE;
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" "<<"REJECT BY POST PSEUDORANGE RESIDUAL res="<<*max_v_code<<" el="<<el*R2D<<" thres="<<3.0/sin(el);
        }
        else if(code_flag&&!norm_code_flag){
            epoch_sat_info_collect_[j1].stat=SAT_NO_USE;
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" "<<"REJECT BY POST PSEUDORANGE RESIDUAL res="<<*max_v_code<<" el="<<el*R2D<<" thres="<<3.0/sin(el);
        }
        else if(!code_flag&&norm_code_flag){
            epoch_sat_info_collect_[j2].stat=SAT_NO_USE;
            int aa=code_idx[idx_max_norm_v_code];
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" "<<"REJECT BY POST STANDARDIZED PSEUDORANGE RESIDUAL norm_res="<<*max_norm_v_code<<" res="<<v[aa]<<" sig="<<SQRT(fabs(Qvv(aa,aa)))<<" el="<<el*R2D<<" thres="<<2.0;
        }

        if(code_flag||norm_code_flag) need_iter=true;

        if(!need_iter){
            double thres_phase=0.03,thres_norm_phase=2.0;
            // step-3: test phase residual
            bool phase_flag=false;
            auto max_v_phase=max_element(begin(v_phase),end(v_phase));
            int idx_max_v_phase=distance(begin(v_phase),max_v_phase);

            sat1=(vflag_[phase_idx[idx_max_v_phase]]>>16)&0xFF;
            sat2=(vflag_[phase_idx[idx_max_v_phase]]>>8)&0xFF;
            f1=(vflag_[phase_idx[idx_max_v_phase]])&0xF;

            for(j1=0;j1<sat_no.size();j1++){
                if(sat_no[j1]==sat2) break;
            }
            el=epoch_sat_info_collect_[ir[j1]].el_az[0];

            if(*max_v_phase>thres_phase/sin(el)){
                phase_flag=true;
            }

            // step-4: test standardized phase residual
            bool norm_phase_flag=false;
            auto max_norm_v_phase=max_element(begin(norm_v_phase),end(norm_v_phase));
            int idx_max_norm_v_phase=distance(begin(norm_v_phase),max_norm_v_phase);
            sat1=(vflag_[phase_idx[idx_max_norm_v_phase]]>>16)&0xFF;
            sat2=(vflag_[phase_idx[idx_max_norm_v_phase]]>>8)&0xFF;
            f2=(vflag_[phase_idx[idx_max_norm_v_phase]])&0xF;

            if(*max_norm_v_phase>thres_norm_phase){
                norm_phase_flag=true;
            }

            if(norm_phase_flag){
                for(j2=0;j2<sat_no.size();j2++){
                    if(sat_no[j2]==sat2) break;
                }
            }

            int ia;
            if(phase_flag&&norm_phase_flag){
                ia=para_.IndexAmb(f1,sat_no[j1]);
                InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" L"<<f1+1<<" REINITIALIZE AMBIGUITY BY POST PHASE RESIDUAL res="<<*max_v_phase<<" el="<<el*R2D<<" thres="<<0.03/sin(el);
            }
            else if(phase_flag&&!norm_phase_flag){
//                ia=para_.IndexAmb(f2,sat_no[j1]);
//                InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
//                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j1]].sat.sat_.id<<" L"<<f1+1<<" REINITIALIZE AMBIGUITY BY POST PHASE RESIDUAL res="<<*max_v_phase<<" el="<<el*R2D<<" thres="<<0.03/sin(el);
            }
            else if(!phase_flag&&norm_phase_flag){
//                int aa=phase_idx[idx_max_norm_v_phase];
//                ia=para_.IndexAmb(f2,sat_no[j2]);
//                el=epoch_sat_info_collect_[ir[j2]].el_az[0];
//                InitX(full_x_[ia],SQR(60.0),ia,full_x_.data(),full_Px_.data());
//                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<epoch_sat_info_collect_[ir[j2]].sat.sat_.id<<" L"<<f2+1<<" REINITIALIZE AMBIGUITY BY POST STANDARDIZE PHASE RESIDUAL norm_res="<<*max_norm_v_phase<<" res="<<v_phase[aa]<<" sig="<<SQRT(fabs(Qvv(aa,aa)))<<" "<<" el="<<el*R2D<<" thres="<<2.0;
            }

            if(phase_flag||norm_phase_flag) need_iter=true;
        }


        v_code.clear();v_phase.clear();norm_v_code.clear();norm_v_phase.clear();code_idx.clear();phase_idx.clear();
        return need_iter;
    }

    void cPpkSolver::ReSetAmb(double *bias, double *xa,int nb) {
        int i,n,m,f,index[MAX_SAT_NUM]={0},nv=0,nf=para_.GetGnssUsedFrqs();
        tSatInfoUnit *sat_info= nullptr;

        for(i=0;i<num_real_x_fix_;i++) xa[i]=real_x_fix_[i];

#if 0
        for(m=0;m<NSYS;m++){
            for(f=0;f<nf;f++){
                for(n=i=0;i<MAX_SAT_NUM;i++){
                    sat_info=&previous_sat_info_[i];
                    if(sat_info->sat.sat_.sys_idx!=m||sat_info->fix[f]!=2) continue;
                    index[n++]=para_.IndexAmb(f,sat_info->sat.sat_.no);
                }

                if(n<2) continue;
                xa[index[0]]=full_x_[index[0]];

                for(i=1;i<n;i++){
                    xa[index[i]]=xa[index[0]]-bias[nv++];
                }
            }
        }
#endif
        int sat1,sat2,ib1,ib2;
        for(i=0;i<nb;i++){
            sat1=ddambs_[i].ref_sat;
            sat2=ddambs_[i].sat;
            f=ddambs_[i].f;

            ib1=para_.IndexAmb(f,sat1);
            ib2=para_.IndexAmb(f,sat2);
            // 基准卫星仍是浮点解
            xa[ib2]=xa[ib1]-bias[i];

            previous_sat_info_[sat1-1].fix_amb[f]=xa[ib1];
            previous_sat_info_[sat2-1].fix_amb[f]=-bias[i];
        }
    }

    int cPpkSolver::DdMat(double *D, int gps, int glo) {
        int i,j,k,m,f,nb=0,na=num_real_x_fix_,nf=para_.GetGnssUsedFrqs(),no_fix;
        double fix[MAX_SAT_NUM],ref[MAX_SAT_NUM],s[2];
        tSatInfoUnit *ref_sat_info= nullptr;
        tSatInfoUnit *sat_info= nullptr;

        for(i=0;i<MAX_SAT_NUM;i++) for(j=0;j<MAX_GNSS_USED_FRQ_NUM;j++){
            previous_sat_info_[i].fix[j]=0;
        }

        for(i=0;i<na;i++) D[i+i*num_full_x_]=1.0;

#if 0
        for(m=0;m<NSYS;m++){

            no_fix=(m==0&&gps==0)||(m==1&&ppk_conf_.gnssC.bds_ar_mode==0)||(m==3&&glo==0);

            for(f=0,k=na;f<nf;f++,k+=MAX_SAT_NUM){

                for(i=k;i<k+MAX_SAT_NUM;i++){
                    sat_info=&previous_sat_info_[i-k];
                    if(full_x_[i]==0.0||sat_info->sat.sat_.sys_idx!=m||!sat_info->vsat[f]){
                        continue;
                    }
                    if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&sat_info->el_az[0]*R2D>ppk_conf_.gnssC.ar_el_mask&&!no_fix){
                        sat_info->fix[f]=2;
                        break;
                    }
                    else sat_info->fix[f]=1;
                }

                if(sat_info->fix[f]!=2) continue;

                for(j=k;j<k+MAX_SAT_NUM;j++){
                    sat_info=&previous_sat_info_[j-k];
                    if(i==j||full_x_[j]==0.0||sat_info->sat.sat_.sys_idx!=m||!sat_info->vsat[f]) continue;
                    if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&previous_sat_info_[i-k].vsat[f]&&sat_info->el_az[0]*R2D>=ppk_conf_.gnssC.ar_el_mask&&!no_fix){
                        D[i+(na+nb)*num_full_x_]=1.0;
                        D[j+(na+nb)*num_full_x_]=-1.0;
                        ref[nb]=i-k+1;
                        fix[nb++]=j-k+1;
                        sat_info->fix[f]=2;
                    }
                    else sat_info->fix[f]=1;
                }

            }
        }
#endif
        int sat1,sat2;
        if(ddambs_.size()) ddambs_.clear();
        tDdAmb dd_amb={0};
        for(i=0;i<vflag_.size();i++){
            if(((vflag_[i]>>4)&0xF)==1) continue;
            sat1=(vflag_[i]>>16)&0xFF;
            sat2=(vflag_[i]>> 8)&0xFF;
            f=vflag_[i]&0xF;

            no_fix=(m==0&&gps==0)||(m==1&&ppk_conf_.gnssC.bds_ar_mode==0)||(m==3&&glo==0);
            ref_sat_info=&previous_sat_info_[sat1-1];
            sat_info=&previous_sat_info_[sat2-1];
            k=para_.IndexAmb(f,sat1);
            j=para_.IndexAmb(f,sat2);

            ref_sat_info->fix[f]=2;
            if(full_x_[j]==0.0||full_x_[k]==0.0) continue;

            if(sat_info->slip[f]) continue;

            if(sat_info->lock[f]>=0&&!(sat_info->slip[f]&2)&&sat_info->vsat[f]&&sat_info->el_az[0]*R2D>=ppk_conf_.gnssC.ar_el_mask&&!no_fix){
                D[k+(na+nb)*num_full_x_]= 1.0;
                D[j+(na+nb)*num_full_x_]=-1.0;

                dd_amb.ref_sat=ref_sat_info->sat.sat_.no;
                dd_amb.sat=sat_info->sat.sat_.no;
                dd_amb.f=f;
                ddambs_.push_back(dd_amb);
                ref[nb]=ref_sat_info->sat.sat_.no;
                fix[nb++]=sat_info->sat.sat_.no;
                sat_info->fix[f]=2;
            }
            else sat_info->fix[f]=1;
        }

        if(nb>0){
            VectorXd ref_sats,fix_sats;
            ref_sats=Map<VectorXd>(ref,nb,1);
            fix_sats=Map<VectorXd>(fix,nb,1);

            LOG(DEBUG)<<"REF SATELLITES: "<<ref_sats.transpose();
            LOG(DEBUG)<<"FIX SATELLITES: "<<fix_sats.transpose();
        }
        return nb;
    }

    int cPpkSolver::ResolveAmbLambda(double *xa, int gps, int glo) {
        int i,j,ny,nb,info,na=num_real_x_fix_;
        VectorXd QQb(MAX_SAT_NUM);
        double var=0.0;
        VectorXd s(2);
        string buff;

        ppplib_sol_.ratio=0.0;
        if(ppk_conf_.gnssC.ar_mode==AR_OFF||ppk_conf_.gnssC.ar_thres[0]<1.0){
            ppplib_sol_.num_ar_sat=0;
            return 0;
        }

        // skip AR if position variance too lager to avoid false fix
#if 1
        for(i=0;i<3;i++) var+=full_Px_(i,i);
        var=var/3.0;
        if(var>ppk_conf_.gnssC.ar_thres[1]){ //0.25
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"POSITIONING VARIANCE TOO LARGE FOR PPK AR var="<<var;
            ppplib_sol_.num_ar_sat=0;
            return 0;
        }
#endif
        // create single to DD transformation matrix,used to translate SD phase biases to DD
        MatrixXd D(num_full_x_,num_full_x_);
        D=Eigen::MatrixXd::Zero(num_full_x_,num_full_x_);

        if((nb=DdMat(D.data(),gps,glo))<(ppk_conf_.gnssC.min_sat_num2fix-1)){
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"NO ENOUGH VALID SATELLITE MAKING DOUBLE DIFFERENCE MATRIX";
            return -1;
        }

        ppplib_sol_.num_ar_sat=nb;

        ny=na+nb;
        VectorXd y(ny);
        MatrixXd Qy(ny,ny),DP(ny,num_full_x_);
        VectorXd db(nb,1);
        MatrixXd b(nb,2),Qb(nb,nb),Qab(na,nb),QQ(na,nb);

        y=D.transpose()*full_x_;
        MatMul("TN",ny,num_full_x_,num_full_x_,1.0,D.data() ,full_Px_.data(),0.0,DP.data());   /* DP=D'*P */
        MatMul("NN",ny,ny,num_full_x_,1.0,DP.data(),D.data(),0.0,Qy.data());

        //phase-bias covariance(Qb) and real-parameters to bias covariance(Qab)
        for(i=0;i<nb;i++){
            QQb[i]=Qy.data()[(na+i)+(na+i)*ny];
            for(j=0;j<nb;j++){
                Qb.data()[i+j*nb]=Qy.data()[(na+i)+(na+j)*ny];
            }
        }
        for(i=0;i<na;i++) for(j=0;j<nb;j++) Qab.data()[i+j*na]=Qy.data()[i+(na+j)*ny];

        VectorXd float_amb(nb),fix_float_diff(nb);
        float_amb=Map<VectorXd>(y.data()+na,nb,1); //DD ambiguity

        char ss[20];
        for(i=0;i<float_amb.size();i++){
            sprintf(ss,"%5.2f  ",float_amb[i]);
            buff+=ss;
        }
        LOG(DEBUG)<<"DD_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        for(i=0;i<float_amb.size();i++){
            sprintf(ss,"%7.5f  ",QQb[i]);
            buff+=ss;
        }
        LOG(DEBUG)<<" Q_AMB(0): "<<buff;
        buff.clear();
        ss[0]='\0';

        if(!(info=lambda_.IntegerAmb(float_amb,Qb,b,nb,2,s))){
            ppplib_sol_.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
            if(ppplib_sol_.ratio>999.9) ppplib_sol_.ratio=999.9f;

            for(i=0;i<nb;i++){
                sprintf(ss,"%5.2f  ",b.col(0)[i]);
                buff+=ss;
            }
            LOG(DEBUG)<<"DD_AMB(1): "<<buff;
            buff.clear();
            ss[0]='\0';

            for(i=0;i<nb;i++){
                sprintf(ss,"%5.2f  ",b.col(1)[i]);
                buff+=ss;
            }
            LOG(DEBUG)<<"DD_AMB(2): "<<buff;
            buff.clear();
            ss[0]='\0';

            // validation by popular ratio-test of residuals
            if(s[0]<=0.0||s[1]/s[0]>=ppk_conf_.gnssC.ar_thres[0]){
                // init non phase-bias states and covariance with float solution values
                for(i=0;i<na;i++){
                    real_x_fix_[i]=full_x_[i];
                    for(j=0;j<na;j++) real_Px_fix_.data()[i+j*na]=full_Px_.data()[i+j*num_full_x_];
                }

                for(i=0;i<nb;i++){
                    fix_float_diff[i]=y.data()[na+i]-b.col(0).data()[i];
                }

                MatrixXd Qb_=Qb.inverse();
                db=Qb_*fix_float_diff; //db=Qb^-1*(b0-b)
                real_x_fix_=real_x_fix_-Qab*db;   //xa=x-Qab-db
                QQ=Qab*Qb_;            // QQ=Qab*Qb^-1
                real_Px_fix_=real_Px_fix_-QQ*Qab.transpose();  //Pa=P-QQ*Qab^T

                ReSetAmb(b.col(0).data(),xa,nb);
                LOG(DEBUG)<<" AR VALIDATION OK, ratio="<<s[1]/s[0];
            }
            else{
                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<" RATIO VALIDATION FAILED ratio="<<s[1]/s[0];
                nb=0;
            }
        }
        else{
            LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<"PPK LAMBDA ERROR";
            nb=0;
        }

        return nb;
    }

    // function from rtklib
    bool cPpkSolver::ResolvePpkAmb(vector<int> cmn_sat, int nf, double *xa) {
        int nb=0,f,i,ar=0,lockc[MAX_GNSS_USED_FRQ_NUM],arsats[64]={0};
        bool exc_flag=false;
        tSatInfoUnit *sat_info= nullptr;

        LOG(DEBUG)<<"PREVIOUS EPOCH AR: ratio1="<<pre_epoch_ar_ratio1<<" ratio2="<<pre_epoch_ar_ratio2<<" ar_sat_num="<<ppplib_sol_.num_ar_sat;

        // 若前一历元没有ＡＲ成功
        if(pre_epoch_ar_ratio2<ppk_conf_.gnssC.ar_thres[0]&&ppplib_sol_.num_ar_sat>=ppk_conf_.gnssC.min_sat_num2drop){
            // previous epoch ar failed and satellite enough to drop
            for(f=0;f<nf;f++){
                for(i=0;i<cmn_sat.size();i++){
                    sat_info=&previous_sat_info_[cmn_sat[i]-1];
                    if(sat_info->vsat[f]&&sat_info->lock[f]>=0&&sat_info->el_az[0]*R2D>=ppk_conf_.gnssC.ele_min){
                        arsats[ar++]=i;
                    }
                }
            }

            if(exc_sat_index<ar){
                i=cmn_sat[arsats[exc_sat_index]];
                for(f=0;f<nf;f++){
                    lockc[f]=previous_sat_info_[i-1].lock[f];
                    previous_sat_info_[i-1].lock[f]=-ppplib_sol_.num_ar_sat;
                }
                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<previous_sat_info_[i-1].sat.sat_.id<<" EXCLUDE BY AR";
                exc_flag=true;
            }
            else exc_sat_index=0;
        }

        double ratio1;
        int dly;
        bool rerun=false;
        if(ppk_conf_.gnssC.glo_ar_mode!=GLO_AR_FIXHOLD){

            nb=ResolveAmbLambda(xa,1,0);

            ratio1=ppplib_sol_.ratio;
            if(ppk_conf_.gnssC.partial_ar){
                if(nb>=0&&pre_epoch_ar_ratio2>=ppk_conf_.gnssC.ar_thres[0]&&(ppplib_sol_.ratio<ppk_conf_.gnssC.ar_thres[0])||
                        (ppplib_sol_.ratio<ppk_conf_.gnssC.ar_thres[0]*1.1&&ppplib_sol_.ratio<pre_epoch_ar_ratio1/2.0)){
                    LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1);
                    dly=2;
                    for(i=0;i<cmn_sat.size();i++){
                        for(f=0;f<nf;f++){
                            if(previous_sat_info_[cmn_sat[i]-1].fix[f]!=2) continue;
                            if(previous_sat_info_[cmn_sat[i]-1].lock[f]==0){
                                LOG(WARNING)<<epoch_sat_info_collect_[0].t_tag.GetTimeStr(1)<<" "<<previous_sat_info_[cmn_sat[i]-1].sat.sat_.id<<" L"<<f+1<<" IS NEW OBSERVED";
                                previous_sat_info_[cmn_sat[i]-1].lock[f]=-ppk_conf_.gnssC.min_lock2fix-dly;
                                dly+=2;
                                rerun=true;
                            }
                        }
                    }
                }

                if(rerun){
                    if((nb=ResolveAmbLambda(xa,1,0))<2){
                        pre_epoch_ar_ratio1=pre_epoch_ar_ratio2=ratio1;
                        return false;
                    }
                }
            }
//            pre_epoch_ar_ratio1=ratio1;
        }
        else{
            ratio1=0.0;
            nb=0;
        }

        if(exc_flag&&(ppplib_sol_.ratio<ppk_conf_.gnssC.ar_thres[0])&&ppplib_sol_.ratio<1.5*pre_epoch_ar_ratio2){
            i=cmn_sat[arsats[exc_sat_index++]];
            for(f=0;f<nf;f++) previous_sat_info_[i-1].lock[f]=lockc[f];
        }

        pre_epoch_ar_ratio1=ratio1>0?ratio1:ppplib_sol_.ratio;
        pre_epoch_ar_ratio2=ppplib_sol_.ratio;

        if(nb>0) return true;
        else return false;
    }

    cFusionSolver::cFusionSolver() {}

    cFusionSolver::cFusionSolver(tPPPLibConf C) {
        fs_conf_=C;
        num_full_x_=para_.GetPPPLibPar(C);
        full_x_=VectorXd::Zero(num_full_x_);
    }

    cFusionSolver::~cFusionSolver() {}

    bool cFusionSolver::InputImuData(int ws) {
        int i,week;

        if(imu_index_++<0||imu_index_>=imu_data_.data_.size()) return false;

        // prepare imu data for static detect
//        for(i=imu_index_;i>=0&&i<ws&&i<imu_data_.data_.size();i++){
//            imu_data_zd_.push_back(imu_data_.data_.at(i));
//        }

        cur_imu_info_.t_tag=imu_data_.data_[imu_index_].t_tag;

        if(fs_conf_.insC.imu_type==IMU_M39){
            gnss_sols_[0].t_tag.Time2Gpst(&week, nullptr,SYS_GPS);
            cur_imu_info_.t_tag+=week*604800.0;
        }

        cur_imu_info_.raw_gyro=imu_data_.data_[imu_index_].gyro;
        cur_imu_info_.raw_acce=imu_data_.data_[imu_index_].acce;
        return true;
    }

    bool cFusionSolver::MatchGnssObs() {
        double sow1,sow2;
        int i,week=0,wod=0,info=false;

        for(i=rover_idx_-100<0?0:rover_idx_-10;rover_idx_<rover_obs_.GetGnssObs().size();i++){
            sow1=rover_obs_.GetGnssObs()[i].obs_time.Time2Gpst(&week,&wod,SYS_GPS);
            sow2=cur_imu_info_.t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<DTTOL){
                rover_idx_=i;info=true;

                break;
            }
            else if((sow1-sow2)>2.0*DTTOL){
                info=false;break;
            }
        }

        if(info&&fs_conf_.mode_opt==MODE_OPT_PPK){
            for(i=base_idx_-100<0?0:base_idx_-10;base_idx_<gnss_solver_->base_obs_.GetGnssObs().size();i++){
                sow1=gnss_solver_->base_obs_.GetGnssObs().at(i).obs_time.Time2Gpst(&week,&wod,SYS_GPS);
                sow2=rover_obs_.GetGnssObs()[rover_idx_].obs_time.Time2Gpst(nullptr, nullptr,SYS_GPS);
                if(fabs(sow1-sow2)<DTTOL){
                    base_idx_=i;info=true;
                    break;
                }
                else if((sow1-sow2)>2.0*DTTOL){
                    info=false;break;
                }
            }
        }

        return info;
    }

    bool cFusionSolver::MatchGnssSol() {
        double sow1,sow2;
        bool info=false;
        int i;

        for(i=gnss_sol_idx>3?gnss_sol_idx-3:0;i<gnss_sols_.size()&&i>=0;i++){
            sow1=gnss_sols_[i].t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            sow2=cur_imu_info_.t_tag.Time2Gpst(nullptr, nullptr,SYS_GPS);
            if(fabs(sow1-sow2)<0.5/fs_conf_.insC.sample_rate){
                info=true;gnss_sol_idx=i;
                break;
            }
            else if((sow1-sow2)>1.0/fs_conf_.insC.sample_rate){
                info=false;break;
            }
        }
        return info;
    }

    void cFusionSolver::InsSol2PpplibSol(tImuInfoUnit &imu_sol, tSolInfoUnit &ppplib_sol) {
        ppplib_sol_.t_tag=imu_sol.t_tag;
        ppplib_sol.pos=imu_sol.re;
        ppplib_sol.vel=imu_sol.ve;
        ppplib_sol.att=imu_sol.rpy;

        ppplib_sol.accl_bias=imu_sol.ba;
        ppplib_sol.gyro_bias=imu_sol.bg;

        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        Matrix3d Pp,Pv,Pa;
        Pp=full_Px_.block<3,3>(ip,ip);
        Pv=full_Px_.block<3,3>(iv,iv);
        Pa=full_Px_.block<3,3>(ia,ia);

        for(int i=0;i<3;i++) ppplib_sol.q_pos[i]=Pp(i,i);
        ppplib_sol.q_pos[3]=Pp(0,1);
        ppplib_sol.q_pos[4]=Pp(1,2);
        ppplib_sol.q_pos[5]=Pp(2,0);

        for(int i=0;i<3;i++) ppplib_sol.q_vel[i]=Pv(i,i);
        ppplib_sol.q_vel[3]=Pv(0,1);
        ppplib_sol.q_vel[4]=Pv(1,2);
        ppplib_sol.q_vel[5]=Pv(2,0);

        for(int i=0;i<3;i++) ppplib_sol.q_att[i]=Pa(i,i);
        ppplib_sol.q_att[3]=Pa(0,1);
        ppplib_sol.q_att[4]=Pa(1,2);
        ppplib_sol.q_att[5]=Pa(2,0);
    }

    double cFusionSolver::Vel2Yaw(Vector3d vn) {
        return atan2(vn[1],fabs(vn[0])<1E-4?1E-4:vn[0]);
    }

    Vector3d cFusionSolver::Pos2Vel(tSolInfoUnit &sol1, tSolInfoUnit &sol2) {
        Vector3d vel={0,0,0};
        double t=sol1.t_tag.TimeDiff(sol2.t_tag.t_);

        vel=(sol1.pos-sol2.pos)/t;
        return vel;
    }

    bool cFusionSolver::GnssSol2Ins(Vector3d re,Vector3d ve) {
        if(re.norm()==0.0||ve.norm()==0.0) return false;
        Vector3d blh=Xyz2Blh(re),rn;
        Vector3d wiee(0,0,OMGE_GPS);
        Matrix3d Cne;
        Cne=CalcCen(blh,COORD_NED).transpose();
        pre_imu_info_.rn=blh;
        pre_imu_info_.vn=Cne.transpose()*ve;

        pre_imu_info_.rpy[2]=Vel2Yaw(pre_imu_info_.vn);
        pre_imu_info_.Cbn=Euler2RotationMatrix(pre_imu_info_.rpy);
        pre_imu_info_.Cbe=Cne*pre_imu_info_.Cbn;

        // lever correct
        pre_imu_info_.re=re-pre_imu_info_.Cbe*fs_conf_.insC.lever;
        Matrix3d T=VectorSkew(cur_imu_info_.raw_gyro*cur_imu_info_.dt);
        Matrix3d Omge=VectorSkew(wiee);
        pre_imu_info_.ve=ve-pre_imu_info_.Cbe*T*fs_conf_.insC.lever+Omge*pre_imu_info_.Cbe*fs_conf_.insC.lever;
        pre_imu_info_.rn=Xyz2Blh(pre_imu_info_.re);
        pre_imu_info_.vn=Cne.transpose()*pre_imu_info_.ve;

        return true;
    }

    bool cFusionSolver::InsAlign(int use_raw_gnss_obs) {
        int num_sols=5;
        static vector<tSolInfoUnit>sols;

        if(use_raw_gnss_obs){

            if(gnss_solver_->SolverProcess(gnss_conf_,rover_idx_)) {
                LOG(INFO)<<"USING GNSS TO ALIGN INS("<<sols.size()<<"): "<<gnss_solver_->ppplib_sol_.pos.transpose();
                if(gnss_solver_->ppplib_sol_.stat==SOL_FIX){
                    sols.push_back(gnss_solver_->ppplib_sol_);
                }
            }
            else{
                sols.clear();
                return false;
            }

        }
        else{
            sols.push_back(gnss_sols_.at(gnss_sol_idx));
        }

        if(sols.size()>=num_sols){
            if(use_raw_gnss_obs){
                for(int i=0;i<sols.size()-1;i++){
                    if(sols[i+1].t_tag.TimeDiff(sols[i].t_tag.t_)>fs_conf_.gnssC.sample_rate){
                        sols.clear();
                        return false;
                    }
                }
            }
            else{
                for(int i=0;i<sols.size();i++){
                    if(sols[i].stat>SOL_FIX||sols[i].stat==SOL_NONE){
                        sols.clear();
                        return false;
                    }
                }
            }

            Vector3d ve;
            if(sols.back().vel.norm()==0.0){
                ve=Pos2Vel(sols.back(),*(sols.end()-2));
            }
            else ve=sols.back().vel;

            if(ve.norm()<5.0||cur_imu_info_.raw_gyro.norm()>30.0*D2R){
                sols.clear();
                return false;
            }


            if(GnssSol2Ins(sols.back().pos,ve)){
                LOG(INFO)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" INS INITIALIZATION OK";
                LOG(INFO)<<"INIT POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pre_imu_info_.re.transpose()<<" m";
                LOG(INFO)<<"INIT VELOCITY(n): "<<setw(13)<<setprecision(3)<<pre_imu_info_.vn.transpose()<<" m/s";
                LOG(INFO)<<"INIT ATTITUTE(n): "<<setw(13)<<setprecision(3)<<pre_imu_info_.rpy.transpose()*R2D<<" deg";
                ppplib_sol_.ins_stat=SOL_INS_INIT;
                InsSol2PpplibSol(pre_imu_info_,ppplib_sol_);
                return true;
            }
            else sols.clear();
        }
        return false;
    }

    void cFusionSolver::InitSolver(tPPPLibConf C) {
        fs_conf_=C;
        para_=cParSetting(C);
        gnss_err_corr_.InitGnssErrCorr(C,&nav_);
        out_=new cOutSol(C);
        out_->InitOutSol(C,C.fileC.sol);
        out_->WriteHead();
        InitFullPx(C);

        if(fs_conf_.mode_opt==MODE_OPT_SPP){
            tPPPLibConf spp_conf=C;
            if(C.mode==MODE_IGLC) spp_conf.mode=MODE_SPP;
            spp_conf.solC.out_sol=false;
            gnss_solver_=new cSppSolver(spp_conf);
            gnss_solver_->nav_=nav_;
            gnss_solver_->InitSolver(spp_conf);
            gnss_conf_=spp_conf;
        }
        else if(fs_conf_.mode_opt==MODE_OPT_PPP){
            tPPPLibConf ppp_conf=C;
            if(C.mode==MODE_IGLC) ppp_conf.mode=MODE_PPP;
            ppp_conf.gnssC.frq_opt=FRQ_DUAL;
            ppp_conf.solC.out_sol=false;
            gnss_solver_=new cPppSolver(ppp_conf);
            gnss_solver_->nav_=nav_;
            gnss_solver_->InitSolver(C);

            gnss_conf_=ppp_conf;
        }
        else if(fs_conf_.mode_opt==MODE_OPT_PPK){
            tPPPLibConf ppk_conf=C;
            if(C.mode==MODE_IGLC) ppk_conf.mode=MODE_PPK;
            ppk_conf.solC.out_sol=false;
            ppk_conf.solC.out_head=false;
            gnss_solver_=new cPpkSolver(ppk_conf);
            gnss_solver_->nav_=nav_;
            gnss_solver_->InitSolver(ppk_conf);
            gnss_conf_=ppk_conf;
        }
        else if(fs_conf_.mode_opt==MODE_OPT_GSOF){
            cDecodeGsof gsof_decoder;
            gsof_decoder.DecodeGsof(C.fileC.gsof,gnss_sols_);
        }

        if(C.insC.imu_type==IMU_M39){
            cDecodeImuM39 m39_decoder;
            m39_decoder.DecodeM39(C.fileC.imu,C.insC,imu_data_.data_);
        }
        else if(C.insC.imu_type==IMU_NOVTEL_A1||C.insC.imu_type==IMU_NOVTEL_CPT){
            cReadImu imu_reader(C.fileC.imu);
            imu_reader.SetImu(C);
            imu_reader.Reading();
            imu_data_=*imu_reader.GetImus();
        }

        gnss_solver_->rover_obs_=rover_obs_;
    }

    bool cFusionSolver::SolverProcess(tPPPLibConf C,int idx) {

        InitSolver(C);
        tSatInfoUnit sat_info;

        int gnss_obs_flag=false,ins_align=false,gnss_sol_flag=false;
        while(InputImuData(5)){
            if(fs_conf_.mode>MODE_INS){
                if(fs_conf_.mode_opt==MODE_OPT_GSOF){
                    gnss_sol_flag=MatchGnssSol();
                }
                else{
                    gnss_obs_flag=MatchGnssObs();
                }

                if(gnss_sol_flag||gnss_obs_flag){

                    if(!ins_align){
                        if(gnss_obs_flag) ins_align=InsAlign(1);
                        else {
                            if(gnss_sols_[gnss_sol_idx].pos.norm()==0.0) continue;

                            // initial alignment
                            ins_align=InsAlign(0);
                        }
                        pre_imu_info_.t_tag=cur_imu_info_.t_tag;
                        pre_imu_info_.cor_gyro=pre_imu_info_.raw_gyro=cur_imu_info_.raw_gyro;
                        pre_imu_info_.cor_acce=pre_imu_info_.raw_acce=cur_imu_info_.raw_acce;
                        epoch_sat_obs_.epoch_data.clear();
                        base_epoch_sat_obs_.epoch_data.clear();
                        continue;
                    }

                    // measurement update
                    if(C.mode==MODE_IGLC){
                        LooseCouple(C);
                    }
                    else if(C.mode==MODE_IGTC) TightCouple(C);
                }
                else{
                    if(!ins_align){
                        pre_imu_info_=cur_imu_info_;
                        continue;
                    }
                    // time update
                    StateTimeUpdate();
                    ppplib_sol_.stat=SOL_NONE;
                    if(ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
                        ppplib_sol_.ins_stat=SOL_INS_MECH;
                        InsSol2PpplibSol(cur_imu_info_,ppplib_sol_);
                    }
                }
            }
            else{
                // pure ins mech
                ppplib_sol_.stat=SOL_NONE;
                if(ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
                    ppplib_sol_.ins_stat=SOL_INS_MECH;
                    InsSol2PpplibSol(cur_imu_info_,ppplib_sol_);
                }
            }

            pre_imu_info_=cur_imu_info_;
            if(C.solC.out_ins_mech_frq!=0&&ins_mech_idx%C.solC.out_ins_mech_frq==0&&ppplib_sol_.ins_stat==SOL_INS_MECH){
                ppplib_sol_.valid_sat_num=0;
                out_->WriteSol(ppplib_sol_,epoch_idx_);
            }
            else if(ppplib_sol_.ins_stat!=SOL_INS_MECH) out_->WriteSol(ppplib_sol_,epoch_idx_);
        }
    }

    bool cFusionSolver::SolverEpoch() {
        gnss_solver_->epoch_idx_=epoch_idx_;
        gnss_solver_->epoch_sat_obs_=epoch_sat_obs_;
        if(gnss_solver_->SolverEpoch()){
            gnss_solver_->out_->WriteSol(gnss_solver_->ppplib_sol_,epoch_idx_);
        }
    }

    bool cFusionSolver::SolutionUpdate() {

    }

    void cFusionSolver::CloseLoopState(VectorXd& x,tImuInfoUnit* imu_info_corr) {
        int ip=para_.IndexPos();
        imu_info_corr->re[0]-=x[ip+0];
        imu_info_corr->re[1]-=x[ip+1];
        imu_info_corr->re[2]-=x[ip+2];
        imu_info_corr->rn=Xyz2Blh(imu_info_corr->re);

        int iv=para_.IndexVel();
        imu_info_corr->ve[0]-=x[iv+0];
        imu_info_corr->ve[1]-=x[iv+1];
        imu_info_corr->ve[2]-=x[iv+2];
        Matrix3d Cne=CalcCen(imu_info_corr->rn,COORD_NED).transpose();
        imu_info_corr->vn=Cne.transpose()*imu_info_corr->ve;

        int ia=para_.IndexAtt();
        if(x[ia]!=DIS_FLAG){
            Vector3d att(x[ia],x[ia+1],x[ia+2]);
            Matrix3d T=Matrix3d::Identity()-VectorSkew(att);
            imu_info_corr->Cbe=T*imu_info_corr->Cbe;
        }

        int iba=para_.IndexBa();
        if(x[iba]!=DIS_FLAG){
            imu_info_corr->ba[0]+=x[iba+0];
            imu_info_corr->ba[1]+=x[iba+1];
            imu_info_corr->ba[2]+=x[iba+2];
        }

        int ibg=para_.IndexBg();
        if(x[ibg]!=DIS_FLAG){
            imu_info_corr->bg[0]+=x[ibg+0];
            imu_info_corr->bg[1]+=x[ibg+1];
            imu_info_corr->bg[2]+=x[ibg+2];
        }
    }

    void cFusionSolver::DisableX(VectorXd &x, int idx) {
        if(idx>para_.NumClPar()) return;
        else{
            x[idx]=DIS_FLAG;
        }
    }

    void cFusionSolver::StateTimeUpdate() {
        double dt=cur_imu_info_.dt;
        int nx=para_.GetInsTransParNum(fs_conf_);
        MatrixXd F,P,Q;
        F=MatrixXd::Zero(nx,nx);

        F=ins_mech_.StateTransferMat(fs_conf_,pre_imu_info_,cur_imu_info_,nx,dt);
        Q=InitQ(fs_conf_,dt);

        if(fabs(dt)>60.0){

        }
        else{
            PropVariance(F,Q,nx);
        }

        for(int i=0;i<nx;i++) full_x_[i]=1E-20;
    }

    void cFusionSolver::PropVariance(MatrixXd& F,MatrixXd& Q,int nx) {

        MatrixXd prior_P=full_Px_.block<15,15>(0,0);
        MatrixXd PQ(nx,nx),FPF(nx,nx);

        PQ=MatrixXd::Zero(nx,nx);FPF=MatrixXd::Zero(nx,nx);

        int i,j;
        for(i=0;i<nx;i++){
            for(j=0;j<nx;j++){
                PQ.data()[i+j*nx]=prior_P.data()[i+j*nx]+0.5*Q.data()[i+j*nx];
            }
        }

        FPF=F*PQ*F.transpose();

#if 0
        cout<<std::fixed<<setprecision(12)<<prior_P<<endl<<endl;

        cout<<std::fixed<<setprecision(12)<<Q<<endl<<endl;

        cout<<std::fixed<<setprecision(10)<<PQ<<endl<<endl;

        cout<<std::fixed<<setprecision(10)<<F<<endl<<endl;

        cout<<std::fixed<<setprecision(10)<<FPF<<endl<<endl;
#endif

        for(i=0;i<nx;i++){
            for(j=0;j<nx;j++){
                full_Px_.data()[i+j*nx]=FPF.data()[i+j*nx]+0.5*Q.data()[i+j*nx];
            }
        }
    }

    void cFusionSolver::RemoveLever(const tImuInfoUnit &imu_info, Vector3d &lever, Vector3d &gnss_re,
                                    Vector3d &gnss_ve) {
        Matrix3d Cbe=imu_info.Cbe;
        Vector3d wiee(0.0,0.0,OMGE_GPS);
        Vector3d T=Cbe*lever;

        // position correction
        gnss_re=imu_info.re+T;

        //velocity correction
        Vector3d omge=imu_info.cor_gyro,W;
        T=Cbe*VectorSkew(omge)*lever;
        W=VectorSkew(wiee)*Cbe*lever;
        gnss_ve=imu_info.ve+T-W;
    }

    int cFusionSolver::BuildLcHVR(int post,tPPPLibConf C,tImuInfoUnit& imu_info,double *meas_pos,double *meas_vel,Vector3d& q_pos,Vector3d& q_vel) {
        Vector3d gnss_re,gnss_ve;
        Matrix3d Cbe=imu_info.Cbe;
        double omc=0.0;
        vector<double>omcs;
        num_L_=0;

        if(!meas_pos&&meas_vel) return 0;
        if(meas_pos) num_L_+=3;
        if(meas_vel) num_L_+=3;
        H_=MatrixXd::Zero(num_L_,num_full_x_);
        R_=MatrixXd::Zero(num_L_,num_L_);

        RemoveLever(imu_info,C.insC.lever,gnss_re,gnss_ve);

        int ip=para_.IndexPos();
        int iv=para_.IndexVel();
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        if(meas_pos){
            for(int i=0;i<3;i++){
                omc=meas_pos[i]-gnss_re[i];
                omcs.push_back(omc);
            }

            if(!post){
                // position to position jacobian
                H_.block<3,3>(ip,ip)=-1.0*Matrix3d::Identity();

                // position to attitude jacobian
                H_.block<3,3>(ip,ia)=VectorSkew(Cbe*C.insC.lever);
            }
            R_.block<3,3>(ip,ip)=q_pos.asDiagonal();
        }

        if(meas_vel){
            for(int i=0;i<3;i++){
                omc=meas_vel[i]-gnss_ve[i];
                omcs.push_back(omc);
            }

            if(!post){
                // velocity to velocity jacobian
                H_.block<3,3>(iv,iv)=-1.0*Matrix3d::Identity();

                // velocity to attitude jacobian
                Vector3d wiee(0.0,0.0,OMGE_GPS);
                Vector3d T,W;
                T=Cbe*VectorSkew(imu_info.cor_gyro)*C.insC.lever;
                W=VectorSkew(wiee)*Cbe*C.insC.lever;
                H_.block<3,3>(iv,ia)=VectorSkew(T-W);

                // velocity to ba jacobian
                H_.block<3,3>(iv,iba)=Cbe*(VectorSkew(C.insC.lever));

            }
            R_.block<3,3>(3,3)=q_vel.asDiagonal();
        }

        omc_L_=Map<VectorXd>(omcs.data(),num_L_);
        omcs.clear();

        return num_L_;
    }

    bool cFusionSolver::LcFilter(tPPPLibConf C) {
        double *meas_pos= nullptr,*meas_vel= nullptr;
        Vector3d q_pos,q_vel;

        if(C.mode_opt>MODE_OPT_GSOF){
            ppplib_sol_.stat=gnss_solver_->ppplib_sol_.stat;
            ppplib_sol_.pos=gnss_solver_->ppplib_sol_.pos;
            ppplib_sol_.vel=gnss_solver_->ppplib_sol_.vel;
            if(ppplib_sol_.pos.norm()!=0.0){
                meas_pos=ppplib_sol_.pos.data();
                q_pos<<ppplib_sol_.q_pos[0],ppplib_sol_.q_pos[1],ppplib_sol_.q_pos[2];
            }
            if(ppplib_sol_.vel.norm()!=0.0){
                meas_vel=ppplib_sol_.vel.data();
                q_vel<<ppplib_sol_.q_vel[0],ppplib_sol_.q_vel[1],ppplib_sol_.q_vel[2];
                if(q_vel.norm()==0.0){
                    q_vel<<0.1,0.1,0.1;
                }
            }
        }
        else{
            ppplib_sol_.stat=gnss_sols_[gnss_sol_idx].stat;
            ppplib_sol_.valid_sat_num=gnss_sols_[gnss_sol_idx].valid_sat_num;
            if(gnss_sols_[gnss_sol_idx].pos.norm()!=0.0){
                meas_pos=gnss_sols_[gnss_sol_idx].pos.data();
                q_pos<<(gnss_sols_[gnss_sol_idx].q_pos[0]),(gnss_sols_[gnss_sol_idx].q_pos[1]),(gnss_sols_[gnss_sol_idx].q_pos[2]);
            }
            if(gnss_sols_[gnss_sol_idx].vel.norm()!=0.0){
                meas_vel=gnss_sols_[gnss_sol_idx].vel.data();
                q_vel<<(gnss_sols_[gnss_sol_idx].q_vel[0]),(gnss_sols_[gnss_sol_idx].q_vel[1]),(gnss_sols_[gnss_sol_idx].q_vel[2]);
                if(q_vel.norm()==0.0){
                    q_vel<<0.1,0.1,0.1;
                }
            }
        }

        Vector3d ins_re=cur_imu_info_.re,ins_ve=cur_imu_info_.ve;

        Vector3d pos,vel;
        for(int i=0;i<3;i++){
            if(meas_pos) pos[i]=meas_pos[i];
            if(meas_vel) vel[i]=meas_vel[i];
        }
        LOG(DEBUG)<<"INS GNSS   FUSION(-): "<<cur_imu_info_.t_tag.GetTimeStr(3);
        LOG(DEBUG)<<"INS MECH POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.re.transpose();
        LOG(DEBUG)<<"INS MECH VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.ve.transpose();
        LOG(DEBUG)<<"GNSS     POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pos.transpose()<<" "<<q_pos.transpose();
        LOG(DEBUG)<<"GNSS     VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<vel.transpose()<<" "<<q_vel.transpose();

        VectorXd x=full_x_;
        MatrixXd Px=full_Px_;
        tImuInfoUnit imu_corr=cur_imu_info_;
        bool stat=false;
        if(BuildLcHVR(0,C,cur_imu_info_,meas_pos,meas_vel,q_pos,q_vel)){

            if(R_.trace()>=10.0){
                for(int i=0;i<para_.IndexAtt()+para_.NumAtt();i++){
                    DisableX(x,i);
                }
                for(int i=0;i<para_.IndexBa()+para_.NumBa();i++){
                    DisableX(x,i);
                }
                for(int i=0;i<para_.IndexBg()+para_.NumBg();i++){
                    DisableX(x,i);
                }
            }

#if 0
            cout<<omc_L_.transpose()<<endl<<endl;
            cout<<R_<<endl<<endl;
            cout<<H_.transpose()<<endl;
            cout<<full_Px_<<endl;
#endif
            kf_.Adjustment(omc_L_,H_.transpose(),R_,x,Px,num_L_,num_full_x_);
            cout<<H_.transpose()<<endl;

            CloseLoopState(x,&imu_corr);

            if(BuildLcHVR(1,C,imu_corr,meas_pos,meas_vel,q_pos,q_vel)){

                stat=ValidSol(x, 10.0);
                if(stat){
                    full_x_=x;
                    full_Px_=Px;
                    cur_imu_info_=imu_corr;
                    LOG(DEBUG)<<"INS GNSS   FUSION(+): "<<cur_imu_info_.t_tag.GetTimeStr(3);
                    LOG(DEBUG)<<"INS MECH POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.re.transpose();
                    LOG(DEBUG)<<"INS MECH VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<cur_imu_info_.ve.transpose();
                    LOG(DEBUG)<<"GNSS     POSITION(e): "<<setw(13)<<std::fixed<<setprecision(3)<<pos.transpose()<<" "<<q_pos.transpose();
                    LOG(DEBUG)<<"GNSS     VELOCITY(e): "<<setw(13)<<std::fixed<<setprecision(3)<<vel.transpose()<<" "<<q_vel.transpose();
                }
            }
        }
        else{
            LOG(WARNING)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" MAKE OMC ERROR";
            stat=false;
        }

        return stat;
    }

    bool cFusionSolver::ValidSol(VectorXd& x, double thres) {
        int flag=0;
        double fact=SQR(thres);
        int ia=para_.IndexAtt();
        int iba=para_.IndexBa();
        int ibg=para_.IndexBg();
        Vector3d att(x[ia],x[ia+1],x[ia+2]);
        Vector3d ba(x[iba],x[iba+1],x[iba+2]);
        Vector3d bg(x[ibg],x[ibg+1],x[ibg+2]);
        if(x[ia]==DIS_FLAG&&att.norm()>5.0*D2R){
            flag|=1;
            LOG(WARNING)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" "<<"NON-ESTIMATE ATTITUDE AND OVERRUN";
        }
        if(x[iba]==DIS_FLAG&&ba.norm()>1E5*MG2MS2){
            LOG(WARNING)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" "<<"NON-ESTIMATE ACCE BIAS AND OVERRUN";
            flag|=1;
        }
        if(x[ibg]==DIS_FLAG&&bg.norm()>10.0*D2R){
            LOG(WARNING)<<cur_imu_info_.t_tag.GetTimeStr(3)<<" "<<"NON-ESTIMATE GYRO BIAS AND OVERRUN";
            flag|=1;
        }

        if(flag){
            return false;
        }

        int i,j;
        for(j=0,i=0;i<num_L_&&thres;i++){
            if(omc_L_[i]*omc_L_[i]<fact*R_.data()[i+i*num_L_]) continue;
            j++;
        }
        return j<num_L_;
    }

    bool cFusionSolver::LooseCouple(tPPPLibConf C) {

        StateTimeUpdate();
        ppplib_sol_.stat=SOL_NONE;
        ppplib_sol_.valid_sat_num=0;
        if(!ins_mech_.InsMechanization(C.insC.err_model,pre_imu_info_,cur_imu_info_,++ins_mech_idx)){
            return false;
        }

        ppplib_sol_.ins_stat=SOL_INS_MECH;
        epoch_idx_++;

        if(fs_conf_.mode_opt>MODE_OPT_GSOF){
            if(gnss_solver_->SolverProcess(gnss_conf_,rover_idx_)){
                ppplib_sol_=gnss_solver_->ppplib_sol_;
                if(LcFilter(C)){
                    ppplib_sol_.ins_stat=SOL_IG_LC;
                    InsSol2PpplibSol(cur_imu_info_,ppplib_sol_);
                }
            }
        }
        else{
            if(LcFilter(C)){
                ppplib_sol_.ins_stat=SOL_IG_LC;
                InsSol2PpplibSol(cur_imu_info_,ppplib_sol_);
            }
        }


    }

    bool cFusionSolver::TightCouple(tPPPLibConf C) {

    }
}