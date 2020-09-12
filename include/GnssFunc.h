//
// Created by cc on 7/16/20.
//

#ifndef PPPLIB_GNSSFUNC_H
#define PPPLIB_GNSSFUNC_H

#include "CmnFunc.h"

namespace PPPLib{

    typedef struct{
        int no;
        int sys;
        int prn;
        string id;
        int sys_idx;
    }tSat;

    class cSat{
    public:
        cSat();
        cSat(int sat_no);
        cSat(int sat_sys,int sat_prn);
        cSat(string sat_id);
        ~cSat();

    public:
        void SatPrn2No();
        void SatNo2Prn();
        void SatId2No();
        void SatNo2Id();

    public:
        tSat sat_={0};
    };

    typedef struct{
        cSat sat;
        int iode,iodc;
        int sva,svh;
        int week,code;
        cTime toe,toc,ttr;
        double A,e,i0,Omg0,omg,M0,deln,Omgd,idot;
        double crc,crs,cuc,cus,cic,cis;
        double toes;
        double f0,f1,f2;
        Vector4d tgd;
        double Adot,ndot;
    }tBrdEphUnit;

    typedef struct{
        cSat sat;
        int iode,frq;
        int svh,sva,age;
        cTime toe,tof;
        double taun,gamn,dtaun;
        Vector3d pos;
        Vector3d vel;
        Vector3d acc;
    }tBrdGloEphUnit;

    typedef struct{
        cTime t_tag;
        Vector4d pos[MAX_SAT_NUM];
        Vector4d vel[MAX_SAT_NUM];
        Vector4d std_pos[MAX_SAT_NUM];
        Vector4d std_vel[MAX_SAT_NUM];
        double clk[MAX_SAT_NUM];
        double std_clk[MAX_SAT_NUM];
    }tPreOrbUnit;

    typedef struct{
        cTime t_tag;
        double clk[MAX_SAT_NUM];
        float  std[MAX_SAT_NUM];
    }tPreClkUnit;

    typedef struct{
        double mjd;
        double xp,yp;
        double xpr,ypr;
        double ut1_utc;
        double lod;
    }tErpUnit;

    typedef struct{
        cSat  sat;
        cTime ts,te;
        string ant_type;
        string ser_code;
        Vector3d pco[MAX_GNSS_FRQ_NUM*NSYS];
        double   pcv[MAX_GNSS_FRQ_NUM*NSYS][80*30];
        double dazi;
        double zen1,zen2,dzen;
        Vector3d rec_ant_del[2];
    }tAntUnit;

    typedef struct{
        cTime t_tag;
        int ndata[3]{};
        double re;
        double lats[3]{};
        double lons[3]{};
        double hgts[3]{};
        vector<double> data;
        vector<float>  rms;
    }tTecUnit;

    typedef struct {
        string name;
        string marker;
        string ant_desc;
        string ant_seri;
        string rec_type;
        string firm_ver;
        Vector3d del;
        Vector3d apr_pos;
        double ant_hgt;
    }tStaInfoUnit;

    typedef struct{
        vector<tBrdEphUnit>    brd_eph;
        vector<tBrdGloEphUnit> brd_glo_eph;
        vector<tPreOrbUnit> pre_eph;
        vector<tPreClkUnit> pre_clk;
        vector<tErpUnit>erp_paras;
//        vector<tAntUnit>ant_paras;
        vector<tTecUnit>tec_paras;

        tStaInfoUnit sta_paras[2];
        tAntUnit sat_ant[MAX_SAT_NUM];
        tAntUnit rec_ant[2];

        int glo_frq_num[GLO_MAX_PRN+1];
        double glo_cp_bias[4];
        int leaps;
        double ion_para[NSYS][8];
        double utc_para[NSYS][4];
        double code_bias[MAX_SAT_NUM][MAX_GNSS_CODE_BIAS_PAIRS];
        double ocean_paras[2][6*11];

        double wide_line_bias[MAX_SAT_NUM];
    }tNav;

    typedef struct{
        int n;
        int frq[MAX_GNSS_OBS_TYPE];
        int pos[MAX_GNSS_OBS_TYPE];
        unsigned char pri[MAX_GNSS_OBS_TYPE];
        unsigned char type[MAX_GNSS_OBS_TYPE];
        unsigned char code[MAX_GNSS_OBS_TYPE];
        double shift[MAX_GNSS_OBS_TYPE];
    }tGnssSignal;

    typedef struct{
        cSat sat;
        double P[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        unsigned char code[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        double L[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        float D[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        unsigned char SNR[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        unsigned char LLI[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        double frq[MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS];
        double CSC_P[MAX_GNSS_FRQ_NUM];
    }tSatObsUnit;

    typedef struct{
        cTime obs_time;
        int sat_num;
        vector<tSatObsUnit> epoch_data;
    }tEpochSatUnit;

    class cGnssObs{
    public:
        cGnssObs();
        ~cGnssObs();

    public:
        void SetTimeSpan(cTime* ts, cTime* te);
        cTime* GetStartTime();
        cTime* GetEndTime();
        void SetRcvIdx(RECEIVER_INDEX rcv);
        vector<tEpochSatUnit>& GetGnssObs();
        tStaInfoUnit* GetStation();

        void CarrierSmoothCode(int smooth_length);

    public:
        int epoch_num;
        double sample_;
        RECEIVER_INDEX rcv_idx_;
        tGnssSignal *signal_[NSYS];

    private:
        cTime ts_,te_;
        vector<tEpochSatUnit> obs_;
        tStaInfoUnit station_;
    };

    typedef struct {
        cTime t_tag[4];          // last epoch time
        int n[4];                // number of epochs
        double lc_amb[4];        // linear combination average */
        double var_amb[4];       // linear combination variance */
        int fix_cnt;             // fix count */
        char flags[MAX_SAT_NUM]; // fix flags */
    }tLcAmb;

    typedef struct {
        int ref_sat,sat;
        int f;
        int inherit_flag;
    }tDdAmb;

    typedef struct{
        cSat sat;
        cTime t_tag;
        cTime t_trans;
        GNSS_SAT_STAT stat;
        int brd_eph_index;

        unsigned char P_code[MAX_GNSS_USED_FRQ_NUM];
        double raw_P[MAX_GNSS_USED_FRQ_NUM];        // L1 L2 L5
        double raw_L[MAX_GNSS_USED_FRQ_NUM];
        double raw_D[MAX_GNSS_USED_FRQ_NUM];
        unsigned char raw_S[MAX_GNSS_USED_FRQ_NUM];
        unsigned char LLI[MAX_GNSS_USED_FRQ_NUM];
        double csc_P[MAX_GNSS_USED_FRQ_NUM];
        double cor_P[MAX_GNSS_USED_FRQ_NUM];        // corrected code bias BDS satellite-specific multipath
        double cor_L[MAX_GNSS_USED_FRQ_NUM];        // corrected phase bias for PPP-AR
        double cor_D[MAX_GNSS_USED_FRQ_NUM];
        double lam[MAX_GNSS_USED_FRQ_NUM];
        double frq[MAX_GNSS_USED_FRQ_NUM];

        Vector3d brd_pos;
        Vector3d brd_vel;
        Vector2d brd_clk;  // clk, clk-rate
        Vector3d pre_pos;
        Vector3d pre_vel;
        Vector2d pre_clk;

        Vector2d trp_dry_delay; // slant_dry,map_dry
        Vector4d trp_wet_delay; // slant_wet,map_wet,grand_e,grand_n
        Vector2d ion_delay; // L1_slant_ion, map_ion;
        double clk_rel;
        double sagnac;
        double shapiro;

        double code_bias[MAX_GNSS_USED_FRQ_NUM];
        double bd2_mp[3];
        double phase_wp;
        double float_amb[MAX_GNSS_USED_FRQ_NUM]; //L1,L2,L5
        double fix_amb[MAX_GNSS_USED_FRQ_NUM];

        double omc_P[MAX_GNSS_USED_FRQ_NUM];
        double omc_L[MAX_GNSS_USED_FRQ_NUM];

        double prior_res_P[MAX_GNSS_USED_FRQ_NUM];
        double prior_res_L[MAX_GNSS_USED_FRQ_NUM];
        double post_res_P[MAX_GNSS_USED_FRQ_NUM];
        double post_res_L[MAX_GNSS_USED_FRQ_NUM];

        Vector3d sig_vec;
        Vector2d el_az;
        Vector2d raw_mw;  //L1_L2, L1_L5
        Vector2d sm_mw;   //L1_L2, L1_L5
        Vector2d var_mw;
        int mw_idx[2];
        Vector2d gf;          //L1_L2, L1_L5
        Vector2d multipath_comb;
        double tdcp;
        double cmc_P;             // code minus carrier
        double cor_if_P[2];     // SF-IF or DF-IF or TF-IF1/TF-IF2
        double cor_if_L[2];
        int tf_if_idx[3];
        tLcAmb lc_amb;

        int svh;
        double brd_eph_var;
        double pre_eph_var;
        double ion_var;
        double trp_var;
        double c_var_factor[MAX_GNSS_USED_FRQ_NUM];
        double p_var_factor[MAX_GNSS_USED_FRQ_NUM];
        int outc[MAX_GNSS_USED_FRQ_NUM];
        int lock[MAX_GNSS_USED_FRQ_NUM];
        unsigned char rejc[MAX_GNSS_USED_FRQ_NUM]; // reject flag for residual
        unsigned char slip[MAX_GNSS_USED_FRQ_NUM];
        unsigned char fix[MAX_GNSS_USED_FRQ_NUM];
        unsigned char vsat[MAX_GNSS_USED_FRQ_NUM];
        unsigned char res_idx[MAX_GNSS_USED_FRQ_NUM];
    }tSatInfoUnit;

    class cGnssObsOperator {
    public:
        cGnssObsOperator();
        ~cGnssObsOperator();


    public:
        void ReAlignObs(tPPPLibConf C,tSatInfoUnit& sat_info, tSatObsUnit sat_obs,int f,int raw_idx,int frq_idx,int *glo_frq);

        double GnssObsIfComb(double obs1,double obs2,double f1,double f2);
        double GnssObsMwComb(double obs_P1,double obs_P2,double obs_L1,double obs_L2,double lam1,double lam2);
        double GnssObsGfComb(double obs1,double obs2,double lam1,double lam2);
        double GnssObsCmcComb(double obs_P,double obs_L1,double obs_L2,double f1,double f2);
        double GnssObsTdComb(double cur_obs,double pre_obs);
        void MakeGnssObsComb(tPPPLibConf C,GNSS_OBS_COMB type,tSatInfoUnit* sat_info,const tSatInfoUnit previous_sat_info);
        double GnssSdObs(const tSatInfoUnit& sat1,const tSatInfoUnit& sat2,int f,GNSS_OBS type);
        double GnssObsLinearComb(tPPPLibConf C,int i,int j,int k,tSatInfoUnit& sat_info, GNSS_OBS type,double* var);
        double LinearCombLam(int i,int j,int k,tSatInfoUnit& sat_info);

        void MwCycleSlip(tPPPLibConf C,double sample_dt,double dt,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat,tTime last_time);
        void GfCycleSlip(tPPPLibConf C,double sample_dt,double dt,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat);
        bool LliCycleSlip(tPPPLibConf C, tSatInfoUnit& sat_info,int nf,double tt,RECEIVER_INDEX rcv);
        void SmoothMw(tPPPLibConf C,tSatInfoUnit* sat_info,tSatInfoUnit* base_sat);

    };

    typedef struct {
        cTime t_tag;
        int num_use_sat;
        int rcv_status;
        int no;

        int sol_level;
        int vel_flag;
        Vector3d pos;
        Vector3d llh;
        Vector3d vel;
        Vector3d rb_delta;
        double dop[4];
        float cov[8];
        float sig[6]; //enu
    }tGsofUnit;

    typedef struct {
        int n,nmax;
        tGsofUnit *data;
    }tGsofs;

    double GeoDist(Vector3d sat_pos,Vector3d rec_pos,Vector3d& sig_vec);
    double SatElAz(Vector3d rec_blh,Vector3d sig_vec,Vector2d& el_az);
    double GnssMeasVar(tPPPLibConf C, GNSS_OBS obs_type,tSatInfoUnit sat_info);
    double SunMoonPos(cTime ut1t,const double *erp_val,Vector3d& sun_pos, Vector3d& moon_pos);
}



#endif //PPPLIB_GNSSFUNC_H
