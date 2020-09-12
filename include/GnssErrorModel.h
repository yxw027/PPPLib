//
// Created by cc on 7/19/20.
//

#ifndef PPPLIB_GNSSERRORMODEL_H
#define PPPLIB_GNSSERRORMODEL_H

#include "GnssFunc.h"

using namespace PPPLib;

namespace PPPLib{

    class cGnssModel{
    public:
        cGnssModel();
        virtual ~cGnssModel();

    public:
        virtual void InitErrModel(tPPPLibConf C);
        void InitSatInfo(tSatInfoUnit* sat_info,Vector3d* blh);
        virtual void UpdateSatInfo();

    public:
//        tNav nav_;
        tSatInfoUnit* sat_info_;
        Vector3d blh_;
        tPPPLibConf PPPLibC_;
    };

    class cTrpDelayModel: public cGnssModel {
    public:
        cTrpDelayModel();
        cTrpDelayModel(Vector3d blh,tSatInfoUnit& sat_info);
        ~cTrpDelayModel();

    public:
        Vector2d GetTrpError(double humi,double *x,int it);
        Vector2d GetSaasTrp(double humi,Vector2d* sat_trp_dry, Vector4d* sat_trp_wet);
        void UpdateSatInfo() override;
        void InitTrpModel(Vector3d& blh);

    private:
        bool SaasModel(double humi);
        Vector4d EstTrpWet(double humi,double *x,int it);
        void TrpMapNeil(cTime t,double el);

    private:
        Vector2d zenith_trp_={0,0};                   //dry, wet
        Vector2d slant_trp_dry_={0,0};                //dry, map_dry
        Vector4d slant_trp_wet_={0,0,0,0};      //wet, map_wet,grand_e,grand_n
    };

    class cIonDelayModel: public cGnssModel {
    public:
        cIonDelayModel();
        cIonDelayModel(Vector3d blh,tSatInfoUnit& sat_info,tNav nav);
        ~cIonDelayModel();

    public:
        Vector2d GetIonError();
        Vector2d GetKlobIon();
        void UpdateSatInfo() override;
        void InitIondelayModel(double ion_para[NSYS][8]);

    private:
        bool KlobModel();
        bool GimModel();
        bool IonFreeModel();
        void IonMapFunc();
        bool IonEstModel();

    private:
        double ion_para_[NSYS][8]={{0}};
        Vector2d ion_delay_;     // L1_ion/ion_map
    };

    class cCbiasModel:public cGnssModel {
    public:
        cCbiasModel();
        cCbiasModel(tNav nav,tSatInfoUnit& sat_info);
        ~cCbiasModel();

    public:
        void GetCodeBias();
        void InitCbiasModel(double cbias[MAX_SAT_NUM][MAX_GNSS_CODE_BIAS_PAIRS],vector<tBrdEphUnit>& brd_eph);
        void UpdateSatInfo() override;

    private:
        void TgdModel();
        void BsxModel();

    private:
        double code_bias_[MAX_SAT_NUM][MAX_GNSS_CODE_BIAS_PAIRS]={{0}};
        vector<tBrdEphUnit> brd_eph_;

    private:
        double cbias_[NSYS+1][MAX_GNSS_FRQ_NUM]={{0}};
    };

    class cAntModel:public cGnssModel {
    public:
        cAntModel();
        ~cAntModel();

    private:
        double InterpPcv(double ang,const double *var);
        double InterpAziPcv(const tAntUnit& ant,double az,double ze,int f);
        void SatPcvModel(tSatInfoUnit* sat_info,double nadir,double* dant);

    public:
        void InitAntModel(tPPPLibConf C,tAntUnit sat_ant[MAX_SAT_NUM],tAntUnit rec_ant[2]);
        void SatPcvCorr(tSatInfoUnit* sat_info,Vector3d rec_pos,double* pcv_dants);
        void RecAntCorr(tSatInfoUnit* sat_info,double* dantr,RECEIVER_INDEX rec_idx);
    private:
        tAntUnit *sat_ant_;
        tAntUnit *rec_ant_;
    };

    class cEphModel:public cGnssModel {
    public:
        cEphModel();
        ~cEphModel();

    private:

        // broadcast
        int SelectBrdEph(tSatInfoUnit* sat_info,int iode);
        int SelectGloBrdEph(tSatInfoUnit* sat_info,int iode);
        void BrdSatClkErr(tSatInfoUnit* sat_info,tBrdEphUnit eph);
        void BrdGloSatClkErr(tSatInfoUnit* sat_info,tBrdGloEphUnit glo_eph);
        bool BrdClkError(tSatInfoUnit* sat_info);

        double EphUra(int ura);
        double ClkRelErr(double mu,double sinE,tBrdEphUnit eph);
        void BrdSatPos(tSatInfoUnit& sat_info,tBrdEphUnit eph);
        void GloDefEq(const VectorXd x,VectorXd& x_dot,const Vector3d acc);

        void GloOrbit(double t,VectorXd& x,const Vector3d acc);
        void BrdGloSatPos(tSatInfoUnit& sat_info,tBrdGloEphUnit glo_eph);
        bool BrdEphModel(tSatInfoUnit* sat_info);

        // precise
        void SatPcoCorr(tSatInfoUnit* sat_info,Vector3d sat_pos,Vector3d& dant);
        bool PreSatClkCorr(tSatInfoUnit* sat_info);
        bool PreSatPos(tSatInfoUnit* sat_info);
        bool PreEphModel(tSatInfoUnit* sat_info);

    public:
        bool EphCorr(vector<tSatInfoUnit>& epoch_sat_info);
        void InitEphModel(vector<tBrdEphUnit>& brd_eph,vector<tBrdGloEphUnit>& brd_glo_eph,vector<tPreOrbUnit>& pre_orb,
                vector<tPreClkUnit>& pre_clk,tAntUnit* sat_ant);

    private:
        vector<tPreClkUnit> pre_clk_;
        vector<tPreOrbUnit> pre_orb_;
        vector<tBrdEphUnit> brd_eph_;
        vector<tBrdGloEphUnit> brd_glo_eph_;
        tAntUnit *sat_ant_;
    };

    class cTidModel:public cGnssModel {
    public:
        cTidModel();
        ~cTidModel();

    private:
        int GetErpVal(cTime t);
        void TideSl(const double* eu,const double* rp,double GMp,Vector3d blh,double* dr);
        void IersMeanPole(double *xp,double *yp);
        void TidSolid();
        void TidOcean();
        void TidPole();

    public:
        void InitTidModel(tPPPLibConf C,vector<tErpUnit> erp,double ocean_par[2][66]);
        void TidCorr(cTime t,Vector3d rec_pos,Vector3d& dr);

    private:
        cTime tut_;
        double re_;
        Vector3d blh_;
        Matrix3d E_;
        Vector3d tid_dr_[3];
        Vector3d denu_[2];
        Vector3d sun_pos_;
        Vector3d moon_pos_;
        double gmst_;
        double erp_val_[5];

    private:
        vector<tErpUnit> erp_;
        double ocean_par_[2][6*11];
        RECEIVER_INDEX rec_idx_;
    };

    class cGnssErrCorr {
    public:
        cGnssErrCorr();
        ~cGnssErrCorr();


    private:
        bool SatYaw(tSatInfoUnit& sat_info,Vector3d& exs,Vector3d& eys);

    public:
        void InitGnssErrCorr(tPPPLibConf C, tNav* nav);
        void BD2MultipathModel(tSatInfoUnit* sat_info);
        double SagnacCorr(Vector3d sat_xyz,Vector3d rec_xyz);
        double ShapiroCorr(int sys,Vector3d sat_xyz,Vector3d rec_xyz);
        void PhaseWindUp(tSatInfoUnit& sat_info,Vector3d rec_xyz);

    public:
        cEphModel      eph_model_;
        cCbiasModel    cbia_model_;
        cTrpDelayModel trp_model_;
        cIonDelayModel ion_model_;
        cAntModel      ant_model_;
        cTidModel      tid_model_;

    private:
        double bd2_mp_[3];
    };


}



#endif //PPPLIB_GNSSERRORMODEL_H
