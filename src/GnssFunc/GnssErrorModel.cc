//
// Created by cc on 7/19/20.
//

#include "GnssErrorModel.h"

namespace PPPLib{

    cGnssModel::cGnssModel() {}

    cGnssModel::~cGnssModel() {}

    void cGnssModel::UpdateSatInfo() {}

    void cGnssModel::InitErrModel(tPPPLibConf C) {
        PPPLibC_=C;
    }

    void cGnssModel::InitSatInfo(tSatInfoUnit *sat_info, Vector3d* blh) {
        sat_info_=sat_info;
        if(blh) blh_=*blh;
    }

    cTrpDelayModel::cTrpDelayModel() {}

    cTrpDelayModel::cTrpDelayModel(Vector3d blh, tSatInfoUnit &sat_info) {
        blh_=blh;
        sat_info_=&sat_info;
    }

    cTrpDelayModel::~cTrpDelayModel() {}

    Vector2d cTrpDelayModel::GetSaasTrp(double humi,Vector2d* sat_trp_dry,Vector4d* sat_trp_wet){
        SaasModel(humi);
        if(sat_trp_dry) sat_trp_dry=&slant_trp_dry_;
        if(sat_trp_wet) sat_trp_wet=&slant_trp_wet_;
        return zenith_trp_;
    }

    Vector2d cTrpDelayModel::GetTrpError(double humi,double* x,int it) {
        for(int i=0;i<4;i++) sat_info_->trp_wet_delay[i]=0.0;
        for(int i=0;i<2;i++) sat_info_->trp_dry_delay[i]=0.0;
        if(PPPLibC_.gnssC.trp_opt==TRP_SAAS) GetSaasTrp(humi, nullptr,nullptr);
        else if(PPPLibC_.gnssC.trp_opt>=TRP_EST_WET){
            EstTrpWet(humi,x,it);
            sat_info_->trp_var=0.0001;
        }
        return zenith_trp_;
    }

    void cTrpDelayModel::UpdateSatInfo() {
        sat_info_->trp_dry_delay=slant_trp_dry_;
        sat_info_->trp_wet_delay=slant_trp_wet_;
        for(int i=0;i<4;i++) slant_trp_wet_[i]=0.0;
        slant_trp_dry_[0]=slant_trp_dry_[1]=0.0;
    }

    void cTrpDelayModel::InitTrpModel(Vector3d &blh) {
        blh_=blh;
    }

    bool cTrpDelayModel::SaasModel(double humi) {
        const double temp0=15.0;
        double hgt,pres,temp,e,z;
        double el=sat_info_->el_az[0];

        if(blh_[2]<-100.0||1E6<blh_[2]||el<=0) return false;
        if(blh_[2]>=1.0/2.2557E-5) return false;
        hgt=blh_[2]<0.0?0.0:blh_[2];
        if(hgt>15000.0) hgt=15000.0;

        pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
        temp=temp0-6.5E-3*hgt+273.16;
        e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));
        z=PI/2.0-el;
        zenith_trp_[0]=0.0022768*pres/(1.0-0.00266*cos(2.0*blh_[0])-0.00028*hgt/1E3);
        zenith_trp_[1]=0.0022770*(1255.0/temp+0.05)*e;
        double map=1.0/cos(z);
        slant_trp_dry_[0]=zenith_trp_[0]*map;
        slant_trp_dry_[1]=map;
        slant_trp_wet_[0]=zenith_trp_[1]*map;
        slant_trp_wet_[1]=map;
        sat_info_->trp_var=SQR(0.3);
    }

    Vector4d cTrpDelayModel::EstTrpWet(double humi,double *x,int it) {
        GetSaasTrp(humi, nullptr, nullptr);
        TrpMapNeil(sat_info_->t_tag,sat_info_->el_az[0]);

        double grad_n=0.0,grad_e=0.0;
        if(PPPLibC_.gnssC.trp_opt==TRP_EST_GRAD&&sat_info_->el_az[0]>0){
            double cotz=1.0/tan(sat_info_->el_az[0]);
            grad_n=slant_trp_wet_[1]*cotz*cos(sat_info_->el_az[1]);
            grad_e=slant_trp_wet_[1]*cotz*sin(sat_info_->el_az[1]);
            slant_trp_wet_[1]+=grad_n*x[it+1]+grad_e*x[it+2];
            slant_trp_wet_[2]=grad_n*x[it+1];
            slant_trp_wet_[3]=grad_e*x[it+2];
        }

        slant_trp_dry_[0]=zenith_trp_[0]*slant_trp_dry_[1];
        slant_trp_wet_[0]=slant_trp_wet_[1]*x[it];
    }

    static double Interpc(const double coef[],double lat) {
        int i=(int)(lat/15.0);

        if(i<1) return coef[0];else if(i>4) return coef[4];
        return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
    }

    static double Mapf(double el,double a,double b,double c) {
        double sinel=sin(el);
        return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
    }

    void cTrpDelayModel::TrpMapNeil(cTime t, double el) {
        const double coef[][5]={
                { 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
                { 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
                { 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},

                { 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
                { 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
                { 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},

                { 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
                { 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
                { 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
        };
        const double aht[]={ 2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */

        double y,cosy,ah[3],aw[3],dm,lat=blh_[0]*R2D,hgt=blh_[2];
        int i;

        if(el<=0.0) {return;}
        y=(t.Time2Doy()-28.0)/365.25+(lat<0.0?0.5:0.0);
        cosy=cos(2.0*PI*y);
        lat=fabs(lat);

        for(i=0;i<3;i++){
            ah[i]=Interpc(coef[i],lat)-Interpc(coef[i+3],lat)*cosy;
            aw[i]=Interpc(coef[i+6],lat);
        }
        dm=(1.0/sin(el)-Mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;
        slant_trp_wet_[1]=Mapf(el,aw[0],aw[1],aw[2]);
        slant_trp_dry_[1]=Mapf(el,ah[0],ah[1],ah[2])+dm;
    }

    cIonDelayModel::cIonDelayModel() {}

    cIonDelayModel::cIonDelayModel(Vector3d blh, tSatInfoUnit &sat_info, tNav nav) {
        blh_=blh;
        sat_info_=&sat_info;
    }

    cIonDelayModel::~cIonDelayModel() {}

    Vector2d cIonDelayModel::GetIonError() {
        sat_info_->ion_delay[0]=sat_info_->ion_delay[1]=0.0;
        if(PPPLibC_.gnssC.ion_opt==ION_IF||PPPLibC_.gnssC.ion_opt==ION_IF_DUAL) IonFreeModel();
        else if(PPPLibC_.gnssC.ion_opt==ION_KLB) KlobModel();
        else if(PPPLibC_.gnssC.ion_opt>ION_EST){

        }
        return ion_delay_;
    }

    Vector2d cIonDelayModel::GetKlobIon() {
        KlobModel();
        return ion_delay_;
    }

    void cIonDelayModel::UpdateSatInfo() {
        sat_info_->ion_delay=ion_delay_;
        ion_delay_[0]=ion_delay_[1]=0.0;
    }

    void cIonDelayModel::InitIondelayModel(double ion_para[NSYS][8]) {
        for(int i=0;i<NSYS;i++){
            for(int j=0;j<8;j++)
                ion_para_[i][j]=ion_para[i][j];
        }
    }

    bool cIonDelayModel::KlobModel() {
        const double ion_default[]={ /* 2004/1/1 */
                0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
                0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
        };
        double tt,psi,phi,lam,amp,per,x;
        double el=sat_info_->el_az[0];
        double az=sat_info_->el_az[1];
        int week;
        cTime obs_time=sat_info_->t_tag;

        if (blh_[2]<-1E3||el<=0) return false;

        VectorXd ion_par=VectorXd::Map(ion_para_[SYS_INDEX_GPS],8);
        if(ion_par.norm()<=0.0) ion_par=VectorXd::Map(ion_default,8);

        /* earth centered angle (semi-circle) */
        psi=0.0137/(el/PI+0.11)-0.022;

        /* subionospheric latitude/longitude (semi-circle) */
        phi=blh_[0]/PI+psi*cos(az);
        // latitude boundary protection
        if (phi> 0.416) phi= 0.416;
        else if (phi<-0.416) phi=-0.416;
        lam=blh_[1]/PI+psi*sin(az)/cos(phi*PI);

        /* geomagnetic latitude (semi-circle) */
        phi+=0.064*cos((lam-1.617)*PI);

        /* local time (s) */
        tt=43200.0*lam+obs_time.Time2Gpst(&week, nullptr,SYS_GPS);
        tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */

        /* slant factor */
        double map_ion=0.0;
        map_ion=sat_info_->ion_delay[1]=1.0+16.0*pow(0.53-el/PI,3.0);

        /* ionospheric delay */
        amp=ion_par[0]+phi*(ion_par[1]+phi*(ion_par[2]+phi*ion_par[3]));
        per=ion_par[4]+phi*(ion_par[5]+phi*(ion_par[6]+phi*ion_par[7]));
        amp=amp<0.0?0.0:amp;
        per=per<72000.0?72000.0:per;
        x=2.0*PI*(tt-50400.0)/per;

        // GPS L1
        double sion=CLIGHT*map_ion*(fabs(x)<1.57 ? 5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
        double ion_var=SQR(sion*0.5);

        double lam_G1=CLIGHT/FREQ_GPS_L1;
        ion_delay_[0]=sion*SQR(sat_info_->lam[0]/lam_G1);
        ion_delay_[1]=map_ion;
        sat_info_->ion_var=ion_var*SQR(sat_info_->lam[0]/lam_G1);

        return true;
    }

    bool cIonDelayModel::IonFreeModel() {
        ion_delay_[0]=ion_delay_[1]=0.0;
    }

    void cIonDelayModel::IonMapFunc() {
        double el=sat_info_->el_az[0];
        if(blh_[2]>=350000.0) ion_delay_[1]=1.0;
        else{
            ion_delay_[1]=1.0/cos(asin((WGS84_EARTH_LONG_RADIUS+blh_[2])/(WGS84_EARTH_LONG_RADIUS+350000.0)*sin(PI/2.0-el)));
        }
    }

    bool cIonDelayModel::IonEstModel() {
        IonMapFunc();
    }

    cCbiasModel::cCbiasModel() {}

    cCbiasModel::cCbiasModel(tNav nav,tSatInfoUnit& sat_info) {
//        nav_=nav;
        sat_info_=&sat_info;
    }

    cCbiasModel::~cCbiasModel() {}

    void cCbiasModel::GetCodeBias() {
//        if(PPPLibC_.fileC.cbias.empty()){
//            TgdModel();
//            return;
//        }
        BsxModel();
    }

    void cCbiasModel::InitCbiasModel(double cbias[MAX_SAT_NUM][MAX_GNSS_CODE_BIAS_PAIRS],vector<tBrdEphUnit>& brd_eph) {
        for(int i=0;i<MAX_SAT_NUM;i++){
            for(int j=0;j<MAX_GNSS_CODE_BIAS_PAIRS;j++)
                code_bias_[i][j]=cbias[i][j];
        }
        brd_eph_=brd_eph;
    }

    void cCbiasModel::UpdateSatInfo() {
        int sys=sat_info_->sat.sat_.sys;
        int sys_idx=sat_info_->sat.sat_.sys_idx;
        int *frqs=PPPLibC_.gnssC.gnss_frq[sys_idx];
        double *bias=cbias_[sys_idx];

        for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
            sat_info_->code_bias[i]=0.0;
        }

        if(sys==SYS_BDS&&sat_info_->sat.sat_.prn>18){
            frqs=PPPLibC_.gnssC.gnss_frq[NSYS];
            bias=cbias_[NSYS];
        }

        sat_info_->code_bias[0]=bias[0];
        if(PPPLibC_.gnssC.frq_opt==FRQ_DUAL){
            sat_info_->code_bias[0]=bias[0];
            sat_info_->code_bias[1]=bias[1];
        }
        else if(PPPLibC_.gnssC.frq_opt==FRQ_TRIPLE){
            sat_info_->code_bias[0]=bias[0];
            sat_info_->code_bias[1]=bias[1];
            sat_info_->code_bias[2]=bias[2];
        }
    }

    void cCbiasModel::TgdModel() {
        int i;
        double gamma=0.0;
        double tgd1=0.0,tgd2=0.0;

        for(i=0;i<NSYS+1;i++){
            for(int j=0;j<MAX_GNSS_FRQ_NUM;j++){
                cbias_[i][j]=0.0;
            }
        }

        for(i=0;i;i++){
            if(brd_eph_.at(i).sat.sat_.no!=sat_info_->sat.sat_.no) continue;
            tgd1=CLIGHT*(brd_eph_.at(i).tgd[0]);
            tgd2=CLIGHT*(brd_eph_.at(i).tgd[1]);
            break;
        }
        if(sat_info_->sat.sat_.sys==SYS_GPS||sat_info_->sat.sat_.sys==SYS_QZS){
            gamma=SQR(FREQ_GPS_L1/FREQ_GPS_L2);
            cbias_[SYS_INDEX_GPS][0]=tgd1;                        // L1   tgd=dcb/(1.0-gamma)
            cbias_[SYS_INDEX_GPS][1]=tgd1*gamma;                  // L2
        }
        else if(sat_info_->sat.sat_.sys==SYS_GAL){
            gamma=SQR(FREQ_GAL_E5A/FREQ_GAL_E1);
            cbias_[SYS_INDEX_GAL][0]=gamma*tgd1;                   // E1
            cbias_[SYS_INDEX_GAL][1]=tgd1;                         // E5a
            cbias_[SYS_INDEX_GAL][2]=gamma*tgd1+(1.0-gamma)*tgd2;  // E5b
        }
        else if(sat_info_->sat.sat_.sys==SYS_BDS){
            if(sat_info_->sat.sat_.prn>18){
                cbias_[NSYS][0]=-tgd1;                             // BD3-B1I
                cbias_[NSYS][1]=0.0;
                cbias_[NSYS][2]=0.0;                               // BD3-B3I
            }
            else{
                cbias_[SYS_INDEX_BDS][0]=-tgd1;                    // BD2-B1I
                cbias_[SYS_INDEX_BDS][1]=-tgd2;                    // BD2-B2I
                cbias_[SYS_INDEX_BDS][2]=0.0;                      // BD2-B3I
            }
        }
    }

    void cCbiasModel::BsxModel() {
        int sys=sat_info_->sat.sat_.sys;
        int sat_no=sat_info_->sat.sat_.no;
        int sys_idx=sat_info_->sat.sat_.sys_idx;

        if(sys==SYS_BDS&&sat_info_->sat.sat_.prn>18) sys_idx=NSYS;

        for(int j=0;j<MAX_GNSS_FRQ_NUM;j++)
            cbias_[sys_idx][j]=0.0;

        if(sys==SYS_GPS){
            double DCB_p1p2=code_bias_[sat_no-1][GPS_C1WC2W];
            double a=(SQR(FREQ_GPS_L1)-SQR(FREQ_GPS_L2));
            double alpha=SQR(FREQ_GPS_L1)/a;
            double beta=-SQR(FREQ_GPS_L2)/a;

            for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                    cbias_[SYS_INDEX_GPS][i]=0.0;
                    continue;
                }
                if(i==0){
                    cbias_[SYS_INDEX_GPS][0]=beta*DCB_p1p2;     // L1
                    if(sat_info_->P_code[i]==GNSS_CODE_L1C) cbias_[SYS_INDEX_GPS][i]+=code_bias_[sat_no-1][GPS_C1CC1W];
                }
                else if(i==1){
                    cbias_[SYS_INDEX_GPS][1]=-alpha*DCB_p1p2;   // L2
                    if(sat_info_->P_code[i]==GNSS_CODE_L2C) cbias_[SYS_INDEX_GPS][i]+=code_bias_[sat_no-1][GPS_C2CC2W];
                    else if(sat_info_->P_code[i]==GNSS_CODE_L2S) cbias_[SYS_INDEX_GPS][i]-=code_bias_[sat_no-1][GPS_C2WC2S];
                    else if(sat_info_->P_code[i]==GNSS_CODE_L2L) cbias_[SYS_INDEX_GPS][i]-=code_bias_[sat_no-1][GPS_C2WC2L];
                    else if(sat_info_->P_code[i]==GNSS_CODE_L2X) cbias_[SYS_INDEX_GPS][i]-=code_bias_[sat_no-1][GPS_C2WC2X];
                }
                else if(i==2){
                    double DCB_p1p3=0.0;
                    double beta_13=-SQR(FREQ_GPS_L5)/(SQR(FREQ_GPS_L1)-SQR(FREQ_GPS_L5));
                    if(sat_info_->P_code[i]==GNSS_CODE_L5Q) DCB_p1p3=code_bias_[sat_no-1][GPS_C1CC5Q]-code_bias_[sat_no-1][GPS_C1CC1W];
                    else if(sat_info_->P_code[i]==GNSS_CODE_L5X) DCB_p1p3=code_bias_[sat_no-1][GPS_C1CC5X]-code_bias_[sat_no-1][GPS_C1CC1W];
                    cbias_[SYS_INDEX_GPS][2]=beta_13*DCB_p1p2-DCB_p1p3;
                }
            }
        }
        else if(sys==SYS_BDS){

            if(PPPLibC_.gnssC.eph_opt==EPH_BRD){
                // base on b3
                if(sat_info_->sat.sat_.prn>18){
                    //BD3
                    double DCB_b1b3=code_bias_[sat_no-1][BD3_C2IC6I];
                    for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                        if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                            cbias_[SYS_INDEX_BDS][i]=0.0;
                            continue;
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L2I){
                            cbias_[NSYS][i]=DCB_b1b3;
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L7I){
                            cbias_[NSYS][i]=0.0;
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L6I){
                            cbias_[NSYS][i]=0.0;
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L1X){
                            cbias_[NSYS][i]=code_bias_[sat_no-1][BD3_C1XC6I];
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L1P){
                            cbias_[NSYS][i]=code_bias_[sat_no-1][BD3_C1PC6I];
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L1D){
                            cbias_[NSYS][i]=code_bias_[sat_no-1][BD3_C1DC6I];
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L5X){
                            cbias_[NSYS][i]=code_bias_[sat_no-1][BD3_C1XC6I]-code_bias_[sat_no-1][BD3_C1XC5X];
                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L5P){
                            cbias_[NSYS][i]=code_bias_[sat_no-1][BD3_C1PC6I]-code_bias_[sat_no-1][BD3_C1PC5P];

                        }
                        else if(sat_info_->P_code[i]==GNSS_CODE_L5D){
                            cbias_[NSYS][i]=code_bias_[sat_no-1][BD3_C1DC6I]-code_bias_[sat_no-1][BD3_C1DC5D];
                        }
                    }
                }
                else{
                    //BD2
                    double DCB_b1b2=code_bias_[sat_no-1][BD2_C2IC7I];
                    double DCB_b1b3=code_bias_[sat_no-1][BD2_C2IC6I];
                    for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                        if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                            cbias_[SYS_INDEX_BDS][i]=0.0;
                            continue;
                        }
                        if(i==0){
                            cbias_[SYS_INDEX_BDS][i]=DCB_b1b3;
                        }
                        else if(i==1){
                            cbias_[SYS_INDEX_BDS][i]=DCB_b1b3-DCB_b1b2;
                        }
                        else if(i==2){
                            cbias_[SYS_INDEX_BDS][i]=0.0;
                        }
                    }
                }
            }
            else if(PPPLibC_.gnssC.eph_opt==EPH_PRE){
                if(PPPLibC_.gnssC.ac_opt==AC_COM){
                    // base on b1b2
                    for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                        if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                            cbias_[SYS_INDEX_BDS][i]=0.0;
                            continue;
                        }
                        if(sat_info_->sat.sat_.prn>18){
                            cbias_[SYS_INDEX_BDS][i]=0.0;
                            continue;
                        }
                        double DCB_b1b2=code_bias_[sat_no-1][BD2_C2IC7I];
                        double DCB_b1b3=code_bias_[sat_no-1][BD2_C2IC6I];
                        double a=(SQR(FREQ_BDS_B1)-SQR(FREQ_BDS_B2));
                        double alpha=SQR(FREQ_BDS_B1)/a;
                        double beta=-SQR(FREQ_BDS_B2)/a;
                        if(i==0){
                            cbias_[SYS_INDEX_BDS][i]=beta*DCB_b1b2;
                        }
                        else if(i==1){
                            cbias_[SYS_INDEX_BDS][i]=-alpha*DCB_b1b2;
                        }
                        else if(i==2){
                            double beta_13=-SQR(FREQ_BDS_B3)/(SQR(FREQ_BDS_B1)-SQR(FREQ_BDS_B3));
                            cbias_[SYS_INDEX_BDS][i]=beta_13*DCB_b1b2-DCB_b1b3;
                        }
                    }
                }
                else if(PPPLibC_.gnssC.ac_opt==AC_WUM||PPPLibC_.gnssC.ac_opt==AC_GBM){
                    // base on b1b3
                    double a=(SQR(FREQ_BDS_B1)-SQR(FREQ_BDS_B3));
                    double alpha=SQR(FREQ_BDS_B1)/a;
                    double beta=-SQR(FREQ_BDS_B3)/a;
                    for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                        if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                            cbias_[SYS_INDEX_BDS][i]=0.0;
                            continue;
                        }
                        if(sat_info_->sat.sat_.prn>18){
                            double DCB_b1b3=code_bias_[sat_no-1][BD3_C2IC6I];
                            if(i==0){
                                cbias_[SYS_INDEX_BDS][i]=beta*DCB_b1b3;
                            }
                            else if(i==1){
                                cbias_[SYS_INDEX_BDS][i]=0.0;
                            }
                            else if(i==2){
                                cbias_[SYS_INDEX_BDS][i]=-alpha*DCB_b1b3;
                            }
                            else if(i==3){
                                //B1C
                                double beta_13=-SQR(FREQ_BDS_B1C)/(SQR(FREQ_BDS_B1)-SQR(FREQ_BDS_B3));
                                double DCB_b1b2=0.0;
                                if(sat_info_->P_code[i]==GNSS_CODE_L1X) DCB_b1b2=code_bias_[sat_no-1][BD3_C1XC6I];
                                else if(sat_info_->P_code[i]==GNSS_CODE_L1P) DCB_b1b2=code_bias_[sat_no-1][BD3_C1PC6I];
                                else if(sat_info_->P_code[i]==GNSS_CODE_L1D) DCB_b1b2=code_bias_[sat_no-1][BD3_C1DC6I];
                                double dcb=DCB_b1b3-DCB_b1b2;
                                cbias_[SYS_INDEX_BDS][i]=beta_13*DCB_b1b3-dcb;
                            }
                            else if(i==4){
                                //B2a
                                double beta_13=-SQR(FREQ_BDS_B2A)/(SQR(FREQ_BDS_B1)-SQR(FREQ_BDS_B2A));
                                double DCB_b1b2=0.0;
                                if(sat_info_->P_code[i]==GNSS_CODE_L5X) DCB_b1b2=code_bias_[sat_no-1][BD3_C1XC6I]-code_bias_[sat_no-1][BD3_C1XC5X];
                                else if(sat_info_->P_code[i]==GNSS_CODE_L5P) DCB_b1b2=code_bias_[sat_no-1][BD3_C1PC6I]-code_bias_[sat_no-1][BD3_C1PC5P];
                                else if(sat_info_->P_code[i]==GNSS_CODE_L5D) DCB_b1b2=code_bias_[sat_no-1][BD3_C1DC6I]-code_bias_[sat_no-1][BD3_C1DC5D];
                                double dcb=DCB_b1b3-DCB_b1b2;
                                cbias_[SYS_INDEX_BDS][i]=beta_13*DCB_b1b3-dcb;
                            }
                            else if(i==5){
                                //B2b
                                cbias_[SYS_INDEX_BDS][i]=0.0;
                            }

                        }
                        else{
                            double DCB_b1b2=code_bias_[sat_no-1][BD2_C2IC7I];
                            double DCB_b1b3=code_bias_[sat_no-1][BD2_C2IC6I];
                            if(i==0){
                                cbias_[SYS_INDEX_BDS][i]=beta*DCB_b1b3;
                            }
                            else if(i==1){
                                double beta_13=-SQR(FREQ_BDS_B2)/(SQR(FREQ_BDS_B1)-SQR(FREQ_BDS_B2));
                                cbias_[SYS_INDEX_BDS][i]=beta_13*DCB_b1b3-DCB_b1b2;
                            }
                            else if(i==2){
                                cbias_[SYS_INDEX_BDS][i]=-alpha*DCB_b1b3;
                            }
                        }
                    }
                }
            }
        }
        else if(sys==SYS_GAL){
            double DCB_p1p2=code_bias_[sat_no-1][GAL_C1CC5Q];
            double a=(SQR(FREQ_GAL_E1)-SQR(FREQ_GAL_E5A));
            double alpha=SQR(FREQ_GAL_E1)/a;
            double beta=-SQR(FREQ_GAL_E5A)/a;

            for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                    cbias_[SYS_INDEX_GAL][i]=0.0;
                    continue;
                }
                if(i==0){
                    if(sat_info_->P_code[i]==GNSS_CODE_L1X) DCB_p1p2=code_bias_[sat_no-1][GAL_C1XC5X];
                    cbias_[SYS_INDEX_GAL][i]=beta*DCB_p1p2;
                }
                else if(i==1){
                    if(sat_info_->P_code[i]==GNSS_CODE_L5X) DCB_p1p2=code_bias_[sat_no-1][GAL_C1XC5X];
                    cbias_[SYS_INDEX_GAL][i]=-alpha*DCB_p1p2;
                }
                else if(i==2){
                    double DCB_p1p3=0.0;
                    double beta_13=-SQR(FREQ_GAL_E5B)/(SQR(FREQ_GAL_E1)-SQR(FREQ_GAL_E5B));
                    if(sat_info_->P_code[i]==GNSS_CODE_L7X){
                        DCB_p1p2=code_bias_[sat_no-1][GAL_C1XC5X];
                        DCB_p1p3=code_bias_[sat_no-1][GAL_C1XC7X];
                    }
                    else if(sat_info_->P_code[i]==GNSS_CODE_L7Q) DCB_p1p3=code_bias_[sat_no-1][GAL_C1CC7Q];
                    cbias_[SYS_INDEX_GAL][i]=beta_13*DCB_p1p2-DCB_p1p3;
                }
            }
        }
        else if(sys==SYS_GLO){
            double DCB_p1p2=code_bias_[sat_no-1][GLO_C1PC2P];
            double a=(SQR(FREQ_GLO_G1)-SQR(FREQ_GLO_G2));
            double alpha=SQR(FREQ_GLO_G1)/a;
            double beta=-SQR(FREQ_GLO_G2)/a;

            for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                    cbias_[SYS_INDEX_GLO][i]=0.0;
                    continue;
                }
                if(i==0){
                    cbias_[SYS_INDEX_GLO][i]=beta*DCB_p1p2;
                    if(sat_info_->P_code[i]==GNSS_CODE_L1C) cbias_[SYS_INDEX_GLO][i]+=code_bias_[sat_no-1][GLO_C1CC1P];
                }
                else if(i==1){
                    cbias_[SYS_INDEX_GLO][i]=-alpha*DCB_p1p2;
                    if(sat_info_->P_code[i]==GNSS_CODE_L2C) cbias_[SYS_INDEX_GLO][i]+=code_bias_[sat_no-1][GLO_C2CC2P];
                }
            }
        }
        else if(sys==SYS_QZS){
            double DCB_p1p2=code_bias_[sat_no-1][QZS_C1CC2L];
            double a=(SQR(FREQ_QZS_L1)-SQR(FREQ_QZS_L2));
            double alpha=SQR(FREQ_QZS_L1)/a;
            double beta=-SQR(FREQ_QZS_L2)/a;

            for(int i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
                if(sat_info_->P_code[i]==GNSS_CODE_NONE){
                    cbias_[SYS_INDEX_QZS][i]=0.0;
                    continue;
                }
                if(i==0){
                    if(sat_info_->P_code[i]==GNSS_CODE_L1X) DCB_p1p2=code_bias_[sat_no-1][QZS_C1XC2X];
                    cbias_[SYS_INDEX_QZS][i]=beta*DCB_p1p2;
                }
                else if(i==1){
                    if(sat_info_->P_code[i]==GNSS_CODE_L2X) DCB_p1p2=code_bias_[sat_no-1][QZS_C1XC2X];
                    cbias_[SYS_INDEX_QZS][i]=-alpha*DCB_p1p2;
                }
                else if(i==2){
                    if(sat_info_->P_code[0]==GNSS_CODE_L1X) DCB_p1p2=code_bias_[sat_no-1][QZS_C1XC2X];
                    double beta_13=-SQR(FREQ_QZS_L1)/(SQR(FREQ_QZS_L1)-SQR(FREQ_QZS_L5));
                    double DCB_p1p3=0.0;
                    if(sat_info_->P_code[i]==GNSS_CODE_L5Q) DCB_p1p3=code_bias_[sat_no-1][QZS_C1CC5Q];
                    else if(sat_info_->P_code[i]==GNSS_CODE_L5X) DCB_p1p3=code_bias_[sat_no-1][QZS_C1CC5X];
                    cbias_[SYS_INDEX_QZS][i]=beta_13*DCB_p1p2-DCB_p1p3;
                }
            }
        }
    }

    cAntModel::cAntModel() {}

    cAntModel::~cAntModel() {}

    void cAntModel::InitAntModel(tPPPLibConf C, tAntUnit *sat_ant, tAntUnit *rec_ant) {
        PPPLibC_=C;
        sat_ant_=sat_ant;
        rec_ant_=rec_ant;
    }

    double cAntModel::InterpAziPcv(const tAntUnit &ant, double az, double ze, int f) {
        double p,q,rpcv=0.0;
        int ize,iaz;
        int i=(int)((ant.zen2-ant.zen1)/ant.dzen)+1;

        ize=(int)((ze-ant.zen1)/(ant.dzen));
        iaz=(int)(az/ant.dazi);

        p=ze/ant.dzen-ize;
        q=az/ant.dazi-iaz;

        rpcv=(1.0-p)*(1.0-q)*ant.pcv[f][(iaz+0)*i+(ize+0)]
             +p*(1.0-q)*ant.pcv[f][(iaz+0)*i+(ize+1)]
             +q*(1.0-p)*ant.pcv[f][(iaz+1)*i+(ize+0)]
             +p*q*ant.pcv[f][(iaz+1)*i+(ize+1)];
        return rpcv;
    }

    double cAntModel::InterpPcv(double ang, const double *var) {
        double a=ang/5.0;
        int i=(int)a;

        if(i<0) return var[0];
        else if(i>=18) return var[18];

        return var[i]*(1.0-a+i)+var[i+1]*(a-i);
    }

    void cAntModel::SatPcvModel(tSatInfoUnit *sat_info, double nadir, double *dant) {
        int i,ii,sys=sat_info->sat.sat_.sys;
        tAntUnit ant=sat_ant_[sat_info->sat.sat_.no-1];

        for(i=0;i<3;i++){
            switch (sys) {
                case SYS_GPS:
                    ii=i;
                    if(i==2) ii=1;
                    break;
                case SYS_BDS:
                    ii=i+MAX_GNSS_FRQ_NUM;
                    if(i==2) {ii=1+MAX_GNSS_FRQ_NUM;}
                    break;
                case SYS_GAL:
                    ii=i+2*MAX_GNSS_FRQ_NUM;
                    if(i==2) {ii=1+2*MAX_GNSS_FRQ_NUM;}
                    break;
                case SYS_GLO:
                    ii=i+3*MAX_GNSS_FRQ_NUM;
                    if(i==2) {ii=1+3*MAX_GNSS_FRQ_NUM;}
                    break;
                case SYS_QZS:
                    ii=i+4*MAX_GNSS_FRQ_NUM;
                    if(i==2) ii=1+4*MAX_GNSS_FRQ_NUM;
                    break;
            }
            dant[i]=InterpPcv(nadir*R2D*5.0,ant.pcv[ii]);
        }
    }

    void cAntModel::SatPcvCorr(tSatInfoUnit *sat_info, Vector3d rec_pos, double *pcv_dants) {
        double ru[3],rz[3],eu[3],ez[3],nadir,cosa;
        int i;

        double rs[3]={sat_info->pre_pos[0],sat_info->pre_pos[1],sat_info->pre_pos[2]};
        for(i=0;i<3;i++){
            ru[i]=rec_pos[i]-rs[i];
            rz[i]=-rs[i];
        }
        if(!NormV3(ru,eu)||!NormV3(rz,ez)) return;

        cosa=Dot(eu,ez,3);
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        nadir=acos(cosa);

        SatPcvModel(sat_info,nadir,pcv_dants);
    }

    void cAntModel::RecAntCorr(tSatInfoUnit *sat_info, double *dantr,RECEIVER_INDEX rec_idx) {
        double e[3],off[3],cosel=cos(sat_info->el_az[0]);
        int i,j,k=0;

        e[0]=sin(sat_info->el_az[1])*cosel;
        e[1]=cos(sat_info->el_az[1])*cosel;
        e[2]=sin(sat_info->el_az[0]);

        int sys=sat_info->sat.sat_.sys,ii;
        for(i=0;i<MAX_GNSS_USED_FRQ_NUM;i++){
            if(sys==SYS_GPS||sys==SYS_BDS||sys==SYS_GAL||sys==SYS_QZS){
                //使用GPS的接收机天线PCV近似替代
                ii=i;
                if(i==2) ii=1; //L5-->L2
            }
            else if(sys==SYS_GLO){
                ii=i+3*MAX_GNSS_FRQ_NUM;
                if(ii==2) ii=1+3*MAX_GNSS_FRQ_NUM;
            }
//            for(j=0;j<3;j++) off[j]=rec_ant_[rec_idx_].pco[ii][j]+ant_del[rov_bas][j];
            for(j=0;j<3;j++) off[j]=rec_ant_[rec_idx].pco[ii][j]+rec_ant_[rec_idx].rec_ant_del[rec_idx][j];
            if(rec_ant_[rec_idx].dazi!=0.0){
                dantr[i]=-Dot(off,e,3)+InterpAziPcv(rec_ant_[rec_idx],sat_info->el_az[1]*R2D,90-sat_info->el_az[0]*R2D,ii);
            }
            else{
                dantr[i]=-Dot(off,e,3)+InterpPcv(90-sat_info->el_az[1]*R2D,rec_ant_[rec_idx].pcv[ii]);
            }
        }
    }

    cEphModel::cEphModel() {}

    cEphModel::~cEphModel() {}

    int cEphModel::SelectBrdEph(tSatInfoUnit* sat_info,int iode) {
        double t,tmax,tmin;
        int neph=-1;
        int sys=sat_info->sat.sat_.sys;
        switch(sys){
            case SYS_BDS: tmax=MAX_DTOE_BDS+1.0;break;
            case SYS_GAL: tmax=MAX_DTOE_GAL+1.0;break;
            case SYS_GLO: tmax=MAX_DTOE_GLO+1.0;break;
            case SYS_QZS: tmax=MAX_DTOE_QZS+1.0;break;
            default:
                tmax=MAX_DTOE_GPS+1.0;break;
        }
        tmin=tmax+1.0;

        for(int i=0;i<brd_eph_.size();i++){
            if(brd_eph_[i].sat.sat_.no!=sat_info->sat.sat_.no) continue;
            if(iode>=0&&brd_eph_[i].iode!=iode) continue;
            if((t=fabs(brd_eph_[i].toe.TimeDiff(sat_info->t_tag.t_)))>tmax) continue;
            if(iode>0) return i;
            if(t<=tmin){
                neph=i;tmin=t;
            }
        }
        if(iode>=0||neph<0) return -2;

        return neph;
    }

    int cEphModel::SelectGloBrdEph(tSatInfoUnit* sat_info,int iode) {
        double t,tmax=MAX_DTOE_GLO,tmin=tmax+1.0;
        int ngeph=-1;

        for(int i=0;i<brd_glo_eph_.size();i++){
            if(brd_glo_eph_[i].sat.sat_.no!=sat_info->sat.sat_.no) continue;
            if(iode>=0&&brd_glo_eph_[i].iode!=iode) continue;
            if((t=fabs(brd_glo_eph_[i].toe.TimeDiff(sat_info->t_tag.t_)))>tmax) continue;
            if(iode>=0) return i;
            if(t<=tmin) {ngeph=i;tmin=t;}
        }

        if(iode>=0||ngeph<0) return -1;

        return ngeph;
    }


    void cEphModel::BrdSatClkErr(tSatInfoUnit* sat_info,tBrdEphUnit eph) {
        double t;

        t=sat_info->t_trans.TimeDiff(eph.toc.t_);
        for(int i=0;i<2;i++){
            t-=eph.f0+eph.f1*t+eph.f2*t*t;
        }
        sat_info->brd_clk[0]=eph.f0+eph.f1*t+eph.f2*t*t;
    }

    void cEphModel::BrdGloSatClkErr(tSatInfoUnit* sat_info,tBrdGloEphUnit glo_eph) {
        double t;

        t=sat_info->t_trans.TimeDiff(glo_eph.toe.t_);

        for(int i=0;i<2;i++){
            t-=-glo_eph.taun+glo_eph.gamn*t;
        }
        sat_info->brd_clk[0]=-glo_eph.taun+glo_eph.gamn*t;
    }

    bool cEphModel::BrdClkError(tSatInfoUnit* sat_info) {
        int sys=sat_info->sat.sat_.sys,idx_eph=0,idx_glo_eph=0;

        if(sys==SYS_GPS||sys==SYS_BDS||sys==SYS_GAL||sys==SYS_QZS){
            if((sat_info->brd_eph_index=idx_eph=SelectBrdEph(sat_info,-1))<0){
                sat_info->stat=SAT_NO_PROD;
//                return false;
            }
            else BrdSatClkErr(sat_info,brd_eph_[idx_eph]);
        }
        else if(sys==SYS_GLO){
            if((sat_info->brd_eph_index=idx_glo_eph=SelectGloBrdEph(sat_info,-1))<0){
                sat_info->stat=SAT_NO_PROD;
//                return false;
            }
            else BrdGloSatClkErr(sat_info,brd_glo_eph_[idx_glo_eph]);
        }

        return true;
    }

    double cEphModel::EphUra(int ura) {
        const double ura_value[]={
                2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
                3072.0,6144.0
        };
        return ura<0||15<ura?SQR(6144.0):SQR(ura_value[ura]);
    }

    double cEphModel::ClkRelErr(double mu, double sinE,tBrdEphUnit eph) {
        return 2.0*sqrt(mu*eph.A)*eph.e*sinE/SQR(CLIGHT);
    }

    void cEphModel::BrdSatPos(tSatInfoUnit& sat_info,tBrdEphUnit eph) {
        const double RTOL_KEPLER=1E-13;
        const int MAX_ITER_KEPLER=30;
        const double SIN_5=-0.0871557427476582;
        const double COS_5=0.9961946980917456;
        double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
        double xg,yg,zg,sino,coso;
        cSat sat=sat_info.sat;
        int sys=sat.sat_.sys,n;

        if (eph.A<=0.0) {
            return;
        }
        tk=sat_info.t_trans.TimeDiff(eph.toe.t_);

        switch (sys) {
            case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
            case SYS_BDS: mu=MU_BDS; omge=OMGE_BDS; break;
            default:  mu=MU_GPS; omge=OMGE_GPS; break;
        }
        M=eph.M0+(sqrt(mu/(eph.A*eph.A*eph.A))+eph.deln)*tk;

        for (n=0,E=M,Ek=0.0; fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER; n++) {
            Ek=E; E-=(E-eph.e*sin(E)-M)/(1.0-eph.e*cos(E));
        }
        if (n>=MAX_ITER_KEPLER) {
            return;
        }
        sinE=sin(E); cosE=cos(E);

        u=atan2(sqrt(1.0-eph.e*eph.e)*sinE,cosE-eph.e)+eph.omg;
        r=eph.A*(1.0-eph.e*cosE);
        i=eph.i0+eph.idot*tk;
        sin2u=sin(2.0*u); cos2u=cos(2.0*u);
        u+=eph.cus*sin2u+eph.cuc*cos2u;
        r+=eph.crs*sin2u+eph.crc*cos2u;
        i+=eph.cis*sin2u+eph.cic*cos2u;
        x=r*cos(u); y=r*sin(u); cosi=cos(i);

        /* beidou geo satellite (ref [9]) */
        if (sys==SYS_BDS&&sat.sat_.prn<=5) {
            O=eph.Omg0+eph.Omgd*tk-omge*eph.toes;
            sinO=sin(O); cosO=cos(O);
            xg=x*cosO-y*cosi*sinO;
            yg=x*sinO+y*cosi*cosO;
            zg=y*sin(i);
            sino=sin(omge*tk); coso=cos(omge*tk);
            sat_info.brd_pos[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
            sat_info.brd_pos[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
            sat_info.brd_pos[2]=-yg*SIN_5+zg*COS_5;
        }
        else {
            O=eph.Omg0+(eph.Omgd-omge)*tk-omge*eph.toes;
            sinO=sin(O); cosO=cos(O);
            sat_info.brd_pos[0]=x*cosO-y*cosi*sinO;
            sat_info.brd_pos[1]=x*sinO+y*cosi*cosO;
            sat_info.brd_pos[2]=y*sin(i);
        }
        tk=sat_info.t_trans.TimeDiff(eph.toc.t_);
        sat_info.brd_clk[0]=eph.f0+eph.f1*tk+eph.f2*tk*tk;

        // relativity correction
        // the onboad clock will be affected by a relativistic clock correction
        // caused by the motion of the satellite as well as the change in the
        // gravitational potential, the effect due to the orbit eccentricity can
        // be computed by
        sat_info.clk_rel=ClkRelErr(mu,sinE,eph);
        sat_info.brd_clk[0]-=sat_info.clk_rel;

        /* position and clock error variance */
        sat_info.brd_eph_var=EphUra(eph.sva);
    }

    void cEphModel::GloDefEq(const VectorXd x, VectorXd &x_dot, const Vector3d acc) {
        const double J2_GLO=1.0826257E-3;
        const double RE_GLO=6378136.0;
        double a,b,c,r2=x.dot(x),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);

        if (r2<=0.0) {
            x_dot[0]=x_dot[1]=x_dot[2]=x_dot[3]=x_dot[4]=x_dot[5]=0.0;
            return;
        }
        /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
        a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
        b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
        c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
        x_dot[0]=x[3]; x_dot[1]=x[4]; x_dot[2]=x[5];
        x_dot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
        x_dot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
        x_dot[5]=(c-2.0*a)*x[2]+acc[2];
    }

    void cEphModel::GloOrbit(double t, VectorXd& x, const Vector3d acc) {
        VectorXd k1(6),k2(6),k3(6),k4(6),w(6);
        int i;

        GloDefEq(x,k1,acc); for(i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
        GloDefEq(w,k2,acc); for(i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
        GloDefEq(w,k3,acc); for(i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
        GloDefEq(w,k4,acc);
        for(i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
    }

    void cEphModel::BrdGloSatPos(tSatInfoUnit& sat_info,tBrdGloEphUnit glo_eph) {
        const double TSTEP=60.0;
        double t,tt;
        VectorXd x(6);
        int i;

        x<<glo_eph.pos[0],glo_eph.pos[1],glo_eph.pos[2],glo_eph.vel[0],glo_eph.vel[1],glo_eph.vel[2];

        t=sat_info.t_trans.TimeDiff(glo_eph.toe.t_);

        sat_info.brd_clk[0]=-glo_eph.taun+glo_eph.gamn*t;

        for (tt=t<0.0 ? -TSTEP : TSTEP; fabs(t)>1E-9; t-=tt) {
            if (fabs(t)<TSTEP) tt=t;
            GloOrbit(tt,x,glo_eph.acc);
        }
        for (i=0; i<3; i++) sat_info.brd_pos[i]=x[i];

        sat_info.brd_eph_var=SQR(5.0);
    }

    bool cEphModel::
    BrdEphModel(tSatInfoUnit* sat_info) {
        sat_info->svh=-1;
        tSatInfoUnit d1E_3;

        int sys=sat_info->sat.sat_.sys;
        if(sys==SYS_GPS||sys==SYS_BDS||sys==SYS_GAL||sys==SYS_QZS){
            BrdSatPos(*sat_info,brd_eph_[sat_info->brd_eph_index]);
            if(sat_info->brd_eph_var<129.0) sat_info->svh=brd_eph_[sat_info->brd_eph_index].svh;

            d1E_3.t_trans=sat_info->t_trans;d1E_3.t_trans+=1E-3;
            BrdSatPos(d1E_3,brd_eph_[sat_info->brd_eph_index]);
        }
        else if(sys==SYS_GLO){
            BrdGloSatPos(*sat_info,brd_glo_eph_[sat_info->brd_eph_index]);
            sat_info->svh=brd_glo_eph_[sat_info->brd_eph_index].svh;

            d1E_3.t_trans=sat_info->t_trans;d1E_3.t_trans+=1E-3;
            BrdGloSatPos(d1E_3,brd_glo_eph_[sat_info->brd_eph_index]);
        }
        else return false;

        sat_info->brd_vel=(sat_info->brd_pos-d1E_3.brd_pos)/1E-3;
        sat_info->brd_clk[1]=(sat_info->brd_clk[0]-d1E_3.brd_clk[0])/1E-3;

        return true;
    }

    void cEphModel::SatPcoCorr(tSatInfoUnit *sat_info, Vector3d sat_pos, Vector3d& dant) {
        double erp_val[5]={0};
        Vector3d sun_pos={0,0,0},moon_pos={0,0,0};
        const tAntUnit *sat_ant=&sat_ant_[sat_info->sat.sat_.no-1];

        cTime t=sat_info->t_trans.Gpst2Utc();
        SunMoonPos(sat_info->t_trans.Gpst2Utc(),erp_val,sun_pos,moon_pos);

        Vector3d r,ez,es,ey,ex;
        for(int i=0;i<3;i++) r[i]=-sat_pos[i];
        if(!NormV3(r,ez.data())) return;
        for(int i=0;i<3;i++) r[i]=sun_pos[i]-sat_pos[i];
        if(!NormV3(r,es.data())) return;
        CrossVec3(ez,es,r.data());
        if(!NormV3(r,ey.data())) return;
        CrossVec3(ey,ez,ex.data());

        int sys=sat_info->sat.sat_.sys;
        int i=0,j=1;
        //pco[MAX_GNSS_FREQ]
        //GPS L1 L2 L5
        //BDS B1I(2) B2I(7) B3I(6) B1C(1) B2A(5) B2B(7)
        //GAL E1(1) E5a(5) E5b
        //GLO G1 G2
        //QZS L1 L2 L5
        double gamma,C1,C2;
        if(sys==SYS_GPS){
            i=0;j=1;
            gamma=SQR(FREQ_GPS_L1)/SQR(FREQ_GPS_L2);
        }
        else if(sys==SYS_BDS){
            i=0+MAX_GNSS_FRQ_NUM;j=2+MAX_GNSS_FRQ_NUM;
            gamma=SQR(FREQ_BDS_B1)/SQR(FREQ_BDS_B3);
            if(PPPLibC_.gnssC.ac_opt==AC_COM){
                j=1+MAX_GNSS_FRQ_NUM;
                gamma=SQR(FREQ_BDS_B1)/SQR(FREQ_BDS_B2);
            }
        }
        else if(sys==SYS_GAL){
            gamma=SQR(FREQ_GAL_E1)/SQR(FREQ_GAL_E5A);
            i=0+2*MAX_GNSS_FRQ_NUM;j=1+2*MAX_GNSS_FRQ_NUM;
        }
        else if(sys==SYS_GLO){
            gamma=SQR(FREQ_GLO_G1)/SQR(FREQ_GLO_G2);
            i=0+3*MAX_GNSS_FRQ_NUM;j=1+3*MAX_GNSS_FRQ_NUM;
        }
        else if(sys==SYS_QZS){
            gamma=SQR(FREQ_QZS_L1)/SQR(FREQ_QZS_L2);
            i=0+4*MAX_GNSS_FRQ_NUM;j=1+4*MAX_GNSS_FRQ_NUM;
        }
        C1=gamma/(gamma-1.0);C2=-1.0/(gamma-1.0);

        for(int k=0;k<3;k++){
            double dant1=sat_ant->pco[i][0]*ex[k]+sat_ant->pco[i][1]*ey[k]+sat_ant->pco[i][2]*ez[k];
            double dant2=sat_ant->pco[j][0]*ex[k]+sat_ant->pco[j][1]*ey[k]+sat_ant->pco[j][2]*ez[k];
            dant[k]=C1*dant1+C2*dant2;
        }
    }

    const int NMAX=10;
    const double MAXDTE=900.0;
    const double EXTERR_CLK=1E-3;
    const double EXTERR_EPH=5E-7;

    bool cEphModel::PreSatClkCorr(tSatInfoUnit *sat_info) {
        int i,j,k,idx;
        double t[2],c[2],std;

        if(pre_clk_.size()<2||
           sat_info->t_trans.TimeDiff(pre_clk_[0].t_tag.t_)<-MAXDTE||
           sat_info->t_trans.TimeDiff(pre_clk_[pre_clk_.size()-1].t_tag.t_)>MAXDTE){
            return false;
        }
        for(i=0,j=pre_clk_.size()-1;i<j;){
            k=(i+j)/2;
            if(pre_clk_[k].t_tag.TimeDiff(sat_info->t_trans.t_)<0.0) i=k+1;
            else j=k;
        }
        idx=i<=0?0:i-1;

        t[0]=sat_info->t_trans.TimeDiff(pre_clk_[idx].t_tag.t_);
        t[1]=sat_info->t_trans.TimeDiff(pre_clk_[idx+1].t_tag.t_);
        c[0]=pre_clk_[idx].clk[sat_info->sat.sat_.no-1];
        c[1]=pre_clk_[idx+1].clk[sat_info->sat.sat_.no-1];

        if(t[0]<=0.0){
            if((sat_info->pre_clk[0]=c[0])==0.0) return 0;
            std=pre_clk_[idx].std[sat_info->sat.sat_.no-1]*CLIGHT-EXTERR_CLK*t[0];
        }
        else if(t[1]>=0.0){
            if((sat_info->pre_clk[0]=c[1])==0.0) return 0;
            std=pre_clk_[idx+1].std[sat_info->sat.sat_.no-1]*CLIGHT+EXTERR_CLK*t[1];
        }
        else if(c[0]!=0.0&&c[1]!=0.0){
            sat_info->pre_clk[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
            i=t[0]<-t[1]?0:1;
            std=pre_clk_[idx+i].std[sat_info->sat.sat_.no-1];

            if(std*CLIGHT>0.05) std+=EXTERR_CLK*fabs(t[i]);
            else                std=std*CLIGHT+EXTERR_CLK*fabs(t[i]);
        }
        else{
            return false;
        }
        sat_info->pre_eph_var+=SQR(std);

        return true;
    }

    static void EarthRotCorr(int k,Vector3d pos,double p[3][NMAX+1],double dt){
        double sinl,cosl;
        sinl=sin(OMGE_GPS*dt);
        cosl=cos(OMGE_GPS*dt);
        p[0][k]=cosl*pos[0]-sinl*pos[1];
        p[1][k]=sinl*pos[0]+cosl*pos[1];
        p[2][k]=pos[2];
    }

    static double InterpolNev(const double*x ,double* y,int n){
        int i,j;
        for(j=1;j<n;j++){
            for(i=0;i<n-j;i++){
                y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
            }
        }
        return y[0];
    }

    bool cEphModel::PreSatPos(tSatInfoUnit *sat_info) {
        double t[NMAX+1],p[3][NMAX+1],c[2],std=0.0,sinl,cosl;
        Vector3d s;
        Vector4d pos;
        int i,j,k,index;
        int id[NMAX+1],idx;

        sat_info->pre_pos<<0.0,0.0,0.0;
        sat_info->pre_clk<<0.0,0.0;
        if (pre_orb_.size()<NMAX||
            sat_info->t_trans.TimeDiff(pre_orb_[0].t_tag.t_)<-MAXDTE||
            sat_info->t_trans.TimeDiff(pre_orb_[pre_orb_.size()-1].t_tag.t_)>MAXDTE) {
            return false;
        }
        /* binary search */
        for (i=0,j=pre_orb_.size()-1; i<j;) {
            k=(i+j)/2;
            if (pre_orb_[k].t_tag.TimeDiff(sat_info->t_trans.t_)<0.0) i=k+1;
            else j=k;
        }
        index=i<=0?0:i-1;

        /* polynomial interpolation for orbit */
        i=index-(NMAX+1)/2;
        if (i<0) i=0; else if (i+NMAX>pre_orb_.size()) i=pre_orb_.size()-NMAX-1;

        for(j=k=0;j<NMAX*50;j++){
            if(index+j>=0&&index+j<pre_orb_.size()&&k<=NMAX){
                id[k]=index+j;
                t[k]=pre_orb_[id[k]].t_tag.TimeDiff(sat_info->t_trans.t_);
                pos=pre_orb_[id[k]].pos[sat_info->sat.sat_.no-1];
                Vector3d pp(pos[0],pos[1],pos[2]);
                if(pp.norm()>0.0){
                    EarthRotCorr(k,pp,p,t[k]);
                    k++;
                }
            }
            if(k==NMAX+1) break;
            if(index-j>=0&&index-j<pre_orb_.size()&&k<=NMAX&&j!=0){
                id[k]=index-j;
                t[k]=pre_orb_[id[k]].t_tag.TimeDiff(sat_info->t_trans.t_);
                if(t[k]==10){
                    int a=1;
                }
                pos=pre_orb_[id[k]].pos[sat_info->sat.sat_.no-1];
                Vector3d pp(pos[0],pos[1],pos[2]);
                if(pp.norm()>0.0){
                    EarthRotCorr(k,pp,p,t[k]);
                    k++;
                }
            }
            if(k==NMAX+1) break;
        }
        if(k<=NMAX) return false;

        for (i=0;i<NMAX+1;i++) {
            if(i==9){
                int a=1;
            }
            for(j=i+1;j<NMAX+1;j++){
                if(t[i]<=t[j]) continue;
                sinl=t[j];t[j]=t[i];t[i]=sinl;
                k=id[j];id[j]=id[i];id[i]=k;
                for(k=0;k<3;k++){
                    sinl=p[k][j];
                    p[k][j]=p[k][i];
                    p[k][i]=sinl;
                }
            }
        }

        idx=0;
        for(i=0;i<=NMAX;i++){
            if(fabs(t[idx])<=fabs(t[i])) idx=i;
        }
        index=id[idx];
        if(t[0]>900.0||t[NMAX]<-900.0) return false;

        sat_info->pre_pos[0]=InterpolNev(t,p[0],NMAX+1);
        sat_info->pre_pos[1]=InterpolNev(t,p[1],NMAX+1);
        sat_info->pre_pos[2]=InterpolNev(t,p[2],NMAX+1);

        /* satellite position variance */
        for(i=0;i<3;i++) s[i]=pre_orb_[index].std_pos[sat_info->sat.sat_.no-1][i];
        std=s.norm();
        if(t[0]>0.0)         std+=EXTERR_EPH*SQR(t[0])/2.0;
        else if(t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        sat_info->pre_eph_var=SQR(std);
        return true;
    }

    bool cEphModel::PreEphModel(tSatInfoUnit *sat_info) {
        tSatInfoUnit d1E_3=*sat_info;
        d1E_3.t_trans+=(1E-3);

        if(!PreSatPos(sat_info)||!PreSatClkCorr(sat_info)) return false;
        if(!PreSatPos(&d1E_3)||!PreSatClkCorr(&d1E_3)) return false;

        sat_info->pre_vel=(d1E_3.pre_pos-sat_info->pre_pos)/1E-3;

        Vector3d dant(0,0,0);
        SatPcoCorr(sat_info,sat_info->pre_pos,dant);
        sat_info->pre_pos+=dant;


        if(sat_info->pre_clk[0]!=0.0){
            sat_info->pre_clk[0]=sat_info->pre_clk[0]-2.0*sat_info->pre_pos.dot(sat_info->pre_vel)/CLIGHT/CLIGHT;
            sat_info->pre_clk[1]=(d1E_3.pre_clk[0]-sat_info->pre_clk[0])/1E-3;
        }

        return true;
    }

    void cEphModel::InitEphModel(vector<tBrdEphUnit> &brd_eph, vector<tBrdGloEphUnit> &brd_glo_eph,
                                 vector<tPreOrbUnit> &pre_orb, vector<tPreClkUnit> &pre_clk,tAntUnit* sat_ant) {
        if(brd_eph.size()>0) brd_eph_=brd_eph;
        if(brd_glo_eph.size()>0) brd_glo_eph_=brd_glo_eph;
        if(pre_orb.size()>0) pre_orb_=pre_orb;
        if(pre_clk.size()>0) pre_clk_=pre_clk;
        if(sat_ant) sat_ant_=sat_ant;
    }

    bool cEphModel::EphCorr(vector<tSatInfoUnit>& epoch_sat_info) {
        int i;
        tSatInfoUnit* sat_info= nullptr;
        char buff[MAX_BUFF]={'\0'};

        for(int j=0;j<epoch_sat_info.size();j++){
            sat_info=&epoch_sat_info.at(j);
            // signal transmit time
            double pr;
            sat_info->t_trans=sat_info->t_tag;
            for(i=0,pr=0.0;i<MAX_GNSS_USED_FRQ_NUM;i++) if((pr=sat_info->raw_P[i])!=0.0) break;
            if(i<MAX_GNSS_USED_FRQ_NUM) sat_info->t_trans+=(-pr/CLIGHT);

            // broadcast clock error
            if(!BrdClkError(sat_info)) return false;

            if(sat_info->stat!=SAT_USED) continue;
            // broadcast clock error correction
            sat_info->t_trans+=(-sat_info->brd_clk[0]);

            if(PPPLibC_.gnssC.eph_opt==EPH_BRD){
                LOG_IF(j==0,DEBUG)<<"SATELLITE BROADCAST EPHEMERIS: "<<sat_info->t_tag.GetTimeStr(1);
                if(!BrdEphModel(sat_info)){
                    sat_info->stat=SAT_NO_USE;
                    continue;
                }
                sprintf(buff,"%s %14.3f %14.3f %14.3f %14.3f %10.3f %10.3f %10.3f",sat_info->sat.sat_.id.c_str(),
                        sat_info->brd_pos[0],sat_info->brd_pos[1],sat_info->brd_pos[2],sat_info->brd_clk[0]*CLIGHT,
                        sat_info->brd_vel[0],sat_info->brd_vel[1],sat_info->brd_vel[2]);
                LOG(DEBUG)<<buff;
            }
            else if(PPPLibC_.gnssC.eph_opt==EPH_PRE){
                LOG_IF(j==0,DEBUG)<<"SATELLITE PRECISE EPHEMERIS: "<<sat_info->t_tag.GetTimeStr(1);
                if(!PreEphModel(sat_info)){
                    sat_info->stat=SAT_NO_USE;
                    continue;
                }
                sprintf(buff,"%s %14.3f %14.3f %14.3f %14.3f %10.3f %10.3f %10.3f",sat_info->sat.sat_.id.c_str(),
                        sat_info->pre_pos[0],sat_info->pre_pos[1],sat_info->pre_pos[2],sat_info->pre_clk[0]*CLIGHT,
                        sat_info->pre_vel[0],sat_info->pre_vel[1],sat_info->pre_vel[2]);
                LOG(DEBUG)<<buff;
            }
        }
    }

    cTidModel::cTidModel() {}

    cTidModel::~cTidModel() {}

    void cTidModel::InitTidModel(tPPPLibConf C, vector<tErpUnit> erp, double ocean_par[2][66]) {
        PPPLibC_=C;
        erp_=erp;
        for(int i=0;i<2;i++){
            for(int j=0;j<66;j++){
                ocean_par_[i][j]=ocean_par[i][j];
            }
        }
    }

    int cTidModel::GetErpVal(cTime time) {
        const double ep[]={2000,1,1,12,0,0};
        const cTime t(ep);
        double mjd,day,a;
        int i,j,k;

        if(erp_.size()<=0) return 0;

        mjd=51544.5+time.Gpst2Utc().TimeDiff(t.t_)/86400.0;

        if(mjd<=erp_[0].mjd){
            day=mjd-erp_[0].mjd;
            erp_val_[0]=erp_[0].xp+erp_[0].xpr*day;
            erp_val_[1]=erp_[0].yp+erp_[0].ypr*day;
            erp_val_[2]=erp_[0].ut1_utc-erp_[0].lod*day;
            erp_val_[3]=erp_[0].lod;
            return 1;
        }
        int n=erp_.size()-1;
        if(mjd>=erp_[n].mjd){
            day=mjd-erp_[0].mjd;
            erp_val_[0]=erp_[n].xp+erp_[n].xpr*day;
            erp_val_[1]=erp_[n].yp+erp_[n].ypr*day;
            erp_val_[2]=erp_[n].ut1_utc-erp_[n].lod*day;
            erp_val_[3]=erp_[n].lod;
            return 1;
        }
        for(j=0,k=n;j<k-1;){
            i=(j+k)/2;
            if(mjd<erp_[i].mjd) k=i;else j=i;
        }
        if(erp_[j].mjd==erp_[j+1].mjd) a=0.5;
        else{
            a=(mjd-erp_[j].mjd)/(erp_[j+1].mjd-erp_[j].mjd);
        }
        erp_val_[0]=(1.0-a)*erp_[j].xp+a*erp_[j+1].xp;
        erp_val_[1]=(1.0-a)*erp_[j].yp+a*erp_[j+1].yp;
        erp_val_[2]=(1.0-a)*erp_[j].ut1_utc+a*erp_[j+1].ut1_utc;
        erp_val_[3]=(1.0-a)*erp_[j].lod+a*erp_[j+1].lod;
        return 1;
    }

    #define GME         3.986004415E+14 /* earth gravitational constant */
    void cTidModel::TideSl(const double* eu,const double* rp,double GMp,Vector3d blh,double* dr) {
        const double H3=0.292,L3=0.015;
        double r,ep[3],latp,lonp,p,K2,K3,a,H2,L2,dp,du,cosp,sinl,cosl;
        int i;

        if ((r=Norm(rp,3))<=0.0) return;

        for (i=0; i<3; i++) ep[i]=rp[i]/r;

        K2=GMp/GME*SQR(WGS84_EARTH_LONG_RADIUS)*SQR(WGS84_EARTH_LONG_RADIUS)/(r*r*r);
        K3=K2*WGS84_EARTH_LONG_RADIUS/r;
        latp=asin(ep[2]); lonp=atan2(ep[1],ep[0]);
        cosp=cos(latp); sinl=sin(blh[0]); cosl=cos(blh[0]);

        /* step1 in phase (degree 2) */
        p=(3.0*sinl*sinl-1.0)/2.0;
        H2=0.6078-0.0006*p;
        L2=0.0847+0.0002*p;
        a=Dot(ep,eu,3);
        dp=K2*3.0*L2*a;
        du=K2*(H2*(1.5*a*a-0.5)-3.0*L2*a*a);

        /* step1 in phase (degree 3) */
        dp+=K3*L3*(7.5*a*a-1.5);
        du+=K3*(H3*(2.5*a*a*a-1.5*a)-L3*(7.5*a*a-1.5)*a);

        /* step1 out-of-phase (only radial) */
        du+=3.0/4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*blh[0])*sin(blh[1]-lonp);
        du+=3.0/4.0*0.0022*K2*cosp*cosp*cosl*cosl*sin(2.0*(blh[1]-lonp));

        dr[0]=dp*ep[0]+du*eu[0];
        dr[1]=dp*ep[1]+du*eu[1];
        dr[2]=dp*ep[2]+du*eu[2];
    }

    #define GMS         1.327124E+20    /* sun gravitational constant */
    #define GMM         4.902801E+12    /* moon gravitational constant */
    void cTidModel::TidSolid() {
        double dr1[3],dr2[3],du,dn,sinl,sin2l;
        Vector3d eu;

        /* step1: time domain */
        eu[0]=E_(2,0); eu[1]=E_(2,1); eu[2]=E_(2,2);

        /* step2: frequency domain, only K1 radial */
        sin2l=sin(2.0*blh_[0]);
        du=-0.012*sin2l*sin(gmst_+blh_[1]);
        TideSl(eu.data(),sun_pos_.data(),GMS,blh_,dr1);
        TideSl(eu.data(),moon_pos_.data(),GMM,blh_,dr2);

        tid_dr_[0][0]=dr1[0]+dr2[0]+du*E_(2,0);
        tid_dr_[0][1]=dr1[1]+dr2[1]+du*E_(2,1);
        tid_dr_[0][2]=dr1[2]+dr2[2]+du*E_(2,2);

        /* eliminate permanent deformation */
        sinl=sin(blh_[0]);
        du=0.1196*(1.5*sinl*sinl-0.5);
        dn=0.0247*sin2l;
        tid_dr_[0][0]+=du*E_(2,0)+dn*E_(1,0);
        tid_dr_[0][1]+=du*E_(2,1)+dn*E_(1,1);
        tid_dr_[0][2]+=du*E_(2,2)+dn*E_(1,2);
    }

    void cTidModel::TidOcean() {
        if (!ocean_par_[rec_idx_]) return ;
        const double args[][5]={
                { 1.40519E-4, 2.0,-2.0, 0.0, 0.00 },  /* M2 */
                { 1.45444E-4, 0.0, 0.0, 0.0, 0.00 },  /* S2 */
                { 1.37880E-4, 2.0,-3.0, 1.0, 0.00 },  /* N2 */
                { 1.45842E-4, 2.0, 0.0, 0.0, 0.00 },  /* K2 */
                { 0.72921E-4, 1.0, 0.0, 0.0, 0.25 },  /* K1 */
                { 0.67598E-4, 1.0,-2.0, 0.0,-0.25 },  /* O1 */
                { 0.72523E-4,-1.0, 0.0, 0.0,-0.25 },  /* P1 */
                { 0.64959E-4, 1.0,-3.0, 1.0,-0.25 },  /* Q1 */
                { 0.53234E-5, 0.0, 2.0, 0.0, 0.00 },  /* Mf */
                { 0.26392E-5, 0.0, 1.0,-1.0, 0.00 },  /* Mm */
                { 0.03982E-5, 2.0, 0.0, 0.0, 0.00 }   /* Ssa */
        };
        const double ep1975[]={ 1975,1,1,0,0,0 };
        double fday,days,t,t2,t3,a[5],ang,dp[3]={ 0 };
        int i,j;

        /* angular argument: see subroutine arg.f for reference [1] */
        cTime time=tut_,t1975(ep1975);
        time.Time2Epoch();
        fday=time.GetEpoch()[3]*3600.0+time.GetEpoch()[4]*60.0+time.GetEpoch()[5];
        time.GetEpoch()[3]=time.GetEpoch()[4]=time.GetEpoch()[5]=0.0;
        days=(time.Epoch2Time(time.GetEpoch())->TimeDiff(t1975.t_))/86400.0+1.0;
        t=(27392.500528+1.000000035*days)/36525.0;
        t2=t*t; t3=t2*t;

        a[0]=fday;
        a[1]=(279.69668+36000.768930485*t+3.03E-4*t2)*D2R; /* H0 */
        a[2]=(270.434358+481267.88314137*t-0.001133*t2+1.9E-6*t3)*D2R; /* S0 */
        a[3]=(334.329653+4069.0340329577*t-0.010325*t2-1.2E-5*t3)*D2R; /* P0 */
        a[4]=2.0*PI;

        /* displacements by 11 constituents */
        for (i=0; i<11; i++) {
            ang=0.0;
            for (j=0; j<5; j++) ang+=a[j]*args[i][j];
            for (j=0; j<3; j++) dp[j]+=ocean_par_[rec_idx_][j+i*6]*cos(ang-ocean_par_[rec_idx_][j+3+i*6]*D2R);
        }
        denu_[0][0]=-dp[1];
        denu_[0][1]=-dp[2];
        denu_[0][2]= dp[0];
    }

    void cTidModel::IersMeanPole(double *xp,double *yp) {
        const double ep2000[]={ 2000,1,1,0,0,0 };
        double y,y2,y3;

        cTime t2000(ep2000);
        y=tut_.TimeDiff(t2000.t_)/86400.0/365.25;

        if (y<3653.0/365.25) { /* until 2010.0 */
            y2=y*y; y3=y2*y;
            *xp= 55.974+1.8243*y+0.18413*y2+0.007024*y3; /* (mas) */
            *yp=346.346+1.7896*y-0.10729*y2-0.000908*y3;
        }
        else { /* after 2010.0 */
            *xp= 23.513+7.6141*y; /* (mas) */
            *yp=358.891-0.6287*y;
        }
    }

    void cTidModel::TidPole() {
        double xp_bar,yp_bar,m1,m2,cosl,sinl;

        /* iers mean pole (mas) */
        IersMeanPole(&xp_bar,&yp_bar);

        /* ref [7] eq.7.24 */
        m1= erp_val_[0]/AS2R-xp_bar*1E-3; /* (as) */
        m2=-erp_val_[1]/AS2R+yp_bar*1E-3;

        /* sin(2*theta) = sin(2*phi), cos(2*theta)=-cos(2*phi) */
        cosl=cos(blh_[1]);
        sinl=sin(blh_[1]);
        denu_[1][0]=  9E-3*cos(blh_[0])    *(m1*sinl-m2*cosl); /* de= Slambda (m) */
        denu_[1][1]= -9E-3*cos(2.0*blh_[0])*(m1*cosl+m2*sinl); /* dn=-Stheta  (m) */
        denu_[1][2]=-32E-3*sin(2.0*blh_[0])*(m1*cosl+m2*sinl); /* du= Sr      (m) */
    }

    void cTidModel::TidCorr(cTime t, Vector3d rec_xyz, Vector3d& dr) {
        dr<<0,0,0;
        if(erp_.size()) GetErpVal(t);
        tut_=t;tut_+=erp_val_[2];
        if((re_=rec_xyz.norm())<=0.0) return;

        blh_[0]=asin(rec_xyz[2]/re_);
        blh_[1]=atan2(rec_xyz[1],rec_xyz[0]);
        blh_[2]=0.0;

        //TODO check;
        Vector3d rec_blh=Xyz2Blh(rec_xyz);
        E_=CalcCen(blh_,COORD_ENU);//Cen

        if(PPPLibC_.gnssC.tid_opt&TID_SOLID){
            gmst_=SunMoonPos(t,erp_val_,sun_pos_,moon_pos_);
            TidSolid();
            dr+=tid_dr_[0];

            TidOcean();
            tid_dr_[1]=E_.transpose()*denu_[0];
            dr+=tid_dr_[1];

            TidPole();
            tid_dr_[2]=E_.transpose()*denu_[1];
            dr+=tid_dr_[2];
        }
        else if(PPPLibC_.gnssC.tid_opt&TID_OCEAN){

        }
        else if(PPPLibC_.gnssC.tid_opt&TID_POLE){

        }
    }

    cGnssErrCorr::cGnssErrCorr() {}

    cGnssErrCorr::~cGnssErrCorr() {}

    void cGnssErrCorr::InitGnssErrCorr(tPPPLibConf C, tNav* nav) {
        eph_model_.PPPLibC_=cbia_model_.PPPLibC_=ion_model_.PPPLibC_=trp_model_.PPPLibC_=ant_model_.PPPLibC_=tid_model_.PPPLibC_=C;
        eph_model_.InitEphModel(nav->brd_eph,nav->brd_glo_eph,nav->pre_eph,nav->pre_clk,nav->sat_ant);
        cbia_model_.InitCbiasModel(nav->code_bias,nav->brd_eph);
        ion_model_.InitIondelayModel(nav->ion_para);
        ant_model_.InitAntModel(C,nav->sat_ant,nav->rec_ant);
        tid_model_.InitTidModel(C,nav->erp_paras,nav->ocean_paras);
    }


    void cGnssErrCorr::BD2MultipathModel(tSatInfoUnit* sat_info) {
        const static double kBDIGSOCoef[3][10]={
                {-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},	//B1
                {-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},	//B2
                {-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32},	//B3
        };
        const static double kBDMEOCoef[3][10]={
                {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},	//B1
                {-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},	//B2
                {-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47},	//B3
        };

        double el=sat_info->el_az[0]*R2D*0.1;
        int int_el=(int)el;
        int prn=sat_info->sat.sat_.prn;
        int BD2_IGSO,BD2_MEO;
        double mp[3]={0};

        BD2_IGSO=std::binary_search(kBD2_IGSO,kBD2_IGSO+NUM_BD2_IGSO,prn);
        BD2_MEO=std::binary_search(kBD2_MEO,kBD2_MEO+NUM_BD2_MEO,prn);

        if(!BD2_IGSO&&!BD2_MEO) return;

        if(BD2_IGSO){
            if(el<0) for(int i=0;i<3;i++) mp[i]=kBDIGSOCoef[i][0];
            else if(el>=9) for(int i=0;i<3;i++) mp[i]=kBDIGSOCoef[i][9];
            else for(int i=0;i<3;i++) mp[i]=kBDIGSOCoef[i][int_el]*(1.0-el+int_el)+kBDIGSOCoef[i][int_el+1]*(el-int_el);
        }
        if(BD2_MEO){
            if(el<0) for(int i=0;i<3;i++) mp[i]=kBDMEOCoef[i][0];
            else if(el>=9) for(int i=0;i<3;i++) mp[i]=kBDMEOCoef[i][9];
            else for(int i=0;i<3;i++) mp[i]=kBDMEOCoef[i][int_el]*(1.0-el+int_el)+kBDMEOCoef[i][int_el+1]*(el-int_el);
        }

        for(int i=0;i<3;i++){
            if(sat_info->frq[i]==0.0) continue;
            if(sat_info->frq[i]==FREQ_BDS_B1) sat_info->bd2_mp[i]=mp[i];
            if(sat_info->frq[i]==FREQ_BDS_B2) sat_info->bd2_mp[i]=mp[i];
            if(sat_info->frq[i]==FREQ_BDS_B3) sat_info->bd2_mp[i]=mp[i];
        }

        for(int i=0;i<3;i++) bd2_mp_[i]=mp[i];
    }

    double cGnssErrCorr::SagnacCorr(Vector3d sat_xyz,Vector3d rec_xyz){
        return OMGE_GPS*(sat_xyz[0]*rec_xyz[1]-sat_xyz[1]*rec_xyz[0])/CLIGHT;
    }

    double cGnssErrCorr::ShapiroCorr(int sys, Vector3d sat_xyz, Vector3d rec_xyz) {
        double rr=rec_xyz.norm();
        double rs=sat_xyz.norm();
        double l=(sat_xyz-rec_xyz).norm();
        double mu=0.0;

        switch(sys){
            case SYS_BDS: mu=MU_BDS;break;
            case SYS_GAL: mu=MU_GAL;break;
            case SYS_GLO: mu=MU_GLO;break;
            default: mu=MU_GPS;
                break;
        }
        return -2.0*mu/SQR(CLIGHT)*log((rr+rs+l)/(rr+rs-l));
    }

    bool cGnssErrCorr::SatYaw(tSatInfoUnit &sat_info,Vector3d& exs,Vector3d& eys) {
        Vector3d sun_pos(0,0,0),moon_pos;
        double erp[5]={0};
        cTime t=sat_info.t_tag;
        SunMoonPos(t.Gpst2Utc(),erp,sun_pos,moon_pos);

        Vector3d sat_pos,sat_vel,n,p;
        sat_pos=sat_info.pre_pos;sat_vel=sat_info.pre_vel;
        if(sat_vel.norm()==0.0) sat_vel=sat_info.brd_vel;
        sat_vel[0]-=OMGE_GPS*sat_vel[1];sat_vel[1]+=OMGE_GPS*sat_vel[0];
        CrossVec3(sat_pos,sat_vel,n.data());
        CrossVec3(sun_pos,n,p.data());

        Vector3d es,esun,en,ep;
        if(!NormV3(sat_pos,es.data())||(!NormV3(sun_pos,esun.data()))||(!NormV3(n,en.data()))||(!NormV3(p,ep.data()))) return false;
        double beta=PI/2.0-acos(esun.dot(en));
        double E=acos(es.dot(ep));
        double mu=PI/2.0+(es.dot(esun)<=0?-E:E);
        if(mu<-PI/2.0) mu+=2.0*PI;
        else if(mu>=PI/2.0) mu-=2.0*PI;

        double yaw=0.0;
        if(fabs(beta)<1E-12&&fabs(mu)<1E-12) yaw=PI;
        else yaw=atan2(-tan(beta),sin(mu))+PI;

        Vector3d ex;
        CrossVec3(en,es,ex.data());
        double cosy=cos(yaw);
        double siny=sin(yaw);
        for(int i=0;i<3;i++){
            exs[i]=-siny*en[i]+cosy*ex[i];
            eys[i]=-cosy*en[i]-siny*ex[i];
        }
    }

    void cGnssErrCorr::PhaseWindUp(tSatInfoUnit &sat_info, Vector3d rec_xyz) {
        Vector3d exs,eys;
        SatYaw(sat_info,exs,eys);

        Vector3d r=rec_xyz-sat_info.pre_pos;
        Vector3d ek;
        if(!NormV3(r,ek.data())) return;

        Vector3d rec_blh=Xyz2Blh(rec_xyz);
        Matrix3d E;
        E=CalcCen(rec_blh,COORD_ENU);
        Vector3d exr,eyr;
        exr[0]=E(1,0);exr[1]=E(1,1);exr[2]=E(1,2);
        eyr[0]=-E(0,0);eyr[1]=-E(0,1);eyr[2]=-E(0,2);

        Vector3d eks,ekr;
        CrossVec3(ek,eys,eks.data());
        CrossVec3(ek,eyr,ekr.data());

        Vector3d ds,dr;
        for(int i=0;i<3;i++){
            ds[i]=exs[i]-ek[i]*ek.dot(exs)-eks[i];
            dr[i]=exr[i]-ek[i]*ek.dot(exr)+ekr[i];
        }

        double a=Norm(ds,3);
        double b=Norm(dr,3);
        double c=ds.dot(dr);
        double cosp=ds.dot(dr)/(Norm(ds,3))/(Norm(dr,3));
        if(cosp<-1.0) cosp=-1.0;
        double ph=acos(cosp)/2.0/PI;
        Vector3d drs;
        CrossVec3(ds,dr,drs.data());
        if(isnan(ph)) {
            ph=0.0;
            LOG(WARNING)<<sat_info.t_tag.GetTimeStr(1)<<" "<<sat_info.sat.sat_.id<<" phase wind up error";
        }
        if(ek.dot(drs)<0.0) ph=-ph;
        sat_info.phase_wp=ph+floor(sat_info.phase_wp-ph+0.5);

    }

}