//
// Created by cc on 7/16/20.
//

#include "InsFunc.h"

extern const double Crf[9]={0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0};

namespace PPPLib {
    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec){
        Eigen::Matrix3d dcm=Matrix3d::Zero();

        dcm(0,1)=-vec(2);
        dcm(0,2)=vec(1);

        dcm(1,0)=vec(2);
        dcm(1,2)=-vec(0);

        dcm(2,0)=-vec(1);
        dcm(2,1)=vec(0);
        return dcm;
    }

    Eigen::Matrix3d Quaternion2RotationMatrix(const Eigen::Quaterniond& q){
        return q.toRotationMatrix();
    }

    Eigen::Quaterniond RotationMatrix2Quaternion(const Eigen::Matrix3d& m){
        Eigen::Quaterniond q(m);
        return q;
    }

    Eigen::Quaterniond Euler2Quaternion(const Vector3d& rpy){
        Eigen::AngleAxisd roll_angle(Eigen::AngleAxisd(rpy(2),Eigen::Vector3d::UnitZ()));
        Eigen::AngleAxisd pitch_angle(Eigen::AngleAxisd(rpy(1),Eigen::Vector3d::UnitY()));
        Eigen::AngleAxisd yaw_angle(Eigen::AngleAxisd(rpy(0),Eigen::Vector3d::UnitX()));

        return roll_angle*pitch_angle*yaw_angle;
    }

    Eigen::Matrix3d Euler2RotationMatrix(const Vector3d& rpy){
        return Quaternion2RotationMatrix(Euler2Quaternion(rpy));
    }

    template <typename Derived>
    static Eigen::Matrix<Derived,3,1> EulerAngles(const Eigen::Matrix<Derived,3,3>& m){
        Eigen::Matrix<Derived, 3, 1> res;

        const size_t i = 2;
        const size_t j = 1;
        const size_t k = 2;
        typedef Eigen::Matrix<Derived, 2, 1> Vector2;

        res[2]=atan2(m(j,k),m(k,k));
        Derived c2=Vector2(m(i,i),m(i,j)).norm();
//        if(res[2]<Derived(0))
//        {
//            res[2]+=Derived(EIGEN_PI)*2;
//        }
        res[1]=atan2(-m(i,k),c2);
        Derived s1=sin(res[2]);
        Derived c1=cos(res[2]);
        res[0]=atan2(s1*m(k,i)-c1*m(j,i),c1*m(j,j)-s1*m(k,j));

        return res;
    }

    Eigen::Vector3d RotationMatrix2Euler(const Matrix3d &m){
//        return m.eulerAngles(2,1,0);
//        return EulerAngles(m);
        Vector3d rpy;
        rpy[0]=atan2(m(1,2),m(2,2));
        rpy[1]=-asin(m(0,2));
        rpy[2]=atan2(m(0,1),m(0,0));
        return rpy;
    }

    Eigen::Vector3d Quaternion2Euler(const Quaterniond& q){
        return RotationMatrix2Euler(q.toRotationMatrix());
    }

    Eigen::Quaterniond RotationVector2Quaternion(const Vector3d& rv){
        Eigen::Quaterniond qfromrv;
        Vector3d rv_2=rv*0.5;
        double norm=rv_2.norm();
        qfromrv.w()=cos(norm);
        qfromrv.vec()=norm<1E-8?rv_2:(sin(norm)/norm)*rv_2;
        return qfromrv;
    }

    cImuData::cImuData(){}

    cImuData::cImuData(PPPLib::cTime *ts, PPPLib::cTime *te){
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cImuData::~cImuData() {data_.clear();}

    void cImuData::SetImu(tInsConf C){
        imu_type_=C.imu_type;
        imu_coord_type_=C.coord_type;
        data_format_=C.data_format;
        gyro_format_=C.gyro_val_format;
        hz_=C.sample_rate;
    }

    void cImuData::SetImuType(PPPLib::IMU_TYPE type) {imu_type_=type;}

    void cImuData::SetImuCoordType(PPPLib::IMU_COORD_TYPE type) {imu_coord_type_=type;}

    void cImuData::SetTimeSpan(cTime *ts, cTime *te) {
        if(ts) ts_=*ts;
        if(te) te_=*te;
    }

    cInsMech::cInsMech() {}

    cInsMech::~cInsMech() {}

    void cInsMech::RotScullCorr(PPPLib::tImuInfoUnit &pre_imu_info, PPPLib::tImuInfoUnit &cur_imu_info, double dt,double *da,double *dv) {
        Vector3d pre_da,pre_dv,cur_da,cur_dv;
        int i;

        pre_da=pre_imu_info.cor_gyro*cur_imu_info.dt;
        pre_dv=pre_imu_info.cor_acce*cur_imu_info.dt;

        cur_da=cur_imu_info.cor_gyro*dt;
        cur_dv=cur_imu_info.cor_acce*dt;

        Vector3d t1,t2,t3,t4;
        CrossVec3(cur_da.data(),cur_dv.data(),t1.data());
        CrossVec3(pre_da.data(),cur_dv.data(),t2.data());
        CrossVec3(pre_dv.data(),cur_da.data(),t3.data());
        CrossVec3(cur_da.data(),t1.data(),t4.data());

        double a1,a2;
        double b=cur_da.norm();
        if(fabs(b)>1E-6){
            a1=(1.0-cos(b))/SQR(b);
            a2=1.0/SQR(b)*(1.0-sin(b)/b);
        }
        else{
            a1=0.5-SQR(b)/24.0+SQR(SQR(b))/720.0;
            a2=1.0/6.0-SQR(b)/120.0+SQR(SQR(b))/5040.0;
        }

        for(i=0;i<3&&dv;i++){
            dv[i]=a1*t1[i]+a2*t4[i]+1.0/12.0*(t2[i]+t3[i]);
        }

        if(da){
            CrossVec3(pre_da,cur_da,da);
            for(i=0;i<3;i++) da[i]*=1.0/12.0;
        }
    }

    Eigen::Quaterniond cInsMech::AttitudeUpdate(PPPLib::tImuInfoUnit &pre_imu_info,
                                                PPPLib::tImuInfoUnit &cur_imu_info,double dt,Vector3d da) {
        Vector3d theta_k(cur_imu_info.cor_gyro*dt),theta_k_1(pre_imu_info.cor_gyro*dt);

        //等效旋转矢量
        Vector3d cur_phi=theta_k+da; //单子样+前一周期
        Quaterniond quat_bb=RotationVector2Quaternion(cur_phi);

        Vector3d wiee(0,0,-OMGE_GPS);
        Vector3d zeta=wiee*dt;
        Quaterniond quat_ee=RotationVector2Quaternion(zeta);
        Quaterniond quat_k_1=RotationMatrix2Quaternion(pre_imu_info.Cbe.transpose()).conjugate();
        Quaterniond qbn_k=quat_ee*quat_k_1*quat_bb;

        return qbn_k.normalized();
    }

    Eigen::Vector3d cInsMech::VelocityUpdate(PPPLib::tImuInfoUnit &pre_imu_info,
                                             PPPLib::tImuInfoUnit &cur_imu_info,double dt,Vector3d dv) {
        Vector3d pos=pre_imu_info.re, vel=pre_imu_info.ve;
        Vector3d wiee(0,0,OMGE_GPS);
        Vector3d theta_k(cur_imu_info.cor_gyro*dt),theta_k_1(pre_imu_info.cor_gyro*dt);
        Vector3d vb_k(cur_imu_info.cor_acce*dt+dv),vb_k_1(pre_imu_info.cor_acce*dt);

        Vector3d coord_blh=Xyz2Blh(pos);
        Matrix3d Cen=CalcCen(coord_blh,COORD_NED);
        Vector3d ge=Cen.transpose()*CalculateGravity(coord_blh,false);
        Vector3d omgea_n=wiee*2.0;
        Vector3d delta_gcor=(ge-omgea_n.cross(vel))*dt;

        Matrix3d Cee=Matrix3d::Identity()-VectorSkew(wiee*0.5*dt);

//        Vector3d vrot=theta_k.cross(vb_k)*0.5;
//        Vector3d vscul=(theta_k_1.cross(vb_k)+vb_k_1.cross(theta_k))/12.0;

        Quaterniond pre_quat=RotationMatrix2Quaternion(pre_imu_info.Cbe);
        Matrix3d Cbe=pre_quat.toRotationMatrix();
        Vector3d delta_ve=Cee*Cbe*(vb_k);

        return (vel+delta_gcor+delta_ve);
    }

    Eigen::Vector3d cInsMech::PositionUpdate(const tImuInfoUnit& pre_imu_info,const Vector3d& cur_vel, double dt) {
        Eigen::Vector3d pos;
        pos=(pre_imu_info.ve+cur_vel)*0.5*dt+pre_imu_info.re;
        return pos;
    }

    bool cInsMech::InsMechanization(bool err_model,tImuInfoUnit &pre_imu_info, tImuInfoUnit &cur_imu_info,int idx) {

        for(int i=0;i<3;i++){
            if(isnan(pre_imu_info.re[i])||isnan(pre_imu_info.ve[i])||
               isinf(pre_imu_info.re[i])||isinf(pre_imu_info.ve[i])){
                LOG(ERROR)<<cur_imu_info.t_tag.GetTimeStr(4)<<" "<<idx<<" NUMERIC ERROR";
                return false;
            }
        }
        TraceInsMechInfo(pre_imu_info,true,idx);

        double dt=cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_);
        if(dt>60.0||fabs(dt)<1E-6){
            cur_imu_info.dt=dt;
            cur_imu_info.pt=pre_imu_info.pt;
            LOG(WARNING)<<"TIME DIFFERENCE TOO LARGER";
            return false;
        }

        if(err_model){

        }else{
            cur_imu_info.cor_gyro=cur_imu_info.raw_gyro-pre_imu_info.bg;
            cur_imu_info.cor_acce=cur_imu_info.raw_acce-pre_imu_info.ba;
        }

        Vector3d da,dv;
        RotScullCorr(pre_imu_info,cur_imu_info,dt,da.data(),dv.data());
        cur_imu_info.Cbe=Quaternion2RotationMatrix(AttitudeUpdate(pre_imu_info,cur_imu_info,dt,da));
        cur_imu_info.ve=VelocityUpdate(pre_imu_info,cur_imu_info,dt,dv);
        cur_imu_info.re=PositionUpdate(pre_imu_info,cur_imu_info.ve,cur_imu_info.t_tag.TimeDiff(pre_imu_info.t_tag.t_));

        Vector3d blh=Xyz2Blh(cur_imu_info.re);
        Matrix3d Cen=CalcCen(blh,COORD_NED);
        Matrix3d Cnb=cur_imu_info.Cbe.transpose()*Cen.transpose();
        cur_imu_info.rpy=RotationMatrix2Euler(Cnb);
        cur_imu_info.vn=Cen*cur_imu_info.ve;
        cur_imu_info.rn=Cen*cur_imu_info.re;

        TraceInsMechInfo(cur_imu_info,false,idx);
        cur_imu_info.dt=dt;
        cur_imu_info.pt=pre_imu_info.t_tag;

        return true;
    }

    Eigen::MatrixXd cInsMech::StateTransferMat(tPPPLibConf C,PPPLib::tImuInfoUnit &pre_imu_info,
                                               PPPLib::tImuInfoUnit &cur_imu_info,int nx,double dt) {
        using Eigen::Matrix3d;
        using Eigen::MatrixXd;
        using Eigen::Vector3d;

        auto &vel=cur_imu_info.ve;
        auto &fb=cur_imu_info.cor_acce;
        auto &wb=cur_imu_info.cor_gyro;
        Vector3d wiee(0,0,OMGE_GPS);
        auto &Cbe=cur_imu_info.Cbe;

        MatrixXd F=MatrixXd::Zero(nx,nx);

        int ip=0;
        int iv=3;
        int ia=6;
        int iba=9;
        int ibg=12;

        //position-velocity
        F.block<3,3>(ip,iv)=Matrix3d::Identity();

        //velocity-velocity
        F.block<3,3>(iv,iv)=(-2.0*VectorSkew(wiee));
        //velocity-attitude
        F.block<3,3>(iv,ia)=-VectorSkew(Cbe*fb);
        //velocity-ba
        F.block<3,3>(iv,iba)=Cbe;

        //attitude-attitude
        F.block<3,3>(ia,ia)=-1.0*VectorSkew(wiee);
        //attitute-bg
        F.block<3,3>(ia,ibg)=Cbe;

        //ba-ba
        F.block<3,3>(iba,iba)=Matrix3d::Identity()*(-1.0/C.insC.correction_time_ba);
        //bg-bg
        F.block<3,3>(ibg,ibg)=Matrix3d::Identity()*(-1.0/C.insC.correction_time_bg);

        return MatrixXd::Identity(nx,nx)+F*dt;
    }

    Eigen::Vector3d CalculateGravity(const Vector3d coord_blh,bool is_ecef){
        if(is_ecef){
            const double constant_J2=0.00108263;
            const double constant_J4=-2.37091222e-6;
            const double constant_J6=6.08347e-9;
            double p = sqrt(coord_blh(0) * coord_blh(0) + coord_blh(1) * coord_blh(1) + coord_blh(2) * coord_blh(2));
            double t = coord_blh(2) / p;
            double a_p = WGS84_EARTH_LONG_RADIUS/ p;
            double a1 = -WGS84_GM / p / p;
            double a2 = 1 + 1.5 * constant_J2 * a_p * a_p - (15.0 / 8) * constant_J4 * pow(a_p,3) * a_p + (35.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a3 = -4.5 * constant_J2 * a_p * a_p + (75.0 / 4) * constant_J4 * pow(a_p,3) * a_p - (735.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a4 = -(175.0 / 8) * constant_J4 * pow(a_p,3) * a_p + (2205.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double a5 = -(1617.0 / 16) * constant_J6 * pow(a_p,3) * pow(a_p,3);

            double b1 = 3 * constant_J2 * a_p * a_p - (15.0 / 2) * constant_J4 * pow(a_p,3) * a_p + (105.0 / 8) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double b2 = (35.0 / 2) * constant_J4 * pow(a_p,3) * a_p - (945.0 / 12) * constant_J6 * pow(a_p,3) * pow(a_p,3);
            double b3 = (693.0 / 8) * constant_J6 * pow(a_p,3) * pow(a_p,3);

            double c1 = a2;
            double c2 = a3 - b1;
            double c3 = a4 - b2;
            double c4 = a5 - b3;
            double d1 = a2 + b1;
            double d2 = c2 + b2;
            double d3 = c3 + b3;
            double d4 = c4;
            Vector3d ge_vec;
            ge_vec(0) = (c1 + c2 * t * t + c3 * pow(t,3) * t + c4 * pow(t,3) * pow(t,3)) * coord_blh(0) * a1 / p + OMGE_GPS * OMGE_GPS * coord_blh(0);
            ge_vec(1) = (c1 + c2 * t * t + c3 * pow(t,3) * t + c4 * pow(t,3) * pow(t,3)) * coord_blh(1) * a1 / p + OMGE_GPS * OMGE_GPS * coord_blh(1);
            ge_vec(2) = (d1 + d2 * t * t + d3 * pow(t,3) * t + d4 * pow(t,3) * pow(t,3)) * coord_blh(2) * a1 / p;
            return ge_vec;
        }
        else{
            double gn = 9.7803267715 * (1 + 0.0052790414 * sin(coord_blh(0)) * sin(coord_blh(0)) + 0.0000232719 * pow(sin(coord_blh(0)),3) * sin(coord_blh(0)));
            gn += (-0.0000030876910891 + 0.0000000043977311 * sin(coord_blh(0)) * sin(coord_blh(0))) * coord_blh(2);
            gn += 0.0000000000007211 * coord_blh(2) * coord_blh(2);
            Vector3d gn_vec{0, 0, gn};
            return gn_vec;
        }
    }

    void cInsMech::TraceInsMechInfo(PPPLib::tImuInfoUnit &imu_info,bool prior,int idx) {
        LOG(DEBUG)<<"INS MECHANIZATION"<<(prior?"- ":"+ ")<< "("<<idx<<"): "<<imu_info.t_tag.GetTimeStr(4);
        LOG(DEBUG)<<"   "<<"GYRO VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_gyro.transpose()<<" rad/s";
        LOG(DEBUG)<<"   "<<"ACCE VALUE: "<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.cor_acce.transpose()<<" m/s^2";
        LOG(DEBUG)<<"   "<<"ATTITUDE:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.rpy.transpose()*R2D<<" deg";
        LOG(DEBUG)<<"   "<<"VELOCITY:(n)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.vn.transpose()<<" m/s";
        LOG(DEBUG)<<"   "<<"POSITION:(e)"<<setw(13)<<std::fixed<<setprecision(4)<<imu_info.re.transpose()<<" m";
        LOG(DEBUG)<<"   "<<"GYRO BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.bg.transpose()<<" ";
        LOG(DEBUG)<<"   "<<"ACCE BIAS:  "<<setw(13)<<std::fixed<<setprecision(6)<<imu_info.ba.transpose()<<" ";
    }

    void AdjustImuData(tImuDataUnit& imu_data,IMU_COORD_TYPE coord_type,IMU_DATA_FORMAT data_format,GYRO_DATA_FORMAT gyro_val_format,double dt) {
        Vector3d gyro,acce;

        if(coord_type==IMU_COORD_RFU){
            gyro=imu_data.gyro;
            acce=imu_data.acce;
            MatMul("NN",3,1,3,1.0,Crf,gyro.data(),0.0,imu_data.gyro.data());
            MatMul("NN",3,1,3,1.0,Crf,acce.data(),0.0,imu_data.acce.data());
        }
        if(data_format==IMU_FORMAT_INCR){
            for(int j=0;j<3;j++) imu_data.acce[j]/=dt;
            for(int j=0;j<3;j++) imu_data.gyro[j]/=dt;
        }
        if(gyro_val_format==GYRO_FORMAT_DEG){
            for(int j=0;j<3;j++) imu_data.gyro[j]*=D2R;
        }
    }

}