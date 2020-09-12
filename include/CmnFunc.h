//
// Created by cc on 7/9/20.
//

#ifndef PPPLIB_CMNFUNC_H
#define PPPLIB_CMNFUNC_H

#include <Eigen/Dense>
#include <iomanip>
#include <regex>
#include "LogInfo.h"
#include "Constant.h"

#ifdef WIN32
#define FILEPATHSEP '\\'
#else
#define FILEPATHSEP '/'
#include <dirent.h>
#endif

using namespace std;
using namespace Eigen;

namespace PPPLib{
    void CoutVector(VectorXd& vec,int p,int q,string s,bool trans);
    void CoutMatrix(MatrixXd& mat,int p,int q,string s,bool trans);

    int Round(double d);
    double VectorMean(vector<double>& seri);
    double NormDistribution(const double u);
    double ReNorm(double p);

    template <typename Iter1,typename Iter2>
    double Dot(const Iter1 VecA,const Iter2 VecB,int SizeVec){
        double dInn=0.0;

        while (--SizeVec>=0){
            dInn+=VecA[SizeVec]*VecB[SizeVec];
        }
        return dInn;
    }
    template <typename Iter>
    double Norm(const Iter VecA,int SizeVec){
        return sqrt(Dot(VecA,VecA,SizeVec));
    }

    template <typename Iter1,typename Iter2>
    int NormV3(const Iter1 vec1,Iter2 vec2){
        double r;
        if((r=Norm(vec1,3))<=0.0) return 0;
        vec2[0]=vec1[0]/r;
        vec2[1]=vec1[1]/r;
        vec2[2]=vec1[2]/r;
        return 1;
    }

    template <typename Iter1,typename Iter2, typename Iter3>
    void CrossVec3(const Iter1 vec1,const Iter2 vec2,Iter3 vec3){
        vec3[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
        vec3[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
        vec3[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
    }

    Eigen::Matrix3d VectorSkew(const Eigen::Vector3d& vec);
    extern void MatMul(const char *tr, int n, int k, int m, double alpha,
                       const double *A, const double *B, double beta, double *C);
    extern int MatInv(double *A, int n);

    string Doul2Str(int str_len, int dec_len, const string str_filler, const double src_num, string &dst_str);
    string Int2Str(int str_len, const string str_filler, const int src_num, string &dst_str);
    int Str2Double(string src_str,double &dst_num);
    int Str2Int(const string src_str,int &dst_num);
    string StrTrim(string s);
    void SplitString(const string& s, vector<string>& v,string c);
    vector<string> MultiSplitStr(const string &s, const string &seperator);
    vector<string> TextSplit(const string &in, const string &delim);
    extern void CreateDir(const char *path);

    typedef struct {
        time_t long_time;
        double sec;
    }tTime;

    typedef struct {
        long sn;
        double tos;
    }tSod;

    typedef struct {
        long day;
        tSod sod;
    }tMjd;

    class cTime {
    public:
        cTime();
        cTime(const double *ep);
        cTime(string str_time);
        cTime operator=(const tTime t);
        cTime operator+(double sec);
        void operator+=(double sec);
        ~cTime();

    public:
        double* GetEpoch();
        string GetTimeStr(int n);
        tMjd* GetMjd();
        int GetDoy();
        double TimeDiff(tTime t1);
        int Str2Time(string s);

        cTime *Epoch2Time(const double *ep);
        void Time2Epoch();
        double Time2Gpst(int* week,int* day, int sys);
        cTime* Gpst2Time(int week, double wos, int sys);
        cTime Utc2Gpst();
        cTime Gpst2Utc();
        cTime Gpst2Bdst();
        cTime Bdst2Gpst();
        double Utc2Gmst(double ut1_utc);
        cTime* AdjWeek(cTime t);
        cTime* AdjDay(cTime t);
        double Time2Doy();

    private:
        string Time2Str(int n);
        void Time2Mjd();
        double Time2Sec(cTime& day);

    private:
        double epoch_[6];
        string time_str_;
        int doy_;
        tMjd  mjd_;

    public:
        tTime t_;
    };

#if 0
    class cCoord{
    public:
        cCoord();
        cCoord(const Vector3d& coord, const COORDINATE_TYPE coord_type);
        cCoord(const double *coord, const COORDINATE_TYPE coord_type);
        ~cCoord();

    public:
        Vector3d GetCoordEnu(cCoord ref_coord);
        Vector3d GetCoordXyz();
        Vector3d GetCoordBlh();
        Vector3d GetCoordNed(cCoord ref_coord);
        Matrix3d GetCne(const COORDINATE_TYPE type);

    private:
        void Xyz2Blh();
        void Blh2Xyz();
        void Xyz2Enu(cCoord ref_coord);
        void Enu2Xyz(cCoord ref_coord);
        void Enu2Ned();
        void Ned2Enu();
        void CalcCne(const COORDINATE_TYPE type);
        double CalcLat();

    private:
        Vector3d coord_XYZ_;
        Vector3d coord_BLH_;
        Vector3d coord_ENU_;
        Vector3d coord_NED_;
        Matrix3d Cne_;
    };
#endif
    Vector3d Xyz2Blh(Vector3d& coord_xyz);
    Vector3d Blh2Xyz(Vector3d& coord_blh);
    Vector3d Xyz2Enu(Vector3d& coord_blh,Vector3d& coord_xyz);
    Vector3d Enu2Xyz(Vector3d& coord_blh,Vector3d& coord_enu);
    Vector3d Enu2Ned(Vector3d& coord_enu);
    Vector3d Ned2Enu(Vector3d& coord_ned);
    Matrix3d CalcCen(Vector3d& coord_blh,COORDINATE_TYPE lf_type);

    typedef struct{
        int nav_sys;
        double sample_rate;
        double ele_min;
        GNSS_FRQ_OPT frq_opt;
        int gnss_frq[NSYS+1][MAX_GNSS_USED_FRQ_NUM];
        bool adj_obs;
        bool csc;
        bool use_doppler;
        double code_phase_ratio;
        double meas_err_factor[5]; // [1-3] error factor a/b/c/of phase
        GNSS_AC_OPT  ac_opt;
        GNSS_EPH_OPT eph_opt;
        GNSS_ION_OPT ion_opt;
        GNSS_TRP_OPT trp_opt;
        GNSS_TID_OPT tid_opt;
        GLO_IFCB_OPT glo_ifcb_opt;
        bool sat_pcv;
        bool rec_ant;
        double cs_thres[2];   //mw and gf
        double max_pdop;
        double max_prior;
        int max_out;
        double max_inno;
        double ait_psd[3];
        bool check_dual_phase;
        GNSS_AR_MODE ar_mode;
        GLO_AR_MODE glo_ar_mode;
        bool bds_ar_mode;
        double ar_thres[8];
        double ar_el_mask;
        int min_sat_num2fix;
        int min_sat_num2drop;
        int min_lock2fix;
        bool partial_ar;
        bool res_qc;
        Vector3d rb;
    }tGnssConf;

    typedef struct{
        IMU_TYPE imu_type;
        IMU_COORD_TYPE coord_type;
        IMU_DATA_FORMAT data_format;
        GYRO_DATA_FORMAT gyro_val_format;
        double sample_rate;
        Vector3d lever;
        double correction_time_ba;
        double correction_time_bg;
        double init_pos_unc;
        double init_vel_unc;
        double init_att_unc;
        double init_ba_unc;
        double init_bg_unc;
        double psd_ba;
        double psd_bg;
        double psd_acce;
        double psd_gyro;
        bool err_model;
    }tInsConf;

    typedef struct{
        string rover;
        string base;
        string brd;
        string cbias;
        string clk;
        string sp3[3];
        string erp;
        string atx;
        string gim;
        string blq;
        string imu;
        string gsof;
        string ref;
        string sol;
        string sol_stat;
    }tFileConf;

    typedef struct{
        bool out_sol;
        bool out_head;
        bool out_vel;
        bool out_att;
        bool out_ba;
        bool out_bg;
        bool out_err;
        bool out_stat;
        int out_ins_mech_frq;
        COORDINATE_TYPE sol_coord;
    }tSolConf;

    typedef struct{
        cTime prc_date;
        string data_dir;
        string site_name;
        bool use_custom_dir;
        PPPLIB_MODE mode;
        PPPLIB_MODE_OPT mode_opt;
        int dynamic;
        SOLVE_ESTIMATOR estimator;
        tGnssConf gnssC;
        tInsConf  insC;
        tFileConf fileC;
        tSolConf  solC;
    }tPPPLibConf;

    typedef struct {
        cTime t_tag;
        SOL_STAT stat;
        SOL_INS_STAT ins_stat;
        int valid_sat_num;
        Vector3d pos;
        Vector3d vel;
        double q_pos[6];
        double q_vel[6];
        double clk_error[NSYS];
        double rec_dcb[NSYS];
        double rec_ifb[NSYS];
        Vector2d zenith_trp_delay; // dry and wet
        double dops[4];
        float ratio;
        double age;
        int num_ar_sat;
        double sigma;

        Vector3d att{0,0,0};  //roll pitch yaw
        double q_att[6];
        Vector3d gyro_bias{0,0,0};
        Vector3d accl_bias{0,0,0};
    }tSolInfoUnit;

    class cParSetting{
    public:
        cParSetting();
        cParSetting(tPPPLibConf conf);
        ~cParSetting();

    public:
        int GetGnssUsedFrqs();
        int GetNumObsType();
        int GetInsTransParNum(tPPPLibConf C);
        int GetPPPLibPar(tPPPLibConf C);
        int GetRealFixParNum(tPPPLibConf C);

        int NumPos();
        int NumVel();
        int NumAtt();
        int NumBa();
        int NumBg();
        int NumClPar();
        int NumClk();
        int NumClkDrift();
        int NumDcb();
        int NumIfb();
        int NumGloIfcb();
        int NumTrp();
        int NumIon();
        int NumAmb();

        int IndexPos();
        int IndexVel();
        int IndexAtt();
        int IndexBa();
        int IndexBg();
        int IndexClk(int sys_index);
        int IndexClkDritf(int sys_index);
        int IndexDcb(int sys_index);
        int IndexIfb(int sys_index);
        int IndexGloIfcb(int i);
        int IndexTrp();
        int IndexIon(int sat_no);
        int IndexAmb(int f,int sat_no);

    public:
        tPPPLibConf PPPLibC_;
    };

    class Config {
    public:
        using Ptr_=std::shared_ptr<Config>;

    private:
        Config() = default;
        Config(Config &&) = delete;
        Config(const Config &)= delete;
        Config &operator=(Config &&)= delete;
        Config &operator=(const Config &)= delete;
        static Ptr_ config_info_;
        std::map<std::string, std::string> storage_;

    public:
        ~Config() = default;

    public:
        static Ptr_ GetInstance();
        bool Open(std::string config_file);
        template <typename T>
        T Get(std::string key){
            transform(key.begin(),key.end(),key.begin(),::tolower);
            if(storage_.count(key)>0){
                try{
                    double value=stod(storage_[key]);
                    return static_cast<T>(value);
                }
                catch (const std::exception &e){
                    std::cerr<<e.what()<<'\n';
                }
            }
            else{
//            LOG(ERROR)<<"The key of "<<key<<" does not exist";
//            getchar();
                return T(0x0);
            }
        }
        template <typename  T>
        std::vector<T> GetArray(std::string key){
            std::vector<T> data;
            transform(key.begin(),key.end(),key.begin(),::tolower);
            if(storage_.count(key)>0){
                try{
                    auto text=TextSplit(storage_[key],",");
                    for(auto index:text){
                        double value=stod(index);
                        data.emplace_back(static_cast<T>(value));
                    }
                }
                catch(const std::exception &e){
                    std::cerr<<e.what()<<'\n';
                }
            }
            else{
//            LOG(ERROR)<<"The key of "<<key<<" does not exist";
//            getchar();
            }
            return data;
        }
    };

    template <>
    inline std::string Config::Get<std::string>(std::string key){
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(storage_.count(key)>0){
            try{
                return std::string(storage_[key]);
            }
            catch (const std::exception &e){
                std::cerr<<e.what()<<'\n';
            }
        }
        else{
//        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
//            getchar();
            return "";
        }
    }

    template <>
    inline std::vector<std::string> Config::GetArray<std::string>(std::string key){
        std::vector<std::string> data;
        transform(key.begin(),key.end(),key.begin(),::tolower);
        if(storage_.count(key)>0){
            try{
                data=TextSplit(storage_[key],",");
            }
            catch(const std::exception &e){
                std::cerr<<e.what()<<'\n';
            }
        }
        else{
//        LOG(ERROR)<<"The key of "<<key<<" does not exist"<<endl;
            getchar();
        }
        return data;
    }




}


#endif //PPPLIB_CMNFUNC_H
