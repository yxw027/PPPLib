//
// Created by cc on 8/22/20.
//

#include "DecodeRaw.h"

#define GSOFHEADLEN   4                       /* number bytes of gsof message header */
#define GSOFMIDLEN    3                       /* number bytes of gsof message information */
#define GSOFMMLEN     2                       /* number bytes of gsof message record header */
#define GSOFTYPE01    1                       /* type of gsof message */
#define GSOFTYPE02    2
#define GSOFTYPE03    3
#define GSOFTYPE06    6
#define GSOFTYPE08    8
#define GSOFTYPE09    9
#define GSOFTYPE10    10
#define GSOFTYPE11    11
#define GSOFTYPE12    12
#define GSOFTYPE38    38
#define GSOFTYPE41    41

#define U1(p)       (*((unsigned char *)(p)))
#define I1(p)       (*((signed char *)(p)))

#define FREQOCXO    1E8                 /* crystal frequency (100Mhz) */
#define NUMBYTES_GI310  43                      /* numbers of bytes of gi310 imu raw data */
#define MAXDIFFTIME     10.0                    /* max time difference to reset  */
const unsigned char gi310_head[2]={0x55,0xAA};  /* imu message header */
static unsigned short U2(unsigned char *p) {unsigned short u; memcpy(&u,p,2); return u;}
static unsigned int   U4(unsigned char *p) {unsigned int   u; memcpy(&u,p,4); return u;}
static short          I2(unsigned char *p) {short          i; memcpy(&i,p,2); return i;}
static int            I4(unsigned char *p) {int            i; memcpy(&i,p,4); return i;}
static float          R4(unsigned char *p) {float          r; memcpy(&r,p,4); return r;}
static double         R8(unsigned char *p) {double         r; memcpy(&r,p,8); return r;}

namespace PPPLib{

    cDecodeGsof::cDecodeGsof() {}

    cDecodeGsof::~cDecodeGsof() {}

    bool cDecodeGsof::DecodeGsof(const string file, vector<tSolInfoUnit>& gsofs) {
        int data;
        tRaw raw={{0}};
        tSolInfoUnit sol;
        FILE *fp= nullptr;

        if(!(fp=fopen(file.c_str(),"r"))){
            LOG(WARNING)<<"FILE "<<file<<" OPEN ERROR";
            return false;
        }

        while(true) {
            if((data=fgetc(fp))==EOF) break;
            if((InputGsofMessage(&raw,(unsigned char)data))){
                GsofUnit2ppplibsol(raw.gsof_unit,sol);
                if(gsofs.size()>=3){
                    if(fabs(gsofs.back().t_tag.TimeDiff((gsofs.end()-2)->t_tag.t_))<DTTOL) continue;
                    else gsofs.push_back(sol);
                }
                else gsofs.push_back(sol);
            }
        }

        sort(gsofs.begin(),gsofs.end(),CmpSolInfo);
        fclose(fp);
    }

    bool cDecodeGsof::InputGsofMessage(PPPLib::tRaw *raw, unsigned char data) {
        if(raw->dire) return DecodeMessageBackward(raw,data);
        else          return DecodeMessageForward(raw,data);
    }

    int cDecodeGsof::DecodeMessageForward(PPPLib::tRaw *raw, unsigned char data) {
        raw->buff[raw->nbyte++]=data;

        if(raw->nbyte<GSOFHEADLEN) return false;

        if(CheckHead(raw)){
            raw->len=U1(raw->buff+3);
        }
        else{
            raw->nbyte=0;return false;
        }
        if(raw->nbyte<raw->len+6) return false;

        if(!CheckSum(raw)){
            raw->buff[0]='\0';
            raw->nbyte=0;
            return -1;
        }

        DecodeMessage(raw);

        raw->buff[0]='\0';raw->nbyte=0;
        return 6;
    }

    int cDecodeGsof::DecodeMessageBackward(PPPLib::tRaw *raw, unsigned char data) {

    }

    void cDecodeGsof::DecodeMessage(PPPLib::tRaw *raw) {
        int type,len,cb=0,nt=GSOFHEADLEN+GSOFMIDLEN;

        raw->gsof_unit.rcv_status=U1(raw->buff+1);
        raw->gsof_unit.no=U1(raw->buff+4);

        while(cb<raw->len-GSOFMIDLEN&&raw->len>0){
            type=U1(raw->buff+nt+cb);
            len=U1(raw->buff+nt+cb+1);
            switch(type){
                case GSOFTYPE01 : DecodeType01(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE02 : DecodeType02(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE03 : DecodeType03(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE06 : DecodeType06(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE08 : DecodeType08(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE09 : DecodeType09(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE10 : DecodeType10(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE11 : DecodeType11(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE12 : DecodeType12(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE38 : DecodeType38(raw,raw->buff+GSOFMMLEN+nt+cb); break;
                case GSOFTYPE41 : DecodeType41(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            }
            cb+=(len+2);
        }
    }

    bool cDecodeGsof::CheckHead(PPPLib::tRaw *raw) {
        return raw->buff[0]==0x02&&raw->buff[2]==0x40;
    }

    bool cDecodeGsof::CheckSum(PPPLib::tRaw *raw) {
        unsigned char cs=0;
        int i;
        for(i=1;i<raw->len+4;i++) cs+=raw->buff[i];
        return cs%256==raw->buff[raw->len+4]&&raw->buff[raw->len+5]==0x03;
    }

    unsigned short cDecodeGsof::GetU16(unsigned char **data) {
        unsigned short ret;
        unsigned char *pb;

        pb=(unsigned char*)(&ret)+1;
        *pb--=*(*data)++;*pb=*(*data)++;
        return ret;
    }

    unsigned long cDecodeGsof::GetU32(unsigned char **data) {
        unsigned int ret;
        unsigned char *pb;

        pb=(unsigned char*)(&ret)+3;
        *pb--=*(*data)++; *pb--=*(*data)++;
        *pb--=*(*data)++; *pb  =*(*data)++;
        return ret;
    }

    float cDecodeGsof::GetFloat(unsigned char **data) {
        float ret;
        unsigned char *pb;

        pb=(unsigned char*)(&ret)+3;
        *pb--=*(*data)++; *pb--=*(*data)++;
        *pb--=*(*data)++; *pb  =*(*data)++;
        return ret;
    }

    double cDecodeGsof::GetDouble(unsigned char **data) {
        double ret;
        unsigned char *pb;

        pb=(unsigned char *)(&ret)+7;
        *pb--=*(*data)++; *pb--=*(*data)++;
        *pb--=*(*data)++; *pb--=*(*data)++;
        *pb--=*(*data)++; *pb--=*(*data)++;
        *pb--=*(*data)++; *pb  =*(*data)++;
        return ret;
    }

    void cDecodeGsof::DecodeType01(PPPLib::tRaw *raw, unsigned char *data) {
        int week;double sow;

        raw->gsof_unit.num_use_sat=data[6];
        if((data[8]&0x00)==0x00) raw->gsof_unit.sol_level=SOL_SPP;
        if((data[8]&0x01)==0x01) raw->gsof_unit.sol_level=SOL_DGPS;
        if((data[8]&0x03)==0x03) raw->gsof_unit.sol_level=SOL_FLOAT;
        if((data[8]&0x07)==0x07) raw->gsof_unit.sol_level=SOL_FIX;

        sow=GetU32(&data)/1000.0; week=GetU16(&data);
        raw->gsof_unit.t_tag.Gpst2Time(week,sow,SYS_GPS);
    }

    void cDecodeGsof::DecodeType02(PPPLib::tRaw *raw, unsigned char *data) {
        raw->gsof_unit.llh[0]=GetDouble(&data);
        raw->gsof_unit.llh[1]=GetDouble(&data);
        raw->gsof_unit.llh[2]=GetDouble(&data);

        raw->gsof_unit.pos=Blh2Xyz(raw->gsof_unit.llh);
    }

    void cDecodeGsof::DecodeType03(PPPLib::tRaw *raw, unsigned char *data) {
        raw->gsof_unit.pos[0]=GetDouble(&data);
        raw->gsof_unit.pos[1]=GetDouble(&data);
        raw->gsof_unit.pos[2]=GetDouble(&data);

        raw->gsof_unit.llh=Xyz2Blh(raw->gsof_unit.pos);
    }

    void cDecodeGsof::DecodeType06(PPPLib::tRaw *raw, unsigned char *data) {
        raw->gsof_unit.rb_delta[0]=GetDouble(&data);
        raw->gsof_unit.rb_delta[1]=GetDouble(&data);
        raw->gsof_unit.rb_delta[2]=GetDouble(&data);
    }

    void cDecodeGsof::DecodeType08(PPPLib::tRaw *raw, unsigned char *data) {
        double s,head;
        raw->gsof_unit.vel_flag=data[0]; data++;
        s   =GetFloat(&data);
        head=GetFloat(&data);
        raw->gsof_unit.vel[0]=s*sin(head); /* east */
        raw->gsof_unit.vel[1]=s*cos(head); /* north */
        raw->gsof_unit.vel[2]=GetFloat(&data);  /* up */
    }

    void cDecodeGsof::DecodeType09(PPPLib::tRaw *raw, unsigned char *data) {
        raw->gsof_unit.dop[0]=GetDouble(&data);
        raw->gsof_unit.dop[1]=GetDouble(&data);
        raw->gsof_unit.dop[2]=GetDouble(&data);
        raw->gsof_unit.dop[3]=GetDouble(&data);
    }

    void cDecodeGsof::DecodeType10(PPPLib::tRaw *raw, unsigned char *data) {
//        raw->gsof_unit.clk.flags=data[0]; data++;
//        raw->gsof_unit.clk.off  =GetDouble(&data)/1000.0;
//        raw->gsof_unit.clk.foff =GetDouble(&data)/1E6;
    }

    void cDecodeGsof::DecodeType11(PPPLib::tRaw *raw, unsigned char *data) {
        raw->gsof_unit.cov[0]=GetFloat(&data);
        raw->gsof_unit.cov[1]=GetFloat(&data);
        raw->gsof_unit.cov[2]=GetFloat(&data);
        raw->gsof_unit.cov[3]=GetFloat(&data);
        raw->gsof_unit.cov[4]=GetFloat(&data);
        raw->gsof_unit.cov[5]=GetFloat(&data);
        raw->gsof_unit.cov[6]=GetFloat(&data);
        raw->gsof_unit.cov[7]=GetFloat(&data);
    }

    void cDecodeGsof::DecodeType12(PPPLib::tRaw *raw, unsigned char *data) {
        raw->gsof_unit.sig[0]=GetFloat(&data);
        raw->gsof_unit.sig[1]=GetFloat(&data); /* east */
        raw->gsof_unit.sig[2]=GetFloat(&data); /* north */
        raw->gsof_unit.sig[3]=GetFloat(&data); /* east-north */
        raw->gsof_unit.sig[4]=GetFloat(&data); /* up */
        GetFloat(&data); GetFloat(&data); GetFloat(&data);
        raw->gsof_unit.sig[5]=GetFloat(&data); /* unit variance */
    }

    void cDecodeGsof::DecodeType38(PPPLib::tRaw *raw, unsigned char *data) {

    }

    void cDecodeGsof::DecodeType41(PPPLib::tRaw *raw, unsigned char *data) {
        int week; double sow; sow=GetU32(&data)/1000.0; week=GetU16(&data);
//        raw->gsof_unit.base.t      =gpst2time(week,sow);
//        raw->gsof_unit.base.pos[0] =GetDouble(&data);
//        raw->gsof_unit.base.pos[1] =GetDouble(&data);
//        raw->gsof_unit.base.pos[2] =GetDouble(&data);
//        raw->gsof_unit.base.quality=data[30];
    }

    void cDecodeGsof::GsofUnit2ppplibsol(PPPLib::tGsofUnit &gsof, PPPLib::tSolInfoUnit &ppplib_sol) {
        Matrix3d Cen,covn,cove;
        Vector3d std,e;
        covn=Matrix3d::Zero();
        cove=Matrix3d::Zero();

        ppplib_sol.t_tag=gsof.t_tag;

        ppplib_sol.pos=gsof.pos;
        Cen=CalcCen(gsof.llh,COORD_ENU);
        ppplib_sol.vel=Cen.transpose()*gsof.vel;

        covn(0,0)=SQR(gsof.sig[2]==0.0?30.0:gsof.sig[2]);
        covn(1,1)=SQR(gsof.sig[1]==0.0?30.0:gsof.sig[1]);
        covn(2,2)=SQR(gsof.sig[4]==0.0?30.0:gsof.sig[4]);
        covn(0,1)=covn(1,0)=(gsof.sig[3]<0.0?-SQR(gsof.sig[3]):SQR(gsof.sig[3]));

        Cen=CalcCen(gsof.llh,COORD_NED);
        cove=Cen.transpose()*covn*Cen;

        for(int i=0;i<3;i++) ppplib_sol.q_pos[i]=cove(i,i);
        for(int i=0;i<3;i++) ppplib_sol.q_vel[i]=0.0;

        ppplib_sol.valid_sat_num=gsof.num_use_sat;
        ppplib_sol.stat= static_cast<SOL_STAT>(gsof.sol_level);
    }

    bool cDecodeGsof::CmpSolInfo(const PPPLib::tSolInfoUnit &p1, const PPPLib::tSolInfoUnit &p2) {
        cTime t1=p1.t_tag,t2=p2.t_tag;
        return t1.TimeDiff(t2.t_)<0.0;
    }

    cDecodeImuM39::cDecodeImuM39() {}

    cDecodeImuM39::~cDecodeImuM39() {}

    int cDecodeImuM39::DecodeM39(string file,tInsConf insC, vector<tImuDataUnit> &imus) {
        int data;
        FILE *fp= nullptr;
        tRaw raw={{0}};

        if(!(fp=fopen(file.c_str(),"r"))){
            LOG(WARNING)<<"FILE "<<file<<" OPEN ERROR";
            return false;
        }

        while(true){
            if((data=fgetc(fp))==EOF) break;
            if((InputM39Message(&raw,(unsigned char)data))){
                if(imus.size()>=3){
                    if(imus.back().t_tag.TimeDiff((imus.end()-2)->t_tag.t_)!=0.0){
                        AdjustImuData(raw.imu_unit,insC.coord_type,insC.data_format,insC.gyro_val_format,insC.sample_rate);
                        imus.push_back(raw.imu_unit);
                    }
                    else continue;
                }
                else{
                    AdjustImuData(raw.imu_unit,insC.coord_type,insC.data_format,insC.gyro_val_format,insC.sample_rate);
                    imus.push_back(raw.imu_unit);
                }
            }
        }

//        sort(imus.begin(),imus.end(),CmpImuData);
        qsort(imus.data(), imus.size(), sizeof(tImuDataUnit), CmpImuData);
        fclose(fp);
        return imus.size();
    }

    int cDecodeImuM39::InputM39Message(PPPLib::tRaw *raw, unsigned char data) {
        raw->len=NUMBYTES_GI310;
        if(raw->dire) return DecodeM39Backward(raw,data);
        else          return DecodeM39Forward(raw,data);
    }

    bool cDecodeImuM39::DecodeM39Forward(PPPLib::tRaw *raw, unsigned char data) {
        raw->buff[raw->nbyte++]=data;

        if(raw->nbyte<2) return false;

        if(!CheckHead(raw->buff)){
            raw->nbyte=0;
            return false;
        }
        if(raw->nbyte<NUMBYTES_GI310) return false;

        if(CheckSum(raw->buff+2,40)==data){
            raw->nbyte=0;
            return DecodeMessage(raw);
        }
        else{
            raw->nbyte=0;
            return false;
        }
    }

    bool cDecodeImuM39::DecodeM39Backward(PPPLib::tRaw *raw, unsigned char data) {

    }

    int cDecodeImuM39::DecodeMessage(PPPLib::tRaw *raw) {
        int i,week=0;
        unsigned int dc=0;
        static int start=0;
        static double sow=0.0,timeu=0.0;
        cTime t0;

        DecodeSowTime(raw,&sow,&start);

        if(start){
            AdjustTime(raw,sow,&sow,&timeu,&week,&dc);
        }
        else return 0;

        if(dc*1.0/FREQOCXO>MAXDIFFTIME) return 0;

        t0.Gpst2Time(week,timeu,SYS_GPS);

        raw->imu_unit.t_tag=t0;
        raw->imu_unit.pps=U4(raw->buff+06);
        raw->imu_unit.imu_cnt=U4(raw->buff+10);
        for(int i=0;i<3;i++) {
            raw->imu_unit.gyro[i]=R4(raw->buff+14+i*4);
            raw->imu_unit.acce[i]=R4(raw->buff+26+i*4);
        }

        return timeu>0.0?4:0;
    }

    void cDecodeImuM39::DecodeSowTime(PPPLib::tRaw *raw, double *sow, int *start) {
        static unsigned int pps=0;

        if(pps==0){
            pps=U4(raw->buff+6);
            *sow=U4(raw->buff+2);
            return;
        }
        if(pps!=U4(raw->buff+6)){
            (*sow)++;*start=1;
        }
        pps=U4(raw->buff+6);
    }

    void cDecodeImuM39::AdjustTime(PPPLib::tRaw *raw, const double sowi, double *sowo, double *timu, int *week,
                                   unsigned int *dcc) {
        static unsigned int imuc=0,dc=0;
        int d=0;

        if (imuc==0) {
            imuc=U4(raw->buff+10); *sowo=sowi;
            return;
        }
        *dcc=dc=U4(raw->buff+10)-imuc<0?UINT_MAX+U4(raw->buff+10)-imuc:U4(raw->buff+10)-imuc;
        imuc=U4(raw->buff+10);

        /* increase week */
        if ((*sowo=(int)(1.0/FREQOCXO*dc+sowi))>=604800.0) {
            *sowo-=604800.0; (*week)++;
        }
        d=U4(raw->buff+10)-U4(raw->buff+06);
        *timu=*sowo+1.0/FREQOCXO*d;
    }

    static double res=2048; /* odometry resolution */
    static double d=0.73;   /* odometry wheel diameter (m) */
    bool cDecodeImuM39::DecodeM39Odo(PPPLib::tRaw *raw, double dt) {
        static int dc;

        if (raw->imu_unit.t_tag.t_.long_time==0||dt==0.0) {
            raw->imu_unit.odo_cnt=I2(raw->buff+38); return 0;
        }
        dc=I2(raw->buff+38)-raw->imu_unit.odo_cnt<=SHRT_MIN?
           I2(raw->buff+38)-raw->imu_unit.odo_cnt+USHRT_MAX:
           I2(raw->buff+38)-raw->imu_unit.odo_cnt>=SHRT_MAX?
           I2(raw->buff+38)-raw->imu_unit.odo_cnt-USHRT_MAX:
           I2(raw->buff+38)-raw->imu_unit.odo_cnt;

        raw->imu_unit.odo_cnt=I2(raw->buff+38);
        raw->imu_unit.odo.t_tag=raw->imu_unit.t_tag;
        raw->imu_unit.odo.dt=dt;
        raw->imu_unit.odo.dr=dc/res*PI*d;
        return true;
    }

    bool cDecodeImuM39::CheckHead(const unsigned char *buff) {
        return (buff[0]==gi310_head[0])&&(buff[1]==gi310_head[1]);
    }

    unsigned char cDecodeImuM39::CheckSum(const unsigned char *buff, int len) {
        int i;
        unsigned char sum=0;
        for(i=0;i<len;i++) sum+=buff[i];return sum;
    }

    int cDecodeImuM39::CmpImuData(const void *p1, const void *p2) {
        tImuDataUnit *q1=(tImuDataUnit *)p1, *q2=(tImuDataUnit *)p2;
        cTime t1=q1->t_tag,t2=q2->t_tag;
        return t1.TimeDiff(t2.t_)<0.0?-1:1;
    }
}