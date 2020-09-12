//
// Created by cc on 7/9/20.
//

#include "CmnFunc.h"

namespace PPPLib{

    void CoutVector(VectorXd& vec,int p,int q,string s,bool trans){
        if(s.length()) cout<<s<<endl;
        cout<<"------------------------------------------------------------"<<endl;
        cout<<std::fixed<<setprecision(q)<<setw(p)<<(trans?vec.transpose():vec)<<endl;
        cout<<"------------------------------------------------------------"<<endl;
    }

    void CoutMatrix(MatrixXd& mat,int p,int q,string s,bool trans){
        if(s.length()) cout<<s<<endl;
        cout<<"------------------------------------------------------------"<<endl;
        cout<<setw(p)<<std::fixed<<setprecision(q)<<(trans?mat.transpose():mat)<<endl;
        cout<<"------------------------------------------------------------"<<endl;
    }

    int Round(double d){
        int i;
        if(d>=0) i=(int)(d+0.5);
        else i=(int)(d-0.5);
        return i;
    }

    double VectorMean(vector<double>& seri){
        if(seri.size()<=0) return 0.0;

        double total=0.0;
        for(int i=0;i<seri.size();i++){
            total+=seri[i];
        }

        return total/seri.size();
    }

    double NormDistribution(const double u){
        if(u<-5.0) return 0.0;
        if(u>5.0) return 1.0;

        double y=fabs(u)/sqrt(2.0);

        double p=1.0+y*(0.0705230784+y*(0.0422820123+y*(0.0092705272+
                                                        y*(0.0001520143+y*(0.0002765672+y*0.0000430638)))));
        double er=1.0-pow(p,-16.0);
        p=(u<0.0)?0.5-0.5*er:0.5+0.5*er;
        return p;
    }

    double ReNorm(double p){
        if(p==0.5) return 0.0;
        if(p>0.9999997) return 5.0;
        if(p<0.0000003) return -5.0;
        if(p<0.5) return -ReNorm(1.0-p);

        double y=-log(4.0*p*(1.0-p));
        y=y*(1.570796288+y*(0.3706987906e-1
                            +y*(-0.8364353589e-3+y*(-0.2250947176e-3
                                                    +y*(0.6841218299e-5+y*(0.5824238515e-5
                                                                           +y*(-0.1045274970e-5+y*(0.8360937017e-7
                                                                                                   +y*(-0.3231081277e-8+y*(0.3657763036e-10
                                                                                                                           +y*0.6936233982e-12))))))))));
        return sqrt(y);
    }

    extern void MatMul(const char *tr, int n, int k, int m, double alpha,
                       const double *A, const double *B, double beta, double *C)
    {
        double d;
        int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);

        for (i=0;i<n;i++) for (j=0;j<k;j++) {
                d=0.0;
                switch (f) {
                    case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
                    case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
                    case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
                    case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
                }
                if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
            }
    }

    extern double *mat(int n, int m)
    {
        double *p;

        if (n<=0||m<=0) return NULL;
        if (!(p=(double *)malloc(sizeof(double)*n*m))) {
            LOG(ERROR)<<"matrix memory allocation error";
        }
        return p;
    }

    extern int *imat(int n, int m)
    {
        int *p;

        if (n<=0||m<=0) return NULL;
        if (!(p=(int *)malloc(sizeof(int)*n*m))) {
            LOG(ERROR)<<"integer matrix memory allocation error";
        }
        return p;
    }

    extern void matcpy(double *A, const double *B, int n, int m)
    {
        memcpy(A,B,sizeof(double)*n*m);
    }

    /* LU decomposition ----------------------------------------------------------*/
    static int ludcmp(double *A, int n, int *indx, double *d)
    {
        double big,s,tmp,*vv=mat(n,1);
        int i,imax=0,j,k;

        *d=1.0;
        for (i=0;i<n;i++) {
            big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
            if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
        }
        for (j=0;j<n;j++) {
            for (i=0;i<j;i++) {
                s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
            }
            big=0.0;
            for (i=j;i<n;i++) {
                s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
                if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
            }
            if (j!=imax) {
                for (k=0;k<n;k++) {
                    tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
                }
                *d=-(*d); vv[imax]=vv[j];
            }
            indx[j]=imax;
            if (A[j+j*n]==0.0) {free(vv); return -1;}
            if (j!=n-1) {
                tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
            }
        }
        free(vv);
        return 0;
    }
/* LU back-substitution ------------------------------------------------------*/
    static void lubksb(const double *A, int n, const int *indx, double *b)
    {
        double s;
        int i,ii=-1,ip,j;

        for (i=0;i<n;i++) {
            ip=indx[i]; s=b[ip]; b[ip]=b[i];
            if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
            b[i]=s;
        }
        for (i=n-1;i>=0;i--) {
            s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
        }
    }

    extern int MatInv(double *A, int n)
    {
        double d,*B;
        int i,j,*indx;

        indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
        if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
        for (j=0;j<n;j++) {
            for (i=0;i<n;i++) A[i+j*n]=0.0;
            A[j+j*n]=1.0;
            lubksb(B,n,indx,A+j*n);
        }
        free(indx); free(B);
        return 0;
    }

    string Doul2Str(int str_len, int dec_len, const string str_filler, const double src_num, string &dst_str){
        dst_str=to_string(src_num); /* with 6 decimal digit */
        if(dec_len>0)dst_str=dst_str.substr(0,dst_str.length()-6+dec_len); /* decimal digits */
        else dst_str=dst_str.substr(0,dst_str.length()-7);
        while (str_len>dst_str.length())
            dst_str=str_filler+dst_str;
        return dst_str;
    }

    string Int2Str(int str_len, const string str_filler, const int src_num, string &dst_str){
        dst_str=to_string(src_num);
        while(str_len>dst_str.length())
            dst_str=str_filler+dst_str;
        return dst_str;
    }

    int Str2Double(string src_str,double &dst_num){
        int i,fnum;

        if (src_str.length()<=0) return 0;

        for (i=0,fnum=1; i<src_str.length(); i++){
            if (fnum && src_str[i]<='9' && src_str[i]>='0') fnum=0;
            if (fnum==0 && src_str[i]=='D') { src_str[i]='E'; break; }
        }
        if (fnum) return 0;

        dst_num=stod(src_str);

        return 1;
    }

    int Str2Int(const string src_str,int &dst_num){
        int i;
        if(src_str.empty()) dst_num=0;
        for (i=0; i<src_str.length(); i++){
            if (src_str[i]<='9'&&src_str[i]>='0') break;
        }
        if (i>=src_str.length()) return 0;

        dst_num=stoi(src_str);

        return 1;
    }

    string StrTrim(string s){
        if(s.empty()) return s;
        s.erase(0,s.find_first_not_of(" "));
        s.erase(s.find_last_not_of(" ")+1);
        return s;
    }

    void SplitString(const string& s, vector<string>& v,string c)
    {
        string::size_type pos1,pos2;
        pos2=s.find(c);
        pos1=0;
        while(string::npos!=pos2){
            v.push_back(s.substr(pos1,pos2-pos1));
            pos1=pos2+c.size();
            pos2=s.find(c,pos1);
        }
        if(pos1!=s.length())
            v.push_back(s.substr(pos1));
    }

    vector<string> MultiSplitStr(const string &s, const string &seperator) {
        vector<string> result;
        typedef string::size_type string_size;
        string_size i = 0;

        while(i!=s.size()){
            int flag = 0;
            while(i!=s.size()&&flag==0){
                flag = 1;
                for(string_size x =0;x<seperator.size();++x){
                    if(s[i]==seperator[x]){
                        ++i;flag=0;break;
                    }
                }
            }
            flag=0;string_size j=i;
            while(j!=s.size()&&flag==0){
                for(string_size x=0;x<seperator.size();++x){
                    if(s[j]==seperator[x]){
                        flag=1;break;
                    }
                }
                if(flag==0) ++j;
            }
            if(i!=j){
                result.push_back(s.substr(i,j-i));i=j;
            }
        }
        return result;
    }

    vector<string> TextSplit(const string &in, const string &delim) {
        std::vector<std::string> ret;
        try
        {
            std::regex re{delim};
            return std::vector<std::string>{std::sregex_token_iterator(in.begin(), in.end(), re, -1), std::sregex_token_iterator()};
        }
        catch (const std::exception &e)
        {
            std::cout << "error:" << e.what() << std::endl;
        }
        return ret;
    }

    extern void CreateDir(const char *path) {
        char buff[1024],*p;

        strcpy(buff,path);
        if (!(p=strrchr(buff,FILEPATHSEP))) return;
        *p='\0';

#ifdef WIN32
        CreateDirectory(buff,NULL);
#else
        mkdir(buff,0777);
#endif
    }


    cTime::cTime(){
        for(double & i : epoch_) i=0.0;
        epoch_[0]=0;epoch_[1]=0;epoch_[2]=0;
        Epoch2Time(epoch_);
    }

    cTime::cTime(const double *ep){
        if(ep){
            for(int i=0;i<6;i++) epoch_[i]=ep[i];
            Epoch2Time(ep);
        }
        else{
            for(int i=0;i<6;i++) epoch_[i]=0.0;
            Epoch2Time(epoch_);
        }

    }

    cTime::cTime(string s):time_str_(s){
        Str2Time(s);
    }

    cTime cTime::operator=(const tTime t) {
        t_=t;
    }

    cTime cTime::operator+(double sec){
#if 1
        double tt;
        this->t_.sec+=sec;
        tt=floor(this->t_.sec);
        this->t_.long_time+=(int)tt;
        this->t_.sec-=tt;
        return *this;
#endif

    }

    void cTime::operator+=(double sec) {
        double tt;
        t_.sec+=sec;
        tt=floor(t_.sec);
        t_.long_time+=(int)tt;
        t_.sec-=tt;
    }

    cTime::~cTime() {}

    double* cTime::GetEpoch() {
        return epoch_;
    }

    string cTime::GetTimeStr(int n) {
        Time2Str(n);
        return time_str_;
    }

    int cTime::GetDoy() {
        Time2Doy();
        return doy_;
    }

    tMjd* cTime::GetMjd() {
        Time2Mjd();
        return &mjd_;
    }

    string cTime::Time2Str(int n){
        if (n<0) n=0; else if (n>12) n=12;
        string str;
        if (1.0-t_.sec<0.5/pow(10.0,n)) { t_.long_time++; t_.sec=0.0; };
        Time2Epoch();
        time_str_=Int2Str(4,"0",(int)epoch_[0],str)+"/"+Int2Str(2,"0",(int)epoch_[1],str)+"/"+
                  Int2Str(2,"0",(int)epoch_[2],str)+" "+Int2Str(2,"0",(int)epoch_[3],str)+":"+
                  Int2Str(2,"0",(int)epoch_[4],str)+":"+Doul2Str(2+n+1,n,"0",epoch_[5],str);
        return time_str_;
    }

    cTime* cTime::Epoch2Time(const double *ep){
        const int doy[]={ 1,32,60,91,121,152,182,213,244,274,305,335 };
        int days, dsec, year=int(ep[0]), mon=(int)ep[1], day=(int)ep[2];

        if (year<1970||2099<year||mon<1||12<mon) {
            t_.long_time=0;t_.sec=0; return this;
        }

        /* leap year if year%4==0 in 1901-2099 */
        days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3 ? 1 : 0);
        dsec=(int)floor(ep[5]);
        t_.long_time=(time_t)days*86400+(time_t)ep[3]*3600+(time_t)ep[4]*60+dsec;
        t_.sec=ep[5]-dsec;

        return this;
    }

    void cTime::Time2Epoch() {
        const int mday[]={ /* # of days in a month */
                31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
                31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
        };
        int days, dsec, mon, day;

        /* leap year if year%4==0 in 1901-2099 */
        days=(int)(t_.long_time/86400);
        dsec=(int)(t_.long_time-(time_t)days*86400);
        for (day=days%1461, mon=0; mon<48; mon++) {
            if (day>=mday[mon]) day-=mday[mon]; else break;
        }
        epoch_[0]=1970+days/1461*4+mon/12; epoch_[1]=mon%12+1; epoch_[2]=day+1;
        epoch_[3]=dsec/3600; epoch_[4]=dsec%3600/60; epoch_[5]=dsec%60+t_.sec;
    }

    int cTime::Str2Time(string s) {
        if (sscanf(s.c_str(),"%lf %lf %lf %lf %lf %lf",epoch_,epoch_+1,epoch_+2,epoch_+3,epoch_+4,epoch_+5)<6){
            LOG(ERROR)<<"time string error";
            return -1;
        }

        if (epoch_[0]<100) epoch_[0]+=2000;
        if (epoch_[0]<=1900||epoch_[1]==0||epoch_[2]==0){
            LOG(ERROR)<<"time format error";
            return -1;
        }

        Epoch2Time(epoch_);

        time_str_=s;

        return 0;
    }

    double cTime::TimeDiff(PPPLib::tTime t1) {
        return difftime(t_.long_time,t1.long_time)+t_.sec-t1.sec;
    }

    void cTime::Time2Mjd() {
        this->Time2Epoch();
        int year=(int)Round(epoch_[0]);
        int mon=(int)Round(epoch_[1]);
        if(mon<=2){
            year=year-1;
            mon=mon+12;
        }

        int a=(int)(365.25*year);
        int b=(int)(30.6001)*(mon+1);
        mjd_.day=a+b+Round(epoch_[2])-679019;
        mjd_.sod.sn=Round(epoch_[3])*3600+Round(epoch_[4])*60+Round(epoch_[5]);
        mjd_.sod.tos=epoch_[5]-Round(epoch_[5]);
    }

    double cTime::Time2Doy() {
        double ep0[6]={0};

        Time2Epoch();
        ep0[0]=epoch_[0];ep0[1]=ep0[2]=1.0;
        cTime t0(ep0);
        doy_=(int) (TimeDiff(t0.t_))/86400.0+1.0;
        return doy_;
    }

    double cTime::Time2Sec(cTime& day){
        double sec;
        double ep0[6]={0};
        int i;
        Time2Epoch();

        sec=epoch_[3]*3600.0+epoch_[4]*60.0+epoch_[5];
        for(i=0;i<3;i++) ep0[i]=epoch_[i];
        day.Epoch2Time(ep0);
        return sec;
    }


    double cTime::Time2Gpst(int* week, int* day, int sys) {
        const double *ep=kGpsTimeStart;
        switch(sys){
            case SYS_BDS: ep=kBdsTimeStart;break;
            case SYS_GAL: ep=kGalTimeStart;break;
        }
        cTime t0(ep);
        time_t sec;

        sec=t_.long_time-t0.t_.long_time;

        int w=(int)(sec/(86400*7));
        if(week) *week=w;

        double s=(double)(sec-(double)w*86400*7)+t_.sec;
        int d=(int)s/86400;
        if(day) *day=d;

        return s;
    }

    cTime* cTime::Gpst2Time(int week, double wos, int sys) {
        const double *ep=kGpsTimeStart;
        switch(sys){
            case SYS_BDS: ep=kBdsTimeStart;break;
            case SYS_GAL: ep=kGalTimeStart;break;
        }
        Epoch2Time(ep);

        if(wos<-1E9||1E9<wos) wos=0.0;
        t_.long_time+=(time_t)86400*7*week+(int)wos;
        t_.sec=wos-(int)wos;
        return this;
    }

    cTime cTime::Utc2Gpst() {
        int i;
        cTime t,t0,tg;
        t=*this;

        for(i=0;kUtcLeapSec[i][0]>0;i++){
            if((t.TimeDiff(t0.Epoch2Time(kUtcLeapSec[i])->t_))>=0.0){
                tg=(t+(-kUtcLeapSec[i][6]));
                return tg;
            }
        }
    }

    cTime cTime::Gpst2Utc() {
        int i;
        cTime tu,t0;

        for(i=0; kUtcLeapSec[i][0]>0;i++){
            tu=*this;
            tu+=(kUtcLeapSec[i][6]);
            t0.Epoch2Time(kUtcLeapSec[i]);
            if((TimeDiff(t0.t_))>=0.0){
                *this=tu;
                return *this;
            }
        }
    }

    cTime cTime::Gpst2Bdst(){
        return *this+(-14.0);
    }

    cTime cTime::Bdst2Gpst(){
        return *this+14.0;
    }

    double cTime::Utc2Gmst(double ut1_utc) {
        const double epoch_2000[]={2000,1,1,12,0,0};
        cTime tut,tut0,t2000;
        double ut,t1,t2,t3,gmst0,gmst;

        tut=*this;tut+=ut1_utc;
        ut=tut.Time2Sec(tut0);
        t1=tut0.TimeDiff(t2000.Epoch2Time(epoch_2000)->t_)/86400.0/36525.0;
        t2=t1*t1; t3=t2*t1;
        gmst0=24110.54841+8640184.812866*t1+0.093104*t2-6.2E-6*t3;
        gmst=gmst0+1.002737909350795*ut;

        return fmod(gmst,86400.0)*PI/43200.0; /* 0 <= gmst <= 2*PI */
    }

    cTime* cTime::AdjWeek(cTime t) {
        double dt=TimeDiff(t.t_);
        if(dt<-302400.0) *this+=604800.0;
        if(dt> 302400.0) *this+=(-604800.0);
        return this;
    }

    cTime* cTime::AdjDay(cTime t) {
        double dt=TimeDiff(t.t_);
//        if(dt<-43200.0) return *this+  86400.0;
//        if(dt> 43200.0) return *this+(-86400.0);
    }
#if 0
    cCoord::cCoord() {
        coord_ENU_={0.0,0.0,0.0};
        coord_BLH_={0.0,0.0,0.0};
        coord_XYZ_={0.0,0.0,0.0};
    }

    cCoord::cCoord(const Vector3d& coord, const PPPLib::COORDINATE_TYPE coord_type) {
        if(coord_type==COORD_ENU){
            coord_ENU_=coord;
            Enu2Ned();
        }
        else if(coord_type==COORD_NED){
            coord_NED_=coord;
            Ned2Enu();
        }
        else if(coord_type==COORD_BLH){
            coord_BLH_=coord;
            Blh2Xyz();
        }
        else if(coord_type==COORD_XYZ){
            coord_XYZ_=coord;
            Xyz2Blh();
        }
        else{
            coord_ENU_={0.0,0.0,0.0};
            coord_BLH_={0.0,0.0,0.0};
            coord_XYZ_={0.0,0.0,0.0};
        }
    }

    cCoord::cCoord(const double* coord, const PPPLib::COORDINATE_TYPE coord_type) {
        if(coord_type==COORD_ENU){
            coord_ENU_[0]=coord[0];
            coord_ENU_[1]=coord[1];
            coord_ENU_[2]=coord[2];
            Enu2Ned();
        }
        else if(coord_type==COORD_BLH){
            coord_BLH_[0]=coord[0];
            coord_BLH_[1]=coord[1];
            coord_BLH_[2]=coord[2];
            Blh2Xyz();
        }
        else if(coord_type==COORD_XYZ){
            coord_XYZ_[0]=coord[0];
            coord_XYZ_[1]=coord[1];
            coord_XYZ_[2]=coord[2];
            Xyz2Blh();
        }
        else{
            coord_ENU_={0.0,0.0,0.0};
            coord_BLH_={0.0,0.0,0.0};
            coord_XYZ_={0.0,0.0,0.0};
        }
    }

    cCoord::~cCoord() {

    }

    Vector3d cCoord::GetCoordEnu(cCoord ref_coord){
        Xyz2Enu(ref_coord);
        return coord_ENU_;
    }

    Vector3d cCoord::GetCoordBlh(){
        Xyz2Blh();
        return coord_BLH_;
    }

    Vector3d cCoord::GetCoordXyz(){
        return coord_XYZ_;
    }

    Vector3d cCoord::GetCoordNed(cCoord ref_coord){
        Xyz2Enu(ref_coord);
        Enu2Ned();
        return coord_NED_;
    }

    Matrix3d cCoord::GetCne(const COORDINATE_TYPE type){
        CalcCne(type);
        return Cne_;
    }

    void cCoord::Xyz2Blh() {
        double lon,lat,hgt;
        if(coord_XYZ_[0]>LAT_ACCURACY){
            lon=atan(coord_XYZ_[1]/coord_XYZ_[0]);
        }
        else if(coord_XYZ_[0]<LAT_ACCURACY){
            lon=atan(coord_XYZ_[1]/coord_XYZ_[0])+PI;
        }
        else{
            if(coord_XYZ_[1]>0)
                lon=PI*0.5;
            else
                lon=PI*1.5;
        }
        lat=CalcLat();
        double N=WGS84_EARTH_LONG_RADIUS/sqrt(1.0-WGS84_FIRST_E2*sin(lat)*sin(lat));
        hgt=sqrt(coord_XYZ_[0]*coord_XYZ_[0]+coord_XYZ_[1]*coord_XYZ_[1])*cos(lat)+coord_XYZ_[2]*sin(lat)
                -N*(1.0-WGS84_FIRST_E2*pow(sin(lat),2));

        coord_BLH_[0]=lat;
        coord_BLH_[1]=lon;
        coord_BLH_[2]=hgt;
    }

    void cCoord::Blh2Xyz() {
        const double sin_lat=sin(coord_BLH_[0]);
        const double cos_lat=cos(coord_BLH_[0]);
        const double sin_lon=sin(coord_BLH_[1]);
        const double cos_lon=cos(coord_BLH_[1]);
        double N=WGS84_EARTH_LONG_RADIUS/sqrt(1.0-WGS84_FIRST_E2*sin_lat*sin_lat);

        coord_XYZ_[0]=(N+coord_BLH_[2])*cos_lat*cos_lon;
        coord_XYZ_[1]=(N+coord_BLH_[2])*cos_lat*sin_lon;
        coord_XYZ_[2]=(N*(1.0-WGS84_FIRST_E2)+coord_BLH_[2])*sin_lat;
    }

    void cCoord::Xyz2Enu(cCoord ref_coord) {
        Vector3d temp_xyz;
        temp_xyz[0]=coord_XYZ_[0]-ref_coord.coord_XYZ_[0];
        temp_xyz[1]=coord_XYZ_[1]-ref_coord.coord_XYZ_[1];
        temp_xyz[2]=coord_XYZ_[2]-ref_coord.coord_XYZ_[2];

        double lat=ref_coord.coord_BLH_[0];
        double lon=ref_coord.coord_BLH_[1];

        coord_ENU_[0]=-sin(lon)*temp_xyz[0]+cos(lon)*temp_xyz[1];
        coord_ENU_[1]=-sin(lat)*cos(lon)*temp_xyz[0]-sin(lat)*sin(lon)*temp_xyz[1]+cos(lat)*temp_xyz[2];
        coord_ENU_[2]=cos(lat)*cos(lon)*temp_xyz[0]+cos(lat)*sin(lon)*temp_xyz[1]+sin(lat)*temp_xyz[2];
    }

    void cCoord::Enu2Xyz(cCoord ref_coord) {
        double lat=ref_coord.coord_BLH_[0],lon=ref_coord.coord_BLH_[1];
        double sin_lat=sin(lat);
        double cos_lat=cos(lat);
        double sin_lon=sin(lon);
        double cos_lon=cos(lon);

        double temp1=sin_lat*cos_lon*coord_ENU_[1];
        double temp2=sin_lon*coord_ENU_[0];
        double temp3=cos_lat*cos_lon*coord_ENU_[2];
        coord_XYZ_[0]=ref_coord.coord_XYZ_[0]-temp1-temp2+temp3;

        temp1=sin_lat*sin_lon*coord_ENU_[1];
        temp2=cos_lon*coord_ENU_[0];
        temp3=cos_lat*sin_lon*coord_ENU_[2];
        coord_XYZ_[1]=ref_coord.coord_XYZ_[1]-temp1+temp2+temp3;

        temp1=cos_lat*coord_ENU_[1];
        temp2=sin_lat*coord_ENU_[2];
        coord_XYZ_[2]=ref_coord.coord_XYZ_[2]+temp1+temp2;
    }

    void cCoord::Enu2Ned() {
        coord_NED_[0]= coord_ENU_[1];
        coord_NED_[1]= coord_ENU_[0];
        coord_NED_[2]=-coord_ENU_[2];
    }

    void cCoord::Ned2Enu() {
        coord_ENU_[0]= coord_NED_[1];
        coord_ENU_[1]= coord_NED_[0];
        coord_ENU_[2]=-coord_NED_[2];
    }

    void cCoord::CalcCne(const COORDINATE_TYPE type) {
        double lat=coord_BLH_[0],lon=coord_BLH_[1];
        double sin_lat=sin(lat),cos_lat=cos(lat);
        double sin_lon=sin(lon),cos_lon=cos(lon);
        if(type==COORD_NED){
            Cne_(0,0)=-sin_lat*cos_lon;
            Cne_(0,1)=-sin_lon;
            Cne_(0,2)=-cos_lat*cos_lon;

            Cne_(1,0)=-sin_lat*sin_lon;
            Cne_(1,1)=cos_lon;
            Cne_(1,2)=-cos_lat*sin_lon;

            Cne_(2,0)=cos_lat;
            Cne_(2,1)=0.0;
            Cne_(2,2)=-sin_lat;
        }
        else if(type==COORD_ENU){
            Cne_(0,0)=-sin_lon;
            Cne_(0,1)=cos_lon;
            Cne_(0,2)=0.0;

            Cne_(1,0)=-sin_lat*cos_lon;
            Cne_(1,1)=-sin_lat*sin_lon;
            Cne_(1,2)=cos_lat;

            Cne_(2,0)=cos_lat*cos_lon;
            Cne_(2,1)=cos_lat*sin_lon;
            Cne_(2,2)=sin_lat;
        }
    }

    double cCoord::CalcLat() {
        double temp_lat1=0.0,temp_lat2=0.0;
        double N=0.0;
        temp_lat2=atan(coord_XYZ_[2]/sqrt(coord_XYZ_[0]*coord_XYZ_[0]+coord_XYZ_[1]*coord_XYZ_[1]));
        while(true){
            temp_lat1=temp_lat2;
            N=WGS84_EARTH_LONG_RADIUS/sqrt(1.0-WGS84_FIRST_E2*sin(temp_lat1)*sin(temp_lat1));
            temp_lat2=atan((coord_XYZ_[2]+N*WGS84_FIRST_E2*sin(temp_lat1))/sqrt(coord_XYZ_[0]*coord_XYZ_[0]+coord_XYZ_[1]*coord_XYZ_[1]));
            if(fabs(temp_lat2-temp_lat1)<LAT_ACCURACY){
                return temp_lat2;
            }
        }
    }

    static double CalcLat(Vector3d coord_xyz) {
        double temp_lat1=0.0,temp_lat2=0.0;
        double N=0.0;
        temp_lat2=atan(coord_xyz[2]/sqrt(coord_xyz[0]*coord_xyz[0]+coord_xyz[1]*coord_xyz[1]));
        while(true){
            temp_lat1=temp_lat2;
            N=WGS84_EARTH_LONG_RADIUS/sqrt(1.0-WGS84_FIRST_E2*sin(temp_lat1)*sin(temp_lat1));
            temp_lat2=atan((coord_xyz[2]+N*WGS84_FIRST_E2*sin(temp_lat1))/sqrt(coord_xyz[0]*coord_xyz[0]+coord_xyz[1]*coord_xyz[1]));
            if(fabs(temp_lat2-temp_lat1)<LAT_ACCURACY){
                return temp_lat2;
            }
        }
    }
#endif

    Vector3d Xyz2Blh(Vector3d& coord_xyz){
#if 0
        double lon,lat,hgt;
        if(coord_xyz[0]>LAT_ACCURACY){
            lon=atan(coord_xyz[1]/coord_xyz[0]);
        }
        else if(coord_xyz[0]<LAT_ACCURACY){
            lon=atan(coord_xyz[1]/coord_xyz[0])+PI;
        }
        else{
            if(coord_xyz[1]>0)
                lon=PI*0.5;
            else
                lon=PI*1.5;
        }
        lat=CalcLat(coord_xyz);
        double N=WGS84_EARTH_LONG_RADIUS/sqrt(1.0-WGS84_FIRST_E2*sin(lat)*sin(lat));
        hgt=sqrt(coord_xyz[0]*coord_xyz[0]+coord_xyz[1]*coord_xyz[1])*cos(lat)+coord_xyz[2]*sin(lat)
            -N*(1.0-WGS84_FIRST_E2*pow(sin(lat),2));

        return Vector3d(lat,lon,hgt);
#endif
        double re=WGS84_EARTH_LONG_RADIUS,fe=WGS84_EARTH_OBLATEO;
        Vector3d blh;
        double r2=SQR(coord_xyz[0])+SQR(coord_xyz[1]),e2=fe*(2.0-fe),z,zk,v=re,sinp;

        for(z=coord_xyz[2],zk=0.0;fabs(z-zk)>=1E-4;){
            zk=z;
            sinp=z/sqrt(r2+z*z);
            v=re/sqrt(1.0-e2*sinp*sinp);
            z=coord_xyz[2]+v*e2*sinp;
        }
        blh[0]=r2>1E-12?atan(z/sqrt(r2)):(coord_xyz[2]>0.0?PI/2.0:-PI/2.0);
        blh[1]=r2>1E-12?atan2(coord_xyz[1],coord_xyz[0]):0.0;
        blh[2]=sqrt(r2+z*z)-v;
        return blh;
    }

    Vector3d Blh2Xyz(Vector3d& coord_blh){
        const double sin_lat=sin(coord_blh[0]);
        const double cos_lat=cos(coord_blh[0]);
        const double sin_lon=sin(coord_blh[1]);
        const double cos_lon=cos(coord_blh[1]);
        double N=WGS84_EARTH_LONG_RADIUS/sqrt(1.0-WGS84_FIRST_E2*sin_lat*sin_lat);

        double x=(N+coord_blh[2])*cos_lat*cos_lon;
        double y=(N+coord_blh[2])*cos_lat*sin_lon;
        double z=(N*(1.0-WGS84_FIRST_E2)+coord_blh[2])*sin_lat;

        return Vector3d(x,y,z);
    }

    Vector3d Xyz2Enu(Vector3d& coord_blh,Vector3d& coord_xyz){
#if 0
        Vector3d temp_xyz;
        temp_xyz[0]=coord_xyz[0]-ref_xyz[0];
        temp_xyz[1]=coord_xyz[1]-ref_xyz[1];
        temp_xyz[2]=coord_xyz[2]-ref_xyz[2];

        Vector3d ref_blh=Xyz2Blh(ref_xyz);

        double lat=ref_blh[0];
        double lon=ref_blh[1];

        Vector3d enu(0,0,0);

        enu[0]=-sin(lon)*temp_xyz[0]+cos(lon)*temp_xyz[1];
        enu[1]=-sin(lat)*cos(lon)*temp_xyz[0]-sin(lat)*sin(lon)*temp_xyz[1]+cos(lat)*temp_xyz[2];
        enu[2]=cos(lat)*cos(lon)*temp_xyz[0]+cos(lat)*sin(lon)*temp_xyz[1]+sin(lat)*temp_xyz[2];

        return enu;
#endif
        Matrix3d Cen;
        Cen=CalcCen(coord_blh,COORD_ENU);
        return Cen*coord_xyz;
    }

    // coord in enu frame to ECEF frame coord
    Vector3d Enu2Xyz(Vector3d& coord_blh,Vector3d& coord_enu){
#if 0
        Vector3d ref_blh=Xyz2Blh(ref_xyz);
        double lat=ref_blh[0],lon=ref_blh[1];
        double sin_lat=sin(lat);
        double cos_lat=cos(lat);
        double sin_lon=sin(lon);
        double cos_lon=cos(lon);

        double temp1=sin_lat*cos_lon*coord_enu[1];
        double temp2=sin_lon*coord_enu[0];
        double temp3=cos_lat*cos_lon*coord_enu[2];

        Vector3d coord_xyz;
        coord_xyz[0]=ref_xyz[0]-temp1-temp2+temp3;

        temp1=sin_lat*sin_lon*coord_enu[1];
        temp2=cos_lon*coord_enu[0];
        temp3=cos_lat*sin_lon*coord_enu[2];
        coord_xyz[1]=ref_xyz[1]-temp1+temp2+temp3;

        temp1=cos_lat*coord_enu[1];
        temp2=sin_lat*coord_enu[2];
        coord_xyz[2]=ref_xyz[2]+temp1+temp2;

        return coord_xyz;
#endif
        Matrix3d Cne;
        Cne=CalcCen(coord_blh,COORD_ENU).transpose();
        return Cne*coord_enu;
    }

    Vector3d Enu2Ned(Vector3d& coord_enu){
        Vector3d ned;
        ned[0]=coord_enu[1];
        ned[1]=coord_enu[0];
        ned[2]=-coord_enu[2];
        return ned;
    }

    Vector3d Ned2Enu(Vector3d& coord_ned){
        Vector3d enu;
        enu[0]=coord_ned[1];
        enu[1]=coord_ned[0];
        enu[2]=-coord_ned[2];

        return enu;
    }

    // n refer to enu
    Matrix3d CalcCen(Vector3d& coord_blh,COORDINATE_TYPE lf_type){
        double lat=coord_blh[0],lon=coord_blh[1];
        double sin_lat=sin(lat),cos_lat=cos(lat);
        double sin_lon=sin(lon),cos_lon=cos(lon);
        Matrix3d Cen;
#if 0
        Cne(0,0)=-sin_lat*cos_lon;
        Cne(0,1)=-sin_lon;
        Cne(0,2)=-cos_lat*cos_lon;

        Cne(1,0)=-sin_lat*sin_lon;
        Cne(1,1)=cos_lon;
        Cne(1,2)=-cos_lat*sin_lon;

        Cne(2,0)=cos_lat;
        Cne(2,1)=0.0;
        Cne(2,2)=-sin_lat;
#endif
        if(lf_type==COORD_ENU){
            Cen(0,0)=-sin_lon;
            Cen(0,1)=cos_lon;
            Cen(0,2)=0.0;

            Cen(1,0)=-sin_lat*cos_lon;
            Cen(1,1)=-sin_lat*sin_lon;
            Cen(1,2)=cos_lat;

            Cen(2,0)=cos_lat*cos_lon;
            Cen(2,1)=cos_lat*sin_lon;
            Cen(2,2)=sin_lat;
        }
        else if(lf_type==COORD_NED){
            Cen(0,0)=-sin_lat*cos_lon;
            Cen(0,1)=-sin_lat*sin_lon;
            Cen(0,2)=cos_lat;

            Cen(1,0)=-sin_lon;
            Cen(1,1)=cos_lon;
            Cen(1,2)=0.0;

            Cen(2,0)=-cos_lat*cos_lon;
            Cen(2,1)=-cos_lat*sin_lon;
            Cen(2,2)=-sin_lat;
        }

        return Cen;
    }

    cParSetting::cParSetting(){}

    cParSetting::cParSetting(tPPPLibConf conf) {PPPLibC_=conf;}

    cParSetting::~cParSetting() {}

    int cParSetting::GetGnssUsedFrqs() {
        if(PPPLibC_.gnssC.ion_opt==ION_IF){
            if(PPPLibC_.gnssC.frq_opt==FRQ_SINGLE||PPPLibC_.gnssC.frq_opt==FRQ_DUAL) return 1;
            else if(PPPLibC_.gnssC.frq_opt==FRQ_TRIPLE) return 3;
        }
        else if(PPPLibC_.gnssC.ion_opt==ION_IF_DUAL) return 2;
        else return PPPLibC_.gnssC.frq_opt+1;
    }

    int cParSetting::GetNumObsType() {
        if(PPPLibC_.mode==MODE_SPP||PPPLibC_.mode_opt==MODE_OPT_SPP) return 1;
        else return 1+1;
    }

    int cParSetting::GetInsTransParNum(tPPPLibConf C) {
        PPPLibC_=C;
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg();
    }

    int cParSetting::GetPPPLibPar(tPPPLibConf C) {
        PPPLibC_=C;
        int a=NumPos();
        int b=NumVel();
        int c=NumAtt();
        int d=NumClk();
        int e=NumBa();
        int f=NumBg();
        int g=NumDcb();
        int h=NumIfb();
        int i=NumGloIfcb();
        int j=NumTrp();
        int k=NumIon();
        int l=NumAmb();
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()
                       +NumClk()+NumDcb()+NumIfb()+NumGloIfcb()
                       +NumTrp()+NumIon()+NumAmb();
    }

    int cParSetting::GetRealFixParNum(tPPPLibConf C) {
        PPPLibC_=C;
        int a=NumPos();
        int b=NumVel();
        int c=NumAtt();
        int d=NumClk();
        int e=NumBa();
        int f=NumBg();
        int g=NumDcb();
        int h=NumIfb();
        int i=NumGloIfcb();
        int j=NumTrp();
        int k=NumIon();
        int l=NumAmb();
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()
               +NumClk()+NumDcb()+NumIfb()+NumGloIfcb()
               +NumTrp()+NumIon();
    }

    int cParSetting::NumPos() {
        return 3;
    }

    int cParSetting::NumVel() {
        if(PPPLibC_.dynamic||PPPLibC_.mode>MODE_INS) return 3;
        else return 0;
    }

    int cParSetting::NumAtt() {
        if(PPPLibC_.mode>MODE_INS) return 3;
        else return 0;
    }

    int cParSetting::NumBa() {
        if(PPPLibC_.mode>MODE_INS) return 3;
        else return 0;
    }

    int cParSetting::NumBg() {
        if(PPPLibC_.mode>MODE_INS) return 3;
        else return 0;
    }

    int cParSetting::NumClPar() {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg();
    }

    int cParSetting::NumClk() {
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC||PPPLibC_.mode==MODE_PPK||PPPLibC_.mode_opt==MODE_OPT_PPK) return 0;
        return NSYS;
    }

    int cParSetting::NumClkDrift() {
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC) return 0;
        if(PPPLibC_.gnssC.use_doppler) return NSYS;
    }

    int cParSetting::NumDcb() {
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC) return 0;
        if(PPPLibC_.gnssC.ion_opt==ION_CONST) return NSYS;
        return 0;
    }

    int cParSetting::NumIfb() {
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC) return 0;
        if(PPPLibC_.gnssC.frq_opt>=FRQ_TRIPLE&&PPPLibC_.mode<MODE_INS) return NSYS;
        return 0;
    }

    int cParSetting::NumGloIfcb() {
        if(!(PPPLibC_.gnssC.nav_sys&SYS_GLO)) return 0;
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC||PPPLibC_.mode==MODE_SPP||PPPLibC_.gnssC.glo_ifcb_opt==GLO_IFCB_OFF) return 0;
        if(PPPLibC_.gnssC.glo_ifcb_opt==GLO_IFCB_LNF) return 1;
        else if(PPPLibC_.gnssC.glo_ifcb_opt==GLO_IFCB_QUAD) return 2;
        else if(PPPLibC_.gnssC.glo_ifcb_opt==GLO_IFCB_1SAT) return NUM_GLO_SAT;
        else if(PPPLibC_.gnssC.glo_ifcb_opt==GLO_IFCB_1FRQ) return 13;
        else return 0;
    }

    int cParSetting::NumTrp() {
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC) return 0;
        if(PPPLibC_.gnssC.trp_opt<=TRP_SAAS) return 0;
        else if(PPPLibC_.gnssC.trp_opt==TRP_EST_WET) return 1;
        else if(PPPLibC_.gnssC.trp_opt==TRP_EST_GRAD) return 1+2;
    }

    int cParSetting::NumIon() {
        if(PPPLibC_.mode==MODE_INS||PPPLibC_.mode==MODE_IGLC) return 0;
        if(PPPLibC_.gnssC.ion_opt<=ION_IF_DUAL) return 0;
        else if(PPPLibC_.gnssC.ion_opt<=ION_CONST) return MAX_SAT_NUM;
    }

    int cParSetting::NumAmb() {
        if(PPPLibC_.mode==MODE_IGLC) return 0;
        if(PPPLibC_.mode==MODE_PPP||PPPLibC_.mode==MODE_PPK||PPPLibC_.mode_opt==MODE_OPT_PPP||PPPLibC_.mode_opt==MODE_OPT_PPK)
            return MAX_SAT_NUM*GetGnssUsedFrqs();
        else return 0;
    }

    int cParSetting::IndexPos() {
        return 0;
    }

    int cParSetting::IndexVel() {
        if(NumVel()>0) return NumPos();
        else return -1;
    }

    int cParSetting::IndexAtt() {
        if(NumAtt()>0) return NumPos()+NumVel();
    }

    int cParSetting::IndexBa() {
        if(NumBa()>0) return NumPos()+NumVel()+NumAtt();
    }

    int cParSetting::IndexBg() {
        if(NumBg()>0) return NumPos()+NumVel()+NumAtt()+NumBa();
    }

    int cParSetting::IndexClk(int sys_index) {
        if(NumClk()>0) return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+sys_index;
        else return -1;
    }

    int cParSetting::IndexClkDritf(int sys_index) {
        if(NumClkDrift()) return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+sys_index;
    }

    int cParSetting::IndexDcb(int sys_index) {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+NumClkDrift()+sys_index;
    }

    int cParSetting::IndexIfb(int sys_index) {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+NumClkDrift()+NumDcb()+sys_index;
    }

    int cParSetting::IndexGloIfcb(int i) {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+NumClkDrift()+NumIfb()+i-1;
    }

    int cParSetting::IndexTrp() {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+NumClkDrift()+NumIfb()+NumGloIfcb();
    }

    int cParSetting::IndexIon(int sat_no) {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+NumClkDrift()+NumIfb()+NumGloIfcb()+NumTrp()+sat_no-1;
    }

    int cParSetting::IndexAmb(int f, int sat_no) {
        return NumPos()+NumVel()+NumAtt()+NumBa()+NumBg()+NumClk()+NumClkDrift()+NumIfb()+NumGloIfcb()+NumTrp()+NumIon()
               +f*MAX_SAT_NUM+sat_no-1;
    }

    Config::Ptr_ Config::config_info_(new Config());

    bool Config::Open(std::string config_file) {
        ifstream inf(config_file,ios::in);
        if(!inf){
            cout<<"Configure file path error";
            return false;
        }
        while(!inf.eof()){
            string line,key;
            getline(inf,line);
            if(line.size()==0) continue;
            else if(line.substr(0,1)=="#"||line.substr(0,1)=="*"||line.substr(0,1)=="!") continue;
            auto data=TextSplit(line,";");
            key=StrTrim(data[0]);
            transform(key.begin(),key.end(),key.begin(),::tolower);
            if(data.size()>1){
                storage_[key]=StrTrim(data[1]);
            }
        }
        return true;
    }

    Config::Ptr_ Config::GetInstance() {
        return config_info_;
    }
}






