//
// Created by chenc on 2020/7/13.
//

#include "ReadFiles.h"

#define MAX_CSV_COLS  30

namespace PPPLib{

    cMatchFile::cMatchFile() {}

    cMatchFile::~cMatchFile() {}

    void cMatchFile::InitMatchFile(tPPPLibConf &C,char sep) {
        sep_=sep;
        C_=&C;

        year_=(int)C_->prc_date.GetEpoch()[0];
        yy_=year_>2000.0?Round(year_-2000.0):Round(year_-1900.0);
        C_->prc_date.Time2Gpst(&week_, nullptr,SYS_GPS);
        doy_=(int)C_->prc_date.Time2Doy();
        wod_=(int)(C_->prc_date.Time2Gpst(&week_, nullptr,SYS_GPS)/86400.0);
    }

    int cMatchFile::MatchFileAuto() {
        if(!MatchProd()){
            return 0;
        }
        if(!MatchPrec()){
            return 0;
        }
        if(!MatchCmn()){
            return 0;
        }
        if(!MatchOut()){
            return 0;
        }

        return 1;
    }

    int cMatchFile::MatchProd() {
        string prods_dir;
        char f[1024]={'\0'};

        if(!C_->use_custom_dir){
            sprintf(f,"%s%c%d%c%d%c%03d%c%s",C_->data_dir.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,"prods");
            prods_dir=f;f[0]='\0';
        }else prods_dir=C_->data_dir;

        LOG(DEBUG)<<"PRODUCT DIR: "<<prods_dir;

        sprintf(f,"%s%cbrdm%03d0.%02dp",prods_dir.c_str(),sep_,doy_,yy_);
        if((access(f,0))==-1){
            LOG(ERROR)<<"BRDM FILE NO FOUND: "<<f;
            return 0;
        }
        C_->fileC.brd=f;
        LOG(DEBUG)<<"BROADCAST EPHEMERICS FILE: "<<C_->fileC.brd;
        f[0]='\0';

        int ion=C_->gnssC.ion_opt==ION_TEC||C_->gnssC.ion_opt==ION_CONST;
        if(ion){
            sprintf(f,"%s%cCODG%03d0.%02dI",prods_dir.c_str(),sep_,doy_,yy_);
            if((access(f,0))==-1){
                LOG(ERROR)<<"IONOSPHERE PRODUCT NO FOUND: "<<f;
                return 0;
            }
            C_->fileC.gim=f;f[0]='\0';
            LOG(DEBUG)<<"GIM FILE: "<<C_->fileC.gim;
        }

        sprintf(f,"%s%cCAS0MGXRAP_%d%03d0000_01D_01D_DCB.BSX",prods_dir.c_str(),sep_,year_,doy_);
        if(access(f,0)==-1){
            LOG(ERROR)<<"CAS DCB NO FOUND: "<<f;
            return 0;
        }
        C_->fileC.cbias=f;f[0]='\0';
        LOG(DEBUG)<<"CAS DCB FILE: "<<C_->fileC.cbias;

        return 1;
    }

    int cMatchFile::MatchPrec() {
        string pres_dir;
        int wod1,wk1;
        int wod2,wk2;
        char f[1024]={'\0'};

        if(!C_->use_custom_dir){
            sprintf(f,"%s%c%d%c%d%c%s",C_->data_dir.c_str(),sep_,year_,sep_,week_,sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str());
            pres_dir=f;f[0]='\0';
        }else pres_dir=C_->data_dir;

        LOG(DEBUG)<<"PRECISE PRODUCTS DIR: "<<pres_dir;
        if(C_->gnssC.eph_opt==EPH_PRE){
            sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),week_,wod_,".sp3");
            if((access(f,0))==-1){
                LOG(ERROR)<<" PRECISE ORBIT FILE NO FOUND: "<<f;
                return 0;
            }
            C_->fileC.sp3[1]=f;f[0]='\0';
            LOG(DEBUG)<<"ON THE DAY SP3 FILE: "<<C_->fileC.sp3[1];

            if(wod_==0) {
                wod1=6;wk1=week_-1;
            }
            else {
                wod1=wod_-1;wk1=week_;
            }
            sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),wk1,wod1,".sp3");
            if((access(f,0))==-1){
                LOG(DEBUG)<<"BEFORE DAY SP3 FILE NO FOUND: "<<f;
            }
            else{
                C_->fileC.sp3[0]=f;f[0]='\0';
                LOG(DEBUG)<<"BEFORE DAY SP3 FILE: "<<C_->fileC.sp3[0];
            }

            if(wod_==6) {
                wod2=0;wk2=week_+1;
            }
            else {
                wod2=wod_+1;wk2=week_;
            }
            sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),wk2,wod2,".sp3");
            if((access(f,0))==-1){
                LOG(DEBUG)<<"AFTER DAY SP3 NO FOUND: "<<f;
            }
            else{
                C_->fileC.sp3[2]=f;f[0]='\0';
                LOG(DEBUG)<<"AFTER DAY SP3 FILE: "<<C_->fileC.sp3[2];
            }

            sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),week_,wod_,".clk");
            if((access(f,0))==-1){
                LOG(ERROR)<<"PRECISE CLOCK FILE NO FOUND: "<<f;
                return 0;
            }
            C_->fileC.clk=f;f[0]='\0';
            LOG(DEBUG)<<"PRECISE CLOCK FILE: "<<C_->fileC.clk;
        }
        if(C_->gnssC.tid_opt>TID_OFF){
            sprintf(f,"%s%c%s%d%d%s",pres_dir.c_str(),sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),week_,wod_,".erp");
            if((access(f,0))!=-1){
                C_->fileC.erp=f;f[0]='\0';
                LOG(DEBUG)<<"ERP FILE: "<<C_->fileC.erp;
            }
            else{
                sprintf(f,"%s%c%d%c%digs%2dP%d.erp",C_->data_dir.c_str(),sep_,year_,sep_,week_,yy_,week_);
                if((access(f,0))==-1){
                    LOG(DEBUG)<<"ERP FILE NO FOUND "<<f;
                }
                C_->fileC.erp=f;
                LOG(DEBUG)<<"ERP FILE: "<<C_->fileC.erp;
            }
        }
        return 1;
    }

    int cMatchFile::MatchCmn() {
        string cmn_prods;
        char f[1024]={'\0'};

        if(!C_->use_custom_dir){
            sprintf(f,"%s%ccmnprods",C_->data_dir.c_str(),sep_);
            cmn_prods=f;f[0]='\0';
        }else cmn_prods=C_->data_dir;

        LOG(DEBUG)<<"COMMON PRODUCTS DIR: "<<cmn_prods;
        if(C_->gnssC.rec_ant||C_->gnssC.sat_pcv)
        {
            sprintf(f,"%s%cigs14_2097.atx",cmn_prods.c_str(),sep_);
            if((access(f,0))==-1){
                LOG(DEBUG)<<"ATX FILE NO FOUND: "<<f;
            }
            C_->fileC.atx=f;f[0]='\0';
            LOG(DEBUG)<<"ATX FILE: "<<C_->fileC.atx;
        }
        if(C_->gnssC.tid_opt){
            sprintf(f,"%s%cocnload.blq",cmn_prods.c_str(),sep_);
            if((access(f,0))==-1){
                LOG(DEBUG)<<"BLQ FILE NO FOUND: "<<f;
            }
            C_->fileC.blq=f;
            LOG(DEBUG)<<"BLQ FILE: "<<C_->fileC.blq;
        }
        return 1;
    }

    int cMatchFile::MatchOut() {
        string sys_str;
        string out_dir,mode_dir,par_dir,name;
        char f[1024]={'\0'};

        if(C_->gnssC.nav_sys&SYS_GPS){
            sys_str+="G";
        }
        if(C_->gnssC.nav_sys&SYS_BDS){
            sys_str+="B";
        }
        if(C_->gnssC.nav_sys&SYS_GAL){
            sys_str+="E";
        }
        if(C_->gnssC.nav_sys&SYS_GLO){
            sys_str+="R";
        }
        if(C_->gnssC.nav_sys&SYS_QZS){
            sys_str+="J";
        }
        if(C_->use_custom_dir){
            sprintf(f,"%s%cresult_%s%c",C_->data_dir.c_str(),sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),sep_);
        }else sprintf(f,"%s%c%d%c%d%c%03d%cresult_%s%c",C_->data_dir.c_str(),sep_,year_,sep_,week_,sep_,doy_,sep_,kGnssAcStr[C_->gnssC.ac_opt].c_str(),sep_);
        out_dir=f;
        CreateDir(out_dir.c_str());
        sprintf(f,"%s%c%s_%s%c",out_dir.c_str(),sep_,kPpplibModeStr[C_->mode].c_str(),kPpplibModeOptStr[C_->mode_opt].c_str(),sep_);
        mode_dir=f;f[0]='\0';
        CreateDir(mode_dir.c_str());
        sprintf(f,"%s%c%s_%s_%s%c",mode_dir.c_str(),sep_,sys_str.c_str(),kGnssFrqStr[C_->gnssC.frq_opt].c_str(),kIonOptStr[C_->gnssC.ion_opt].c_str(),sep_);
        par_dir=f;f[0]='\0';
        CreateDir(par_dir.c_str());
        sprintf(f,"%s%s.pos",par_dir.c_str(),C_->site_name.c_str());
        C_->fileC.sol=f;
        f[0]='\0';
        if(C_->solC.out_stat){
            sprintf(f,"%s%s.pos.stat",par_dir.c_str(),C_->site_name.c_str());
        }
        C_->fileC.sol_stat=f;
        LOG(DEBUG)<<"SOLUTION FILE: "<<C_->fileC.sol;

        return 1;
    }

    cReadFile::cReadFile() { }

    cReadFile::cReadFile(string file_path):file_(file_path) {}

    cReadFile::~cReadFile() {}

    bool cReadFile::OpenFile() {
        inf_.open(file_);
        if(!inf_.is_open()){
            LOG(WARNING)<<"File open error: "<<file_;
            return false;
        }
        return true;
    }

    void cReadFile::CloseFile() {
        if(inf_.is_open()){
            inf_.close();
        }
    }

    bool cReadFile::ReadHead() {}

    bool cReadFile::Reading() {return true;}

    cReadImu::cReadImu() {}

    cReadImu::cReadImu(string file_path){file_=file_path;}

    cReadImu::~cReadImu() {}

    bool cReadImu::DecodeImu() {
        if(imu_data_.imu_type_==IMU_NOVTEL_CPT){
            double a_scale=1.52587890625E-06;
            double g_scale=1.085069444444444E-07;
            DecodeNovatel(a_scale, g_scale);
        }else if(imu_data_.imu_type_==IMU_NOVTEL_A1){
            double a_scale=1.52587890625E-06;
            double g_scale=1.085069444444444E-07;
            DecodeNovatel(a_scale, g_scale);
        }
        else if(imu_data_.imu_type_==IMU_MTI_CSV){
            DecodeCsv(TIME_FMT_WS);
        }
        else{
            LOG(ERROR)<<"Unsupport imu type";
            return false;
        }
        return true;
    }

    bool cReadImu::DecodeNovatel(double a_scale, double g_scale) {
        char *seps=(char *)",;*";
        char *token;
        int i,gps_week;
        tImuDataUnit imu_data={0};
        Matrix3d Crf;
        Crf<<0,1,0,1,0,0,0,0,-1;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            token=strtok((char *)line_str_.c_str(),seps);
            if(!strcmp(token,"%RAWIMUSA")){
                for(i=0;i<3;i++) token=strtok(NULL,seps);
                gps_week=atoi(token);
                token=strtok(NULL,seps);
                imu_data.t_tag.Gpst2Time(gps_week, atof(token),SYS_GPS);

                if(imu_data_.ts_.t_.long_time!=0.0){
                    if(imu_data.t_tag.TimeDiff(imu_data_.ts_.t_)<0) continue;
                }
                if(imu_data_.te_.t_.long_time!=0.0){
                    if(imu_data.t_tag.TimeDiff(imu_data_.te_.t_)>0) continue;
                }

                for(i=0;i<2;i++) token=strtok(NULL,seps);
                imu_data.acce[2]=atof(token)*a_scale;
                token=strtok(NULL,seps);
                imu_data.acce[1]=-atof(token)*a_scale;
                token=strtok(NULL,seps);
                imu_data.acce[0]=atof(token)*a_scale;

                token=strtok(NULL,seps);
                imu_data.gyro[2]=atof(token)*g_scale;
                token=strtok(NULL,seps);
                imu_data.gyro[1]=atof(token)*g_scale;
                token=strtok(NULL,seps);
                imu_data.gyro[0]=atof(token)*g_scale; //RFU,increment

                AdjustImuData(imu_data,imu_data_.imu_coord_type_,imu_data_.data_format_,imu_data_.gyro_format_,1.0/imu_data_.hz_);
#if 0
                imu_data.gyro=Crf*imu_data.gyro;   //RFU->FRD
                imu_data.acce=Crf*imu_data.acce;
#endif
                imu_data_.data_.push_back(imu_data);
            }
        }
        if(imu_data_.data_.empty()) return false;
        return true;
    }

    bool cReadImu::DecodeCsv(TIME_FORMAT time_f) {
        double data[MAX_CSV_COLS]={0};
        string part_str;
        tImuDataUnit imu_data={{nullptr}};
        int i,week;
        double sow;
        string t_str;
        cTime t_tag;

        getline(inf_,line_str_);
        while(getline(inf_,line_str_)&&!inf_.eof()){

            if(!line_str_.substr(0,3).compare("EOF")) break;

            istringstream read_str(line_str_);

            for(auto &j:data) j=0.0;
            for(i=0;i<cols_&&cols_<MAX_CSV_COLS;i++){

                getline(read_str,part_str,sep_);
                if(i==idx_t_){
                    if(time_f==TIME_FMT_STR){
                        t_str=part_str;
                        t_tag.Str2Time(t_str);
                        continue;
                    }
                    else{
                        week=atof(part_str.c_str());
                        getline(read_str,part_str,sep_);
                        sow=atof(part_str.c_str());
                        t_tag.Gpst2Time(week,sow,SYS_GPS);
                        i++;
                        continue;
                    }
                }
                else if(i>idx_t_){
                    data[i]=atof(part_str.c_str());
                }
            }

            imu_data.t_tag=t_tag;
            imu_data.acce[0]=data[idx_ax_];
            imu_data.acce[1]=data[idx_ay_];
            imu_data.acce[2]=data[idx_az_];
            imu_data.gyro[0]=data[idx_gx_];
            imu_data.gyro[1]=data[idx_gy_];
            imu_data.gyro[2]=data[idx_gz_];
            imu_data_.data_.push_back(imu_data);
        }
    }

    cImuData* cReadImu::GetImus() {return &imu_data_;}

    bool cReadImu::SetImu(tPPPLibConf C) {
        if(C.insC.imu_type==IMU_UNKNOW){
            LOG(ERROR)<<"Unknow imu type";
            return false;
        }
        imu_data_.SetImu(C.insC);
    }

    void cReadImu::SetDataIdx(int t, int ax, int ay, int az, int gx, int gy, int gz,const char sep, int cols) {
        idx_t_=t;
        idx_ax_=ax;
        idx_ay_=ay;
        idx_az_=az;
        idx_gx_=gx;
        idx_gy_=gy;
        idx_gz_=gz;
        sep_=sep;
        cols_=cols;
    }

    void cReadImu::SetImuTimeSpan(PPPLib::cTime *ts, PPPLib::cTime *te) {
        if(ts->TimeDiff(te->t_)>=0){
            LOG(WARNING)<<"Wrong time interval setting";
            return;
        }
        else{
            imu_data_.SetTimeSpan(ts,te);
        }
    }

    bool cReadImu::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open imu file error: "<<file_;
            return false;
        }

        DecodeImu();

        if(OpenFile()) CloseFile();
        return true;
    }

    void cReadImu::OutImu(string out_path) {
        ofstream fout(out_path);
        string sep="   ";
        tImuDataUnit imu_data={0};

        for(int i=0;i<imu_data_.data_.size();i++){
            imu_data=imu_data_.data_[i];
            fout<<imu_data.t_tag.GetTimeStr(3)<<sep<<imu_data.acce[0]<<sep<<imu_data.acce[1]<<sep<<imu_data.acce[2];
            fout<<sep<<imu_data.gyro[0]<<sep<<imu_data.gyro[1]<<sep<<imu_data.gyro[2]<<endl;
        }
    }

    cReadPos::cReadPos() {}

    cReadPos::cReadPos(string file_path) {file_=file_path;}

    cReadPos::~cReadPos() {}

    bool cReadPos::DecodePos() {
        string line_str;

        while(getline(inf_,line_str)&&!inf_.eof()){
            if(line_str.size()==0) continue;
        }
    }

    bool cReadPos::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open imu file error: "<<file_;
            return false;
        }

        DecodePos();

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadRnx::cReadRnx(){sys_mask_=SYS_ALL;}

    cReadRnx::cReadRnx(string file_path) {file_=file_path;sys_mask_=SYS_ALL;}

    cReadRnx::~cReadRnx() {}

    void cReadRnx::SetGnssSysMask(int mask) {
        sys_mask_=mask;
    }

    bool cReadRnx::ReadRnxHead() {

        while(getline(inf_,line_str_)&&!inf_.eof()){

            if(line_str_.find("RINEX VERSION / TYPE")!=string::npos){
                Str2Double(line_str_.substr(0,9),rnx_ver_);
                rnx_type_=line_str_.substr(20,1);
                switch(line_str_[40]){
                    case ' ':
                    case 'G': sat_sys_=SYS_GPS; time_sys_=TIME_GPS;break;
                    case 'C': sat_sys_=SYS_BDS; time_sys_=TIME_BDS;break;
                    case 'E': sat_sys_=SYS_GAL; time_sys_=TIME_GAL;break;
                    case 'R': sat_sys_=SYS_GLO; time_sys_=TIME_utc;break;
                    case 'J': sat_sys_=SYS_QZS; time_sys_=TIME_QZS;break;
                    case 'M': sat_sys_=SYS_NONE;time_sys_=TIME_GPS;break;
                    default: break;
                }
                break;
            }
        }

        if(rnx_type_.empty()){
            LOG(ERROR)<<"Rinex file type error";
            return false;
        }
        return true;
    }

    cGnssSignal::cGnssSignal() {}

    cGnssSignal::~cGnssSignal() {}

    tGnssSignal * cGnssSignal::GetGnssSignal() {return &signal_;}

    unsigned char cGnssSignal::Signal2Code(string signal, int *frq, int sat_sys) {
        int i;
        if(frq) *frq=0;
        for(i=0;i<MAX_GNSS_CODE_TYPE;i++){
            if(kGnssSignalCodes[i]!=signal) continue;
            if(frq){
                switch(sat_sys){
                    case SYS_GPS: *frq=kGpsFreqBand[i];break;
                    case SYS_BDS: *frq=kBdsFreqBand[i];break;
                    case SYS_GAL: *frq=kGalFreqBand[i];break;
                    case SYS_GLO: *frq=kGloFreqBand[i];break;
                    case SYS_QZS: *frq=kQzsFreqBand[i];break;
                    case SYS_IRN: *frq=0;break;
                    default:
                        *frq=0;break;
                }
            }
            return (unsigned char)i;
        }
        return GNSS_CODE_NONE;
    }

    string cGnssSignal::Code2Signal(unsigned char code, int *frq, int sat_sys) {
        if(frq) *frq=0;
        if(code<=GNSS_CODE_NONE||MAX_GNSS_CODE_TYPE<code) return "";
        if(frq){
            switch(sat_sys){
                case SYS_GPS: *frq=kGpsFreqBand[code];break;
                case SYS_BDS: *frq=kBdsFreqBand[code];break;
                case SYS_GAL: *frq=kGalFreqBand[code];break;
                case SYS_GLO: *frq=kGloFreqBand[code];break;
                case SYS_QZS: *frq=kQzsFreqBand[code];break;
                case SYS_IRN: *frq=0;break;
                default:
                    *frq=0;break;
            }
        }
        return kGnssSignalCodes[code];
    }

    int cGnssSignal::GetCodePri(int sat_sys, unsigned char code) {
        if(code==GNSS_CODE_NONE) return 0;
        size_t str_pos;
        string signal;
        int i,f;

        signal=Code2Signal(code,&f,sat_sys);
        switch(sat_sys){
            case SYS_GPS: i=SYS_INDEX_GPS;break;
            case SYS_BDS: i=SYS_INDEX_BDS;break;
            case SYS_GAL: i=SYS_INDEX_GAL;break;
            case SYS_GLO: i=SYS_INDEX_GLO;break;
            case SYS_QZS: i=SYS_INDEX_QZS;break;
            case SYS_IRN: i=0;break;
            default:
                i=0;break;
        }
        return ((str_pos=kGnssSignalPriors[i][f-1].find(signal[1]))!=string::npos)?14-(int)str_pos:0;
    }

    void cGnssSignal::GnssSignalIndex(int sat_sys, string *obs_code_type) {
        size_t p;
        string s;
        int i,j,num_signal;
        int k_code=-1,k_phase=-1,k_doppler=-1,k_snr=-1;

        for(i=num_signal=0;obs_code_type[i][0];i++,num_signal++){
            signal_.code[i]=Signal2Code(obs_code_type[i].substr(1),signal_.frq+i,sat_sys);
            signal_.type[i]=((p=kGnssObsCode.find(obs_code_type[i][0]))!=string::npos)?(int)p:0;
            signal_.pri[i]=GetCodePri(sat_sys,signal_.code[i]);
            signal_.pos[i]=-1;
        }

        for(i=0;i<MAX_GNSS_FRQ_NUM;i++){
            for(j=0;j<num_signal;j++){
                if(signal_.type[j]==GNSS_OBS_CODE){
                    if(signal_.frq[j]==i+1&&signal_.pri[j]&&
                            (k_code<0||signal_.pri[j]>=signal_.pri[k_code])){
                        k_code=j;
                    }
                }
                else if(signal_.type[j]==GNSS_OBS_PHASE){
                    if(signal_.frq[j]==i+1&&signal_.pri[j]&&
                            (k_phase<0||signal_.pri[j]>=signal_.pri[k_phase])){
                        k_phase=j;
                    }
                }
                else if(signal_.type[j]==GNSS_OBS_DOPPLER){
                    if(signal_.frq[j]==i+1&&signal_.pri[j]&&
                            (k_doppler<0||signal_.pri[j]>=signal_.pri[k_doppler])){
                        k_doppler=j;
                    }
                }
                else if(signal_.type[j]==GNSS_OBS_SNR){
                    if(signal_.frq[j]==i+1&&signal_.pri[j]&&
                            (k_snr<0||signal_.pri[j]>=signal_.pri[k_snr])){
                        k_snr=j;
                    }
                }
            }

            if(k_code>=0)    signal_.pos[k_code]=i;
            if(k_phase>=0)   signal_.pos[k_phase]=i;
            if(k_doppler>=0) signal_.pos[k_doppler]=i;
            if(k_snr>=0)     signal_.pos[k_snr]=i;

            k_code=k_phase=k_doppler=k_snr=-1;
        }

        int k;
        for(i=0;i<GNSS_NUM_EXOBS;i++){
            for(j=0;j<num_signal;j++){
                if(signal_.code[j]&&signal_.pri[j]&&signal_.pos[j]<0) break;
             }
            if(j>=num_signal) break;
            for(k=0;k<num_signal;k++){
                if(signal_.code[k]==signal_.code[j]) signal_.pos[k]=MAX_GNSS_FRQ_NUM+i;
            }
        }

        signal_.n=num_signal;
    }

    cReadGnssObs::cReadGnssObs() {}

    cReadGnssObs::cReadGnssObs(string file_path,PPPLib::tNav& nav,cGnssObs& obss,RECEIVER_INDEX rcv) {
        file_=file_path;
        gnss_data_=&obss;
        gnss_data_->SetRcvIdx(rcv);
        nav_=&nav;
    }

    cReadGnssObs::~cReadGnssObs() {}

    bool cReadGnssObs::CmpEpochSatData(const PPPLib::tSatObsUnit &p1, const PPPLib::tSatObsUnit &p2) {
        return p1.sat.sat_.no<p2.sat.sat_.no;
    }

    bool cReadGnssObs::SortEpochSatData(PPPLib::tEpochSatUnit& epoch_sat_data) {
        if(epoch_sat_data.sat_num<=0) return false;
        sort(epoch_sat_data.epoch_data.begin(),epoch_sat_data.epoch_data.end(),CmpEpochSatData);
    }

    int cReadGnssObs::DecodeEpoch(cTime& t, int& obs_flag) {
        int n=0;
        Str2Int(line_str_.substr(32,3),n);
        if(n<=0) return 0;

        Str2Int(line_str_.substr(31,1),obs_flag);
        if(3<obs_flag&&obs_flag<5) return n;

        if(line_str_[0]!='>'||t.Str2Time(line_str_.substr(1,28))!=0)
            return 0;

        return n;
    }

    bool cReadGnssObs::DecodeEpochSatObs(tSatObsUnit& obs) {
        cGnssSignal *gnss_signal={0};
        double val[MAX_GNSS_OBS_TYPE]={0};
        unsigned char lli[MAX_GNSS_OBS_TYPE]={0};
        string sat_id;
        int i,j,num=0;
        bool stat=true;

        if(rnx_ver_>2.99){
            sat_id=line_str_.substr(0,3);
            if(!sat_id.compare(0,1,"S")) return false;
            if(!sat_id.compare(0,1,"I")) return false;
            if(!sat_id.compare(1,1," ")) sat_id[1]='0';
            obs.sat=cSat(sat_id);
            obs.sat.SatId2No();
            if(!(obs.sat.sat_.sys&sys_mask_)) return false;
        }
        if(!obs.sat.sat_.no) stat=false;

        switch(obs.sat.sat_.sys){
            case SYS_BDS: gnss_signal=signal_index+SYS_INDEX_BDS;break;
            case SYS_GAL: gnss_signal=signal_index+SYS_INDEX_GAL;break;
            case SYS_GLO: gnss_signal=signal_index+SYS_INDEX_GLO;break;
            case SYS_QZS: gnss_signal=signal_index+SYS_INDEX_QZS;break;
            case SYS_IRN: gnss_signal=signal_index+SYS_INDEX_IRN;break;
            default:
                gnss_signal=signal_index;break;
        }

        if(!stat) return false;


        for(i=0;i<MAX_GNSS_FRQ_NUM+GNSS_NUM_EXOBS;i++){
            obs.P[i]=obs.L[i]=0.0;obs.D[i]=0.0f;
            obs.SNR[i]=obs.LLI[i]=obs.code[i]=0;
        }

        for(i=0,j=rnx_ver_<=2.99?0:3;i<gnss_signal->GetGnssSignal()->n&&j+15<line_str_.length();i++,j+=16){
            Str2Double(line_str_.substr(j,14),val[i]);
            val[i]+=gnss_signal->GetGnssSignal()->shift[i];

            Str2Int(line_str_.substr(j+14,1),num);
            lli[i]=(unsigned char)num&3;
        }

        int n,m,p[MAX_GNSS_OBS_TYPE],k[16],l[16];
        for(i=n=m=0;i<gnss_signal->GetGnssSignal()->n;i++){
            p[i]=rnx_ver_<=2.11?gnss_signal->GetGnssSignal()->frq[i]:gnss_signal->GetGnssSignal()->pos[i];
            if(gnss_signal->GetGnssSignal()->type[i]==0&&p[i]==0) k[n++]=i;
            if(gnss_signal->GetGnssSignal()->type[i]==0&&p[i]==1) l[m++]=i;
        }

        if(n>=2){
            if(val[k[0]]==0.0&&val[k[1]]==0.0){
                p[k[0]]=-1;p[k[1]]=-1;
            }
            else if(val[k[0]]!=0.0&&val[k[1]]==0.0){
                p[k[0]]=0;p[k[1]]=-1;
            }
            else if(val[k[0]]==0.0&&val[k[1]]!=0.0){
                p[k[0]]=-1;p[k[1]]=0;
            }
            else if(gnss_signal->GetGnssSignal()->pri[k[1]]>gnss_signal->GetGnssSignal()->pri[k[0]]){
                p[k[1]]=0;p[k[0]]=-1;
            }
            else{
                p[k[0]]=0;p[k[1]]=-1;
            }

        }
        if(m>=2){
            if(val[l[0]]==0.0&&val[l[1]]==0.0){
                p[l[0]]=-1;p[l[1]]=-1;
            }
            else if(val[l[0]]!=0.0&&val[l[1]]==0.0){
                p[l[0]]=1;p[l[1]]=-1;
            }
            else if(val[l[0]]==0.0&&val[l[1]]!=0.0){
                p[l[0]]=-1;p[l[1]]=1;
            }
            else if(gnss_signal->GetGnssSignal()->pri[l[1]]>gnss_signal->GetGnssSignal()->pri[l[0]]){
                p[l[1]]=1;p[l[0]]=-1;
            }
            else{
                p[l[0]]=1;p[l[1]]=-1;
            }
        }

        for(i=0;i<gnss_signal->GetGnssSignal()->n;i++){
            if(p[i]<0||val[i]==0.0) continue;
            switch(gnss_signal->GetGnssSignal()->type[i]){
                case GNSS_OBS_CODE:    obs.P[p[i]]=val[i];obs.code[p[i]]=gnss_signal->GetGnssSignal()->code[i];break;
                case GNSS_OBS_PHASE:   obs.L[p[i]]=val[i];obs.LLI[p[i]]=lli[i];break;
                case GNSS_OBS_DOPPLER: obs.D[p[i]]=(float)val[i];break;
                case GNSS_OBS_SNR:     obs.SNR[p[i]]=(unsigned char)(val[i]*4.0+0.5); break;
            }
        }
        return true;
    }

    int cReadGnssObs::ReadObsBody(tEpochSatUnit& epoch_sat_data) {
        int line_idx=0,num_sat=0,n,obs_flag=0;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            tSatObsUnit sat_data={0};
            if(line_idx==0){
                if((num_sat=DecodeEpoch(epoch_sat_data.obs_time,obs_flag))<=0) continue;
                if(GetGnssData()->GetStartTime()->t_.long_time!=0.0){
                    if(GetGnssData()->GetStartTime()->TimeDiff(epoch_sat_data.obs_time.t_)>0) continue;
                }
                if(GetGnssData()->GetEndTime()->t_.long_time!=0.0){
                    if(GetGnssData()->GetEndTime()->TimeDiff(epoch_sat_data.obs_time.t_)<0) continue;
                }
            }
            else if(obs_flag<=2||obs_flag==6){
                if(DecodeEpochSatObs(sat_data)){
                    epoch_sat_data.epoch_data.push_back(sat_data);
                }
            }
            if(++line_idx>num_sat){
                return num_sat;
            }
        }
        return 0;
    }

    void cReadGnssObs::SetGnssTimeSpan(cTime* ts,cTime* te) {
        if(ts->TimeDiff(te->t_)>=0){
            LOG(WARNING)<<"Wrong time interval setting";
            return;
        }
        else{
            gnss_data_->SetTimeSpan(ts,te);
        }
    }


    cGnssObs* cReadGnssObs::GetGnssData() {
        return gnss_data_;
    }

    tNav * cReadGnssObs::GetGnssNav() {
        return nav_;
    }

    bool cReadGnssObs::ReadHead() {
        int i,j,k,num_signal,idx_signal,num_line=0,prn,fcn;
        tStaInfoUnit *sta=&nav_->sta_paras[GetGnssData()->rcv_idx_];

        if(!ReadRnxHead()) return false;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if (line_str_.find("MARKER NAME")!=string::npos && sta)
                sta->name=StrTrim(line_str_.substr(0,10));
            else if (line_str_.find("MARKER NUMBER")!=string::npos && sta)
                sta->marker=StrTrim(line_str_.substr(0,20));
            else if (line_str_.find("MARKER TYPE")!=string::npos) continue;
            else if (line_str_.find("OBSERVER / AGENCY")!=string::npos) continue;
            else if (line_str_.find("REC # / TYPE / VERS")!=string::npos && sta){
                sta->ant_seri=StrTrim(line_str_.substr(0,20));
                sta->rec_type=StrTrim(line_str_.substr(20,20));
                sta->firm_ver=StrTrim(line_str_.substr(40,20));
            }
            else if (line_str_.find("ANT # / TYPE")!=string::npos && sta){
                sta->ant_seri=StrTrim(line_str_.substr(0,20));
                sta->ant_desc=StrTrim(line_str_.substr(20,20));
            }
            else if (line_str_.find("APPROX POSITION XYZ")!=string::npos && sta)
                for (i=0;i<3;i++){
                    Str2Double(line_str_.substr(i*14,14),sta->apr_pos[i]);
                }
            else if (line_str_.find("ANTENNA: DELTA H/E/N")!=string::npos && sta){
                Str2Double(line_str_.substr(0,14),sta->del[2]);  /* h */
                Str2Double(line_str_.substr(14,14),sta->del[0]); /* e */
                Str2Double(line_str_.substr(28,14),sta->del[1]); /* n */
            }
            else if (line_str_.find("ANTENNA: DELTA X/Y/Z")!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("ANTENNA: PHASECENTER")!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("ANTENNA: B.SIGHT XYZ")!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("ANTENNA: ZERODIR AZI")!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("ANTENNA: ZERODIR XYZ")!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("CENTER OF MASS: XYZ" )!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("SYS / # / OBS TYPES" )!=string::npos) { /* ver.3 */
                if((i=kSatSysCode.find(line_str_[0]))==string::npos){
                    LOG(WARNING)<<"Invalid satellite system: "<<line_str_[0];
                    continue;
                }
                Str2Int(line_str_.substr(3,3),num_signal);
                for(j=idx_signal=0,k=7;j<num_signal;j++,k+=4){
                    if(k>58){
                        if(!getline(inf_,line_str_)) break;
                        k=7;
                    }
                    if(idx_signal<MAX_GNSS_OBS_TYPE-1) obs_type_code_[i][idx_signal++]=line_str_.substr(k,3);
                }
                if(i==1){
                    if(rnx_ver_<3.04){
                        for(j=0;j<num_signal;j++){
                            if(obs_type_code_[i][j][1]=='1'){
                                obs_type_code_[i][j][1]='2';
                                LOG(DEBUG)<<"BD2 change C1x to C2x";
                            }
                        }
                    }
                }
            }
            else if (line_str_.find("WAVELENGTH FACT L1/2")!=string::npos) continue; /* opt ver.2 */
            else if (line_str_.find("# / TYPES OF OBSERV")!=string::npos) { /* ver.2 */
            }
            else if (line_str_.find("SIGNAL STRENGTH UNIT")!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("INTERVAL"			 )!=string::npos) continue; /* opt */
            else if (line_str_.find("TIME OF FIRST OBS"   )!=string::npos) {
                if      (!line_str_.compare(48,3,"GPS")) time_sys_=TIME_GPS;
                else if (!line_str_.compare(48,3,"BDT")) time_sys_=TIME_BDS; /* ver.3.02 */
                else if (!line_str_.compare(48,3,"GAL")) time_sys_=TIME_GAL;
                else if (!line_str_.compare(48,3,"GLO")) time_sys_=TIME_utc;
                else if (!line_str_.compare(48,3,"QZS")) time_sys_=TIME_QZS; /* ver.3.02 */
            }
            else if (line_str_.find("TIME OF LAST OBS"    )!=string::npos) continue; /* opt */
            else if (line_str_.find("RCV CLOCK OFFS APPL" )!=string::npos) continue; /* opt */
            else if (line_str_.find("SYS / DCBS APPLIED"  )!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("SYS / PCVS APPLIED"  )!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("SYS / SCALE FACTOR"  )!=string::npos) continue; /* opt ver.3 */
            else if (line_str_.find("SYS / PHASE SHIFTS"  )!=string::npos) continue; /* ver.3.01 */
            else if (line_str_.find("GLONASS SLOT / FRQ #")!=string::npos && &nav_) { /* ver.3.02 */
                for (i=0;i<8;i++) {
                    string aa=line_str_.substr(7*i+8,2);
                    if (line_str_.compare(7*i+4,1,"R")!=0||!line_str_.compare(7*i+8,2,"  ")) continue;
                    string a=line_str_.substr(7*i+5,2);
                    string b=line_str_.substr(7*i+8,2);
                    Str2Int(line_str_.substr(7*i+5,2),prn);
                    Str2Int(line_str_.substr(7*i+8,2),fcn);
                    if (1<=prn&&prn<=GLO_MAX_PRN) nav_->glo_frq_num[prn-1]=fcn+8;
                }
            }
            else if (line_str_.find("GLONASS COD/PHS/BIS" )!=string::npos && &nav_) {  /* ver.3.02 */
                for (i=0;i<4;i++) {
                    if      (line_str_.compare(13*i+1,3,"C1C"))
                        Str2Double(line_str_.substr(13*i+5,8), nav_->glo_cp_bias[0]);
                    else if (line_str_.compare(13*i+1,3,"C1P"))
                        Str2Double(line_str_.substr(13*i+5,8),nav_->glo_cp_bias[1]);
                    else if (line_str_.compare(13*i+1,3,"C2C"))
                        Str2Double(line_str_.substr(13*i+5,8),nav_->glo_cp_bias[2]);
                    else if (line_str_.compare(13*i+1,3,"C2P"))
                        Str2Double(line_str_.substr(13*i+5,8),nav_->glo_cp_bias[3]);
                }
            }
            else if (line_str_.find("LEAP SECONDS")!=string::npos && &nav_) {/* opt */
                Str2Int(line_str_.substr(0,6),nav_->leaps);
            }
            else if (line_str_.find("# OF SALTELLITES")!=string::npos) continue;/* opt */
            else if (line_str_.find("PRN / # OF OBS"  )!=string::npos) continue;/* opt */
            else if (line_str_.find("PGM / RUN BY / DATE")!=string::npos) continue;
            else if (line_str_.find("COMMENT" )!=string::npos) continue;
            if (line_str_.find("END OF HEADER")!=string::npos)
                break;
            if (++num_line>=1024 && rnx_type_.compare(" ")==0) return false; /* no rinex file */
        }

        signal_index[SYS_INDEX_GPS].GnssSignalIndex(SYS_GPS,obs_type_code_[SYS_INDEX_GPS]);
        signal_index[SYS_INDEX_BDS].GnssSignalIndex(SYS_BDS,obs_type_code_[SYS_INDEX_BDS]);
        signal_index[SYS_INDEX_GAL].GnssSignalIndex(SYS_GAL,obs_type_code_[SYS_INDEX_GAL]);
        signal_index[SYS_INDEX_GLO].GnssSignalIndex(SYS_GLO,obs_type_code_[SYS_INDEX_GLO]);
        signal_index[SYS_INDEX_QZS].GnssSignalIndex(SYS_QZS,obs_type_code_[SYS_INDEX_QZS]);

        return true;

    }

    bool cReadGnssObs::Reading() {
        unsigned char slips[MAX_SAT_NUM][MAX_GNSS_FRQ_NUM]={{0}};

        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss obs file error: "<<file_;
            return false;
        }

        if(!ReadHead()) return false;

        int epoch_sat_num=0;
        tEpochSatUnit epoch_sat_data={0};
        int kk=0;
        while((epoch_sat_num=ReadObsBody(epoch_sat_data))>=0){
            if(inf_.eof()) break;
            kk++;
            for(int i=0;i<epoch_sat_data.sat_num;i++){
                if(time_sys_==TIME_utc) epoch_sat_data.obs_time.Utc2Gpst();

                for(int f=0;f<MAX_GNSS_FRQ_NUM;f++){
                    if(epoch_sat_data.epoch_data[i].LLI[f]&1) slips[epoch_sat_data.epoch_data[i].sat.sat_.no-1][f]|=0x01;
                }
            }

            for(int i=0;i<epoch_sat_data.sat_num;i++){
                for(int f=0;f<MAX_GNSS_FRQ_NUM;f++){
                    if(slips[epoch_sat_data.epoch_data[i].sat.sat_.no-1][f]&1) epoch_sat_data.epoch_data[i].LLI[f]|=0x01;
                    slips[epoch_sat_data.epoch_data[i].sat.sat_.no-1][f]=0;
                }
            }

            epoch_sat_data.sat_num=epoch_sat_data.epoch_data.size();
            SortEpochSatData(epoch_sat_data);
            gnss_data_->GetGnssObs().push_back(epoch_sat_data);

            epoch_sat_data.sat_num=0;
            epoch_sat_data.epoch_data.clear();
        }

        for(int i=0;i<NSYS;i++){
            gnss_data_->signal_[i]= reinterpret_cast<tGnssSignal *>(&signal_index[i]);
        }

        gnss_data_->epoch_num=gnss_data_->GetGnssObs().size();
        if(gnss_data_->GetStartTime()->t_.long_time==0.0){
            gnss_data_->GetStartTime()->t_=gnss_data_->GetGnssObs()[0].obs_time.t_;
        }
        if(gnss_data_->GetEndTime()->t_.long_time==0.0){
            gnss_data_->GetEndTime()->t_=gnss_data_->GetGnssObs().back().obs_time.t_;
        }

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssBrdEph::cReadGnssBrdEph() {}

    cReadGnssBrdEph::cReadGnssBrdEph(string file_path, PPPLib::tNav &nav) {
        file_=file_path;
        nav_=&nav;
    }

    cReadGnssBrdEph::~cReadGnssBrdEph() {}

    bool cReadGnssBrdEph::CmpBrdEph(const tBrdEphUnit& p1, const tBrdEphUnit& p2){
        return p1.ttr.t_.long_time!=p2.ttr.t_.long_time?(int)(p1.ttr.t_.long_time<p2.ttr.t_.long_time):
               (p1.toe.t_.long_time!=p2.toe.t_.long_time?(int)(p1.toe.t_.long_time<p2.toe.t_.long_time):
               (p1.sat.sat_.no!=p2.sat.sat_.no?(int)(p1.sat.sat_.no<p2.sat.sat_.no):(int)(p1.code<p2.code)));
    }

    bool cReadGnssBrdEph::CmpBrdGloEph(const tBrdGloEphUnit &p1, const tBrdGloEphUnit &p2) {
        return p1.tof.t_.long_time!=p2.tof.t_.long_time?(int)(p1.tof.t_.long_time<p2.tof.t_.long_time):
               (p1.toe.t_.long_time!=p2.toe.t_.long_time?(int)(p1.toe.t_.long_time<p2.toe.t_.long_time):
               (p1.sat.sat_.no!=p2.sat.sat_.no?(p1.sat.sat_.no<p2.sat.sat_.no):false));

    }

    bool cReadGnssBrdEph::SortBrdEph() {
        vector<tBrdEphUnit>& ephs=nav_->brd_eph;
        if(nav_->brd_eph.size()<=0) return false;

        sort(ephs.begin(),ephs.end(),CmpBrdEph);

        int i,j;
        for(i=1,j=0;i<ephs.size();i++){
            if((ephs[i].sat.sat_.no==ephs[j].sat.sat_.no)&&(ephs[i].iode==ephs[j].iode)){
                __gnu_cxx::__normal_iterator<tBrdEphUnit *, vector<tBrdEphUnit>> iter= ephs.begin()+j;
                ephs.erase(iter);
            }
        }
        return true;
    }

    bool cReadGnssBrdEph::SortBrdGloEph() {
        if(nav_->brd_glo_eph.size()<=0) return false;

        sort(nav_->brd_glo_eph.begin(),nav_->brd_glo_eph.end(),CmpBrdGloEph);

        int i,j;
        for(i=1,j=0;i<nav_->brd_glo_eph.size();i++){
            if((nav_->brd_glo_eph[i].sat.sat_.no==nav_->brd_glo_eph[j].sat.sat_.no)&&(nav_->brd_glo_eph[i].iode==nav_->brd_glo_eph[j].iode)){
                __gnu_cxx::__normal_iterator<tBrdGloEphUnit *, vector<tBrdGloEphUnit>> iter= nav_->brd_glo_eph.begin()+j;
                nav_->brd_glo_eph.erase(iter);
            }
        }
        return true;
    }

    void cReadGnssBrdEph::ClearEphData() {
        for(auto &i:eph_data_) i=0.0;
    }

    void cReadGnssBrdEph::DecodeEph(PPPLib::cTime toc, PPPLib::cSat sat, PPPLib::tBrdEphUnit &brd_eph) {
        brd_eph.sat=sat;
        brd_eph.toc=toc;

        brd_eph.f0=eph_data_[0];
        brd_eph.f1=eph_data_[1];
        brd_eph.f2=eph_data_[2];

        brd_eph.A=SQR(eph_data_[10]); brd_eph.e=eph_data_[ 8]; brd_eph.i0  =eph_data_[15]; brd_eph.Omg0=eph_data_[13];
        brd_eph.omg =eph_data_[17]; brd_eph.M0 =eph_data_[ 6]; brd_eph.deln=eph_data_[ 5]; brd_eph.Omgd=eph_data_[18];
        brd_eph.idot=eph_data_[19]; brd_eph.crc=eph_data_[16]; brd_eph.crs =eph_data_[ 4]; brd_eph.cuc =eph_data_[ 7];
        brd_eph.cus =eph_data_[ 9]; brd_eph.cic=eph_data_[12]; brd_eph.cis =eph_data_[14];

        if(sat.sat_.sys==SYS_GPS||sat.sat_.sys==SYS_QZS){
            brd_eph.iode=(int)eph_data_[3];
            brd_eph.iodc=(int)eph_data_[26];
            brd_eph.toes=     eph_data_[11];
            brd_eph.week=     eph_data_[21];
            brd_eph.toe.Gpst2Time(brd_eph.week,eph_data_[11],SYS_GPS)->AdjWeek(toc);
            brd_eph.ttr.Gpst2Time(brd_eph.week,eph_data_[27],SYS_GPS)->AdjWeek(toc);

            brd_eph.code=(int)eph_data_[20];
            brd_eph.svh =(int)eph_data_[24];
            brd_eph.sva =0;
            brd_eph.tgd[0]=eph_data_[25]; //L1_L2
        }
        else if(sat.sat_.sys==SYS_GAL){
            brd_eph.iode=(int)eph_data_[3];
            brd_eph.toes=     eph_data_[11];
            brd_eph.week=(int)eph_data_[21];
            brd_eph.toe.Gpst2Time(brd_eph.week,eph_data_[11],SYS_GPS)->AdjWeek(toc);
            brd_eph.ttr.Gpst2Time(brd_eph.week,eph_data_[27],SYS_GPS)->AdjWeek(toc);
            brd_eph.code=(int)eph_data_[20];
            brd_eph.svh =(int)eph_data_[24];
            brd_eph.sva=0;
            brd_eph.tgd[0]=eph_data_[25]; //E1_E5a
            brd_eph.tgd[1]=eph_data_[26]; //E1_E5b
        }
        else if(sat.sat_.sys==SYS_BDS){
            cTime toc_temp=toc;
            brd_eph.toc=toc_temp.Bdst2Gpst();
            brd_eph.iode=(int)eph_data_[3];
            brd_eph.iodc=(int)eph_data_[28];
            brd_eph.toes=     eph_data_[11];
            brd_eph.week=(int)eph_data_[21];
            brd_eph.toe.Gpst2Time(brd_eph.week,eph_data_[11],SYS_BDS)->Bdst2Gpst();
            brd_eph.ttr.Gpst2Time(brd_eph.week,eph_data_[27],SYS_BDS)->Bdst2Gpst();
            brd_eph.toe.AdjWeek(toc);
            brd_eph.ttr.AdjWeek(toc);
            brd_eph.svh=(int)eph_data_[24];
            brd_eph.sva=0;
            brd_eph.tgd[0]=eph_data_[25]; //B1I_B3I
            brd_eph.tgd[1]=eph_data_[26]; //B2I_B3I
        }
        else if(sat.sat_.sys==SYS_IRN){
            brd_eph.iode=(int)eph_data_[3];
            brd_eph.toes=     eph_data_[11];
            brd_eph.week=(int)eph_data_[21];
            brd_eph.toe.Gpst2Time(brd_eph.week,eph_data_[11],SYS_GPS)->AdjWeek(toc);
            brd_eph.ttr.Gpst2Time(brd_eph.week,eph_data_[27],SYS_GPS)->AdjWeek(toc);
            brd_eph.svh=(int)eph_data_[24];
            brd_eph.sva=0;
            brd_eph.tgd[0]=eph_data_[25];
        }

        if(brd_eph.iode<0||1023<brd_eph.iode) brd_eph.svh=-1;
        if(brd_eph.iodc<0||1023<brd_eph.iodc) brd_eph.svh=-1;
    }

    void cReadGnssBrdEph::DecodeGloEph(cTime toc,const cSat& sat,PPPLib::tBrdGloEphUnit &glo_eph) {
        double tow,tod;
        int week,dow;
        cTime tof;

        glo_eph.sat=sat;
        tow=toc.Time2Gpst(&week,&dow,SYS_GPS);
        toc.Gpst2Time(week,floor((tow+450.0)/900.0)*900,SYS_GPS);
        dow=(int)floor(tow/86400.0);
        tod=rnx_ver_<=2.99?eph_data_[2]:fmod(eph_data_[2],86400.0);
        tof.Gpst2Time(week,tod+dow*86400.0,SYS_GPS);
        tof.AdjDay(toc);
        cTime t1=toc,t2=toc;
        glo_eph.toe=t1.Utc2Gpst();
        glo_eph.tof=t2.Utc2Gpst();

        glo_eph.iode=(int)(fmod(tow+10800.0,86400.0)/900.0+0.5);
        glo_eph.taun=-eph_data_[0];
        glo_eph.gamn= eph_data_[1];
        glo_eph.pos[0]=eph_data_[3]*1E3; glo_eph.pos[1]=eph_data_[7]*1E3; glo_eph.pos[2]=eph_data_[11]*1E3;
        glo_eph.vel[0]=eph_data_[4]*1E3; glo_eph.vel[1]=eph_data_[8]*1E3; glo_eph.vel[2]=eph_data_[12]*1E3;
        glo_eph.acc[0]=eph_data_[5]*1E3; glo_eph.acc[1]=eph_data_[9]*1E3; glo_eph.acc[2]=eph_data_[13]*1E3;
        glo_eph.svh=(int)eph_data_[ 6];
        glo_eph.frq=(int)eph_data_[10];
        glo_eph.age=(int)eph_data_[14];
        if (glo_eph.frq>128) glo_eph.frq-=256;
        if(glo_eph.frq<GLO_MIN_FREQ||GLO_MAX_FREQ<glo_eph.frq){
            glo_eph.svh=-1;
        }
    }

    bool cReadGnssBrdEph::ReadBrdBody() {
        int i=0,j,sp=3,prn=0,flag=0,sys=SYS_GPS;
        string id;
        cSat sat;
        cTime toc;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(line_str_.compare(0,3,"   ")!=0) i=0;
            if(i==0){
                if(rnx_ver_>=3.0||sat_sys_==SYS_GAL||sat_sys_==SYS_QZS||sat_sys_==SYS_NONE){
                    id=line_str_.substr(0,3);
                    if(id=="R01"){
                        int a=1;
                    }
                    sat=cSat(id);
                    sat.SatId2No();
                    sp=4;
                    if(rnx_ver_>=3.0) sys=sat.sat_.sys;
                }
                else{
                    Str2Int(line_str_.substr(0,2),prn);
                    if(sat_sys_==SYS_GLO) sat=cSat(SYS_GLO,prn);
                    else sat=cSat(SYS_GPS,prn);
                }
                if(toc.Str2Time(line_str_.substr(sp,19))) flag=0;
                else flag=1;

                for(j=0;j<3;j++){
                    Str2Double(line_str_.substr(sp+19*(j+1),19), eph_data_[i++]);
                }
            }
            else if(flag==1){
                for(j=0;j<4;j++){
                    if(sp+19*(j+1)<=line_str_.size()){
                        Str2Double(line_str_.substr(sp+19*j,19),eph_data_[i++]);
                    }
                    else eph_data_[i++]=0.0;
                }

                if(sys==SYS_GLO&&i>=15){
                    if(!(sys_mask_&sys)) continue;
                    tBrdGloEphUnit glo_eph={0};
                    DecodeGloEph(toc,sat,glo_eph);
                    nav_->brd_glo_eph.push_back(glo_eph);
                    ClearEphData();
                }
                else if(i>=31){
                    if(!(sys_mask_&sys)) continue;
                    tBrdEphUnit eph={0};
                    DecodeEph(toc,sat,eph);
                    nav_->brd_eph.push_back(eph);
                    ClearEphData();
                }
            }
        }
        return true;
    }

    tNav* cReadGnssBrdEph::GetGnssNav() {return nav_;}

    bool cReadGnssBrdEph::ReadHead() {
        int nline=0,i,j;

        if(!ReadRnxHead()) return false;

        while(getline(inf_,line_str_)&&!inf_.eof()){

            if (line_str_.find("IONOSPHERIC CORR",60)!=string::npos) { /* opt ver.3 */
                if (line_str_.compare(0,4,"GPSA")==0) {
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_GPS][i]);
                }
                else if (line_str_.compare(0,4,"GPSB")==0) {
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_GPS][i+4]);
                }
                else if (line_str_.compare(0,3,"GAL")==0) {
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_GAL][i]);
                }
                else if (line_str_.compare(0,4,"QZSA")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_QZS][i]);
                }
                else if (line_str_.compare(0,4,"QZSB")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_QZS][i+4]);
                }
                else if (line_str_.compare(0,4,"BDSA")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_BDS][i]);
                }
                else if (line_str_.compare(0,4,"BDSB")==0) { /* v.3.02 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_BDS][i+4]);
                }
                else if (line_str_.compare(0,4,"IRNA")==0) { /* v.3.03 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_IRN][i]);
                }
                else if (line_str_.compare(0,4,"IRNB")==0) { /* v.3.03 */
                    for (i=0,j=5;i<4;i++,j+=12) Str2Double(line_str_.substr(j,12),nav_->ion_para[SYS_INDEX_IRN][i+4]);
                }
            }
            else if (line_str_.find("TIME SYSTEM CORR",60)!=string::npos) { /* opt ver.3 */
                if (line_str_.compare(0,4,"GPUT")==0) {
                    Str2Double(line_str_.substr( 5,17),nav_->utc_para[SYS_INDEX_GPS][0]);
                    Str2Double(line_str_.substr(22,16),nav_->utc_para[SYS_INDEX_GPS][1]);
                    Str2Double(line_str_.substr(38, 7),nav_->utc_para[SYS_INDEX_GPS][2]);
                    Str2Double(line_str_.substr(45, 5),nav_->utc_para[SYS_INDEX_GPS][3]);
                }
                else if (line_str_.compare(0,4,"GLUT")==0) {
                    Str2Double(line_str_.substr( 5,17),nav_->utc_para[SYS_INDEX_GLO][0]);
                    Str2Double(line_str_.substr(22,16),nav_->utc_para[SYS_INDEX_GLO][1]);
                    Str2Double(line_str_.substr(38, 7),nav_->utc_para[SYS_INDEX_GLO][2]);
                    Str2Double(line_str_.substr(45, 5),nav_->utc_para[SYS_INDEX_GLO][3]);
                }
                else if (line_str_.compare(0,4,"GAUT")==0) { /* v.3.02 */
                    Str2Double(line_str_.substr( 5,17),nav_->utc_para[SYS_INDEX_GAL][0]);
                    Str2Double(line_str_.substr(22,16),nav_->utc_para[SYS_INDEX_GAL][1]);
                    Str2Double(line_str_.substr(38, 7),nav_->utc_para[SYS_INDEX_GAL][2]);
                    Str2Double(line_str_.substr(45, 5),nav_->utc_para[SYS_INDEX_GAL][3]);
                }
                else if (line_str_.compare(0,4,"QZUT")==0) { /* v.3.02 */
                    Str2Double(line_str_.substr( 5,17),nav_->utc_para[SYS_INDEX_GAL][0]);
                    Str2Double(line_str_.substr(22,16),nav_->utc_para[SYS_INDEX_GAL][1]);
                    Str2Double(line_str_.substr(38, 7),nav_->utc_para[SYS_INDEX_GAL][2]);
                    Str2Double(line_str_.substr(45, 5),nav_->utc_para[SYS_INDEX_GAL][3]);
                }
                else if (line_str_.compare(0,4,"BDUT")==0) { /* v.3.02 */
                    Str2Double(line_str_.substr( 5,17),nav_->utc_para[SYS_INDEX_BDS][2]);
                    Str2Double(line_str_.substr(22,16),nav_->utc_para[SYS_INDEX_BDS][2]);
                    Str2Double(line_str_.substr(38, 7),nav_->utc_para[SYS_INDEX_BDS][2]);
                    Str2Double(line_str_.substr(45, 5),nav_->utc_para[SYS_INDEX_BDS][3]);
                }
                else if (line_str_.compare(0,4,"IRUT")==0) { /* v.3.03 */
                    Str2Double(line_str_.substr( 5,17),nav_->utc_para[SYS_INDEX_IRN][0]);
                    Str2Double(line_str_.substr(22,16),nav_->utc_para[SYS_INDEX_IRN][1]);
                    Str2Double(line_str_.substr(38, 7),nav_->utc_para[SYS_INDEX_IRN][2]);
                    Str2Double(line_str_.substr(45, 5),nav_->utc_para[SYS_INDEX_IRN][3]);
                }
            }
            else if (line_str_.find("LEAP SECONDS",60)!=string::npos) { /* opt */
                Str2Int(line_str_.substr(0,6),nav_->leaps);
            }
            if (line_str_.find("END OF HEADER")!=string::npos) return true;
            if (++nline>=1024 && rnx_type_.compare(" ")==0) break; /* no rinex file */
        }
        return false;
    }

    bool cReadGnssBrdEph::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss broadcast ephemeris file error: "<<file_;
            return false;
        }

        if(!ReadHead()) return false;

        if(!ReadBrdBody()) return false;

        SortBrdEph();
        SortBrdGloEph();

        int i,j;
        for(i=1;i<=NUM_GLO_SAT;i++){
            nav_->glo_frq_num[i-1]=0;
        }
        if(nav_->brd_glo_eph.size()){
            for(i=1;i<=NUM_GLO_SAT;i++){
                for(j=0;j<nav_->brd_glo_eph.size();j++){
                    if(i==nav_->brd_glo_eph[j].sat.sat_.prn) break;
                }
                if(j<nav_->brd_glo_eph.size()) nav_->glo_frq_num[i-1]=nav_->brd_glo_eph[j].frq;
            }
        }

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssPreEph::cReadGnssPreEph(){}

    cReadGnssPreEph::cReadGnssPreEph(string file_path, PPPLib::tNav &nav) {
        file_=file_path;
        nav_=&nav;
        sys_mask_=SYS_ALL;
    }

    cReadGnssPreEph::~cReadGnssPreEph() {}

    void cReadGnssPreEph::ReadPreOrbHead() {
        int clm=0;
        string last_block;
        num_sat_=0;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(!line_str_.substr(0,1).compare("*")) break;
            if(clm==0){
                pre_eph_time_.Str2Time(line_str_.substr(3,28));
            }
            else if(!line_str_.substr(0,2).compare("+ ")){
                if(!last_block.compare("##")||clm==2) Str2Int(line_str_.substr(3,3),num_sat_);
            }
            else if(!line_str_.substr(0,2).compare("%c")&&!last_block.compare("++")){

            }
            else if(!line_str_.substr(0,2).compare("%f")&&!last_block.compare("%c")){

            }

            last_block=line_str_.substr(0,2);
            clm++;
        }
    }

    void cReadGnssPreEph::ReadPreOrbBody() {
        cTime epoch_t;
        cSat sat;
        tPreOrbUnit eph={nullptr};

        while(!inf_.eof()){
            if(!line_str_.substr(0,3).compare("EOF")) break;
            if(line_str_[0]!='*'||epoch_t.Str2Time(line_str_.substr(3,28))) continue;
            nav_->pre_eph.push_back(eph);
            nav_->pre_eph.back().t_tag=epoch_t;
            for(int i=0;i<num_sat_&&getline(inf_,line_str_);i++){
                if(line_str_.length()<4||line_str_[0]!='P') continue;
                sat=cSat(line_str_.substr(1,3));
                sat.SatId2No();
                if(!(sat.sat_.no)) continue;
                if(!(sat.sat_.sys&sys_mask_)) continue;
                for(int j=0;j<4;j++){
                    double std=0.0,val=0.0;
                    Str2Double(line_str_.substr(4+j*14,14),val);
                    if(line_str_.length()>=80) Str2Double(line_str_.substr(61+j*3,j<3?2:3),std);
                    if(line_str_[0]=='P'){
                        if(val!=0.0&&fabs(val-999999.999999)>=1E-6){
                            nav_->pre_eph.back().pos[sat.sat_.no-1][j]=val*(j<3?1000.0:1E-6);
                        }
                    }
                }
            }
            getline(inf_,line_str_);
        }
    }

    void cReadGnssPreEph::ReadPreClkHead() {
        int block=0;
        double bias;
        cSat sat;
        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(strstr(line_str_.c_str(),"COMMENT")){
                if (strstr(line_str_.c_str(),"WIDELANE SATELLITE FRACTIONNAL BIASES FOR GALILEO")||
                    strstr(line_str_.c_str(),"WIDELANE SATELLITE FRACTIONNAL BIASES USED IN THIS SOLUTION")) {
                    block=1;
                }
                else if(block){
                    if(!strncmp(line_str_.c_str(),"WL",2)){
                        sat=cSat(line_str_.substr(3,3));
                        sat.SatId2No();
                        if(sat.sat_.no>0&&sscanf(line_str_.c_str()+40,"%lf",&bias)==1){
                            nav_->wide_line_bias[sat.sat_.no-1]=bias;
                        }
                    }
                }
            }
            else if(!line_str_.substr(60,13).compare("END OF HEADER")) break;
            else continue;
        }
    }

    void cReadGnssPreEph::SetGnssSatMask(int mask) {sys_mask_=mask;}

    void cReadGnssPreEph::ReadHead(int type) {
        if(type==0) ReadPreOrbHead();
        else if(type==1) ReadPreClkHead();
    }

    void cReadGnssPreEph::ReadPreClkBody() {
        int k=1;
        double data[2]={0};
        cTime epoch_t;
        cSat sat;
        tPreClkUnit clk={nullptr};

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(line_str_.substr(0,2).compare("AS")||epoch_t.Str2Time(line_str_.substr(8,26))) continue;
            if(nav_->pre_clk.empty()||fabs(epoch_t.TimeDiff(nav_->pre_clk.back().t_tag.t_)>1E-9)){
                nav_->pre_clk.push_back(clk);
            }
            nav_->pre_clk.back().t_tag=epoch_t;
            sat=cSat(line_str_.substr(3,3));
            sat.SatId2No();
            if(!(sat.sat_.sys&sys_mask_)) continue;
            if(!(sat.sat_.no)) continue;
            if(line_str_.length()>60) k=2;
            for(int i=0,j=40;i<k;i++,j+=20) Str2Double(line_str_.substr(j,19),data[i]);
            nav_->pre_clk.back().clk[sat.sat_.no-1]=data[0];
            nav_->pre_clk.back().std[sat.sat_.no-1]=(float)data[1];
        }
    }

    bool cReadGnssPreEph::Reading(int type) {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss precise products file error: "<<file_;
            return false;
        }

        ReadHead(type);
        if(type==0){
            ReadPreOrbBody();
        }
        else if(type==1){
            ReadPreClkBody();
        }

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssCodeBias::cReadGnssCodeBias() {}

    cReadGnssCodeBias::cReadGnssCodeBias(string file_path, PPPLib::tNav &nav) {
        file_=file_path;
        nav_=&nav;
    }

    cReadGnssCodeBias::~cReadGnssCodeBias() {}

    void cReadGnssCodeBias::DecodeCasMgexDcb() {
        cSat sat;
        string code_pair;
        double cbias=0.0;
        int i,j=0;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if((!line_str_.compare(1,3,"DSB"))&&(!line_str_.compare(65,2,"ns"))&&
               (!line_str_.compare(15,4,"    "))){
                sat=cSat(line_str_.substr(11,3));
                sat.SatId2No();
                code_pair=line_str_.substr(25,3)+"-"+line_str_.substr(30,3);
                Str2Double(line_str_.substr(80,10),cbias);
                if(sat.sat_.sys==SYS_GPS){
                    for(i=0;i<MAX_GNSS_CODE_BIAS_PAIRS;i++){
                        if(code_pair==kGnssCodeBiasPairs[SYS_INDEX_GPS][i]) break;
                    }
                    nav_->code_bias[sat.sat_.no-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_.sys==SYS_BDS){
                    for(i=0;i<MAX_GNSS_CODE_BIAS_PAIRS;i++){
                        if(sat.sat_.prn>18){
                            if(code_pair==kGnssCodeBiasPairs[SYS_INDEX_BDS][i]){
                                if(i==1) i=BD3_C2IC6I;
                                break;
                            }
                        }
                        else if(code_pair==kGnssCodeBiasPairs[SYS_INDEX_BDS][i]) break;
                    }
                    nav_->code_bias[sat.sat_.no-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_.sys==SYS_GAL){
                    for(i=0;i<MAX_GNSS_CODE_BIAS_PAIRS;i++){
                        if(code_pair==kGnssCodeBiasPairs[SYS_INDEX_GAL][i]) break;
                    }
                    nav_->code_bias[sat.sat_.no-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_.sys==SYS_GLO){
                    for(i=0;i<MAX_GNSS_CODE_BIAS_PAIRS;i++){
                        if(code_pair==kGnssCodeBiasPairs[SYS_INDEX_GLO][i]) break;
                    }
                    nav_->code_bias[sat.sat_.no-1][i]=cbias*1E-9*CLIGHT;
                }
                else if(sat.sat_.sys==SYS_QZS){
                    for(i=0;i<MAX_GNSS_CODE_BIAS_PAIRS;i++){
                        if(code_pair==kGnssCodeBiasPairs[SYS_INDEX_QZS][i]) break;
                    }
                    nav_->code_bias[sat.sat_.no-1][i]=cbias*1E-9*CLIGHT;
                }
            }
            else continue;
        }
    }

    bool cReadGnssCodeBias::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss code bias file error: "<<file_;
            return false;
        }

        DecodeCasMgexDcb();

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssErp::cReadGnssErp() {}

    cReadGnssErp::cReadGnssErp(string file_path, PPPLib::tNav &nav) {
        file_=file_path;
        nav_=&nav;
    }

    cReadGnssErp::~cReadGnssErp() {}

    void cReadGnssErp::DecodeErpPara() {
        double v[14]={0};
        tErpUnit erp0={0};

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if (sscanf(line_str_.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12,v+13)<5)
                continue;
            erp0.mjd=v[0];
            erp0.xp=v[1]*1E-6*AS2R;
            erp0.yp=v[2]*1E-6*AS2R;
            erp0.ut1_utc=v[3]*1E-7;
            erp0.lod=v[4]*1E-7;
            erp0.xpr=v[12]*1E-6*AS2R;
            erp0.ypr=v[13]*1E-6*AS2R;
            nav_->erp_paras.push_back(erp0);
        }
    }

    bool cReadGnssErp::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss erp file error: "<<file_;
            return false;
        }

        DecodeErpPara();

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssOcean::cReadGnssOcean() {}

    cReadGnssOcean::cReadGnssOcean(string file_path, PPPLib::tNav &nav,string site,RECEIVER_INDEX idx) {
        file_=file_path;
        nav_=&nav;
        site_name_=site;
        index_=idx;
    }

    cReadGnssOcean::~cReadGnssOcean() {}

    void cReadGnssOcean::DecodeOceanPara() {
        string name="    ";
        transform(site_name_.begin(),site_name_.begin()+4,name.begin(),::toupper);

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(line_str_.size()<2||line_str_.compare(0,2,"$$")==0) continue;
            if(line_str_.compare(2,4,name)==0){
                double v[11]={0};
                int n=0;
                while(getline(inf_,line_str_)&&!inf_.eof()){
                    if(line_str_.compare(0,2,"$$")==0) continue;
                    if (sscanf(line_str_.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                               v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10)<11) continue;
                    for(int i=0;i<11;i++) nav_->ocean_paras[index_][n+i*6]=v[i];
                    if(++n==6){return;}
                }
            }
        }
    }

    bool cReadGnssOcean::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss ocean_load file error: "<<file_;
            return false;
        }

        DecodeOceanPara();

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssAntex::cReadGnssAntex() {}

    cReadGnssAntex::cReadGnssAntex(string file_path, PPPLib::tNav &nav) {
        file_=file_path;
    }

    cReadGnssAntex::~cReadGnssAntex() {}

    tAntUnit* cReadGnssAntex::SearchAntPar(cTime t,int sat, const string &type) {
        int i,j;
        tAntUnit *ant=nullptr;
        if(sat){
            for(i=0;i<ant_paras_.size();i++){
                ant=&ant_paras_.at(i);
                if(ant->sat.sat_.no!=sat) continue;
                if(ant->ts.t_.long_time!=0.0&&(ant->ts.TimeDiff(t.t_)>0.0)) continue;
                if(ant->te.t_.long_time!=0.0&&(ant->te.TimeDiff(t.t_)<0.0)) continue;
                return ant;
            }
        }
        else{
            int n=0;
            char *p,*types[2];
            const string& buff=type;
            for(p=strtok((char*)buff.c_str()," ");p&&n<2;p=strtok(nullptr," ")) types[n++]=p;
            if(n<=0) return nullptr;

            for(i=0;i<ant_paras_.size();i++){
                ant=&ant_paras_.at(i);
                for(j=0;j<n;j++) if(!strstr(ant->ant_type.c_str(),types[j])) break;
                if(j>=n) return ant;
            }

            for(i=0;i<ant_paras_.size();i++){
                ant=&ant_paras_.at(i);
                if(strstr(ant->ant_type.c_str(),types[0])!=ant->ant_type) continue;
                return ant;
            }
        }
        return nullptr;
    }

    int cReadGnssAntex::DecodeAntPcv(char* p,int n,double* v) {
        int i;
        for(i=0;i<n;i++) v[i]=0.0;
        for(i=0,p=strtok(p," ");p&&i<n;p=strtok(NULL," ")){
            v[i++]=atof(p)*1E-3;
        }
        return i;
    }

    void cReadGnssAntex::ReadAntBody() {
        int stat=0,frq=0,f,i;
        tAntUnit ant0={0};
        string sys;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(line_str_.length()<60||line_str_.find("COMMENT",60)!=string::npos) continue;
            if(line_str_.find("START OF ANTENNA",60)!=string::npos){
                ant_paras_.push_back(ant0);stat=1;
            }
            if(line_str_.find("END OF ANTENNA",60)!=string::npos) stat=0;
            if(!stat) continue;

            if(line_str_.find("TYPE / SERIAL NO",60)!=string::npos){
                ant_paras_.back().ant_type=line_str_.substr(0,20);
                ant_paras_.back().ser_code=line_str_.substr(20,20);
                ant_paras_.back().ant_type=StrTrim(ant_paras_.back().ant_type);
                if(ant_paras_.back().ser_code.compare(3,8,"        ")==0){
                    if(!ant_paras_.back().ser_code.substr(0,3).compare(0,1,"I")) continue;
                    ant_paras_.back().sat=cSat(ant_paras_.back().ser_code.substr(0,3));
                    ant_paras_.back().sat.SatId2No();
                }
                ant_paras_.back().ser_code=StrTrim(ant_paras_.back().ser_code);
            }
            else if(line_str_.find("VALID FROM",60)!=string::npos){
                if(!ant_paras_.back().ts.Str2Time(line_str_.substr(0,43))) continue;
            }
            else if(line_str_.find("VALID UNTIL",60)!=string::npos){
                if(!ant_paras_.back().te.Str2Time(line_str_.substr(0,43))) continue;
            }
            else if(line_str_.find("DAZI",60)!=string::npos){
                Str2Double(line_str_.substr(2,6),ant_paras_.back().dazi); continue;
            }
            else if(line_str_.find("ZEN1 / ZEN2 / DZEN")!=string::npos) {
                Str2Double(line_str_.substr(2,6),ant_paras_.back().zen1);
                Str2Double(line_str_.substr(8,6),ant_paras_.back().zen2);
                Str2Double(line_str_.substr(14,6),ant_paras_.back().dzen);
                continue;
            }
            else if(line_str_.find("START OF FREQUENCY",60)!=string::npos){
                if(sscanf(line_str_.c_str()+4,"%d",&f)<1) continue;
                sys=line_str_[3];
                if(sys=="G") frq=f;
                else if(sys=="C"){
                    if(f==1) frq=f+1*MAX_GNSS_FRQ_NUM;
                    else if(f==7) frq=2+1*MAX_GNSS_FRQ_NUM;
                    else if(f==6) frq=3+1*MAX_GNSS_FRQ_NUM;
                    //TODO 不同版本的atx文件对BDS定义的频率不一样,注意这里的频率匹配
                }
                else if(sys=="E"){
                    if(f==1) frq=f+2*MAX_GNSS_FRQ_NUM;
                    else if(f==5) frq=2+2*MAX_GNSS_FRQ_NUM;
                    else if(f==6) frq=3+2*MAX_GNSS_FRQ_NUM;
                }
                else if(sys=="R"){
                    frq=f+3*MAX_GNSS_FRQ_NUM;
                }
                else if(sys=="J"){
                    if(f==1) frq=f+4*MAX_GNSS_FRQ_NUM;
                    else if(f==2) frq=2+4*MAX_GNSS_FRQ_NUM;
                    else if(f==5) frq=3+4*MAX_GNSS_FRQ_NUM;
                }
                else frq=0;
            }
            else if(line_str_.find("END OF FREQUENCY",60)!=string::npos){
                frq=0;
            }
            else if(line_str_.find("NORTH / EAST / UP",60)!=string::npos){
                double neu[3]={0};
                if(frq<1) continue;
                if(sscanf(line_str_.c_str(),"%lf %lf %lf",neu,neu+1,neu+2)<3) continue;
                ant_paras_.back().pco[frq-1][0]=1E-3*neu[ant_paras_.back().sat.sat_.no?0:1];
                ant_paras_.back().pco[frq-1][1]=1E-3*neu[ant_paras_.back().sat.sat_.no?1:0];
                ant_paras_.back().pco[frq-1][2]=1E-3*neu[2];
            }
            else if(line_str_.find("NOAZI")!=string::npos){
                if(frq<1) continue;
                double dd=(ant_paras_.back().zen2-ant_paras_.back().zen1)/ant_paras_.back().dzen+1;
                if(dd!=round(dd)||dd<=1){
                    LOG_N_TIMES(1,WARNING)<<"Number of PCV NOAZI parameter decode error";
                    continue;
                }
                if(ant_paras_.back().dazi==0.0){
                    i=DecodeAntPcv(const_cast<char*>(line_str_.c_str()+8),(int)dd,ant_paras_.back().pcv[frq-1]);
                    if(i<=0) {
                        LOG_N_TIMES(1,WARNING)<<"Number of PCV NOAZI parameter decode error";
                        continue;
                    }
                    else if(i!=(int)dd){
                        LOG_N_TIMES(1,WARNING)<<"Number of PCV NOAZI parameter decode error";
                        continue;
                    }
                }
                else{
                    int id=(int)(360-0)/ant_paras_.back().dazi+1;
                    for(i=0;i<id;i++){
                        getline(inf_,line_str_);
                        int j=DecodeAntPcv(const_cast<char*>(line_str_.c_str()+8),(int)dd,&ant_paras_.back().pcv[frq-1][i*(int)dd]);
                        if(j<=0){
                            LOG_N_TIMES(1,WARNING)<<"Number of PCV AZI parameter decode error";
                            continue;
                        }
                        else if(j!=(int)dd){
                            LOG_N_TIMES(1,WARNING)<<"Number of PCV AZI parameter decode error";
                            continue;
                        }
                    }
                }
            }
        }
    }

    void cReadGnssAntex::AlignAntPar2Sat(tPPPLibConf C,cTime t, tStaInfoUnit* sta,tAntUnit *sat_ant, tAntUnit *rec_ant) {
        int i;
        tAntUnit *ant;
        cSat sat;

        for(i=0;i<MAX_SAT_NUM;i++){
            sat=cSat(i+1);
            if(sat.sat_.sys&C.gnssC.nav_sys){
                if(!(ant=SearchAntPar(t,sat.sat_.no,""))){
                    continue;
                }
                sat_ant[i]=*ant;
            }
        }

        int dgnss=(C.mode==MODE_DGNSS||C.mode==MODE_PPK);
        for(i=0;i<dgnss?2:1;i++){
            if(!(ant=SearchAntPar(t,0,sta[i].ant_desc))){
                return;
            }
            else {
                rec_ant[i]=*ant;
                rec_ant[i].rec_ant_del[i]=sta->del;
            }
        }
    }

    bool cReadGnssAntex::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss antex file error: "<<file_;
            return false;
        }

        ReadAntBody();

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadGnssIonex::cReadGnssIonex() {}

    cReadGnssIonex::cReadGnssIonex(string file_path, PPPLib::tNav &nav) {
        file_=file_path;
        nav_=nav;
    }

    cReadGnssIonex::~cReadGnssIonex() {}

    bool cReadGnssIonex::ReadHead() {
        double ver=0;
        if(!inf_.is_open()) return 0;
        while(getline(inf_,line_str_)&&!inf_.eof()){
            if (line_str_.length()<60) continue;

            if (line_str_.find("IONEX VERSION / TYPE")!=string::npos){
                if (line_str_[20]=='I') Str2Double(line_str_.substr(0,8),ver);
            }
            else if (line_str_.find("BASE RADIUS")!=string::npos) {
                Str2Double(line_str_.substr(0,8),re_);
            }
            else if (line_str_.find("HGT1 / HGT2 / DHGT")!=string::npos){
                Str2Double(line_str_.substr(2, 6),hgts_[0]);
                Str2Double(line_str_.substr(8, 6),hgts_[1]);
                Str2Double(line_str_.substr(14,6),hgts_[2]);
            }
            else if (line_str_.find("LAT1 / LAT2 / DLAT")!=string::npos){
                Str2Double(line_str_.substr(2, 6),lats_[0]);
                Str2Double(line_str_.substr(8, 6),lats_[1]);
                Str2Double(line_str_.substr(14,6),lats_[2]);
            }
            else if (line_str_.find("LON1 / LON2 / DLON")!=string::npos) {
                Str2Double(line_str_.substr(2, 6),lons_[0]);
                Str2Double(line_str_.substr(8, 6),lons_[1]);
                Str2Double(line_str_.substr(14,6),lons_[2]);
            }
            else if (line_str_.find("EXPONENT")!=string::npos){
                Str2Double(line_str_.substr(0,6),factor_);
            }
            else if (line_str_.find("START OF AUX DATA")!=string::npos&&
                    line_str_.find("DIFFERENTIAL CODE BIASES")!=string::npos){
                continue;
            }
            else if (line_str_.find("PRN / BIAS / RMS")!=string::npos){
                continue;
            }
            else if (line_str_.find("END OF HEADER")!=string::npos){
                return ver;
            }
        }
    }

    int cReadGnssIonex::DataIndex(int i, int j, int k, const int *ndata) {
        if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
        return i+ndata[0]*(j+ndata[1]*k);
    }

    int cReadGnssIonex::GetIndex(double val, const double *range) {
        if(range[2]==0.0) return 0;
        if(range[2]>0.0&&(val<range[0]||range[1]<val)) return -1;
        if(range[2]<0.0&&(val<range[1]||range[0]<val)) return -1;
        return (int)floor((val-range[0])/range[2]+0.5);
    }

    int cReadGnssIonex::GetNumItems(const double *range) {
        return GetIndex(range[1],range)+1;
    }

    tTecUnit* cReadGnssIonex::AddTec() {
        tTecUnit tec0;
        int num_data[3];

        num_data[0]=GetNumItems(lats_);
        num_data[1]=GetNumItems(lons_);
        num_data[2]=GetNumItems(hgts_);
        if(num_data[0]<=1||num_data[1]<=1||num_data[2]<=0) return NULL;

        nav_.tec_paras.push_back(tec0);
        nav_.tec_paras.back().re=re_;
        for(int i=0;i<3;i++){
            nav_.tec_paras.back().ndata[i]=num_data[i];
            nav_.tec_paras.back().lats[i]=lats_[i];
            nav_.tec_paras.back().lons[i]=lons_[i];
            nav_.tec_paras.back().hgts[i]=hgts_[i];
        }
        int n=num_data[0]*num_data[1]*num_data[2];

        nav_.tec_paras.back().data.assign(n,0.0);
        nav_.tec_paras.back().rms.assign(n,0.0);

        for(int i=0;i<n;i++){
            nav_.tec_paras.back().data[i]=0.0;
            nav_.tec_paras.back().rms[i]=0.0f;
        }
        return &nav_.tec_paras.back();
    }

    void cReadGnssIonex::ReadIonBody() {
        tTecUnit *tec=nullptr;
        int flag=0;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(line_str_.length()<60) continue;
            if(line_str_.find("START OF TEC MAP")!=string::npos){
                tec = AddTec();
                if(tec) flag=1;
            }
            else if(line_str_.find("END OF TEC MAP")!=string::npos){
                flag=0;tec=nullptr;
            }
            else if(line_str_.find("START OF RMS MAP")!=string::npos){
                flag=2;tec=nullptr;
            }
            else if(line_str_.find("END OF RMS MAP")!=string::npos){
                flag=0;tec=nullptr;
            }
            else if(line_str_.find("EPOCH OF CURRENT MAP")!=string::npos){
                if(ion_time_.Str2Time(line_str_.substr(0,36))) continue;
                if(flag==2){
                    for(int i=nav_.tec_paras.size()-1;i>=0;i--){
                        if(fabs(ion_time_.TimeDiff(nav_.tec_paras[i].t_tag.t_))>=1.0) continue;
                        tec=&nav_.tec_paras[i];
                        break;
                    }
                }
                else if(tec) tec->t_tag=ion_time_;
            }
            else if(line_str_.find("LAT/LON1/LON2/DLON/H")!=string::npos&&tec){
                double lat,lon[3],hgt,x;
                Str2Double(line_str_.substr(2,6), lat);
                Str2Double(line_str_.substr(8,6), lon[0]);
                Str2Double(line_str_.substr(14,6),lon[1]);
                Str2Double(line_str_.substr(20,6),lon[2]);
                Str2Double(line_str_.substr(26,6),hgt);
                int i=GetIndex(lat,tec->lats);
                int k=GetIndex(hgt,tec->hgts);
                int n=GetNumItems(lon);

                for (int m=0;m<n;m++) {
                    int index;
                    if (m%16==0&&!getline(inf_,line_str_)) break;
                    int j=GetIndex(lon[0]+lon[2]*m,tec->lons);
                    if ((index=DataIndex(i,j,k,tec->ndata))<0) continue;
                    Str2Double(line_str_.substr(m%16*5,5),x);
                    if (x==9999.0) continue;
                    if (flag==1) tec->data[index]=x*pow(10.0,factor_);
                    else tec->rms[index]=(float)(x*pow(10.0,factor_));
                }
            }
        }
    }

    bool cReadGnssIonex::Reading() {
        if(!OpenFile()){
            LOG(ERROR)<<"Open gnss ionex file error: "<<file_;
            return false;
        }

        if(!ReadHead()) return false;

        ReadIonBody();

        if(OpenFile()) CloseFile();
        return true;
    }

    cReadRefSol::cReadRefSol() {}

    cReadRefSol::cReadRefSol(string file_path,vector<tSolInfoUnit>& ref_sol) {file_=file_path;ref_sols_=&ref_sol;}

    cReadRefSol::~cReadRefSol() {}

    void cReadRefSol::ReadRefIe() {
        int i,gpsw;
        double sec;
        char *seps=(char *) ", ";
        char *token;
        tSolInfoUnit sol;

        while(getline(inf_,line_str_)&&!inf_.eof()){
            if(48>=line_str_[1]||line_str_[1]>=57) continue;
            token=strtok((char *)line_str_.c_str(),seps);
            gpsw=atoi(token);
            token=strtok(NULL,seps);
            sec=atof(token);
            sol.t_tag.Gpst2Time(gpsw,sec,SYS_GPS);

            for(i=0;i<3;i++){
                token=strtok(NULL,seps);
                sol.pos[i]=atof(token);
            }
            for(i=0;i<3;i++){
                token=strtok(NULL,seps);
                sol.vel[i]=atof(token);
            }
            for(i=0;i<3;i++){
                token=strtok(NULL,seps);
                sol.att[i]=atof(token)*D2R;
            }
            for(i=0;i<3;i++){
                if(sol.att[i]<0.0) sol.att[i]+=360*D2R;
            }
            token=strtok(NULL,seps);
            sol.accl_bias[0]=atof(token);
            token=strtok(NULL,seps);
            sol.accl_bias[1]=atof(token);
            token=strtok(NULL,seps);
            sol.accl_bias[2]=atof(token);
            token=strtok(NULL,seps);
            sol.gyro_bias[0]=atof(token);
            token=strtok(NULL,seps);
            sol.gyro_bias[1]=atof(token);
            token=strtok(NULL,seps);
            sol.gyro_bias[2]=atof(token);
            sol.valid_sat_num=0;
            sol.stat=SOL_REF;
            sol.ins_stat=SOL_INS_REF;
            sol.valid_sat_num=0;
            sol.ratio=0.0;
            sol.sigma=0.0;
            sol.age=0.0;
            ref_sols_->push_back(sol);
        }

    }

    void cReadRefSol::ReadRefCsv(TIME_FORMAT time_f) {
        double data[MAX_CSV_COLS]={0};
        string part_str;
        tSolInfoUnit sol;
        int i,week;
        double sow;
        string t_str;
        cTime t_tag;

        getline(inf_,line_str_);
        while(getline(inf_,line_str_)&&!inf_.eof()){

            if(!line_str_.substr(0,3).compare("EOF")) break;

            istringstream read_str(line_str_);

            for(auto &j:data) j=0.0;
            for(i=0;i<cols_&&cols_<MAX_CSV_COLS;i++){

                getline(read_str,part_str,sep_);
                if(i==idx_t_){
                    if(time_f==TIME_FMT_STR){
                        t_str=part_str;
                        t_tag.Str2Time(t_str);
                        continue;
                    }
                    else{
                        week=atof(part_str.c_str());
                        getline(read_str,part_str,sep_);
                        sow=atof(part_str.c_str());
                        t_tag.Gpst2Time(week,sow,SYS_GPS);
                        i++;
                        continue;
                    }
                }
                else if(i>idx_t_){
                    data[i]=atof(part_str.c_str());
                }
            }

            sol.t_tag=t_tag;
            Vector3d blh(data[px_]*D2R,data[py_]*D2R,data[pz_]*D2R);
            sol.pos=Blh2Xyz(blh);
            sol.vel[0]=data[vx_];
            sol.vel[1]=data[vy_];
            sol.vel[2]=data[vz_];
            sol.att[0]=data[ax_]*D2R;
            sol.att[1]=data[ay_]*D2R;
            sol.att[2]=data[az_]*D2R;
            sol.stat=SOL_REF;
            sol.ins_stat=SOL_INS_REF;
            sol.valid_sat_num=0;
            sol.ratio=0.0;
            sol.sigma=0.0;
            sol.age=0.0;
            ref_sols_->push_back(sol);
        }
    }

    vector<tSolInfoUnit> cReadRefSol::GetRefSols() {return *ref_sols_;}

    void cReadRefSol::SetDataIdx(int t, int px, int py, int pz, int vx, int vy, int vz, int ax, int ay, int az, const char sep, int cols) {
        idx_t_=t;
        px_=px;
        py_=py;
        pz_=pz;
        vx_=vx;
        vy_=vy;
        vz_=vz;
        ax_=ax;
        ay_=ay;
        az_=az;
        sep_=sep;
        cols_=cols;
    }

    bool cReadRefSol::Reading(int type) {
        if(!OpenFile()){
            LOG(ERROR)<<"Open reference solution file error: "<<file_;
            return false;
        }

        if(type){
            ReadRefIe();
        }
        else{
            ReadRefCsv(TIME_FMT_WS);
        }

        if(OpenFile()) CloseFile();
        return true;
    }

}







