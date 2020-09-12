//
// Created by cc on 8/19/20.
//

#include "GnssAR.h"

#define LOOPMAX     10000           /* maximum count of search loop */
#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

namespace PPPLib {

    static double p_gamma(double a, double x, double log_gamma_a);
    static double q_gamma(double a, double x, double log_gamma_a){
        double y,w,la=1.0,lb=x+1.0-a,lc;
        int i;

        if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
        w=exp(-x+a*log(x)-log_gamma_a);
        y=w/lb;
        for (i=2;i<100;i++) {
            lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
            la=lb; lb=lc;
            w*=(i-1-a)/i;
            y+=w/la/lb;
            if (fabs(w/la/lb)<1E-15) break;
        }
        return y;
    }

    static double p_gamma(double a, double x, double log_gamma_a){
        double y,w;
        int i;

        if (x==0.0) return 0.0;
        if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);

        y=w=exp(a*log(x)-x-log_gamma_a)/a;

        for (i=1;i<100;i++) {
            w*=x/(a+i);
            y+=w;
            if (fabs(w)<1E-15) break;
        }
        return y;
    }

    static double f_erfc(double x){
        return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
    }

    double IntAmbConfFunc(int fix_amb,double float_amb,double sig){
        double x,p=1.0;
        int i;

        x=fabs(float_amb-fix_amb);
        for (i=1;i<8;i++) {
            p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
        }
        return p;
    }

    cLambda::cLambda() {}

    cLambda::~cLambda() {}

    int cLambda::FactorizationLD(int num_amb,MatrixXd &var_amb,VectorXd &D,MatrixXd &L) {
        int i,j,k,info=0;
        double a;
        Eigen::MatrixXd A=var_amb;

        for(i=num_amb-1;i>=0;i--){
            if((D[i]=A.data()[i+i*num_amb])<=0.0){info=-1;break;}
            a=sqrt(D[i]);
            for(j=0;j<=i;j++) L.data()[i+j*num_amb]=A.data()[i+j*num_amb]/a;
            for(j=0;j<=i-1;j++) for(k=0;k<=j;k++) A.data()[j+k*num_amb]-=L.data()[i+k*num_amb]*L.data()[i+j*num_amb];
            for(j=0;j<=i;j++) L.data()[i+j*num_amb]/=L.data()[i+i*num_amb];
        }
        return info;
    }

    void cLambda::Permutations(int num_amb, int num_j, double del,MatrixXd& Z_trans,MatrixXd& L,VectorXd& D) {
        int k;
        double eta,lam,a0,a1;

        eta=D[num_j]/del;
        lam=D[num_j+1]*L(num_j+1,num_j)/del;
        D[num_j]=eta*D[num_j+1];D[num_j+1]=del;
        for(k=0;k<=num_j-1;k++){
            a0=L(num_j,k);a1=L(num_j+1,k);
            L(num_j,k)=-L(num_j+1,num_j)*a0+a1;
            L(num_j+1,k)=eta*a0+lam*a1;
        }
        L(num_j+1,num_j)=lam;
        for(k=num_j+2;k<num_amb;k++) SWAP(L(k,num_j),L(k,num_j+1));
        for(k=0;k<num_amb;k++) SWAP(Z_trans(k,num_j),Z_trans(k,num_j+1));

    }

    void cLambda::Gauss(int num_amb,int num_i, int num_j,MatrixXd& Z_trans,MatrixXd& L) {
        int k,mu;

        if((mu=(int)ROUND(L(num_i,num_j)))!=0){
            for(k=num_i;k<num_amb;k++){
                L(k,num_j)-=(double)mu*L(k,num_i);
            }
            for(k=0;k<num_amb;k++){
                Z_trans(k,num_j)-=(double)mu*Z_trans(k,num_i);
            }
        }
    }

    void cLambda::ReductionZ(int num_amb,VectorXd& D,MatrixXd& L,MatrixXd& Z_trans) {
        int i,j,k;
        double del;

        j=num_amb-2;k=num_amb-2;
        while(j>=0){
            if(j<=k) for(i=j+1;i<num_amb;i++) Gauss(num_amb,i,j,Z_trans,L);
            del=D[j]+L(j+1,j)*L(j+1,j)*D[j+1];
            if(del+1E-6<D[j+1]){
                Permutations(num_amb,j,del,Z_trans,L,D);
                k=j;j=num_amb-2;
            }
            else j--;
        }
    }

    int cLambda::SearchFix(int num_amb,int num_fix_sols, VectorXd &sum_res,VectorXd &D,MatrixXd &L,VectorXd& z_float,MatrixXd& z_fix) {
        int i,j,k,c,nn=0,imax=0;
        double new_dist,max_dist=1E99,y;
        Eigen::MatrixXd S;
        Eigen::VectorXd dist,zb,z,step;

        S=Eigen::MatrixXd::Zero(num_amb,num_amb);
        zb=z=step=dist=Eigen::VectorXd::Zero(num_amb);

        k=num_amb-1; dist[k]=0.0;
        zb[k]=z_float[k];
        z[k]=ROUND(zb[k]);y=zb[k]-z[k];step[k]=SGN(y); // steps towards closest integer

        for(c=0;c<LOOPMAX;c++){
            new_dist=dist[k]+y*y/D[k];
            if(new_dist<max_dist){
                if(k!=0){
                    dist[--k]=new_dist;
                    for(i=0;i<=k;i++){
                        S(k,i)=S(k+1,i)+(z[k+1]-zb[k+1])*L(k+1,i);
                    }
                    zb[k]=z_float[k]+S(k,k);
                    z[k]=ROUND(zb[k]);y=zb[k]-z[k];step[k]=SGN(y);
                }
                else{
                    if(nn<num_fix_sols){
                        if(nn==0||new_dist>sum_res[imax]) imax=nn;
                        for(i=0;i<num_amb;i++) z_fix(i,nn)=z[i];
                        sum_res[nn++]=new_dist;
                    }
                    else{
                        if(new_dist<sum_res[imax]){
                            for(i=0;i<num_amb;i++) z_fix(i,imax)=z[i];
                            sum_res[imax]=new_dist;
                            for(i=imax=0;i<num_fix_sols;i++) if(sum_res[imax]<sum_res[i]) imax=i;
                        }
                        max_dist=sum_res[imax];
                    }
                    z[0]+=step[0];y=zb[0]-z[0];step[0]=-step[0]-SGN(step[0]);
                }
            }
            else{
                if(k==num_amb-1) break;
                else{
                    k++;
                    z[k]+=step[k];y=zb[k]-z[k];step[k]=-step[k]-SGN(step[k]);
                }
            }
        }

        for(i=0;i<num_fix_sols-1;i++){
            for(j=i+1;j<num_fix_sols;j++){
                if(sum_res[i]<sum_res[j]) continue;
                SWAP(sum_res[i],sum_res[j]);
                for(k=0;k<num_amb;k++) SWAP(z_fix(k,i),z_fix(k,j));
            }
        }

        if(c>=LOOPMAX){
            LOG(WARNING)<<"SEARCH LOOP COUNT OVERFLOW";
            return -1;
        }
        return 0;
    }

    // float_amb:   float DD_ambiguity vector
    // var_amb:     covariance matrix of float ambiguity
    // fix_amb:     fixed ambiguity (num_amb,num_fix_sols)
    // num_amb:     number of float parameters
    // num_fix_sols:number of fixed solutions
    // sum_res:     sum of squared residuals of fixed solutions
    int cLambda::IntegerAmb(VectorXd &float_amb, MatrixXd &var_amb, MatrixXd &fix_amb, int num_amb, int num_fix_sols, VectorXd &sum_res) {
        int state=-1;

        if(num_amb<=0||num_fix_sols<=0) return -1;
        MatrixXd L= MatrixXd::Zero(num_amb,num_amb);
        VectorXd D= VectorXd::Zero(num_amb);
        MatrixXd Z_trans= MatrixXd::Identity(num_amb,num_amb);
        VectorXd z_float= VectorXd::Zero(num_amb);
        MatrixXd z_fix= MatrixXd::Zero(num_amb,num_fix_sols);

        // LD(lower diaganol) factorization(Q=L'*diag(D)*L)
        if(!(state=FactorizationLD(num_amb,var_amb,D,L))){
//            cout<<L<<endl<<endl;
//            cout<<D<<endl<<endl;
            ReductionZ(num_amb,D,L,Z_trans);
            z_float=Z_trans.transpose()*float_amb;

            if(!(state=SearchFix(num_amb,num_fix_sols,sum_res,D,L,z_float,z_fix))){
                fix_amb=Z_trans.transpose().inverse()*z_fix;
            }
        }

        return state;
    }

}