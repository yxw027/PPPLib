//
// Created by cc on 7/26/20.
//

#ifndef PPPLIB_ADJFUNC_H
#define PPPLIB_ADJFUNC_H

#include "CmnFunc.h"

namespace PPPLib{

    class cAdjuster {
    public:
        cAdjuster();
        virtual ~cAdjuster();

    public:
        virtual int Adjustment(VectorXd L,const MatrixXd H,const MatrixXd R,VectorXd& X, MatrixXd& Px,int nl,int nx);

    public:
        VectorXd dx_;
        VectorXd v_;
        MatrixXd Qvv_;
        double unit_weight_STD_=0.0;
    };

    class cLsqAdjuster:public cAdjuster{
    public:
        cLsqAdjuster();
        ~cLsqAdjuster();

    public:
        int Adjustment(VectorXd L,const MatrixXd H,const MatrixXd R,VectorXd& X, MatrixXd& Px,int nl,int nx) override;
    };

    class cKfAdjuster:public cAdjuster{
    public:
        cKfAdjuster();
        ~cKfAdjuster();

    public:
        int Adjustment(VectorXd L,const MatrixXd H,const MatrixXd R,VectorXd& X, MatrixXd& Px,int nl,int nx) override;
    };
}




#endif //PPPLIB_ADJFUNC_H
