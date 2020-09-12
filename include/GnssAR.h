//
// Created by cc on 8/19/20.
//

#ifndef PPPLIB_GNSSAR_H
#define PPPLIB_GNSSAR_H

#include "CmnFunc.h"

#define SQRT2           1.41421356237309510 /* sqrt(2) */
#define LOG_PI          1.14472988584940017 /* log(pi) */


namespace PPPLib {

    double IntAmbConfFunc(int fix_amb,double float_amb,double sig);

    class cLambda{
    public:
        cLambda();
        ~cLambda();

    private:
        void Gauss(int num_amb,int num_i,int num_j,MatrixXd& Z_trans,MatrixXd& L);
        void Permutations(int num_amb,int num_j,double del,MatrixXd& Z_trans,MatrixXd& L,VectorXd& D);
        int FactorizationLD(int num_amb,MatrixXd &var_amb,VectorXd &D,MatrixXd &L);
        void ReductionZ(int num_amb,VectorXd& D,MatrixXd& L,MatrixXd& Z_trans);
        int SearchFix(int num_amb,int fix_num,VectorXd& sum_var,VectorXd &D,MatrixXd &L,VectorXd& z_float,MatrixXd& z_fix);

    public:
        int IntegerAmb(VectorXd& float_amb,MatrixXd& var_amb,MatrixXd& fix_amb,int num_amb,int num_fix_sols,VectorXd& sum_res);

//    private:
//        Eigen::MatrixXd L_,Z_trans_,z_fix_;
//        Eigen::VectorXd D_,z_float_;
    };

}

#endif //PPPLIB_GNSSAR_H
