#include <iostream>
#include "Kalman.h"
#include "FileOperator.h"
#include <chrono>
#include <string>
#include <vector>
#include "quill/Quill.h"
#include "Kalman.h"
#include <eigen3/Eigen/Dense>
#include "cppTypes.h"
#include "FootPlace.h"

double acc[3], eul[3], omegaL[3];
double s, swing, legSptInd[2];
double pB_fk[6];
double pas_delta[4];
int count = 0;
double pCoM[3], vCoM[3], pfW[3], peW[6], pf2cen_ar[3];
double L[3], Lm[3];
double yaw0;

FileOperator fileRW("../test_in.txt");

Angular_Momentum_class AngularMomentum;
FootPlace_class     foot;

Eigen::Matrix<double,3,3> Rx3(double theta);
Eigen::Matrix<double,3,3> Ry3(double theta);
Eigen::Matrix<double,3,3> Rz3(double theta);

double ini[10] = {0.0585,  0.0427-0.125, -0.6, 0.0, 0.0,
                  0.0585, -0.0427+0.125, -0.6, 0.0, 0.0};

int main() {
    std::vector<double> tmpValue;
    std::string tmpStr;
    // Get a pointer to the default logger
    quill::Handler *file_handler = quill::file_handler("../test_out.txt", "w");
    file_handler->set_pattern(QUILL_STRING("%(message)"));
    quill::Logger *dl = quill::create_logger("logger", file_handler);
    quill::start();

    // main loop
    for (int i = 0; i < fileRW.getTotalLine(); i++) {

        // read data
        fileRW.getNewLine();
        fileRW.getNumsInLine();
        eul[0] = fileRW.values[0];
        eul[1] = fileRW.values[1];
        eul[2] = fileRW.values[2];
        acc[0] = fileRW.values[3];
        acc[1] = fileRW.values[4];
        acc[2] = fileRW.values[5];
        omegaL[0] = fileRW.values[6];
        omegaL[1] = fileRW.values[7];
        omegaL[2] = fileRW.values[8];
           s    = fileRW.values[9];
        swing = fileRW.values[10];
        pB_fk[0] = fileRW.values[11];
        pB_fk[1] = fileRW.values[12];
        pB_fk[2] = fileRW.values[13];
        pB_fk[3] = fileRW.values[14];
        pB_fk[4] = fileRW.values[15];
        pB_fk[5] = fileRW.values[16];
        pas_delta[0] = fileRW.values[17];
        pas_delta[1] = fileRW.values[18];
        pas_delta[2] = fileRW.values[19];
        pas_delta[3] = fileRW.values[20];
        yaw0 = fileRW.values[21];

        pCoM[0] = fileRW.values[22];
        pCoM[1] = fileRW.values[23];
        pCoM[2] = fileRW.values[24];
        vCoM[0] = fileRW.values[25];
        vCoM[1] = fileRW.values[26];
        vCoM[2] = fileRW.values[27];
        peW[0] = fileRW.values[28];
        peW[1] = fileRW.values[29];
        peW[2] = fileRW.values[30];
        peW[3] = fileRW.values[31];
        peW[4] = fileRW.values[32];
        peW[5] = fileRW.values[33];
        legSptInd[0] = fileRW.values[34];
        legSptInd[1] = fileRW.values[35];

        // kalman filter
        Eigen::Matrix<double,3,3> R, I, Iw;
        Eigen::Matrix<double,3,1> pf2cen, ww, vw, p_cen_b, p_cen_w, p_com_w;
        R = Rz3(eul[2]-yaw0)*Ry3(eul[1])*Rx3(eul[0]);

        p_cen_b << 0.0455025, 0.000187389, -0.0576864;
        p_com_w << pCoM[0], pCoM[1], pCoM[2];
        p_cen_w = R*p_cen_b  + p_com_w;

        if (swing < 1.1 && swing > 0.9){//((legSwingInd <1.5 & legr < 0.5) | (legSwingInd > 1.5 & legl > 0.5)){
            pfW[0] = peW[3];
            pfW[1] = peW[4];
            pfW[2] = peW[5];
        }
        else if (swing < 2.1 && swing > 1.9){
            pfW[0] = peW[0];
            pfW[1] = peW[1];
            pfW[2] = peW[2];
        }

        pf2cen(0) = p_cen_w(0) - pfW[0];
        pf2cen(1) = p_cen_w(1) - pfW[1];
        pf2cen(2) = p_cen_w(2) - pfW[2];

//        printf("eul[2] = %.6f,  eul[1] = %.6f,  eul[0] = %.6f,  yaw0 = %.6f\n", eul[2], eul[1], eul[0], yaw0);
//        printf("p_com_w = %.6f, %.6f, %.6f\n", pCoM[0], pCoM[1], pCoM[2]);
//        printf("p_cen_w = %.6f, %.6f, %.6f\n", p_cen_w(0), p_cen_w(1), p_cen_w(2));
//        printf("pfW = %.6f, %.6f, %.6f\n", pfW[0], pfW[1], pfW[2]);
//        printf("pf2cen = %.6f, %.6f, %.6f\n", pf2cen(0), pf2cen(1), pf2cen(2));
//        printf("%d\n",i);

        ww = R* (Eigen::Matrix<double, 3, 1>(omegaL[0], omegaL[1], omegaL[2]));
        vw << vCoM[0], vCoM[1], 0;//vCoM[2];
        I << 0.523524, 3.19194e-5, -0.0279889,
                3.19194e-5, 0.393347, -0.000262413,
                -0.0279889, -0.000262413, 0.21856;
        Iw = R*I*R.transpose();
        if (count > 46742 && count < 46744) // need to set the ini state when the robot is steady
            AngularMomentum.init();
        AngularMomentum.cal(pf2cen, ww, vw, Iw);

        //foot
        L[0] = AngularMomentum.L(0);
        L[1] = AngularMomentum.L(1);
        L[2] = AngularMomentum.L(2);
        Lm[0] = AngularMomentum.Lm(0);
        Lm[1] = AngularMomentum.Lm(1);
        Lm[2] = AngularMomentum.Lm(2);
        pf2cen_ar[0] = pf2cen(0);
        pf2cen_ar[1] = pf2cen(1);
        pf2cen_ar[2] = pf2cen(2);

        double   v_des[3] {0}, p_des[3] {0};
        if (swing == 1){
            for(int j = 0; j<3; j++)
                p_des[j] = ini[j];
        }
        else if (swing == 2){
            for(int j = 0; j<3; j++)
                p_des[j] = ini[5+j];
        }

        foot.cal(L, pf2cen_ar, s, v_des, p_des);

        count++;
        // output data
        tmpValue.clear();

        tmpValue.push_back(pCoM[0]);//1
        tmpValue.push_back(pCoM[1]);//2
        tmpValue.push_back(pCoM[2]);//3
        tmpValue.push_back(vCoM[0]);//4
        tmpValue.push_back(vCoM[1]);//5
        tmpValue.push_back(vCoM[2]);//6
        tmpValue.push_back(peW[0]);//7
        tmpValue.push_back(peW[1]);//8
        tmpValue.push_back(peW[2]);//9
        tmpValue.push_back(peW[3]);//10
        tmpValue.push_back(peW[4]);//11
        tmpValue.push_back(peW[5]);//12
        tmpValue.push_back(legSptInd[0]);//13
        tmpValue.push_back(legSptInd[1]);//14

        tmpValue.push_back(AngularMomentum.L[0]);//15
        tmpValue.push_back(AngularMomentum.L[1]);//16
        tmpValue.push_back(AngularMomentum.L[2]);//17
        tmpValue.push_back(pf2cen(0));//18
        tmpValue.push_back(pf2cen(1));//19
        tmpValue.push_back(pf2cen(2));//20
        tmpValue.push_back(vw(0));//21
        tmpValue.push_back(vw(1));//22
        tmpValue.push_back(vw(2));//23
        tmpValue.push_back(AngularMomentum.Lm[0]);//24
        tmpValue.push_back(AngularMomentum.Lm[1]);//25
        tmpValue.push_back(AngularMomentum.Lm[2]);//26

        tmpValue.push_back(foot.psw_Tkf_t[0]);//27
        tmpValue.push_back(foot.psw_Tkf_t[1]);//28
        tmpValue.push_back(foot.psw_Tkf_t[2]);//29
        tmpValue.push_back(foot.L_Tkf_t[0]);//30
        tmpValue.push_back(foot.L_Tkf_t[1]);//31
        tmpValue.push_back(foot.L_Tkf_t[2]);//32
        tmpValue.push_back(foot.L_des_Tk1f[0]);//33
        tmpValue.push_back(foot.L_des_Tk1f[1]);//34
        tmpValue.push_back(foot.L_des_Tk1f[2]);//35

        tmpValue.push_back(foot.L_Tkf_t_1[0]);//36
        tmpValue.push_back(foot.L_Tkf_t_2[0]);//37
        tmpValue.push_back(foot.L_Tkf_t_1[1]);//38
        tmpValue.push_back(foot.L_Tkf_t_2[1]);//39

        tmpValue.push_back(foot.t);//40
        tmpValue.push_back(swing);//41

        tmpValue.push_back(p_cen_w(0));//42
        tmpValue.push_back(p_cen_w(1));//43
        tmpValue.push_back(p_cen_w(2));//44

        tmpValue.push_back(ww(0));//45
        tmpValue.push_back(ww(1));//46
        tmpValue.push_back(ww(2));//47

        tmpValue.push_back(AngularMomentum.Lcom(0));//48
        tmpValue.push_back(AngularMomentum.Lcom(1));//49
        tmpValue.push_back(AngularMomentum.Lcom(2));//50

        tmpStr = fmt::format("{:.5f}", fmt::join(tmpValue, " "));
        LOG_INFO(dl, "{}", tmpStr);
    }

    return 0;
}

Eigen::Matrix<double,3,3> Rx3(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Eigen::Matrix<double,3,3> M;
    M << 1, 0, 0,
            0, cos(theta), -sin(theta),
            0, sin(theta),  cos(theta);
    return M;
}

Eigen::Matrix<double,3,3> Ry3(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Eigen::Matrix<double,3,3> M;
    M << cos(theta), 0, sin(theta),
            0,             1,               0,
            -sin(theta), 0, cos(theta);
    return M;
}

Eigen::Matrix<double,3,3> Rz3(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Eigen::Matrix<double,3,3> M;
    M << cos(theta), -sin(theta), 0,
            sin(theta), cos(theta),  0,
            0,             0,              1;
    return M;
}












