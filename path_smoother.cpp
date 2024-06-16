#include <bits/stdc++.h>
#include <sys/time.h>

#include "BlockMatrix.h"
#include "MPC.h"
#include "FunctionG.h"
#include <cmath>

#include "json.hpp"

using namespace std;
//double shift test
const int SizeX = 6, SizeU = 3;
const int SizeYx = 6, SizeYu = 1;
const int HorizonNum = 49;
const double pi = M_PI;

typedef Eigen::Matrix<double, SizeX, SizeX> MatrixX;
typedef Eigen::Matrix<double, SizeU, SizeU> MatrixU;
typedef Eigen::Matrix<double, SizeX, SizeU> MatrixB;
typedef Eigen::Matrix<double, SizeX, 1> VectorX;
typedef Eigen::Matrix<double, SizeU, 1> VectorU;
typedef Eigen::Matrix<double, SizeYx, SizeX> MatrixHx;
VectorX x_init;
MatrixX Ai, Qi;
MatrixB Bi;
MatrixHx Hxi;
MatrixU Ri, Hui;
VectorX ci, Li;
VectorU Wi;

const double inf = 1e5;

// std::array<double, 3> aaa, bbb, ccc;
// std::array<double, 4> ppp;
std::vector<double> aaa, bbb, ccc;
std::vector<double> ppp;

// void solveunit3D(vector<vector<double>> center, vector<double> dt, vector<double> px, vector<double> Vx, vector<double> Amin, vector<double> Amax, vector<double> Xsafe1, vector<double> Xsafe2, vector<double> Vlim, vector<double> Vlaw, vector<double> B1, vector<double> Theta, vector<double> Ella, vector<double> Ellb, vector<double> lamb1, vector<double> lamb2, vector<double> lamb3, vector<double> lamb4, vector<double> lamb5, vector<double> lamb6, vector<double> lamb7, vector<double> lamb8, vector<double> lamb9, vector<double> lamb10, int K = 250) {
void solveunit3D(vector<double> dt, vector<double> Px, vector<double> Py, vector<double> Pz, vector<double> Vx, vector<double> Vy, vector<double> Vz, vector<double> lamb1, vector<double> lamb2, vector<double> lamb3, vector<double> lamb4, vector<double> lamb5, vector<vector<Hyperplane>> CorridorP, vector<vector<Ellipsoid>> CorridorE, int K = 250) {
    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    std::array<MatrixHx, HorizonNum + 1> Hx;
    std::array<MatrixU, HorizonNum + 1> Hu;

    std::array<MatrixU, HorizonNum> Rk;
    std::array<double, HorizonNum+1> Vnorm;
    // 状态变量(px, py, pz, vx, vy, vz)^T
    // 控制输入(mu1, mu2, mu3)^T
    // std::array<Eigen::Matrix<double, SizeYx - SizeX, 1>, HorizonNum + 1> new_center;
    VectorX x_init;
    std::array<std::array<FunctionG<double>, SizeYx + SizeYu>, HorizonNum + 1> g;
    x_init << 0, 0, 0, 0;
    
    for(int i = 0; i <= HorizonNum; ++i) {
        // 此处为m个切平面约束+一个椭球约束
        // 计算优化路点到切平面之间的距离有个常数跟着怎么处理
        Hxi << 1, 0, 0, 0, 0, 0,
                      0, 1, 0, 0, 0, 0,
                      0, 0, 1, 0, 0, 0,
                      0, 0, 0, 1, 0, 0;
               // ...m行到切平面距离
               // 1行椭球到圆心约束
        
            //    ,
            //    cos(Theta[i]) / Ella[i], sin(Theta[i]) / Ella[i], 0, 0,
            //    -sin(Theta[i]) / Ellb[i], cos(Theta[i]) / Ellb[i], 0, 0;
        // 此处为一个曲率平方约束+一个曲率罚函数   ck = sqrt(mu2^2+mu3^2)
        Hui << 0, 1, 1,
                      0, 1, 1;
        Hx[i] = Hxi;
        Hu[i] = Hui;
        // new_center[i] << cos(Theta[i]) / Ella[i] * center[i][0] + sin(Theta[i]) / Ella[i] * center[i][1],
        //                  -sin(Theta[i]) / Ellb[i] * center[i][0] + cos(Theta[i]) / Ellb[i] * center[i][1];
    }
    // 根据采样参考点计算n-1个Rk，Rk第一列单位向量与vk平行，且为正交矩阵

    // 根据Vx， Vy， Vz计算Vnorm[i]

    for(int i = 0; i < HorizonNum; ++i) {
        // 6*6
        Ai << 1, 0, 0, dt[i], 0, 0,
                   0, 1, 0, 0,  dt[i], 0,
                   0, 0, 1 , 0, 0, dt[i],
                   0, 0, 0 , 1, 0, 0,
                   0, 0, 0 , 0, 1, 0,
                   0, 0, 0 , 0, 0, 1;
        // 6*3
        Bi << 0.5*Vnorm[i][i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][0][0], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][0][1], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][0][2],
                   0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][1][0], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][1][1], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][1][2],
                   0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][2][0], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][2][1], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*dt[i]*Rk[i][2][2],
                   Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][0][0], Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][0][1], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][0][2],
                   Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][1][0], Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][1][1], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][1][2],
                   Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][2][0], Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][2][1], 0.5*Vnorm[i]*Vnorm[i]*dt[i]*Rk[i][2][2];
        // 6*1
        ci << 0, 0, 0, 0, 0, 0;
        A[i] = Ai; B[i] = Bi; c[i] = ci;
    }

    for(int i = 0; i <= HorizonNum; ++i) {
        // 追参考位置
        Qi << lamb1[i], 0, 0, 0, 0, 0,
                   0, lamb1[i], 0, 0, 0, 0,
                   0, 0, lamb1[i], 0, 0, 0,
                   0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0;
        Li << -Px[i] * lamb1[i],
                   -Py[i] * lamb1[i],
                   -Pz[i] * lamb1[i],
                   0,
                   0,
                   0;
        // 限制纵向加速度
        Ri << 0, 0, 0,
                   0, 0, 0,
                   0, 0, lamb2[i];
        Wi << 0,
                    0, 
                    0;
        Q[i] = Qi; R[i] = Ri; L[i] = Li; W[i] = Wi;
    }

    // 设置状态变量相关的罚函数形状
    for(int i = 0;i <= 3; ++i) {
        aaa.push_back(0);
        bbb.push_back(0);
        ccc.push_back(0);
        ppp.push_back(0);
    }
    for(int i = 0; i <= HorizonNum; ++i) {s
        g[i][2].AddIndicator(-inf, inf);
        g[i][3]. (Amin[i], Amax[i]);
        
        double lamb = lamb5[i] / 100 * 0;
        double Xsafem = (Xsafe1[i] + Xsafe2[i]) / 2;
        ppp[0] = -inf;
        aaa[0] = (lamb5[i] + lamb); bbb[0] = (lamb5[i] + lamb) * (-2) * Xsafe1[i]; ccc[0] = (lamb5[i] + lamb) * Xsafe1[i] * Xsafe1[i] + lamb * (Xsafem - Xsafe1[i]) * (Xsafem - Xsafe1[i]);
        ppp[1] = Xsafe1[i];
        aaa[1] = lamb; bbb[1] = lamb * (-2) * Xsafem; ccc[1] = lamb * Xsafem * Xsafem; 
        ppp[2] = Xsafe2[i];
        aaa[2] = (lamb5[i] + lamb); bbb[2] = (lamb5[i] + lamb) * (-2) * Xsafe2[i]; ccc[2] = (lamb5[i] + lamb) * Xsafe2[i] * Xsafe2[i] + lamb * (Xsafem - Xsafe2[i]) * (Xsafem - Xsafe2[i]);
        ppp[3] = inf;
        g[i][0].AddQuadratic(3, aaa, bbb, ccc, ppp);

        double lim1, lim2, weig1, weig2;
        if(Vlim[i] < Vlaw[i]) {
            lim1 = Vlim[i]; lim2 = Vlaw[i];
            weig1 = lamb6[i]; weig2 = lamb8[i];
        } else {
            lim2 = Vlim[i]; lim1 = Vlaw[i];
            weig2 = lamb6[i]; weig1 = lamb8[i];
        }
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0; ppp[0] = -inf;
        aaa[1] = weig1; bbb[1] = weig1 * (-2) * lim1; ccc[1] = weig1 * lim1 * lim1; ppp[1] = lim1;
        aaa[2] = weig1 + weig2;
        bbb[2] = weig1 * (-2) * lim1 + weig2 * (-2) * lim2;
        ccc[2] = weig1 * lim1 * lim1 + weig2 * lim2 * lim2;
        ppp[2] = lim2;
        ppp[3] = inf;
        g[i][2].AddQuadratic(3, aaa, bbb, ccc, ppp);

        ppp[0] = -inf;
        aaa[0] = lamb7[i]; bbb[0] = lamb7[i] * (-2) * B1[i]; ccc[0] = lamb7[i] * B1[i] * B1[i];
        aaa[1] = 0; bbb[1] = 0; ccc[1] = 0; ppp[1] = B1[i]; ppp[2] = inf;
        g[i][3].AddQuadratic(2, aaa, bbb, ccc, ppp);

        ppp[0] = 0;
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0; ppp[1] = 0.5;
        aaa[1] = lamb9[i]; bbb[1] = 0; ccc[1] = (0 - lamb9[i]) * 0.25; ppp[2] = 0.8;
        aaa[2] = lamb10[i]; bbb[2] = 0; ccc[2] = (lamb9[i] - lamb10[i]) * 0.64 + (0 - lamb9[i]) * 0.25; ppp[3] = inf;
        g[i][4].AddQuadratic(3, aaa, bbb, ccc, ppp);
    }


    MPC_ADMMSolver<double, HorizonNum, SizeX, SizeU, SizeYx, SizeYu> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, new_center, g, 100, K);
    BlockVector<double, HorizonNum + 1, SizeX + SizeU> res = mpc.solve();
    

    using json = nlohmann::json;
    json j;
    json jx = json::array();
    json jy = json::array();
    json jv = json::array();
    json ja = json::array();
    json ju = json::array();
    // vector<double> t(HorizonNum), x(HorizonNum), v(HorizonNum), a(HorizonNum), u(HorizonNum);
    for(int i = 0; i <= HorizonNum; ++i) {
        jx.push_back(res.v[i](0, 0));
        jy.push_back(res.v[i](1, 0));
        jv.push_back(res.v[i](2, 0));
        ja.push_back(res.v[i](3, 0));
        ju.push_back(res.v[i](4, 0));
    }
    j["x"] = jx;
    j["y"] = jy;
    j["v"] = jv;
    j["a"] = ja;
    j["u"] = ju;
    j["Px"] = Px;
    ofstream out("test1.out");
    if(out.is_open()) {
        out << j.dump(4) << endl;
        out.close();
    }
    //res.print("RESULT");
    return;
}

int main() {
    double q=0.35, st=0, wei=10, weig=50; int K=250;
    cin >> K;
    struct timeval T1,T2;
    double timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<double> dt, Px, Vx, Amax, Amin, Vcon, Vlaw, Xsafe1, Xsafe2, B1, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, lamb7, lamb8, lamb9, lamb10;
    vector<vector<double>> center;
    vector<vector<Hyperplane>> CorridorP;//分段切平面约束
    vector<vector<Ellipsoid>> CorridorE;//分段

    //0.获取凸走廊结果，一个二维vector，分段的表示全部切平面约束
    
    //1.基于曲率采样HorizonNum段弧长，得到HorizonNum+1个采样路点

    //2.根据弧长vector与对应的V求解非均衡dt

    //3.将采样路点与分段凸走廊匹配，在第i段有m个凸走廊约束就有g[i-1][m-1],让Hxi正确匹配到对应的凸走廊约束



    

    for(int i = 0;i <= HorizonNum; ++i) {
        // dt.push_back(0.2);
        // Px.push_back(0);
        // Vx.push_back(0);
        // Amax.push_back(10);
        // Amin.push_back(-10);

        // 4.将采样路点的位置作为Px、Py、Pz

        // 5.将采样路点的速度作为Vx、Vy、Vz

        // 
        // Vcon.push_back(25);
        // Vlaw.push_back(30);
        // B1.push_back(-inf);
        // Theta.push_back(pi / 3);
        // center.push_back({2., 1.});
        // Ella.push_back(1);
        // Ellb.push_back(0.5);
        lamb1.push_back(1);
        lamb2.push_back(1);
        // lamb3.push_back(3);
        // lamb4.push_back(2);
        // lamb5.push_back(3);
        // lamb6.push_back(100);
        // lamb7.push_back(0.6);
        // lamb8.push_back(3000);
        lamb3.push_back(1);
        lamb4.push_back(1);
        lamb5.push_back(1);
        // lamb6.push_back(1);
        // lamb7.push_back(0.6);
        // lamb8.push_back(3000);
        // if (i < 20 && i >= 11 ){
        //     Xsafe1.push_back(1.7);
        //     Xsafe2.push_back(2.2);
        // } else if (i < 40 && i >= 31) {
        //     Xsafe1.push_back(-1.2);
        //     Xsafe2.push_back(-0.5);
        // } else {
        //     Xsafe1.push_back(-inf);
        //     Xsafe2.push_back(inf);
        // }
        // if (i < 14 && i >= 9) {
        //     lamb9.push_back(10);
        //     lamb10.push_back(1000);
        // } else {
        //     lamb9.push_back(0);
        //     lamb10.push_back(0);
        // }
    }
    solveunit3D(center, dt, Px, Vx, Amin, Amax, Xsafe1, Xsafe2, Vcon, Vlaw, B1, Theta, Ella, Ellb, lamb1, lamb2, lamb3, lamb4, lamb5, CorridorP, CorridorE, K);
    gettimeofday(&T2,NULL);
    timeuse = (T2.tv_sec - T1.tv_sec) + (double)(T2.tv_usec - T1.tv_usec)/1000000.0;
    cout<<"time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
    return 0;
}