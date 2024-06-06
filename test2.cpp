#include <bits/stdc++.h>
#include <sys/time.h>

#include "BlockMatrix.h"
#include "MPC.h"
#include "FunctionG.h"

#include "json.hpp"

using namespace std;
const int SizeX = 3, SizeU = 1;
const int HorizonNum = 29;

typedef Eigen::Matrix<double, SizeX, SizeX> MatrixX;
typedef Eigen::Matrix<double, SizeU, SizeU> MatrixU;
typedef Eigen::Matrix<double, SizeX, SizeU> MatrixB;
typedef Eigen::Matrix<double, SizeX, 1> VectorX;
typedef Eigen::Matrix<double, SizeU, 1> VectorU;
VectorX x_init;
MatrixX Ai, Qi, Hxi;
MatrixB Bi;
MatrixU Ri, Hui;
VectorX ci, Li;
VectorU Wi;
const double inf = 1e5;

// std::array<double, 3> aaa, bbb, ccc;
// std::array<double, 4> ppp;
std::vector<double> aaa, bbb, ccc;
std::vector<double> ppp;

void solveunit(vector<double> dt, vector<double> Xref, vector<double> Vref, vector<double> Amin, vector<double> Amax, vector<double> Xsafe, vector<double> Vlim, vector<double> Vlaw, vector<double> B1, vector<double> lamb1, vector<double> lamb2, vector<double> lamb3, vector<double> lamb4, vector<double> lamb5, vector<double> lamb6, vector<double> lamb7, vector<double> lamb8, int K = 250) {
    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    std::array<MatrixX, HorizonNum + 1> Hx;
    std::array<MatrixU, HorizonNum + 1> Hu;
    VectorX x_init;
    std::array<std::array<FunctionG<double>, SizeX + SizeU>, HorizonNum + 1> g;
    x_init << 0, 20, 0;
    Hxi << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;
    Hui << 1;
    for(int i = 0; i <= HorizonNum; ++i) {
        Hx[i] = Hxi; Hu[i] = Hui;
    }
    for(int i = 0; i < HorizonNum; ++i) {
        Ai << 1, dt[i], dt[i] * dt[i] * 0.5, 
              0, 1 ,    dt[i],
              0, 0 ,    1;
        Bi << 0,
              dt[i] * dt[i] * 0.5,
              dt[i];
        ci << 0, 0, 0;
        A[i] = Ai; B[i] = Bi; c[i] = ci;
    }
    for(int i = 0; i <= HorizonNum; ++i) {
        Qi << lamb1[i], 0, 0,
              0, lamb2[i], 0,
              0, 0, lamb3[i];
        Ri << lamb4[i];
        Li << -Xref[i] * lamb1[i],
              -Vref[i] * lamb2[i],
              0;
        Wi << 0;
        Q[i] = Qi; R[i] = Ri; L[i] = Li; W[i] = Wi;
    }
    for(int i = 0;i < 3; ++i) {
        aaa.push_back(0);
        bbb.push_back(0);
        ccc.push_back(0);
        ppp.push_back(0);
    }
    for(int i = 0; i <= HorizonNum; ++i) {
        g[i][1].AddIndicator(0, inf);
        g[i][2].AddIndicator(Amin[i], Amax[i]);
        
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0;
        aaa[1] = lamb5[i]; bbb[1] = lamb5[i] * (-2) * Xsafe[i]; ccc[1] = lamb5[i] * Xsafe[i] * Xsafe[i];
        ppp[0] = 0; ppp[1] = Xsafe[i]; ppp[2] = inf;
        // for(int i = 0;i < 3; ++i) cout << ppp[i] <<' ';
        g[i][0].AddQuadratic(2, aaa, bbb, ccc, ppp);


        double lim1, lim2, weig1, weig2;
        if(Vlim[i] < Vlaw[i]) {
            lim1 = Vlim[i]; lim2 = Vlaw[i];
            weig1 = lamb6[i]; weig2 = lamb8[i];
        } else {
            lim2 = Vlim[i]; lim1 = Vlaw[i];
            weig2 = lamb6[i]; weig1 = lamb8[i];
        }
        aaa[0] = 0; bbb[0] = 0; ccc[0] = 0; ppp[0] = 0;
        aaa[1] = weig1; bbb[1] = weig1 * (-2) * lim1; ccc[1] = weig1 * lim1 * lim1; ppp[1] = lim1;
        aaa[2] = weig1 + weig2;
        bbb[2] = weig1 * (-2) * lim1 + weig2 * (-2) * lim2;
        ccc[2] = weig1 * lim1 * lim1 + weig2 * lim2 * lim2;
        ppp[2] = lim2;
        ppp[3] = inf;
        g[i][1].AddQuadratic(3, aaa, bbb, ccc, ppp);

        ppp[0] = -inf;
        aaa[0] = lamb7[i]; bbb[0] = lamb7[i] * (-2) * B1[i]; ccc[0] = lamb7[i] * B1[i] * B1[i];
        aaa[1] = 0; bbb[1] = 0; ccc[1] = 0; ppp[1] = B1[i]; ppp[2] = inf;
        g[i][2].AddQuadratic(2, aaa, bbb, ccc, ppp);
    }


    
    MPC_ADMMSolver<double, HorizonNum, SizeX, SizeU, SizeX, SizeU> mpc(Q, R, L, W, A, B, c, x_init, Hx, Hu, g, 100, K);
    BlockVector<double, HorizonNum + 1, SizeX + SizeU> res = mpc.solve();
    

    using json = nlohmann::json;
    json j;
    json jx = json::array();
    json jv = json::array();
    json ja = json::array();
    json ju = json::array();
    // vector<double> t(HorizonNum), x(HorizonNum), v(HorizonNum), a(HorizonNum), u(HorizonNum);
    for(int i = 0; i <= HorizonNum; ++i) {
        jx.push_back(res.v[i](0, 0));
        jv.push_back(res.v[i](1, 0));
        ja.push_back(res.v[i](2, 0));
        ju.push_back(res.v[i](3, 0));
    }
    j["x"] = jx;
    j["v"] = jv;
    j["a"] = ja;
    j["u"] = ju;
    j["Xref"] = Xref;
    ofstream out("test.out");
    if(out.is_open()) {
        out << j.dump(4) << endl;
        out.close();
    }
    // res.print("RESULT");
    return ;
}

int main() {
    double q=0.35, st=0, wei=10, weig=50; int K=250;
    cin >> K;
    struct timeval T1,T2;
    double timeuse;
    gettimeofday(&T1,NULL);
    // For test: 10, 10, 250
    // solve(0.2, 0.2, 2.0, 2.0/3, q, 0.7, 1.2, -0.5, -1.2, 10, 18, 30, 38, st, wei, weig, K);
    vector<double> dt, Xref, Vref, Amax, Amin, Vcon, Vlaw, Xsafe, B1, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, lamb7, lamb8;
    for(int i = 0; i < 25; ++i) {
        dt.push_back(0.2);
        Xref.push_back(20 + i * 0.2 * 23);
    }
    Xref.push_back(Xref[24] + 0.2 * 23);
    for(int i = 25; i < HorizonNum; ++i) {
        dt.push_back(0.5);
        Xref.push_back(Xref[i] + 0.5 * 23);
    }
    for(int i = 0;i <= HorizonNum; ++i) {
        Vref.push_back(23);
        Amax.push_back(10);
        Amin.push_back(-2);
        Vcon.push_back(25);
        Vlaw.push_back(30);
        B1.push_back(-0.4);
        lamb1.push_back(1.0);
        lamb2.push_back(5);
        // lamb3.push_back(3);
        // lamb4.push_back(2);
        // lamb5.push_back(3);
        // lamb6.push_back(100);
        // lamb7.push_back(0.6);
        // lamb8.push_back(3000);
        lamb3.push_back(0.2);
        lamb4.push_back(0.2);
        lamb5.push_back(0);
        lamb6.push_back(1);
        lamb7.push_back(0.6);
        lamb8.push_back(3000);
        Xsafe.push_back(inf);
    }
    Xsafe[1] = 20; Xsafe[2] = 40; Xsafe[3] = 60; Xsafe[4] = 80; Xsafe[5] = 100;
    solveunit(dt, Xref, Vref, Amin, Amax, Xsafe, Vcon, Vlaw, B1, lamb1, lamb2, lamb3, lamb4, lamb5, lamb6, lamb7, lamb8, K);
    gettimeofday(&T2,NULL);
    timeuse = (T2.tv_sec - T1.tv_sec) + (double)(T2.tv_usec - T1.tv_usec)/1000000.0;
    cout<<"time = "<<timeuse<<endl;  //输出时间（单位：ｓ）
    return 0;
}