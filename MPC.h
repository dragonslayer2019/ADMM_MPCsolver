#pragma once

#include "BlockMatrix.h"
#include "FunctionG.h"
#include <Eigen/Dense>
#include "json.hpp"

template <typename T, int HorizonNum, int SizeX, int SizeU, int SizeYx, int SizeYu>
class MPC_ADMMSolver{
private:
    typedef Eigen::Matrix<T, SizeX, SizeX> MatrixX;
    typedef Eigen::Matrix<T, SizeU, SizeU> MatrixU;
    typedef Eigen::Matrix<T, SizeX, SizeU> MatrixB;
    typedef Eigen::Matrix<T, SizeX, 1> VectorX;
    typedef Eigen::Matrix<T, SizeU, 1> VectorU;
    typedef Eigen::Matrix<T, SizeYx, SizeX> MatrixHx;
    typedef Eigen::Matrix<T, SizeYu, SizeU> MatrixHu;
    typedef BlockVector<T, HorizonNum + 1, SizeX + SizeU> BlockVectorComb;
    typedef BlockVector<T, HorizonNum + 1, SizeX> BlockVectorX;
    typedef BlockVector<T, HorizonNum + 1, SizeU> BlockVectorU;
    const double inf = 1e6;

public:
    /**
     * MPC problem coefficients
     */
    std::array<MatrixX, HorizonNum + 1> Q;
    std::array<MatrixU, HorizonNum + 1> R;
    std::array<VectorX, HorizonNum + 1> L;
    std::array<VectorU, HorizonNum + 1> W;
    std::array<MatrixX, HorizonNum> A;
    std::array<MatrixB, HorizonNum> B;
    std::array<VectorX, HorizonNum> c;
    VectorX x_init;
    std::array<MatrixHx, HorizonNum + 1> Hx;
    std::array<MatrixHu, HorizonNum + 1> Hu;
    std::array<std::array<FunctionG<T>, SizeX + SizeU>, HorizonNum + 1> g;
    int K;
    T costweight;

    /**
     * ADMM coefficients
     */
    BlockLowerBidiagonalMatrix<T, HorizonNum + 1, SizeX, SizeX + SizeU> bA;
    BlockLowerBidiagonalMatrix<T, HorizonNum + 1, SizeX, SizeX> bL;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeYx + SizeYu, SizeX + SizeU> bE;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeX + SizeU, SizeYx + SizeYu> bET;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeX + SizeU, SizeX + SizeU> bG, biG;
    BlockSymmetricTridiagonalMatrix<T, HorizonNum + 1, SizeX> bAiGAT;
    BlockSymmetricDiagonalMatrix<T, HorizonNum + 1, SizeX + SizeU, SizeX + SizeU> barE;
    BlockVectorX bb;
    BlockVectorComb bF;
    T rho = 1;

    /**
     * ADMM new coefficients
     */
    std::array<MatrixX, HorizonNum> barA;
    std::array<MatrixB, HorizonNum> barB;
    std::array<Eigen::Matrix<T, SizeU, SizeX>, HorizonNum> barBT;
    std::array<VectorX, HorizonNum> barc;
    std::array<MatrixHx, HorizonNum + 1> bEx;
    std::array<MatrixHu, HorizonNum + 1> bEu;
    /**
     * LQR backward matrix
     */
    std::array<MatrixX, HorizonNum + 1> Q_lqr;
    std::array<MatrixU, HorizonNum + 1> R_lqr;
    std::array<MatrixU, HorizonNum + 1> G_lqr;
    std::array<MatrixB, HorizonNum + 1> K_lqr;
    std::array<MatrixX, HorizonNum + 1> Fx_lqr;
    std::array<MatrixX, HorizonNum + 1> P_lqr;
    std::array<VectorU, HorizonNum + 1> BTPc;
    std::array<MatrixX, HorizonNum + 1> FTP;
    std::array<MatrixB, HorizonNum + 1> KRT;
    VectorX bar_x_init;

    MPC_ADMMSolver(
        std::array<MatrixX, HorizonNum + 1> _Q,
        std::array<MatrixU, HorizonNum + 1> _R,
        std::array<VectorX, HorizonNum + 1> _L,
        std::array<VectorU, HorizonNum + 1> _W,
        std::array<MatrixX, HorizonNum> _A,
        std::array<MatrixB, HorizonNum> _B,
        std::array<VectorX, HorizonNum> _c,
        VectorX _x_init,
        std::array<MatrixHx, HorizonNum + 1> _Hx,
        std::array<MatrixHu, HorizonNum + 1> _Hu,
        std::array<std::array<FunctionG<T>, SizeX + SizeU>, HorizonNum + 1> _g,
        T _costweight,
        int _K = 250
    ) : Q(_Q), R(_R), L(_L), W(_W), A(_A), B(_B), c(_c), x_init(_x_init), Hx(_Hx), Hu(_Hu), g(_g), K(_K) {}

    void PreScaling1() {
        // P_x, P_u
        std::array<MatrixX, HorizonNum + 1> Px, iPx, PxT, iPxT;
        std::array<MatrixU, HorizonNum + 1> Pu, iPu, PuT, iPuT;
        for(int i = 0; i <= HorizonNum; ++i) {
            Px[i] = Q[i].llt().matrixL();
            Pu[i] = R[i].llt().matrixL();
            iPx[i] = Px[i].inverse();
            iPu[i] = Pu[i].inverse();
            PxT[i] = Px[i].transpose();
            PuT[i] = Pu[i].transpose();
            iPxT[i] = PxT[i].inverse();
            iPuT[i] = PuT[i].inverse();
        }
        bar_x_init = PxT[0] * x_init;

        std::array<Eigen::Matrix<T, SizeX, SizeX>, HorizonNum + 1> LambX;
        std::array<Eigen::Matrix<T, SizeU, SizeU>, HorizonNum + 1> LambU;
        for(int i = 0; i <= HorizonNum; ++i) {
            // \bar{A}
            if (i < HorizonNum) {
                barA[i] = PxT[i + 1] * A[i] * iPxT[i];
            }
            // \bar{b}
            if (i < HorizonNum) {
                barB[i] = PxT[i + 1] * B[i] * iPuT[i];
                barBT[i] = barB[i].transpose();
            }
            // \bar{c}
            if (i < HorizonNum) {
                barc[i] = PxT[i + 1] * c[i];
            }
            // \bar{F}
            bF.v[i] << iPx[i] * L[i], iPu[i] *  W[i];

            // \Lambda_x, \Lambda_u
            LambX[i] = (Hx[i] * iPxT[i]).rowwise().norm().asDiagonal().inverse();
            LambU[i] = (Hu[i] * iPuT[i]).rowwise().norm().asDiagonal().inverse();

            for(int j = 0; j < SizeYx; ++j) g[i][j].Prescaling(LambX[i].diagonal()(j));
            for(int j = SizeX; j < SizeYx + SizeYu; ++j) g[i][j].Prescaling(LambU[i].diagonal()(j - SizeX));

            // (x,u) = barE * barz; w = bE * barz;
            barE.d[i] << iPxT[i], Eigen::Matrix<T, SizeX, SizeU>::Zero(SizeX, SizeU),
                      Eigen::Matrix<T, SizeU, SizeX>::Zero(SizeU, SizeX), iPuT[i];
            bE.d[i] << LambX[i] * Hx[i] * iPxT[i], Eigen::Matrix<T, SizeYx, SizeU>::Zero(SizeX, SizeU),
                      Eigen::Matrix<T, SizeYu, SizeX>::Zero(SizeU, SizeX), LambU[i] * Hu[i] * iPuT[i];
            bET.d[i] = bE.d[i].transpose();
            bEx[i] = LambX[i] * Hx[i] * iPxT[i];
            bEu[i] = LambU[i] * Hu[i] * iPuT[i];
        }
    }

    void ADMMPrework1() {
        Eigen::Matrix<T, SizeX, SizeX> Ix;
        Ix.setIdentity(SizeX, SizeX);
        Eigen::Matrix<T, SizeU, SizeU> Iu;
        Iu.setIdentity(SizeU, SizeU);

        for (int i = HorizonNum; i >= 0; --i) {
            Q_lqr[i] = (rho * bEx[i].transpose() * bEx[i] + Ix) * 0.5;
            R_lqr[i] = (rho * bEu[i].transpose() * bEu[i] + Iu) * 0.5;
        }
        P_lqr[HorizonNum] = Q_lqr[HorizonNum];
        K_lqr[HorizonNum] = MatrixB::Zero(SizeX, SizeU);
        for (int i = HorizonNum - 1; i >= 0; --i) {
            G_lqr[i] = (barBT[i] * P_lqr[i + 1] * barB[i] + R_lqr[i]).inverse();
            K_lqr[i] = -barA[i].transpose() * P_lqr[i + 1] * barB[i] * G_lqr[i].transpose();
            Fx_lqr[i] = barA[i] + barB[i] * K_lqr[i].transpose();
            P_lqr[i] = Fx_lqr[i].transpose() * P_lqr[i + 1] * Fx_lqr[i] + K_lqr[i] * R_lqr[i] * K_lqr[i].transpose() + Q_lqr[i];
            BTPc[i] = -2 * barBT[i] * P_lqr[i + 1] * barc[i];
            FTP[i] = Fx_lqr[i].transpose() * P_lqr[i + 1];
            KRT[i] = K_lqr[i] * R_lqr[i].transpose();
        }
    }

    BlockVectorComb LQR_Solver1(BlockVectorComb f) {
        std::array<VectorU, HorizonNum + 1> r_lqr;
        std::array<VectorX, HorizonNum + 1> E_lqr;
        std::array<VectorX, HorizonNum + 1> Fr_lqr;
        E_lqr[HorizonNum] = f.v[HorizonNum].block(0, 0, SizeX, 1);
        for (int i = HorizonNum - 1; i >= 0; --i) {
            r_lqr[i] = 0.5 * G_lqr[i] * (BTPc[i] - barBT[i] * E_lqr[i + 1] - f.v[i].block(SizeX, 0, SizeU, 1));
            Fr_lqr[i] = barB[i] * r_lqr[i];
            E_lqr[i] = 2 * FTP[i] * (Fr_lqr[i] + barc[i]) + Fx_lqr[i].transpose() * E_lqr[i + 1] + 2 * KRT[i] * r_lqr[i] + f.v[i].block(0, 0, SizeX, 1) + K_lqr[i] * f.v[i].block(SizeX, 0, SizeU, 1);
        }
        BlockVectorComb res;
        res.v[0].block(0, 0, SizeX, 1) = bar_x_init;
        res.v[HorizonNum].block(SizeX, 0, SizeU, 1) = -0.5 * R_lqr[HorizonNum].inverse() * f.v[HorizonNum].block(SizeX, 0, SizeU, 1);
        for (int i = 0; i < HorizonNum; ++i) {
            res.v[i + 1].block(0, 0, SizeX, 1) = Fx_lqr[i] * res.v[i].block(0, 0, SizeX, 1) + Fr_lqr[i] + barc[i];
            res.v[i].block(SizeX, 0, SizeU, 1) = K_lqr[i].transpose() * res.v[i].block(0, 0, SizeX, 1) + r_lqr[i];
        }
        return res;
    }
    
    BlockVectorComb G_Solver(BlockVectorComb barz, BlockVectorComb nu, T rho) {
        BlockVectorComb u = bE * barz - nu;
        for(int i = 0; i <= HorizonNum; ++i) {
            for(int j = 0; j < SizeX + SizeU; ++j) {
                T ui = u.v[i](j);
                u.v[i](j) = g[i][j].Minimizer(ui, rho);
            }
        }
        return u;
    }

    T CalculateCost(BlockVector<T, HorizonNum + 1, SizeX + SizeU> barz) {
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> z = barE * barz, w = bE * barz;
        T res1 = 0, res2 = 0, res3 = 0;
        for(int i = 0; i <= HorizonNum; ++i) {
            VectorX x = z.v[i].block(0, 0, SizeX, 1);
            VectorU u = z.v[i].block(SizeX, 0, SizeU, 1);
            res1 += 0.5 * x.transpose() * (Q[i] * x) + L[i].dot(x);
            VectorU v = R[i] * u;
            res1 += 0.5 * u.transpose() * v + W[i].dot(u);
            // check this line, not equal to the real cost, don't effect result
            // res1 += 0.5 * (L[i][0] * L[i][0] / 1 + L[i][1] * L[i][1] / 1);
            //cout<< i<< ' '<< u.transpose() << ' ' << 0.5 * u.transpose() * v + W[i].dot(u) << ' ' << 0.5 * x.transpose() * (Q[i] * x) <<' '<< L[i].dot(x)<<' '<<0.5 * (L[i][0] * L[i][0] / 1 + L[i][1] * L[i][1] / 1)<<' '<<x[0]<<' '<<x[1]<<' '<<x[2]<<endl;
        }
        for(int j = 0; j < SizeX + SizeU; ++j) {
            std::array<T, HorizonNum + 1> d;
            for(int i = 0; i <= HorizonNum; ++i) {
                T x = w.v[i](j);
                res2 += g[i][j].CostOfQuadraticPart(x);
                d[i] = g[i][j].DistanceOfIndicatorPart(x);
            }
            for(int i = 0; i < HorizonNum; ++i) {
                res3 += costweight * (d[i + 1] - d[i] > 0 ? d[i + 1] - d[i] : 0);
            }
            // cout << res3 << endl;
        }
        return res1 + res2 + res3;
    }

    BlockVectorComb ADMMIteration() {
        int k = 1; // iteration num
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> barz, res;
        BlockVector<T, HorizonNum + 1, SizeYx + SizeYu> w, nu;
        barz.setZero(); w.setZero(); nu.setZero(); 
        res = barz;
        T min_cost = inf;
        std::vector<T> cvec;
        while(k <= K) {
            barz = LQR_Solver1(bF - (bET * (w + nu)) * (rho));
            w = G_Solver(barz, nu, rho);
            nu = nu + w - bE * barz;
            k++;
            T cost = CalculateCost(barz);
            cvec.push_back(cost);
            if(!(k % 10)) {
                if(cost < min_cost) {
                    res = barz;
                    min_cost = cost;
                }
                cout << "Episodes: " << k << ' ' << cost << endl;
            }
        }

        using json = nlohmann::json;
        json j;
        json jcost = json::array();
        for(int i = 0; i < K; ++i) {
            jcost.push_back(cvec[i]);
        }
        j["cost"] = jcost;
        ofstream out("testori.out");
        if(out.is_open()) {
            out << j.dump(4) << endl;
            out.close();
        }

        cout << "Final Cost: " << min_cost << endl;
        return barE * res;
    }

    BlockVector<T, HorizonNum + 1, SizeX + SizeU> solve() {
        PreScaling1();
        ADMMPrework1();
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> res = ADMMIteration();
        return res;
    }

};
