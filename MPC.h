#pragma once

#include "BlockMatrix.h"
#include "FunctionG.h"
#include <Eigen/Dense>

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

    void PreScaling() {
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

        typedef Eigen::DiagonalMatrix<T, SizeX> DiagMatrixX;
        typedef Eigen::DiagonalMatrix<T, SizeU> DiagMatrixU;
        
        std::array<Eigen::Matrix<T, SizeX, SizeX>, HorizonNum + 1> LambX;
        std::array<Eigen::Matrix<T, SizeU, SizeU>, HorizonNum + 1> LambU;
        std::array<T, SizeX + SizeU> Lamb;
        MatrixX IA; IA.setIdentity(SizeX, SizeX);
        for(int i = 0; i <= HorizonNum; ++i) {
            // \bar{A}
            bA.d[i] << IA, Eigen::Matrix<T, SizeX, SizeU>::Zero(SizeX, SizeU);
            if(i < HorizonNum) bA.sd[i] << -PxT[i] * A[i] * iPxT[i], -PxT[i] * B[i] * iPuT[i];
            
            // \bar{b}
            if(i == 0) bb.v[i] = PxT[0] * x_init;
            else bb.v[i] = PxT[i - 1] * c[i - 1];

            // \bar{F}
            bF.v[i] << iPx[i] * L[i], iPu[i] *  W[i];

            // \Lambda_x, \Lambda_u
            LambX[i] = (Hx[i] * iPxT[i]).rowwise().norm().asDiagonal().inverse();
            LambU[i] = (Hu[i] * iPuT[i]).rowwise().norm().asDiagonal().inverse();

            for(int j = 0; j < SizeYx; ++j) g[i][j].Prescaling(LambX[i].diagonal()(j));
            for(int j = SizeX; j < SizeYx + SizeYu; ++j) g[i][j].Prescaling(LambU[i].diagonal()(j - SizeX));
            

            // (x,u) = barE * barz; Calculate bE = \bar{E}
            barE.d[i] << iPxT[i], Eigen::Matrix<T, SizeX, SizeU>::Zero(SizeX, SizeU),
                      Eigen::Matrix<T, SizeU, SizeX>::Zero(SizeU, SizeX), iPuT[i];
            bE.d[i] << LambX[i] * Hx[i] * iPxT[i], Eigen::Matrix<T, SizeYx, SizeU>::Zero(SizeX, SizeU),
                      Eigen::Matrix<T, SizeYu, SizeX>::Zero(SizeU, SizeX), LambU[i] * Hu[i] * iPuT[i];
        }
    }

    void ADMMPrework() {
        Eigen::Matrix<T, SizeX + SizeU, SizeX + SizeU> I;
        I.setIdentity(SizeX + SizeU, SizeX + SizeU);
        for(int i = 0; i <= HorizonNum; ++i) {
            bET.d[i] = bE.d[i].transpose();
            bG.d[i] = rho * bET.d[i] * bE.d[i] + I;
            biG.d[i] = bG.d[i].inverse();
        }
        
        BlockLowerBidiagonalMatrix<T, HorizonNum + 1, SizeX, SizeX + SizeU> bAiG;
        for(int i = 0; i <= HorizonNum; ++i) {
            bAiG.d[i] = bA.d[i] * biG.d[i];
            if(i < HorizonNum) bAiG.sd[i] = bA.sd[i] * biG.d[i];
        }
        bAiGAT.d[0] = bAiG.d[0] * bA.d[0].transpose();
        bAiGAT.sd[0] = bAiG.sd[0] * bA.d[0].transpose();
        for(int i = 1; i <= HorizonNum; ++i) {
            bAiGAT.d[i] = bAiG.sd[i - 1] * bA.sd[i - 1].transpose() + bAiG.d[i] * bA.d[i].transpose();
            if(i < HorizonNum) bAiGAT.sd[i] = bAiG.sd[i] * bA.d[i].transpose();
        }
        bL = bAiGAT.Cholesky();
    }
    
    BlockVectorComb LQR_Solver(BlockVectorComb f) {
        BlockVectorX y, lamb;
        BlockVectorComb iGf = biG * f;
        y = bL.LinearSolver(bb * (-1) - bA.Multiply(iGf));
        lamb = bL.TransposeLinearSolver(y);
        return (biG * bA.TransposeMultiply(lamb) + iGf) * (-1);
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
            // cout << u.transpose() << ' ' << 0.5 * u.transpose() * v + W[i].dot(u) << ' ' << 0.5 * x.transpose() * (Q[i] * x) + L[i].dot(x) << endl;
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
        while(k <= K) {
            barz = LQR_Solver(bF - (bET * (w + nu)) * (rho));
            w = G_Solver(barz, nu, rho);
            nu = nu + w - bE * barz;
            k++;
            if(k > 100 && !(k % 10)) {
                T cost = CalculateCost(barz);
                if(cost < min_cost) {
                    res = barz;
                    min_cost = cost;
                }
                cout << "Episodes: " << k << ' ' << cost << endl;
            }
        }
        cout << "Final Cost: " << min_cost << endl;
        return barE * res;
    }

    BlockVector<T, HorizonNum + 1, SizeX + SizeU> solve() {
        PreScaling();
        ADMMPrework();
        BlockVector<T, HorizonNum + 1, SizeX + SizeU> res = ADMMIteration();
        return res;
    }

};
