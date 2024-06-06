#pragma once

#include <Eigen/Dense>
#include <cstring>

template <typename T, int HorizonNum, int SizeNum>
class BlockVector{
public:
    typedef Eigen::Matrix<T, SizeNum, 1> UnitBlockVector;
    std::array <UnitBlockVector, HorizonNum> v;

    BlockVector<T, HorizonNum, SizeNum> operator + (BlockVector<T, HorizonNum, SizeNum> other) {
        for(int i = 0; i < HorizonNum; ++i) {
            other.v[i] += v[i];
        }
        return other;
    }

    BlockVector<T, HorizonNum, SizeNum> operator - (BlockVector<T, HorizonNum, SizeNum> other) {
        for(int i = 0; i < HorizonNum; ++i) {
            other.v[i] = v[i] - other.v[i];
        }
        return other;
    }

    BlockVector<T, HorizonNum, SizeNum> operator * (T val) {
        BlockVector<T, HorizonNum, SizeNum> res;
        for(int i = 0; i < HorizonNum; ++i) {
            res.v[i] = v[i] * val;
        }
        return res;
    }

    std::array<T, SizeNum> ChangeToArray(int i) {
        std::array<T, SizeNum> res;
        for(int j = 0; j < SizeNum; ++j) {
            res[j] = v[i](j);
        }
        return res;
    }

    void AssignFromArray(int i, std::array<T, SizeNum> val) {
        for(int j = 0; j < SizeNum; ++j) {
            v[i](j) = val[j];
        }
        return ;
    }

    void setZero() {
        for(int i = 0; i < HorizonNum; ++i) {
            v[i].setZero(SizeNum, 1);
        }
        return;
    }

    // This function is for debug
    void print(std::string name = "default name") {
        puts("==================================");
        printf("Print of Block Vector: "); std::cout << name << '\n';
        printf("Its size is %d\n, with each block of %d\n", HorizonNum, SizeNum);
        for(int i = 0; i < HorizonNum; ++i) {
            std::cout << "Block " << i << ": \n";
            std::cout << v[i] << '\n';
        }
        puts("==================================");
    }
};

template <typename T, int HorizonNum, int SizeNum>
class BlockSymmetricTridiagonalMatrix;

template <typename T, int HorizonNum, int RowNum, int ColNum>
class BlockLowerBidiagonalMatrix{
public:
    typedef Eigen::Matrix<T, RowNum, ColNum> UnitBlockMatrix;
    std::array<UnitBlockMatrix, HorizonNum> d;
    std::array<UnitBlockMatrix, HorizonNum - 1> sd;

    bool operator == (const BlockLowerBidiagonalMatrix &other) {
        for(int i = 0; i < HorizonNum; ++i) {
            if(!d[i].isApprox(other.d[i], 1e-6)) {return false;}
            if(i < HorizonNum - 1 && !sd[i].isApprox(other.sd[i], 1e-6)) {return false;}
        }
        return true;
    }

    BlockSymmetricTridiagonalMatrix<T, HorizonNum, RowNum> MultiplyWithTranspose() {
        BlockSymmetricTridiagonalMatrix<T, HorizonNum, RowNum> A;
        A.d[0] = d[0] * d[0].transpose();
        if(HorizonNum == 1) return A;
        A.sd[0] = sd[0] * d[0].transpose();
        for(int i = 1; i < HorizonNum; ++i) {
            A.d[i] = sd[i - 1] * sd[i - 1].transpose() + d[i] * d[i].transpose();
            if(i < HorizonNum - 1) A.sd[i] = sd[i] * d[i].transpose();
        }
        return A;
    }

    BlockVector<T, HorizonNum, ColNum> LinearSolver(const BlockVector<T, HorizonNum, ColNum> &b) {
        BlockVector<T, HorizonNum, ColNum> x;
        x.v[0] = d[0].inverse() * b.v[0];
        for(int i = 1; i < HorizonNum; ++i) {
            x.v[i] = d[i].inverse() * (b.v[i] - sd[i - 1] * x.v[i - 1]);
        }
        return x;
    }

    BlockVector<T, HorizonNum, ColNum> TransposeLinearSolver(const BlockVector<T, HorizonNum, ColNum> &b) {
        BlockVector<T, HorizonNum, ColNum> x;
        x.v[HorizonNum - 1] = d[HorizonNum - 1].transpose().inverse() * b.v[HorizonNum - 1];
        for(int i = HorizonNum - 2; i >= 0 ; --i) {
            x.v[i] = d[i].transpose().inverse() * (b.v[i] - sd[i].transpose() * x.v[i + 1]);
        }
        return x;
    }

    BlockVector<T, HorizonNum, RowNum> Multiply(const BlockVector<T, HorizonNum, ColNum> &x) {
        BlockVector<T, HorizonNum, RowNum> res;
        for(int i = 0; i < HorizonNum; ++i) {
            if(i == 0) {
                res.v[i] = d[i] * x.v[i];
            } else {
                res.v[i] = d[i] * x.v[i] + sd[i - 1] *  x.v[i - 1];
            }
        }
        return res;
    }

    BlockVector<T, HorizonNum, ColNum> TransposeMultiply(const BlockVector<T, HorizonNum, RowNum> &x) {
        BlockVector<T, HorizonNum, ColNum> res;
        for(int i = 0; i < HorizonNum; ++i) {
            if(i == HorizonNum - 1) {
                res.v[i] = d[i].transpose() * x.v[i];
            } else {
                res.v[i] = d[i].transpose() * x.v[i] + sd[i].transpose() *  x.v[i + 1];
            }
        }
        return res;
    }

    // This function is for debug
    void print(std::string name = "default name") {
        puts("==================================");
        printf("Print of Block Lower Bidiagonal Matrix: "); std::cout << name << '\n';
        printf("Its size is %d\n, with each block of %d * %d\n", HorizonNum, RowNum, ColNum);
        for(int i = 0; i < HorizonNum; ++i) {
            std::cout << "Diagonal Block " << i << ":\n";
            std::cout << d[i] << '\n';
            if(i < HorizonNum - 1) {
                std::cout << "Sub Diagonal Block " << i << ":\n";
                std::cout << sd[i] << '\n';
            }
        }
        puts("==================================");
    }
};
using namespace std;
template <typename T, int HorizonNum, int SizeNum>
class BlockSymmetricTridiagonalMatrix{
public:
    typedef Eigen::Matrix<T, SizeNum, SizeNum> UnitBlockMatrix;
    std::array<UnitBlockMatrix, HorizonNum> d;  //主对角线
    std::array<UnitBlockMatrix, HorizonNum - 1> sd; //主对角线下方的次对角线
    
    BlockLowerBidiagonalMatrix<T, HorizonNum, SizeNum, SizeNum> Cholesky() {
        BlockLowerBidiagonalMatrix<T, HorizonNum, SizeNum, SizeNum> L;
        
        L.d[0] = d[0].llt().matrixL();
        if(HorizonNum == 1) return L;
        L.sd[0] = sd[0] * L.d[0].transpose().inverse();
        for(int i = 1; i < HorizonNum; ++i) {
            L.d[i] = (d[i] - L.sd[i - 1] * L.sd[i - 1].transpose()).llt().matrixL();
            if(i < HorizonNum - 1) L.sd[i] = sd[i] * L.d[i].transpose().inverse();
        }
        return L;
    }
    
    BlockVector<T, HorizonNum, SizeNum> operator * (const BlockVector<T, HorizonNum, SizeNum> &x) {
        BlockVector<T, HorizonNum, SizeNum> res;
        Eigen::Matrix<T, SizeNum, 1> b;
        for(int i = 0; i < HorizonNum; ++i) {
            b = d[i] * x.v[i];
            if(i >= 1) b += sd[i - 1] * x.v[i - 1];
            if(i < HorizonNum - 1) b += sd[i].transpose() * x.v[i + 1];
            res.v[i] = b;
        }
        return res;
    }
    
    BlockSymmetricTridiagonalMatrix<T, HorizonNum, SizeNum> operator + \
                (BlockSymmetricTridiagonalMatrix<T, HorizonNum, SizeNum> other) {
        for(int i = 0; i < HorizonNum; ++i) {
            other.d[i] += d[i];
            if(i < HorizonNum - 1) other.sd[i] += sd[i];
        }
        return other;
    }

    // This function is for debug
    void print(std::string name = "default name") {
        puts("==================================");
        printf("Print of Block Symmetric Tridiagonal Matrix: "); std::cout << name << '\n';
        printf("Its size is %d\n, with each block of %d * %d\n", HorizonNum, SizeNum, SizeNum);
        for(int i = 0; i < HorizonNum; ++i) {
            std::cout << "Diagonal Block " << i << ":\n";
            std::cout << d[i] << '\n';
            if(i < HorizonNum - 1) {
                std::cout << "Sub Diagonal Block " << i << ":\n";
                std::cout << sd[i] << '\n';
            }
        }
        puts("==================================");
    }
};

template <typename T, int HorizonNum, int RowNum, int ColNum>
class BlockSymmetricDiagonalMatrix{
public:
    typedef Eigen::Matrix<T, RowNum, ColNum> UnitBlockMatrix;
    std::array<UnitBlockMatrix, HorizonNum> d;  //主对角线

    BlockVector<T, HorizonNum, RowNum> operator * (const BlockVector<T, HorizonNum, ColNum> &other) {
        BlockVector<T, HorizonNum, RowNum> res;
        for(int i = 0; i < HorizonNum; ++i) {
            res.v[i] = d[i] * other.v[i];
        }
        return res;
    }
    // This function is for debug
    void print(std::string name = "default name") {
        puts("==================================");
        printf("Print of Block Symmetric Tridiagonal Matrix: "); std::cout << name << '\n';
        printf("Its size is %d\n, with each block of %d * %d\n", HorizonNum, RowNum, ColNum);
        for(int i = 0; i < HorizonNum; ++i) {
            std::cout << "Diagonal Block " << i << ":\n";
            std::cout << d[i] << '\n';
        }
        puts("==================================");
    }
};