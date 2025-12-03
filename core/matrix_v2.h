#ifndef GONS_CORE_MATRIX_H_
#define GONS_CORE_MATRIX_H_

#include "config.h"
#include "utilites.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace gons {
using namespace utilites::LOG_MSG;

template <typename T, GONS_UINT R> class Vector;

template <typename T, GONS_UINT R, GONS_UINT C> class Matrix {
private:
  static_assert(R > 0 && C > 0, "Matrix dimensions must be positive");
  static_assert(std::is_arithmetic<T>::value, "Matrix type must be arithmetic");

  T data_[R][C] = {{0}};

public:
  // =============== 构造和赋值 ===============

  // 默认构造函数 - 零初始化
  constexpr Matrix() = default;

  // 拷贝构造和赋值
  Matrix(const Matrix &other) = default;
  Matrix &operator=(const Matrix &other) = default;

  // 移动构造和赋值
  Matrix(Matrix &&other) = default;
  Matrix &operator=(Matrix &&other) = default;

  // 初始化列表构造 - 行优先
  Matrix(const std::initializer_list<std::initializer_list<T>> &init_list) {
    GONS_UINT i = 0;
    for (const auto &row : init_list) {
      CHECK(i >= R, "Initializer list has more rows than matrix");
      GONS_UINT j = 0;
      for (const auto &elem : row) {
        CHECK(j >= C, "Initializer list has more columns than matrix");
        data_[i][j] = elem;
        ++j;
      }
      // 填充剩余列为0
      for (; j < C; ++j) {
        data_[i][j] = T(0);
      }
      ++i;
    }
    // 填充剩余行为0
    for (; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        data_[i][j] = T(0);
      }
    }
  }

  // 扁平化初始化列表构造
  explicit Matrix(const std::initializer_list<T> &init_list) {
    auto it = init_list.begin();
    for (GONS_UINT i = 0; i < R && it != init_list.end(); ++i) {
      for (GONS_UINT j = 0; j < C && it != init_list.end(); ++j, ++it) {
        data_[i][j] = *it;
      }
    }
  }

  // =============== 数据访问 ===============

  // 获取底层数据指针
  constexpr T *data() noexcept { return &data_[0][0]; }
  constexpr const T *data() const noexcept { return &data_[0][0]; }

  // 获取行列数
  static constexpr GONS_UINT rows() noexcept { return R; }
  static constexpr GONS_UINT cols() noexcept { return C; }
  static constexpr std::pair<GONS_UINT, GONS_UINT> shape() noexcept {
    return {R, C};
  }

  // 元素访问 - 带边界检查
  T &operator()(GONS_UINT row, GONS_UINT col) {
    CHECK(row >= R || col >= C, "Matrix index out of bounds");
    return data_[row][col];
  }

  const T &operator()(GONS_UINT row, GONS_UINT col) const {
    CHECK(row >= R || col >= C, "Matrix index out of bounds");
    return data_[row][col];
  }

  // =============== 算术运算 ===============

  // 矩阵加法
  Matrix operator+(const Matrix &other) const {
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result.data_[i][j] = data_[i][j] + other.data_[i][j];
      }
    }
    return result;
  }

  Matrix &operator+=(const Matrix &other) {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        data_[i][j] += other.data_[i][j];
      }
    }
    return *this;
  }

  // 矩阵减法
  Matrix operator-(const Matrix &other) const {
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result.data_[i][j] = data_[i][j] - other.data_[i][j];
      }
    }
    return result;
  }

  Matrix &operator-=(const Matrix &other) {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        data_[i][j] -= other.data_[i][j];
      }
    }
    return *this;
  }

  // 矩阵乘法
  template <GONS_UINT C2>
  Matrix<T, R, C2> operator*(const Matrix<T, C, C2> &other) const {
    Matrix<T, R, C2> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C2; ++j) {
        T sum = 0;
        for (GONS_UINT k = 0; k < C; ++k) {
          sum += data_[i][k] * other(k, j);
        }
        result(i, j) = sum;
      }
    }
    return result;
  }

  // 标量乘法
  Matrix operator*(T scalar) const {
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result.data_[i][j] = data_[i][j] * scalar;
      }
    }
    return result;
  }

  Matrix &operator*=(T scalar) {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        data_[i][j] *= scalar;
      }
    }
    return *this;
  }

  friend Matrix operator*(T scalar, const Matrix &matrix) {
    return matrix * scalar;
  }

  // 标量除法
  Matrix operator/(T scalar) const {
    CHECK(std::abs(scalar) < GONS_FLT_EPSILON, "Division by zero");
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result.data_[i][j] = data_[i][j] / scalar;
      }
    }
    return result;
  }

  Matrix &operator/=(T scalar) {
    CHECK(std::abs(scalar) < GONS_FLT_EPSILON, "Division by zero");
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        data_[i][j] /= scalar;
      }
    }
    return *this;
  }

  // =============== 矩阵操作 ===============

  // 转置
  Matrix<T, C, R> transpose() const {
    Matrix<T, C, R> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(j, i) = data_[i][j];
      }
    }
    return result;
  }

  // 单位矩阵
  static Matrix identity() {
    static_assert(R == C, "Identity matrix must be square");
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      result(i, i) = T(1);
    }
    return result;
  }

  // 零矩阵
  static Matrix zeros() { return Matrix(); }

  // 增广矩阵
  Matrix<T, R, 2 * C> augment() const {
    Matrix<T, R, 2 * C> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(i, j) = data_[i][j];
      }
      result(i, C + i) = T(1);
    }
    return result;
  }

  // 矩阵求逆 - 使用高斯-约当消元法
  Matrix inverse() const {
    static_assert(R == C, "Inverse requires square matrix");
    CHECK(R == 0, "Matrix is empty");

    // 使用增广矩阵方法
    auto augmented = this->augment();

    for (GONS_UINT i = 0; i < R; ++i) {
      // 寻找主元
      GONS_UINT pivot_row = i;
      T max_val = std::abs(augmented(i, i));

      for (GONS_UINT k = i + 1; k < R; ++k) {
        T current_val = std::abs(augmented(k, i));
        if (current_val > max_val) {
          max_val = current_val;
          pivot_row = k;
        }
      }

      CHECK(max_val < GONS_FLT_EPSILON,
            "Matrix is singular or nearly singular");

      // 交换行
      if (pivot_row != i) {
        for (GONS_UINT j = 0; j < 2 * C; ++j) {
          std::swap(augmented(i, j), augmented(pivot_row, j));
        }
      }

      // 主元归一化
      T pivot = augmented(i, i);
      for (GONS_UINT j = 0; j < 2 * C; ++j) {
        augmented(i, j) /= pivot;
      }

      // 消元
      for (GONS_UINT k = 0; k < R; ++k) {
        if (k != i) {
          T factor = augmented(k, i);
          for (GONS_UINT j = 0; j < 2 * C; ++j) {
            augmented(k, j) -= factor * augmented(i, j);
          }
        }
      }
    }

    // 提取逆矩阵
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(i, j) = augmented(i, C + j);
      }
    }
    return result;
  }

  // =============== 范数计算 ===============

  // Frobenius 范数
  T frobeniusNorm() const {
    T sum = 0;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        sum += data_[i][j] * data_[i][j];
      }
    }
    return std::sqrt(sum);
  }

  // 2-范数（对于矩阵来说通常是谱范数，这里简化为Frobenius范数）
  T norm() const { return frobeniusNorm(); }

  // 最大绝对值
  T maxAbs() const {
    T max_val = std::abs(data_[0][0]);
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        T current = std::abs(data_[i][j]);
        if (current > max_val) {
          max_val = current;
        }
      }
    }
    return max_val;
  }

  // =============== 行列操作 ===============

  // 获取行
  Vector<T, C> getRow(GONS_UINT row) const {
    CHECK(row >= R, "Row index out of range");
    Vector<T, C> result;
    for (GONS_UINT j = 0; j < C; ++j) {
      result(j) = data_[row][j];
    }
    return result;
  }

  // 获取列
  Vector<T, R> getCol(GONS_UINT col) const {
    CHECK(col >= C, "Column index out of range");
    Vector<T, R> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      result(i) = data_[i][col];
    }
    return result;
  }

  // 设置行
  void setRow(GONS_UINT row, const Vector<T, C> &vec) {
    CHECK(row >= R, "Row index out of range");
    for (GONS_UINT j = 0; j < C; ++j) {
      data_[row][j] = vec(j);
    }
  }

  // 设置列
  void setCol(GONS_UINT col, const Vector<T, R> &vec) {
    CHECK(col >= C, "Column index out of range");
    for (GONS_UINT i = 0; i < R; ++i) {
      data_[i][col] = vec(i);
    }
  }

  // =============== 矩阵-向量运算 ===============

  Vector<T, R> operator*(const Vector<T, C> &vec) const {
    Vector<T, R> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      T sum = 0;
      for (GONS_UINT j = 0; j < C; ++j) {
        sum += data_[i][j] * vec(j);
      }
      result(i) = sum;
    }
    return result;
  }

  friend Vector<T, C> operator*(const Vector<T, R> &vec, const Matrix &mat) {
    Vector<T, C> result;
    for (GONS_UINT j = 0; j < C; ++j) {
      T sum = 0;
      for (GONS_UINT i = 0; i < R; ++i) {
        sum += vec(i) * mat(i, j);
      }
      result(j) = sum;
    }
    return result;
  }

  // =============== 输出和调试 ===============

#if (ENABLE_PRINT_MATRIX == 1)
  void print(const std::string &prefix = "") const {
    if (!prefix.empty()) {
      std::cout << prefix << "\n";
    }
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        std::cout << data_[i][j] << " ";
      }
      std::cout << "\n";
    }
  }
#endif

#if (EMBEDDED_MODE == 0)
  friend std::ostream &operator<<(std::ostream &os, const Matrix &mat) {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        os << mat.data_[i][j];
        if (j < C - 1)
          os << ' ';
      }
      if (i < R - 1)
        os << '\n';
    }
    return os;
  }
#endif

  // =============== 比较操作 ===============

  bool operator==(const Matrix &other) const {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        if (std::abs(data_[i][j] - other.data_[i][j]) > GONS_FLT_EPSILON) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const Matrix &other) const { return !(*this == other); }

}; // class Matrix

} // namespace gons

#endif // GONS_CORE_MATRIX_H_