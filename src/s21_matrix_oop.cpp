#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() {
  rows_ = cols_ = 0;
  matrix_ = nullptr;
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::invalid_argument("Invalid matrix size 1");
  } else {
    rows_ = rows, cols_ = cols;
    matrix_ = new double* [rows_] {};
    for (int i{}; i < rows_; i++) {
      matrix_[i] = new double[cols_]{};
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  this->rows_ = other.rows_, this->cols_ = other.cols_;

  this->matrix_ = new double* [rows_] {};
  for (int i{}; i < rows_; i++) this->matrix_[i] = new double[cols_]{};

  for (int i{}; i < rows_; i++)
    for (int j{}; j < cols_; j++) this->matrix_[i][j] = other.matrix_[i][j];
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  this->rows_ = other.rows_, this->cols_ = other.cols_,
  this->matrix_ = other.matrix_;
  other.matrix_ = nullptr, other.rows_ = other.cols_ = 0;
}

S21Matrix::~S21Matrix() noexcept {
  if (matrix_) {
    for (int i{}; i < rows_; ++i) {
      if (matrix_[i]) {
        delete[] matrix_[i];
        matrix_[i] = nullptr;
      }
    }
    delete[] matrix_;
    matrix_ = nullptr;
    rows_ = cols_ = 0;
  }
}

int S21Matrix::GetRows() { return rows_; }

int S21Matrix::GetCols() { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (!matrix_ || rows <= 0) {
    throw std::invalid_argument("Invalid matrix size 2");
  }
  int i_rows = (rows_ < rows) ? rows_ : rows;
  S21Matrix res(rows, cols_);
  for (int i{}; i < i_rows; i++)
    for (int j{}; j < cols_; j++) res.matrix_[i][j] = matrix_[i][j];
  *this = res;
}

void S21Matrix::SetCols(int cols) {
  if (!matrix_ || cols <= 0)
    throw std::invalid_argument("Invalid matrix size 3");

  int i_cols = (cols_ < cols) ? cols_ : cols;
  S21Matrix res(rows_, cols);
  for (int i{}; i < rows_; i++)
    for (int j{}; j < i_cols; j++) res.matrix_[i][j] = matrix_[i][j];
  *this = res;
}

bool S21Matrix::CheckRowsAndCols(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_) return false;
  return true;
}

bool S21Matrix::CheckRowsAndCols() {
  if (this->rows_ <= 0 || this->cols_ <= 0) return false;
  return true;
}

bool S21Matrix::CheckSquareMatrix() {
  if (this->rows_ != this->cols_) return false;
  return true;
}

bool S21Matrix::CheckMatrix(const S21Matrix& other) {
  if (this->rows_ <= 0 || other.rows_ <= 0 || this->cols_ <= 0 ||
      other.cols_ <= 0 || !this->matrix_ || !other.matrix_)
    return false;
  return true;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (!this->CheckMatrix(other))
    throw std::invalid_argument("Invalid matrix size");
  if (!this->CheckRowsAndCols(other)) return false;
  for (int i{}; i < rows_; i++)
    for (int j{}; j < cols_; j++)
      if (this->matrix_[i][j] != other.matrix_[i][j]) return false;
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!this->CheckRowsAndCols(other))
    throw std::invalid_argument("Exceptional situation: Invalid matrix size");
  if (!this->CheckMatrix(other))
    throw std::invalid_argument("Invalid matrix size");
  for (int i{}; i < rows_; i++)
    for (int j{}; j < cols_; j++) this->matrix_[i][j] += other.matrix_[i][j];
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!this->CheckRowsAndCols(other))
    throw std::invalid_argument("Exceptional situation: Invalid matrix size");
  if (!this->CheckMatrix(other))
    throw std::invalid_argument("Invalid matrix size");
  for (int i{}; i < rows_; i++)
    for (int j{}; j < cols_; j++) this->matrix_[i][j] -= other.matrix_[i][j];
}

void S21Matrix::MulNumber(const double num) {
  if (!this->CheckRowsAndCols() || !this->matrix_)
    throw std::invalid_argument("Invalid matrix size");
  for (int i{}; i < rows_; i++)
    for (int j{}; j < cols_; j++) matrix_[i][j] *= num;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->rows_ != other.cols_)
    throw std::invalid_argument("Exceptional situation: Invalid matrix size");
  if (!this->CheckMatrix(other))
    throw std::invalid_argument("Invalid matrix size");
  S21Matrix res(this->rows_, other.cols_);
  for (int i{}; i < this->rows_; i++)
    for (int j{}; j < other.cols_; j++)
      for (int z{}; z < other.rows_; z++)
        res.matrix_[i][j] += this->matrix_[i][z] * other.matrix_[z][j];
  (*this) = res;
}

S21Matrix S21Matrix::Transpose() {
  if (!this->CheckRowsAndCols() || !matrix_) {
    throw std::invalid_argument("Invalid matrix");
  }
  S21Matrix res(this->cols_, this->rows_);
  for (int i{}; i < this->rows_; i++)
    for (int j{}; j < this->cols_; j++) res.matrix_[j][i] = this->matrix_[i][j];
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (!this->CheckSquareMatrix())
    throw std::invalid_argument("Exceptional situation: Invalid matrix");
  if (!this->CheckRowsAndCols() || !matrix_)
    throw std::invalid_argument("Invalid matrix");

  S21Matrix res(rows_, cols_);
  if (rows_ == 1) {
    res.matrix_[0][0] = this->matrix_[0][0];
    return res;
  }
  for (int i{}; i < rows_; i++)
    for (int j{}; j < cols_; j++) {
      double det = 0;
      S21Matrix minor(rows_ - 1, cols_ - 1);
      Minor(this, &minor, i, j);
      det = minor.Determinant();
      res.matrix_[i][j] = det * pow(-1, i + j);  // corrected formula
    }
  return res;
}

double S21Matrix::Determinant() {
  if (!this->CheckSquareMatrix()) {
    throw std::invalid_argument("Exceptional situation: Invalid matrix");
  }
  if (!this->CheckRowsAndCols() || !this->matrix_)
    throw std::invalid_argument("Invalid matrix");
  double res = 0;
  if (cols_ == 1)
    return matrix_[0][0];
  else if (cols_ == 2)
    return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  else {
    bool isZero = true;
    for (int i{}; i < rows_; i++) {
      for (int j{}; j < cols_; j++) {
        if (matrix_[i][j] != 0) {
          isZero = false;
          break;
        }
      }
      if (!isZero) break;
    }
    if (isZero) return 0;
    for (int i{}; i < rows_; i++) {
      S21Matrix minor(rows_ - 1, cols_ - 1);
      double det;
      Minor(this, &minor, 0, i);
      det = minor.Determinant();
      res += pow(-1, i) * matrix_[0][i] * det;
    }
  }

  return res;
}

void S21Matrix::Minor(const S21Matrix* A, S21Matrix* minor, int rows, int col) {
  if (A->rows_ - 1 <= 0 || A->cols_ - 1 <= 0 || !this->CheckMatrix(*A) ||
      !this->CheckMatrix(*minor)) {
    throw std::invalid_argument("Invalid matrix size");
  }
  if (A->rows_ == 1 && A->cols_ == 1) {
    minor->rows_ = 1;
    minor->cols_ = 1;
    minor->matrix_[0][0] = 1;
    return;
  }
  minor->SetRows(A->rows_ - 1);
  minor->SetCols(A->cols_ - 1);
  int x = 0, flag = 0;
  for (int i{}; i < A->rows_; i++) {
    if (flag && i != rows && x < minor->GetRows()) {
      x++;
    }
    int y = 0;
    for (int j{}; j < A->cols_; j++)
      if (i != rows && j != col && y < minor->GetCols()) {
        minor->matrix_[x][y] = A->matrix_[i][j];
        flag = 1;
        y++;
      }
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  if (!this->CheckRowsAndCols() || !this->matrix_ || !this->CheckSquareMatrix())
    throw std::invalid_argument("Invalid matrix size");

  double inv_det = this->Determinant();
  if (inv_det == 0)
    throw std::invalid_argument("Exceptional situation: Invalid determinant");
  else
    inv_det = 1 / inv_det;

  S21Matrix res = this->CalcComplements();
  if (this->rows_ == 1) {
    res.matrix_[0][0] = inv_det;
    return res;
  }
  S21Matrix res2 = res.Transpose();
  res2.MulNumber(inv_det);
  return res2;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) return *this;

  for (int i{}; i < rows_; ++i) delete[] matrix_[i];
  delete[] matrix_;

  this->rows_ = other.rows_;
  this->cols_ = other.cols_;

  matrix_ = new double*[rows_]();
  for (int i{}; i < rows_; ++i) matrix_[i] = new double[cols_]();

  for (int i{}; i < rows_; ++i)
    for (int j{}; j < cols_; ++j) matrix_[i][j] = other.matrix_[i][j];

  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this == &other) return *this;

  for (int i{}; i < rows_; ++i) delete[] matrix_[i];
  delete[] matrix_;

  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;

  other.matrix_ = nullptr;
  other.rows_ = other.cols_ = 0;

  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
    throw std::out_of_range("Matrix index out of range");
  }
  return matrix_[i][j];
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

S21Matrix operator*(double num, S21Matrix& other) {
  S21Matrix res(other);
  res.MulNumber(num);
  return res;
}

S21Matrix operator*(S21Matrix& other, double num) { return num * other; }

bool S21Matrix::operator==(const S21Matrix& other) {
  bool eq;
  eq = this->EqMatrix(other);
  return eq;
}
S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}
S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}
S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}
S21Matrix& S21Matrix::operator*=(const double& num) {
  this->MulNumber(num);
  return *this;
}
