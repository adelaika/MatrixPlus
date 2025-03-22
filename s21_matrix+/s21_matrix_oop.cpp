#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr){}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
    if (rows < 1 || cols < 1) 
        throw std::length_error("Invalid matrix size");
    CreateMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other) : rows_(other.rows_), cols_(other.cols_){
    CreateMatrix();
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < cols_; j++){
            matrix_[i][j] = other.matrix_[i][j];
        }
    }
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_){
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
    RemoveMatrix();
}

void S21Matrix::CreateMatrix(){
    matrix_ = new double*[rows_];
    for (int i = 0; i<rows_; i++){
        matrix_[i] = new double[cols_]();
    }
}

void S21Matrix::RemoveMatrix(){
    if (matrix_){
        for (int i = 0; i < rows_; i++){
            delete[] matrix_[i];
        }
        delete[] matrix_;
        matrix_ = nullptr;
    }
}

int S21Matrix::getRow() {
    return rows_;
}

int S21Matrix::getCol() {
    return cols_;
}

// void S21Matrix::setRow(int rows){
//     if (rows <= 0) {
//         throw std::length_error("matrix rows should be > 1");
//     }
//     if (rows != rows_){
//         S21Matrix temp(rows, cols_);
//         for (int i = 0; i < std::min(rows_, rows); i++){
//             for (int j = 0; j < cols_; j++){
//                 temp.matrix_[i][j] = matrix_[i][j];
//             }
//         }
//         *this = temp;
//         rows_ = rows;
//     }
// }

void S21Matrix::setRow(int rows) {
  if (rows < 1) {
    throw std::length_error("matrix rows should be > 1");
  }
  if (rows != rows_) {
    S21Matrix tmp(rows, cols_);
    int min = std::min(rows_, rows);
    for (int i = 0; i < min; ++i) {
      for (int j = 0; j < cols_; ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }
    *this = std::move(tmp);
  }
}

void S21Matrix::setCol(int cols){
    if (cols < 1) {
        throw std::length_error("matrix cols should be > 1");
    }
    if (cols != cols_){
        S21Matrix temp(rows_, cols);
        for (int i = 0; i < rows_; i++){
            for (int j = 0; j < std::min(cols_, cols); j++){
                temp.matrix_[i][j] = matrix_[i][j];
            }
        }
        *this = std::move(temp);
        cols_ = cols;
    }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const{
    if (!this->isValid() || !other.isValid()) {
        throw std::logic_error("invalid matrix");
    }
    bool res = true;
    if (rows_ != other.rows_ || cols_ != other.cols_) res = false;
    const double epsilon = 1e-7;
    if (res){
        for (int i = 0; i < rows_ && res; i++){
            for (int j = 0; j < cols_ && res; j++){
                if (fabs(matrix_[i][j] - other.matrix_[i][j]) > epsilon) res = false;
            }
        }
    }
    
    return res;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
    if (!this->isValid() || !other.isValid()) {
        throw std::logic_error("invalid matrix");
    }
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::logic_error("Matrix dimensions do not match for addition");
    }
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < cols_; j++){
            matrix_[i][j] += other.matrix_[i][j];
        }
    }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
    if (!this->isValid() || !other.isValid()) {
        throw std::logic_error("invalid matrix");
    }
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::logic_error("Matrix dimensions do not match for addition");
    }
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < cols_; j++){
            matrix_[i][j] -= other.matrix_[i][j];
        }
    }
}

void S21Matrix::MulNumber(const double num) {
    if (!this->isValid()) {
        throw std::logic_error("invalid matrix");
    }
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < cols_; j++){
            matrix_[i][j] *= num;
        }
    }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
    if (!this->isValid() || !other.isValid()) {
        throw std::logic_error("invalid matrix");
    }
    if (cols_ != other.rows_){
        throw std::logic_error("MatrixNoneEq");
    }
    S21Matrix result(rows_, other.cols_);
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < other.cols_; j++){
            for (int k = 0; k < cols_; k++){
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
            }
        }
    }
    *this = result;
}

S21Matrix S21Matrix::Transpose() {
    if (!this->isValid()) {
        throw std::logic_error("invalid matrix");
    }
    S21Matrix result(rows_, cols_);
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < cols_; j++){
            result.matrix_[j][i] = matrix_[i][j];
        }
    }
    return result;
}

S21Matrix S21Matrix::CalcComplements() {
    if (!this->isValid()) {
        throw std::logic_error("invalid matrix");
    }
    if (rows_ != cols_) {
        throw std::logic_error(
            "Matrix dimensions must be equal");
    }
    S21Matrix calc(rows_, cols_);
    for (int i = 0; i < rows_; i++){
        for (int j = 0; j < cols_; j++){
            double min_det = 0;
            S21Matrix minor = this->GetMinor(i, j);
            min_det = minor.Determinant();
            calc(i, j) = min_det * pow(-1, i + j);
        }
    }
    return calc;
}

double S21Matrix::Determinant() {
    if (!this->isValid()) {
        throw std::logic_error("invalid matrix");
    }
    if (rows_ != cols_){
        if (rows_ != cols_) {
            throw std::logic_error(
                "To get Determinant matrix dimensions must be equal");
        }
    }
    double det = 0;
    if (rows_ == 1){
        det = matrix_[0][0];
    } else if (rows_ == 2){
        det = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    } else {
        int sign = 1;
        for (int i = 0; i < cols_; i++){
            S21Matrix minor = this->GetMinor(0, i);
            double det_temp = minor.Determinant();
            det += sign * matrix_[0][i] * det_temp;
            sign *= -1;
        }
    }
    return det;
}

S21Matrix S21Matrix::InverseMatrix(){
    if (!this->isValid()) {
        throw std::logic_error("invalid matrix");
    }
    if (rows_ != cols_) {
        throw std::logic_error("Matrix dimensions must be equal");
    }
    double det = Determinant();
    if (det == 0) {
        throw std::logic_error("Determinant must be non-zero");
    }

    S21Matrix calc = CalcComplements();
    S21Matrix result = calc.Transpose();
    result.MulNumber(1.0/det);
    return result;
}

S21Matrix S21Matrix::GetMinor(const int row, const int col){
    S21Matrix minor(rows_ - 1, cols_ - 1);
    int m_r = 0;
    for (int i = 0; i < rows_; i++){
        if (i != row){
            int m_c = 0;
            for (int j = 0; j < cols_; j++){
                if (j != col){
                    minor(m_r, m_c) = matrix_[i][j];
                    m_c++;
                }
            }
            m_r++;
        }
    }
    return minor;
}

bool S21Matrix::isValid() const{
    return matrix_ != nullptr && rows_ > 0 && cols_ > 0;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other){
    S21Matrix result(*this);
    result.SumMatrix(other);
    return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other){
    S21Matrix result(*this);
    result.SubMatrix(other);
    return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other){
    S21Matrix result(*this);
    result.MulMatrix(other);
    return result;
}

S21Matrix S21Matrix::operator*(const double num){
    S21Matrix result(*this);
    result.MulNumber(num);
    return result;
}

bool S21Matrix::operator==(const S21Matrix &other) const{
    return EqMatrix(other);
}

S21Matrix& S21Matrix::operator=(const S21Matrix &other){
    if (this != &other){
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::logic_error("Matrix dimensions must match for assignment");
        }
        S21Matrix temp(other);
        std::swap(rows_, temp.rows_);
        std::swap(cols_, temp.cols_);
        std::swap(matrix_, temp.matrix_);
    }
    return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix &other){
    SumMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix &other){
    SubMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix &other){
    MulMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
    MulNumber(num);
    return *this;
}

double& S21Matrix::operator()(int row, int col){
    if (row > rows_ || col > cols_ || row < 0 || col < 0){
        throw std::out_of_range("Index out of range");
    }
    return matrix_[row][col];
}

// int main(){
//     S21Matrix test(2, 2);
//     test.setRow(5);
//     printf("%d\n", test.getRow());
// }
