#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H
#pragma once
#include <stdexcept>
#include <algorithm>
#include <math.h>

class S21Matrix {
    private:
        int rows_, cols_;
        double **matrix_;
        void CreateMatrix();
        void RemoveMatrix();

    public:
        S21Matrix();
        S21Matrix(int rows, int cols);
        S21Matrix(const S21Matrix& other);
        S21Matrix(S21Matrix&& other) noexcept;
        ~S21Matrix();
        int getRow();
        int getCol();
        void setRow(int rows);
        void setCol(int cols);
        bool isValid() const;

        bool EqMatrix(const S21Matrix& other) const;
        void SumMatrix(const S21Matrix& other);
        void SubMatrix(const S21Matrix& other); 
        void MulNumber(const double num);
        void MulMatrix(const S21Matrix& other);
        S21Matrix Transpose();
        S21Matrix CalcComplements();
        double Determinant();
        S21Matrix InverseMatrix();
        S21Matrix GetMinor(const int row, const int col);
        S21Matrix operator+(const S21Matrix &other);
        S21Matrix operator-(const S21Matrix &other);
        S21Matrix operator*(const S21Matrix &other);
        S21Matrix operator*(const double num);
        bool operator==(const S21Matrix &other) const;
        S21Matrix& operator=(const S21Matrix &other);
        S21Matrix& operator+=(const S21Matrix &other);
        S21Matrix& operator-=(const S21Matrix &other);
        S21Matrix& operator*=(const S21Matrix &other);
        S21Matrix& operator*=(const double num);
        double &operator()(int row, int col);
};

#endif
