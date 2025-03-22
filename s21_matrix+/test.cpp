#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <iostream>

#include "s21_matrix_oop.h"

TEST(CreateMatrix1, StandartConstructor) {
  S21Matrix matrix_a;
  ASSERT_TRUE(matrix_a.getCol() == 0);
  ASSERT_TRUE(matrix_a.getRow() == 0);
}

TEST(CreateMatrix2, ParamConstr) { ASSERT_NO_THROW(S21Matrix matrix(2, 2)); }

TEST(CreateMatrix3, ParamConstWithWrongSizes) {
  EXPECT_THROW(S21Matrix matrix(-1, 2), std::length_error);
}

TEST(CreateMatrix4, MoveConstr) {
  S21Matrix matrix_a(2, 2);
  ASSERT_NO_THROW(S21Matrix matrix_b = std::move(matrix_a));
}

TEST(CreateMatrix5, CopyConstr) {
  S21Matrix matrix_a(2, 2);
  ASSERT_NO_THROW(S21Matrix matrix_b(matrix_a));
}

TEST(OperatorParentheses3, MatrixsetterWrongIndexes) {
  S21Matrix matrix_a(2, 2);
  EXPECT_THROW(matrix_a(3, 0) = 1, std::logic_error);
}

TEST(eqMatrix3, NonEqMatrix) {
  S21Matrix matrix_a(1, 3);
  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;
  matrix_a(0, 2) = 3;
  S21Matrix matrix_b(2, 2);
  matrix_b(0, 0) = 1;
  matrix_b(0, 1) = 2;
  matrix_b(1, 0) = 3;
  matrix_b(1, 1) = 4;
  ASSERT_FALSE(matrix_a == matrix_b);
}

TEST(SumMatrix2, SumNonEqMatrix) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;

  matrix_b(0, 0) = 1;
  matrix_b(0, 1) = 2;
  matrix_b(1, 0) = 3;
  matrix_b(1, 1) = 4;

  EXPECT_THROW(matrix_a.SumMatrix(matrix_b), std::logic_error);
}

TEST(SubMatrix2, SubNonEqMatrix) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a.SubMatrix(matrix_b), std::logic_error);
}

TEST(MulMatrix2, MulNonEqMatrix) {
  S21Matrix matrix_a(2, 1);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(1, 0) = -6.6;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a.MulMatrix(matrix_b), std::logic_error);
}

TEST(OperatorSumMatrix2, OperatorSumWithNoEqSizesMatrix) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;

  matrix_b(0, 0) = 1;
  matrix_b(0, 1) = 2;
  matrix_b(1, 0) = 3;
  matrix_b(1, 1) = 4;

  EXPECT_THROW(matrix_a += matrix_b, std::logic_error);
}

TEST(OperatorSubMatrix2, False) {
  S21Matrix matrix_a(1, 2);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(0, 1) = 2;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a -= matrix_b, std::logic_error);
}

TEST(OperatorMulMatrix2, False) {
  S21Matrix matrix_a(2, 1);
  S21Matrix matrix_b(2, 2);

  matrix_a(0, 0) = 3;
  matrix_a(1, 0) = -6.6;

  matrix_b(0, 0) = -7;
  matrix_b(0, 1) = 0;
  matrix_b(1, 0) = -3.5;
  matrix_b(1, 1) = 2;

  EXPECT_THROW(matrix_a *= matrix_b, std::logic_error);
}

TEST(OperatorEqal, Invalid) {
  S21Matrix matrix_a(2, 1);
  S21Matrix matrix_b;

  EXPECT_THROW(matrix_a = matrix_b, std::logic_error);
}

TEST(invalidmatrix, False) {
  S21Matrix matrix_a(2, 1);
  S21Matrix matrix_b;

  matrix_a(0, 0) = 3;
  matrix_a(1, 0) = -6.6;

  EXPECT_THROW(matrix_b *= matrix_a, std::logic_error);
  EXPECT_THROW(matrix_b.EqMatrix(matrix_a), std::logic_error);
  EXPECT_THROW(matrix_b.Transpose(), std::logic_error);
  EXPECT_THROW(matrix_b.InverseMatrix(), std::logic_error);
  EXPECT_THROW(matrix_b.Determinant(), std::logic_error);
  EXPECT_THROW(matrix_b.CalcComplements(), std::logic_error);
  EXPECT_THROW(matrix_b *= 2, std::logic_error);
  EXPECT_THROW(matrix_b += matrix_a, std::logic_error);
  EXPECT_THROW(matrix_b -= matrix_a, std::logic_error);
}

TEST(Inverse1, False) {
  S21Matrix matrix_a(4, 3);

  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;
  matrix_a(0, 2) = 3;
  matrix_a(1, 0) = 4;
  matrix_a(1, 1) = 5;
  matrix_a(1, 2) = 6;
  matrix_a(2, 0) = 7;
  matrix_a(2, 1) = 8;
  matrix_a(2, 2) = 9;
  matrix_a(3, 0) = 10;
  matrix_a(3, 1) = 11;
  matrix_a(3, 2) = 12;

  EXPECT_THROW(matrix_a.InverseMatrix(), std::logic_error);
}

TEST(Det52, False) {
  S21Matrix matrix_a(4, 3);

  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;
  matrix_a(0, 2) = 3;
  matrix_a(1, 0) = 4;
  matrix_a(1, 1) = 5;
  matrix_a(1, 2) = 6;
  matrix_a(2, 0) = 7;
  matrix_a(2, 1) = 8;
  matrix_a(2, 2) = 9;
  matrix_a(3, 0) = 10;
  matrix_a(3, 1) = 11;
  matrix_a(3, 2) = 12;

  EXPECT_THROW(matrix_a.Determinant(), std::logic_error);
}

TEST(calc52, False) {
  S21Matrix matrix_a(4, 3);

  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;
  matrix_a(0, 2) = 3;
  matrix_a(1, 0) = 4;
  matrix_a(1, 1) = 5;
  matrix_a(1, 2) = 6;
  matrix_a(2, 0) = 7;
  matrix_a(2, 1) = 8;
  matrix_a(2, 2) = 9;
  matrix_a(3, 0) = 10;
  matrix_a(3, 1) = 11;
  matrix_a(3, 2) = 12;

  EXPECT_THROW(matrix_a.CalcComplements(), std::logic_error);
}

TEST(Inverse3, False) {
  S21Matrix matrix_a(3, 3);
  matrix_a(0, 0) = 1;
  matrix_a(0, 1) = 2;
  matrix_a(0, 2) = 3;
  matrix_a(1, 0) = 4;
  matrix_a(1, 1) = 5;
  matrix_a(1, 2) = 6;
  matrix_a(2, 0) = 7;
  matrix_a(2, 1) = 8;
  matrix_a(2, 2) = 9;

  EXPECT_THROW(matrix_a.InverseMatrix(), std::logic_error);
}

TEST(Constructor, Default) {
  auto test = S21Matrix();
  EXPECT_EQ(test.getRow(), 0);
  EXPECT_EQ(test.getCol(), 0);
}

TEST(Constructor, By2Args) {
  auto test = S21Matrix(3, 3);
  EXPECT_EQ(test.getRow(), 3);
  EXPECT_EQ(test.getCol(), 3);
}

TEST(Constructor, Copy) {
  auto test = S21Matrix(3, 3);
  test(0, 0) = 1;
  auto temp = S21Matrix(test);
  EXPECT_EQ(test, temp);
}

TEST(Constructor, Move) {
  auto test1 = S21Matrix(3, 3);
  test1(1, 1) = 5;
  auto test2 = S21Matrix(test1);
  auto test3 = S21Matrix(std::move(test1));
  EXPECT_EQ(test2, test3);
}

TEST(SimpleMathOperations, EqMatrix) {
  auto test1 = S21Matrix(3, 3);
  test1(1, 2) = 5;
  auto test2 = test1;
  EXPECT_TRUE(test1.EqMatrix(test2));
  test2(1, 1) = 1;
  EXPECT_FALSE(test1.EqMatrix(test2));
  auto test3 = S21Matrix(5, 5);
  EXPECT_FALSE(test1.EqMatrix(test3));
}

TEST(SimpleMathOperations, SumMatrix) {
  auto test1 = S21Matrix(1, 1);
  test1(0, 0) = 1;
  auto test2 = test1;
  test1.SumMatrix(test2);
  EXPECT_EQ(test1(0, 0), 2);
}

TEST(SimpleMathOperations, SubMatrix) {
  auto test1 = S21Matrix(1, 1);
  test1(0, 0) = 1;
  auto test2 = test1;
  test1.SubMatrix(test2);
  EXPECT_EQ(test1(0, 0), 0);
}

TEST(SimpleMathOperations, SubMatrix1) {
  auto test1 = S21Matrix(1, 1);
  test1(0, 0) = 1;
  auto test2 = test1;
  test1 = test1 - test2;
  EXPECT_EQ(test1(0, 0), 0);
}

TEST(SimpleMathOperations, SumMatrix1) {
  auto test1 = S21Matrix(2, 2);
  auto test_res = S21Matrix(2, 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      test1(i, j) = 1;
      test_res(i, j) = 2;
    }
  }
  auto test2 = test1;
  test1 = (test1 + test2);
  EXPECT_EQ(test1, test_res);
}

TEST(SimpleMathOperations, MulMatrix2) {
  auto test1 = S21Matrix(1, 1);
  test1(0, 0) = 1;
  auto test2 = test1;
  test1 = (test1 * 0);
  EXPECT_EQ(test1(0, 0), 0);
}

TEST(SimpleMathOperations, MulNumber) {
  auto test1 = S21Matrix(1, 1);
  test1(0, 0) = 5;
  test1.MulNumber(3);
  EXPECT_EQ(test1(0, 0), 15);
}

TEST(SimpleMathOperations, MulMatrix) {
  auto test1 = S21Matrix(3, 3);
  for (int i = 0; i < test1.getRow(); i++) {
    for (int j = 0; j < test1.getCol(); j++) {
      test1(i, j) = test1.getRow() * i + j + 1;
    }
  }
  auto test2 = test1;
  test1.MulMatrix(test2);
  test2(0, 0) = 30;
  test2(0, 1) = 36;
  test2(0, 2) = 42;
  test2(1, 0) = 66;
  test2(1, 1) = 81;
  test2(1, 2) = 96;
  test2(2, 0) = 102;
  test2(2, 1) = 126;
  test2(2, 2) = 150;
  EXPECT_EQ(test1, test2);
}

TEST(SimpleMathOperations, MulMatrix1) {
  auto test1 = S21Matrix(3, 3);
  for (int i = 0; i < test1.getRow(); i++) {
    for (int j = 0; j < test1.getCol(); j++) {
      test1(i, j) = test1.getRow() * i + j + 1;
    }
  }
  auto test2 = test1;
  test1 *= test2;

  test2(0, 0) = 30;
  test2(0, 1) = 36;
  test2(0, 2) = 42;
  test2(1, 0) = 66;
  test2(1, 1) = 81;
  test2(1, 2) = 96;
  test2(2, 0) = 102;
  test2(2, 1) = 126;
  test2(2, 2) = 150;
  EXPECT_EQ(test1, test2);
}

TEST(LinearOperations, Transpose) {
  auto test1 = S21Matrix(3, 3);
  for (int i = 0; i < test1.getRow(); i++) {
    for (int j = 0; j < test1.getCol(); j++) {
      test1(i, j) = test1.getRow() * i + j + 1;
    }
  }
  test1 = test1.Transpose();
  auto test2 = S21Matrix(3, 3);
  test2(0, 0) = 1;
  test2(0, 1) = 4;
  test2(0, 2) = 7;
  test2(1, 0) = 2;
  test2(1, 1) = 5;
  test2(1, 2) = 8;
  test2(2, 0) = 3;
  test2(2, 1) = 6;
  test2(2, 2) = 9;
  EXPECT_EQ(test1, test2);
}

TEST(LinearOperations, CalcComplements) {
  auto test1 = S21Matrix(3, 3);
  test1(0, 0) = 1;
  test1(0, 1) = 2;
  test1(0, 2) = 3;
  test1(1, 0) = 0;
  test1(1, 1) = 4;
  test1(1, 2) = 2;
  test1(2, 0) = 5;
  test1(2, 1) = 2;
  test1(2, 2) = 1;
  test1 = test1.CalcComplements();
  auto test2 = S21Matrix(3, 3);
  test2(0, 0) = 0;
  test2(0, 1) = 10;
  test2(0, 2) = -20;
  test2(1, 0) = 4;
  test2(1, 1) = -14;
  test2(1, 2) = 8;
  test2(2, 0) = -8;
  test2(2, 1) = -2;
  test2(2, 2) = 4;
  EXPECT_EQ(test1, test2);
}

TEST(LinearOperations, CalcComplements1) {
  auto test1 = S21Matrix(2, 3);
  test1(0, 0) = 1;
  test1(0, 1) = 2;
  test1(0, 2) = 3;
  test1(1, 0) = 0;
  test1(1, 1) = 4;
  test1(1, 2) = 2;
  EXPECT_THROW(test1.CalcComplements(), std::logic_error);
}

TEST(LinearOperations, Determinant) {
  auto test1 = S21Matrix(3, 3);
  test1(0, 0) = 1;
  test1(0, 1) = 2;
  test1(0, 2) = 3;
  test1(1, 0) = 0;
  test1(1, 1) = 4;
  test1(1, 2) = 2;
  test1(2, 0) = 5;
  test1(2, 1) = 2;
  test1(2, 2) = 1;
  EXPECT_DOUBLE_EQ(test1.Determinant(), -40);
}

TEST(LinearOperations, Determinant1) {
  auto test1 = S21Matrix(1, 1);
  test1(0, 0) = -40;
  EXPECT_DOUBLE_EQ(test1.Determinant(), -40);
}

TEST(LinearOperations, InverseMatrix) {
  auto test1 = S21Matrix(3, 3);
  test1(0, 0) = 2;
  test1(0, 1) = 5;
  test1(0, 2) = 7;
  test1(1, 0) = 6;
  test1(1, 1) = 3;
  test1(1, 2) = 4;
  test1(2, 0) = 5;
  test1(2, 1) = -2;
  test1(2, 2) = -3;
  test1 = test1.InverseMatrix();
  auto test2 = S21Matrix(3, 3);
  test2(0, 0) = 1;
  test2(0, 1) = -1;
  test2(0, 2) = 1;
  test2(1, 0) = -38;
  test2(1, 1) = 41;
  test2(1, 2) = -34;
  test2(2, 0) = 27;
  test2(2, 1) = -29;
  test2(2, 2) = 24;
  EXPECT_EQ(test1, test2);
}

TEST(Technical, getRow) {
  auto test1 = S21Matrix();
  EXPECT_EQ(test1.getRow(), 0);
  auto test2 = S21Matrix(2, 2);
  EXPECT_EQ(test2.getRow(), 2);
}

TEST(Technical, getCol) {
  auto test1 = S21Matrix();
  EXPECT_EQ(test1.getCol(), 0);
  auto test2 = S21Matrix(2, 2);
  EXPECT_EQ(test2.getCol(), 2);
}

TEST(Technical, setRows) {
  auto test1 = S21Matrix(1, 1);
  test1.setRow(5);
  EXPECT_EQ(test1.getRow(), 5);
}

TEST(Technical, setRows1) {
  auto test1 = S21Matrix(1, 1);
  EXPECT_THROW(test1.setRow(-1), std::length_error);
}

TEST(Technical, setCols) {
  auto test1 = S21Matrix(1, 1);
  test1.setCol(5);
  EXPECT_EQ(test1.getCol(), 5);
}

TEST(Technical, setCols1) {
  auto test1 = S21Matrix(1, 1);
  EXPECT_THROW(test1.setCol(0), std::length_error);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

