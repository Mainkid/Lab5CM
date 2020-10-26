#include "Matrix.h"
#include <iostream>
#include<iomanip>
#include<math.h>

void Matrix::Clear(double** A, int n)
{
	for (int i = 0; i < n; i++)
		delete[] A[i];
	delete[] A;
}

void Matrix::Show(double** matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << " ";
			std::cout << std::setw(12) << std::left << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

double** Matrix::GetEMatrix(int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
	{
		res[i] = new double[n]();
		res[i][i] = 1;
	}
	return res;
}

double** Matrix::MultMatrix(double** A, double** B, int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
		res[i] = new double[n]();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				res[i][j] += A[i][k] * B[k][j];

	return res;
}

double** Matrix::Clone(double** A, int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
		res[i] = new double[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			res[i][j] = A[i][j];

	return res;
}

double** Matrix::GetZeroMatrix(int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
		res[i] = new double[n]();
	return res;
}

double* Matrix::MultMatrixByVector(double** A, double* x, int n)
{
	double* b = new double[n]();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			b[i] += A[i][j] * x[j];

	return b;
}

void Matrix::ShowVector(double* x, int n)
{
	for (int i = 0; i < n; i++)
	{
		std::cout << " ";
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;
}

double* Matrix::SubVectors(double* a, double* b, int n)
{
	double* res = new double[n];
	for (int i = 0; i < n; i++)
	{
		res[i] = a[i] - b[i];
	}
	return res;
}

double Matrix::Norm1(double** Matrix, int N)
{
	double max = 0;
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		sum = 0;

		for (int j = 0; j < N; j++)
		{
			sum += abs(Matrix[i][j]);


		}
		if (sum > max)
			max = sum;
	}
	return max;
}

double Matrix::Norm2(double** Matrix, int N)
{
	double max = 0;
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		sum = 0;

		for (int j = 0; j < N; j++)
		{
			sum += abs(Matrix[j][i]);
		}
		if (sum > max)
			max = sum;
	}
	return max;
}

double** Matrix::GetTransposedMatrix(double** Matrix, int N)
{
	double** TransposedMatrix = new double* [N];

	for (int i = 0; i < N; i++)
		TransposedMatrix[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			TransposedMatrix[i][j] = Matrix[j][i];
		}
	return TransposedMatrix;
}

double Matrix::Norm3(double** Matrix, int N)
{
	double** TransposedMatrix = GetTransposedMatrix(Matrix, N);
	double** MultipliedMatrix = MultMatrix(TransposedMatrix, Matrix, N);
	double** resultedMatrix = Clone(MultipliedMatrix, N);
	double* result = new double[N];

	double tg_2_alpha = 0;
	double sin_alpha = 0;
	double cos_alpha = 0;
	double maxElement = 1;
	double epsilon = 0.001;
	double maxValue = MultipliedMatrix[0][1];
	int iRotation = 0;
	int jRotation = 0;
	double** E_Matrix = GetEMatrix(N);
	double** E_Matrix_Transposed = NULL;
	int k = 0;


	while (maxElement > epsilon)
	{
		k++;

		maxValue = 0;
		tg_2_alpha = 0;
		sin_alpha = 0;
		cos_alpha = 0;
		maxElement = 0;
		iRotation = 0;
		jRotation = 0;

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				if (i == j)
					E_Matrix[i][j] = 1;
				else
					E_Matrix[i][j] = 0;

		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				if (abs(resultedMatrix[i][j]) > maxValue)
				{
					tg_2_alpha = 2 * resultedMatrix[i][j] / (resultedMatrix[i][i] - resultedMatrix[j][j]);
					iRotation = i;
					jRotation = j;
					maxValue = abs(resultedMatrix[i][j]);
				}

		sin_alpha = sin(0.5 * atan(tg_2_alpha));
		cos_alpha = cos(0.5 * atan(tg_2_alpha));

		E_Matrix[iRotation][iRotation] = cos_alpha;
		E_Matrix[jRotation][jRotation] = cos_alpha;
		E_Matrix[iRotation][jRotation] = -sin_alpha;
		E_Matrix[jRotation][iRotation] = sin_alpha;

		E_Matrix_Transposed = GetTransposedMatrix(E_Matrix, N);

		double** tmp;
		double** tmp2;

		tmp = MultMatrix(E_Matrix_Transposed, resultedMatrix, N);
		Clear(resultedMatrix, N);
		tmp2 = MultMatrix(tmp, E_Matrix, N);

		resultedMatrix = Clone(tmp2, N);

		Clear(tmp, N);
		Clear(tmp2, N);

		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				if (abs(resultedMatrix[i][j]) > maxElement)
				{
					maxElement = abs(resultedMatrix[i][j]);
				}

		Clear(E_Matrix_Transposed, N);
	}

	for (int i = 0; i < N; i++)
		result[i] = resultedMatrix[i][i];

	Clear(resultedMatrix, N);
	Clear(MultipliedMatrix, N);
	Clear(TransposedMatrix, N);
	Clear(E_Matrix, N);

	double max = abs(result[0]);
	for (int i = 1; i < N; i++)
		if (abs(result[i]) > max)
			max = abs(result[i]);
	delete[] result;

	return sqrt(max);
}

double** Matrix::SubMatrixes(double** A, double** B, int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
		res[i] = new double[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			res[i][j] = A[i][j] - B[i][j];

	return res;
}

double Matrix::VectorNorm1(double* x, int n)
{
	double max = abs(x[0]);

	for (int i = 1; i < n; i++)
		if (abs(x[i]) > max)
			max = abs(x[i]);

	return max;
}

double Matrix::VectorNorm2(double* x, int n)
{
	double sum = 0.0;

	for (int i = 0; i < n; i++)
		sum += abs(x[i]);

	return sum;
}

double Matrix::VectorNorm3(double* x, int n)
{
	double res = 0.0;

	for (int i = 0; i < n; i++)
		res += x[i] * x[i];

	return sqrt(res);
}

double Matrix::MinVectorNorm(double* x, int n)
{
	double norm1 = VectorNorm1(x, n);
	double norm2 = VectorNorm2(x, n);
	double norm3 = VectorNorm3(x, n);

	double min = norm1;
	if (norm2 < min)
		min = norm2;
	if (norm3 < min)
		min = norm3;

	return min;
}

double* Matrix::FindVectorX(double** L, double** U, double* b, int n)
{
	double* x = new double[n]();
	double* y = new double[n]();

	for (int i = 0; i < n; i++)
	{
		y[i] = b[i];
		for (int j = 0; j < i; j++)
		{
			y[i] -= y[j] * L[i][j];
		}
		y[i] /= L[i][i];
	}

	for (int i = n - 1; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= U[i][j] * x[j];
		}
	}

	delete[] y;
	return x;
}

double* Matrix::LU(double** A, double** L, double** U, double** P, double* b, int n)
{
	double** A_Clone = Clone(A, n);

	for (int i = 0; i < n; i++)
	{
		int indexOfMaxElem = i;
		double max = abs(U[i][i]);

		for (int j = i + 1; j < n; j++)
			if (abs(U[j][i]) > max)
			{
				max = abs(U[j][i]);
				indexOfMaxElem = j;
			}

		if (max < 0.00001)
		{
			return NULL;
		}

		if (i != indexOfMaxElem)
		{
			std::swap(U[i], U[indexOfMaxElem]);
			std::swap(A_Clone[i], A_Clone[indexOfMaxElem]);
			std::swap(P[i], P[indexOfMaxElem]);
			std::swap(b[i], b[indexOfMaxElem]);
		}

		double diagElem = U[i][i];
		for (int j = 0; j < n; j++)
		{
			U[i][j] /= diagElem;
		}

		for (int j = i + 1; j < n; j++)
		{
			double d = U[j][i];
			for (int k = 0; k < n; k++)
			{
				U[j][k] -= d * U[i][k];
			}
		}

		L[i][0] = A_Clone[i][0];
		for (int j = 1; j < n; j++)
		{
			if (j > i)
				break;
			L[i][j] += A_Clone[i][j];
			for (int k = 0; k < j; k++)
				L[i][j] -= L[i][k] * U[k][j];
		}

	}
	Clear(A_Clone, n);

	return FindVectorX(L, U, b, n);
}

double** Matrix::FindInverseMatrix(double** L, double** U, double** P, int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
		res[i] = new double[n];
	double* e = new double[n]();
	double* x;

	for (int i = 0; i < n; i++)
	{
		if (i > 0)
			e[i - 1] = 0;

		e[i] = 1;

		x = FindVectorX(L, U, e, n);
		for (int j = 0; j < n; j++)
			res[j][i] = x[j];

		delete[] x;
	}

	double** final_res = Matrix::MultMatrix(res, P, n);
	Clear(res, n);
	delete[] e;

	return final_res;
}

double Matrix::GetCond(double** A, double** AInversed, int FormNumber, int N)
{
	if (FormNumber == 1) 
	{
		return Norm1(A, N) * Norm1(AInversed, N);
	}
	else if (FormNumber == 2)
	{
		return Norm2(A, N) * Norm2(AInversed, N);
	}
	else if (FormNumber == 3)
	{
		double maxValue1 = Norm3(A, N);
		double maxValue2 = Norm3(AInversed, N);

		return maxValue1 * maxValue2;
	}
}