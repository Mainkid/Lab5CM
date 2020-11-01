#pragma once
class Matrix
{
public:
	//������ ������� ����������� n �� �������
	static void Show(double** matrix, int n);

	//��������� ������� ����������� n � ��������� �� ������� ���������
	static double** GetEMatrix(int n);

	//������������ ������ A � B � ������������� n
	static double** MultMatrix(double** A, double** B, int n);

	//��������� ����� ������� A ������� n
	static double** Clone(double** A, int n);

	//��������� ������� ������� ������� n
	static double** GetZeroMatrix(int n);

	//������� ������ ��-��� ���������� �������
	static void Clear(double** A, int n);

	//��������� ������� A ������� n �� ������ x ������� n
	static double* MultMatrixByVector(double** A, double* x, int n);

	//������ ������� x ������� n � �������
	static void ShowVector(double* x, int n);

	//�������� �������� a � b � ������������� n
	static double* SubVectors(double* a, double* b, int n);

	//�������� ������ ����������� n
	static double** SubMatrixes(double** A, double** B, int n);

	//������ ����� �������
	static double Norm1(double** Matrix, int N);

	//������ ����� �������
	static double Norm2(double** Matrix, int N);

	//���������� ����������������� �������
	static double** GetTransposedMatrix(double** Matrix, int N);

	//������ ����� �������
	static double Norm3(double** Matrix, int N);

	//���������� �����
	static double VectorNorm1(double* x, int n);

	//�������������� �����
	static double VectorNorm2(double* x, int n);

	//��������� �����
	static double VectorNorm3(double* x, int n);

	//����������� �� 3-� ����
	static double MinVectorNorm(double* x, int n);

	//LU ����������
	static double* LU(double** A, double* b, int n);

	//���������� ����� ��������������� ������� A
	static double GetCond(double** A, double** AInversed, int FormNumber, int N);

	//���������� �������� ������� A ����� LU ����������
	static double** FindInverseMatrix(double** L, double** U, double** P, int n);
private:
	//���������� ������� x ����� ���������� ������� U, L � ������ b
	//�� �������� Ly = b � Ux = y
	static double* FindVectorX(double** L, double** U, double* b, int n);
};

