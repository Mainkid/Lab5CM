#pragma once
class Matrix
{
public:
	//ѕечать матрицы размерности n на консоль
	static void Show(double** matrix, int n);

	//ѕолучение матрицы размерности n с единицами на главной диагонали
	static double** GetEMatrix(int n);

	//ѕроизведение матриц A и B с размерност€ми n
	static double** MultMatrix(double** A, double** B, int n);

	//ѕолучение копии матрицы A размера n
	static double** Clone(double** A, int n);

	//ѕолучение нулевой матрицы размера n
	static double** GetZeroMatrix(int n);

	//ќчистка пам€ти из-под двумерного массива
	static void Clear(double** A, int n);

	//”множение матрицы A размера n на вектор x размера n
	static double* MultMatrixByVector(double** A, double* x, int n);

	//ѕечать вектора x размера n в консоли
	static void ShowVector(double* x, int n);

	//–азность векторор a и b с размерност€ми n
	static double* SubVectors(double* a, double* b, int n);

	//–азность матриц размерности n
	static double** SubMatrixes(double** A, double** B, int n);

	//ѕерва€ норма матрицы
	static double Norm1(double** Matrix, int N);

	//¬тора€ норма матрицы
	static double Norm2(double** Matrix, int N);

	//Ќахождение транспонированной матрицы
	static double** GetTransposedMatrix(double** Matrix, int N);

	//“реть€ норма матрицы
	static double Norm3(double** Matrix, int N);

	//кубическа€ норма
	static double VectorNorm1(double* x, int n);

	//октаэдрическа€ норма
	static double VectorNorm2(double* x, int n);

	//эвклидова норма
	static double VectorNorm3(double* x, int n);

	//минимальна€ из 3-х норм
	static double MinVectorNorm(double* x, int n);

	//LU разложение
	static double* LU(double** A, double* b, int n);

	//Ќахождение числа обусловленности матрицы A
	static double GetCond(double** A, double** AInversed, int FormNumber, int N);

	//Ќахождение обратной матрицы A через LU разложение
	static double** FindInverseMatrix(double** L, double** U, double** P, int n);
private:
	//Ќахождение вектора x через полученные матрицы U, L и вектор b
	//по формулам Ly = b и Ux = y
	static double* FindVectorX(double** L, double** U, double* b, int n);
};

