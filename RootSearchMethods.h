#pragma once
class RootSearchMethods
{
private:
	static void PrintTheHead(int i);
	static void PrintIntermediateResults(int Itr, double q, double normaNevyazki, double ErrorEstimate, double x, double y,double F1,double F2);
	static void PrintIntermediateResults(int Itr, double normaNevyazki, double x, double y, double F1, double F2);
	static void PrintIntermediateResults1(int itr, double x, double y, double alpha, double normaNevyazki, double f1, double f2, double ff);
	static double FunctionCount(double x, double y, int variant, int equationPos);
	static void Jacobian(double**J, double x, double y, int variant);
	static double ErrorNorm(double x1, double x2, double y1, double y2);
	static double F1(double x, double y, double variant);
	static double F2(double x, double y, double variant);
	static double** CalculateMatrixF(double x, double y,int var);
	static void CalculateReverseMatrixF(double** F);
	static double CalculateTheError(double x, double y, double xNext, double yNext, double q);
	static double FunctionFF(double x, double y,int var);
	static double* CalculateVectorDF(double x, double y,int var);

public:
	static void NewtonMethod(double x, double y, double eps,int var);
	static void SimpleIterationMethod(int variant, double epsilon, double x, double y);
	static void GradientDescentMethods(double x, double y,int var,double eps);
};