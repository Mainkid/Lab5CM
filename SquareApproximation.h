#pragma once

enum class variants { var20, var24 };

struct Table
{	
public:
	Table();
	~Table();
	double** m_table;
	int n = 0;
	void Create(variants variant, double a, double b, int n);
	void Print();
};

struct Functions
{
private:
	variants variant;
public:
	Functions(variants variant);
	double F(double x);
	double FG1(double x);
	double FG2(double x);
	double FG3(double x);
	double FF(double x);
};

class SquareApproximation
{
private:
	void PrintEquation(double* C, int n);
	Table* table = new Table();
	variants variant;
public:
	SquareApproximation(variants variant);
	~SquareApproximation();
	void DiscreteMethod();
	void IntegrationMethod();
	double ScalarProduct(double* x, double* y, int n);
};

