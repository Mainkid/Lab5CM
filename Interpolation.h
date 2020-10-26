#pragma once

class Interpolation {

public:
	static double F20(double x);
	static double F24(double x);
	static void MakeDifTable(double a, double b, double step, int variant);
	static double RecF(int start, int finish, double* Arr,int variant);



};