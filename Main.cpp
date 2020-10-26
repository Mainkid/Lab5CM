#include<iostream>
#include <iomanip>
#include<Windows.h>
#include <fstream>
#include <ctime>
#include "Matrix.h"
#include "RootSearchMethods.h"
#include "Interpolation.h"

using namespace std;

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	
	Interpolation::MakeDifTable(1, 2, 0.2, 20);
	
	return 0;
}