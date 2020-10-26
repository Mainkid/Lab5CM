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
	

	Interpolation::MakeDifTable(1, 2, 0.2,19);
	
	return 0;
}