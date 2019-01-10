#include<iostream>
#include<fstream>
using namespace std;


int main()
{
	ifstream infile, infile_o, infile_p;
	infile.open("serial_size_time.txt",ios::in);
	infile_p.open("pthread_size_time.txt",ios::in);
	infile_o.open("openmp_size_time.txt",ios::in);
	ofstream outfile_p, outfile_o;
	outfile_p.open("efficiency_pthread.txt",ios::out);
	outfile_o.open("efficiency_openmp.txt",ios::out);
	double serial,pthread,openmp;
	for(int i = 0; i < 46; i++)
	{
		infile>>serial;
		infile_p>>pthread;
		infile_o>>openmp;
		pthread = (double)serial/pthread;
		pthread = pthread / 4.0;
		openmp = (double)serial/openmp;
		openmp = openmp / 4.0;
		outfile_p<<pthread<<endl;
		outfile_o<<openmp<<endl;
	}
	infile.close();
	infile_p.close();
	infile_o.close();
	outfile_p.close();	
	outfile_o.close();
	return 0;
 } 
