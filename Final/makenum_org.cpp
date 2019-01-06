#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

int main(int argc, char *argv[]){
	ofstream fout;
	fout.open("s.txt");
	srand(time(NULL));
	int n;
	/*cout << "Enter the number of digits of s (n >=2): " << endl;
	cin >> n;*/

	n = atoi(argv[1]);

	fout << (rand()%9)+1;
	for(int i = 1; i <= n-1; i++){
		fout << rand()%10;
	}

	fout.close();

	fout.open("t.txt");
	/*cout << "Enter the number of digits of t (n >= 2): " << endl;
	cin >> n;*/
	n = atoi(argv[2]);
	fout << (rand()%9)+1;
	for(int i = 1; i <= n-1; i++){
		fout << rand()%10;
	}

	fout.close();
	
	return 0;
}