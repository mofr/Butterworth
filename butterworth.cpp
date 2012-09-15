#include <stdio.h>
#include <complex>
#include <cmath>
#include <iostream>

typedef double var_type;
typedef std::complex<var_type> Complex;

void zeros2coeffs(Complex * zeros, Complex * coeffs, int size)//symmetric polynoms
{

	int * b = new int[size+1];
	memset(b, 0, sizeof(int)*(size+1));
	int i = 0, j = 1;
	while( !b[size] ){
		i = 0;
		while( b[i]) b[i++] = 0; // моделируем перенос в следующий разряд, возникающий при сложении
		b[i] = 1;
 
		int subset_size = 0;
		Complex summand = 1.0;
		for(i = 0; i < size; i++){
			if(b[i]){
				summand *= zeros[i];
				++subset_size;
			}
		}
		coeffs[size-subset_size] += summand;
	}

	delete [] b;

	for(int i=(size+1)%2; i<size; i+=2){
		coeffs[i] = -coeffs[i];
	}

	coeffs[size]=Complex(1.0, 0.0);
}

void butterworth_poles(Complex * poles, int order, var_type cutoff)//cutoff Hz
{
	using namespace std;
	var_type arg;
	for(int i=0; i<order; ++i){
		//arg = (M_PI*(2*i+1))/order;
		arg = M_PI*(2*i+order+1)/2/order;
		poles[i] = Complex(std::cos(arg), std::sin(arg));
		poles[i] *= 2*M_PI*cutoff;
	}
}

void butterworth_z_zeros(Complex * zeros, int order)
{
	for(int i=0; i<order; ++i){
		zeros[i] = Complex(-1.0, 0.0);
	}
}

var_type warp_freq(var_type freq, var_type Fs)
{
	return Fs*tan(M_PI*freq/Fs)/M_PI;
}

void p2z(Complex * p, Complex * z, int size, var_type Fs)
{
	for(int i=0; i<size; ++i){
		z[i] = (Complex(2*Fs, 0.0)+p[i])/(Complex(2*Fs, 0.0)-p[i]);
	}
}


//recalc coeffs from x to 1/x polynomial variable
void inverse_poly(Complex * coeffs, int size)
{
	Complex denom = coeffs[size-1];
	for(int i=0; i<size; ++i){
		coeffs[i] /= denom;
	}
	reverse(coeffs, coeffs+size);
}

void print_line(int len=80)
{
	printf("<");
	for(int i=0; i<len-2; ++i)
		printf("=");
	printf(">\n");
}


int main()
{
	using namespace std;
	#define MAX_SIZE 30
	Complex p_poles[MAX_SIZE];
	Complex z_poles[MAX_SIZE];
	Complex b_coeffs[MAX_SIZE];

	Complex z_zeros[MAX_SIZE];
	Complex a_coeffs[MAX_SIZE];

	int order = 2;
	var_type Fs = 44100;
	var_type analog_cutoff = Fs/2*0.7;//analog freq
	var_type cutoff = warp_freq(analog_cutoff, Fs);//corresponding digital freq

	print_line();
	cout << "Butterworth low-pass filter design" << endl;
	cout << "Order: " << order << endl;
	cout << "Cutoff freq: " << analog_cutoff << " Hz" << endl;
	cout << "Digital cutoff freq: " << cutoff << " Hz" << endl;
	print_line();

	butterworth_poles(p_poles, order, cutoff);
	cout << "poles in p-domain:" << endl;
	for(int i=0; i<order; ++i){
		cout << "\t" << p_poles[i] << endl;
	}

	p2z(p_poles, z_poles, order, Fs);
	cout << "poles in z-domain:" << endl;
	for(int i=0; i<order; ++i){
		cout << "\t" << z_poles[i] << endl;
	}
	
	zeros2coeffs(z_poles, b_coeffs, order);
	cout << "z denominator coeffs:" << endl;
	for(int i=0; i<order+1; ++i){
		cout << "\t" << b_coeffs[i].real() << endl;
	}

//	Complex denom = coeffs[order];

//	inverse_poly(b_coeffs, order+1);
//	cout << "1/z denominator coeffs:" << endl;
//	for(int i=0; i<order+1; ++i){
//		cout << "\t" << b_coeffs[i].real() << endl;
//	}
	print_line();

	butterworth_z_zeros(z_zeros, order);
	cout << "zeros in z-domain:" << endl;
	for(int i=0; i<order; ++i){
		cout << "\t" << z_zeros[i] << endl;
	}
	
	zeros2coeffs(z_zeros, a_coeffs, order);
	cout << "z numerator coeffs:" << endl;
	for(int i=0; i<order+1; ++i){
		cout << "\t" << a_coeffs[i].real() << endl;
	}

//	inverse_poly(a_coeffs, order+1);
//	cout << "1/z numerator coeffs:" << endl;
//	for(int i=0; i<order+1; ++i){
//		cout << "\t" << a_coeffs[i].real() << endl;
//	}
	print_line();
	
	var_type k0_numer=0.0;
	var_type k0_denom=0.0;
	for(int i=0; i<order+1; ++i){
		k0_numer += a_coeffs[i].real();
		k0_denom += b_coeffs[i].real();
	}
	cout << "k0: " << k0_numer/k0_denom << endl;
	print_line();
	
	return 0;
}
