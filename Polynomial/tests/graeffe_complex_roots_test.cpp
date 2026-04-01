// Тестирование Хосмана (Грефе для комплексных корней)
// Степени 5-30, вещ./компл./кластериз./кратные, double/long double/float_precision
// Результаты -> .txt

#ifdef _WIN32
#include <windows.h>
#endif

#include <algorithm>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "Graeffe_complex_roots.h"
#include "NumericConstants.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

template <typename T> double to_dbl(const T &v) { return static_cast<double>(v); }
template <typename T> string tname() {
	if constexpr (is_same_v<T,double>) return "double";
	else if constexpr (is_same_v<T,long double>) return "long double";
	else return "float_precision";
}

template <typename T>
double residual_check(const vector<complex<T>> &roots, const vector<T> &coef_asc) {
	double mx=0;
	for (auto &z : roots) {
		complex<T> val(T(0),T(0));
		complex<T> pw(T(1),T(0));
		for (size_t i=0; i<coef_asc.size(); ++i) {
			val += pw * complex<T>(coef_asc[i], T(0));
			pw *= z;
		}
		double r = to_dbl(sqrt_val(val.real()*val.real()+val.imag()*val.imag()));
		mx = max(mx, r);
	}
	return mx;
}

struct TR { string name,dtype; int deg,found,expected; double residual,ms; bool ok; };

template <typename T>
TR run_test(const string &nm, unsigned P, unsigned ncp, unsigned ncl,
            const vector<unsigned>&cc, const vector<T>&crad,
            const vector<pair<unsigned,unsigned>>&mg, T dcr, uint64_t seed,
            double res_tol, int maxIt, ostream &log) {
	TR r; r.name=nm; r.dtype=tname<T>(); r.deg=P;
	vector<T> coef,rrep,uq; vector<unsigned> rm; vector<complex<T>> cx;
	generate_high_degree_polynomial(P,ncp,ncl,cc,crad,mg,dcr,true,seed,coef,rrep,uq,rm,cx);
	vector<T> desc(coef.rbegin(),coef.rend());
	T eps=numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_PRECISE);
	r.expected=int(rrep.size()+cx.size());
	auto t0=chrono::high_resolution_clock::now();
	auto roots=hosman_modification_graeffe(desc,eps,maxIt);
	r.ms=chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();
	r.found=int(roots.size());
	r.residual=residual_check(roots,coef);
	r.ok=(r.residual<res_tol)&&r.found>0;
	log<<"  "<<nm<<" ["<<tname<T>()<<"] deg="<<P
	   <<"  "<<r.found<<"/"<<r.expected<<"  |P(z)|="<<scientific<<setprecision(3)<<r.residual
	   <<"  "<<fixed<<setprecision(1)<<r.ms<<"мс  "<<(r.ok?"PASSED":"FAILED")<<endl;
	return r;
}

template <typename T> vector<TR> run_all(ostream &log) {
	vector<TR> res; int mi=50; double rtol=1.0;
	if constexpr(is_same_v<T,float_precision>) mi=20;
	log<<"\n=== ХОСМАН (компл.): "<<tname<T>()<<" ===\n";
	for(unsigned d:{5u,10u,15u,20u}){
		if constexpr(is_same_v<T,float_precision>) if(d>15) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr; res.push_back(run_test<T>("Веществ.",d,0,0,{},cr,{},dc,12345+d,rtol,mi,log));
	}
	for(unsigned d:{6u,10u,20u}){
		if constexpr(is_same_v<T,float_precision>) if(d>10) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr; res.push_back(run_test<T>("Компл.",d,d/4,0,{},cr,{},dc,54321+d,rtol,mi,log));
	}
	for(unsigned d:{10u,20u}){
		if constexpr(is_same_v<T,float_precision>) if(d>10) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr; vector<pair<unsigned,unsigned>> mg={{2,2}};
		res.push_back(run_test<T>("Кратные",d,0,0,{},cr,mg,dc,77777+d,rtol,mi,log));
	}
	return res;
}

void save(const vector<TR>&r,const string&fn){
	ofstream o(fn); if(!o) return;
	o<<"========================================================================\n";
	o<<"ХОСМАН (комплексные корни)\nТип: "<<(r.empty()?"":r[0].dtype)<<"\n";
	o<<"========================================================================\n\n";
	o<<left<<setw(16)<<"Тест"<<setw(8)<<"Степ."<<setw(12)<<"Найд."<<setw(12)<<"Ожид."
	 <<setw(16)<<"Макс.|P(z)|"<<setw(12)<<"Время,мс"<<setw(10)<<"Статус\n"<<string(86,'-')<<"\n";
	int p=0,t=0;
	for(auto&x:r){o<<left<<setw(16)<<x.name<<setw(8)<<x.deg<<setw(12)<<x.found<<setw(12)<<x.expected
	  <<setw(16)<<scientific<<setprecision(3)<<x.residual<<setw(12)<<fixed<<setprecision(1)<<x.ms
	  <<setw(10)<<(x.ok?"PASSED":"FAILED")<<"\n"; ++t; if(x.ok)++p;}
	o<<"\nИтого: "<<p<<" из "<<t<<" пройдено.\n"; o.close();
	cout<<"Сохранено: "<<fn<<endl;
}

int main(){
#ifdef _WIN32
	SetConsoleOutputCP(65001); SetConsoleCP(65001);
#endif
	cout<<"=== ТЕСТ: ХОСМАН (компл.) ===\n";
	{stringstream l; auto r=run_all<double>(l); cout<<l.str(); save(r,"graeffe_complex_double.txt");}
	{stringstream l; auto r=run_all<long double>(l); cout<<l.str(); save(r,"graeffe_complex_long_double.txt");}
	{stringstream l; auto r=run_all<float_precision>(l); cout<<l.str(); save(r,"graeffe_complex_float_precision.txt");}
	cout<<"\nХосман — все тесты завершены.\n";
	return 0;
}
