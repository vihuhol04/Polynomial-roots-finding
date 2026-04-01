// Тестирование метода Вилфа
// Степени 5-50, вещ./компл./кластериз./кратные, только double/long double
// (Wilf не поддерживает float_precision — static_assert is_floating_point)
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

#include "Wilf.h"
#include "NumericConstants.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

template <typename T> double to_dbl(const T &v) { return static_cast<double>(v); }
template <typename T> string tname() {
	if constexpr (is_same_v<T,double>) return "double";
	else return "long double";
}

template <typename T>
double residual(const vector<complex<T>> &roots, const vector<T> &coef_asc) {
	double mx=0;
	for (auto &z : roots) {
		complex<T> val(T(0),T(0)), pw(T(1),T(0));
		for (size_t i=0; i<coef_asc.size(); ++i) { val+=pw*complex<T>(coef_asc[i],T(0)); pw*=z; }
		double r=to_dbl(sqrt(val.real()*val.real()+val.imag()*val.imag()));
		mx=max(mx,r);
	}
	return mx;
}

struct TR { string name,dtype; int deg,found,expected; double res,ms; bool ok; };

template <typename T>
TR run_test(const string &nm, unsigned P, unsigned ncp, unsigned ncl,
            const vector<unsigned>&cc, const vector<T>&crad,
            const vector<pair<unsigned,unsigned>>&mg, T dcr, uint64_t seed,
            double rtol, ostream &log) {
	TR r; r.name=nm; r.dtype=tname<T>(); r.deg=P;
	vector<T> coef,rrep,uq; vector<unsigned> rm; vector<complex<T>> cx;
	generate_high_degree_polynomial(P,ncp,ncl,cc,crad,mg,dcr,true,seed,coef,rrep,uq,rm,cx);
	// Wilf takes ascending + explicit eps
	T eps=T(1e-10);
	r.expected=int(P);
	auto t0=chrono::high_resolution_clock::now();
	auto roots=find_roots_by_Wilf(coef,eps);
	r.ms=chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();
	r.found=int(roots.size());
	r.res=residual(roots,coef);
	r.ok=(r.res<rtol)&&r.found>0;
	log<<"  "<<nm<<" ["<<tname<T>()<<"] deg="<<P
	   <<"  "<<r.found<<"/"<<r.expected<<"  |P(z)|="<<scientific<<setprecision(3)<<r.res
	   <<"  "<<fixed<<setprecision(1)<<r.ms<<"мс  "<<(r.ok?"PASSED":"FAILED")<<endl;
	return r;
}

template <typename T> vector<TR> run_all(ostream &log) {
	vector<TR> res; double rtol=1.0;
	log<<"\n=== ВИЛФ: "<<tname<T>()<<" ===\n";
	for(unsigned d:{5u,8u,10u,12u}){
		T dc=T(0.1); vector<T> cr;
		res.push_back(run_test<T>("Веществ.",d,0,0,{},cr,{},dc,12345+d,rtol,log));
	}
	for(unsigned d:{6u,8u,10u}){
		T dc=T(0.1); vector<T> cr;
		res.push_back(run_test<T>("Компл.",d,d/4,0,{},cr,{},dc,54321+d,rtol,log));
	}
	for(unsigned d:{8u,10u}){
		T dc=T(0.1); vector<T> cr={T(0.01),T(0.01)};
		res.push_back(run_test<T>("Кластер.",d,0,2,{3,3},cr,{},dc,99999+d,rtol,log));
	}
	for(unsigned d:{8u,10u}){
		T dc=T(0.1); vector<T> cr;
		vector<pair<unsigned,unsigned>> mg={{2,2},{3,1}};
		res.push_back(run_test<T>("Кратные",d,0,0,{},cr,mg,dc,77777+d,rtol,log));
	}
	for(unsigned d:{15u,20u}){
		T dc=T(0.1); vector<T> cr;
		res.push_back(run_test<T>("Высок.",d,d/5,0,{},cr,{},dc,33333+d,rtol,log));
	}
	return res;
}

void save(const vector<TR>&r,const string&fn){
	ofstream o(fn); if(!o) return;
	o<<"========================================================================\n";
	o<<"ВИЛФ\nТип: "<<(r.empty()?"":r[0].dtype)<<"\n";
	o<<"========================================================================\n\n";
	o<<left<<setw(16)<<"Тест"<<setw(8)<<"Степ."<<setw(12)<<"Найд."<<setw(12)<<"Ожид."
	 <<setw(16)<<"Макс.|P(z)|"<<setw(12)<<"Время,мс"<<setw(10)<<"Статус\n"<<string(86,'-')<<"\n";
	int p=0,t=0;
	for(auto&x:r){o<<left<<setw(16)<<x.name<<setw(8)<<x.deg<<setw(12)<<x.found<<setw(12)<<x.expected
	  <<setw(16)<<scientific<<setprecision(3)<<x.res<<setw(12)<<fixed<<setprecision(1)<<x.ms
	  <<setw(10)<<(x.ok?"PASSED":"FAILED")<<"\n"; ++t; if(x.ok)++p;}
	o<<"\nИтого: "<<p<<" из "<<t<<" пройдено.\n"; o.close();
	cout<<"Сохранено: "<<fn<<endl;
}

int main(){
#ifdef _WIN32
	SetConsoleOutputCP(65001); SetConsoleCP(65001);
#endif
	cout<<"=== ТЕСТ: ВИЛФ ===\n";
	// Wilf supports only float/double/long double (no float_precision)
	{stringstream l; auto r=run_all<double>(l); cout<<l.str(); save(r,"wilf_double.txt");}
	{stringstream l; auto r=run_all<long double>(l); cout<<l.str(); save(r,"wilf_long_double.txt");}
	cout<<"\nВилф — все тесты завершены.\n";
	return 0;
}
