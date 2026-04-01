// Тестирование метода непрерывных дробей (вещественные корни)
// Степени 5-50, только вещественные корни, double/long double/float_precision
// Результаты -> .txt

// Реализация: Павлова Анастасия, КМБО-01-22 vihuhol04@mail.ru

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

#include "Con_Frac.h"
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
double max_root_err(const vector<T> &found, const vector<T> &expected) {
	int ne=int(expected.size()), nc=int(found.size());
	if(!ne) return 0;
	double mx=0;
	vector<bool> used(nc,false);
	for(int i=0;i<ne;++i){
		double e=to_dbl(expected[i]), best=1e30; int bj=-1;
		for(int j=0;j<nc;++j){if(used[j]) continue; double d=abs(to_dbl(found[j])-e); if(d<best){best=d;bj=j;}}
		if(bj>=0){used[bj]=true; mx=max(mx, best/max(abs(e),0.01));}
		else mx=max(mx,1.0);
	}
	return mx;
}

struct TR { string name,dtype; int deg,found,expected; double err,ms; bool ok; };

template <typename T>
TR run_test(const string &nm, unsigned P, unsigned ncl,
            const vector<unsigned>&cc, const vector<T>&crad,
            const vector<pair<unsigned,unsigned>>&mg, T dcr, uint64_t seed,
            T tol, ostream &log) {
	TR r; r.name=nm; r.dtype=tname<T>(); r.deg=P;
	vector<T> coef,rrep,uq; vector<unsigned> rm; vector<complex<T>> cx;
	generate_high_degree_polynomial(P,0,ncl,cc,crad,mg,dcr,true,seed,coef,rrep,uq,rm,cx);
	// Con_Frac takes ascending coefficients (same as generator output)
	T eps=numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_PRECISE);
	r.expected=int(rrep.size());
	auto t0=chrono::high_resolution_clock::now();
	auto roots=find_roots_by_Continued_Fractions(coef,eps);
	r.ms=chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();
	r.found=int(roots.size());
	// Flatten expected: unique real roots only
	vector<T> sorted_exp=rrep;
	sort(sorted_exp.begin(),sorted_exp.end());
	vector<T> sorted_found=roots;
	sort(sorted_found.begin(),sorted_found.end());
	r.err=max_root_err(sorted_found,sorted_exp);
	r.ok=(r.err<to_dbl(tol))&&r.found>0;
	log<<"  "<<nm<<" ["<<tname<T>()<<"] deg="<<P
	   <<"  "<<r.found<<"/"<<r.expected<<"  err="<<scientific<<setprecision(3)<<r.err
	   <<"  "<<fixed<<setprecision(1)<<r.ms<<"мс  "<<(r.ok?"PASSED":"FAILED")<<endl;
	return r;
}

template <typename T> vector<TR> run_all(ostream &log) {
	vector<TR> res; T tol;
	if constexpr(is_same_v<T,float_precision>) tol=T("0.5"); else tol=T(0.5);
	log<<"\n=== НЕПРЕРЫВНЫЕ ДРОБИ: "<<tname<T>()<<" ===\n";
	// Only real roots (num_complex_pairs=0)
	for(unsigned d:{5u,10u,15u,20u,30u}){
		if constexpr(is_same_v<T,float_precision>) if(d>15) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr;
		res.push_back(run_test<T>("Веществ.",d,0,{},cr,{},dc,12345+d,tol,log));
	}
	for(unsigned d:{10u,20u}){
		if constexpr(is_same_v<T,float_precision>) if(d>10) continue;
		T r1,r2,dc; if constexpr(is_same_v<T,float_precision>){r1=T("0.01");r2=T("0.01");dc=T("0.1");}
		else{r1=T(0.01);r2=T(0.01);dc=T(0.1);}
		vector<T> cr={r1,r2};
		res.push_back(run_test<T>("Кластер.",d,2,{3,3},cr,{},dc,99999+d,tol,log));
	}
	for(unsigned d:{10u,15u}){
		if constexpr(is_same_v<T,float_precision>) if(d>10) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr; vector<pair<unsigned,unsigned>> mg={{2,2},{3,1}};
		res.push_back(run_test<T>("Кратные",d,0,{},cr,mg,dc,77777+d,tol,log));
	}
	return res;
}

void save(const vector<TR>&r,const string&fn){
	ofstream o(fn); if(!o) return;
	o<<"========================================================================\n";
	o<<"НЕПРЕРЫВНЫЕ ДРОБИ (вещественные корни)\nТип: "<<(r.empty()?"":r[0].dtype)<<"\n";
	o<<"========================================================================\n\n";
	o<<left<<setw(16)<<"Тест"<<setw(8)<<"Степ."<<setw(12)<<"Найд."<<setw(12)<<"Ожид."
	 <<setw(16)<<"Макс.ошиб."<<setw(12)<<"Время,мс"<<setw(10)<<"Статус\n"<<string(86,'-')<<"\n";
	int p=0,t=0;
	for(auto&x:r){o<<left<<setw(16)<<x.name<<setw(8)<<x.deg<<setw(12)<<x.found<<setw(12)<<x.expected
	  <<setw(16)<<scientific<<setprecision(3)<<x.err<<setw(12)<<fixed<<setprecision(1)<<x.ms
	  <<setw(10)<<(x.ok?"PASSED":"FAILED")<<"\n"; ++t; if(x.ok)++p;}
	o<<"\nИтого: "<<p<<" из "<<t<<" пройдено.\n"; o.close();
	cout<<"Сохранено: "<<fn<<endl;
}

int main(){
#ifdef _WIN32
	SetConsoleOutputCP(65001); SetConsoleCP(65001);
#endif
	cout<<"=== ТЕСТ: НЕПРЕРЫВНЫЕ ДРОБИ ===\n";
	{stringstream l; auto r=run_all<double>(l); cout<<l.str(); save(r,"cf_double.txt");}
	{stringstream l; auto r=run_all<long double>(l); cout<<l.str(); save(r,"cf_long_double.txt");}
	// float_precision не поддерживается методом непрерывных дробей (stack overflow)
	cout<<"float_precision: пропущен (не поддерживается)\n";
	cout<<"\nНепрерывные дроби — все тесты завершены.\n";
	return 0;
}
