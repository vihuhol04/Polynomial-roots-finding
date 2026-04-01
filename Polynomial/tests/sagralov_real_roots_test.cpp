// Тестирование метода Саградова (вещественные корни - изолирующие интервалы)
// Степени 5-50, только вещественные корни, double/long double/float_precision
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

#include "Sagralov_real_roots.h"
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
int count_covered(const vector<T> &expected_roots,
                  const vector<RootWithMultiplicity<T>> &intervals) {
	int covered=0;
	for (auto &root : expected_roots) {
		for (auto &iv : intervals) {
			T lo=min(iv.interval.first,iv.interval.second);
			T hi=max(iv.interval.first,iv.interval.second);
			if (root>=lo && root<=hi) { ++covered; break; }
		}
	}
	return covered;
}

struct TR { string name,dtype; int deg,found,expected,covered; double ms; bool ok; };

template <typename T>
TR run_test(const string &nm, unsigned P, unsigned ncl,
            const vector<unsigned>&cc, const vector<T>&crad,
            const vector<pair<unsigned,unsigned>>&mg, T dcr, uint64_t seed,
            ostream &log) {
	TR r; r.name=nm; r.dtype=tname<T>(); r.deg=P;
	vector<T> coef,rrep,uq; vector<unsigned> rm; vector<complex<T>> cx;
	generate_high_degree_polynomial(P,0,ncl,cc,crad,mg,dcr,true,seed,coef,rrep,uq,rm,cx);
	T eps=numeric_constants::adaptive_epsilon<T>(numeric_constants::EPSILON_SCALE_STANDARD);
	r.expected=int(rrep.size());
	auto t0=chrono::high_resolution_clock::now();
	auto intervals=find_real_roots_by_Sagralov(coef,eps);
	r.ms=chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();
	r.found=int(intervals.size());
	r.covered=count_covered(rrep,intervals);
	r.ok=(r.covered>=r.expected*0.7)&&r.found>0;
	log<<"  "<<nm<<" ["<<tname<T>()<<"] deg="<<P
	   <<"  интерв.="<<r.found<<"  корн.ожид.="<<r.expected<<"  покрыто="<<r.covered
	   <<"  "<<fixed<<setprecision(1)<<r.ms<<"мс  "<<(r.ok?"PASSED":"FAILED")<<endl;
	return r;
}

template <typename T> vector<TR> run_all(ostream &log) {
	vector<TR> res;
	log<<"\n=== САГРАДОВ (вещ.): "<<tname<T>()<<" ===\n";
	for(unsigned d:{5u,10u,15u,20u,30u}){
		if constexpr(is_same_v<T,float_precision>) if(d>15) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr;
		res.push_back(run_test<T>("Веществ.",d,0,{},cr,{},dc,12345+d,log));
	}
	for(unsigned d:{10u,20u}){
		if constexpr(is_same_v<T,float_precision>) if(d>10) continue;
		T r1,r2,dc; if constexpr(is_same_v<T,float_precision>){r1=T("0.01");r2=T("0.01");dc=T("0.1");}
		else{r1=T(0.01);r2=T(0.01);dc=T(0.1);}
		vector<T> cr={r1,r2};
		res.push_back(run_test<T>("Кластер.",d,2,{3,3},cr,{},dc,99999+d,log));
	}
	for(unsigned d:{10u,15u}){
		if constexpr(is_same_v<T,float_precision>) if(d>10) continue;
		T dc; if constexpr(is_same_v<T,float_precision>) dc=T("0.1"); else dc=T(0.1);
		vector<T> cr; vector<pair<unsigned,unsigned>> mg={{2,2},{3,1}};
		res.push_back(run_test<T>("Кратные",d,0,{},cr,mg,dc,77777+d,log));
	}
	if constexpr(!is_same_v<T,float_precision>)
		for(unsigned d:{40u,50u}){
			T dc=T(0.1); vector<T> cr;
			res.push_back(run_test<T>("Высок.",d,0,{},cr,{},dc,33333+d,log));
		}
	return res;
}

void save(const vector<TR>&r,const string&fn){
	ofstream o(fn); if(!o) return;
	o<<"========================================================================\n";
	o<<"САГРАДОВ (вещественные корни, изолирующие интервалы)\nТип: "<<(r.empty()?"":r[0].dtype)<<"\n";
	o<<"========================================================================\n\n";
	o<<left<<setw(16)<<"Тест"<<setw(8)<<"Степ."<<setw(10)<<"Интерв."<<setw(10)<<"Ожид."
	 <<setw(10)<<"Покрыто"<<setw(12)<<"Время,мс"<<setw(10)<<"Статус\n"<<string(76,'-')<<"\n";
	int p=0,t=0;
	for(auto&x:r){o<<left<<setw(16)<<x.name<<setw(8)<<x.deg<<setw(10)<<x.found<<setw(10)<<x.expected
	  <<setw(10)<<x.covered<<setw(12)<<fixed<<setprecision(1)<<x.ms
	  <<setw(10)<<(x.ok?"PASSED":"FAILED")<<"\n"; ++t; if(x.ok)++p;}
	o<<"\nИтого: "<<p<<" из "<<t<<" пройдено.\n"; o.close();
	cout<<"Сохранено: "<<fn<<endl;
}

int main(){
#ifdef _WIN32
	SetConsoleOutputCP(65001); SetConsoleCP(65001);
#endif
	cout<<"=== ТЕСТ: САГРАДОВ (вещ.) ===\n";
	{stringstream l; auto r=run_all<double>(l); cout<<l.str(); save(r,"sagralov_real_double.txt");}
	{stringstream l; auto r=run_all<long double>(l); cout<<l.str(); save(r,"sagralov_real_long_double.txt");}
	{stringstream l; auto r=run_all<float_precision>(l); cout<<l.str(); save(r,"sagralov_real_float_precision.txt");}
	cout<<"\nСаградов (вещ.) — все тесты завершены.\n";
	return 0;
}
