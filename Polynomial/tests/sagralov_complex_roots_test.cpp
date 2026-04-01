// Тестирование CIsolate (комплексная изоляция корней - диски)
// Степени 5-30, вещ./компл./кластериз./кратные, double/long double
// CIsolate работает с long double по умолчанию; float_precision не поддерживается
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

#include "Sagralov_complex_roots.h"
#include "NumericConstants.h"
#include "generate_high_degree_polynomial.h"

using namespace std;

template <typename T> double to_dbl(const T &v) { return static_cast<double>(v); }

template <typename Tin>
string tname_in() {
	if constexpr (is_same_v<Tin,double>) return "double";
	else return "long double";
}

template <typename Tin, typename T>
int count_covered(const vector<Tin> &real_roots,
                  const vector<complex<Tin>> &complex_roots,
                  const vector<Disk<T>> &disks) {
	int covered=0;
	auto in_disk = [&](double re, double im) -> bool {
		for (auto &d : disks) {
			double dx=re-to_dbl(d.center.real());
			double dy=im-to_dbl(d.center.imag());
			if (sqrt(dx*dx+dy*dy) <= to_dbl(d.radius)*1.1+1e-8)
				return true;
		}
		return false;
	};
	for (auto &r : real_roots)
		if (in_disk(to_dbl(r), 0.0)) ++covered;
	for (auto &c : complex_roots)
		if (in_disk(to_dbl(c.real()), to_dbl(c.imag()))) ++covered;
	return covered;
}

struct TR { string name,dtype; int deg,disks,expected,covered,total_roots_in_disks; double ms; bool ok; };

template <typename Tin, typename T=long double>
TR run_test(const string &nm, unsigned P, unsigned ncp, unsigned ncl,
            const vector<unsigned>&cc, const vector<Tin>&crad,
            const vector<pair<unsigned,unsigned>>&mg, Tin dcr, uint64_t seed,
            ostream &log) {
	TR r; r.name=nm; r.dtype=tname_in<Tin>(); r.deg=P;
	vector<Tin> coef,rrep,uq; vector<unsigned> rm; vector<complex<Tin>> cx;
	generate_high_degree_polynomial(P,ncp,ncl,cc,crad,mg,dcr,true,seed,coef,rrep,uq,rm,cx);
	r.expected=int(rrep.size()+cx.size());
	auto t0=chrono::high_resolution_clock::now();
	auto disks=CIsolate<Tin,T>(coef);
	r.ms=chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();
	r.disks=int(disks.size());
	r.total_roots_in_disks=0;
	for(auto &d : disks) r.total_roots_in_disks+=d.num_roots;
	r.covered=count_covered(rrep,cx,disks);
	r.ok=(r.covered>=int(r.expected*0.5))&&r.disks>0;
	log<<"  "<<nm<<" ["<<tname_in<Tin>()<<"] deg="<<P
	   <<"  дисков="<<r.disks<<"  корн.в дисках="<<r.total_roots_in_disks
	   <<"  ожид.="<<r.expected<<"  покрыто="<<r.covered
	   <<"  "<<fixed<<setprecision(1)<<r.ms<<"мс  "<<(r.ok?"PASSED":"FAILED")<<endl;
	return r;
}

template <typename Tin> vector<TR> run_all(ostream &log) {
	vector<TR> res;
	log<<"\n=== CIsolate (компл.диски): "<<tname_in<Tin>()<<" ===\n";
	for(unsigned d:{5u,10u,15u,20u}){
		Tin dc=Tin(0.1); vector<Tin> cr;
		res.push_back(run_test<Tin>("Веществ.",d,0,0,{},cr,{},dc,12345+d,log));
	}
	for(unsigned d:{6u,10u,20u}){
		Tin dc=Tin(0.1); vector<Tin> cr;
		res.push_back(run_test<Tin>("Компл.",d,d/4,0,{},cr,{},dc,54321+d,log));
	}
	for(unsigned d:{10u,15u}){
		Tin dc=Tin(0.1); vector<Tin> cr={Tin(0.01),Tin(0.01)};
		res.push_back(run_test<Tin>("Кластер.",d,0,2,{3,3},cr,{},dc,99999+d,log));
	}
	for(unsigned d:{10u,15u}){
		Tin dc=Tin(0.1); vector<Tin> cr;
		vector<pair<unsigned,unsigned>> mg={{2,2}};
		res.push_back(run_test<Tin>("Кратные",d,0,0,{},cr,mg,dc,77777+d,log));
	}
	for(unsigned d:{25u,30u}){
		Tin dc=Tin(0.1); vector<Tin> cr;
		res.push_back(run_test<Tin>("Высок.",d,d/5,0,{},cr,{},dc,33333+d,log));
	}
	return res;
}

void save(const vector<TR>&r,const string&fn){
	ofstream o(fn); if(!o) return;
	o<<"========================================================================\n";
	o<<"CIsolate (комплексные корни, изолирующие диски)\nТип: "<<(r.empty()?"":r[0].dtype)<<"\n";
	o<<"========================================================================\n\n";
	o<<left<<setw(16)<<"Тест"<<setw(8)<<"Степ."<<setw(10)<<"Дисков"<<setw(12)<<"Корн.дисков"
	 <<setw(10)<<"Ожид."<<setw(10)<<"Покрыто"<<setw(12)<<"Время,мс"<<setw(10)<<"Статус\n"<<string(88,'-')<<"\n";
	int p=0,t=0;
	for(auto&x:r){o<<left<<setw(16)<<x.name<<setw(8)<<x.deg<<setw(10)<<x.disks
	  <<setw(12)<<x.total_roots_in_disks<<setw(10)<<x.expected<<setw(10)<<x.covered
	  <<setw(12)<<fixed<<setprecision(1)<<x.ms
	  <<setw(10)<<(x.ok?"PASSED":"FAILED")<<"\n"; ++t; if(x.ok)++p;}
	o<<"\nИтого: "<<p<<" из "<<t<<" пройдено.\n"; o.close();
	cout<<"Сохранено: "<<fn<<endl;
}

int main(){
#ifdef _WIN32
	SetConsoleOutputCP(65001); SetConsoleCP(65001);
#endif
	cout<<"=== ТЕСТ: CIsolate (компл.) ===\n";
	{stringstream l; auto r=run_all<double>(l); cout<<l.str(); save(r,"cisolate_double.txt");}
	{stringstream l; auto r=run_all<long double>(l); cout<<l.str(); save(r,"cisolate_long_double.txt");}
	cout<<"\nCIsolate — все тесты завершены.\n";
	return 0;
}
