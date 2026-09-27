////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: Linear_Basis<mint> B; B.setN(n); B.insert(x); B.contains(x);
//   T must be an exact field type; coordinates are 0..n-1
//
////////////////////////////////////////////////////////////////
template<typename T>
struct Linear_Basis{
static_assert(!is_integral_v<T>&&!is_floating_point_v<T>,
	"use an exact field type");
int n=0,rank=0;
vector<vector<T>> a;
void setN(int _n){assert(_n>=0);n=_n;rank=0;a.assign(n,{});}
bool insert(vector<T> x){
	assert(x.size()==(size_t)n);
	for(int i=n-1;i>=0;i--)if(x[i]){
		if(a[i].empty()){
			T v=T(1)/x[i];
			for(int j=0;j<i;j++)x[j]*=v;
			x[i]=T(1);
			x.resize(i+1);
			a[i]=move(x);
			rank++;
			return true;
		}
		T v=x[i];
		for(int j=0;j<i;j++)x[j]-=v*a[i][j];
		x[i]=T(0);
	}
	return false;
}
bool contains(vector<T> x)const{
	assert(x.size()==(size_t)n);
	for(int i=n-1;i>=0;i--)if(x[i]){
		if(a[i].empty())return false;
		T v=x[i];
		for(int j=0;j<i;j++)x[j]-=v*a[i][j];
		x[i]=T(0);
	}
	return true;
}
};
// end for math/linear-basis/field.cpp
/////////////////////////
// !!!!! T must be an exact field type, such as mint with a prime modulus. !!!!
