////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/math/linear-basis/field-prefix.cpp
//
// usage: Prefix_Linear_Basis<mint> B; B.setN(n); int removed=B.insert(x,id);
//   id must increase; 0: rank increased, otherwise removed id
//
////////////////////////////////////////////////////////////////
template<typename T>
struct Prefix_Linear_Basis{
static_assert(!is_integral_v<T>&&!is_floating_point_v<T>,
	"use an exact field type");
int n=0,rank=0,last=0;
vector<vector<T>> a;
vector<int> pos;
void setN(int _n){
	assert(_n>=0);
	n=_n;
	rank=last=0;
	a.assign(n,{});
	pos.assign(n,0);
}
// Ids must be positive and strictly increasing. Return 0 if rank grows;
// otherwise return the removed id, possibly the new element itself.
int insert(vector<T> x,int id){
	assert(x.size()==(size_t)n&&id>last);last=id;
	for(int i=n-1;i>=0;i--)if(x[i]){
		if(a[i].empty()){
			T v=T(1)/x[i];
			for(int j=0;j<i;j++)x[j]*=v;
			x[i]=T(1);
			x.resize(i+1);
			a[i]=move(x);
			pos[i]=id;
			rank++;
			return 0;
		}
		if(pos[i]<id){
			x.resize(i+1);
			swap(a[i],x);
			swap(pos[i],id);
			T v=T(1)/a[i][i];
			for(int j=0;j<i;j++)a[i][j]*=v;
			a[i][i]=T(1);
		}
		T v=x[i];
		for(int j=0;j<i;j++)x[j]-=v*a[i][j];
		x[i]=T(0);
	}
	return id;
}
int get_rank(int l=1)const{
	assert(l>=1);
	int ans=0;
	for(int p:pos)ans+=p>=l;
	return ans;
}
bool contains(vector<T> x,int l=1)const{
	assert(x.size()==(size_t)n&&l>=1);
	for(int i=n-1;i>=0;i--)if(x[i]){
		if(pos[i]<l)return false;
		T v=x[i];
		for(int j=0;j<i;j++)x[j]-=v*a[i][j];
		x[i]=T(0);
	}
	return true;
}
};
// end for math/linear-basis/field-prefix.cpp
/////////////////////////
// !!!!! T must be an exact field; ids must be positive and increasing. !!!!
