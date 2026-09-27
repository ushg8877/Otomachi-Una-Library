////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/math/gaussian-elimination.cpp
//
// usage: Gauss<mint> g; int res=g.solve(a,m);
//   g.x; g.rank; g.det(a);
//
////////////////////////////////////////////////////////////////
template<typename T>
struct Gauss{
static_assert(!is_integral_v<T>,"use a field type");
long double eps=1e-12L;
int rank=0;
vector<int> where;
vector<T> x;
bool zero(const T &v)const{
	if constexpr(is_floating_point_v<T>)return abs(v)<=eps;
	else return !v;
}
int pivot(const vector<vector<T>> &a,int r,int c)const{
	int p=-1;
	for(int i=r;i<(int)a.size();i++)if(!zero(a[i][c])){
		if constexpr(is_floating_point_v<T>){
			if(p==-1||abs(a[i][c])>abs(a[p][c]))p=i;
		}
		else return i;
	}
	return p;
}
// Augmented matrix with m variables; T is a field or floating point.
// Return -1 / 0 / 1 for none / multiple / unique; x stores one solution.
// Free variables are zero; rank / where store elimination information.
int solve(vector<vector<T>> a,int m){
	assert(m>=0&&m<INT_MAX&&a.size()<INT_MAX&&eps>=0);
	int n=a.size();
	for(auto &v:a)assert(v.size()==(size_t)m+1);
	rank=0;
	where.assign(m,-1);
	x.assign(m,T(0));
	for(int c=0;c<m&&rank<n;c++){
		int p=pivot(a,rank,c);if(p==-1)continue;
		swap(a[p],a[rank]);where[c]=rank;
		T z=T(1)/a[rank][c];a[rank][c]=T(1);
		for(int j=c+1;j<=m;j++)a[rank][j]*=z;
		for(int i=rank+1;i<n;i++){
			T v=a[i][c];
			a[i][c]=T(0);
			if(zero(v))continue;
			for(int j=c+1;j<=m;j++)a[i][j]-=v*a[rank][j];
		}
		rank++;
	}
	for(int i=rank;i<n;i++)if(!zero(a[i][m])){x.clear();return -1;}
	for(int c=m-1;c>=0;c--)if(where[c]!=-1){
		int r=where[c];x[c]=a[r][m];
		for(int j=c+1;j<m;j++)x[c]-=a[r][j]*x[j];
	}
	return rank==m?1:0;
}
// Determinant of a square matrix; floating-point zero tests use eps.
T det(vector<vector<T>> a)const{
	assert(a.size()<INT_MAX&&eps>=0);int n=a.size();
	for(auto &v:a)assert(v.size()==(size_t)n);
	T ans=1;
	for(int c=0;c<n;c++){
		int p=pivot(a,c,c);if(p==-1)return T(0);
		if(p!=c){swap(a[p],a[c]);ans=-ans;}
		ans*=a[c][c];T z=T(1)/a[c][c];
		for(int i=c+1;i<n;i++){
			T v=a[i][c];
			a[i][c]=T(0);
			if(zero(v))continue;
			v*=z;
			for(int j=c+1;j<n;j++)a[i][j]-=v*a[c][j];
		}
	}
	return ans;
}
};
// end for math/gaussian-elimination.cpp
/////////////////////////
