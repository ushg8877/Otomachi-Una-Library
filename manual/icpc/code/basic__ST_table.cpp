struct ST_table{
int n=0;
vector<vector<int>> a;
void setN(int _n){
	assert(1<=_n&&_n<INT_MAX);n=_n;
	a.assign(__lg(n)+1,vector<int>(n+1));
}
inline int chk(const int &x,const int &y){return max(x,y);} // compare
void build(){
	assert(n>0);
	for(int i=1;i<(int)a.size();i++) for(int j=1;j<=n-(1<<i)+1;j++) 
		a[i][j]=chk(a[i-1][j],a[i-1][j+(1<<(i-1))]);
}
int ask(int l,int r){
	assert(1<=l&&l<=r&&r<=n);
	int k=__lg(r-l+1);
	return chk(a[k][l],a[k][r-(1<<k)+1]);
}
}S;
