template<typename T,bool MAX=false> struct Linear_RMQ{
private:
int n=0,B=1;
vector<T> a;
vector<unsigned> mask;
vector<vector<int>> st;
bool better(int x,int y)const{if constexpr(MAX)return a[y]<a[x];else return
	a[x]<a[y];}
int chk(int x,int y)const{return better(y,x)?y:x;}
int small(int l,int r)const{
	unsigned s=mask[r]&(~0u<<(l%B));
	return r/B*B+__builtin_ctz(s);
}
public:
void set(){n=0;B=1;a.clear();mask.clear();st.clear();}
void build(const vector<T> &v){
	assert(v.size()<(size_t)INT_MAX);a=v;n=a.size();mask.resize(n);st.clear();
	if(!n){B=1;return;}
	// B grows with log(n); the block sparse table has only O(n) entries.
	B=max(1,31-__builtin_clz((unsigned)n));
	int m=(n-1)/B+1;st.resize(32-__builtin_clz((unsigned)m));st[0].resize(m);
	unsigned s=0;
	for(int i=0;i<n;i++){
		if(i%B==0)s=0;
		while(s){
			int j=31-__builtin_clz(s);
			if(!better(i,i/B*B+j))break;
			s^=1u<<j;
		}
		mask[i]=s|=1u<<(i%B);
		st[0][i/B]=i/B*B+__builtin_ctz(s);
	}
	for(int k=1;k<(int)st.size();k++){
		st[k].resize(m-(1<<k)+1);
		for(int i=0;i<(int)st[k].size();i++)st[k][i]=chk(st[k-1][i],
			st[k-1][i+(1<<(k-1))]);
	}
}
int pos(int l,int r)const{
	assert(0<=l&&l<=r&&r<n);
	int x=l/B,y=r/B;if(x==y)return small(l,r);
	int ans=small(l,(x+1)*B-1);
	if(x+1<y){
		int k=31-__builtin_clz((unsigned)(y-x-1));
		ans=chk(ans,chk(st[k][x+1],st[k][y-(1<<k)]));
	}
	return chk(ans,small(y*B,r));
}
T ask(int l,int r)const{return a[pos(l,r)];}
};
