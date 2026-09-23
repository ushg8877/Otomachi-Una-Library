struct range_max{
int n=0;vector<vector<int>> a;
void setN(int _n){assert(1<=_n&&_n<INT_MAX);n=_n;a.assign(__lg(n)+1,vector<int>(n+1));}
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
};
struct range_min{
int n=0;vector<vector<int>> a;
void setN(int _n){assert(1<=_n&&_n<INT_MAX);n=_n;a.assign(__lg(n)+1,vector<int>(n+1));}
inline int chk(const int &x,const int &y){return min(x,y);} // compare
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
};
// Segment Tree for Xihe Tree
struct segt_minadd{
private:
//////////////////////////////////////////////////////////////
// basic operation
struct lazt{
	// modify here
	int a;
	bool empty(){return a==0;}
	void clear(){a=0;}
	lazt():a(0){}
	lazt(int _a):a(_a){}
	inline lazt& operator +=(const lazt &y){
		a+=y.a;
		return *this;
	}
	inline lazt operator +(const lazt &y)const{
		return lazt(*this)+=y;
	}
};
struct info{
	// modify here
	int s;
	info():s(0){}
	info(int _s):s(_s){}
	inline info& operator +=(const info &y){
		chkmin(s,y.s);
		return *this;
	}
	inline info operator +(const info &y)const{
		return info(*this)+=y;
	}
	inline info& operator +=(const lazt &y){
		s+=y.a;
		return *this;
	}
	inline info operator +(const lazt &y)const{
		return info(*this)+=y;
	}
	// for find_first or find_last
	inline bool operator <= (const info &a)const{
		return s<=a.s;
	}
};
/////////////////////////////////////////////////////////////////
// basic segment tree template, private operations
int n=0;
vector<info> val;
vector<lazt> laz;
inline void push(int id,const lazt &x){
	laz[id]+=x;val[id]+=x;
}
inline void pushdown(int id){
	if(!laz[id].empty()) push(id<<1,laz[id]),push(id<<1|1,laz[id]),laz[id].clear();
}
inline void pushup(int id){
	val[id]=val[id<<1]+val[id<<1|1];
}
void build(int id,int l,int r){
	laz[id]=lazt();
	if(l==r){
		val[id]=info();
		return;
	}
	int mid=(l+r)>>1;
	build(id<<1,l,mid);build(id<<1|1,mid+1,r);
	pushup(id);
}
void modify(int ql,int qr,const lazt &v,int id,int l,int r){
	if(ql<=l&&r<=qr){push(id,v);return;}
	int mid=(l+r)>>1;
	pushdown(id);
	if(ql<=mid) modify(ql,qr,v,id<<1,l,mid);
	if(mid<qr) modify(ql,qr,v,id<<1|1,mid+1,r);
	pushup(id);
}
int find_first(int id,int l,int r){
	if(l==r) return l;
	int mid=(l+r)>>1;
	pushdown(id);
	// cerr<<l<<' '<<r<<' '<<val[id<<1].s<<endl;
	return val[id<<1].s>0?find_first(id<<1|1,mid+1,r):find_first(id<<1,l,mid);
}
/////////////////////////////////////////////////////////////
// public operations, remember to init()
// please, transfer (int,int) to info here !!!!
public:
void setN(int _n){
	assert(1<=_n&&_n<=INT_MAX/4);n=_n;
	val.assign(4*n,info());laz.assign(4*n,lazt());
	build(1,1,n);
}
void add(int ql,int qr,ll a){
	// make a[x] += v, forall x in [ql,qr]
	assert(1<=ql&&ql<=qr&&qr<=n&&INT_MIN<=a&&a<=INT_MAX); 
	modify(ql,qr,lazt(a),1,1,n);
}
int find_first(){
	return find_first(1,1,n);
}
};

struct xihe_tree{
int n=0,rt=0,cnt=0,tp1=0,tp2=0,tp=0;
vector<int> a,id,L,R,M,st1,st2,st;
vector<char> vis;
enum Type{LEAF,INC,DEC,PRIME};
vector<Type> typ;
vector<vector<int>> son;
segt_minadd T;
void setN(int _n){
	assert(1<=_n&&_n<(INT_MAX-2)/2);n=_n;
	a.assign(n+2,0);id.assign(n+2,0);vis.assign(n+2,0);
	st1.assign(n+2,0);st2.assign(n+2,0);st.assign(n+2,0);
	L.assign(2*n+2,0);R.assign(2*n+2,0);M.assign(2*n+2,0);
	typ.assign(2*n+2,LEAF);son.assign(2*n+2,{});tp1=tp2=tp=rt=cnt=0;
}
range_max mx;range_min mn;
bool judge(int l,int r){
	// [l,r] 是否是一个连续段
	assert(1<=l&&l<=r&&r<=n);
	return mx.ask(l,r)-mn.ask(l,r)==r-l;
}
void build(){
	assert(1<=n);
	for(int i=1;i<=n;i++) vis[i]=false;
	for(int i=1;i<=n;i++) assert(1<=a[i]&&a[i]<=n&&!vis[a[i]]),vis[a[i]]=true;
	for(int i=0;i<=2*n;i++) son[i].clear();
	mx.setN(n);mn.setN(n);
	for(int i=1;i<=n;i++) mn.a[0][i]=mx.a[0][i]=a[i];
	mx.build();mn.build();tp1=tp2=tp=0;
	T.setN(n);cnt=0;
	for(int i=1;i<=n;i++){
		while(tp1&&a[i]<=a[st1[tp1]]) T.add(st1[tp1-1]+1,st1[tp1],a[st1[tp1]]),tp1--;
		while(tp2&&a[i]>=a[st2[tp2]]) T.add(st2[tp2-1]+1,st2[tp2],-a[st2[tp2]]),tp2--;
		T.add(st1[tp1]+1,i,-a[i]);
		T.add(st2[tp2]+1,i,a[i]);
		st1[++tp1]=st2[++tp2]=i;
		id[i]=++cnt;L[cnt]=R[cnt]=M[cnt]=i;typ[cnt]=LEAF;
		int le=T.find_first(),now=cnt;
		while(tp&&L[st[tp]]>=le){
			if((typ[st[tp]]==INC||typ[st[tp]]==DEC)&&judge(M[st[tp]],i)){
				// 当儿子
				R[st[tp]]=i;M[st[tp]]=L[now];son[st[tp]].push_back(now);now=st[tp--];
			}else if(judge(L[st[tp]],i)){
				// 新合点
				++cnt;
				if(a[L[st[tp]]]<a[i]) typ[cnt]=INC;
				else typ[cnt]=DEC;
				L[cnt]=L[st[tp]];R[cnt]=i;M[cnt]=L[now];
				son[cnt].push_back(st[tp]);
				son[cnt].push_back(now);
				now=cnt;tp--;
			}else{
				// 新析点
				++cnt;typ[cnt]=PRIME;
				son[cnt].push_back(now);
				do son[cnt].push_back(st[tp--]);
				while(tp&&!judge(L[st[tp]],i));
				L[cnt]=(tp?L[st[tp]]:L[now]);R[cnt]=i;
				if(tp) son[cnt].push_back(st[tp--]);
				reverse(son[cnt].begin(),son[cnt].end());
				now=cnt;
			}
		}
		st[++tp]=now;
		T.add(1,i,-1);
	}
	rt=st[1];
	// cerr<<"root = "<<rt<<endl;
	// for(int i=1;i<=cnt;i++){
	// 	cerr<<"vertex #"<<i<<' '<<typ[i]<<endl;
	// 	for(int j:son[i]) cerr<<i<<' '<<j<<endl;
	// }
}
void solve(){
	// solving problem from here!!!
}
};
///////////////////// END for xihe_tree.cpp ////////////////////////////////////
