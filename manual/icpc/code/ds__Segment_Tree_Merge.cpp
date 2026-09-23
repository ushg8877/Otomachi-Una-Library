struct segt_merge{
private:
struct info{
	// modify here
	ll s,s1;
	info():s(0),s1(0){}
	info(ll _s,ll _s1):s(_s),s1(_s1){}
	inline info& operator +=(const info &y){
		s+=y.s;s1+=y.s1;
		return *this;
	}
	inline info operator +(const info &y)const{
		return info(*this)+=y;
	}
};
vector<info> val=vector<info>(1);
vector<int> ls{0},rs{0};int n=0,tot=0;
inline void pushup(int id){val[id]=val[ls[id]]+val[rs[id]];}
int add(int x,info v,int id,int l,int r){
	if(!id){
		assert(tot<INT_MAX-1);id=++tot;
		val.emplace_back();ls.push_back(0);rs.push_back(0);
	}
	if(l==r){val[id]+=v;return id;}
	int mid=l+(r-l)/2;
	if(x<=mid){int t=add(x,v,ls[id],l,mid);ls[id]=t;}
	else {int t=add(x,v,rs[id],mid+1,r);rs[id]=t;}
	pushup(id);return id;
}
info query(int ql,int qr,int id,int l,int r){
	if(!id||max(l,ql)>min(r,qr)) return info(); 
	if(ql<=l&&r<=qr) return val[id];
	int mid=l+(r-l)/2;
	return query(ql,qr,ls[id],l,mid)+query(ql,qr,rs[id],mid+1,r);
}
void merge(int &u,int &v,int l,int r){
	if(!u||!v){u|=v;return;}
	if(l==r){val[u]+=val[v];return;}
	int mid=l+(r-l)/2;
	merge(ls[u],ls[v],l,mid);
	merge(rs[u],rs[v],mid+1,r);
	pushup(u);
}
public:
void setN(int _n){
	assert(_n>=0);n=_n;tot=0;
	val.assign(1,info());ls.assign(1,0);rs.assign(1,0);
}
void add(int &id,int x,int o){
	assert(0<=x&&x<=n);
	assert(0<=id&&id<=tot);
	id=add(x,info(o,1ll*x*o),id,0,n);
}
ll ask(int &id,int x){
	assert(0<=id&&id<=tot&&0<=x&&x<=n);
	auto it=query(x,n,id,0,n);
	return it.s1-it.s*x;
}
void merge(int &x,int &y){
	assert(0<=x&&x<=tot&&0<=y&&y<=tot&&(&x!=&y)&&(!x||x!=y));
	merge(x,y,0,n);
	y=0;
}
};
