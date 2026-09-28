////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/data-structure/segment-tree/merge.cpp
//
// usage: segt_merge sm; sm.setN(n); int rt=0;
//   sm.add(rt,i,v); sm.ask(rt,x); sm.merge(a,b);
//
////////////////////////////////////////////////////////////////
struct segt_merge{
public:
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
private:
vector<info> val=vector<info>(1);
vector<int> ls{0},rs{0};int n=0,tot=0;
inline void pushup(int id){val[id]=val[ls[id]]+val[rs[id]];}
int add(int x,info v,int id,int l,int r){
	if(!id){
		assert(tot<INT_MAX-1);id=++tot;
		val.emplace_back();
		ls.push_back(0);
		rs.push_back(0);
	}
	if(l==r){val[id]+=v;return id;}
	int mid=l+(r-l)/2;
	if(x<=mid){int t=add(x,v,ls[id],l,mid);ls[id]=t;}
	else {int t=add(x,v,rs[id],mid+1,r);rs[id]=t;}
	pushup(id);
	return id;
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
template<class F>
ll search(int q,bool rev,F &f,info &s,int id,int l,int r)const{
	if(!id||(rev?l>q:r<q)) return rev?(ll)0-1:(ll)n+1;
	if(rev?r<=q:q<=l){
		info t=rev?val[id]+s:s+val[id];
		if(f(t)){s=t;return rev?(ll)0-1:(ll)n+1;}
		if(l==r) return l;
	}
	int mid=l+((ll)r-l)/2;
	ll p=rev?search(q,rev,f,s,rs[id],mid+1,r):
		search(q,rev,f,s,ls[id],l,mid);
	if(p!=(rev?(ll)0-1:(ll)n+1)) return p;
	return rev?search(q,rev,f,s,ls[id],l,mid):
		search(q,rev,f,s,rs[id],mid+1,r);
}

public:
void setN(int _n){
	assert(_n>=0);
	n=_n;
	tot=0;
	val.assign(1,info());
	ls.assign(1,0);
	rs.assign(1,0);
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
// Largest r with f(info of [l,r]) true; return {r,info}, empty r=l-1.
// O(log n). f(empty) must be true; extension may only change true to false.
template<class F>
pair<ll,info> max_right(int rt,ll l,F f)const{
	assert(0<=rt&&rt<=tot&&0<=l&&l<=(ll)n+1);
	info s=info();
	bool ok=f(s);assert(ok);(void)ok;
	if(l==(ll)n+1) return {l-1,s};
	ll p=search((int)l,false,f,s,rt,0,n);
	return {p-1,s};
}
// Smallest l with f(info of [l,r]) true; return {l,info}, empty l=r+1.
// Same predicate requirements; the information keeps left-to-right order.
template<class F>
pair<ll,info> min_left(int rt,ll r,F f)const{
	assert(0<=rt&&rt<=tot&&(ll)0-1<=r&&r<=n);
	info s=info();
	bool ok=f(s);assert(ok);(void)ok;
	if(r==(ll)0-1) return {r+1,s};
	ll p=search((int)r,true,f,s,rt,0,n);
	return {p+1,s};
}

};
// end for data-structure/segment-tree/merge.cpp
/////////////////////////
// !!!!! merge consumes roots; do not reuse the old separate trees. !!!!
// !!!!! Search predicates must be monotone and accept empty info. !!!!
