////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/data-structure/segment-tree/point.cpp
//
// usage: segt st; st.setN(n);  // or set(l,r), 1-indexed
//   upd/add point; ask(l,r); find_first/last
//
////////////////////////////////////////////////////////////////
// remove this line when info has no subtraction / comparison
#define USE_COMPARE
struct segt{
private:
//////////////////////////////////////////////////////////////
// basic operation
// Edit info for the problem; find_first/last require monotonicity.
struct info{
	// modify here
	ll s,len;
	info():s(0),len(0){}
	info(ll _s,ll _len):s(_s),len(_len){}
	inline info& operator +=(const info &y){
		s+=y.s;len+=y.len;
		return *this;
	}
	inline info operator +(const info &y)const{
		return info(*this)+=y;
	}
#ifdef USE_COMPARE
	// for find_first or find_last
	inline info& operator -=(const info &y){
		s-=y.s;len-=y.len;
		return *this;
	}
	inline info operator -(const info &y)const{
		return info(*this)-=y;
	}
	inline bool operator <= (const info &a)const{
		return s<=a.s;
	}
#endif
};
/////////////////////////////////////////////////////////////////
// basic segment tree template, private operations
int vl=1,vr=0;
vector<info> val;
inline void pushup(int id){
	val[id]=val[id<<1]+val[id<<1|1];
}
void build(int id,int l,int r){
	if(l==r){
		val[id]=info(0,1); // please modify this init value
		return;
	}
	int mid=l+((ll)r-l)/2;
	build(id<<1,l,mid);build(id<<1|1,mid+1,r);
	pushup(id);
}
void update(int x,const info &v,int id,int l,int r){
	if(l==r){
		assert(x==l);
		val[id]=v;
		return;
	}
	int mid=l+((ll)r-l)/2;
	if(x<=mid) update(x,v,id<<1,l,mid);
	else update(x,v,id<<1|1,mid+1,r);
	pushup(id);
}
void modify(int x,const info &v,int id,int l,int r){
	if(l==r){
		assert(x==l);
		val[id]+=v;
		return;
	}
	int mid=l+((ll)r-l)/2;
	if(x<=mid) modify(x,v,id<<1,l,mid);
	else modify(x,v,id<<1|1,mid+1,r);
	pushup(id);
}
info query(int ql,int qr,int id,int l,int r){
	if(ql<=l&&r<=qr) return val[id];
	int mid=l+((ll)r-l)/2;
	if(qr<=mid) return query(ql,qr,id<<1,l,mid);
	if(mid<ql) return query(ql,qr,id<<1|1,mid+1,r);
	return query(ql,qr,id<<1,l,mid)+query(ql,qr,id<<1|1,mid+1,r);
}
#ifdef USE_COMPARE
int find_first(info x,int id,int l,int r){
	// if no compare betweeen infos, delete this
	if(l==r) return l;
	int mid=l+((ll)r-l)/2;
	return x<=val[id<<1]?find_first(x,id<<1,l,mid):find_first(x-val[id<<1],
		id<<1|1,mid+1,r);
}
int find_last(info x,int id,int l,int r){
	// if no compare betweeen infos, delete this
	if(l==r) return l;
	int mid=l+((ll)r-l)/2;
	return x<=val[id<<1|1]?find_last(x,id<<1|1,mid+1,r):
		find_last(x-val[id<<1|1],id<<1,l,mid);
}
#endif
/////////////////////////////////////////////////////////////
// public operations, remember to setN()/set()
// convert arguments to info / lazt, and extract the answer here
public:
void set(){vl=1;vr=0;val.clear();}
void set(int _l,int _r){
	ll n=(ll)_r-_l+1;assert(1<=n&&n<=INT_MAX/4);
	vl=_l;
	vr=_r;
	val.assign(4*n,info());
	build(1,vl,vr);
}
void setN(int _n){
	set(1,_n);
}
void upd(int x,ll a){
	assert(vl<=x&&x<=vr);
	update(x,info(a,1),1,vl,vr);
}
void add(int x,ll a){
	assert(vl<=x&&x<=vr);
	modify(x,info(a,0),1,vl,vr);
}
ll ask(int ql,int qr){
	assert(vl<=ql&&ql<=qr&&qr<=vr);
	return query(ql,qr,1,vl,vr).s;
}
ll ask(int x){
	assert(vl<=x&&x<=vr);
	// query a[x]
	return query(x,x,1,vl,vr).s;
}
#ifdef USE_COMPARE
int find_first(ll x){
	assert(!val.empty()&&x>0&&info(x,0)<=val[1]);
	return find_first(info(x,0),1,vl,vr);
}
int find_last(ll x){
	assert(!val.empty()&&x>0&&info(x,0)<=val[1]);
	return find_last(info(x,0),1,vl,vr);
}
#endif
};

#undef USE_COMPARE
// end for data-structure/segment-tree/point.cpp
/////////////////////////
