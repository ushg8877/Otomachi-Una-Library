struct MaxAdd_segt{
private:
//////////////////////////////////////////////////////////////
// basic operation
struct lazt{
	// modify here
	ll a;
	bool empty(){return a==0;}
	void clear(){a=0;}
	lazt():a(0){}
	lazt(ll _a):a(_a){}
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
	ll s;
	info():s(0){}
	info(ll _s):s(_s){}
	inline info& operator +=(const info &y){
		chkmax(s,y.s);
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
};
/////////////////////////////////////////////////////////////////
// basic segment tree template, private operations
int vl=1,vr=0;
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
		val[id]=info(); // please modify this init value
		return;
	}
	int mid=l+((ll)r-l)/2;
	build(id<<1,l,mid);build(id<<1|1,mid+1,r);
	pushup(id);
}
void update(int x,const info &v,int id,int l,int r){
	if(l==r){
		assert(x==l);
		val[id]=v;return;
	}
	int mid=l+((ll)r-l)/2;
	pushdown(id);
	if(x<=mid) update(x,v,id<<1,l,mid);
	else update(x,v,id<<1|1,mid+1,r);
	pushup(id);
}
void modify(int x,const lazt &v,int id,int l,int r){
	if(l==r){
		assert(x==l);
		val[id]+=v;return;
	}
	int mid=l+((ll)r-l)/2;
	pushdown(id);
	if(x<=mid) modify(x,v,id<<1,l,mid);
	else modify(x,v,id<<1|1,mid+1,r);
	pushup(id);
}
void modify(int ql,int qr,const lazt &v,int id,int l,int r){
	if(ql<=l&&r<=qr){push(id,v);return;}
	int mid=l+((ll)r-l)/2;
	pushdown(id);
	if(ql<=mid) modify(ql,qr,v,id<<1,l,mid);
	if(mid<qr) modify(ql,qr,v,id<<1|1,mid+1,r);
	pushup(id);
}
info query(int ql,int qr,int id,int l,int r){
	if(ql<=l&&r<=qr) return val[id];
	int mid=l+((ll)r-l)/2;
	pushdown(id);
	if(qr<=mid) return query(ql,qr,id<<1,l,mid);
	if(mid<ql) return query(ql,qr,id<<1|1,mid+1,r);
	return query(ql,qr,id<<1,l,mid)+query(ql,qr,id<<1|1,mid+1,r);
}
/////////////////////////////////////////////////////////////
// public operations, remember to setN()/setRange()
public:
void setRange(int _l,int _r){
	ll n=(ll)_r-_l+1;assert(1<=n&&n<=INT_MAX/4);
	vl=_l;vr=_r;val.assign(4*n,info());
	laz.assign(4*n,lazt());
	build(1,vl,vr);
}
void setN(int _n){
	setRange(1,_n);
}
void upd(int x,ll a){
	assert(vl<=x&&x<=vr);
	update(x,info(a),1,vl,vr);
}
void add(int x,ll a){
	assert(vl<=x&&x<=vr);
	modify(x,lazt(a),1,vl,vr);
}
void add(int ql,int qr,ll a){
	assert(vl<=ql&&(ll)ql<=1ll+qr&&qr<=vr);
	if(ql<=qr) modify(ql,qr,lazt(a),1,vl,vr);
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
};
