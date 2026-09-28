////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/data-structure/segment-tree/range-add-sum.cpp
//
// usage: AddSum_segt st; st.setN(n);  // or set(l,r), 1-indexed
//   upd/add point; add(l,r,v) range; ask; find_first/last
//
////////////////////////////////////////////////////////////////
struct AddSum_segt{
private:
struct node{ll s=0,a=0;};
int vl=1,vr=0;
vector<node> tr;
inline void pushup(int id){tr[id].s=tr[id<<1].s+tr[id<<1|1].s;}
inline void push(int id,ll x,int len){
	tr[id].s+=x*len;
	if(len>1)tr[id].a+=x;
}
inline void pushdown(int id,int l,int r){
	ll x=tr[id].a;if(!x)return;
	int mid=l+((ll)r-l)/2;
	push(id<<1,x,mid-l+1);push(id<<1|1,x,r-mid);
	tr[id].a=0;
}
void update(int x,ll v,bool assign,int id,int l,int r){
	if(l==r){if(assign)tr[id].s=v;else tr[id].s+=v;tr[id].a=0;return;}
	int mid=l+((ll)r-l)/2;pushdown(id,l,r);
	if(x<=mid)update(x,v,assign,id<<1,l,mid);
	else update(x,v,assign,id<<1|1,mid+1,r);
	pushup(id);
}
void modify(int ql,int qr,ll x,int id,int l,int r){
	if(ql<=l&&r<=qr){push(id,x,r-l+1);return;}
	int mid=l+((ll)r-l)/2;pushdown(id,l,r);
	if(ql<=mid)modify(ql,qr,x,id<<1,l,mid);
	if(mid<qr)modify(ql,qr,x,id<<1|1,mid+1,r);
	pushup(id);
}
ll query(int ql,int qr,int id,int l,int r,ll a)const{
	// carry ancestor tags; queries do not change the tree
	if(ql<=l&&r<=qr)return tr[id].s+a*(r-l+1);
	int mid=l+((ll)r-l)/2;a+=tr[id].a;
	if(qr<=mid)return query(ql,qr,id<<1,l,mid,a);
	if(mid<ql)return query(ql,qr,id<<1|1,mid+1,r,a);
	ll x=query(ql,qr,id<<1,l,mid,a),y=query(ql,qr,id<<1|1,mid+1,r,a);
	return x+y;
}
int find(ll x,bool rev)const{
	int id=1,l=vl,r=vr;ll a=0;
	while(l<r){
		int mid=l+((ll)r-l)/2;a+=tr[id].a;
		ll s=tr[id<<1|rev].s+a*(rev?r-mid:mid-l+1);
		bool right=rev;
		if(x>s){x-=s;right=!right;}
		id=id<<1|right;
		if(right)l=mid+1;else r=mid;
	}
	return l;
}
template<class F>
ll search(int q,bool rev,F &f,ll &s,int id,int l,int r,ll a)const{
	if((rev?l>q:r<q)) return rev?(ll)vl-1:(ll)vr+1;
	if(rev?r<=q:q<=l){
		ll t=s+tr[id].s+a*(r-l+1);
		if(f(t)){s=t;return rev?(ll)vl-1:(ll)vr+1;}
		if(l==r) return l;
	}
	int mid=l+((ll)r-l)/2;a+=tr[id].a;
	ll p=rev?search(q,rev,f,s,id<<1|1,mid+1,r,a):
		search(q,rev,f,s,id<<1,l,mid,a);
	if(p!=(rev?(ll)vl-1:(ll)vr+1)) return p;
	return rev?search(q,rev,f,s,id<<1,l,mid,a):
		search(q,rev,f,s,id<<1|1,mid+1,r,a);
}

public:
void set(){vl=1;vr=0;tr.clear();}
void set(int l,int r){
	ll n=(ll)r-l+1;assert(1<=n&&n<=INT_MAX/4);
	vl=l;
	vr=r;
	tr.assign(4*n,node());
}
void setN(int n){set(1,n);}
void upd(int x,ll v){assert(vl<=x&&x<=vr);update(x,v,true,1,vl,vr);}
void add(int x,ll v){assert(vl<=x&&x<=vr);update(x,v,false,1,vl,vr);}
void add(int l,int r,ll v){
	assert(vl<=l&&(ll)l<=1ll+r&&r<=vr);
	if(l<=r)modify(l,r,v,1,vl,vr);
}
ll ask(int l,int r)const{
	assert(vl<=l&&l<=r&&r<=vr);
	return query(l,r,1,vl,vr,0);
}
ll ask(int x)const{return ask(x,x);}
// all leaf values must be non-negative
int find_first(ll x)const{
	assert(!tr.empty()&&0<x&&x<=tr[1].s);
	return find(x,false);
}
int find_last(ll x)const{
	assert(!tr.empty()&&0<x&&x<=tr[1].s);
	return find(x,true);
}
// Largest r with f(info of [l,r]) true; return {r,info}, empty r=l-1.
// O(log n). f(empty) must be true; extension may only change true to false.
template<class F>
pair<ll,ll> max_right(ll l,F f)const{
	assert(vl<=l&&l<=(ll)vr+1);
	ll s=0;
	bool ok=f(s);assert(ok);(void)ok;
	if(l==(ll)vr+1) return {l-1,s};
	ll p=search((int)l,false,f,s,1,vl,vr,0);
	return {p-1,s};
}
// Smallest l with f(info of [l,r]) true; return {l,info}, empty l=r+1.
// Same predicate requirements; the information keeps left-to-right order.
template<class F>
pair<ll,ll> min_left(ll r,F f)const{
	assert((ll)vl-1<=r&&r<=vr);
	ll s=0;
	bool ok=f(s);assert(ok);(void)ok;
	if(r==(ll)vl-1) return {r+1,s};
	ll p=search((int)r,true,f,s,1,vl,vr,0);
	return {p+1,s};
}

};
// end for data-structure/segment-tree/range-add-sum.cpp
/////////////////////////
// !!!!! Prefix-sum searches require nonnegative leaf values. !!!!
// !!!!! Search predicates must be monotone and accept empty info. !!!!
