////////////////////////////////////////////////////////////////
//
// template for Segment Tree — 区间加 / 区间最大
// version 1.0 (Last Update Jun 27th, 2026)
//
// usage:
//   MaxAdd_segt st; st.setN(n); upd/add; ask(l,r) for range max
//   default leaf value 0
//
////////////////////////////////////////////////////////////////
struct MaxAdd_segt{
private:
struct node{ll s=0,a=0;};
int vl=1,vr=0;
vector<node> tr;
inline void pushup(int id){tr[id].s=max(tr[id<<1].s,tr[id<<1|1].s);}
inline void push(int id,ll x){
	tr[id].s+=x;
	tr[id].a+=x;
}
inline void pushdown(int id){
	ll x=tr[id].a;if(!x)return;
	push(id<<1,x);push(id<<1|1,x);
	tr[id].a=0;
}
void update(int x,ll v,bool assign,int id,int l,int r){
	if(l==r){if(assign)tr[id].s=v;else tr[id].s+=v;tr[id].a=0;return;}
	int mid=l+((ll)r-l)/2;pushdown(id);
	if(x<=mid)update(x,v,assign,id<<1,l,mid);
	else update(x,v,assign,id<<1|1,mid+1,r);
	pushup(id);
}
void modify(int ql,int qr,ll x,int id,int l,int r){
	if(ql<=l&&r<=qr){push(id,x);return;}
	int mid=l+((ll)r-l)/2;pushdown(id);
	if(ql<=mid)modify(ql,qr,x,id<<1,l,mid);
	if(mid<qr)modify(ql,qr,x,id<<1|1,mid+1,r);
	pushup(id);
}
ll query(int ql,int qr,int id,int l,int r){
	if(ql<=l&&r<=qr)return tr[id].s;
	int mid=l+((ll)r-l)/2;pushdown(id);
	if(qr<=mid)return query(ql,qr,id<<1,l,mid);
	if(mid<ql)return query(ql,qr,id<<1|1,mid+1,r);
	ll x=query(ql,qr,id<<1,l,mid),y=query(ql,qr,id<<1|1,mid+1,r);
	return max(x,y);
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
ll ask(int l,int r){assert(vl<=l&&l<=r&&r<=vr);return query(l,r,1,vl,vr);}
ll ask(int x){return ask(x,x);}

};
// end for data-structure/segment-tree/range-add-max.cpp
/////////////////////////
