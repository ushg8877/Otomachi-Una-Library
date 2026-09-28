////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/data-structure/segment-tree/dynamic.cpp
//
// usage: dynamic_segt T; T.set(); T.add(pos, val); T.ask(l, r);
//   sparse segt on value domain [L,R]
//
////////////////////////////////////////////////////////////////
struct dynamic_segt{
private:
static const int L=-1e9,R=1e9;
int tot=0;
public:
struct info{
	ll s;
	info():s(0){}
	info(ll _s):s(_s){}
	inline info& operator +=(const info &x){s+=x.s;return *this;}
	inline info operator +(const info &x)const{return info(*this)+=x;}
};
private:
struct node{
	int p=0,ls=0,rs=0;
	info self,val;
};
vector<node> tr=vector<node>(1);
int rt=0;
int new_node(int p,info x){
	assert(tr.size()<INT_MAX);tr.push_back({p,0,0,x,x});
	return ++tot;
}
inline void pushup(int id){
	tr[id].val=tr[tr[id].ls].val+tr[id].self+tr[tr[id].rs].val;
}
int add(int p,info v,int id,int l,int r){
	if(!id)return new_node(p,v);
	if(tr[id].p==p){tr[id].self+=v;pushup(id);return id;}
	int mid=l+((ll)r-l)/2;
	if(p<=mid){
		if(p>tr[id].p)swap(tr[id].p,p),swap(tr[id].self,v);
		int x=add(p,v,tr[id].ls,l,mid);tr[id].ls=x;
	}else{
		if(p<tr[id].p)swap(tr[id].p,p),swap(tr[id].self,v);
		int x=add(p,v,tr[id].rs,mid+1,r);tr[id].rs=x;
	}
	pushup(id);
	return id;
}
void query(int ql,int qr,int id,int l,int r,info &X){
	if(max(ql,l)>min(r,qr)||!id) return;
	if(ql<=l&&r<=qr){X+=tr[id].val;return;}
	int mid=l+((ll)r-l)/2;
	query(ql,qr,tr[id].ls,l,mid,X);
	if(ql<=tr[id].p&&tr[id].p<=qr) X+=tr[id].self;
	query(ql,qr,tr[id].rs,mid+1,r,X);
}
template<class F>
ll search(int q,bool rev,F &f,info &s,int id,int l,int r)const{
	if(!id||(rev?l>q:r<q)) return rev?(ll)L-1:(ll)R+1;
	if(rev?r<=q:q<=l){
		info t=rev?tr[id].val+s:s+tr[id].val;
		if(f(t)){s=t;return rev?(ll)L-1:(ll)R+1;}
	}
	int mid=l+((ll)r-l)/2;
	ll p=rev?search(q,rev,f,s,tr[id].rs,mid+1,r):
		search(q,rev,f,s,tr[id].ls,l,mid);
	if(p!=(rev?(ll)L-1:(ll)R+1)) return p;
	if(rev?tr[id].p<=q:q<=tr[id].p){
		info t=rev?tr[id].self+s:s+tr[id].self;
		if(!f(t)) return tr[id].p;
		s=t;
	}
	return rev?search(q,rev,f,s,tr[id].ls,l,mid):
		search(q,rev,f,s,tr[id].rs,mid+1,r);
}

public:
void set(int m=0){
	assert(m>=0);
	tot=rt=0;
	tr.assign(1,node());
	tr.reserve((size_t)m+1);
}
void add(int p,ll s){
	assert(L<=p&&p<=R);
	rt=add(p,info(s),rt,L,R);
}
ll ask(int l,int r){
	assert(L<=l&&l<=r&&r<=R);
	info X;
	query(l,r,rt,L,R,X);
	return X.s;
}
// Largest r with f(info of [l,r]) true; return {r,info}, empty r=l-1.
// O(log(R-L+1)); f(empty) is true, then at most one true-to-false change.
template<class F>
pair<ll,info> max_right(ll l,F f)const{
	assert(L<=l&&l<=(ll)R+1);
	info s=info();
	bool ok=f(s);assert(ok);(void)ok;
	if(l==(ll)R+1) return {l-1,s};
	ll p=search((int)l,false,f,s,rt,L,R);
	return {p-1,s};
}
// Smallest l with f(info of [l,r]) true; return {l,info}, empty l=r+1.
// Same predicate requirements; the information keeps left-to-right order.
template<class F>
pair<ll,info> min_left(ll r,F f)const{
	assert((ll)L-1<=r&&r<=R);
	info s=info();
	bool ok=f(s);assert(ok);(void)ok;
	if(r==(ll)L-1) return {r+1,s};
	ll p=search((int)r,true,f,s,rt,L,R);
	return {p+1,s};
}

}T;
// end for data-structure/segment-tree/dynamic.cpp
/////////////////////////
// !!!!! Search predicates must be monotone and accept empty info. !!!!
