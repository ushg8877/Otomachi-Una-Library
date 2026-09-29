struct dynamic_segt{
private:
static const int L=-1e9,R=1e9;
int tot=0;
struct info{
	ll s;
	info():s(0){}
	info(ll _s):s(_s){}
	inline info& operator+=(const info &x){s+=x.s;return *this;}
	inline info operator+(const info &x)const{return info(*this)+=x;}
};
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
	pushup(id);return id;
}
void query(int ql,int qr,int id,int l,int r,info &X){
	if(max(ql,l)>min(r,qr)||!id) return;
	if(ql<=l&&r<=qr){X+=tr[id].val;return;}
	int mid=l+((ll)r-l)/2;
	query(ql,qr,tr[id].ls,l,mid,X);
	if(ql<=tr[id].p&&tr[id].p<=qr) X+=tr[id].self;
	query(ql,qr,tr[id].rs,mid+1,r,X);
}
public:
void set(int m=0){assert(m>=0);tot=rt=0;tr.assign(1,node());
	tr.reserve((size_t)m+1);}
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
}T;
