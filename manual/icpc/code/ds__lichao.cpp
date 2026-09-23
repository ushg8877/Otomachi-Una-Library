struct line{
	ll k,b;
	ll f(ll x){return k*x+b;}
	// please, check whether use __int128
};
struct lichao{
// you may add y=kb+x where x \in [l,r]
// you can query x, answer the **min** y (if no such, return inf = 2e18)
static const ll inf=2e18,V=1e9;
vector<int> ls{0},rs{0},p{0};int tot=0,cnt=0,rt=0;
vector<line> L{{0,inf}};
// if you don't use range add, this can up to n
inline ll gv(const int &id,const ll &v){return L[id].f(v);}
void init(){
	ls.assign(1,0);rs.assign(1,0);p.assign(1,0);
	L.assign(1,{0,inf});tot=cnt=rt=0;
}
int new_node(){
	assert(tot<INT_MAX-1);ls.push_back(0);rs.push_back(0);p.push_back(0);
	return ++tot;
}
int add(int x,int id,int l,int r){
	if(!id)id=new_node();
	if(l==r){
		if(gv(x,l)<gv(p[id],l))p[id]=x;
		return id;
	}
	int mid=l+((ll)r-l)/2;
	if(gv(x,mid)<gv(p[id],mid))swap(x,p[id]);
	if(gv(x,l)<gv(p[id],l)){int t=add(x,ls[id],l,mid);ls[id]=t;}
	if(gv(x,r)<gv(p[id],r)){int t=add(x,rs[id],mid+1,r);rs[id]=t;}
	return id;
}
int upd(int x,int ql,int qr,int id,int l,int r){
	if(max(ql,l)>min(r,qr))return id;
	if(!id)id=new_node();
	if(ql<=l&&r<=qr)return add(x,id,l,r);
	int mid=l+((ll)r-l)/2;
	int a=upd(x,ql,qr,ls[id],l,mid),b=upd(x,ql,qr,rs[id],mid+1,r);
	ls[id]=a;rs[id]=b;return id;
}
ll ask(int x,int id,int l,int r){
	if(!id) return inf;
	ll ans=gv(p[id],x);
	if(l==r) return ans;
	int mid=l+((ll)r-l)/2;
	if(x<=mid) ans=min(ans,ask(x,ls[id],l,mid));
	else ans=min(ans,ask(x,rs[id],mid+1,r));
	return ans;
}

//////////////////////////////////////////////////

void add(line x,int l=-V,int r=V){
	assert(-V<=l&&l<=r&&r<=V);
	assert(-V<=x.k&&x.k<=V);
	assert((__int128)x.k*l+x.b>=-inf&&(__int128)x.k*l+x.b<=inf);
	assert((__int128)x.k*r+x.b>=-inf&&(__int128)x.k*r+x.b<=inf);
	assert(cnt<INT_MAX-1);L.push_back(x);++cnt;
	rt=upd(cnt,l,r,rt,-V,V);
}
ll ask(int x){
	assert(-V<=x&&x<=V);
	return ask(x,rt,-V,V);
}
};
