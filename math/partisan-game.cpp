////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/math/partisan-game.cpp
//
// usage: Partisan_Game g; g.add("one",{"zero"},{}); g.print();
//   auto x=g.get("one"); cout<<x+x; cout<<x.winner('R');
//
////////////////////////////////////////////////////////////////
struct Dyadic{
BigInt p=0;
int k=0; // p/2^k, reduced; do not modify p or k directly.
Dyadic(BigInt _p=0,int _k=0):p(move(_p)),k(_k){
	assert(k>=0);
	if(!p){k=0;return;}
	while(k&&p.a[0]%2==0){p/=2;k--;}
}
static BigInt mul2(BigInt x,int k){
	assert(k>=0);
	while(k>=29){x=x.mul(1<<29);k-=29;}
	return x.mul(1<<k);
}
Dyadic operator-()const{return Dyadic(-p,k);}
Dyadic operator+(const Dyadic &b)const{
	int t=max(k,b.k);
	return Dyadic(mul2(p,t-k)+mul2(b.p,t-b.k),t);
}
Dyadic operator-(const Dyadic &b)const{return *this+(-b);}
Dyadic& operator+=(const Dyadic &b){return *this=*this+b;}
Dyadic& operator-=(const Dyadic &b){return *this=*this-b;}
int cmp(const Dyadic &b)const{
	int t=max(k,b.k);
	return mul2(p,t-k).cmp(mul2(b.p,t-b.k));
}
bool operator==(const Dyadic &b)const{return k==b.k&&p==b.p;}
bool operator!=(const Dyadic &b)const{return !(*this==b);}
bool operator<(const Dyadic &b)const{return cmp(b)<0;}
bool operator>(const Dyadic &b)const{return cmp(b)>0;}
bool operator<=(const Dyadic &b)const{return cmp(b)<=0;}
bool operator>=(const Dyadic &b)const{return cmp(b)>=0;}
// Simplest dyadic strictly between l and r; nullopt means no bound.
static Dyadic between(optional<Dyadic> l,optional<Dyadic> r){
	assert(!l||!r||*l<*r);
	Dyadic zero;
	if((!l||*l<zero)&&(!r||zero<*r)) return zero;
	if(r&&*r<=zero) return -between(-*r,l?optional<Dyadic>(-*l):nullopt);
	// Here l>=0. Find the first denominator with a point in (l,r).
	auto next=[&](int t){
		BigInt p=t>=l->k?mul2(l->p,t-l->k):
			l->p/mul2(1,l->k-t);
		return Dyadic(p+1,t);
	};
	if(!r) return next(0);
	assert(max(l->k,r->k)<INT_MAX);
	int lo=0,hi=max(l->k,r->k)+1;
	while(lo<hi){
		int mid=lo+(hi-lo)/2;
		if(next(mid)<*r) hi=mid;
		else lo=mid+1;
	}
	return next(lo);
}
friend ostream& operator<<(ostream &o,const Dyadic &x){
	o<<x.p;
	if(x.k) o<<'/'<<mul2(1,x.k);
	return o;
}
};
struct Partisan_Game{
private:
struct node{
	vector<int> l,r;
	optional<Dyadic> num;
	int nim=-1;
};
vector<node> a;
map<pair<vector<int>,vector<int>>,int> pool;
map<pair<int,int>,bool> cmp;
map<pair<int,int>,int> sum;
map<int,int> neg;
map<string,int> id;
vector<string> name;
vector<array<vector<int>,2>> edg;
vector<int> val;
bool ready=false;
int intern(vector<int> l,vector<int> r){
	for(auto p:{&l,&r}){
		sort(p->begin(),p->end());
		p->erase(unique(p->begin(),p->end()),p->end());
	}
	auto key=make_pair(l,r);
	auto it=pool.find(key);
	if(it!=pool.end()) return it->second;
	assert(a.size()<INT_MAX);
	int u=a.size();
	pool.emplace(move(key),u);a.push_back({move(l),move(r),nullopt,-1});
	return u;
}
// G<=H iff no G^L>=H and no H^R<=G.
bool le(int u,int v){
	if(u==v) return true;
	if(a[u].num&&a[v].num) return *a[u].num<=*a[v].num;
	auto key=make_pair(u,v);
	auto it=cmp.find(key);
	if(it!=cmp.end()) return it->second;
	for(int x:a[u].l) if(le(v,x)) return cmp[key]=false;
	for(int x:a[v].r) if(le(x,u)) return cmp[key]=false;
	return cmp[key]=true;
}
void prune(vector<int> &v,bool right){
	vector<int> b;
	for(int x:v){
		bool bad=false;
		for(int y:b) if(right?le(y,x):le(x,y)){bad=true;break;}
		if(bad) continue;
		int t=0;
		for(int y:b) if(!(right?le(x,y):le(y,x))) b[t++]=y;
		b.resize(t);b.push_back(x);
	}
	v=move(b);
}
int make(vector<int> l,vector<int> r){
	while(true){
		prune(l,false);prune(r,true);
		int u=intern(l,r);
		bool changed=false;
		for(int side=0;side<2&&!changed;side++){
			auto &v=side?r:l;
			for(int i=0;i<(int)v.size()&&!changed;i++){
				auto &reply=side?a[v[i]].l:a[v[i]].r;
				for(int x:reply) if(side?le(u,x):le(x,u)){
					// Bypass the reversible option through this reply.
					auto &to=side?a[x].r:a[x].l;
					v.erase(v.begin()+i);v.insert(v.end(),to.begin(),to.end());
					changed=true;break;
				}
			}
		}
		if(changed) continue;
		optional<Dyadic> lo,hi;
		bool number=true;
		for(int x:l){
			if(!a[x].num){number=false;break;}
			if(!lo||*lo<*a[x].num) lo=a[x].num;
		}
		for(int x:r){
			if(!a[x].num){number=false;break;}
			if(!hi||*a[x].num<*hi) hi=a[x].num;
		}
		if(number&&(!lo||!hi||*lo<*hi)) a[u].num=Dyadic::between(lo,hi);
		if(a[u].l==a[u].r){
			vector<int> vis(a[u].l.size()+1);
			bool impartial=true;
			for(int x:a[u].l){
				if(a[x].nim<0){impartial=false;break;}
				if(a[x].nim<(int)vis.size()) vis[a[x].nim]=1;
			}
			if(impartial){
				int k=0;
				while(vis[k]) k++;
				a[u].nim=k;
			}
		}
		return u;
	}
}
int add(int u,int v){
	if(u>v) swap(u,v);
	if(!u) return v;
	auto key=make_pair(u,v);
	auto it=sum.find(key);
	if(it!=sum.end()) return it->second;
	// Copies are needed: recursive calls can reallocate a.
	auto x=a[u],y=a[v];
	vector<int> l,r;
	for(int w:x.l) l.push_back(add(w,v));
	for(int w:y.l) l.push_back(add(u,w));
	for(int w:x.r) r.push_back(add(w,v));
	for(int w:y.r) r.push_back(add(u,w));
	return sum[key]=make(move(l),move(r));
}
int minus(int u){
	auto it=neg.find(u);
	if(it!=neg.end()) return it->second;
	auto x=a[u];
	vector<int> l,r;
	for(int v:x.r) l.push_back(minus(v));
	for(int v:x.l) r.push_back(minus(v));
	int v=make(move(l),move(r));
	neg[u]=v;neg[v]=u;
	return v;
}
int node_id(const string &s){
	auto it=id.find(s);
	if(it!=id.end()) return it->second;
	assert(name.size()<INT_MAX);
	int u=name.size();
	id[s]=u;name.push_back(s);edg.push_back({});ready=false;
	return u;
}
void print_value(int u,ostream &o)const{
	if(a[u].num){o<<*a[u].num;return;}
	if(a[u].nim>=0){o<<'*';if(a[u].nim!=1)o<<a[u].nim;return;}
	o<<'{';
	for(int side=0;side<2;side++){
		if(side) o<<'|';
		auto &v=side?a[u].r:a[u].l;
		for(int i=0;i<(int)v.size();i++){
			if(i) o<<',';
			print_value(v[i],o);
		}
	}
	o<<'}';
}
public:
struct Game{
private:
	Partisan_Game *g;
	int u;
	friend struct Partisan_Game;
	Game(Partisan_Game *_g,int _u):g(_g),u(_u){}
public:
	// Disjunctive sum: each turn moves in exactly one component.
	Game operator+(const Game &b)const{
		assert(g==b.g);return Game(g,g->add(u,b.u));
	}
	Game operator-()const{return Game(g,g->minus(u));}
	Game operator-(const Game &b)const{return *this+(-b);}
	Game& operator+=(const Game &b){return *this=*this+b;}
	Game& operator-=(const Game &b){return *this=*this-b;}
	// Partial order: e.g. * and 0 are incomparable (both <= are false).
	bool operator<=(const Game &b)const{
		assert(g==b.g);return g->le(u,b.u);
	}
	bool operator>=(const Game &b)const{return b<=*this;}
	bool operator==(const Game &b)const{return *this<=b&&b<=*this;}
	bool operator!=(const Game &b)const{return !(*this==b);}
	bool operator<(const Game &b)const{return *this<=b&&!(b<=*this);}
	bool operator>(const Game &b)const{return b<*this;}
	// Return an exact dyadic, or nullopt if this game is not a number.
	optional<Dyadic> number()const{return g->a[u].num;}
	// Normal play: the player with no move loses. Return 'L' or 'R'.
	char winner(char first='L')const{
		assert(first=='L'||first=='R');
		return first=='L'?(g->le(u,0)?'R':'L'):(g->le(0,u)?'L':'R');
	}
	void print(ostream &o)const{g->print_value(u,o);}
	friend ostream& operator<<(ostream &o,const Game &x){x.print(o);return o;}
};
Partisan_Game(){set();}
// Clear everything. Previously returned Game objects become invalid.
void set(){
	a.clear();pool.clear();cmp.clear();sum.clear();neg.clear();id.clear();
	name.clear();edg.clear();val.clear();ready=false;
	int u=intern({},{});a[u].num=Dyadic();a[u].nim=0;neg[0]=0;
}
// Set/replace both option lists. Unspecified labels are terminal nodes.
void add(const string &s,const vector<string> &l,const vector<string> &r){
	int u=node_id(s);
	vector<int> x,y;
	for(auto &t:l) x.push_back(node_id(t));
	for(auto &t:r) y.push_back(node_id(t));
	edg[u]={move(x),move(y)};ready=false;
}
// Evaluate the DAG from sinks to sources; called automatically by get/print.
// General canonicalization and sums may take exponential time and space.
void solve(){
	if(ready) return;
	int n=name.size();
	vector<int> deg(n),q;
	vector<vector<int>> rev(n);
	for(int u=0;u<n;u++) for(int side=0;side<2;side++){
		for(int v:edg[u][side]){deg[u]++;rev[v].push_back(u);}
	}
	for(int u=0;u<n;u++) if(!deg[u]) q.push_back(u);
	val.resize(n);
	for(int i=0;i<(int)q.size();i++){
		int u=q[i];
		vector<int> l,r;
		for(int v:edg[u][0]) l.push_back(val[v]);
		for(int v:edg[u][1]) r.push_back(val[v]);
		val[u]=make(move(l),move(r));
		for(int v:rev[u]) if(!--deg[v]) q.push_back(v);
	}
	assert((int)q.size()==n);ready=true;
}
// Query an existing label. Game objects retain their values after add().
Game get(const string &s){
	auto it=id.find(s);assert(it!=id.end());
	solve();return Game(this,val[it->second]);
}
// Print every label, canonical value, and winner for each starting player.
void print(ostream &o=cout){
	solve();
	for(int i=0;i<(int)name.size();i++){
		Game x(this,val[i]);
		o<<name[i]<<" = "<<x<<"; L: "<<x.winner('L');
		o<<", R: "<<x.winner('R')<<'\n';
	}
}
};
// end for math/partisan-game.cpp
/////////////////////////
// !!!!! Paste math/big-int.cpp first. Input must be a DAG. !!!!
// !!!!! Small games only; keep the owning Partisan_Game alive and unmoved. !!!!
