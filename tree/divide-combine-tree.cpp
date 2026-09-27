////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: xihe_tree tr; auto t=tr.build(vector<int>{3,1,2});
//
////////////////////////////////////////////////////////////////
struct xihe_tree{
enum Type{LEAF,INC,DEC,PRIME};
struct result{
	int n=0,rt=0,cnt=0;
	vector<vector<int>> son;
	vector<Type> typ;
	vector<int> L,R,mn,mx,id;
};
int n=0;vector<int> a;
void set(){n=0;a.clear();}
void setN(int _n){assert(0<=_n&&_n<=INT_MAX/4);n=_n;a.assign(n+1,0);}
private:
struct segt{
	int n;vector<int> val,laz;
	segt(int _n):n(_n),val(4*n),laz(4*n){}
	void push(int id,int x){val[id]+=x;laz[id]+=x;}
	void pushdown(int id){
		if(laz[id]){
			push(id<<1,laz[id]);
			push(id<<1|1,laz[id]);
			laz[id]=0;
		}
	}
	void add(int ql,int qr,int x,int id,int l,int r){
		if(ql<=l&&r<=qr){push(id,x);return;}
		int mid=(l+r)>>1;pushdown(id);
		if(ql<=mid)add(ql,qr,x,id<<1,l,mid);
		if(mid<qr)add(ql,qr,x,id<<1|1,mid+1,r);
		val[id]=min(val[id<<1],val[id<<1|1]);
	}
	void add(int l,int r,int x){add(l,r,x,1,1,n);}
	int first(){
		int id=1,l=1,r=n;assert(val[1]<=0);
		while(l<r){
			int mid=(l+r)>>1;pushdown(id);
			if(val[id<<1]<=0){id<<=1;r=mid;}
			else{id=id<<1|1;l=mid+1;}
		}
		return l;
	}
};
public:
// Input a[1..n] must be a permutation of 1..n.
// Return rt, son, typ, L/R, mn/mx and leaf id; children follow position order.
// O(n log n) time and O(n) space, with no problem-specific solver.
result build()const{
	assert(a.size()==(size_t)n+1||!n);
	result res;res.n=n;
	auto &son=res.son;
	auto &typ=res.typ;
	auto &L=res.L;
	auto &R=res.R;
	auto &id=res.id;
	int &cnt=res.cnt;
	son.resize(2*n+1);
	typ.resize(2*n+1,LEAF);
	L.resize(2*n+1);
	R.resize(2*n+1);
	id.resize(n+1);
	if(!n){res.mn.resize(1);res.mx.resize(1);return res;}
	vector<char> vis(n+1);
	for(int i=1;i<=n;i++){assert(1<=a[i]&&a[i]<=n&&!vis[a[i]]);vis[a[i]]=true;}
	vector<int> v(a.begin()+1,a.end());
	Linear_RMQ<int> mn;
	Linear_RMQ<int,true> mx;
	mn.build(v);
	mx.build(v);
	auto judge=[&](int l,int r){return mx.ask(l-1,r-1)-mn.ask(l-1,r-1)==r-l;};
	vector<int> M(2*n+1),st1(n+1),st2(n+1),st(n+1);
	int tp1=0,tp2=0,tp=0;segt T(n);
	for(int i=1;i<=n;i++){
		while(tp1&&a[i]<=a[st1[tp1]]) T.add(st1[tp1-1]+1,st1[tp1],a[st1[tp1]]),
			tp1--;
		while(tp2&&a[i]>=a[st2[tp2]]) T.add(st2[tp2-1]+1,st2[tp2],-a[st2[tp2]]),
			tp2--;
		T.add(st1[tp1]+1,i,-a[i]);
		T.add(st2[tp2]+1,i,a[i]);
		st1[++tp1]=st2[++tp2]=i;
		id[i]=++cnt;
		L[cnt]=R[cnt]=M[cnt]=i;
		typ[cnt]=LEAF;
		int le=T.first(),now=cnt;
		while(tp&&L[st[tp]]>=le){
			if((typ[st[tp]]==INC||typ[st[tp]]==DEC)&&judge(M[st[tp]],i)){
				// 当儿子
				R[st[tp]]=i;
				M[st[tp]]=L[now];
				son[st[tp]].push_back(now);
				now=st[tp--];
			}else if(judge(L[st[tp]],i)){
				// 新合点
				++cnt;
				if(a[L[st[tp]]]<a[i]) typ[cnt]=INC;
				else typ[cnt]=DEC;
				L[cnt]=L[st[tp]];
				R[cnt]=i;
				M[cnt]=L[now];
				son[cnt].push_back(st[tp]);
				son[cnt].push_back(now);
				now=cnt;tp--;
			}else{
				// 新析点
				++cnt;typ[cnt]=PRIME;
				son[cnt].push_back(now);
				do son[cnt].push_back(st[tp--]);
				while(tp&&!judge(L[st[tp]],i));
				L[cnt]=(tp?L[st[tp]]:L[now]);R[cnt]=i;
				if(tp) son[cnt].push_back(st[tp--]);
				reverse(son[cnt].begin(),son[cnt].end());
				now=cnt;
			}
		}
		st[++tp]=now;
		T.add(1,i,-1);
	}
	assert(tp==1);res.rt=st[1];
	son.resize(cnt+1);
	typ.resize(cnt+1);
	L.resize(cnt+1);
	R.resize(cnt+1);
	res.mn.resize(cnt+1);res.mx.resize(cnt+1);
	for(int i=1;i<=cnt;i++){
		res.mn[i]=mn.ask(L[i]-1,R[i]-1);
		res.mx[i]=mx.ask(L[i]-1,R[i]-1);
	}
	return res;
}
result build(const vector<int> &v){
	assert(v.size()<=(size_t)INT_MAX/4);
	setN(v.size());
	copy(v.begin(),v.end(),a.begin()+1);
	return build();
}
};
// end for tree/divide-combine-tree.cpp
/////////////////////////
// !!!!! Paste data-structure/linear-rmq.cpp first. !!!!
