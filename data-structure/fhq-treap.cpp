////////////////////////////////////////////////////////////////
//
// template for FHQ Treap
//
// usage:
//   FHQ t; t.setN(n); split_by_size / split_by_key / merge / merge_it
//   t.push(id, x) for lazy add on key[]
//
////////////////////////////////////////////////////////////////
mt19937_64 rnd(time(0));
struct FHQ{
int tot=0;
vector<int> ls{0},rs{0},key{0},siz{0},tag{0};
vector<ll> val{0};
inline void push(const int &id,const int &x){
	assert(1<=id&&id<=tot);
	key[id]+=x;tag[id]+=x;
}
inline void pushup(const int &id){
	siz[id]=siz[ls[id]]+siz[rs[id]]+1;
}
inline void pushdown(const int &id){
	if(tag[id]){
		if(ls[id]) push(ls[id],tag[id]);
		if(rs[id]) push(rs[id],tag[id]);
		tag[id]=0;
	}
}
void set(){setN(0);}
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);tot=_n;
	ls.assign(tot+1,0);
	rs.assign(tot+1,0);
	tag.assign(tot+1,0);
	key.assign(tot+1,0);
	siz.assign(tot+1,1);
	siz[0]=0;
	val.resize(tot+1);
	for(int i=1;i<=tot;i++)val[i]=rnd();
}
inline int newnode(int x){
	assert(tot<INT_MAX-1);++tot;
	ls.push_back(0);
	rs.push_back(0);
	tag.push_back(0);
	siz.push_back(1);
	key.push_back(x);
	val.push_back(rnd());
	return tot;
}
int leftmost(int x){
	while(ls[x]) pushdown(x),x=ls[x];
	return x;
}
int rightmost(int x){
	while(rs[x]) pushdown(x),x=rs[x];
	return x;
}
int merge(int p,int q){
	if(!p||!q) return p|q;
	if(val[p]>val[q]) {pushdown(p);rs[p]=merge(rs[p],q);pushup(p);return p;}
	else {pushdown(q);ls[q]=merge(p,ls[q]);pushup(q);return q;}
}
void split_by_key(int rt,int k,int &p,int &q){
	if(!rt){p=q=0;return;}
	pushdown(rt);
	if(k<key[rt]) q=rt,split_by_key(ls[rt],k,p,ls[rt]);
	else p=rt,split_by_key(rs[rt],k,rs[rt],q);
	pushup(rt);
}
void split_by_size(int rt,int k,int &p,int &q){
	assert(0<=k&&k<=siz[rt]);
	if(k==0){p=0;q=rt;return;}
	if(k==siz[rt]){p=rt;q=0;return;}
	pushdown(rt);
	if(k<=siz[ls[rt]]) q=rt,split_by_size(ls[rt],k,p,ls[rt]),pushup(rt);
	else p=rt,split_by_size(rs[rt],k-1-siz[ls[rt]],rs[rt],q),pushup(rt);
}
int merge_it(int p,int q){
	// ~O(logn*logV) to merge two set (maybe intersect)
	if(!p||!q) return p|q;
	if(key[rightmost(p)]<=key[leftmost(q)]) return merge(p,q);
	if(key[rightmost(q)]<=key[leftmost(p)]) return merge(q,p);
	if(val[p]<val[q]) swap(p,q);
	int r,s;
	pushdown(p);
	split_by_key(q,key[p],r,s);
	ls[p]=merge_it(ls[p],r);
	rs[p]=merge_it(rs[p],s);
	pushup(p);
	return p;
}
void alldown(int id){
	if(!id) return;
	pushdown(id);
	alldown(ls[id]);
	alldown(rs[id]);
}
};
// end for data-structure/fhq-treap.cpp
/////////////////////////
