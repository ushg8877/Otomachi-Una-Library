////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/data-structure/priority-queue-two-stacks.cpp
//
// usage: pq2stack<T> pq; pq.push(x,w); pq.prep(); pq.pop();
//   pop removes max weight; needs chkmin; call prep() before each query
//
////////////////////////////////////////////////////////////////
template<typename T>
struct pq2stack{
private:
struct info{
	T e;
	int w;
	bool tag;
	inline bool operator < (const info &x)const{return w<x.w;}
};
int n=0,r=0,d=0;
vector<info> I;
public:
void set(){n=r=d=0;I.clear();}
vector<T> prep(){
	vector<T> ans;
	while(d<I.size()) ans.push_back(I[d++].e);
	return ans;
}
void push(T x,int w){
	I.push_back({x,w,0});
}
void imerge(int i,int j,int k){
	if(i==j||j==k) return;
	static vector<info> y;
	y.resize(k-i);
	merge(I.begin()+i,I.begin()+j,I.begin()+j,I.begin()+k,y.begin());
	for(int t=i;t<k;t++) I[t]=y[t-i];
}
int pop(){
	// pq2stack core
	assert(!I.empty());
	int tar=(r?r&-r:inf<int>),mn=I.size()-1,hd=mn;
	for(int i=I.size()-1;i>=0;i--){
		if(I[i].tag){
			I[i].tag=false;
			tar--;
		}
		if(I[i]<I[mn]) mn=i;
		else hd=i;
		if(tar==0) break;
	}
	sort(I.begin()+n,I.end());
	int cut=I.size();
	for(int i=I.size()-1;i>hd;i--) if(I[i]<I[i-1]){
		imerge(i,cut,I.size());
		cut=i;
	}
	imerge(hd,cut,I.size());
	int ans=max(0,d-hd);
	int mx=(int)I.size()-1;
	for(int i=n;i<(int)I.size();i++) if(I[i].w>I[mx].w) mx=i;
	if(mx!=(int)I.size()-1) swap(I[mx],I.back());
	I.pop_back();
	chkmin(d,hd);
	n=I.size();
	if(r==0){
		for(int i=0;i<n;i++) I[i].tag=true;
		r=n;
	}else{
		for(int i=n-(r&-r)+1;i<n;i++) I[i].tag=true;
		r--;
	}
	return ans;
}
};
// end for data-structure/priority-queue-two-stacks.cpp
/////////////////////////
// !!!!! Requires chkmin; call prep() before each query. !!!!
