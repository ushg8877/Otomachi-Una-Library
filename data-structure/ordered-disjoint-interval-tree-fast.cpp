////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/data-structure/ordered-disjoint-interval-tree-fast.cpp
//
// usage: Fast_ODT fot; fot.setN(n); fot.insert(l,r,c); fot.extract(l,r);
//   requires data-structure/fast-set.cpp
//
////////////////////////////////////////////////////////////////
struct Fast_ODT{
private:
vector<int> a;int n=0;
FastSet S;
void cut(int y){
	if(S[y]) return;
	a[y]=a[S.next(y)];
	S.insert(y);
}
public:
void setN(int _n){
	assert(1<=_n&&_n<INT_MAX);
	n=_n;
	a.assign(n+1,0);
	S.setN(n);
	S.insert(n);
	a[n]=-1;
}

vector<array<int,3>> extract(int l,int r){
	// extract segments from [l,r]
	// before you use this function, pay attention is the TC correct?
	assert(1<=l&&l<=r&&r<=n);
	cut(l-1);cut(r);
	vector<array<int,3>> I;
	for(int i=l,j;i<=r;){
		j=S.next(i);
		I.push_back({i,j,a[j]});
		i=j+1;
	}
	return I;
}
vector<array<int,3>> insert(int l,int r,int c){
	// erase segments from [l,r]
	// return vector for segs, [l,r,c]
	assert(1<=l&&l<=r&&r<=n);
	cut(l-1);cut(r);
	vector<array<int,3>> I;
	for(int i=l,j;i<=r;){
		j=S.next(i);
		S.erase(j);
		I.push_back({i,j,a[j]});
		i=j+1;
	}
	a[r]=c;S.insert(r);
	return I;
}
};
// end for data-structure/ordered-disjoint-interval-tree-fast.cpp
/////////////////////////
// !!!!! Paste data-structure/fast-set.cpp first. !!!!
