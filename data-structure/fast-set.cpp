////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/data-structure/fast-set.cpp
//
// usage: FastSet fs; fs.setN(n); fs.insert(x); fs.erase(x);
//   fs[x]; fs.prev/next; needs chkmin/chkmax/first_bit/last_bit
//
////////////////////////////////////////////////////////////////
struct FastSet{
vector<vector<ull>> a;int V=0,B=0;
void setN(int _V){
	assert(0<=_V&&_V<INT_MAX);
	V=_V+1;
	a.clear();
	for(int n=V;;){
		n=(n-1)/64+1;a.emplace_back(n,0);
		if(n==1)break;
	}
	B=(int)a.size()-1;
}
inline bool operator [](int x)const{
	assert(0<=x&&x<V);
	return a[0][x>>6]>>(x&63)&1;
}
vector<int> to_vector()const{
	vector<int> I;
	if(a.empty())return I;
	for(int i=0;i<(int)a[0].size();i++)
		for(ull x=a[0][i];x;x&=x-1)I.push_back(i*64+first_bit(x));
	return I;
}
void insert(int x){
	assert(0<=x&&x<V);
	if((*this)[x]) return;
	for(int i=0;i<=B;i++){
		a[i][x>>6]|=(1ull<<(x&63));
		x>>=6;
	}
}
void erase(int x){
	assert(0<=x&&x<V);
	if(!(*this)[x]) return;
	for(int i=0;i<=B;i++){
		a[i][x>>6]&=~(1ull<<(x&63));
		if(a[i][x>>6]) break;
		x>>=6;
	}
}
int prev(int x){
	// find the first element <= x, if no, return -1
	if(x<0||a.empty())return -1;
	chkmin(x,V-1);
	for(int h=0;h<=B;h++){
		if(x<0||(x>>6)>=(int)a[h].size()) break;
		ull f=a[h][x>>6]<<(63-(x&63));
		if(!f){
			x=(x>>6)-1;
			continue;
		}
		x-=63-last_bit(f);
		for(int g=h-1;g>=0;g--){
			x<<=6;
			x+=last_bit(a[g][x>>6]);
		}
		return x;
	}
	return -1;
}
int next(int x){
	// find the first element >=x, if no, return -1
	if(x>=V||a.empty())return -1;
	chkmax(x,0);
	for(int h=0;h<=B;h++){
		if(x<0||(x>>6)>=(int)a[h].size()) break;
		ull f=a[h][x>>6]>>(x&63);
		if(!f){
			x=(x>>6)+1;
			continue;
		}
		x+=first_bit(f);
		for(int g=h-1;g>=0;g--){
			x<<=6;
			x+=first_bit(a[g][x>>6]);
		}
		return x;
	}
	return -1;
}
inline int front(){return next(0);}
inline int back(){return prev(V-1);}
};
// end for data-structure/fast-set.cpp
/////////////////////////
// !!!!! Requires chkmin/chkmax/first_bit/last_bit from basic/template.cpp. !!!!
