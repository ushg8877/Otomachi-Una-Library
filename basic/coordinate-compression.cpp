////////////////////////////////////////////////////////////////
//
// template for Coordinate Compression
// version 2.0 (Last Update Jul 3rd, 2026)
//
// usage: compress<int> cc; cc.set(); cc.pb(x); cc.build();
//   cc.id(x); cc[i]; cc.range(l,r);  // call set() each problem
//
////////////////////////////////////////////////////////////////
template<typename T>
struct compress{
vector<T> I;
vector<T*> II;
// make compress values to 1 ~ size
void set(){I.clear();II.clear();}
void pb(T &x){I.push_back(x);II.push_back(&x);}
void pb(const T &x){I.push_back(x);}
int id(T x)const{
	auto it=lower_bound(I.begin(),I.end(),x);
	assert(it!=I.end()&&*it==x);
	return it-I.begin()+1;
}
// Inserting new values and rebuilding may change every compressed id.
void build(){
	sort(I.begin(),I.end());
	I.erase(unique(I.begin(),I.end()),I.end());
	vector<int> a; a.reserve(II.size());
	for(auto p:II)a.push_back(id(*p));
	for(int i=0;i<(int)II.size();i++)*II[i]=a[i];
	II.clear();
}
int size()const{return I.size();}
array<int,2> range(const T &l,const T &r)const{
	int a=lower_bound(I.begin(),I.end(),l)-I.begin();
	int b=upper_bound(I.begin(),I.end(),r)-I.begin();
	return {a+1,b};
}
T operator [](int x)const{
	assert(1<=x&&x<=I.size());
	return I[x-1];
}
};
// end for basic/coordinate-compression.cpp
/////////////////////////
