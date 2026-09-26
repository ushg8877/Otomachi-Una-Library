////////////////////////////////////////////////////////////////
//
// template for vector operations
// usage:
//   auto a=vrange(1,6); // {1,2,3,4,5,6}, [l,r]
//   auto b=vcat(a,a); // concatenate
//   a+x, a-x, a*x, a/x, a%x; a&x, a|x, a^x, a<<x, a>>x
//   a==x, a!=x, a<x, a<=x, a>x, a>=x, a&&x, a||x -> vector<bool>
//   +a, -a, ~a, !a; a+=x, a-=x, ...; ++a, a++, --a, a--
//   vadd(a,2); vsub(a,1); // named versions
//   auto e=vmap(a,[](int x){return 1ll*x*x;}); // vector<ll>
//   auto f=vfilter(a,[](int x){return x&1;}); // keep odd elements
//   auto g=vunique(b); // sort and remove duplicates
//   a&b / a|b: set intersection / union; a&=b / a|=b modify a.
//   vintersection(a,b); vunion(a,b); // sorted and unique results
//   vsort(a); vreverse(a); vslice(a,1,3); // return new vectors
//   Compound assignments and ++/-- modify a; other operations return copies.
//   O(n), except vsort / vunique: O(n log n); vslice: O(r-l+1).
//   Set intersection / union: O(n log n+m log m).
//
////////////////////////////////////////////////////////////////
template<typename T=int>
vector<T> vrange(T l,T r){
	static_assert(is_integral_v<T>&&sizeof(T)<=8&&!is_same_v<T,bool>);
	assert(l<=r);
	__int128 n=(__int128)r-l+1;
	assert(n<=vector<T>().max_size());
	vector<T> a((size_t)n);
	a[0]=l;
	for(size_t i=1;i<a.size();i++)a[i]=a[i-1]+1;
	return a;
}
template<typename T>
vector<T> vcat(vector<T> a,const vector<T> &b){
	assert(b.size()<=a.max_size()-a.size());
	a.reserve(a.size()+b.size());
	a.insert(a.end(),b.begin(),b.end());
	return a;
}
template<typename T,typename F,typename Alloc>
auto vmap(const vector<T,Alloc> &a,F f){
	using U=decay_t<decltype(f((const T&)a[0]))>;
	vector<U> b;
	b.reserve(a.size());
	for(const T &x:a)b.push_back(f(x));
	return b;
}
template<typename T,typename U>
auto vadd(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y+x;});
}
template<typename T,typename U>
auto vsub(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y-x;});
}
template<typename T,typename U>
auto operator+(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y+x;});
}
template<typename T,typename U>
auto operator-(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y-x;});
}
template<typename T,typename U>
auto operator*(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y*x;});
}
template<typename T,typename U>
auto operator/(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y/x;});
}
template<typename T,typename U>
auto operator%(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y%x;});
}
template<typename T,typename U>
auto operator&(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y&x;});
}
template<typename T,typename U>
auto operator|(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y|x;});
}
template<typename T,typename U>
auto operator^(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y^x;});
}
template<typename T,typename U>
auto operator<<(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y<<x;});
}
template<typename T,typename U>
auto operator>>(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y>>x;});
}
// Deducing Alloc preserves the STL vector-vector comparison overloads.
template<typename T,typename U,typename Alloc>
auto operator==(const vector<T,Alloc> &a,const U &x){
	return vmap(a,[&](const T &y){return y==x;});
}
template<typename T,typename U,typename Alloc>
auto operator!=(const vector<T,Alloc> &a,const U &x){
	return vmap(a,[&](const T &y){return y!=x;});
}
template<typename T,typename U,typename Alloc>
auto operator<(const vector<T,Alloc> &a,const U &x){
	return vmap(a,[&](const T &y){return y<x;});
}
template<typename T,typename U,typename Alloc>
auto operator<=(const vector<T,Alloc> &a,const U &x){
	return vmap(a,[&](const T &y){return y<=x;});
}
template<typename T,typename U,typename Alloc>
auto operator>(const vector<T,Alloc> &a,const U &x){
	return vmap(a,[&](const T &y){return y>x;});
}
template<typename T,typename U,typename Alloc>
auto operator>=(const vector<T,Alloc> &a,const U &x){
	return vmap(a,[&](const T &y){return y>=x;});
}
template<typename T,typename U>
auto operator&&(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y&&x;});
}
template<typename T,typename U>
auto operator||(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y||x;});
}
// Compound assignments modify a; copy x in case it refers to a[i].
template<typename T,typename U>
vector<T>& operator+=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y+=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator-=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y-=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator*=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y*=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator/=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y/=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator%=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y%=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator&=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y&=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator|=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y|=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator^=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y^=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator<<=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y<<=x;a[i]=y;}
	return a;
}
template<typename T,typename U>
vector<T>& operator>>=(vector<T> &a,U x){
	for(size_t i=0;i<a.size();i++){T y=a[i];y>>=x;a[i]=y;}
	return a;
}
template<typename T>
auto operator+(const vector<T> &a){
	return vmap(a,[](const T &x){return +x;});
}
template<typename T>
auto operator-(const vector<T> &a){
	return vmap(a,[](const T &x){return -x;});
}
template<typename T>
auto operator~(const vector<T> &a){
	return vmap(a,[](const T &x){return ~x;});
}
template<typename T>
auto operator!(const vector<T> &a){
	return vmap(a,[](const T &x){return !x;});
}
template<typename T>
vector<T>& operator++(vector<T> &a){
	for(T &x:a)++x;
	return a;
}
template<typename T>
vector<T> operator++(vector<T> &a,int){
	vector<T> b=a;++a;
	return b;
}
template<typename T>
vector<T>& operator--(vector<T> &a){
	for(T &x:a)--x;
	return a;
}
template<typename T>
vector<T> operator--(vector<T> &a,int){
	vector<T> b=a;--a;
	return b;
}
template<typename T,typename F>
vector<T> vfilter(const vector<T> &a,F f){
	vector<T> b;
	b.reserve(a.size());
	for(const T &x:a)if(f(x))b.push_back(x);
	return b;
}
template<typename T>
vector<T> vsort(vector<T> a){
	sort(a.begin(),a.end());
	return a;
}
template<typename T>
vector<T> vunique(vector<T> a){
	sort(a.begin(),a.end());
	a.erase(unique(a.begin(),a.end()),a.end());
	return a;
}
// Set operations: input may be unsorted and contain duplicates.
// O(n log n+m log m); the result is sorted and unique.
template<typename T>
vector<T> vintersection(vector<T> a,vector<T> b){
	a=vunique(move(a));b=vunique(move(b));
	vector<T> c;
	c.reserve(min(a.size(),b.size()));
	set_intersection(a.begin(),a.end(),b.begin(),b.end(),back_inserter(c));
	return c;
}
template<typename T>
vector<T> vunion(vector<T> a,vector<T> b){
	a=vunique(move(a));b=vunique(move(b));
	vector<T> c;
	assert(b.size()<=c.max_size()-a.size());
	c.reserve(a.size()+b.size());
	set_union(a.begin(),a.end(),b.begin(),b.end(),back_inserter(c));
	return c;
}
template<typename T>
vector<T> operator&(const vector<T> &a,const vector<T> &b){
	return vintersection(a,b);
}
template<typename T>
vector<T> operator|(const vector<T> &a,const vector<T> &b){
	return vunion(a,b);
}
template<typename T>
vector<T>& operator&=(vector<T> &a,const vector<T> &b){
	return a=vintersection(a,b);
}
template<typename T>
vector<T>& operator|=(vector<T> &a,const vector<T> &b){
	return a=vunion(a,b);
}
template<typename T>
vector<T> vreverse(vector<T> a){
	reverse(a.begin(),a.end());
	return a;
}
template<typename T>
vector<T> vslice(const vector<T> &a,size_t l,size_t r){
	assert(l<=r&&r<a.size());
	return vector<T>(a.begin()+l,a.begin()+r+1);
}
// end for basic/vector.cpp
/////////////////////////
// !!!!! GNU C++17; ranges are [l,r]; vslice uses 0-based indices.
// vunique sorts first. Binary operators take a scalar on the right, except
// vector & vector / vector | vector, which are set intersection / union.
// Vector-vector comparisons keep the STL lexicographical / equality rules.
// ^ is xor; !a is elementwise NOT. && / || do not short-circuit the operands.
// Result types follow element expressions; compound assignments keep T.
// Elements must support the operation. Avoid overflow, zero divisors and
// invalid shifts, just as with scalar arithmetic. !!!!
