////////////////////////////////////////////////////////////////
//
// template for vector operations
// usage: auto a=vrange(1,6); auto b=a+2; auto c=vmap(a,f);
//   vcat(a,b); vunique(a); a&b; a|b;
//
////////////////////////////////////////////////////////////////
// Integers in [l,r], inclusive; O(r-l+1).
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
// Concatenate two vectors; O(a.size()+b.size()).
template<typename T>
vector<T> vcat(vector<T> a,const vector<T> &b){
	assert(b.size()<=a.max_size()-a.size());
	a.reserve(a.size()+b.size());
	a.insert(a.end(),b.begin(),b.end());
	return a;
}
// Return f(x) for each element; result type is inferred. O(n).
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
// Alloc keeps the STL vector-vector comparisons more specific.
// Scalar-right arithmetic / bitwise operations return a new vector.
// Comparisons and && / || return vector<bool>; no operand short-circuit.
// Vector-vector comparisons keep STL rules; ^ is bitwise xor. O(n).
#define VEC_BINARY(op) \
template<typename T,typename U,typename Alloc> \
auto operator op(const vector<T,Alloc> &a,const U &x){ \
	return vmap(a,[&](const T &y){return y op x;}); \
}
VEC_BINARY(+) VEC_BINARY(-) VEC_BINARY(*) VEC_BINARY(/) VEC_BINARY(%)
VEC_BINARY(&) VEC_BINARY(|) VEC_BINARY(^) VEC_BINARY(<<) VEC_BINARY(>>)
VEC_BINARY(==) VEC_BINARY(!=) VEC_BINARY(<) VEC_BINARY(<=)
VEC_BINARY(>) VEC_BINARY(>=) VEC_BINARY(&&) VEC_BINARY(||)
#undef VEC_BINARY

// Copy x in case it refers to a[i]; the local T also supports vector<bool>.
// Compound assignments modify a and keep the element type. O(n).
#define VEC_ASSIGN(op) \
template<typename T,typename U> \
vector<T>& operator op(vector<T> &a,U x){ \
	for(size_t i=0;i<a.size();i++){T y=a[i];y op x;a[i]=y;} \
	return a; \
}
VEC_ASSIGN(+=) VEC_ASSIGN(-=) VEC_ASSIGN(*=) VEC_ASSIGN(/=) VEC_ASSIGN(%=)
VEC_ASSIGN(&=) VEC_ASSIGN(|=) VEC_ASSIGN(^=) VEC_ASSIGN(<<=) VEC_ASSIGN(>>=)
#undef VEC_ASSIGN

// Elementwise +, -, ~ and !; ! returns vector<bool>. O(n).
#define VEC_UNARY(op) \
template<typename T> \
auto operator op(const vector<T> &a){ \
	return vmap(a,[](const T &x){return op x;}); \
}
VEC_UNARY(+) VEC_UNARY(-) VEC_UNARY(~) VEC_UNARY(!)
#undef VEC_UNARY

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
// Keep elements satisfying f(x), preserving their order. O(n).
template<typename T,typename F>
vector<T> vfilter(const vector<T> &a,F f){
	vector<T> b;
	b.reserve(a.size());
	for(const T &x:a)if(f(x))b.push_back(x);
	return b;
}
// Return a sorted copy; O(n log n).
template<typename T>
vector<T> vsort(vector<T> a){
	sort(a.begin(),a.end());
	return a;
}
// Return a sorted copy with duplicates removed; O(n log n).
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
// Return a reversed copy; O(n).
template<typename T>
vector<T> vreverse(vector<T> a){
	reverse(a.begin(),a.end());
	return a;
}
// Copy 0-based indices [l,r], inclusive; O(r-l+1).
template<typename T>
vector<T> vslice(const vector<T> &a,size_t l,size_t r){
	assert(l<=r&&r<a.size());
	return vector<T>(a.begin()+l,a.begin()+r+1);
}
// end for basic/vector.cpp
/////////////////////////
