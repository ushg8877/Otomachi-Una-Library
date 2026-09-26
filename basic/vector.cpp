////////////////////////////////////////////////////////////////
//
// template for vector operations
// usage:
//   auto a=vrange(1,6); // {1,2,3,4,5,6}, [l,r]
//   auto b=vcat(a,a); // concatenate
//   auto c=a+2; auto d=a-1; auto h=a^3; // element + / - / bitwise xor
//   vadd(a,2); vsub(a,1); // named versions
//   auto e=vmap(a,[](int x){return 1ll*x*x;}); // vector<ll>
//   auto f=vfilter(a,[](int x){return x&1;}); // keep odd elements
//   auto g=vunique(b); // sort and remove duplicates
//   vsort(a); vreverse(a); vslice(a,1,3); // return new vectors
//   All operations leave the input unchanged; assign the result if needed.
//   O(n), except vsort / vunique: O(n log n); vslice: O(r-l+1).
//
////////////////////////////////////////////////////////////////
template<typename T=int> vector<T> vrange(T l,T r){
	static_assert(is_integral_v<T>&&sizeof(T)<=8&&!is_same_v<T,bool>);
	assert(l<=r);__int128 n=(__int128)r-l+1;
	assert(n<=vector<T>().max_size());
	vector<T> a((size_t)n);a[0]=l;
	for(size_t i=1;i<a.size();i++)a[i]=a[i-1]+1;return a;
}
template<typename T> vector<T> vcat(vector<T> a,const vector<T> &b){
	assert(b.size()<=a.max_size()-a.size());a.reserve(a.size()+b.size());
	a.insert(a.end(),b.begin(),b.end());return a;
}
template<typename T,typename F> auto vmap(const vector<T> &a,F f){
	using U=decay_t<invoke_result_t<F&,const T&>>;
	vector<U> b;b.reserve(a.size());
	for(const T &x:a)b.push_back(f(x));return b;
}
template<typename T,typename U> auto vadd(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y+x;});
}
template<typename T,typename U> auto vsub(const vector<T> &a,const U &x){
	return vmap(a,[&](const T &y){return y-x;});
}
template<typename T,typename U,enable_if_t<is_arithmetic_v<U>,int> =0>
auto operator+(const vector<T> &a,const U &x)->vector<decay_t<decltype(declval<const T&>()+x)>>{
	return vmap(a,[&](const T &y){return y+x;});
}
template<typename T,typename U,enable_if_t<is_arithmetic_v<U>,int> =0>
auto operator-(const vector<T> &a,const U &x)->vector<decay_t<decltype(declval<const T&>()-x)>>{
	return vmap(a,[&](const T &y){return y-x;});
}
template<typename T,typename U,enable_if_t<is_arithmetic_v<U>,int> =0>
auto operator^(const vector<T> &a,const U &x)->vector<decay_t<decltype(declval<const T&>()^x)>>{
	return vmap(a,[&](const T &y){return y^x;});
}
template<typename T,typename F> vector<T> vfilter(const vector<T> &a,F f){
	vector<T> b;b.reserve(a.size());for(const T &x:a)if(f(x))b.push_back(x);return b;
}
template<typename T> vector<T> vsort(vector<T> a){sort(a.begin(),a.end());return a;}
template<typename T> vector<T> vunique(vector<T> a){
	sort(a.begin(),a.end());a.erase(unique(a.begin(),a.end()),a.end());return a;
}
template<typename T> vector<T> vreverse(vector<T> a){reverse(a.begin(),a.end());return a;}
template<typename T> vector<T> vslice(const vector<T> &a,size_t l,size_t r){
	assert(l<=r&&r<a.size());return vector<T>(a.begin()+l,a.begin()+r+1);
}
// end for basic/vector.cpp
/////////////////////////
// !!!!! GNU C++17; ranges are [l,r]; vslice uses 0-based indices. vunique sorts first. Operators take a scalar on the right; ^ means bitwise xor. vmap / scalar operations infer the result type from the expression; arithmetic must not overflow. !!!!
