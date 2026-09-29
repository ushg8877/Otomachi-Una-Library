#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using ull=unsigned long long;
// head file
#ifndef DEBUG
	#define cerr for(;false;) cerr
#endif
template<typename T>
inline constexpr T inf=T{};
template<>
inline constexpr int inf<int> = 1010000000;
template<>
inline constexpr long long inf<long long> = 2020000000000000000LL;
template<typename T,size_t S>
istream& operator >>(istream &o,array<T,S> &I){
	assert(S>0);
	for(int i=0;i<S;i++) o>>I[i];
	return o;
}
template<typename T,size_t S>
ostream& operator <<(ostream &o,const array<T,S> &I){
	o<<"[";
	for(int i=0;i<S;i++){
		if(i==50){o<<"...(total "<<S<<" elements)";break;}
		o<<I[i];
		if(i+1<I.size()) o<<",";
	}
	o<<"]";
	return o;
}
template<typename T>
istream& operator >>(istream &o,vector<T> &I){
	assert(!I.empty());
	for(int i=0;i<I.size();i++) o>>I[i];
	return o;
}
template<typename T>
ostream& operator <<(ostream &o,const vector<T> &I){
	o<<"[";
	for(int i=0;i<I.size();i++){
		if(i==50){o<<"...(total "<<I.size()<<" elements)";break;}
		o<<I[i];
		if(i+1<I.size()) o<<",";
	}
	o<<"]";
	return o;
}

#define DEBUG_1(x) #x"="<<(x)
#define DEBUG_2(x,y) DEBUG_1(x)<<", "<<DEBUG_1(y)
#define DEBUG_3(x,y,z) DEBUG_2(x,y)<<", "<<DEBUG_1(z)
#define DEBUG_4(x,y,z,w) DEBUG_3(x,y,z)<<", "<<DEBUG_1(w)
#define DEBUG_5(x,y,z,w,a) DEBUG_4(x,y,z,w)<<", "<<DEBUG_1(a)
#define DEBUG_6(x,y,z,w,a,b) DEBUG_5(x,y,z,w,a)<<", "<<DEBUG_1(b)
#define DEBUG_7(x,y,z,w,a,b,c) DEBUG_6(x,y,z,w,a,b)<<", "<<DEBUG_1(c)
#define DEBUG_8(x,y,z,w,a,b,c,d) DEBUG_7(x,y,z,w,a,b,c)<<", "<<DEBUG_1(d)
#define DEBUG_9(x,y,z,w,a,b,c,d,e) DEBUG_8(x,y,z,w,a,b,c,d)<<", "<<DEBUG_1(e)
#define DEBUG_10(x,y,z,w,a,b,c,d,e,f) DEBUG_9(x,y,z,w,a,b,c,d,e)<<", "<<DEBUG_1(f)
#define GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,NAME,...) NAME
#ifdef DEBUG
#define debug(...) do { \
	cerr<<"(line "<< __LINE__ <<"): "; \
	cerr<<GET_MACRO(__VA_ARGS__,DEBUG_10,DEBUG_9,DEBUG_8, \
						   DEBUG_7,DEBUG_6,DEBUG_5,DEBUG_4, \
						   DEBUG_3,DEBUG_2,DEBUG_1)(__VA_ARGS__); \
	cerr<<endl; \
} while(0)
#define debugl(s,...) do { \
	cerr<<"(line "<< __LINE__ <<"): <"<<s<<"> "; \
	cerr<<GET_MACRO(__VA_ARGS__,DEBUG_10,DEBUG_9,DEBUG_8, \
						   DEBUG_7,DEBUG_6,DEBUG_5,DEBUG_4, \
						   DEBUG_3,DEBUG_2,DEBUG_1)(__VA_ARGS__); \
	cerr<<endl; \
} while(0)
#else
#define debug(...) ((void)0)
#define debugl(...) ((void)0)
#endif
inline int last_bit(ull x){assert(x);return 63-__builtin_clzll(x);}
inline int first_bit(ull x){assert(x);return __builtin_ctzll(x);}
inline int popc(ull x){return __builtin_popcountll(x);}
template<typename A,typename B> inline bool chkmax(A &x,const B &y){if(x<y){
	x=y;return true;}return false;}
template<typename A,typename B> inline bool chkmin(A &x,const B &y){if(x>y){
	x=y;return true;}return false;}
int main(){
	#ifdef LOCAL
		if(!freopen("Otomachi_Una.in","r",stdin)) return 1;
		if(!freopen("Otomachi_Una.out","w",stdout)) return 1;
	#endif
	ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
	return 0;
}
