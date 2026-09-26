////////////////////////////////////////////////////////////////
//
// template for inv, CRT, etc.
// version 1.0 (Last Update Jun 23rd, 2026)
//
// usage:
//    inv(x,M); COE a(r,M); COE c=a+b;
//   c.empty() if no solution
//
////////////////////////////////////////////////////////////////
using ull=unsigned long long;
inline ull gcd(ull x,ull y){
	// Stein's Algorithm
	if(!x||!y) return x|y;
	int c=__builtin_ctzll(x|y);
	x>>=__builtin_ctzll(x),y>>=__builtin_ctzll(y);
	while(x!=y){
		if(x<y)swap(x,y);
		x-=y;
		x>>=__builtin_ctzll(x);
	}
	return y<<c;
}
inline ull lcm(ull x,ull y){
	if(!x||!y)return 0;
	x/=gcd(x,y);
	assert(x<=ULLONG_MAX/y);
	return x*y;
}
inline ll inv(ll x,ll M){
	// return x^{-1} mod M
	assert(M>0);
	x%=M;
	if(x<0)x+=M;
	ull a=M,b=x;__int128 y=0,z=1;
	while(b){
		const ull q=a/b,c=a-q*b;
		a=b,b=c;
		const __int128 w=y-static_cast<__int128>(q)*z;
		y=z,z=w;
	}
	assert(a==1U);y%=M;
	return (y>=0?y:M+y);
}
////////////////////////////////
// Congruence equation, make, merge, etc.
struct COE{
ll r,M;
COE():r(0),M(1){}
COE(ll _r,ll _M):r(_r),M(_M){assert(M>0);r%=M;if(r<0)r+=M;}
inline bool empty()const{return M==-1;}
inline bool accept(ll x)const{
	if(empty())return false;
	x%=M;
	return (x<0?x+M:x)==r;
}
inline COE& operator +=(const COE &x){
	if(empty()||x.empty()){M=-1;return *this;}
	ll g=gcd(M,x.M);
	if(r%g!=x.r%g){M=-1;return *this;}
	assert(M/g<=LLONG_MAX/x.M);
	ll i=inv(M/g,x.M/g),new_M=M/g*x.M;
	r=((__int128)r+static_cast<__int128>(x.r-r)/g*i%(x.M/g)*M)%new_M;
	M=new_M;
	if(r<0) r+=M;
	return *this;
}
inline COE operator +(const COE &x)const{
	return COE(*this)+=x;
}
};
// end for number-theory/modular-arithmetic.cpp
/////////////////////////
// !!!!! CRT results may use __int128; basic/fast-io.cpp provides decimal IO if
// needed. !!!!
