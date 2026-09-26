////////////////////////////////////////////////////////////////
//
// template for dynamic bitset
// usage:
//   Bitset a(n),b(n); a.set(i); a.set(i,false); a.flip(i); a[i];
//   a.setN(n); // resize and clear; set() clears, fill(true) sets all bits
//   a&b; a|b; a^b; ~a; a<<k; a>>k; // compound assignments supported
//   a.count(); a.count(l,r); // count 1s in [l,r], inclusive
//   auto v=a.positions(); // increasing indices of all 1s
//   a.any(); a.none(); a.all();
//   a.first(); a.next(i); // first / strictly next 1; n if absent
//   a.last(); a.prev(i); // last / strictly previous 1; -1 if absent
//   a.window_xor(l,r,b,to); // b[to+i] ^= a[l+i], 0<=i<=r-l
//   window_and / window_or are similar; overlapping self updates are safe.
//   String / vector constructors: input[i] is bit i; str() uses this order.
//   Point operations O(1); bitwise / count / find O(ceil(n/64)) worst case.
//   count(l,r): O((r-l+1)/64+1); positions(): O(ceil(n/64)+count()).
//   String / vector conversions O(n).
//   Window operations O((r-l+1)/64+1), O(1) extra space.
//
////////////////////////////////////////////////////////////////
struct Bitset{
using ull=unsigned long long;
int n=0;
vector<ull> a;
Bitset()=default;
explicit Bitset(int n){setN(n);}
template<typename T>
explicit Bitset(const vector<T> &v){
	assert(v.size()<=INT_MAX);setN(v.size());
	for(int i=0;i<n;i++){assert(v[i]==0||v[i]==1);set(i,v[i]);}
}
explicit Bitset(const string &s){
	assert(s.size()<=INT_MAX);setN(s.size());
	for(int i=0;i<n;i++){
		assert(s[i]=='0'||s[i]=='1');set(i,s[i]=='1');
	}
}
void setN(int _n){assert(_n>=0);n=_n;a.assign((n+63ull)/64,0);}
void set(){fill(false);}
void fill(bool v){std::fill(a.begin(),a.end(),v?~0ull:0);trim();}
int size()const{return n;}
void set(int i,bool v=true){
	assert(0<=i&&i<n);
	ull m=1ull<<(i&63);
	if(v)a[i>>6]|=m;else a[i>>6]&=~m;
}
void reset(int i){set(i,false);}
void flip(int i){assert(0<=i&&i<n);a[i>>6]^=1ull<<(i&63);}
void flip(){for(ull &x:a)x=~x;trim();}
bool get(int i)const{assert(0<=i&&i<n);return a[i>>6]>>(i&63)&1;}
bool operator[](int i)const{return get(i);}
int count()const{
	int ans=0;for(ull x:a)ans+=__builtin_popcountll(x);
	return ans;
}
int count(int l,int r)const{
	assert(0<=l&&l<=r&&r<n);
	int L=l>>6,R=r>>6;
	ull x=~0ull<<(l&63),y=mask((r&63)+1);
	if(L==R)return __builtin_popcountll(a[L]&x&y);
	int ans=__builtin_popcountll(a[L]&x)+__builtin_popcountll(a[R]&y);
	for(int i=L+1;i<R;i++)ans+=__builtin_popcountll(a[i]);
	return ans;
}
vector<int> positions()const{
	vector<int> v;v.reserve(count());
	for(int i=0;i<(int)a.size();i++){
		for(ull x=a[i];x;x&=x-1)v.push_back(i*64+__builtin_ctzll(x));
	}
	return v;
}
bool any()const{for(ull x:a)if(x)return true;return false;}
bool none()const{return !any();}
bool all()const{return count()==n;}
int next(int p)const{
	assert(-1<=p&&p<=n);
	if(p>=n-1)return n;
	int i=(p+1)>>6;
	ull x=a[i]&(~0ull<<((p+1)&63));
	while(true){
		if(x)return i*64+__builtin_ctzll(x);
		if(++i==(int)a.size())return n;
		x=a[i];
	}
}
int prev(int p)const{
	assert(0<=p&&p<=n);
	if(!p)return -1;
	int i=(p-1)>>6;
	ull x=a[i]&mask(((p-1)&63)+1);
	while(true){
		if(x)return i*64+63-__builtin_clzll(x);
		if(--i<0)return -1;
		x=a[i];
	}
}
int first()const{return next(-1);}
int last()const{return prev(n);}
int _Find_first()const{return first();}
int _Find_next(int p)const{return next(p);}
#define BITSET_OP(op) \
Bitset& operator op##=(const Bitset &b){ \
	assert(n==b.n); \
	for(size_t i=0;i<a.size();i++)a[i] op##= b.a[i]; \
	return *this; \
} \
friend Bitset operator op(Bitset a,const Bitset &b){return a op##= b;}
BITSET_OP(&) BITSET_OP(|) BITSET_OP(^)
#undef BITSET_OP
Bitset operator~()const{Bitset b=*this;b.flip();return b;}
bool operator==(const Bitset &b)const{return n==b.n&&a==b.a;}
bool operator!=(const Bitset &b)const{return !(*this==b);}
// Compare as nonnegative binary integers; bit 0 is the least significant.
bool operator<(const Bitset &b)const{
	assert(n==b.n);
	for(int i=(int)a.size()-1;i>=0;i--)if(a[i]!=b.a[i])return a[i]<b.a[i];
	return false;
}
Bitset& operator<<=(long long k){
	assert(k>=0);
	if(k>=n){set();return *this;}
	int q=k>>6,s=k&63;
	for(int i=(int)a.size()-1;i>=q;i--){
		ull x=a[i-q]<<s;
		if(s&&i>q)x|=a[i-q-1]>>(64-s);
		a[i]=x;
	}
	std::fill(a.begin(),a.begin()+q,0);trim();
	return *this;
}
Bitset& operator>>=(long long k){
	assert(k>=0);
	if(k>=n){set();return *this;}
	int q=k>>6,s=k&63,m=a.size();
	for(int i=0;i<m-q;i++){
		ull x=a[i+q]>>s;
		if(s&&i+q+1<m)x|=a[i+q+1]<<(64-s);
		a[i]=x;
	}
	std::fill(a.begin()+m-q,a.end(),0);
	return *this;
}
friend Bitset operator<<(Bitset a,long long k){return a<<=k;}
friend Bitset operator>>(Bitset a,long long k){return a>>=k;}
void window_and(int l,int r,Bitset &b,int to)const{window<'&'>(l,r,b,to);}
void window_or(int l,int r,Bitset &b,int to)const{window<'|'>(l,r,b,to);}
void window_xor(int l,int r,Bitset &b,int to)const{window<'^'>(l,r,b,to);}
string str()const{
	string s(n,'0');for(int i=0;i<n;i++)s[i]+=get(i);
	return s;
}
private:
static ull mask(int k){return k==64?~0ull:(1ull<<k)-1;}
void trim(){if(n&63)a.back()&=mask(n&63);}
// Read k bits starting at p, packed into the low bits of a word.
ull read(int p,int k)const{
	int i=p>>6,s=p&63;
	ull x=a[i]>>s;
	if(s&&k>64-s)x|=a[i+1]<<(64-s);
	return x&mask(k);
}
template<char OP>
void window(int l,int r,Bitset &b,int to)const{
	assert(0<=l&&l<=r&&r<n);
	int len=r-l+1;
	assert(0<=to&&len<=b.n&&to<=b.n-len);
	auto apply=[&](int p,int k){
		int i=(to+p)>>6,s=(to+p)&63;
		ull x=read(l+p,k)<<s,m=mask(k)<<s;
		if constexpr(OP=='&')b.a[i]&=x|~m;
		if constexpr(OP=='|')b.a[i]|=x;
		if constexpr(OP=='^')b.a[i]^=x;
	};
	if(this==&b&&to>l){
		for(int p=len;p;){
			int k=min(p,((to+p-1)&63)+1);
			p-=k;apply(p,k);
		}
	}else{
		for(int p=0;p<len;){
			int k=min(len-p,64-((to+p)&63));
			apply(p,k);p+=k;
		}
	}
}
};
using BS=Bitset;
// end for data-structure/bitset.cpp
/////////////////////////
// !!!!! Indices are 0-based; windows use [l,r], inclusive. Binary bitwise
// operations require equal sizes. Shifts discard bits outside [0,n).
// Self windows use the original source bits. Unused high bits stay zero;
// Do not modify n/a directly. str() prints bit 0 first. !!!!
