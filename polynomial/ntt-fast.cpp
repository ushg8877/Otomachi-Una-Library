////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/polynomial/ntt-fast.cpp
//
// usage: poly h=f*g; Inv(f); Ln(f); Exp(f);
//   Sqrt(f); Div(f,g,q,r); Eval(f,x); FSPE(f,g,k);
//
////////////////////////////////////////////////////////////////
template<unsigned _M>
struct ModInt{
static constexpr unsigned MOD=_M;
static_assert(1<MOD&&MOD<=INT_MAX);
unsigned x;
constexpr ModInt():x(0){}
constexpr ModInt(unsigned y):x(y%MOD){}
constexpr ModInt(int y):x((y%=static_cast<int>(MOD))<0?y+MOD:y){}
constexpr ModInt(unsigned long long y):x(y%MOD){}
constexpr ModInt(long long y):x((y%=static_cast<long long>(MOD))<0?y+MOD:y){}
ModInt &operator+=(const ModInt &a){x=(((x+=a.x)>=MOD)?x-MOD:x);return *this;}
ModInt &operator-=(const ModInt &a){x=(((x-=a.x)>=MOD)?x+MOD:x);return *this;}
ModInt &operator*=(const ModInt &a){
	x=static_cast<unsigned long long>(x)*a.x%MOD;
	return *this;
}
ModInt &operator/=(const ModInt &a){return (*this)*=a.inv();}
bool operator==(const ModInt &a)const{return x==a.x;}
bool operator!=(const ModInt &a)const{return x!=a.x;}
explicit operator bool()const{return x!=0;}
bool operator!()const{return x==0;}
ModInt inv()const{
	unsigned a=MOD,b=x;ll y=0,z=1;
	while(b){
		const unsigned q=a/b,c=a-q*b;
		a=b,b=c;
		const ll w=y-(ll)q*z;
		y=z,z=w;
	}
	assert(a==1);
	return ModInt(y);
}
ModInt pow(long long k)const{
	unsigned long long e=k;
	ModInt a=*this,b=1;
	if(k<0){a=a.inv();e=0-e;}
	for(;e;e>>=1){if(e&1)b*=a;a*=a;}
	return b;
}
ModInt &operator++(){*this+=1;return *this;}
ModInt &operator--(){*this-=1;return *this;}
ModInt operator+()const{return *this;}
ModInt operator-()const{ModInt a;a.x=(x?MOD-x:0u);return a;}
ModInt operator+(const ModInt &a)const{return ModInt(*this)+=a;}
ModInt operator-(const ModInt &a)const{return ModInt(*this)-=a;}
ModInt operator*(const ModInt &a)const{return ModInt(*this)*=a;}
ModInt operator/(const ModInt &a)const{return ModInt(*this)/=a;}

template<typename T>
ModInt friend operator+(T a,const ModInt &b){return ModInt(a)+=b;}
template<typename T>
ModInt friend operator-(T a,const ModInt &b){return ModInt(a)-=b;}
template<typename T>
ModInt friend operator*(T a,const ModInt &b){return ModInt(a)*=b;}
template<typename T>
ModInt friend operator/(T a,const ModInt &b){return ModInt(a)/=b;}
friend istream&operator>>(istream&i,ModInt &x){
	ll y;
	if(i>>y)x=ModInt(y);
	return i;
}
friend ostream&operator<<(ostream&o,const ModInt &x){o<<x.x;return o;}
};
//////////////////////////////////////////////////////////////////

const int MOD=998244353;
using mint=ModInt<MOD>;
vector<mint>fac{1},ifac{1},inv{0,1};
void init(int n=0){
	assert(0<=n&&n<MOD);
	int m=fac.size();if(n<m)return;
	n=max(n,(int)min(2ll*m,(ll)MOD-1));
	fac.resize(n+1);
	ifac.resize(n+1);
	inv.resize(max(n+1,2));
	for(int i=m;i<=n;i++){
		fac[i]=fac[i-1]*i;
		if(i>1)inv[i]=-(MOD/i)*inv[MOD%i];
		ifac[i]=ifac[i-1]*inv[i];
	}
}
inline mint C(int x,int y){
	// choose x from y
	if(x<0||y<x)return 0;
	init(y);
	return fac[y]*ifac[x]*ifac[y-x];
}
inline mint binom(int y,int x){return C(x,y);}


using poly=vector<mint>;
using Poly=poly;

constexpr unsigned MO=998244353U;
constexpr unsigned MO2=2U*MO;
constexpr int FFT_MAX=23;
// AVX2 Montgomery / block FPS: QedDust413 & Killer_joke.
#pragma GCC push_options
#pragma GCC target("avx2")
#pragma GCC optimize("O3","unroll-loops")
namespace Poly_Fast{
using i64=int64_t;
using u32=uint32_t;
using u64=uint64_t;
using idt=int;
constexpr u32 M=998244353;
typedef unsigned u32x8 __attribute__((vector_size(32),may_alias));
typedef ull u64x4 __attribute__((vector_size(32),may_alias));
using I256=__m256i;
inline u64x4 fus_mul(u32x8 x,u32x8 y){
	return (u64x4)_mm256_mul_epu32((I256)x, (I256)y);
}
inline u32x8 swaplohi128(u32x8 x){
	return (u32x8)_mm256_permute2x128_si256((I256)x,(I256)x,1);
}
template<int typ>inline u32x8 shuffle(u32x8 x){
	return (u32x8)_mm256_shuffle_epi32((I256)x,typ);
}
template<int typ>inline u32x8 blend(u32x8 x,u32x8 y){
	return (u32x8)_mm256_blend_epi32((I256)x,(I256)y,typ);
}
inline u32x8&x8(u32*p){return*((u32x8*)p);}
inline const u32x8&x8(const u32*p){return*((u32x8*)p);}
inline u32x8 min_u32(u32x8 x,u32x8 y){
	return (u32x8)_mm256_min_epu32((I256)x, (I256)y);
}
constexpr u32x8 padd(u32 x){return (u32x8){x,x,x,x,x,x,x,x};}
inline u32x8 loadu_u32x8(const u32*f){
	return (u32x8)_mm256_loadu_si256((const __m256i_u*)f);
}
inline void storeu(void*f,u32x8 x){_mm256_storeu_si256((__m256i_u*)f,(I256)x);}

constexpr u32 get_nr(u32 M){
	u32 Iv=1;
	for(int i=0;i<5;++i){
		Iv*=2-M*Iv;
	}
	return Iv;
}
constexpr idt bcl(idt x){return ((x<2)?1:idt(2)<<__lg(x-1));}

constexpr u32 R=(-M)%M,E={},nR=M-R,M2=M*2,niv=-get_nr(M),R2=(-u64(M))%M;
constexpr u32 shrk(u32 x){return min(x,x-M);}
constexpr u32 shrk2(u32 x){return min(x,x-M2);}
constexpr u32 dil2(u32 x){return min(x,x+M2);}
constexpr u32 reduce(u64 x){return (x+u64(u32(x)*niv)*M)>>32;}
constexpr u32 reduce_s(u64 x){return shrk(reduce(x));}

constexpr u32 add(u32 x,u32 y){return shrk2(x+y);}
constexpr u32 sub(u32 x,u32 y){return dil2(x-y);}
constexpr u32 mul(u32 x,u32 y){return reduce(u64(x)*y);}
constexpr u32 mul_s(u32 x,u32 y){return reduce_s(u64(x)*y);}
constexpr u32 qpw(u32 a,u32 b,u32 r=R){
	for(;b;b>>=1,a=mul(a,a))if(b&1)r=mul(r,a);
	return r;
}
constexpr u32 inv(u32 x){return qpw(x,M-2);}
constexpr u32 dvs(u32 x,u32 y){return qpw(y,M-2,x);}
constexpr u32 neg(u32 x){return M2-x;}
constexpr u32 in(u32 x){return mul(x,R2);}
constexpr u32 in_s(u32 x){return mul_s(x,R2);}
constexpr u32 out(u32 x){return reduce_s(x);}
constexpr bool equals(u32 x,u32 y){return out(x)==out(y);}
constexpr void clr(u32&x){x=E;}

constexpr auto Rx8=padd(R),Ex8=padd(E),Mx8=padd(M),M2x8=padd(M2),
	nivx8=padd(niv);
inline u32x8 shrk(u32x8 x){return min_u32(x,x-Mx8);}
inline u32x8 dil2(u32x8 x){return min_u32(x,x+M2x8);}
inline u32x8 shrk2(u32x8 x){return min_u32(x,x-M2x8);}
inline u32x8 add(u32x8 x,u32x8 y){return shrk2(x+y);}
inline u32x8 sub(u32x8 x,u32x8 y){return dil2(x-y);}
inline u32x8 reduce(u64x4 a,u64x4 b){
	auto c=fus_mul(u32x8(a),nivx8),d=fus_mul(u32x8(b),nivx8);
	c=fus_mul(u32x8(c),Mx8),d=fus_mul(u32x8(d),Mx8);
	return blend<0xaa>(u32x8((a+c)>>32),u32x8(b+d));
}
inline u32x8 mul(u32x8 x,u32x8 y){
	return reduce(fus_mul(x,y),fus_mul(u32x8(u64x4(x)>>32),
		u32x8(u64x4(y)>>32)));
}
inline u32x8 qpw(u32x8 a,u32 b,u32x8 r=Rx8){
	for(;b;a=mul(a,a),b/=2){
		if(b&1){
			r=mul(r,a);
		}
	}
	return r;
}
inline u32x8 inv(u32x8 x){return qpw(x,M-2);}
inline u32x8 dvs(u32x8 x,u32x8 y){return qpw(y,M-2,x);}
inline u32x8 mul_s(u32x8 x,u32x8 y){return shrk(mul(x,y));}
inline u32x8 neg(u32x8 x){return M2x8-x;}
inline void clr(u32x8&x){x=Ex8;}

constexpr u32 _Amul(u32 a,u32 b,u32 c){return mul(a+b,c);}
constexpr u32 _Smul(u32 a,u32 b,u32 c){return mul(a-b+M2,c);}
inline u32x8 _Amul(u32x8 a,u32x8 b,u32x8 c){return mul(a+b,c);}
inline u32x8 _Smul(u32x8 a,u32x8 b,u32x8 c){return mul(a-b+M2x8,c);}
template<int typ>inline u32x8 Neg(u32x8 x){return blend<typ>(x,M2x8-x);}
constexpr u32x8 powXx8(u32 X){
	auto X2=mul_s(X,X),X3=mul_s(X2,X),X4=mul_s(X2,X2),X5=mul_s(X4,X),
		X6=mul_s(X4,X2),X7=mul_s(X4,X3);
	return (u32x8){R,X,X2,X3,X4,X5,X6,X7};
}
constexpr u32 _ADmul(u32 a,u32 b,u32 c,u32 d){
	return reduce_s(u64(a)*b+u64(c)*d);
}
inline u32x8 _ADmul(u32x8 a,u32x8 b,u32x8 c,u32x8 d){
	return shrk(reduce(fus_mul(a,b)+fus_mul(c,d),fus_mul(u32x8(u64x4(a)>>32),
		u32x8(u64x4(b)>>32))+fus_mul(u32x8(u64x4(c)>>32),u32x8(u64x4(d)>>32))));
}

constexpr auto Half=shrk(inv(in(2))),nHalf=M-Half;
struct vec_v{
	u32x8 v;
	constexpr vec_v(u32 x):v{padd(x)}{}
	constexpr operator u32()const{return v[0];}
	constexpr operator u32x8()const{return v;}
};
template<class F,class Op> inline void vec_op(F f,idt n,Op op){
	idt i=0;
	for(;i+7<n;i+=8){
		op(x8(f+i));
	}
	for(;i<n;++i){
		op(f[i]);
	}
}
template<class F,class G,class Op> inline void vec_op(F f,G g,idt n,Op op){
	idt i=0;
	for(;i+7<n;i+=8){
		op(x8(f+i), x8(g+i));
	}
	for(;i<n;++i){
		op(f[i],g[i]);
	}
}
template<class F,class G,class H,class Op> inline void vec_op(F f,G g,H h,
	idt n,Op op){
	idt i=0;
	for(;i+7<n;i+=8){
		op(x8(f+i),x8(g+i),x8(h+i));
	}
	for(;i<n;++i){
		op(f[i],g[i],h[i]);
	}
}
template<class F,class G,class H,class O,class Op> inline void vec_op(F f,G g,
	H h,O o,idt n,Op op){
	idt i=0;
	for(;i+7<n;i+=8){
		op(x8(f+i),x8(g+i),x8(h+i),x8(o+i));
	}
	for(;i<n;++i){
		op(f[i],g[i],h[i], o[i]);
	}
}
constexpr auto _g=in(3);
constexpr auto lml=__builtin_ctz(M-1);
struct Root{
	u32 t[lml+1];
	constexpr Root(u32 G):t{}{
		t[lml]=qpw(G,M>>lml);
		for(int i=lml;i>0;--i){t[i-1]=mul_s(t[i],t[i]);}
	}
	constexpr u32 operator[](int i)const{return t[i];}
};
struct Info{
	u32 rt3[lml-2],rt3_I[lml-2];
	// The last unused rate update at length 2^23 also needs a slot.
	u32x8 rt4ix8[lml-2],rt4ix8_I[lml-2];
	constexpr Info(const Root&w,const Root&wI):rt3{},rt3_I{},rt4ix8{},
		rt4ix8_I{}{
		auto pr=R,pr_I=R;
		for(int i=0;i<lml-2;pr=mul(pr,wI[i+3]),pr_I=mul(pr_I,w[i+3]),++i){
			rt3[i]=mul_s(pr,w[i+3]),rt3_I[i]=mul_s(pr_I,wI[i+3]);
		}
		pr=R,pr_I=R;
		for(int i=0;i<lml-3;pr=mul(pr,wI[i+4]),pr_I=mul(pr_I,w[i+4]),++i){
			rt4ix8[i]=powXx8(mul_s(pr,w[i+4])),rt4ix8_I[i]=powXx8(mul_s(pr_I,
				wI[i+4]));
		}
	}
};
constexpr Root rt1={_g},rt1_I={inv(_g)};
constexpr Info iab4={rt1,rt1_I};
constexpr auto Img=rt1[2];
constexpr auto Imgx8=padd(Img);
template<bool strict=false>inline void dif_2(u32&x,u32&y){
	auto sum=add(x,y),diff=sub(x,y);
	x=sum,y=diff;
	if constexpr(strict){
		x=shrk(x),y=shrk(y);
	}
}
template<bool strict=false>inline void dif_4(u32&x,u32&y,u32&z,u32&w){
	auto a=sub(x,z),b=_Smul(y,w,Img);
	x=add(x,z),y=add(y,w),z=add(a,b),w=sub(a,b),a=add(x,y),b=sub(x,y),x=a,y=b;
	if constexpr(strict){
		x=shrk(x),y=shrk(y),z=shrk(z),w=shrk(w);
	}
}
template<bool strict=false>inline void vec_dif_base4(u32x8*f,idt n){
	auto L=n/2;
	if(__builtin_ctzll(n)&1){
		for(idt j=0;j<L;++j){
			auto x=f[j],y=f[j+L];
			f[j]=x+y,f[j+L]=x-y+M2x8;
		}
		L/=2;
	}
	L/=2;
	for(idt l=L*4,k;L;l=L,L/=4){
		auto r=R,r2=R,r3=nR;
		k=1;
		for(auto i=f;i!=(f+n);r=mul_s(r,iab4.rt3[__builtin_ctzll(k++)]),
			r2=mul_s(r,r),r3=mul_s(r2,neg(r)),i+=l){
			auto rx8=padd(r),r2x8=padd(r2),r3x8=padd(r3);
			for(auto F0=i,F1=F0+L,F2=F1+L,F3=F2+L;F3!=i+l;++F0,++F1,++F2,++F3){
				auto f0=shrk2(*F0),f1=mul(*F1,rx8),f2=mul(*F2,r2x8),
					f3=mul(*F3,r3x8);
				auto f1f3=_Amul(f1,f3,Imgx8),f02=add(f0,f2),f13=sub(f1,f3),
					f_02=sub(f0,f2);
				*F0=f02+f13,*F1=f02-f13+M2x8,*F2=f_02+f1f3,*F3=f_02-f1f3+M2x8;
			}
		}
	}
	constexpr u32x8 pr2={R,R,R,Img,R,R,R,Img},pr4={R,R,R,R,R,rt1[3],Img,
		mul_s(Img,rt1[3])};
	auto rx8=Rx8;
	for(idt i=0;i<n;++i){
		auto&fi=f[i];
		fi=mul(fi,rx8),rx8=mul_s(rx8,iab4.rt4ix8[__builtin_ctzll(~i)]);
		fi=_Amul(Neg<0xf0>(fi),swaplohi128(fi),pr4);
		fi=_Amul(Neg<0xcc>(fi),shuffle<0x4e>(fi),pr2);
		fi=sub(shuffle<0xb1>(fi),Neg<0x55>(fi));
		if constexpr(strict){
			fi=shrk(fi);
		}
	}
}
template<u32 fx>inline void dit_2(u32&x,u32&y){
	constexpr auto iv2=mul_s(inv(in(2)),fx);
	auto a=_Amul(x,y,iv2),b=_Smul(x,y,iv2);
	x=a,y=b;
}
template<u32 fx>inline void dit_4(u32&x,u32&y,u32&z,u32&w){
	constexpr auto iv4=mul_s(inv(in(4)),fx),Imgi4=mul_s(iv4,Img);
	auto a=_Amul(x,y,iv4),b=_Smul(x,y,iv4);
	x=a,y=b,a=_Amul(z,w,iv4),b=_Smul(w,z,Imgi4),z=sub(x,a),w=sub(y,b),x=add(x,
		a),y=add(y,b);
}
template<u32 fx>inline void vec_dit_base4(u32x8*f,idt n){
	constexpr auto nR2=in_s(nR),M8=(M-1)/8;
	constexpr u32x8 pr2={nR2,nR2,nR2,in(Img),nR2,nR2,nR2,in(Img)},pr4={fx,fx,
		fx,fx,fx,mul_s(fx,rt1_I[3]),mul_s(fx,rt1_I[2]),mul_s(fx,
		mul_s(rt1_I[2],rt1_I[3]))};
	auto rx8=padd(M8>>__builtin_ctzll(n));
	idt L=1;
	for(idt i=0;i<n;++i){
		auto&fi=f[i];
		fi=_Amul(Neg<0xaa>(fi),shuffle<0xb1>(fi),pr2);
		fi=_Amul(Neg<0xcc>(fi),shuffle<0x4e>(fi),pr4);
		fi=_Amul(Neg<0xf0>(fi),swaplohi128(fi),rx8);
		rx8=mul_s(rx8,iab4.rt4ix8_I[__builtin_ctzll(~i)]);
	}
	for(idt l=L*4,k;L<(n/2);L=l,l*=4){
		auto r=R,r2=R,r3=R;
		k=1;
		for(auto i=f;i!=(f+n);r=mul_s(r,iab4.rt3_I[__builtin_ctzll(k++)]),
			r2=mul_s(r,r),r3=mul_s(r2,r),i+=l){
			auto rx8=padd(r),r2x8=padd(r2),r3x8=padd(r3);
			for(auto F0=i,F1=F0+L,F2=F1+L,F3=F2+L;F3!=i+l;++F0,++F1,++F2,++F3){
				auto f0=*F0,f1=*F1,f2=neg(*F2),f3=*F3;
				auto f2f3=_Amul(f3,f2,Imgx8),f01=add(f0,f1),f23=sub(f2,f3),
					f_01=sub(f0,f1);
				*F0=sub(f01,f23),*F1=_Amul(f_01,f2f3,rx8),*F2=_Amul(f01,f23,
					r2x8),*F3=_Smul(f_01,f2f3,r3x8);
			}
		}
	}
	if(__builtin_ctzll(n)&1){
		for(idt j=0;j<L;++j){
			auto x=f[j],y=f[j+L];
			f[j]=add(x,y),f[j+L]=sub(x,y);
		}
	}
}
template<bool strict=false>inline void dif(u32*A,idt lm){
	switch(lm){
		case 1:if constexpr(strict){
			A[0]=shrk(A[0]);
		}
		break;
		case 2:dif_2<strict>(A[0],A[1]);
		break;
		case 4:dif_4<strict>(A[0],A[1],A[2],A[3]);
		break;
		default:vec_dif_base4<strict>((u32x8*)A,lm/8);
	}
}
template<u32 fx=R>inline void dit(u32*A,idt lm){
	switch(lm){
		case 1:if constexpr(!equals(fx,R)){
			A[0]=mul(A[0],fx);
		}
		break;
		case 2:dit_2<fx>(A[0],A[1]);
		break;
		case 4:dit_4<fx>(A[0],A[1],A[2],A[3]);
		break;
		default:vec_dit_base4<fx>((u32x8*)A,lm/8);
	}
}
inline u32*alc(idt n){
	auto p=new(align_val_t(32))u32[n];
	return p;
}
inline void fre(u32*p){
	::operator delete[](p,align_val_t(32));
}
template<class T,idt al>struct Alloc{
	typedef T value_type;
	T*allocate(idt n){return new(align_val_t(al))T[n];}
	template<class Jok>struct rebind{using other=Alloc<Jok,al>;};
	void deallocate(T*p,idt){::operator delete[](p,align_val_t(al));}
};
using vec=vector<u32,Alloc<u32,32> >;
template<class T>inline T*cpy(T*f,const T*g,idt n){
	return (T*)memcpy(f,g, n*sizeof(T));
}

template<class T>inline T*clr(T*f,idt n){return (T*)memset(f,0,n*sizeof(T));}


void dot(u32*f,const u32*g,idt n){
	vec_op(f,g,n,[](auto&fi,auto&gi){
		fi=mul(fi, gi);
	});
}
void dot(u32*f,const u32*g,const u32*h,idt n){
	vec_op(f,g,h,n,[](auto&fi, auto&gi,auto&hi){
		return fi=mul(gi,hi);
	});
}
void add(u32*f,const u32*g,idt n){
	vec_op(f,g,n,[](auto&fi,auto&gi){
		fi=add(fi, gi);
	});
}
void add(u32*f,const u32*g,const u32*h,idt n){
	vec_op(f,g,h,n,[](auto&fi, auto&gi,auto&hi){
		return fi=add(gi,hi);
	});
}
void sub(u32*f,const u32*g,idt n){
	vec_op(f,g,n,[](auto&fi,auto&gi){
		fi=sub(fi, gi);
	});
}
void sub(u32*f,const u32*g,const u32*h,idt n){
	vec_op(f,g,h,n,[](auto&fi, auto&gi,auto&hi){
		return fi=sub(gi,hi);
	});
}

void vec_multi_iv(u32x8*f,const u32x8*g,idt n){
	if(n==0){
		return;
	}
	f[0]=g[0];
	for(idt i=1;i<n;++i){
		f[i]=mul(f[i-1],g[i]);
	}
	f[n-1]=inv(f[n-1]);
	for(auto i=n-1;i;--i){
		auto ivi=f[i];
		f[i]=mul(ivi,f[i-1]),f[i-1]=mul(ivi,g[i]);
	}
}

template<u32 fx=R>inline void conv(u32*f,u32*g,idt lm){
	dif(f,lm),dif(g,lm), dot(f,g,lm),dit<fx>(f,lm);
}

template<class V,idt alz=16>struct Table{
	using T=typename V::value_type;
	function<void(T*,idt,idt)> f;
	mutable V v;
	template<class F> Table(F f):f{f}{}
	const T*raw()const{return v.data();}
	const T*rsv(idt l)const{
		int ol=v.size();
		if(l>ol){l=max((l+alz-1)&-alz,ol*2),v.resize(l),f(v.data(),ol,l);}
		return raw();
	}
	const T&operator[](idt pos)const{return rsv(pos+1)[pos];}
};
Table<vec>
Id=[](u32*f,idt l,idt r){
	for(;l<8;++l){f[l]=in(l);}
	for(auto i=l;i<r;i+=8){x8(f+i)=add(x8(f+i-8),padd(in(8)));}
},
Iv=[](u32*f,idt l,idt r){
	auto id=Id.rsv(r);
	if(l<8){x8(f)=inv(x8(id)),l=8;}
	vec_multi_iv((u32x8*)(f+l),(const u32x8*)(id+l),(r-l)/8);
};
void inv(u32*f,const u32*g,idt n){
	auto lm=bcl(n);
	auto o=alc(lm*2),h=o+lm;
	f[0]=inv(g[0]);
	for(idt t=2,m=1,xl;t<=lm;m=t,t*=2){
		xl=min(n,t),clr(cpy(o,g,xl)+xl,t-xl),clr(cpy(h,f,m)+m,m),conv(o,h,t);
		clr(o,m),dif(o,t),dot(o,h,t),dit<nR>(o,t),cpy(f+m,o+m,xl-m);
	}
	fre(o);
}
void quo(u32*f,const u32*g,const u32*h,idt n){
	if(n<=64){
		auto lm=bcl(n*2);
		auto o=alc(lm*2),s=o+lm;
		inv(o,h,n),clr(o+n,lm-n),cpy(s,g,n),clr(s+n,lm-n),conv(o,s,lm),cpy(f,
			o,n),fre(o);
		return;
	}
	auto bn=bcl(n)/16,bt=(n+bn-1)/bn,bn2=bn*2;
	auto o=alc(bn2),A=alc(bn2);
	inv(o,h,bn),clr(o+bn,bn),clr(cpy(A,g,bn)+bn,bn),conv(A,o,bn2);
	auto Nh=alc(bn2*bt),nh=Nh,Nf=alc(bn2*(bt-1)),nf=Nf;
	cpy(f,A,bn),clr(cpy(nh,h,bn)+bn,bn),dif(nh,bn2);
	for(idt ds=bn,xl;ds<n;ds+=bn){
		xl=min(bn,n-ds),nh+=bn2;
		clr(cpy(nh,h+ds,xl)+xl,bn2-xl),dif(nh,bn2);
		clr(cpy(nf,f+ds-bn,bn)+bn,bn),dif<1>(nf,bn2),clr(A,bn2),nf+=bn2;
		auto nH=nh,nF=Nf,nH1=nH-bn2;
		for(idt dj=0;dj<ds;dj+=bn,nH-=bn2,nH1-=bn2,nF+=bn2){
			for(idt i=0;i<bn;i+=8){
				x8(A+i)=sub(x8(A+i),_Amul(x8(nH+i), x8(nH1+i),x8(nF+i)));
			}
			for(idt i=bn;i<bn2;i+=8){
				x8(A+i)=sub(x8(A+i),_Smul(x8(nH+i), x8(nH1+i),x8(nF+i)));
			}
		}
		dit(A,bn2),clr(A+bn,bn),add(A,g+ds,xl),dif(A,bn2),dot(A,o,bn2),dit(A,
			bn2),cpy(f+ds,A,xl);
	}
	fre(o),fre(A),fre(Nh),fre(Nf);
}
void ln(u32*f,const u32*g,idt n){
	dot(f,Id.rsv(n),g,n),quo(f,f,g,n),dot(f, Iv.rsv(n),n);
}
template<bool c_inv>void __expi(u32*f,u32*h,const u32*g,idt n){
	f[0]=h[0]=R;
	if(n==1){
		return;
	}
	auto lm=bcl(n);
	auto id=Id.rsv(lm),iv=Iv.rsv(lm);
	auto o=alc(lm*3),A=o+lm,B=A+lm;
	clr(A,lm),A[0]=A[1]=R;
	for(idt t=2,m=1,xl;t<=lm;m=t,t*=2){
		xl=min(n,t),dot(o,id,g,m),dif(o,m),dot(o,A,m),dit(o,m);
		dot(o+m,f,id,m);
		vec_op(o+m,o,m,[](auto&fi,auto&gi){
			fi=sub(fi,gi),clr(gi);
		}),dif(o,t);
		clr(cpy(B,h,m)+m,m),dif(B,t),dot(o,B,t),dit(o,t),dot(clr(o,m)+m,iv+m,m);
		sub(o+m,g+m,xl-m),dif(o,t),dot(A,o,t),dit<nR>(A,t),cpy(f+m,A+m,xl-m);
		if(c_inv||(t!=lm)){
			cpy(A,f,m),dif(A,min(t*2,lm)),dot(o,A,B,t),dit(o,t),clr(o,m);
			dif(o,t),dot(o,B,t),dit<nR>(o,t),cpy(h+m,o+m,xl-m);
		}
	}
	fre(o);
}
void exp(u32*f,const u32*g,idt n){
	if(n<=64){
		auto p=alc(n);
		return __expi<false>(f,p,g,n),fre(p);
	}
	auto bn=bcl(n)/16,bt=(n+bn-1)/bn,bn2=bn*2;
	auto o=alc(bn2),h=alc(bn2);
	auto id=Id.rsv(n),iv=Iv.rsv(n);
	__expi<true>(f,h,g,bn),clr(h+bn,bn),dif(h,bn2);
	auto Ng=alc(bn2*bt),ng=Ng,Nf=alc(bn2*(bt-1)),nf=Nf;
	dot(ng,g,id,bn),clr(ng+bn,bn),dif(ng,bn2);
	for(idt ds=bn,xl;ds<n;ds+=bn){
		xl=min(bn,n-ds),ng+=bn2,dot(ng,g+ds,id+ds,xl),clr(ng+xl,bn2-xl),
			dif(ng,bn2);
		clr(cpy(nf,f+ds-bn,bn)+bn,bn),dif<1>(nf,bn2),clr(o,bn2),nf+=bn2;
		auto nG=ng,nF=Nf,nG1=nG-bn2;
		for(idt dj=0;dj<ds;dj+=bn,nG-=bn2,nG1-=bn2,nF+=bn2){
			for(idt i=0;i<bn;i+=8){
				x8(o+i)=sub(x8(o+i),_Amul(x8(nG+i), x8(nG1+i),x8(nF+i)));
			}
			for(idt i=bn;i<bn2;i+=8){
				x8(o+i)=sub(x8(o+i),_Smul(x8(nG+i), x8(nG1+i),x8(nF+i)));
			}
		}
		dit(o,bn2),clr(o+bn,bn),dif(o,bn2),dot(o,h,bn2),dit<nR>(o,bn2);
		dot(o,iv+ds,xl),clr(o+xl,bn2-xl),dif(o,bn2),dot(o,Nf,bn2),dit(o,bn2),
			cpy(f+ds,o,xl);
	}
	fre(o),fre(h),fre(Ng),fre(Nf);
}
template<bool c_inv>void __sqrti(u32*f,u32*h,const u32*g,idt n){
	auto lm=bcl(n);
	auto o=alc(lm*3),H=o+lm,F=H+lm;
	f[0]=h[0]=F[0]=R;
	for(idt t=2,m=1,xl;t<=lm;m=t,t*=2){
		xl=min(t,n),dot(F,F,m),dit(F,m);
		vec_op(F,F+m,g,g+m,m,[](auto&a0,auto&a1,auto&b0,auto&b1){
			a1=sub(sub(a0,b0),b1),clr(a0);
		});
		clr(cpy(H,h,m)+m,m),conv<nHalf>(F,H,t),cpy(f+m,F+m,xl-m);
		if(c_inv||(t!=lm)){
			dif(cpy(o,f,t),t),cpy(F,o,t),dot(o,H,t),dit(o,t),dif(clr(o,m),t),
				dot(o,H,t),dit<nR>(o,t),cpy(h+m,o+m,xl-m);
		}
	}
	fre(o);
}
void sqrt(u32*f,const u32*g,idt n){
	if(n<=64){
		auto p=alc(n);
		return __sqrti<false>(f,p,g,n),fre(p);
	}
	auto bn=bcl(n)/16,bt=(n+bn-1)/bn,bn2=bn*2;
	auto o=alc(bn2),jok=alc(bn2);
	__sqrti<true>(f,o,g,bn),clr(o+bn,bn),dif(o,bn2);
	auto Nf=alc(bn2*(bt-1)),nf=Nf;
	for(idt ds=bn,xl;ds<n;ds+=bn){
		xl=min(bn,n-ds),clr(cpy(nf,f+ds-bn,bn)+bn,bn),dif<1>(nf,bn2),nf+=bn2;
		auto nF=nf,nF1=nf-bn2,NF=Nf;
		for(idt i=0;i<bn;i+=8){
			x8(jok+i)=neg(mul(x8(nF1+i),x8(NF+i)));
		}
		for(idt i=bn;i<bn2;i+=8){
			x8(jok+i)=mul(x8(nF1+i),x8(NF+i));
		}
		for(idt dj=bn;nF-=bn2,nF1-=bn2,NF+=bn2,dj<ds;dj+=bn){
			for(idt i=0;i<bn;i+=8){
				x8(jok+i)=sub(x8(jok+i),_Amul(x8(nF1+i), x8(nF+i),x8(NF+i)));
			}
			for(idt i=bn;i<bn2;i+=8){
				x8(jok+i)=sub(x8(jok+i),_Smul(x8(nF+i), x8(nF1+i),x8(NF+i)));
			}
		}
		dit<nR>(jok,bn2),clr(jok+bn,bn),sub(jok,g+ds,xl),dif(jok,bn2),dot(jok,
			o,bn2),dit<nHalf>(jok,bn2),cpy(f+ds,jok,xl);
	}
	fre(o),fre(jok),fre(Nf);
}

vec read(const mint *f,int n,int l){
	static_assert(sizeof(mint)==sizeof(u32));
	vec a(l);
	int i=0;
	for(;i+7<n;i+=8){
		auto x=(u32x8)_mm256_loadu_si256((const __m256i_u*)(f+i));
		x8(a.data()+i)=mul(x,padd(R2));
	}
	for(;i<n;i++)a[i]=in(f[i].x);
	return a;
}
void write(mint *f,const u32 *a,int n){
	int i=0;
	for(;i+7<n;i+=8)storeu(f+i,shrk(mul(x8(a+i),padd(1))));
	for(;i<n;i++)f[i].x=out(a[i]);
}
}
#pragma GCC pop_options

void fft(mint *a,int n){
	assert(n>0&&!(n&(n-1))&&n<=(1<<FFT_MAX));
	auto b=Poly_Fast::read(a,n,n);
	Poly_Fast::dif(b.data(),n);Poly_Fast::write(a,b.data(),n);
}
void invFft(mint *a,int n){
	assert(n>0&&!(n&(n-1))&&n<=(1<<FFT_MAX));
	auto b=Poly_Fast::read(a,n,n);
	Poly_Fast::dit(b.data(),n);Poly_Fast::write(a,b.data(),n);
}
void fft(poly &a){fft(a.data(),a.size());}
void invFft(poly &a){invFft(a.data(),a.size());}

// AVX-512 incomplete NTT: the last 16 points use local convolution.
#pragma GCC push_options
#pragma GCC target("avx512f")
#pragma GCC optimize("O3")
namespace Poly_Fast512{
using u32=uint32_t;
using u64=uint64_t;
using idt=size_t;
using I128=__m128i;
using I256=__m256i;
using I512=__m512i;
struct Mont{
	u32 M,M2,niv,R,R2;
	Mont()=default;
	Mont(u32 m):M{m},M2{m*2},niv{2+m},R{-m%m},R2((-u64(m))%m){
		for(u32 i=0;i<4;++i){niv*=2+m*niv;}
	}
	template<bool cond=true>u32 shr(u32 x)const{
		if constexpr(cond){
			return min(x,x-M);
		}
		return x;
	}
	u32 shr2(u32 x)const{
		return min(x,x-M2);
	}
	u32 dil(u32 x)const{
		return min(x,x+M);
	}
	u32 dil2(u32 x)const{
		return min(x,x+M2);
	}
	u32 reduce(u64 x)const{
		return (x+u64(u32(x)*niv)*M)>>32;
	}
	template<bool shrk=false>u32 mul(u32 x,u32 y)const{
		return shr<shrk>(reduce(u64(x)*y));
	}
	template<bool shrk=false>u32 qpw(u32 a,u32 b,u32 r)const{
		for(;b;b>>=1,a=mul(a,a)){
			b&1?r=mul(r,a):r;
		}
		return shr<shrk>(r);
	}
	template<bool shrk=false>u32 qpw(u32 a,u32 b)const{
		return qpw<shrk>(a,b,R);
	}
	template<bool shrk=false>u32 inv(u32 a)const{
		return qpw<shrk>(a,M-2);
	}
	template<bool shrk=false>u32 in(u32 x)const{
		return mul<shrk>(x,R2);
	}
	u32 add(u32 x,u32 y)const{
		return shr2(x+y);
	}
	u32 sub(u32 x,u32 y)const{
		return dil2(x-y);
	}
	u32 Ladd(u32 x,u32 y)const{
		return x+y;
	}
	u32 Lsub(u32 x,u32 y)const{
		return x+M2-y;
	}
	u32 neg(u32 x)const{
		return M2-x;
	}
};
struct Mont16{
	I512 M,M2,niv;
	Mont16()=default;
	Mont16(Mont mt):M{_mm512_set1_epi32(mt.M)},M2{_mm512_set1_epi32(mt.M2)},
		niv{_mm512_set1_epi32(mt.niv)}{}
	template<bool cond=true>I512 shr(I512 x)const{
		if constexpr(cond){
			return _mm512_min_epu32(x,_mm512_sub_epi32(x,M));
		}
		return x;
	}
	I512 shr2(I512 x)const{
		return _mm512_min_epu32(x,_mm512_sub_epi32(x,M2));
	}
	I512 dil(I512 x)const{
		return _mm512_min_epu32(x,_mm512_add_epi32(x,M));
	}
	I512 dil2(I512 x)const{
		return _mm512_min_epu32(x,_mm512_add_epi32(x,M2));
	}
	I512 _mul_hi(I512 x,I512 y)const{
		I512 a=_mm512_mul_epu32(x,y);
		return _mm512_add_epi64(a,_mm512_mul_epu32(_mm512_mul_epu32(a,niv),M));
	}
	template<bool shrk=false>I512 mul(I512 x,I512 y)const{
		I512 a=_mm512_mul_epu32(x,y),b=_mm512_mul_epu32(_mm512_srli_epi64(x,
			32),_mm512_srli_epi64(y,32));
		I512 c=_mm512_mul_epu32(a,niv),d=_mm512_mul_epu32(b,niv);
		c=_mm512_mul_epu32(c,M),d=_mm512_mul_epu32(d,M);
		return shr<shrk>(_mm512_mask_blend_epi32(0xaaaa,_mm512_srli_epi64(
			_mm512_add_epi64(a,c),32),_mm512_add_epi64(b,d)));
	}
	template<bool shrk=false>I512 mul_sm(I512 x,I512 y)const{
		I512 a=_mm512_mul_epu32(x,y),b=_mm512_mul_epu32(_mm512_srli_epi64(x,
			32),y);
		I512 c=_mm512_mul_epu32(a,niv),d=_mm512_mul_epu32(b,niv);
		c=_mm512_mul_epu32(c,M),d=_mm512_mul_epu32(d,M);
		return shr<shrk>(_mm512_mask_blend_epi32(0xaaaa,_mm512_srli_epi64(
			_mm512_add_epi64(a,c),32),_mm512_add_epi64(b,d)));
	}
	template<bool shrk=false>I512 mul_hi(I512 x,I512 y)const{
		return shr<shrk>(_mm512_srli_epi64(_mul_hi(x,y),32));
	}
	I512 add(I512 x,I512 y)const{
		return shr2(_mm512_add_epi32(x,y));
	}
	I512 sub(I512 x,I512 y)const{
		return dil2(_mm512_sub_epi32(x,y));
	}
	I512 Ladd(I512 x,I512 y)const{
		return _mm512_add_epi32(x,y);
	}
	I512 Lsub(I512 x,I512 y)const{
		return _mm512_add_epi32(_mm512_sub_epi32(x,y),M2);
	}
	I512 neg(I512 x)const{
		return _mm512_sub_epi32(M2,x);
	}
	template<u32 z>I512 neg_m(I512 x)const{
		return _mm512_mask_sub_epi32(x,z,M2,x);
	}
};
inline u32*alc(idt n){return new(align_val_t(64))u32[n];}
inline void fre(u32*p){::operator delete[](p,align_val_t(64));}
template<class T>inline T*cpy(T*f,const T*g,idt n){return (T*)memcpy(f,g,
	n*sizeof(T));}
template<class T>inline T*clr(T*f,idt n){return (T*)memset(f,0,n*sizeof(T));}
constexpr u32 mxlg=26,_lg_iter_thre=6;
struct NTT{
	Mont mt;
	u32 RT1[mxlg]{};
	Mont16 ms;
	alignas(32) array<u64,4> rt3[mxlg-6]{},rt3i[mxlg-6]{};
	alignas(32) array<u64,4> st[(mxlg-_lg_iter_thre)>>1]{};
	alignas(32) array<u64,4> st2[_lg_iter_thre>>1]{},bwb{};
	I512 imgx16;
	alignas(64) static constexpr u32
	idx[]={0,2,0,4,0,2,0,4,0,2,0,4,0,2,0,4},
	id2x[]={8,0,9,1,10,2,11,3,12,4,13,5,14,6,15,7},
	id2i[]={0,2,4,6,8,10,12,14,1,3,5,7,9,11,13,15};
	NTT()=default;
	NTT(u32 M,u32 k,u32 _g):mt{M},ms{mt}{
		u32 rt1[mxlg-1],rt1i[mxlg-1];
		rt1[k-2]=_g,rt1i[k-2]=mt.inv(_g);
		for(u32 i=k-2;i>0;--i){
			rt1[i-1]=mt.mul(rt1[i],rt1[i]);
			rt1i[i-1]=mt.mul(rt1i[i],rt1i[i]);
		}
		for(u32 i=1;i<=k;++i){
			RT1[i-1]=mt.qpw<true>(_g,3<<(k-i));
		}
		u32 pr=mt.R,pri=mt.R;
		bwb={mt.R,mt.R,mt.M-mt.R};
		for(u32 i=0;i<k-6;pr=mt.mul(pr,rt1i[i+1]),pri=mt.mul(pri,rt1[i+1]),++i){
			u32 r=mt.mul<true>(pr,rt1[i+1]),ri=mt.mul<true>(pri,rt1i[i+1]);
			u32 r2=mt.mul<true>(r,r),r2i=mt.mul<true>(ri,ri);
			u32 r3=mt.mul<true>(r,r2),r3i=mt.mul<true>(ri,r2i);
			rt3[i]={r,r2,r3},rt3i[i]={ri,r2i,r3i};
		}
		u32 w[8],wi[8];
		w[0]=mt.R,wi[0]=mt.R;
		for(u32 i=0;i<3;++i){
			pr=rt1[i],pri=rt1i[i];
			for(u32 j=1<<i,k=0;k<j;++k){
				w[j+k]=mt.mul<true>(w[k],pr);
				wi[j+k]=mt.mul<true>(wi[k],pri);
			}
		}
		imgx16=_mm512_set1_epi32(w[1]);
	}
	template<bool first>static void forward4(I512*p0,I512*p1,I512*p2,I512*p3,
		I512 img,I512 r1,I512 r2,I512 nr3,const Mont16&ms){
		if constexpr(first){
			//in : [0,2mod)
			auto f1=_mm512_load_si512(p1);
			auto f3=_mm512_load_si512(p3);
			auto g1=ms.add(f1,f3);
			auto g3=ms.mul_sm(ms.Lsub(f1,f3),img);
			auto f0=_mm512_load_si512(p0);
			auto f2=_mm512_load_si512(p2);
			auto g0=ms.add(f0,f2);
			auto g2=ms.sub(f0,f2);
			_mm512_store_si512(p0,ms.add(g0,g1));
			_mm512_store_si512(p1,ms.Lsub(g0,g1));
			_mm512_store_si512(p2,ms.Ladd(g2,g3));
			_mm512_store_si512(p3,ms.Lsub(g2,g3));
			return;
		}
		//in : [0,4mod)
		auto f1=ms.mul_sm(_mm512_load_si512(p1),r1);
		auto nf3=ms.mul_sm(_mm512_load_si512(p3),nr3);
		auto g1=ms.sub(f1,nf3);
		auto f0=ms.shr2(_mm512_load_si512(p0));
		auto f2=ms.mul_sm(_mm512_load_si512(p2),r2);
		auto g3=ms.mul_sm(ms.Ladd(f1,nf3),img);
		auto g0=ms.add(f0,f2);
		auto g2=ms.sub(f0,f2);
		_mm512_store_si512(p0,ms.Ladd(g0,g1));
		_mm512_store_si512(p1,ms.Lsub(g0,g1));
		_mm512_store_si512(p2,ms.Ladd(g2,g3));
		_mm512_store_si512(p3,ms.Lsub(g2,g3));
	}
	template<bool first>void dif4(I512*f,idt lm,idt ix){
		const auto ms=this->ms;
		if(lm<=(idt(1)<<_lg_iter_thre)){
			u32 p=0;
			auto id=_mm512_load_si512(idx),img=imgx16;
			idt l=lm,L=lm>>2,yk=ix;
			for(;L;l=L,L>>=2,++p,yk<<=2){
				if constexpr(first){st2[p]=bwb;}
				auto rt=_mm512_zextsi256_si512(_mm256_load_si256((
					const I256*)(st2+p)));
				for(idt i=0,k=yk;i<lm;i+=l,++k){
					auto r1=_mm512_permutexvar_epi32(id,rt);
					auto r2=_mm512_shuffle_epi32(r1,_MM_PERM_BBBB);
					auto nr3=_mm512_shuffle_epi32(r1,_MM_PERM_DDDD);
					I512 tr=_mm512_zextsi256_si512(_mm256_load_si256((
						const I256*)(rt3+__builtin_ctzll(~k))));
					rt=ms.mul_hi<true>(rt,tr);
					for(idt j=0;j<L;++j){
						forward4<false>(f+i+j+L*0,f+i+j+L*1,f+i+j+L*2,
							f+i+j+L*3,img,r1,r2,nr3,ms);
					}
				}
				_mm256_store_si256((I256*)(st2+p),_mm512_castsi512_si256(rt));
			}

			return;
		}
		idt qlm=lm>>2;
		const u32 p=(__builtin_ctzll(lm)-_lg_iter_thre-1)>>1;
		if constexpr(first){
			st[p]=bwb;
		}
		auto rt=_mm512_zextsi256_si512(_mm256_load_si256((const I256*)(st+p)));
		auto r1=_mm512_permutexvar_epi32(_mm512_load_si512(idx),rt),img=imgx16;
		auto r2=_mm512_shuffle_epi32(r1,_MM_PERM_BBBB),nr3=
			_mm512_shuffle_epi32(r1,_MM_PERM_DDDD);
		I512 tr=_mm512_zextsi256_si512(_mm256_load_si256((const I256*)(
			rt3+__builtin_ctzll(~ix))));
		_mm256_store_si256((I256*)(st+p),_mm512_castsi512_si256(
			ms.mul_hi<true>(rt,tr)));
		for(idt j=0;j<qlm;++j){
			forward4<first>(f+j+qlm*0,f+j+qlm*1,f+j+qlm*2,f+j+qlm*3,img,r1,r2,
				nr3,ms);
		}
		dif4<first>(f+qlm*0,qlm,ix<<2|0);
		dif4<false>(f+qlm*1,qlm,ix<<2|1);
		dif4<false>(f+qlm*2,qlm,ix<<2|2);
		dif4<false>(f+qlm*3,qlm,ix<<2|3);
	}
	void dif(I512*f,idt n){
		if(__builtin_ctzll(n)&1){
			const auto ms=this->ms;
			n>>=1;
			for(idt i=0;i<n;++i){
				auto x=_mm512_load_si512(f+i),y=_mm512_load_si512(f+n+i);
				_mm512_store_si512(f+i,ms.add(x,y)),_mm512_store_si512(f+n+i,
					ms.Lsub(x,y));
			}
			dif4<true>(f,n,0),dif4<false>(f+n,n,1);
		}
		else{
			dif4<true>(f,n,0);
		}
	}
	template<bool first,bool shrk=false>static void inverse4(I512*p0,I512*p1,
		I512*p2,I512*p3,I512 img,I512 r1,I512 r2,I512 nr3,const Mont16&ms){
		if constexpr(first){
			auto f2=_mm512_load_si512(p2);
			auto f3=_mm512_load_si512(p3);
			auto g3=ms.mul_sm(ms.Lsub(f3,f2),img);
			auto g2=ms.add(f2,f3);
			auto f0=_mm512_load_si512(p0);
			auto f1=_mm512_load_si512(p1);
			auto g0=ms.add(f0,f1);
			auto g1=ms.sub(f0,f1);
			_mm512_store_si512(p0,ms.shr<shrk>(ms.add(g0,g2)));
			_mm512_store_si512(p1,ms.shr<shrk>(ms.add(g1,g3)));
			_mm512_store_si512(p2,ms.shr<shrk>(ms.sub(g0,g2)));
			_mm512_store_si512(p3,ms.shr<shrk>(ms.sub(g1,g3)));
			return;
		}
		auto f0=_mm512_load_si512(p0);
		auto f1=_mm512_load_si512(p1);
		auto nf2=ms.neg(_mm512_load_si512(p2));
		auto f3=_mm512_load_si512(p3);
		auto g0=ms.add(f0,f1);
		auto ng2=ms.sub(nf2,f3);
		auto g3=ms.mul_sm(ms.Ladd(nf2,f3),img);
		_mm512_store_si512(p2,ms.mul_sm<shrk>(ms.Ladd(g0,ng2),r2));
		_mm512_store_si512(p0,ms.shr<shrk>(ms.sub(g0,ng2)));
		auto g1=ms.sub(f0,f1);
		_mm512_store_si512(p1,ms.mul_sm<shrk>(ms.Ladd(g1,g3),r1));
		_mm512_store_si512(p3,ms.mul_sm<shrk>(ms.Lsub(g3,g1),nr3));
	}
	template<bool first,bool shrk=false>void dit4(I512*f,idt lm,idt ix){
		const auto ms=this->ms;
		if(lm<=(idt(1)<<_lg_iter_thre)){
			idt yk=ix<<__builtin_ctzll(lm);

			auto id=_mm512_load_si512(idx),img=imgx16;
			u32 p=0;
			idt l=4,L=1;
			for(;yk>>=2,L<lm;L=l,l<<=2,++p){
				if constexpr(first){st2[p]=bwb;}
				auto rt=_mm512_zextsi256_si512(_mm256_load_si256((
					const I256*)(st2+p)));
				for(idt i=0,k=yk;i<lm;i+=l,++k){
					auto r1=_mm512_permutexvar_epi32(id,rt);
					auto r2=_mm512_shuffle_epi32(r1,_MM_PERM_BBBB);
					auto nr3=_mm512_shuffle_epi32(r1,_MM_PERM_DDDD);
					I512 tr=_mm512_zextsi256_si512(_mm256_load_si256((
						const I256*)(rt3i+__builtin_ctzll(~k))));
					rt=ms.mul_hi<true>(rt,tr);
					for(idt j=0;j<L;++j){
						inverse4<false>(f+i+j+L*0,f+i+j+L*1,f+i+j+L*2,
							f+i+j+L*3,img,r1,r2,nr3,ms);
					}
				}
				_mm256_store_si256((I256*)(st2+p),_mm512_castsi512_si256(rt));
			}
			if constexpr(shrk){
				for(idt i=0;i<lm;++i){
					_mm512_store_si512(f+i,ms.shr(_mm512_load_si512(f+i)));
				}
			}
			return;
		}
		idt qlm=lm>>2;
		dit4<first>(f+qlm*0,qlm,ix<<2|0);
		dit4<false>(f+qlm*1,qlm,ix<<2|1);
		dit4<false>(f+qlm*2,qlm,ix<<2|2);
		dit4<false>(f+qlm*3,qlm,ix<<2|3);
		const u32 p=(__builtin_ctzll(lm)-_lg_iter_thre-1)>>1;
		if constexpr(first){st[p]=bwb;}
		auto rt=_mm512_zextsi256_si512(_mm256_load_si256((const I256*)(st+p)));
		auto r1=_mm512_permutexvar_epi32(_mm512_load_si512(idx),rt),img=imgx16;
		auto r2=_mm512_shuffle_epi32(r1,_MM_PERM_BBBB),nr3=
			_mm512_shuffle_epi32(r1,_MM_PERM_DDDD);
		I512 tr=_mm512_zextsi256_si512(_mm256_load_si256((const I256*)(
			rt3i+__builtin_ctzll(~ix))));
		_mm256_store_si256((I256*)(st+p),_mm512_castsi512_si256(
			ms.mul_hi<true>(rt,tr)));
		for(idt j=0;j<qlm;++j){
			inverse4<first,shrk>(f+j+qlm*0,f+j+qlm*1,f+j+qlm*2,f+j+qlm*3,img,
				r1,r2,nr3,ms);
		}
	}
	template<bool shrk=false>void dit(I512*f,idt n){

		if(__builtin_ctzll(n)&1){
			n>>=1;
			dit4<true>(f,n,0),dit4<false>(f+n,n,1);
			const auto ms=this->ms;
			for(idt i=0;i<n;++i){
				auto x=_mm512_load_si512(f+i),y=_mm512_load_si512(f+n+i);
				_mm512_store_si512(f+i,ms.shr<shrk>(ms.add(x,y)));
				_mm512_store_si512(f+n+i,ms.shr<shrk>(ms.sub(x,y)));
			}
		}
		else{
			dit4<true,shrk>(f,n,0);
		}
	}
	//mod x^16-w
	static void conv16(I512*a,I512*b,I512 w,I512 fx,const Mont16&ms){
		auto aa=ms.shr(ms.shr2(_mm512_load_si512(a))),bb=ms.mul_sm<true>(
			_mm512_load_si512(b),fx);
		auto ix=_mm512_set1_epi64(u64(8)<<32),al1=_mm512_set1_epi32(1);
		auto a1=_mm512_permutexvar_epi32(_mm512_load_si512(id2x),aa),a0=
			_mm512_shuffle_epi32(a1,_MM_PERM_CDAB);
		auto a1w=ms.mul_sm<true>(a1,w),a0w=_mm512_shuffle_epi32(a1w,
			_MM_PERM_CDAB);
		auto res0=_mm512_setzero_si512(),res1=_mm512_setzero_si512();
		auto res2=_mm512_setzero_si512(),res3=_mm512_setzero_si512();
		// alignr needs an immediate shift; this unroll is required.
		#pragma GCC unroll 8
		for(u32 i=0;i<8;++i){
			auto b0=_mm512_permutexvar_epi32(ix,bb),b1=_mm512_shuffle_epi32(
				b0,_MM_PERM_CDAB);
			//i,i+8 i+8,i
			auto pr1=(i==0)?a1:_mm512_alignr_epi64(a1,a0,8-i);
			auto pr2=(i==0)?a0:_mm512_alignr_epi64(a0,a1w,8-i);
			auto pr3=(i==0)?a1w:_mm512_alignr_epi64(a1w,a0w,8-i);
			res0=_mm512_add_epi64(res0,_mm512_mul_epu32(b0,pr2));
			res1=_mm512_add_epi64(res1,_mm512_mul_epu32(b0,pr1));
			res2=_mm512_add_epi64(res2,_mm512_mul_epu32(b1,pr3));
			res3=_mm512_add_epi64(res3,_mm512_mul_epu32(b1,pr2));
			ix=_mm512_add_epi32(ix,al1);
		}
		res0=_mm512_add_epi64(res0,res2),res1=_mm512_add_epi64(res1,res3);
		res2=_mm512_sub_epi64(res0,ms.M2),res3=_mm512_sub_epi64(res1,ms.M2);
		res0=_mm512_min_epu64(res0,res2),res1=_mm512_min_epu64(res1,res3);
		res2=_mm512_mul_epu32(res0,ms.niv),res3=_mm512_mul_epu32(res1,ms.niv);
		res2=_mm512_mul_epu32(res2,ms.M),res3=_mm512_mul_epu32(res3,ms.M);
		res0=_mm512_add_epi64(res0,res2),res1=_mm512_add_epi64(res1,res3);
		res0=_mm512_mask_blend_epi32(0xaaaa,_mm512_srli_epi64(res0,32),res1);
		_mm512_store_si512(a,_mm512_permutexvar_epi32(_mm512_load_si512(id2i),
			ms.shr2(res0)));
	}
	void dot_I512(I512*a,const I512*b,idt lm){
		const auto ms=this->ms;
		for(idt i=0;i<lm;++i){
			_mm512_store_si512(a+i,ms.mul(_mm512_load_si512(a+i),
				_mm512_load_si512(b+i)));
		}
	}
	void dot(I512*a,I512*b,idt lm){
		u32 R=mt.R;
		const auto ms=this->ms;
		const auto fx=mt.in<true>(mt.in(mt.M-((mt.M-1)>>(__builtin_ctzll(
			lm)))));
		for(idt i=0;i<lm;++i){
			conv16(a+i,b+i,_mm512_set1_epi32(R),_mm512_set1_epi32(fx),ms);
			R=mt.mul(R,RT1[__builtin_ctzll(~i)]);
		}
	}
};

poly mul(const poly &a,const poly &b){
	int m=a.size()+b.size()-1,n=16;
	while(n<m)n<<=1;
	auto f=alc(n),g=alc(n);
	for(int i=0;i<(int)a.size();i++)f[i]=a[i].x;
	for(int i=0;i<(int)b.size();i++)g[i]=b[i].x;
	clr(f+a.size(),n-a.size());clr(g+b.size(),n-b.size());
	static NTT ntt{998244353,23,752838388};
	ntt.dif((I512*)f,n>>4);
	ntt.dif((I512*)g,n>>4);
	ntt.dot((I512*)f,(I512*)g,n>>4);
	ntt.dit<true>((I512*)f,n>>4);
	poly h(m);
	for(int i=0;i<m;i++)h[i].x=f[i];
	fre(f);fre(g);
	return h;
}

}
#pragma GCC pop_options

poly operator*(const poly &f,const poly &g){
	if(f.empty()||g.empty())return{};
	int n=f.size()+g.size()-1;
	assert(n<=(1<<FFT_MAX));
	if(min(f.size(),g.size())<=32){
		poly h(n);
		for(int i=0;i<(int)f.size();i++)
			for(int j=0;j<(int)g.size();j++)h[i+j]+=f[i]*g[j];
		return h;
	}
	static const bool avx512=__builtin_cpu_supports("avx512f");
	if(avx512)return Poly_Fast512::mul(f,g);
	int l=Poly_Fast::bcl(n);
	auto a=Poly_Fast::read(f.data(),f.size(),l);
	auto b=Poly_Fast::read(g.data(),g.size(),l);
	Poly_Fast::conv(a.data(),b.data(),l);
	poly h(n);Poly_Fast::write(h.data(),a.data(),n);
	return h;
}

poly operator/(poly f,poly g){
	if(f.empty()||g.empty())return poly(f.size());
	int m=g.size();
	reverse(g.begin(),g.end());
	g=f*g;
	for(int i=0;i<(int)f.size();i++)f[i]=g[i+m-1];
	return f;
}

void operator+=(poly &f,const poly &g){
	if(f.size()<g.size())f.resize(g.size());
	for(int i=0;i<(int)g.size();i++)f[i]+=g[i];
}

poly operator+(poly f,const poly &g){f+=g;return f;}

void operator-=(poly &f,const poly &g){
	if(f.size()<g.size())f.resize(g.size());
	for(int i=0;i<(int)g.size();i++)f[i]-=g[i];
}

poly operator-(poly f,const poly &g){f-=g;return f;}

poly operator*(poly f,mint x){
	for(mint &i:f)i*=x;
	return f;
}

poly operator*(mint x,poly f){
	for(mint &i:f)i*=x;
	return f;
}

mint value(const poly &f,mint x){
	mint ans=0;
	for(int i=f.size()-1;i>=0;i--)ans=ans*x+f[i];
	return ans;
}

// Formal inverse; constant term must be invertible.
poly Inv(poly f){
	assert(!f.empty()&&f[0]&&f.size()<=(1<<FFT_MAX));
	int n=f.size(),l=Poly_Fast::bcl(n);
	auto a=Poly_Fast::read(f.data(),n,l);
	Poly_Fast::vec b(l);
	Poly_Fast::inv(b.data(),a.data(),n);
	Poly_Fast::write(f.data(),b.data(),n);
	return f;
}

poly integ(poly f){
	int n=f.size();
	::init(n);
	f.resize(n+1);
	for(int i=n;i>=1;i--)f[i]=f[i-1]*::inv[i];
	f[0]=0;
	return f;
}

poly diff(poly f){
	if(f.empty())return{};
	int n=f.size();
	for(int i=0;i<n-1;i++)
		f[i]=f[i+1]*(i+1);
	f.pop_back();
	return f;
}

// Formal logarithm; f[0]=1. Integration denominators must be invertible.
poly Ln(poly f){
	assert(!f.empty()&&f[0].x==1&&f.size()<=(1<<FFT_MAX));
	int n=f.size(),l=Poly_Fast::bcl(n);
	auto a=Poly_Fast::read(f.data(),n,l);
	Poly_Fast::vec b(l);
	Poly_Fast::ln(b.data(),a.data(),n);
	Poly_Fast::write(f.data(),b.data(),n);
	return f;
}

// f[0]=0; block FPS with cached transforms.
poly Exp(poly f){
	assert(!f.empty()&&!f[0]&&f.size()<=(1<<FFT_MAX));
	int n=f.size(),l=Poly_Fast::bcl(n);
	auto a=Poly_Fast::read(f.data(),n,l);
	Poly_Fast::vec b(l);
	Poly_Fast::exp(b.data(),a.data(),n);
	Poly_Fast::write(f.data(),b.data(),n);
	return f;
}

// f[0]=1; choose the square root with constant term 1.
poly Sqrt(poly f){
	assert(!f.empty()&&f[0].x==1&&f.size()<=(1<<FFT_MAX));
	int n=f.size(),l=Poly_Fast::bcl(n);
	auto a=Poly_Fast::read(f.data(),n,l);
	Poly_Fast::vec b(l);
	Poly_Fast::sqrt(b.data(),a.data(),n);
	Poly_Fast::write(f.data(),b.data(),n);
	return f;
}

// f=q*g+r; g.back()!=0. Unlike operator/, this is division.
void Div(poly f,poly g,poly &q,poly &r){
	assert(!g.empty()&&g.back()&&(&q!=&r));
	if(f.size()<g.size()){q={0};r=move(f);return;}
	int n=f.size(),m=g.size(),k=n-m+1;
	poly a(f.rbegin(),f.rbegin()+k),b(g.rbegin(),g.rend());
	b.resize(k);q=a*Inv(b);q.resize(k);
	reverse(q.begin(),q.end());
	b=q*g;r.resize(m-1);
	for(int i=0;i<m-1;i++)r[i]=f[i]-b[i];
}

// Return f(x[i]); repeated evaluation points are allowed.
poly Eval(poly f,poly x){
	int n=x.size();
	if(!n)return{};
	vector<poly>g(n*4);
	auto build=[&](auto &&self,int u,int l,int r)->void{
		if(l==r){g[u]={-x[l],1};return;}
		int m=(l+r)>>1;
		self(self,u*2,l,m);self(self,u*2+1,m+1,r);
		g[u]=g[u*2]*g[u*2+1];
	};
	build(build,1,0,n-1);
	poly ans(n);
	auto solve=[&](auto &&self,poly a,int u,int l,int r)->void{
		if(r-l<32){
			for(int i=l;i<=r;i++)ans[i]=value(a,x[i]);
			return;
		}
		int m=(l+r)>>1;
		poly q,b;
		Div(a,g[u*2],q,b);self(self,move(b),u*2,l,m);
		Div(move(a),g[u*2+1],q,b);self(self,move(b),u*2+1,m+1,r);
	};
	solve(solve,move(f),1,0,n-1);
	return ans;
}

poly BM(poly a){
	poly C{1},B{1};
	int L=0,m=1;
	mint b=1;
	for(int n=0;n<(int)a.size();n++){
		mint d=0;
		for(int i=0;i<=L;i++)d+=C[i]*a[n-i];
		if(!d){
			m++;
		}else{
			poly T=C;
			mint coef=d*b.inv();
			poly xmB(m,0);
			for(mint val:B)xmB.push_back(val);
			if(C.size()<xmB.size())C.resize(xmB.size(),0);
			for(int i=0;i<(int)xmB.size();i++)C[i]-=coef*xmB[i];
			if(2*L<=n){
				L=n+1-L;
				B=T;
				b=d;
				m=1;
			}else{
				m++;
			}
		}
	}
	return C;
}

mint FSPE(poly F,poly G,ll t){
	// find [x^t] F/G
	assert(t>=0&&!G.empty()&&G[0]);
	if(F.empty())return 0;
	while(t){
		poly G1=G;
		for(int i=1;i<(int)G1.size();i+=2)G1[i]=-G1[i];
		F=F*G1;G=G*G1;
		poly f,g;
		for(int i=t&1;i<(int)F.size();i+=2)f.push_back(F[i]);
		for(int i=0;i<(int)G.size();i+=2)g.push_back(G[i]);
		t>>=1;
		F=f;G=g;
		if(F.empty())return 0;
	}
	return F[0]/G[0];
}

mint RSPE(poly f,ll x){
	// find f[x] in O(n^2), note that f must be recursion
	poly g=BM(f);
	poly s(g.size()-1);
	for(int i=0;i<(int)s.size();i++)
		for(int j=0;j<=i&&j<(int)f.size();j++)
			s[i]+=f[j]*g[i-j];
	return FSPE(s,g,x);
}

poly lagrange(vector<mint>x,vector<mint>y){
	int n=x.size();assert(y.size()==x.size());
	if(!n)return{};
	poly g(n+1),q(n),f(n);g[0]=1;
	for(int i=0;i<n;i++){
		for(int j=i+1;j>=1;j--)g[j]=g[j-1]-x[i]*g[j];
		g[0]*=-x[i];
	}
	for(int i=0;i<n;i++){
		q[n-1]=g[n];
		for(int j=n-2;j>=0;j--)q[j]=g[j+1]+x[i]*q[j+1];
		mint d=value(q,x[i]);assert(d);
		mint c=y[i]/d;
		for(int j=0;j<n;j++)f[j]+=q[j]*c;
	}
	return f;
}
poly to_ex(poly f){ // f(e^x)
	int n=f.size();
	if(!n)return{};
	::init(n);
	function<pair<poly,poly>(int,int)>solve=[&](int l,int r){
		if(l==r)return make_pair(poly{f[l]},poly{1,MOD-l});
		int mid=l+r>>1;
		auto ls=solve(l,mid),rs=solve(mid+1,r);
		return make_pair(ls.first*rs.second+ls.second*rs.first,
			ls.second*rs.second);
	};
	auto ans=solve(0,n-1);
	poly g=ans.first*Inv(ans.second);
	g.resize(n);
	for(int i=0;i<n;i++)g[i]=g[i]*ifac[i];
	return g;
}

poly S2line(int n){
	// return S(n,i)
	assert(n>=0);::init(n);
	poly F(n+1),G(n+1);
	// S(n,i) * binom(j,i) -> F(j)
	for(int i=0;i<=n;i++){
		G[i]=ifac[i];F[i]=mint(i).pow(n)*ifac[i];
		if(i&1)F[i]*=-1;
	}
	F=F*G;
	F.resize(n+1);
	for(int i=0;i<=n;i++){
		if(i&1)F[i]*=-1;
	}
	return F;
}
// use init(n) before accessing fac / ifac / inv directly
// end for polynomial/ntt-fast.cpp
/////////////////////////
// !!!!! Requires GNU C++17/20, AVX2 and <immintrin.h>.
// Contains mint already; do not paste another mint or poly. !!!!
