////////////////////////////////////////////////////////////////
//
// template for rational number
// usage:
//   Rat a(1,2),b(2,3); cout<<a+b; // 7/6
//   p/q is reduced, q>0; +, -, *, / and comparisons are supported.
//
////////////////////////////////////////////////////////////////
struct Rat{
ll p=0,q=1;
Rat()=default;
Rat(ll a,ll b=1){*this=make(a,b);}
static Rat make(__int128 a,__int128 b){
	assert(b);if(b<0)a=-a,b=-b;
	__int128 x=a<0?-a:a,y=b;
	while(y){__int128 t=x%y;x=y;y=t;}
	a/=x;b/=x;
	assert(LLONG_MIN<=a&&a<=LLONG_MAX&&b<=LLONG_MAX);
	Rat r;r.p=a;r.q=b;return r;
}
friend Rat operator+(Rat a,Rat b){return make((__int128)a.p*b.q+(__int128)b.p*a.q,(__int128)a.q*b.q);}
friend Rat operator-(Rat a,Rat b){return make((__int128)a.p*b.q-(__int128)b.p*a.q,(__int128)a.q*b.q);}
friend Rat operator*(Rat a,Rat b){return make((__int128)a.p*b.p,(__int128)a.q*b.q);}
friend Rat operator/(Rat a,Rat b){return make((__int128)a.p*b.q,(__int128)a.q*b.p);}
Rat operator-()const{return make(-(__int128)p,q);}
Rat& operator+=(Rat b){return *this=*this+b;}
Rat& operator-=(Rat b){return *this=*this-b;}
Rat& operator*=(Rat b){return *this=*this*b;}
Rat& operator/=(Rat b){return *this=*this/b;}
friend bool operator==(Rat a,Rat b){return a.p==b.p&&a.q==b.q;}
friend bool operator!=(Rat a,Rat b){return !(a==b);}
friend bool operator<(Rat a,Rat b){return (__int128)a.p*b.q<(__int128)b.p*a.q;}
friend bool operator>(Rat a,Rat b){return b<a;}
friend bool operator<=(Rat a,Rat b){return !(b<a);}
friend bool operator>=(Rat a,Rat b){return !(a<b);}
friend ostream& operator<<(ostream &o,Rat a){return o<<a.p<<'/'<<a.q;}
};
// end for math/rational.cpp
/////////////////////////
// !!!!! Requires ll; inputs and reduced results must fit ll, q must be positive
// after reduction. Do not modify p/q directly. !!!!
