struct BIT{
private:
struct info{
// modify this if the info template give is not you wish
ll x;
info():x(0){}
info(ll _x):x(_x){}
inline info& operator +=(const info &y){
	x+=y.x;return *this;
}
inline info operator +(const info &y)const{
	return info(*this)+=y;
}
};
int n=0;vector<info> a;
void modify(int x,info v){
	for(ll i=x;i<=n+1;i+=i&-i)a[i]+=v;
}
info query(int x){
	assert(0<=x&&x<=n);
	info ans;
	while(x) ans+=a[x],x-=(x&-x);
	return ans;
}
public:
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX-1);n=_n;a.assign(n+2,info());
}
void add(int x,ll v){
	assert(1<=x&&x<=n+1);
	modify(x,info(v));
}
ll ask(int x){
	return query(x).x;
}
};
struct rars{
// a implentment for range add range query
int n=0;BIT T0,T1;
void setN(int _n){n=_n;T0.setN(n);T1.setN(n);}
void add(int l,int r,ll x){
	assert(1<=l&&l<=r&&r<=n);
	T0.add(l,x);T0.add(r+1,-x);
	T1.add(l,l*x);T1.add(r+1,-(r+1)*x);
}
ll query(int x){return (x+1)*T0.ask(x)-T1.ask(x);}
ll ask(int l,int r){assert(1<=l&&l<=r&&r<=n);return query(r)-query(l-1);}
}T;
