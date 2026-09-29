const ld pi=acos(-1),dx=1e-9;
struct Point{
	ll x,y;
	friend istream& operator >>(istream &i,Point &P){i>>P.x>>P.y;return i;}
	friend ostream& operator <<(ostream &o,const Point &P){o<<P.x<<' '<<P.y;
		return o;}
	// void read(){cin>>x>>y;}
	// void output(){cout<<x<<' '<<y<<'\n';}
	void output(string s=""){cerr<<"Point "<<s<<": ("<<x<<","<<y<<")"<<endl;}
	double arg(){return atan2(y,x);}
	Point (ll _x=0,ll _y=0){x=_x,y=_y;}
};

inline bool operator==(const Point &x,const Point &y){return
	x.x==y.x&&x.y==y.y;}
inline bool operator!=(const Point &x,const Point &y){return
	x.x!=y.x||x.y!=y.y;}

inline bool operator <(const Point &x,const Point &y){return
	(x.x!=y.x?x.x<y.x:x.y<y.y);}
inline Point operator-(const Point &x,const Point &y){return Point(x.x-y.x,
	x.y-y.y);}
inline Point operator+(const Point &x,const Point &y){return Point(x.x+y.x,
	x.y+y.y);}
inline void operator-=(Point &x,const Point &y){x=x-y;}
inline void operator+=(Point &x,const Point &y){x=x+y;}
inline Point operator*(const Point &x,ll y){return Point(x.x*y,x.y*y);}
inline void operator*=(Point &x,ll y){x=x*y;}
inline Point operator*(ll y,const Point &x){return Point(x.x*y,x.y*y);}
inline void operator*=(ll y,Point &x){x=x*y;}

inline ll operator*(const Point &x,const Point &y){
	// return OX · OY · sin(XOY) = 2 * area(OXY)
	return x.x*y.y-x.y*y.x;
}
inline ll operator ^(const Point &x,const Point &y){
	// return OX · OY · cos(XOY)
	return x.x*y.x+x.y*y.y;
}

int sgn(ll x){return (x==0?0:(x>0?1:-1));}

inline int sgn(const Point &P){
	return P.y<0?-1:P.y>0||P.x<0?1:0;
}
inline bool cmp(const Point &x,const Point &y){
	if(sgn(x)!=sgn(y)) return sgn(x)<sgn(y);
	return x*y>0;
}
inline ll area(const Point &x,const Point &y,const Point &z){
	// return triangle XYZ area **times 2**
	return abs((x-y)*(x-z));
}

inline ld dist(const Point &x){return sqrt((ld)x.x*x.x+(ld)x.y*x.y);}
inline ll dist2(const Point &x){return x.x*x.x+x.y*x.y;}

inline int dir(const Point &x,const Point &y,const Point &z){
	// return is YX clockwise direction to YZ
	return sgn((x-y)*(y-z));
}
inline int suf(int x,int n){return (x+1==n?0:x+1);}
inline int pre(int x,int n){return (x==0?n-1:x-1);}
ll convec_area(const vector<Point> &a){
	// return convex A area **times 2**
	ll s=0;int n=a.size();
	for(int i=0;i<n;i++) s+=area(a[0],a[i],a[suf(i,n)]);
	return s;
}
vector<Point> convex(vector<Point> a){
	sort(a.begin(),a.end());a.erase(unique(a.begin(),a.end()),a.end());
	int n=a.size();
	if(n<=2)return a;
	vector<Point> b;int s=0,t=0;
	for(int i=2;i--;s=t,reverse(a.begin(),a.end())){
		for(Point p:a){
			while(t>=s+2&&dir(b[t-2],b[t-1],p)>=0) t--,b.pop_back();
			b.push_back(p);t++;
		}
		b.pop_back();t--;
	}
	return b;
}
inline bool on_segment(const Point &X,const Point &A,const Point &B){
	// check if X on segment AB
	if((A-B)*(X-B)!=0) return false;
	return min(A.x,B.x)<=X.x&&X.x<=max(A.x,B.x)&&
		min(A.y,B.y)<=X.y&&X.y<=max(A.y,B.y);
}
inline bool strict_intersect(const Point &A,const Point &B,
	const Point &C,const Point &D){
	// check if AB,CD (except endpoint)strictly intersect
	return dir(A,C,D)*dir(B,C,D)==-1&&dir(C,A,B)*dir(D,A,B)==-1;
}
inline bool nonstrict_intersect(const Point &A,const Point &B,
	const Point &C,const Point &D){
	// check if AB,CD (include endpoints) intersect
	return (dir(A,C,D)*dir(B,C,D)==-1&&dir(C,A,B)*dir(D,A,B)==-1)||
	on_segment(A,C,D)||on_segment(B,C,D)||on_segment(C,A,B)||on_segment(D,A,B);
}
inline bool parallel(Point A,Point B,
	Point C,Point D){
	// check if AB // CD
	B-=A,D-=C;A-=C;C=-1*D;
	return sgn(B*C)==0&&sgn(A*B)!=0;
}
inline bool ray_intersect(Point A,Point B,Point C,Point D){
	// closed rays; a zero direction is a single point
	Point u=B-A,v=D-C,w=C-A;
	if(u==Point())return v==Point()?A==C:((A-C)*v==0&&((A-C)^v)>=0);
	if(v==Point())return w*u==0&&(w^u)>=0;
	ll d=u*v;
	if(!d)return w*u==0&&((u^v)>0||(w^u)>=0);
	return sgn(w*v)*sgn(d)>=0&&sgn(w*u)*sgn(d)>=0;
}
inline bool in_triangle(const Point &P,const Point &A,
	const Point &B,const Point &C){
	// check if P in triangle ABC
	// P may lie on bound
	if(!area(A,B,C))return on_segment(P,A,B)||on_segment(P,B,C)||on_segment(P,
		C,A);
	return area(A,B,C)==area(P,A,B)+area(P,B,C)+area(P,C,A);
}
template<typename F> int convex_min(int n,F val){
	assert(n>0);
	auto low=[&](int i){return val(i)<=val((i+n-1)%n)&&val(i)<=val((i+1)%n);};
	if(low(0))return 0;
	int l=0,r=n;
	while(l+1<r){
		int mid=(l+r)>>1;
		if(low(mid))return mid;
		bool a=val((l+1)%n)>=val(l),b=val((mid+1)%n)>=val(mid);
		if(a!=b){if(a)l=mid;else r=mid;}
		else if((val(mid)>val(l))==a)l=mid;
		else r=mid;
	}
	return r%n;
}
inline ll min_cross(const vector<Point> &I,Point X){
	// convex polygon in cyclic order, without collinear intermediate points
	auto val=[&](int i){return I[i]*X;};
	return val(convex_min(I.size(),val));
}
inline ld min_cross(const vector<Point> &I,ld x,ld y){
	auto val=[&](int i){return I[i].x*x+I[i].y*y;};
	return val(convex_min(I.size(),val));
}
