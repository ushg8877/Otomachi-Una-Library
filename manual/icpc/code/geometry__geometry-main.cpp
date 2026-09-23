const ld pi=acos(-1),dx=1e-9;
struct Point{
	ll x,y;
	friend istream& operator >>(istream &i,Point &P){i>>P.x>>P.y;return i;}
	friend ostream& operator <<(ostream &o,Point &P){o<<P.x<<P.y;return o;}
	// void read(){cin>>x>>y;}
	// void output(){cout<<x<<' '<<y<<'\n';}
	void debug(string s=""){cerr<<"Point "<<s<<": ("<<x<<","<<y<<")"<<endl;}
	double arg(){return atan2(y,x);}
	Point (ll _x=0,ll _y=0){x=_x,y=_y;}
};

inline bool operator ==(const Point &x,const Point &y){return x.x==y.x&&x.y==y.y;}
inline bool operator !=(const Point &x,const Point &y){return x.x!=y.x||x.y!=y.y;}

inline bool operator <(const Point &x,const Point &y){return (x.x!=y.x?x.x<y.x:x.y<y.y);}
inline Point operator -(const Point &x,const Point &y){return Point(x.x-y.x,x.y-y.y);}
inline Point operator +(const Point &x,const Point &y){return Point(x.x+y.x,x.y+y.y);}
inline void operator -=(Point &x,const Point &y){x=x-y;}
inline void operator +=(Point &x,const Point &y){x=x+y;}
inline Point operator *(const Point &x,ll y){return Point(x.x*y,x.y*y);}
inline void operator *=(Point &x,ll y){x=x*y;}
inline Point operator *(ll y,const Point &x){return Point(x.x*y,x.y*y);}
inline void operator *=(ll y,Point &x){x=x*y;}

inline ll operator *(const Point &x,const Point &y){
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

inline ld dist(const Point &x){return sqrt(x.x*x.x+x.y*x.y);}
inline ll dist2(const Point &x){return x.x*x.x+x.y*x.y;}

inline int dir(const Point &x,const Point &y,const Point &z){
	// return is YX clockwise direction to YZ 
	return sgn((x-y)*(y-z));
}
inline int suf(int x,int n){return (x+1==n?0:x+1);}
inline int pre(int x,int n){return (x==0?n-1:x-1);}
ll convec_area(vector<Point> a){
	// return convex A area **times 2**
	ll s=0;int n=a.size();
	for(int i=0;i<n;i++) s+=area(a[0],a[i],a[suf(i,n)]);
	return s;
}
vector<Point> convex(vector<Point> a){
	int n=a.size();
	if(n<=2) return a;
	sort(a.begin(),a.end());
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
inline bool ray_intersect(Point A,Point B,
	Point C,Point D){
	// check if ray AB,CD (include endpoints) intersect
	B-=A,D-=C;A-=C;C=-1*D;
	// judge exist non-negtive real number x,y s.t. A+xB+yC=O;
	if(sgn(B*C)==0){
		if(sgn(A*B)==0) return true; // coincide
		return false; // parallel
	}
	return sgn(C*B)*sgn(A*B)<=0&&sgn(B*C)*sgn(A*C)<=0;
}
inline bool in_triangle(const Point &P,const Point &A,
	const Point &B,const Point &C){
	// check if P in triangle ABC
	// P may lie on bound
	return area(A,B,C)==area(P,A,B)+area(P,B,C)+area(P,C,A);
}
inline ll min_cross(const vector<Point> &I,Point X){
	// please make sure I is convex, find min i\in I, i*X
	int n=I.size();assert(n>0);
	auto id=[&](int x){
		if(x<0) x+=n;
		if(x>=n) x-=n;
		return x;
	};
	int p=0;
	for(int i=__lg(n);i>=0;i--){
		if(I[p]*X>=I[id(p+(1<<i))]*X) p=id(p+(1<<i));
		if(I[p]*X>=I[id(p-(1<<i))]*X) p=id(p-(1<<i));
	}
	return I[p]*X;
}
inline ld min_cross(const vector<Point> &I,ld x,ld y){
	// please make sure I is convex, find min i\in I, x(i)*x+y(i)*y
	int n=I.size();assert(n>0);
	auto crs=[&](int p){
		return I[p].x*x+I[p].y*y;
	};
	auto id=[&](int x){
		if(x<0) x+=n;
		if(x>=n) x-=n;
		return x;
	};
	int p=0;
	for(int i=__lg(n);i>=0;i--){
		if(crs(p)>=crs(id(p+(1<<i)))) p=id(p+(1<<i));
		if(crs(p)>=crs(id(p-(1<<i)))) p=id(p-(1<<i));
	}
	return crs(p);
}
