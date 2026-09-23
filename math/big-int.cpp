////////////////////////////////////////////////////////////////
//
// template for big-int
//
// usage:
//   BigInt a,b; cin>>a>>b; cout<<a+b; // also -, *, /, %
//
////////////////////////////////////////////////////////////////
struct BigInt{
	static constexpr int base=1000000000;
	vector<int> a;
	int sign=1;
	BigInt(long long x=0){
		unsigned long long y=x;
		if(x<0)sign=-1,y=0-y;
		for(;y;y/=base)a.push_back(y%base);
	}
	void trim(){while(!a.empty()&&!a.back())a.pop_back();if(a.empty())sign=1;}
	explicit operator bool()const{return !a.empty();}
	int cmp_abs(const BigInt &b)const{
		if(a.size()!=b.a.size())return a.size()<b.a.size()?-1:1;
		for(int i=(int)a.size()-1;i>=0;i--)if(a[i]!=b.a[i])return a[i]<b.a[i]?-1:1;
		return 0;
	}
	int cmp(const BigInt &b)const{return sign!=b.sign?(sign<b.sign?-1:1):sign*cmp_abs(b);}
	bool operator <(const BigInt &b)const{return cmp(b)<0;}
	bool operator >(const BigInt &b)const{return cmp(b)>0;}
	bool operator <=(const BigInt &b)const{return cmp(b)<=0;}
	bool operator >=(const BigInt &b)const{return cmp(b)>=0;}
	bool operator ==(const BigInt &b)const{return sign==b.sign&&a==b.a;}
	bool operator !=(const BigInt &b)const{return !(*this==b);}
	BigInt operator -()const{BigInt b=*this;if(b)b.sign=-b.sign;return b;}
	BigInt operator +()const{return *this;}
	BigInt operator +(const BigInt &b)const{
		if(!b)return *this;
		if(!*this)return b;
		if(sign!=b.sign)return *this-(-b);
		BigInt c;c.sign=sign;c.a.resize(max(a.size(),b.a.size()));
		int carry=0;
		for(int i=0;i<(int)c.a.size();i++){
			long long x=(long long)carry+(i<(int)a.size()?a[i]:0)+(i<(int)b.a.size()?b.a[i]:0);
			c.a[i]=x%base;carry=x/base;
		}
		if(carry)c.a.push_back(carry);
		return c;
	}
	BigInt operator -(const BigInt &b)const{
		if(!b)return *this;
		if(!*this)return -b;
		if(sign!=b.sign)return *this+(-b);
		if(cmp_abs(b)<0)return -(b-*this);
		BigInt c=*this;int carry=0;
		for(int i=0;i<(int)a.size();i++){
			int x=a[i]-carry-(i<(int)b.a.size()?b.a[i]:0);
			carry=x<0;c.a[i]=x+(carry?base:0);
		}
		c.trim();return c;
	}
	BigInt mul(int x)const{
		assert(0<=x&&x<base);
		BigInt c;c.sign=sign;c.a.resize(a.size());long long carry=0;
		for(int i=0;i<(int)a.size();i++){
			long long y=1ll*a[i]*x+carry;c.a[i]=y%base;carry=y/base;
		}
		if(carry)c.a.push_back(carry);
		c.trim();return c;
	}
	BigInt operator *(const BigInt &b)const{
		BigInt c;if(!*this||!b)return c;
		c.sign=sign*b.sign;c.a.resize(a.size()+b.a.size());
		for(int i=0;i<(int)a.size();i++){
			long long carry=0;
			for(int j=0;j<(int)b.a.size();j++){
				long long x=c.a[i+j]+1ll*a[i]*b.a[j]+carry;
				c.a[i+j]=x%base;carry=x/base;
			}
			c.a[i+b.a.size()]=carry;
		}
		c.trim();return c;
	}
	// Division truncates towards zero; remainder has the dividend's sign.
	pair<BigInt,BigInt> divmod(const BigInt &b)const{
		assert(b);
		if(cmp_abs(b)<0)return {0,*this};
		if(b.a.size()==1){
			BigInt q=*this;long long rem=0;
			for(int i=(int)a.size()-1;i>=0;i--){
				long long x=rem*base+a[i];q.a[i]=x/b.a[0];rem=x%b.a[0];
			}
			q.sign=sign*b.sign;q.trim();return {q,BigInt(rem*sign)};
		}
		int norm=base/(b.a.back()+1);
		BigInt x=mul(norm),y=b.mul(norm),q,r;x.sign=y.sign=1;
		q.a.resize(x.a.size());
		for(int i=(int)x.a.size()-1;i>=0;i--){
			r.a.insert(r.a.begin(),x.a[i]);r.trim();
			int n=y.a.size();
			long long hi=r.a.size()>(size_t)n?r.a[n]:0;
			long long lo=r.a.size()>=(size_t)n?r.a[n-1]:0;
			int d=min((hi*base+lo)/y.a.back(),(long long)base-1);
			r=r-y.mul(d);
			while(r.sign<0)r=r+y,--d;
			q.a[i]=d;
		}
		long long carry=0;
		for(int i=(int)r.a.size()-1;i>=0;i--){
			long long x=r.a[i]+carry*base;r.a[i]=x/norm;carry=x%norm;
		}
		q.sign=sign*b.sign;r.sign=sign;q.trim();r.trim();return {q,r};
	}
	BigInt operator /(const BigInt &b)const{return divmod(b).first;}
	BigInt operator %(const BigInt &b)const{return divmod(b).second;}
	BigInt& operator +=(const BigInt &b){return *this=*this+b;}
	BigInt& operator -=(const BigInt &b){return *this=*this-b;}
	BigInt& operator *=(const BigInt &b){return *this=*this*b;}
	BigInt& operator /=(const BigInt &b){return *this=*this/b;}
	BigInt& operator %=(const BigInt &b){return *this=*this%b;}
	friend istream& operator >>(istream &in,BigInt &x){
		string s;if(!(in>>s))return in;
		int l=(s[0]=='+'||s[0]=='-');
		if(l==(int)s.size()){in.setstate(ios::failbit);return in;}
		for(int i=l;i<(int)s.size();i++)if(s[i]<'0'||s[i]>'9'){
			in.setstate(ios::failbit);return in;
		}
		BigInt y;y.sign=s[0]=='-'?-1:1;y.a.reserve((s.size()-l+8)/9);
		for(int r=s.size();r>l;r-=9){
			int v=0;for(int i=max(l,r-9);i<r;i++)v=v*10+s[i]-'0';
			y.a.push_back(v);
		}
		y.trim();x=move(y);return in;
	}
	friend ostream& operator <<(ostream &out,const BigInt &x){
		if(!x)return out<<'0';
		string s=x.sign<0?"-":"";s+=to_string(x.a.back());
		for(int i=(int)x.a.size()-2;i>=0;i--){
			string t=to_string(x.a[i]);s.append(9-t.size(),'0');s+=t;
		}
		return out<<s;
	}
};
// end for math/big-int.cpp
/////////////////////////
// !!!!! Division truncates toward zero; the divisor must be nonzero. !!!!
