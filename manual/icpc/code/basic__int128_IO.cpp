using lll=__int128;
istream& operator >>(istream &o,lll &x){
	string s;if(!(o>>s))return o;
	bool neg=s[0]=='-';size_t i=(neg||s[0]=='+');
	assert(i<s.size());__uint128_t v=0,lim=((__uint128_t)1<<127)-!neg;
	for(;i<s.size();i++){
		assert('0'<=s[i]&&s[i]<='9');unsigned d=s[i]-'0';
		assert(v<=(lim-d)/10);v=v*10+d;
	}
	x=neg?(lll)(__uint128_t(0)-v):(lll)v;
	return o;
}
ostream& operator <<(ostream &o,lll x){
	array<int,40> a;int n=0;
	if(x==0){o<<0;return o;}
	__uint128_t v=x;
	if(x<0){o<<"-";v=__uint128_t(0)-v;}
	while(v) a[++n]=v%10,v/=10;
	for(int i=n;i>=1;i--) o<<a[i];
	return o;
}
