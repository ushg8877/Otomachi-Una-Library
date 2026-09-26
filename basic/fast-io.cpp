////////////////////////////////////////////////////////////////
//
// template for buffered IO / __int128 IO
//
// usage: input.set(stdin); output.set(stdout); // after freopen
//   input >> x; output << x << '\n'; output.flush();
//
////////////////////////////////////////////////////////////////
// cin/cout overloads also support signed / unsigned __int128.
using lll=__int128;
struct Reader{
FILE *f=stdin;
vector<char> buf;
size_t p=0,n=0;
bool ok=true;
Reader(){set(stdin);}
void set(FILE *_f=stdin,size_t sz=1<<16){
	assert(_f&&sz);
	f=_f;
	buf.resize(sz);
	p=n=0;
	ok=true;
}
int get(){
	if(p==n){n=fread(buf.data(),1,buf.size(),f);p=0;if(!n)return EOF;}
	return (unsigned char)buf[p++];
}
int next(){int c;do{c=get();}while(c!=EOF&&c<=32);return c;}
// Integer input must fit T; false / operator bool() reports EOF.
template<typename T>
bool read(T &x){
	int c=next();if(c==EOF)return ok=false;
	bool neg=c=='-';if(neg||c=='+')c=get();
	assert('0'<=c&&c<='9');
	// Accumulate in the unsigned type to handle the signed minimum.
	using U=conditional_t<(sizeof(T)>8),__uint128_t,unsigned long long>;
	U v=0;do{v=v*10+c-'0';c=get();}while('0'<=c&&c<='9');
	x=neg?T(U(0)-v):T(v);
	return ok=true;
}
bool read(char &x){int c=next();if(c==EOF)return ok=false;x=c;return ok=true;}
bool read(string &s){
	s.clear();
	int c=next();
	if(c==EOF)return ok=false;
	do{s.push_back(c);c=get();}while(c!=EOF&&c>32);
	return ok=true;
}
template<typename T>
Reader& operator>>(T &x){read(x);return *this;}
explicit operator bool()const{return ok;}
}input;
struct Writer{
FILE *f=stdout;
vector<char> buf;
size_t p=0;
Writer(){buf.resize(1<<16);}
Writer(const Writer&)=delete;
Writer& operator=(const Writer&)=delete;
// Flush output and fflush(stdout) before an interactive read.
void flush(){if(p){fwrite(buf.data(),1,p,f);p=0;}}
void set(FILE *_f=stdout,size_t sz=1<<16){
	assert(_f&&sz);
	flush();
	f=_f;
	buf.resize(sz);
}
void put(char c){if(p==buf.size())flush();buf[p++]=c;}
Writer& operator<<(char c){put(c);return *this;}
Writer& operator<<(const char *s){while(*s)put(*s++);return *this;}
Writer& operator<<(const string &s){for(char c:s)put(c);return *this;}
template<typename T>
Writer& operator<<(T x){
	using U=conditional_t<(sizeof(T)>8),__uint128_t,unsigned long long>;
	U v=x;if(x<0){put('-');v=U(0)-v;}
	char s[40];
	int n=0;
	do{s[n++]=char('0'+v%10);v/=10;}while(v);
	while(n)put(s[--n]);
	return *this;
}
~Writer(){flush();}
}output;
istream& operator>>(istream &o,__uint128_t &x){
	string s;if(!(o>>s))return o;
	size_t i=s[0]=='+';__uint128_t v=0,lim=~__uint128_t(0);
	if(i==s.size()){o.setstate(ios::failbit);return o;}
	for(;i<s.size();i++){
		unsigned d=s[i]-'0';
		if(d>9||v>(lim-d)/10){o.setstate(ios::failbit);return o;}v=v*10+d;
	}
	x=v;
	return o;
}
istream& operator>>(istream &o,lll &x){
	string s;if(!(o>>s))return o;
	bool neg=s[0]=='-';size_t i=neg||s[0]=='+';
	__uint128_t v=0,lim=(__uint128_t(1)<<127)-!neg;
	if(i==s.size()){o.setstate(ios::failbit);return o;}
	for(;i<s.size();i++){
		unsigned d=s[i]-'0';
		if(d>9||v>(lim-d)/10){o.setstate(ios::failbit);return o;}v=v*10+d;
	}
	x=neg?lll(__uint128_t(0)-v):lll(v);
	return o;
}
ostream& operator<<(ostream &o,__uint128_t x){
	char s[40];
	int n=0;
	do{s[n++]=char('0'+x%10);x/=10;}while(x);
	while(n)o.put(s[--n]);
	return o;
}
ostream& operator<<(ostream &o,lll x){
	__uint128_t v=x;if(x<0){o.put('-');v=__uint128_t(0)-v;}return o<<v;
}
// end for basic/fast-io.cpp
/////////////////////////
// !!!!! Do not mix input/output with cin/cout on the same stream. !!!!
