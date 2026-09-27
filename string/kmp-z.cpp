////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/string/kmp-z.cpp
//
// usage: string s = " ababa"; auto b = border(s); auto z = z_function(s);
//   1-indexed with leading space
//
////////////////////////////////////////////////////////////////
vector<int> border(string s){
	// s is a string with length n+1,
	// this function return a vector with length n+1
	// indication the border array of s
	assert(!s.empty()&&s.size()<=INT_MAX);
	int n=(int)s.length()-1;
	assert(!s.empty()&&s[0]==' ');
	vector<int> a(n+1);
	for(int i=2,j=0;i<=n;i++){
		while(j&&s[j+1]!=s[i]) j=a[j];
		if(s[j+1]==s[i]) j++;
		a[i]=j;
	}
	return a;
}
vector<int> z_function(string s){
	// s is a string with length n+1,
	// this function return a vector with length n+1
	// indication the a[i]=lcp(s[i:],s[])
	assert(!s.empty()&&s.size()<=INT_MAX);
	int n=(int)s.length()-1;
	assert(!s.empty()&&s[0]==' ');
	vector<int> z(n+1);
	for(int i=2,l=1,r=0;i<=n;i++){
		z[i]=r<i?0:min(r-i+1,z[i-l+1]);
		while(i+z[i]<=n&&s[1+z[i]]==s[i+z[i]]) z[i]++;
		if(i+z[i]>r) l=i,r=i+z[i]-1;
	}
	return z;
}
// end for string/kmp-z.cpp
/////////////////////////
