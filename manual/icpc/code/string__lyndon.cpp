vector<pair<int,int>> lyndon(const string &s){
	assert(!s.empty()&&s[0]==' '&&s.size()<INT_MAX);
	int n=(int)s.size()-1;vector<pair<int,int>> ans;
	for(int i=1;i<=n;){
		int j=i+1,k=i;
		while(j<=n&&(unsigned char)s[k]<=(unsigned char)s[j]){
			k=s[k]==s[j]?k+1:i;j++;
		}
		int p=j-k;
		while(i<=k){ans.emplace_back(i,i+p-1);i+=p;}
	}
	return ans;
}
