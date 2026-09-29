struct Runs{
private:
struct LCP{
	int n=0;
	vector<int> rk;
	Linear_RMQ<int> st;
	void build(const string &s){
		SA t;t.build(s,256);n=t.n;rk=move(t.rk);
		vector<int> h(n+1);
		for(int i=1,k=0;i<=n;i++){
			int r=rk[i];if(r==1){k=0;continue;}
			int j=t.sa[r-1];
			while(i+k<=n&&j+k<=n&&s[i+k]==s[j+k])k++;
			h[r]=k;if(k)k--;
		}
		st.build(h);
	}
	int ask(int i,int j)const{
		if(i>n||j>n)return 0;
		if(i==j)return n-i+1;
		int l=rk[i],r=rk[j];if(l>r)swap(l,r);
		return st.ask(l+1,r);
	}
};
public:
vector<array<int,3>> build(const string &s)const{
	assert(!s.empty()&&s[0]==' '&&s.size()<(INT_MAX-1024)/2);
	int n=(int)s.size()-1;
	for(int i=1;i<=n;i++)assert((unsigned char)s[i]>0);
	vector<array<int,3>> ans;if(n<2)return ans;
	LCP a,b;a.build(s);
	string t=s;reverse(t.begin()+1,t.end());b.build(t);
	vector<int> st;st.reserve(n);
	for(int op=0;op<2;op++){
		st.clear();
		for(int i=n;i>=1;i--){
			while(!st.empty()){
				int j=st.back(),k=a.ask(i,j);
				// Reversing the alphabet does not reverse the prefix rule.
				if(j+k>n)break;
				bool less=(unsigned char)s[i+k]<(unsigned char)s[j+k];
				if(less==bool(op))break;
				st.pop_back();
			}
			if(!st.empty()){
				int j=st.back(),p=j-i;
				int l=b.ask(n-i+2,n-j+2),r=a.ask(i,j);
				if(l+r>=p)ans.push_back({i-l,j+r-1,p});
			}
			st.push_back(i);
		}
	}
	sort(ans.begin(),ans.end());ans.erase(unique(ans.begin(),ans.end()),
		ans.end());
	return ans;
}
};
