__int128 floor_sum(ll l,ll r,ll a,ll b,ll c){
	assert(l>=0&&r>=0&&a>=0&&b>=0&&c>0);
	if(l>r)return 0;
	__int128 n=(__int128)r-l+1,A=a,B=(__int128)a*l+b,C=c,ans=0;
	while(1){
		ans+=n*(n-1)/2*(A/C)+n*(B/C);
		A%=C;B%=C;
		__int128 y=A*n+B;
		if(y<C)break;
		n=y/C;B=y%C;swap(A,C);
	}
	return ans;
}
