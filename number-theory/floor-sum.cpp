////////////////////////////////////////////////////////////////
//
// template for Floor Sum, O(log c)
//
// usage:
//   floor_sum(l,r,a,b,c); // sum floor((a*i+b)/c), l<=i<=r
//   0<=l,r,a,b; c>0; return __int128, answer must fit
//
////////////////////////////////////////////////////////////////
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
// end for number-theory/floor-sum.cpp
/////////////////////////
// !!!!! Requires nonnegative l,r,a,b and positive c; result is __int128, basic/fast-io.cpp provides decimal IO. !!!!
