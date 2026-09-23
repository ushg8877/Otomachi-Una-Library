///////////////////////////////////////////////////////////////////////////////////
//
// template for FWT and bitwise convolution
//
// usage:
//   auto c=convolution_and(a,b); // convolution_or / convolution_xor
//   fwt_and(a); fwt_and(a,true); // forward / inverse; also or / xor
//   vector<int>, vector<ll> or vector<mint>; 0-indexed
//   convolution pads to a power of 2; empty input returns {}
//   transform size must be a positive power of 2
//   integer intermediates must fit the type; mint XOR needs inverse of 2
//
///////////////////////////////////////////////////////////////////////////////////
template<typename T>
void fwt_and(vector<T> &a,bool inverse=false){
	assert(!a.empty()&&a.size()<=INT_MAX&&!(a.size()&(a.size()-1)));
	int n=a.size();
	for(int k=1;k<n;k<<=1)
		for(int i=0;i<n;i+=k<<1)
			for(int j=0;j<k;j++){
				if(inverse)a[i+j]-=a[i+j+k];
				else a[i+j]+=a[i+j+k];
			}
}
template<typename T>
void fwt_or(vector<T> &a,bool inverse=false){
	assert(!a.empty()&&a.size()<=INT_MAX&&!(a.size()&(a.size()-1)));
	int n=a.size();
	for(int k=1;k<n;k<<=1)
		for(int i=0;i<n;i+=k<<1)
			for(int j=0;j<k;j++){
				if(inverse)a[i+j+k]-=a[i+j];
				else a[i+j+k]+=a[i+j];
			}
}
template<typename T>
void fwt_xor(vector<T> &a,bool inverse=false){
	assert(!a.empty()&&a.size()<=INT_MAX&&!(a.size()&(a.size()-1)));
	int n=a.size();T inv2=1;
	if constexpr(!is_integral_v<T>)if(inverse)inv2=T(1)/T(2);
	for(int k=1;k<n;k<<=1)
		for(int i=0;i<n;i+=k<<1)
			for(int j=0;j<k;j++){
				if constexpr(is_integral_v<T>){
					using W=conditional_t<(sizeof(T)<=4),long long,__int128>;
					W x=a[i+j],y=a[i+j+k],u=x+y,v=x-y;
					if(inverse){assert(u%2==0&&v%2==0);u/=2;v/=2;}
					assert(numeric_limits<T>::lowest()<=u&&u<=numeric_limits<T>::max());
					assert(numeric_limits<T>::lowest()<=v&&v<=numeric_limits<T>::max());
					a[i+j]=T(u);a[i+j+k]=T(v);
				}else{
					T x=a[i+j],y=a[i+j+k];
					a[i+j]=x+y;a[i+j+k]=x-y;
					if(inverse){a[i+j]*=inv2;a[i+j+k]*=inv2;}
				}
			}
}

template<typename T>
vector<T> convolution_and(vector<T> a,vector<T> b){
	if(a.empty()||b.empty())return {};
	assert(max(a.size(),b.size())<=(1u<<30));
	int n=1;while(n<(int)max(a.size(),b.size()))n<<=1;
	a.resize(n);b.resize(n);
	fwt_and(a);fwt_and(b);
	for(int i=0;i<n;i++)a[i]*=b[i];
	fwt_and(a,true);return a;
}

template<typename T>
vector<T> convolution_or(vector<T> a,vector<T> b){
	if(a.empty()||b.empty())return {};
	assert(max(a.size(),b.size())<=(1u<<30));
	int n=1;while(n<(int)max(a.size(),b.size()))n<<=1;
	a.resize(n);b.resize(n);
	fwt_or(a);fwt_or(b);
	for(int i=0;i<n;i++)a[i]*=b[i];
	fwt_or(a,true);return a;
}

template<typename T>
vector<T> convolution_xor(vector<T> a,vector<T> b){
	if(a.empty()||b.empty())return {};
	assert(max(a.size(),b.size())<=(1u<<30));
	int n=1;while(n<(int)max(a.size(),b.size()))n<<=1;
	a.resize(n);b.resize(n);
	fwt_xor(a);fwt_xor(b);
	for(int i=0;i<n;i++)a[i]*=b[i];
	fwt_xor(a,true);return a;
}
// end for polynomial/fwt.cpp
/////////////////////////
// !!!!! Transform length must be a power of two; integer intermediates must fit T; mint XOR requires 2 to be invertible. !!!!
