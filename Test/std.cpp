#include<bits/stdc++.h>
using namespace std ;
#define rep( i, s, t ) for(int i = s; i <= t; ++ i )
#define re register
#define int long long
const int P = 1004535809 ;
const int Gi = 334845270 ;
const int G = 3 ;
int gi() {
    char cc = getchar() ; int cn = 0, flus = 1 ;
    while( cc < '0' || cc > '9' ) {  if( cc == '-' ) flus = - flus ; cc = getchar() ; }
    while( cc >= '0' && cc <= '9' )  cn = ( cn * 10 + cc - '0' ) % P, cc = getchar() ;
    return cn * flus ;
}
const int N = 1e5 + 5 ;
int n, m, A[N << 2], B[N << 2], limit, L, Inv, R[N << 2] ;
int fpow( int x, int k ) {
    int ans = 1, base = x ;
    while( k ) {
        if( k & 1 ) ans = ( ans * base ) % P ;
        base = ( base * base ) % P, k >>= 1 ;
    } return ans ;
}
void init( int x ) {
    limit = 1, L = 0 ;
    while( limit <= x ) limit <<= 1, ++ L ;
    rep( i, 0, limit ) R[i] = ( R[i >> 1] >> 1 ) | ( ( i & 1 ) << ( L - 1 ) ) ;
    Inv = fpow( limit, P - 2 ) ;
}
void NTT( int *a, int type ) {
    for(int i = 0; i < limit; ++ i ) if( R[i] > i ) swap( a[i], a[R[i]] ) ;
    for(int k = 1; k < limit; k <<= 1 ) {
        int d = fpow( ( type == 1 ) ? G : Gi, ( P - 1 ) / ( k << 1 ) ) ;
        for(int i = 0; i < limit; i += ( k << 1 ) )
            for(int j = i, g = 1 ; j < i + k; ++ j, g = ( g * d ) % P ) {
                int Nx = a[j], Ny = ( a[j + k] * g ) % P ;
                a[j] = ( Nx + Ny ) % P, a[j + k] = ( Nx - Ny + P ) % P ;
            }
    }
    if( type != 1 ) rep( i, 0, limit ) a[i] = a[i] * Inv % P ;
}
signed main()
{
    n = gi(), m = gi() ; int type = gi() ;
    rep( i, 1, n ) A[i - 1] = gi() ; B[0] = 1 ;
    if( type == 0 ) rep( i, 1, n ) B[i] = B[i - 1] * ( m + i - 1 ) % P * fpow( i, P - 2 ) % P ;
    if( type == 1 ) rep( i, 1, n ) B[i] = ( -B[i - 1] * ( m - i + 1 + P ) % P * fpow( i, P - 2 ) % P + P ) % P ;
    init( n + n ), NTT( A, 1 ), NTT( B, 1 ) ;
    rep( i, 0, limit ) A[i] = A[i] * B[i] % P ;
    NTT( A, -1 ) ; rep( i, 1, n ) printf("%lld ", A[i - 1] ) ;
    return 0 ;
}