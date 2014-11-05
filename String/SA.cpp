// http://poj.org/problem?id=3693
/*
最大重复次数子串（字典序最小），如：ccabababc， ans ＝ ababab
*/
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <string>
#include <set>
#include <map>
using namespace std;
typedef long long lint;

const int MAX_N = 1e5 + 1000;
int two[MAX_N];

int sa[MAX_N], sar[MAX_N],
	rk[MAX_N], rkr[MAX_N];
int h[MAX_N], hr[MAX_N];
int dp[22][MAX_N], dpr[22][MAX_N], dprk[22][MAX_N];


const int MAX_M = 128;
int str[MAX_N];

int wa[MAX_N], wb[MAX_N], wv[MAX_N], ww[MAX_N];
int cmp(int *r, int a, int b, int l){
	return r[a] == r[b] && r[a + l] == r[b + l];
}
void da(int *r, int *sa, int n, int m){
	int i, j, p, *x = wa, *y = wb, *t;
	for(i = 0; i < m; i++) ww[i] = 0;
	for(i = 0; i < n; i++) ww[x[i] = r[i]]++;
	for(i = 1; i < m; i++) ww[i] += ww[i - 1];
	for(i = n - 1; i >= 0; i--) sa[--ww[x[i]]] = i;

	for(j = 1, p = 1; p < n; j *= 2, m = p){
		for(p = 0, i = n - j; i < n; i++) y[p++] = i;
		for(i = 0; i < n; i++) if(sa[i] >= j) y[p++] = sa[i] - j;
		for(i = 0; i < n; i++) wv[i] = x[y[i]];
		for(i = 0; i < m; i++) ww[i] = 0;
		for(i = 0; i < n; i++) ww[wv[i]]++;
		for(i = 1; i < m; i++) ww[i] += ww[i - 1];
		for(i = n - 1;  i >= 0; i--) sa[--ww[wv[i]]] = y[i];
		for(swap(x, y), p = 1, x[sa[0]] = 0, i = 1; i < n; i++){
			x[sa[i]] = cmp(y, sa[i - 1], sa[i], j) ? p - 1: p++;
		}
	}
}

void calheight(int * r, int * sa, int n, int *rk, int *h){
	int i, j, k = 0;
	for(i = 1; i <= n; i++) rk[sa[i]] = i;
	for(i = 0; i < n; h[rk[i++]] = k)
	for(k? k--: 0, j = sa[rk[i] - 1]; r[i + k] == r[j + k]; k++);
}

void gao(int * ch, int n, int *sa, int *rk, int *h,int dp[22][MAX_N], bool f = 0){
	da(ch, sa, n, 128);
	calheight(ch, sa, n - 1, rk, h);

	for(int i = 1; i < n; i++)
		dp[0][i - 1] = h[i];

	for(int s = 1; s < 22; s++){
		for(int i = 0; i < n - 1; i++){
			if((i + (1<<s)) > n - 1) break;
			dp[s][i] = min(dp[s - 1][i], dp[s - 1][i + (1<<(s - 1))]);
		}
	}

	if(f == 0) return;
	for(int i = 0; i < n; i++) dprk[0][i] = i;
	for(int s = 1; s < 22; s++){
		for(int i = 0; i < n; i++){
			if((i + (1<<s)) > n) break;
			int l = dprk[s - 1][i], r = dprk[s - 1][i + (1<<(s - 1))];
			dprk[s][i] = rk[l] < rk[r] ? l: r;
		}
	}
}

inline int lca(int a, int b, int rk[], int dp[22][MAX_N]){ // [a, b)
	a = rk[a], b = rk[b];
	if(a > b) swap(a, b);
	int s = two[b - a];
	return min(dp[s][a], dp[s][b - (1<<s)]);
}

inline int getminrk(int a, int b){
	if(a > b) swap(a, b); b++;
	int s = two[b - a];
	int l = dprk[s][a], r = dprk[s][b - (1<<s)];
	return rk[l] < rk[r] ? l: r;
}

char ch[MAX_N];

int main() {
	two[0] = -1;
	for(int i = 1; i < MAX_N; i++){
		two[i] = two[i - 1] + ((!(i & (i - 1))) ? 1: 0);
	}

	int icase = 1;
	while(scanf("%s", ch) == 1){
		if(ch[0] == '#') break;
		int n = strlen(ch); 
		for(int i = 0; i < n; i++) str[i] = ch[n - 1 - i];
		str[n] = 0;
		gao(str, n + 1, sar, rkr, hr, dpr);

		for(int i = 0; i < n; i++) str[i] = ch[i];
		str[n] = 0;
		gao(str, n + 1, sa, rk, h, dp, 1);


		int ansp = 0, ansrep = 1, anslen = 1;
		for(int i = 0; i < n; i++){
			if(ch[i] < ch[ansp]) ansp = i;
		}

		for(int i = 1; i <= n; i++){
			for(int j = i; j < n; j += i){
				// pos j - i && j
				//cerr<<"i = "<<i<<" j = "<<j<<endl;
				int right = lca(j - i, j, rk, dp),
					left = lca(n - 1 - j + i, n - 1 - j, rkr, dpr);
				if(left == 0 || right == 0) continue;

				int len = right + left - 1;
				//[j - i - left + 1, j + right - len] can be
				int rep = (len + i) / i;
				len = i * rep;
				if(rep < ansrep) continue;

				if(j + right - len < j - i - left + 1) continue;
				int p = getminrk(j - i - left + 1, j + right - len);

				if(rep > ansrep || rep == ansrep && rk[ansp] > rk[p]){ 
					ansrep = rep;
					ansp = p;
					anslen = len;
				}
			}
		}

		string s = ch;
		s = s.substr(ansp, anslen);
		printf("Case %d: %s\n", icase++, s.c_str());

	}
	return 0;
}

