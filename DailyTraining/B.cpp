#include<bits/stdc++.h>
using namespace std;
const int MAXN=2e5+10,INF=0x3f3f3f3f;
void file(string s){freopen((s+".in").c_str(),"r",stdin),freopen((s+".out").c_str(),"w",stdout);}
int read(){
    int f=1,a=0;
    char ch=getchar();
    while(ch<'0'||ch>'9'){
        if(ch=='-')f=-f;
        ch=getchar();
    }
    while(ch>='0'&&ch<='9'){
        a=a*10+ch-'0';
        ch=getchar();
    }
    return a*f;
}

string s,t;
int l[MAXN],r[MAXN],pos[35];

void makel(){
    int k=0;
    memset(pos,-1,sizeof(pos));
    for(int i=0;i<s.length();++i){
        if(k<t.length()&&s[i]==t[k]){
            l[i]=k,pos[s[i]-'a']=k++;
        }else{
            l[i]=pos[s[i]-'a'];
        }
    }
}

void maker(){
    int k=t.length()-1;
    memset(pos,INF,sizeof(pos));
    for(int i=s.length()-1;i>=0;--i){
        if(k>=0&&s[i]==t[k]){
            r[i]=k,pos[s[i]-'a']=k--;
        }else{
            r[i]=pos[s[i]-'a'];
        }
    }
}

signed main(){
//  file("");
    cin>>s>>t;
    makel();
    maker();

    for (int i = 0; i < s.size(); i++) {
        cout << l[i] << " ";
    }
    cout << '\n';

    for (int i = 0; i < s.size(); i++) {
        cout << r[i] << " ";
    }
    cout << '\n';


    for(int i=0;i<s.length();++i){
        if(l[i]<r[i]){
            puts("No");
            return 0;
        }
    }
    puts("Yes");
    return 0;
}