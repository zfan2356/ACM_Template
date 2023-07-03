#include<bits/stdc++.h>
using namespace std;
struct edg//边
{
    int to,col,len;
};
struct zxs//主席树
{
    int lc,rc,tot,sum;
}tree[4000010];
vector<edg>g[100010];
int n,cnt,root[100010],dep[100010],fa[100010][20],dis[100010];
inline int read()//快读
{
    char c=getchar();
    int x=0;
    while(c<'0'||c>'9')
        c=getchar();
    while(c>='0'&&c<='9')
    {
        x=(x<<3)+(x<<1)+c-'0';
        c=getchar();
    }
    return x;
}
inline void write(int x)//快写
{
    int sta[10],tp=0;
    while(x)
    {
        sta[++tp]=x%10;
        x/=10;
    }
    while(tp)
        putchar(sta[tp--]+'0');
    putchar('\n');
}
int build(int l,int r)//建树
{
    int p=++cnt;
    tree[p].sum=tree[p].tot=0;
    if(l!=r)
    {
        int mid=(l+r)>>1;
        tree[p].lc=build(l,mid);
        tree[p].rc=build(mid+1,r);
    }
    return p;
}
int change(int q,int l,int r,int x,int y)//修改
{
    int p=++cnt;
    tree[p]=tree[q];
    tree[p].tot++;
    tree[p].sum+=y;
    if(l!=r)
    {
        int mid=(l+r)>>1;
        if(x>mid)
            tree[p].rc=change(tree[q].rc,mid+1,r,x,y);
        else
            tree[p].lc=change(tree[q].lc,l,mid,x,y);
    }
    return p;
}
int ask(int p,int l,int r,int x,int y)//查询
{
    if(l==r)
        return tree[p].tot*y-tree[p].sum;
    int mid=(l+r)>>1;
    if(x>mid)
        return ask(tree[p].rc,mid+1,r,x,y);
    return ask(tree[p].lc,l,mid,x,y);
}
void dfs(int x,int f)
{
    dep[x]=dep[f]+1;
    fa[x][0]=f;
    for(int i=1;i<19;i++)
        fa[x][i]=fa[fa[x][i-1]][i-1];
    for(int i=0;i<g[x].size();i++)
        if(g[x][i].to!=f)
        {
            dis[g[x][i].to]=dis[x]+g[x][i].len;
            root[g[x][i].to]=change(root[x],1,n,g[x][i].col,g[x][i].len);
            dfs(g[x][i].to,x);
        }
}
inline int LCA(int u,int v)//最近公共祖先
{
    if(dep[u]<dep[v])
        swap(u,v);
    for(int i=18;~i;i--)
        if(dep[fa[u][i]]>=dep[v])
            u=fa[u][i];
    if(u==v)
        return u;
    for(int i=18;~i;i--)
        if(fa[u][i]!=fa[v][i])
        {
            u=fa[u][i];
            v=fa[v][i];
        }
    return fa[u][0];
}

int main()
{
    int q,i,x,y,z,Z;
    n=read();
    q=read();
    for(i=1;i<n;i++)
    {
        x=read();
        y=read();
        z=read();
        Z=read();
        g[x].push_back((edg){y,z,Z});
        g[y].push_back((edg){x,z,Z});
    }
    root[1]=build(1,n);
    dfs(1,1);
    while(q--)
    {
        x=read();
        y=read();
        z=read();
        Z=read();
        int l=LCA(z,Z);
        write(dis[z]+dis[Z]-2*dis[l]+ask(root[z],1,n,x,y)+ask(root[Z],1,n,x,y)-2*ask(root[l],1,n,x,y));
    }
    return 0;
}