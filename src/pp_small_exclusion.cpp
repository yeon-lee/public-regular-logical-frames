// Portable copy of ../full_P_PP_scaling/certificates/cluster.cpp.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
using namespace std;
const int MAXN=4096,MAXM=2048,MW=64;
using Word=array<uint64_t,MW>;
int n,m,g,nw,P,limit,maxdeg,side,top_start=0,top_stop=MAXN; vector<vector<int>> hr,vc;vector<Word> Gb;vector<int> Gp;
array<unsigned char,MAXN> forbidden{};array<unsigned char,MAXM> syn{};Word selected{};int sw=0;
uint64_t nodes=0,nodecap;bool aborted=false,found=false; double seconds;chrono::steady_clock::time_point start;
inline void flip(int v){selected[v>>6]^=1ULL<<(v&63);for(int c:vc[v]){sw+=syn[c]?-1:1;syn[c]^=1;}}
bool nontriv(){Word a=selected;for(int i=0;i<(int)Gb.size();i++)if((a[Gp[i]>>6]>>(Gp[i]&63))&1)for(int w=0;w<nw;w++)a[w]^=Gb[i][w];for(int w=0;w<nw;w++)if(a[w])return true;return false;}
bool independent_witness(uint64_t avail,int need,const vector<uint64_t>&adj,int &steps){
 if(need<=0)return true;if(__builtin_popcountll(avail)<need||++steps>120)return false;
 int vertex=-1,min_degree=100;uint64_t temp=avail;
 while(temp){int v=__builtin_ctzll(temp);temp&=temp-1;int d=__builtin_popcountll(adj[v]&avail);if(d<min_degree){min_degree=d;vertex=v;if(d<=1)break;}}
 if(independent_witness(avail&~adj[vertex],need-1,adj,steps))return true;
 return independent_witness(avail&~(1ULL<<vertex),need,adj,steps);
}
void dfs(int wt){
 nodes++;if((nodes&65535)==0){double elapsed=chrono::duration<double>(chrono::steady_clock::now()-start).count();if(nodes>=nodecap||elapsed>=seconds){aborted=true;return;}}
 if(sw==0){if(nontriv()){found=true;cout<<"FOUND weight "<<wt<<" support";for(int v=0;v<n;v++)if((selected[v>>6]>>(v&63))&1)cout<<" "<<v;cout<<"\n"<<flush;}return;}
 int budget=limit-wt;if(budget<=0||sw>budget*maxdeg)return;
 vector<pair<int,int>> checks;checks.reserve(sw);int bestcheck=-1,minavail=MAXN;
 for(int c=0;c<m;c++)if(syn[c]){int count=0;for(int v:hr[c])if(!forbidden[v])count++;if(count==0)return;checks.emplace_back(count,c);if(count<minavail){minavail=count;bestcheck=c;}}
 // Greedy packing of unsatisfied checks with disjoint available neighborhoods.
 sort(checks.begin(),checks.end());Word used{};int packing=0;
 for(auto [cnt,c]:checks){bool overlap=false;for(int v:hr[c])if(!forbidden[v]&&((used[v>>6]>>(v&63))&1)){overlap=true;break;}
  if(!overlap){packing++;if(packing>budget)return;for(int v:hr[c])if(!forbidden[v])used[v>>6]|=1ULL<<(v&63);}}
 if(sw>budget && sw<=63 && packing+3>budget){
  vector<uint64_t> adj(checks.size(),0);array<int,MAXM> idx;fill(idx.begin(),idx.end(),-1);
  for(int i=0;i<(int)checks.size();i++)idx[checks[i].second]=i;
  for(int i=0;i<(int)checks.size();i++){adj[i]=1ULL<<i;for(int v:hr[checks[i].second])if(!forbidden[v])for(int c:vc[v])if(idx[c]>=0)adj[i]|=1ULL<<idx[c];}
  int steps=0;if(independent_witness((1ULL<<checks.size())-1,budget+1,adj,steps))return;
 }
 vector<int> options;for(int v:hr[bestcheck])if(!forbidden[v])options.push_back(v);
 // Put variables repairing most syndrome checks first. This changes no completeness.
 sort(options.begin(),options.end(),[](int a,int b){int ca=0,cb=0;for(int c:vc[a])ca+=syn[c];for(int c:vc[b])cb+=syn[c];return ca!=cb ? ca>cb : a<b;});
 int usedcount=0;
 if(wt==1){cout<<"top_options";for(int v:options)cout<<" "<<v;cout<<" branches "<<top_start<<" "<<top_stop<<"\n"<<flush;}
 for(int j=0;j<(int)options.size();j++){
  if(wt==1 && j>=top_stop)break;
  int v=options[j];forbidden[v]=1;usedcount++;
  if(wt!=1 || j>=top_start){flip(v);dfs(wt+1);flip(v);if(found||aborted)break;}
 }
 for(int j=0;j<usedcount;j++)forbidden[options[j]]=0;
}
int main(int argc,char**argv){
 if(argc<6){cerr<<"cluster input P W seconds nodecap [root_start root_stop]\n";return 2;}ifstream f(argv[1]);f>>n>>m>>g;nw=(n+63)/64;P=stoi(argv[2]);limit=stoi(argv[3]);seconds=stod(argv[4]);nodecap=stoull(argv[5]);
 if(n>MAXN||m>MAXM||nw>MW)return 3;hr.resize(m);vc.resize(n);
 for(int c=0;c<m;c++){int w;f>>w;hr[c].resize(w);for(int &v:hr[c]){f>>v;vc[v].push_back(c);}}
 for(auto v:vc)maxdeg=max(maxdeg,(int)v.size());
 for(int i=0;i<g;i++){Word a{};int w,v;f>>w;while(w--){f>>v;a[v>>6]^=1ULL<<(v&63);}for(int j=0;j<(int)Gb.size();j++)if((a[Gp[j]>>6]>>(Gp[j]&63))&1)for(int w=0;w<nw;w++)a[w]^=Gb[j][w];int p=0;while(p<n&&!((a[p>>6]>>(p&63))&1))p++;if(p<n){Gp.push_back(p);Gb.push_back(a);}}
 int rootstart=argc>6?stoi(argv[6]):0,rootstop=argc>7?stoi(argv[7]):n/P;
 if(argc>8)top_start=stoi(argv[8]);if(argc>9)top_stop=stoi(argv[9]);
 start=chrono::steady_clock::now();
 for(int b=rootstart;b<rootstop;b++){int root=b*P;fill(forbidden.begin(),forbidden.end(),0);fill(syn.begin(),syn.end(),0);selected={};sw=0;for(int v=0;v<=root;v++)forbidden[v]=1;flip(root);dfs(1);flip(root);
  cout<<"root "<<b<<" nodes "<<nodes<<" elapsed "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n"<<flush;if(found||aborted)break;}
 cout<<(found?"COUNTEREXAMPLE":aborted?"INCOMPLETE":"CERTIFIED_NO_LOGICAL")<<" W "<<limit<<" roots "<<rootstart<<" "<<rootstop<<" nodes "<<nodes<<" branches "<<top_start<<" "<<top_stop<<"\n"<<flush;
}
