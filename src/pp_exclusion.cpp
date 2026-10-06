// Exhaustive irreducible-support search. See README.txt for the proof of
// completeness, the coset-leader pruning rule, and the shard partition.
#include "common.hpp"

struct Exact {
    BinaryCode c;int P,W,root,begin,end,maxdeg=0,maxrow=0,sw=0,bad=0;
    double seconds;uint64_t cap,nodes=0,coset_prunes=0,pack_prunes=0;
    bool found=false,incomplete=false,use_coset=true;Word selected{};
    bool hashed=false;int shard=0,shards=1,split_weight=1,roots_completed=0;
    uint64_t frontier_seen=0,frontier_owned=0,frontier_digest=0,common_nodes=0;
    std::vector<uint64_t> available,syndrome;
    std::vector<unsigned char> forbidden;
    std::vector<int> overlap;
    std::vector<std::vector<std::pair<int,uint64_t>>> edges;
    std::array<int,MAXM> index;
    std::chrono::steady_clock::time_point start;
    Exact(const std::string &path,int p,int w,double sec,int r,int b,int e,uint64_t cap,bool prune):
      c(path),P(p),W(w),root(r),begin(b),end(e),seconds(sec),cap(cap),use_coset(prune){
        if(W<1||W>64||root<0||root>=c.n/P||begin<0||end<=begin)throw std::runtime_error("invalid exact search parameters");
        c.validate_translation(P);edges.resize(c.n);available.resize(c.m);syndrome.resize((c.m+63)/64);forbidden.resize(c.n);overlap.resize(c.g);
        for(int j=0;j<c.m;++j){int d=c.hr[j].size();if(d>64)throw std::runtime_error("row degree exceeds 64");maxrow=std::max(maxrow,d);available[j]=d==64?~0ULL:((1ULL<<d)-1);
            for(int k=0;k<d;++k)edges[c.hr[j][k]].push_back({j,1ULL<<k});}
        for(auto &a:c.vc)maxdeg=std::max(maxdeg,int(a.size()));
    }
    bool unsatisfied(int j)const{return (syndrome[j/64]>>(j%64))&1;}
    void block(int v){forbidden[v]=1;for(auto [j,b]:edges[v])available[j]^=b;}
    void unblock(int v){forbidden[v]=0;for(auto [j,b]:edges[v])available[j]^=b;}
    void flip(int v,int delta){
        toggle(selected,v);for(int j:c.vc[v]){sw+=unsatisfied(j)?-1:1;syndrome[j/64]^=1ULL<<(j%64);}
        if(use_coset)for(int j:c.vg[v]){bad-=(2*overlap[j]>int(c.gr[j].size()));overlap[j]+=delta;bad+=(2*overlap[j]>int(c.gr[j].size()));}
    }
    bool independent(uint64_t a,int need,const uint64_t *adj,int &steps){
        if(need<=0)return true;if(__builtin_popcountll(a)<need||++steps>120)return false;
        int v=-1,mind=100;uint64_t tmp=a;
        while(tmp){int j=__builtin_ctzll(tmp);tmp&=tmp-1;int d=__builtin_popcountll(adj[j]&a);if(d<mind){v=j;mind=d;if(d<=1)break;}}
        return independent(a&~adj[v],need-1,adj,steps)||independent(a&~(1ULL<<v),need,adj,steps);
    }
    static uint64_t mix(uint64_t x){
        x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
        x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);
    }
    void dfs(int wt,uint64_t key=0){
        ++nodes;
        if(nodes>cap||((nodes==1||(nodes&16383)==0)&&(stopped||elapsed(start)>=seconds))){incomplete=true;return;}
        // All shards traverse the same short prefix tree. At this frontier,
        // each subtree belongs to exactly one residue class; no pruning rule
        // depends on the shard. See README.txt for the coverage argument.
        if(hashed&&wt<=split_weight)++common_nodes;
        if(hashed&&wt==split_weight){
            ++frontier_seen;frontier_digest^=mix(key);if(key%uint64_t(shards)!=uint64_t(shard))return;++frontier_owned;
        }
        if(bad){++coset_prunes;return;}
        if(!sw){if(c.basis.nontrivial(selected)){
            if(!c.kernel(selected))throw std::runtime_error("internal syndrome mismatch");found=true;
            std::cout<<"{\"event\":\"witness\",\"weight\":"<<wt<<",\"support\":";print_support(selected,c.n);std::cout<<"}\n"<<std::flush;
        }return;}
        int budget=W-wt;if(budget<=0||sw>budget*maxdeg)return;
        // Bucket the active checks by available degree, with no heap allocation
        // and no full scan of the check matrix at each search node.
        int ids[MAXM],next[MAXM],head[65],tail[65],count=0,best=-1,minavail=65;
        std::fill(head,head+65,-1);std::fill(tail,tail+65,-1);
        for(size_t k=0;k<syndrome.size();++k){uint64_t a=syndrome[k];while(a){int j=64*k+__builtin_ctzll(a);a&=a-1;
            int d=__builtin_popcountll(available[j]);if(!d)return;
            ids[count]=j;next[count]=-1;if(head[d]<0)head[d]=count;else next[tail[d]]=count;tail[d]=count;++count;
            if(d<minavail){minavail=d;best=j;}
        }}
        Word used{};int packing=0;
        for(int d=minavail;d<=maxrow;++d)for(int k=head[d];k>=0;k=next[k]){
            int j=ids[k];uint64_t a=available[j];bool intersects=false;
            while(a){int b=__builtin_ctzll(a);a&=a-1;int v=c.hr[j][b];if(bit(used,v)){intersects=true;break;}}
            if(!intersects){if(++packing>budget){++pack_prunes;return;}a=available[j];while(a){int b=__builtin_ctzll(a);a&=a-1;int v=c.hr[j][b];used[v/64]|=1ULL<<(v%64);}}
        }
        if(sw>budget&&sw<=63&&packing+3>budget){
            uint64_t adj[63];for(int i=0;i<count;++i)index[ids[i]]=i;
            for(int i=0;i<count;++i){adj[i]=1ULL<<i;uint64_t a=available[ids[i]];
                while(a){int b=__builtin_ctzll(a);a&=a-1;for(int j:c.vc[c.hr[ids[i]][b]])if(unsatisfied(j))adj[i]|=1ULL<<index[j];}}
            int steps=0;if(independent((1ULL<<count)-1,budget+1,adj,steps)){++pack_prunes;return;}
        }
        struct Option{int v,score;};Option options[64];int size=0;uint64_t a=available[best];
        while(a){int b=__builtin_ctzll(a);a&=a-1;int v=c.hr[best][b],score=0;for(int j:c.vc[v])score+=unsatisfied(j);options[size++]={v,score};}
        std::sort(options,options+size,[](Option a,Option b){return a.score!=b.score?a.score>b.score:a.v<b.v;});
        int blocked=0;
        for(int j=0;j<size;++j){if(wt==1&&j>=end)break;int v=options[j].v;block(v);++blocked;
            if(wt!=1||j>=begin){flip(v,1);dfs(wt+1,hashed&&wt<split_weight?mix(key^uint64_t(v)):key);flip(v,-1);if(found||incomplete)break;}}
        for(int j=0;j<blocked;++j)unblock(options[j].v);
    }
    int run(){
        start=std::chrono::steady_clock::now();int first=hashed?0:root,last=hashed?c.n/P:root+1;
        for(int r=first;r<last;++r){
            selected={};sw=bad=0;std::fill(forbidden.begin(),forbidden.end(),0);std::fill(syndrome.begin(),syndrome.end(),0);std::fill(overlap.begin(),overlap.end(),0);
            for(int j=0;j<c.m;++j){int d=c.hr[j].size();available[j]=d==64?~0ULL:((1ULL<<d)-1);}
            int v=r*P;for(int j=0;j<=v;++j)block(j);flip(v,1);dfs(1,mix(uint64_t(r)));flip(v,-1);
            if(found||incomplete)break;++roots_completed;
        }
        std::cout<<"{\"event\":\"result\",\"status\":\""<<(found?"witness":incomplete?"incomplete":"excluded")
          <<"\",\"W\":"<<W<<",\"root\":"<<root<<",\"branch_begin\":"<<begin<<",\"branch_end\":"<<end
          <<",\"coset_pruning\":"<<(use_coset?"true":"false")<<",\"nodes\":"<<nodes<<",\"coset_prunes\":"<<coset_prunes
          <<",\"packing_prunes\":"<<pack_prunes<<",\"seconds\":"<<elapsed(start)
          <<",\"partition\":\""<<(hashed?"hash_prefix_v1":"root_branch_v1")<<"\",\"shard\":"<<shard<<",\"shards\":"<<shards<<",\"split_weight\":"<<split_weight
          <<",\"roots_completed\":"<<roots_completed<<",\"frontier_seen\":"<<frontier_seen<<",\"frontier_owned\":"<<frontier_owned
          <<",\"frontier_digest\":\""<<frontier_digest<<"\",\"common_nodes\":"<<common_nodes<<"}\n"<<std::flush;
        return found?10:incomplete?75:0;
    }
};
int main(int argc,char **argv){try{
    if(argc<9){std::cerr<<"exact matrix P W seconds root branch_begin branch_end nodecap [no-coset]\n"
      <<"exact matrix P W seconds --shard id count split_weight nodecap [no-coset]\n";return 2;}
    std::signal(SIGTERM,stop_handler);std::signal(SIGINT,stop_handler);
    int p=std::stoi(argv[2]);if(p<1)throw std::runtime_error("invalid P");
    bool hashed=std::string(argv[5])=="--shard";
    if(hashed&&argc<10)throw std::runtime_error("missing shard parameters");
    Exact search(argv[1],p,std::stoi(argv[3]),std::stod(argv[4]),hashed?0:std::stoi(argv[5]),hashed?0:std::stoi(argv[6]),hashed?64:std::stoi(argv[7]),std::stoull(argv[hashed?9:8]),argc<(hashed?11:10));
    if(hashed){search.hashed=true;search.shard=std::stoi(argv[6]);search.shards=std::stoi(argv[7]);search.split_weight=std::stoi(argv[8]);
        if(search.shards<1||search.shard<0||search.shard>=search.shards||search.split_weight<1||search.split_weight>64)throw std::runtime_error("invalid shard partition");}
    return search.run();
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 2;}}
