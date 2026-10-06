#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int MAXN=4096, MAXM=2048, MW=MAXN/64;
using Word=std::array<uint64_t,MW>;
inline int bit(const Word &v,int j){return (v[j/64]>>(j%64))&1;}
inline void toggle(Word &v,int j){v[j/64]^=uint64_t(1)<<(j%64);}
inline int weight(const Word &v,int words){int w=0;for(int j=0;j<words;++j)w+=__builtin_popcountll(v[j]);return w;}
inline void xorin(Word &a,const Word &b,int words){for(int j=0;j<words;++j)a[j]^=b[j];}
inline double elapsed(std::chrono::steady_clock::time_point t){return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();}
inline volatile std::sig_atomic_t stopped=0;
inline void stop_handler(int){stopped=1;}
struct Basis {
    int n,words;std::vector<Word> rows;std::vector<int> pivots;
    explicit Basis(int n):n(n),words((n+63)/64){}
    Word reduce(Word a)const{for(size_t i=0;i<rows.size();++i)if(bit(a,pivots[i]))xorin(a,rows[i],words);return a;}
    bool nontrivial(Word a)const{return weight(reduce(a),words)!=0;}
    void add(Word a){a=reduce(a);for(int j=0;j<n;++j)if(bit(a,j)){rows.push_back(a);pivots.push_back(j);break;}}
};
inline Word readword(std::istream &f,int n){
    Word v{};int w,j;if(!(f>>w)||w<0||w>n)throw std::runtime_error("invalid row length");
    while(w--){if(!(f>>j)||j<0||j>=n||bit(v,j))throw std::runtime_error("invalid/duplicate coordinate");toggle(v,j);}return v;
}
struct BinaryCode {
    int n,m,g,words;std::vector<Word> H,G;std::vector<std::vector<int>> hr,vc,gr,vg;Basis basis;
    explicit BinaryCode(const std::string &path):n(0),m(0),g(0),words(0),basis(0){
        std::ifstream f(path);if(!(f>>n>>m>>g)||n<=0||n>MAXN||m<0||m>MAXM||g<0||g>MAXM)throw std::runtime_error("invalid matrix header");
        words=(n+63)/64;basis=Basis(n);hr.resize(m);vc.resize(n);gr.resize(g);vg.resize(n);
        for(int c=0;c<m;++c){auto v=readword(f,n);H.push_back(v);for(int j=0;j<n;++j)if(bit(v,j)){hr[c].push_back(j);vc[j].push_back(c);}}
        for(int c=0;c<g;++c){auto v=readword(f,n);G.push_back(v);basis.add(v);for(int j=0;j<n;++j)if(bit(v,j)){gr[c].push_back(j);vg[j].push_back(c);}}
        for(auto &a:H)for(auto &b:G){int p=0;for(int w=0;w<words;++w)p^=__builtin_parityll(a[w]&b[w]);if(p)throw std::runtime_error("checks/stabilizers do not commute");}
    }
    bool kernel(const Word &v)const{for(auto &h:H){int p=0;for(int j=0;j<words;++j)p^=__builtin_parityll(h[j]&v[j]);if(p)return false;}return true;}
    void validate_translation(int P)const{
        if(P<1||n%P)throw std::runtime_error("invalid cyclic period");
        Basis h(n);for(auto &v:H)h.add(v);
        auto shifted=[&](const Word &v){Word u{};for(int j=0;j<n;++j)if(bit(v,j))toggle(u,(j/P)*P+(j%P+1)%P);return u;};
        for(auto &v:H)if(h.nontrivial(shifted(v)))throw std::runtime_error("check space is not translation invariant");
        for(auto &v:G)if(basis.nontrivial(shifted(v)))throw std::runtime_error("stabilizer space is not translation invariant");
    }
};
inline void print_support(const Word &v,int n){std::cout<<'[';bool first=true;for(int j=0;j<n;++j)if(bit(v,j)){if(!first)std::cout<<',';std::cout<<j;first=false;}std::cout<<']';}
