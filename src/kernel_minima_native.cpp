// Exhaustive exclusion of nontrivial Z-logicals through weight W for a CSS code with a cyclic block structure.
// Implements the irreducible-support search of the draft's Appendix A:
//   * branch on an unsatisfied check with the fewest available neighbours; k-th branch adds neighbour k and forbids 1..k-1
//   * stop at zero syndrome (test nontriviality against the X-logical basis)
//   * pruning: syndrome weight <= budget * max column degree; disjoint-neighbourhood packing; stabilizer half-overlap rule
//   * translation normalisation: roots (block b, position 0) with all earlier blocks forbidden
// Usage: excl <codefile> <W> <root_lo> <root_hi> [time_limit_sec]
// Exit text: "COMPLETE" (no nontrivial logical of weight <= W with a root in [root_lo, root_hi)) or "FOUND w ..." or "TIMEOUT".
#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
using namespace std;
typedef uint64_t u64;

static int n, P, nblocks, rX, rZ, kX;
static vector<vector<int>> chk_cols, col_chks, stab_cols, col_stabs, xcols; // xcols[j] = indices of X-basis vectors containing j
static vector<int> stab_w;
static int Delta;                       // max column degree in H_X
static const int NW = 16;               // 1024 bits
static int NWcols;                      // words for n bits
static vector<u64> blocked;             // T U forbidden
static vector<int> T;
static vector<u64> synd;                // m bits
static int syndw;
static vector<int> ov;                  // stabilizer overlap counters
static vector<u64> lsyn;                // logical syndrome bits (kX)
static int lsynw;
static long long nodes = 0;
static int Wmax;
static bool found = false;
static vector<int> found_word;
static double tlimit = 1e18;
static chrono::steady_clock::time_point t_start;
static bool timed_out = false;
static int child_lo = 0, child_hi = 1 << 30, depth = 0;
static int modD = -1, modM = 1, modR = 0; static long long modcount = 0;

static inline bool getbit(const vector<u64>& v, int i) { return (v[i >> 6] >> (i & 63)) & 1ULL; }
static inline void flipbit(vector<u64>& v, int i) { v[i >> 6] ^= (1ULL << (i & 63)); }
static inline void setbit(vector<u64>& v, int i) { v[i >> 6] |= (1ULL << (i & 63)); }
static inline void clrbit(vector<u64>& v, int i) { v[i >> 6] &= ~(1ULL << (i & 63)); }
static inline int popc(const vector<u64>& v) { int c = 0; for (u64 x : v) c += __builtin_popcountll(x); return c; }
static inline bool anybit(const vector<u64>& v) { for (u64 x : v) if (x) return true; return false; }

static void add_coord(int j) {
    T.push_back(j);
    setbit(blocked, j);
    for (int c : col_chks[j]) flipbit(synd, c);
    for (int g : col_stabs[j]) ov[g]++;
    for (int x : xcols[j]) flipbit(lsyn, x);
}
static void remove_coord(int j) {
    T.pop_back();
    // blocked bit handled by caller (forbid semantics)
    for (int c : col_chks[j]) flipbit(synd, c);
    for (int g : col_stabs[j]) ov[g]--;
    for (int x : xcols[j]) flipbit(lsyn, x);
}
static bool overlap_ok(int j) {
    for (int g : col_stabs[j]) if (2 * (ov[g] + 1) > stab_w[g]) return false;
    return true;
}

static void dfs(int budget) {
    if (found || timed_out) return;
    if (depth == modD) { long long c = modcount++; if (c % modM != modR) return; }
    nodes++;
    if ((nodes & 0xFFFFF) == 0) {
        double el = chrono::duration<double>(chrono::steady_clock::now() - t_start).count();
        if (el > tlimit) { timed_out = true; return; }
    }
    int sw = popc(synd);
    if (sw == 0) {
        if (anybit(lsyn)) { printf("WORD %d", (int)T.size()); for (int q : T) printf(" %d", q); printf("\n"); }
        return;
    }
    if (budget == 0) return;
    if (sw > budget * Delta) return;
    // collect unsatisfied checks, count available neighbours
    static thread_local vector<int> unsat; unsat.clear();
    for (int w = 0; w < syndw; w++) {
        u64 x = synd[w];
        while (x) { int b = __builtin_ctzll(x); x &= x - 1; unsat.push_back(w * 64 + b); }
    }
    int best_c = -1, best_cnt = INT_MAX;
    static thread_local vector<pair<int,int>> cnts; cnts.clear();
    for (int c : unsat) {
        int cnt = 0;
        for (int j : chk_cols[c]) if (!getbit(blocked, j)) cnt++;
        if (cnt == 0) return;
        cnts.push_back({cnt, c});
        if (cnt < best_cnt) { best_cnt = cnt; best_c = c; }
    }
    // packing bound: greedy disjoint available neighbourhoods, checks with few neighbours first
    if ((int)cnts.size() > budget) {
        sort(cnts.begin(), cnts.end());
        static thread_local vector<u64> covered; covered.assign(NWcols, 0);
        int p = 0;
        for (auto& pc : cnts) {
            int c = pc.second;
            bool disjoint = true;
            for (int j : chk_cols[c]) if (!getbit(blocked, j) && getbit(covered, j)) { disjoint = false; break; }
            if (disjoint) {
                p++;
                if (p > budget) return;
                for (int j : chk_cols[c]) if (!getbit(blocked, j)) setbit(covered, j);
            }
        }
    }
    // branch on best_c's available neighbours, ordered by number of unsatisfied checks repaired (desc)
    static thread_local vector<pair<int,int>> nb;
    vector<pair<int,int>> local;
    for (int j : chk_cols[best_c]) if (!getbit(blocked, j)) {
        int rep = 0;
        for (int c : col_chks[j]) if (getbit(synd, c)) rep++;
        local.push_back({-rep, j});
    }
    sort(local.begin(), local.end());
    vector<int> forb;
    int childidx = 0;
    for (auto& pr : local) {
        int j = pr.second;
        bool explore = true;
        if (depth == 0) { explore = (childidx >= child_lo && childidx < child_hi); childidx++; }
        if (overlap_ok(j) && explore) {
            add_coord(j);
            depth++;
            dfs(budget - 1);
            depth--;
            remove_coord(j);
            // keep blocked (forbidden) for later branches
        } else {
            setbit(blocked, j);
        }
        forb.push_back(j);
        if (depth == 0 && childidx >= child_hi) { for (int jj : forb) clrbit(blocked, jj); return; }
        if (found || timed_out) break;
    }
    for (int j : forb) clrbit(blocked, j);
}

int main(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "usage: excl codefile W root_lo root_hi [tlimit]\n"); return 1; }
    ifstream in(argv[1]);
    Wmax = atoi(argv[2]);
    int root_lo = atoi(argv[3]), root_hi = atoi(argv[4]);
    if (argc > 5) tlimit = atof(argv[5]);
    if (argc > 7) { child_lo = atoi(argv[6]); child_hi = atoi(argv[7]); }
    if (argc > 10) { modD = atoi(argv[8]); modM = atoi(argv[9]); modR = atoi(argv[10]); }
    in >> n >> P >> nblocks;
    in >> rX; chk_cols.resize(rX); col_chks.assign(n, {});
    string line; getline(in, line);
    for (int r = 0; r < rX; r++) { getline(in, line); istringstream ss(line); int j; while (ss >> j) { chk_cols[r].push_back(j); col_chks[j].push_back(r); } }
    in >> rZ; stab_cols.resize(rZ); col_stabs.assign(n, {}); stab_w.assign(rZ, 0); getline(in, line);
    for (int r = 0; r < rZ; r++) { getline(in, line); istringstream ss(line); int j; while (ss >> j) { stab_cols[r].push_back(j); col_stabs[j].push_back(r); stab_w[r]++; } }
    in >> kX; xcols.assign(n, {}); getline(in, line);
    for (int r = 0; r < kX; r++) { getline(in, line); istringstream ss(line); int j; while (ss >> j) xcols[j].push_back(r); }
    Delta = 0; for (int j = 0; j < n; j++) Delta = max(Delta, (int)col_chks[j].size());
    NWcols = (n + 63) / 64; syndw = (rX + 63) / 64; lsynw = (kX + 63) / 64;
    blocked.assign(NWcols, 0); synd.assign(syndw, 0); ov.assign(rZ, 0); lsyn.assign(lsynw, 0);
    t_start = chrono::steady_clock::now();
    fprintf(stderr, "n=%d P=%d blocks=%d rX=%d rZ=%d kX=%d Delta=%d W=%d roots [%d,%d)\n", n, P, nblocks, rX, rZ, kX, Delta, Wmax, root_lo, root_hi);
    for (int b = root_lo; b < root_hi && !found && !timed_out; b++) {
        // forbid earlier blocks
        fill(blocked.begin(), blocked.end(), 0ULL);
        for (int j = 0; j < b * P; j++) setbit(blocked, j);
        T.clear(); fill(synd.begin(), synd.end(), 0ULL); fill(ov.begin(), ov.end(), 0); fill(lsyn.begin(), lsyn.end(), 0ULL);
        int root = b * P;
        if (!overlap_ok(root)) continue;
        add_coord(root);
        long long n0 = nodes;
        dfs(Wmax - 1);
        double el = chrono::duration<double>(chrono::steady_clock::now() - t_start).count();
        fprintf(stderr, "root block %d done: nodes %lld (total %lld) t=%.1fs\n", b, nodes - n0, nodes, el);
    }
    double el = chrono::duration<double>(chrono::steady_clock::now() - t_start).count();
    if (found) {
        printf("FOUND %d", (int)found_word.size());
        for (int j : found_word) printf(" %d", j);
        printf("\n");
    } else if (timed_out) {
        printf("TIMEOUT nodes=%lld t=%.1f\n", nodes, el);
    } else {
        printf("ENUMERATION_COMPLETE W=%d roots=[%d,%d) children=[%d,%d) mod(D=%d,M=%d,r=%d) nodes=%lld t=%.1f\n", Wmax, root_lo, root_hi, child_lo, child_hi, modD, modM, modR, nodes, el);
    }
    return 0;
}
