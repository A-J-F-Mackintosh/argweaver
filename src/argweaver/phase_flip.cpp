// c++ includes
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

// arghmm includes
#include "common.h"
#include "logging.h"
#include "matrices.h"
#include "phase_flip.h"
#include "recomb.h"
#include "sample_thread.h"
#include "states.h"
#include "thread.h"
#include "total_prob.h"

namespace argweaver {

using namespace std;

bool flip_force_accept = false;


// A leaf thread over positions [start, start + path.size()): the index of
// its state among the coalescent states of the partial ARG's local tree at
// each position, and its recombinations (positions after 'start' only).
struct ThreadPath
{
    int start;
    vector<int> path;
    vector<int> recomb_pos;
    vector<Spr> recombs;

    int at(int pos) const { return path[pos - start]; }
    const int *abs_path() const { return &path[0] - start; }
};


static int find_leaf(const LocalTrees *trees, int seqid)
{
    for (int i=0; i<trees->get_num_leaves(); i++)
        if (trees->seqids[i] == seqid)
            return i;
    return -1;
}


static bool same_spr(const Spr &a, const Spr &b)
{
    return a.recomb_node == b.recomb_node && a.recomb_time == b.recomb_time &&
        a.coal_node == b.coal_node && a.coal_time == b.coal_time;
}


// The thread of leaf 'chrom' in the ARG 'full', relative to 'partial' (which
// over [lo, hi) must be 'full' without 'chrom', up to the naming of nodes).
// That is, the thread path and recombinations which, added to 'partial' by
// add_arg_thread, give 'full'.  Returns false if the thread cannot be
// written in that form.
static bool extract_thread(const ArgModel *model, const LocalTrees *full,
                           const LocalTrees *partial, int chrom, int lo,
                           int hi, ThreadPath *out)
{
    out->start = lo;
    out->path.assign(hi - lo, -1);
    out->recomb_pos.clear();
    out->recombs.clear();

    const int leaf = find_leaf(full, chrom);
    assert(leaf != -1);
    StatesModel states_model(model->ntimes, false, 0, model->pop_tree, 0);
    States states;
    vector<Spr> candidates;
    vector<int> recon;

    LocalTrees::const_iterator pit = partial->begin();
    int pstart = partial->start_coord;
    int pend = pstart + pit->blocklen;
    int fend = full->start_coord;
    for (LocalTrees::const_iterator fit=full->begin(); fit!=full->end();
         ++fit) {
        const int fstart = fend;
        fend += fit->blocklen;
        if (fit->blocklen == 0 || fend <= lo)
            continue;
        if (fstart >= hi)
            break;
        while (max(fstart, lo) >= pend) {
            ++pit;
            assert(pit != partial->end());
            pstart = pend;
            pend += pit->blocklen;
        }
        const LocalTree *ftree = fit->tree;
        const LocalTree *ptree = pit->tree;

        // where the leaf's branch joins the rest of the tree
        recon.resize(ftree->nnodes);
        map_congruent_trees(ftree, &full->seqids[0], ptree,
                            &partial->seqids[0], &recon[0]);
        int ptr = leaf;
        while (ptr != -1 && recon[ptr] == -1)
            ptr = ftree->nodes[ptr].parent;
        if (ptr == -1) {
            return false;
        }
        State state(recon[ptr], ftree->nodes[ptr].age);
        states_model.get_coal_states(ptree, states);
        int idx = -1;
        for (unsigned int j=0; j<states.size(); j++)
            if (states[j].node == state.node && states[j].time == state.time)
                idx = j;
        if (idx == -1) {
            return false;
        }
        state = states[idx];

        // a recombination of the thread starts this block if the partial
        // ARG has no break here
        if (fstart > lo && fstart != pstart && !fit->spr.is_null()) {
            const State last_state = states[out->at(fstart - 1)];
            Spr spr;
            if (fit->spr.recomb_node == leaf)
                spr = Spr(-1, fit->spr.recomb_time, state.node, state.time,
                          state.pop_path);
            else
                spr = Spr(state.node, fit->spr.recomb_time, -1, state.time,
                          ptree->nodes[state.node].pop_path);
            candidates.clear();
            get_possible_recomb(model, ptree, last_state, state, false,
                                candidates);
            bool found = false;
            for (unsigned int j=0; j<candidates.size(); j++)
                if (same_spr(candidates[j], spr))
                    found = true;
            if (!found) {
                return false;
            }
            out->recomb_pos.push_back(fstart);
            out->recombs.push_back(spr);
        }

        for (int i=max(fstart, lo); i<min(fend, hi); i++)
            out->path[i - lo] = idx;
    }
    for (unsigned int i=0; i<out->path.size(); i++)
        if (out->path[i] == -1)
            return false;
    return true;
}


// A copy of the local trees of 'trees' over [lo, hi), or NULL if a block
// inside it is empty
static LocalTrees *copy_segment(const LocalTrees *trees, int lo, int hi)
{
    LocalTrees *seg = new LocalTrees(lo, hi, trees->nnodes);
    seg->chrom = trees->chrom;
    seg->seqids = trees->seqids;
    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end();
         ++it) {
        const int start = end;
        end += it->blocklen;
        if (end <= lo)
            continue;
        if (start >= hi)
            break;
        if (it->blocklen == 0) {
            delete seg;
            return NULL;
        }
        LocalTree *tree = new LocalTree();
        tree->copy(*it->tree);
        Spr spr;
        spr.set_null();
        int *mapping = NULL;
        if (seg->get_num_trees() > 0) {
            spr = it->spr;
            if (it->mapping) {
                mapping = new int [trees->nnodes];
                for (int i=0; i<trees->nnodes; i++)
                    mapping[i] = it->mapping[i];
            }
        }
        seg->trees.push_back(LocalTreeSpr(tree, spr,
                                          min(end, hi) - max(start, lo),
                                          mapping));
    }
    return seg;
}


// Rename the leaves of every local tree: leaf i becomes leaf perm[i], and
// takes its sequence id with it.  Internal nodes keep their names.
static void permute_leaves(LocalTrees *trees, const vector<int> &perm)
{
    const int nnodes = trees->nnodes;
    const int nleaves = perm.size();
    vector<int> p(nnodes);
    for (int i=0; i<nnodes; i++)
        p[i] = (i < nleaves) ? perm[i] : i;
    vector<LocalNode> nodes(nnodes);
    vector<int> mapping(nnodes);
    for (LocalTrees::iterator it=trees->begin(); it!=trees->end(); ++it) {
        LocalTree *tree = it->tree;
        for (int i=0; i<tree->nnodes; i++)
            nodes[p[i]] = tree->nodes[i];
        for (int i=0; i<tree->nnodes; i++) {
            LocalNode &node = nodes[i];
            if (node.parent != -1) node.parent = p[node.parent];
            if (node.child[0] != -1) node.child[0] = p[node.child[0]];
            if (node.child[1] != -1) node.child[1] = p[node.child[1]];
            tree->nodes[i] = node;
        }
        if (tree->root != -1)
            tree->root = p[tree->root];
        if (it->mapping) {
            for (int i=0; i<tree->nnodes; i++)
                mapping[p[i]] = (it->mapping[i] == -1) ? -1 :
                    p[it->mapping[i]];
            for (int i=0; i<tree->nnodes; i++)
                it->mapping[i] = mapping[i];
        }
        if (!it->spr.is_null()) {
            it->spr.recomb_node = p[it->spr.recomb_node];
            it->spr.coal_node = p[it->spr.coal_node];
        }
    }
    vector<int> seqids(trees->seqids);
    for (int i=0; i<nleaves; i++)
        trees->seqids[perm[i]] = seqids[i];
}


// Replace the local trees of 'trees' over [lo, hi] with 'window', and swap
// the leaves of haplotypes u and v (structure, not sequence ids) from hi+1
// on.  'window' must match 'trees' at lo, and at hi up to the swap; it is
// left empty.
static void splice_window(const ArgModel *model, LocalTrees *trees,
                          LocalTrees *window, int lo, int hi, int u, int v)
{
    // put the window's leaves in the same order as the ARG's
    vector<int> perm(window->get_num_leaves());
    for (unsigned int i=0; i<perm.size(); i++)
        perm[i] = find_leaf(trees, window->seqids[i]);
    permute_leaves(window, perm);

    LocalTrees *mid = partition_local_trees(trees, lo, true);
    LocalTrees *right = partition_local_trees(mid, hi + 1, true);

    // a recombination exactly at hi+1 is left as an empty block at the end
    // of 'mid'; it moves onto the first tree of 'right'
    Spr stub_spr;
    stub_spr.set_null();
    if (mid->trees.back().blocklen == 0) {
        stub_spr = mid->trees.back().spr;
        mid->trees.back().clear();
        mid->trees.pop_back();
    }
    const LocalTree *old_last = mid->trees.back().tree;

    // swap u and v from hi+1 on
    const int leaf_u = find_leaf(trees, u), leaf_v = find_leaf(trees, v);
    vector<int> swap(trees->get_num_leaves());
    for (unsigned int i=0; i<swap.size(); i++)
        swap[i] = i;
    swap[leaf_u] = leaf_v;
    swap[leaf_v] = leaf_u;
    if (right->get_num_trees() > 0) {
        permute_leaves(right, swap);
        std::swap(right->seqids[leaf_u], right->seqids[leaf_v]);
    }

    if (!stub_spr.is_null()) {
        // the recombination node, renamed: through the swap, then onto the
        // window's last tree (which matches the old one after the swap)
        LocalTrees old_tree(hi, hi + 1, trees->nnodes);
        old_tree.seqids = trees->seqids;
        LocalTree *tree = new LocalTree();
        tree->copy(*old_last);
        old_tree.trees.push_back(LocalTreeSpr(tree, Spr(-1, -1, -1, -1, -1),
                                              1, NULL));
        permute_leaves(&old_tree, swap);
        std::swap(old_tree.seqids[leaf_u], old_tree.seqids[leaf_v]);
        vector<int> recon(trees->nnodes);
        map_congruent_trees(old_tree.front().tree, &old_tree.seqids[0],
                            window->back().tree, &window->seqids[0],
                            &recon[0]);
        Spr spr(stub_spr);
        spr.recomb_node = recon[swap.size() > (unsigned) spr.recomb_node ?
                                swap[spr.recomb_node] : spr.recomb_node];
        spr.coal_node = recon[swap.size() > (unsigned) spr.coal_node ?
                              swap[spr.coal_node] : spr.coal_node];
        right->front().spr = spr;
        right->front().mapping = new int [trees->nnodes];
    }
    delete mid;

    append_local_trees(trees, window, true, model->pop_tree);
    if (right->get_num_trees() > 0)
        append_local_trees(trees, right, true, model->pop_tree);
    delete right;
}


// Forward table for threading 'chrom' (not in 'seg') through 'seg', with its
// state at the first column given
static void pinned_forward(const ArgModel *model, const Sequences *sequences,
                           const LocalTrees *seg, int chrom, int first_state,
                           ArgHmmForwardTable *forward)
{
    ArgHmmMatrixIter matrix_iter(model, sequences, seg, chrom);
    matrix_iter.set_start_pop(sequences->get_pop(chrom));
    States states;
    matrix_iter.begin();
    matrix_iter.get_coal_states(states);
    forward->new_block(matrix_iter.get_block_start(),
                       matrix_iter.get_block_end(), states.size());
    double *col = forward->get_table()[seg->start_coord];
    fill(col, col + states.size(), 0.0);
    col[first_state] = 1.0;
    arghmm_forward_alg(seg, model, sequences, &matrix_iter, forward, NULL,
                       true, false);
}


// Log probability that stochastic traceback (with the last state given)
// draws 'path', mirroring stochastic_traceback
static double traceback_logp(const LocalTrees *trees,
                             ArgHmmMatrixIter *matrix_iter, double **fw,
                             const int *path)
{
    States states;
    vector<double> trans;
    double lnl = 0.0;

    matrix_iter->rbegin();
    int pos = trees->end_coord;
    for (; matrix_iter->more(); matrix_iter->prev()) {
        ArgHmmMatrices &mat = matrix_iter->ref_matrices();
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        mat.states_model.get_coal_states(tree, states);
        const int nstates = max((int) states.size(), 1);
        trans.resize(nstates);
        pos -= mat.blocklen;

        // within the block, then into its first column from the last
        // column of the previous block.  Transition probabilities into
        // state k are recomputed only when k changes.
        const int first = (pos > trees->start_coord) ? pos - 1 : pos;
        int last_k = -1;
        for (int i=pos + mat.blocklen - 2; i>=first; i--) {
            const int k = path[i+1];
            double total = 0.0;
            if (i == pos - 1 && mat.transmat_switch) {
                const int n1 = max(mat.transmat_switch->nstates1, 1);
                for (int j=0; j<n1; j++)
                    total += fw[i][j] * mat.transmat_switch->get(j, k);
                lnl += log(fw[i][path[i]] *
                           mat.transmat_switch->get(path[i], k) / total);
                continue;
            }
            if (k != last_k) {
                for (int j=0; j<nstates; j++)
                    trans[j] = mat.transmat->get(tree, states, j, k);
                last_k = k;
            }
            for (int j=0; j<nstates; j++)
                total += fw[i][j] * trans[j];
            lnl += log(fw[i][path[i]] * trans[path[i]] / total);
        }
    }
    return lnl;
}


// Log probability that sample_recombinations, given thread path 'path',
// chooses the recombinations (recomb_pos, recombs).  Returns -INFINITY if
// it cannot.
static double recombinations_logp(const LocalTrees *trees,
                                  const ArgModel *model,
                                  ArgHmmMatrixIter *matrix_iter,
                                  const int *path,
                                  const vector<int> &recomb_pos,
                                  const vector<Spr> &recombs)
{
    States states;
    LineageCounts lineages(model->ntimes, model->num_pops());
    vector<Spr> candidates;
    double lnl = 0.0;
    unsigned int irecomb = 0;

    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        ArgHmmMatrices &matrices = matrix_iter->ref_matrices();
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        lineages.count(tree, model->pop_tree, false);
        matrices.states_model.get_coal_states(tree, states);

        int start = matrix_iter->get_block_start();
        const int end = matrix_iter->get_block_end();
        if (matrices.transmat_switch || start == trees->start_coord)
            start++;
        // no recombination is allowed where the block starts
        if (irecomb < recomb_pos.size() && recomb_pos[irecomb] < start)
            return -INFINITY;

        // as in sample_recombinations: under the SMC, the position of
        // the next recombination where the state stays the same is drawn
        // as the integer part of an exponential waiting time
        int draw_start = -1;
        double draw_rate = 0.0;
        for (int i=start; i<end; i++) {
            const bool is_recomb = (irecomb < recomb_pos.size() &&
                                    recomb_pos[irecomb] == i);
            if (path[i] == path[i-1]) {
                if (model->smc_prime) {
                    if (is_recomb) return -INFINITY;
                    continue;
                }
                if (draw_start < 0) {
                    const int last_state = path[i-1];
                    TransMatrix *m = matrices.transmat;
                    const int a = states[last_state].time;
                    const double self_trans = m->get(
                        tree, states, last_state, last_state);
                    draw_rate = 1.0 - (m->norecombs[a] / self_trans);
                    draw_start = i;
                }
                if (!is_recomb)
                    continue;
                lnl += -draw_rate * (i - draw_start) +
                    log(1.0 - exp(-draw_rate));
                draw_start = -1;
            } else {
                if (!is_recomb)
                    return -INFINITY;
                if (draw_start >= 0) {
                    lnl += -draw_rate * (i - draw_start);
                    draw_start = -1;
                }
            }

            // which recombination
            const State state = states[path[i]];
            const State last_state = states[path[i-1]];
            candidates.clear();
            get_possible_recomb(model, tree, last_state, state, false,
                                candidates);
            double total = 0.0, chosen = 0.0;
            for (unsigned int j=0; j<candidates.size(); j++) {
                double p = recomb_prob_unnormalized(
                    model, tree, lineages, last_state, state, candidates[j],
                    false);
                total += p;
                if (same_spr(candidates[j], recombs[irecomb]))
                    chosen = p;
            }
            if (chosen == 0.0)
                return -INFINITY;
            lnl += log(chosen / total);
            irecomb++;
        }
        if (draw_start >= 0)
            lnl += -draw_rate * (end - draw_start);
    }
    if (irecomb != recomb_pos.size())
        return -INFINITY;
    return lnl;
}


// Log probability of threading 'chrom' through 'seg' along 'path' with the
// given recombinations, with its states at the first and last columns given
static double pinned_logp(const ArgModel *model, const Sequences *sequences,
                          const LocalTrees *seg, int chrom, const int *path,
                          const vector<int> &recomb_pos,
                          const vector<Spr> &recombs)
{
    ArgHmmForwardTable forward(seg->start_coord, seg->length());
    pinned_forward(model, sequences, seg, chrom, path[seg->start_coord],
                   &forward);
    ArgHmmMatrixIter matrix_iter(model, NULL, seg, chrom);
    matrix_iter.set_start_pop(sequences->get_pop(chrom));
    double lnl = traceback_logp(seg, &matrix_iter, forward.get_table(),
                                path);
    return lnl + recombinations_logp(seg, model, &matrix_iter, path,
                                     recomb_pos, recombs);
}


// Sample a thread for 'chrom' through 'seg', with its states at the first
// and last columns given.  'path' (indexed by position) must hold them.
// Returns the log probability of the thread.
static double sample_pinned(const ArgModel *model, const Sequences *sequences,
                            const LocalTrees *seg, int chrom, int *path,
                            vector<int> &recomb_pos, vector<Spr> &recombs)
{
    ArgHmmForwardTable forward(seg->start_coord, seg->length());
    pinned_forward(model, sequences, seg, chrom, path[seg->start_coord],
                   &forward);
    ArgHmmMatrixIter matrix_iter(model, NULL, seg, chrom);
    matrix_iter.set_start_pop(sequences->get_pop(chrom));
    stochastic_traceback(seg, model, &matrix_iter, forward.get_table(), path,
                         true, false);
    double lnl = traceback_logp(seg, &matrix_iter, forward.get_table(),
                                path);
    recomb_pos.clear();
    recombs.clear();
    sample_recombinations(seg, model, &matrix_iter, path, recomb_pos,
                          recombs, false);
    return lnl + recombinations_logp(seg, model, &matrix_iter, path,
                                     recomb_pos, recombs);
}


// Add the thread of 'chrom' to 'trees': 'left' up to column lo, 'mid'
// (indexed by position) from lo to hi and 'right' from hi on
static void add_composite_thread(const ArgModel *model,
                                 const Sequences *sequences,
                                 LocalTrees *trees, int chrom,
                                 const ThreadPath &left, const int *mid,
                                 const vector<int> &mid_pos,
                                 const vector<Spr> &mid_recombs,
                                 const ThreadPath &right, int lo, int hi)
{
    const int start = trees->start_coord;
    vector<int> path(trees->length());
    for (int i=start; i<trees->end_coord; i++)
        path[i - start] = (i <= lo) ? left.at(i) :
            (i < hi ? mid[i] : right.at(i));
    assert(mid[lo] == left.at(lo) && mid[hi] == right.at(hi));

    vector<int> recomb_pos;
    vector<Spr> recombs;
    for (unsigned int j=0; j<left.recomb_pos.size(); j++)
        if (left.recomb_pos[j] <= lo) {
            recomb_pos.push_back(left.recomb_pos[j]);
            recombs.push_back(left.recombs[j]);
        }
    recomb_pos.insert(recomb_pos.end(), mid_pos.begin(), mid_pos.end());
    recombs.insert(recombs.end(), mid_recombs.begin(), mid_recombs.end());
    for (unsigned int j=0; j<right.recomb_pos.size(); j++)
        if (right.recomb_pos[j] > hi) {
            recomb_pos.push_back(right.recomb_pos[j]);
            recombs.push_back(right.recombs[j]);
        }

    StatesModel states_model(model->ntimes, false, 0, model->pop_tree,
                             sequences->get_pop(chrom));
    add_arg_thread(trees, states_model, model->ntimes, &path[0] - start,
                   chrom, recomb_pos, recombs, model->pop_tree);
}


// recombinations of a thread in (lo, hi]
static void recombs_in(const ThreadPath &thread, int lo, int hi,
                       vector<int> &recomb_pos, vector<Spr> &recombs)
{
    recomb_pos.clear();
    recombs.clear();
    for (unsigned int j=0; j<thread.recomb_pos.size(); j++)
        if (thread.recomb_pos[j] > lo && thread.recomb_pos[j] <= hi) {
            recomb_pos.push_back(thread.recomb_pos[j]);
            recombs.push_back(thread.recombs[j]);
        }
}


void flip_window(const LocalTrees *trees, const PhaseBlocks &blocks, int k,
                 int margin, int *lo, int *hi)
{
    assert(k > 0 && k < blocks.size());
    *lo = max(trees->start_coord, blocks.sites[k-1] - margin);
    *hi = min(trees->end_coord - 1, blocks.sites[k] + margin);
}


bool flip_phase_boundary(const ArgModel *model, Sequences *sequences,
                         LocalTrees *trees, const PhaseBlocks &blocks, int k,
                         int margin, FlipMoveInfo *info, int first)
{
    assert(!blocks.same_block(k));
    FlipMoveInfo local;
    if (!info)
        info = &local;
    int lo, hi;
    flip_window(trees, blocks, k, margin, &lo, &hi);
    const int flip_pos = blocks.sites[k];
    const bool hap1_first = (first == -1) ? (irand(2) == 0) :
        (first == blocks.hap1);
    const int u = hap1_first ? blocks.hap1 : blocks.hap2;
    const int v = hap1_first ? blocks.hap2 : blocks.hap1;
    info->lo = lo;
    info->hi = hi;
    info->first = u;
    info->dprior = info->dlik = info->lq_fwd = info->lq_rev = 0.0;
    info->log_ratio = -INFINITY;
    info->supported = false;

    // Everything is done on the local trees of the window, [lo, hi]: the
    // ARG, the ARG without v, without u, and without both
    LocalTrees *window = copy_segment(trees, lo, hi + 1);
    if (!window)
        return false;
    LocalTrees no_v, no_u, core;
    no_v.copy(*window);
    remove_arg_thread(&no_v, v, model);
    no_u.copy(*window);
    remove_arg_thread(&no_u, u, model);
    core.copy(no_v);
    remove_arg_thread(&core, u, model);

    // old threads: u relative to the core, v relative to the core and u.
    // After the flip, u ends on v's old state relative to the core.
    ThreadPath thread_u, thread_v, thread_v_core;
    if (!extract_thread(model, &no_v, &core, u, lo, hi + 1, &thread_u) ||
        !extract_thread(model, window, &no_v, v, lo, hi + 1, &thread_v) ||
        !extract_thread(model, &no_u, &core, v, hi, hi + 1,
                        &thread_v_core)) {
        delete window;
        return false;
    }

    // probability of the reverse move proposing the old threads
    vector<int> rpos;
    vector<Spr> recs;
    recombs_in(thread_u, lo, hi, rpos, recs);
    double lq_rev = pinned_logp(model, sequences, &core, u,
                                thread_u.abs_path(), rpos, recs);
    recombs_in(thread_v, lo, hi, rpos, recs);
    lq_rev += pinned_logp(model, sequences, &no_v, v, thread_v.abs_path(),
                          rpos, recs);
    if (!isfinite(lq_rev)) {
        delete window;
        return false;
    }
    info->supported = true;
    const double prior_old = calc_arg_prior(model, window);
    const double lik_old = calc_arg_likelihood(model, sequences, window);

    // swap the alleles from the boundary on
    vector<int> swapped;
    for (int i=flip_pos; i<sequences->length(); i++)
        if (sequences->seqs[u][i] != sequences->seqs[v][i]) {
            sequences->switch_alleles(i, u, v);
            swapped.push_back(i);
        }

    // propose u, from its old state at lo to v's old state at hi
    // (removing v then u, as for 'core', so that nodes are named the same)
    vector<int> mid_alloc(hi + 1 - lo);
    int *mid = &mid_alloc[0] - lo;
    mid[lo] = thread_u.at(lo);
    mid[hi] = thread_v_core.at(hi);
    double lq_fwd = sample_pinned(model, sequences, &core, u, mid, rpos,
                                  recs);
    LocalTrees proposal;
    proposal.copy(*window);
    remove_arg_thread(&proposal, v, model);
    remove_arg_thread(&proposal, u, model);
    add_composite_thread(model, sequences, &proposal, u, thread_u, mid,
                         rpos, recs, thread_v_core, lo, hi);

    // propose v, from its old state at lo to u's old state at hi
    ThreadPath left_v, right_v;
    bool ok = extract_thread(model, window, &proposal, v, lo, lo + 1,
                             &left_v);
    const int leaf_u = find_leaf(window, u), leaf_v = find_leaf(window, v);
    window->seqids[leaf_u] = v;
    window->seqids[leaf_v] = u;
    ok = ok && extract_thread(model, window, &proposal, v, hi, hi + 1,
                              &right_v);
    window->seqids[leaf_u] = u;
    window->seqids[leaf_v] = v;
    if (ok) {
        mid[lo] = left_v.at(lo);
        mid[hi] = right_v.at(hi);
        lq_fwd += sample_pinned(model, sequences, &proposal, v, mid, rpos,
                                recs);
        add_composite_thread(model, sequences, &proposal, v, left_v, mid,
                             rpos, recs, right_v, lo, hi);
    } else {
        printError("flip_phase_boundary: cannot rebuild the thread of %d\n",
                   v);
    }

    if (ok) {
        info->dprior = calc_arg_prior(model, &proposal) - prior_old;
        info->dlik = calc_arg_likelihood(model, sequences, &proposal) -
            lik_old;
        info->lq_fwd = lq_fwd;
        info->lq_rev = lq_rev;
        info->log_ratio = info->dprior + info->dlik + lq_rev - lq_fwd;
        if (flip_force_accept || info->log_ratio >= 0.0 ||
            frand() < exp(info->log_ratio)) {
            splice_window(model, trees, &proposal, lo, hi, u, v);
            delete window;
            return true;
        }
    }

    // rejected: restore the alleles
    for (unsigned int j=0; j<swapped.size(); j++)
        sequences->switch_alleles(swapped[j], u, v);
    delete window;
    return false;
}


bool rebuild_leaf_thread(const ArgModel *model, const Sequences *sequences,
                         const LocalTrees *full, int seqid,
                         LocalTrees *rebuilt)
{
    LocalTrees partial;
    partial.copy(*full);
    remove_arg_thread(&partial, seqid, model);
    ThreadPath thread;
    if (!extract_thread(model, full, &partial, seqid, full->start_coord,
                        full->end_coord, &thread))
        return false;
    rebuilt->copy(*full);
    remove_arg_thread(rebuilt, seqid, model);
    StatesModel states_model(model->ntimes, false, 0, model->pop_tree,
                             sequences->get_pop(seqid));
    add_arg_thread(rebuilt, states_model, model->ntimes,
                   &thread.path[0] - thread.start, seqid, thread.recomb_pos,
                   thread.recombs, model->pop_tree);
    return true;
}


int flip_phase_boundaries(const ArgModel *model, Sequences *sequences,
                          LocalTrees *trees, const PhaseBlocks &blocks,
                          int margin, int *ntried, double try_prob)
{
    int naccept = 0, ntry = 0;
    for (int k=1; k<blocks.size(); k++) {
        if (blocks.same_block(k))
            continue;
        if (try_prob < 1.0 && frand() >= try_prob)
            continue;
        ntry++;
        if (flip_phase_boundary(model, sequences, trees, blocks, k, margin))
            naccept++;
    }
    if (ntried)
        *ntried = ntry;
    return naccept;
}

} // namespace argweaver
