// c++ includes
#include <assert.h>
#include <map>
#include <vector>

// arghmm includes
#include "common.h"
#include "logging.h"
#include "matrices.h"
#include "phase_thread.h"
#include "recomb.h"
#include "total_prob.h"
#include "sample_thread.h"
#include "thread.h"
#include "trans.h"

namespace argweaver {

using namespace std;


// probability of moving from phase state a to phase state b at a column
static inline double phase_trans(char mix, int a, int b, double switch_prob)
{
    if (mix == PHASE_MIX_SWITCH)
        return (a == b) ? 1.0 - switch_prob : switch_prob;
    if (mix == PHASE_MIX_REDRAW)
        return 0.5;
    return (a == b) ? 1.0 : 0.0;
}


// Current and swapped orientation emissions of one block when both of the
// individual's haplotypes are in the ARG.  'emit' holds ARGweaver's unphased
// emissions (the average over both orientations at heterozygous sites) and
// 'phase_pr' the fraction of each due to the current orientation.
static void orient_emissions_pair(const ArgHmmMatrices &matrices, int start,
                                  const PhaseProbs &phase_pr,
                                  double **now, double **swap)
{
    const int nstates = max(matrices.nstates2, 1);
    for (int i=0; i<matrices.blocklen; i++) {
        map<int, vector<double> >::const_iterator it =
            phase_pr.probs.find(start + i);
        for (int k=0; k<nstates; k++) {
            double e = matrices.emit[i][k];
            if (it == phase_pr.probs.end()) {
                now[i][k] = swap[i][k] = e;
            } else {
                // un-average
                now[i][k] = 2.0 * e * it->second[k];
                swap[i][k] = 2.0 * e * (1.0 - it->second[k]);
            }
        }
    }
}


// Emissions of one block for the two phase states.  At the individual's
// sites in 'blocks' the two orientations are split between the phase states;
// elsewhere the current orientation is kept.  With point errors
// (point_prob > 0), at sites in a block each phase state's emission is
// averaged over the orientation being its own (1 - point_prob) or the other
// one (point_prob); the emissions of the orientations themselves are then
// saved in orient11[site] and orient12[site] (for each thread state) for
// the traceback.  'given' marks a block whose emissions the forward
// algorithm does not use; their saved emissions are set equal.
static void calc_phase_layer_emissions(
    int blocklen, int nstates, int start, double **now, double **swap,
    const Sequences *sequences, const PhaseBlocks &blocks,
    const vector<int> &block_site, int region_start, double point_prob,
    bool given, double **emit11, double **emit12,
    vector<vector<double> > *orient11, vector<vector<double> > *orient12)
{
    for (int i=0; i<blocklen; i++) {
        const int pos = start + i;
        int site = block_site[pos - region_start];
        for (int k=0; k<nstates; k++) {
            if (site < 0) {
                emit11[i][k] = emit12[i][k] = now[i][k];
            } else if (sequences->seqs[blocks.hap1][pos] ==
                       blocks.allele1[site]) {
                emit11[i][k] = now[i][k];
                emit12[i][k] = swap[i][k];
            } else {
                emit11[i][k] = swap[i][k];
                emit12[i][k] = now[i][k];
            }
        }

        if (site < 0 || point_prob == 0.0 || blocks.block[site] < 0)
            continue;
        vector<double> &o11 = (*orient11)[site];
        vector<double> &o12 = (*orient12)[site];
        o11.resize(nstates);
        o12.resize(nstates);
        for (int k=0; k<nstates; k++) {
            o11[k] = given ? 1.0 : emit11[i][k];
            o12[k] = given ? 1.0 : emit12[i][k];
            emit11[i][k] = (1.0 - point_prob) * o11[k] + point_prob * o12[k];
            emit12[i][k] = (1.0 - point_prob) * o12[k] + point_prob * o11[k];
        }
    }
}


// Forward algorithm over the doubled (thread state, phase state) space.
// See phase_thread.h.
void forward_phase(const LocalTrees *trees, const ArgModel *model,
                   ArgHmmMatrixIter *matrix_iter, PhaseProbs *phase_pr,
                   const Sequences *sequences, const PhaseBlocks &blocks,
                   const vector<int> &block_site,
                   const vector<char> &mix, double switch_prob,
                   double point_prob, ArgHmmForwardTable *forward,
                   ArgHmmForwardTable *forward2,
                   vector<vector<double> > *orient11,
                   vector<vector<double> > *orient12)
{
    vector<vector<double> > own11, own12;
    if (!orient11) {
        orient11 = &own11;
        orient12 = &own12;
    }
    if (point_prob > 0.0) {
        orient11->assign(blocks.size(), vector<double>());
        orient12->assign(blocks.size(), vector<double>());
    }
    const bool internal = true;
    const int start_coord = trees->start_coord;
    LineageCounts lineages(model->ntimes, model->num_pops());
    States states;
    ArgModel local_model;
    int mu_idx=0, rho_idx=0;
    double **fw = forward->get_table();
    double **fw2 = forward2->get_table();
    const char *mixp = &mix[0] - start_coord;  // indexed by position

    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        ArgHmmMatrices &matrices = matrix_iter->ref_matrices(phase_pr);
        int pos = matrix_iter->get_block_start();
        int blocklen = matrices.blocklen;
        int nstates = max(matrices.nstates2, 1);
        model->get_local_model(pos, local_model, &mu_idx, &rho_idx);

        forward->new_block(pos, pos+blocklen, matrices.nstates2);
        forward2->new_block(pos, pos+blocklen, matrices.nstates2);

        // emissions for the current and swapped orientation
        double **now = new_matrix<double>(blocklen, nstates);
        double **swap = new_matrix<double>(blocklen, nstates);
        orient_emissions_pair(matrices, pos, *phase_pr, now, swap);

        double **emit11 = new_matrix<double>(blocklen, nstates);
        double **emit12 = new_matrix<double>(blocklen, nstates);
        // the forward algorithm leaves out the emissions of the first
        // column, and of blocks where the branch has only one place to go
        // (in its internal mode)
        const bool given = (internal && matrices.nstates2 == 0);
        calc_phase_layer_emissions(blocklen, nstates, pos, now, swap,
                                   sequences, blocks, block_site,
                                   start_coord, point_prob, given, emit11,
                                   emit12, orient11, orient12);
        if (pos == start_coord && point_prob > 0.0 && block_site[0] >= 0 &&
            blocks.block[block_site[0]] >= 0) {
            (*orient11)[block_site[0]].assign(nstates, 1.0);
            (*orient12)[block_site[0]].assign(nstates, 1.0);
        }
        delete_matrix<double>(now, blocklen);
        delete_matrix<double>(swap, blocklen);
        double **e11 = emit11, **e12 = emit12;
        double **fw_block = &fw[pos], **fw2_block = &fw2[pos];
        const char *mix_block = &mixp[pos];

        matrices.states_model.get_coal_states(tree, states);
        lineages.count(tree, model->pop_tree, internal);

        if (pos == start_coord) {
            // prior of first state, split evenly between phase states
            int minage = matrices.states_model.minage;
            int subtree_root = tree->nodes[tree->root].child[0];
            if (subtree_root != -1)
                minage = max(minage, tree->nodes[subtree_root].age);
            calc_state_priors(states, &lineages, &local_model, fw[pos],
                              minage);
            for (int k=0; k<nstates; k++) {
                fw[pos][k] *= 0.5;
                fw2[pos][k] = fw[pos][k];
            }
        } else if (matrices.transmat_switch) {
            // first column of a new local tree: use the switch matrix
            arghmm_forward_switch_trans(fw[pos-1], fw[pos],
                                        matrices.transmat_switch);
            arghmm_forward_switch_trans(fw2[pos-1], fw2[pos],
                                        matrices.transmat_switch);
            mix_phase_layers(fw[pos], fw2[pos], nstates, mixp[pos],
                             switch_prob);
            double norm = 0.0;
            for (int k=0; k<nstates; k++) {
                fw[pos][k] *= emit11[0][k];
                fw2[pos][k] *= emit12[0][k];
                norm += fw[pos][k] + fw2[pos][k];
            }
            assert(norm > 0.0 && !isnan(norm) && !isinf(norm));
            for (int k=0; k<nstates; k++) {
                fw[pos][k] /= norm;
                fw2[pos][k] /= norm;
            }
        } else {
            // same state space as the previous block: continue from its
            // last column
            fw_block = &fw[pos-1];
            fw2_block = &fw2[pos-1];
            e11--;
            e12--;
            mix_block--;
            blocklen++;
        }

        arghmm_forward_block(model, tree, blocklen, states, lineages,
                             matrices.transmat, e11, fw_block,
                             fw2_block, e12, mix_block, switch_prob);

        delete_matrix<double>(emit11, matrices.blocklen);
        delete_matrix<double>(emit12, matrices.blocklen);
    }
}




// Sample (thread state, phase state) at column i, given the states already
// drawn at column i+1 and the transition probabilities 'trans[j]' of moving
// from thread state j at column i to the drawn thread state at i+1.
static void sample_column(const double *col1, const double *col2,
                          const double *trans, int nstates,
                          char mix_next, int phase_next, double switch_prob,
                          int *state, char *phase)
{
    double A[2*nstates];
    double p11 = phase_trans(mix_next, 0, phase_next, switch_prob);
    double p12 = phase_trans(mix_next, 1, phase_next, switch_prob);
    for (int j=0; j<nstates; j++) {
        A[j] = col1[j] * trans[j] * p11;
        A[nstates + j] = col2[j] * trans[j] * p12;
    }
    int x = sample(A, 2*nstates);
    *state = x % nstates;
    *phase = x / nstates;
    assert(trans[*state] != 0.0);
}


// Stochastic traceback over the doubled state space
static void traceback_phase(const LocalTrees *trees,
                            ArgHmmMatrixIter *matrix_iter,
                            double **fw, double **fw2,
                            const vector<char> &mix, double switch_prob,
                            int *path, char *phase)
{
    const char *mixp = &mix[0] - trees->start_coord;
    States states;

    // choose last column first
    matrix_iter->rbegin();
    int pos = trees->end_coord;
    {
        ArgHmmMatrices &mat = matrix_iter->ref_matrices();
        const int nstates = max(mat.nstates2, 1);
        double ones[nstates];
        for (int j=0; j<nstates; j++) ones[j] = 1.0;
        // no transition after the last column
        sample_column(fw[pos-1], fw2[pos-1], ones, nstates,
                      PHASE_MIX_REDRAW, 0, switch_prob,
                      &path[pos-1], &phase[pos-1]);
    }

    // iterate backward through blocks
    for (; matrix_iter->more(); matrix_iter->prev()) {
        ArgHmmMatrices &mat = matrix_iter->ref_matrices();
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        mat.states_model.get_coal_states(tree, states);
        const int nstates = max((int) states.size(), 1);
        pos -= mat.blocklen;

        // columns within the block, plus the last column of the previous
        // block if it shares this block's state space
        int first = pos;
        if (pos > trees->start_coord && !mat.transmat_switch)
            first = pos - 1;

        double trans[nstates];
        int last_k = -1;
        for (int i=pos+mat.blocklen-2; i>=first; i--) {
            int k = path[i+1];
            if (k != last_k) {
                for (int j=0; j<nstates; j++)
                    trans[j] = (states.size() == 0) ? 1.0 :
                        mat.transmat->get(tree, states, j, k);
                last_k = k;
            }
            sample_column(fw[i], fw2[i], trans, nstates, mixp[i+1],
                          phase[i+1], switch_prob, &path[i], &phase[i]);
        }

        // last column of the previous block via the switch matrix
        if (pos > trees->start_coord && mat.transmat_switch) {
            int i = pos - 1;
            const int nstates1 = max(mat.transmat_switch->nstates1, 1);
            double strans[nstates1];
            for (int j=0; j<nstates1; j++)
                strans[j] = mat.transmat_switch->get(j, path[i+1]);
            sample_column(fw[i], fw2[i], strans, nstates1, mixp[i+1],
                          phase[i+1], switch_prob, &path[i], &phase[i]);
        }
    }
}


// Sample the orientation of the alleles at each site of 'blocks' given the
// sampled path: the phase state, flipped at sampled point errors.
// orient11/orient12 hold the emissions of each orientation at sites in
// blocks (see calc_phase_layer_emissions).
static void sample_orientations(const PhaseBlocks &blocks, double point_prob,
                                const int *thread_path, const char *phase,
                                const vector<vector<double> > &orient11,
                                const vector<vector<double> > &orient12,
                                vector<int> &orient)
{
    orient.resize(blocks.size());
    for (int k=0; k<blocks.size(); k++) {
        int pos = blocks.sites[k];
        orient[k] = phase[pos];
        if (point_prob == 0.0 || blocks.block[k] < 0)
            continue;
        int state = thread_path[pos];
        double e11 = orient11[k][state], e12 = orient12[k][state];
        double own = (phase[pos] == 0) ? e11 : e12;
        double other = (phase[pos] == 0) ? e12 : e11;
        double err = point_prob * other;
        if (frand() * ((1.0 - point_prob) * own + err) < err)
            orient[k] = 1 - phase[pos];
    }
}


int resample_arg_leaf_phase(const ArgModel *model, Sequences *sequences,
                            LocalTrees *trees, int node,
                            const PhaseBlocks &blocks, double switch_prob,
                            double point_prob)
{
    const bool internal = true;
    const int start_coord = trees->start_coord;
    const int seqlen = trees->length();
    assert(trees->seqids[node] == blocks.hap1 ||
           trees->seqids[node] == blocks.hap2);
    assert(point_prob >= 0.0 && point_prob <= 0.5);

    // remove the leaf's branch (as in resample_arg_leaf)
    const int maxtime = model->get_removed_root_time();
    int *removal_path = new int [trees->get_num_trees()];
    sample_arg_removal_leaf_path(trees, node, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime, model->pop_tree);
    delete [] removal_path;

    int minage = sequences->ages[trees->seqids[node]];
    if (minage > 0) {
        for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it)
            it->tree->nodes[node].age = minage;
    }

    // phase transitions at each column, and which phase block site (if any)
    // each column holds
    vector<char> mix(seqlen, PHASE_MIX_NONE);
    vector<int> block_site(seqlen, -1);
    for (int k=0; k<blocks.size(); k++) {
        int i = blocks.sites[k] - start_coord;
        mix[i] = blocks.same_block(k) ? PHASE_MIX_SWITCH : PHASE_MIX_REDRAW;
        block_site[i] = k;
    }

    // forward pass.  ARGweaver computes emissions for both orientations of
    // the pair only in unphased mode, so switch it on meanwhile.
    ArgHmmForwardTable forward(start_coord, seqlen);
    ArgHmmForwardTable forward2(start_coord, seqlen);
    ArgHmmMatrixIter matrix_iter(model, sequences, trees);
    matrix_iter.set_internal(internal, minage);
    vector<vector<double> > orient11, orient12;
    {
        ArgModel *model_rw = const_cast<ArgModel*>(model);
        const bool unphased = model->unphased;
        model_rw->unphased = true;
        PhaseProbs phase_pr(trees->seqids[node], node, sequences, trees,
                            model);
        assert(phase_pr.treemap2 >= 0);
        forward_phase(trees, model, &matrix_iter, &phase_pr, sequences,
                      blocks, block_site, mix, switch_prob, point_prob,
                      &forward, &forward2, &orient11, &orient12);
        model_rw->unphased = unphased;
    }

    // sample thread and phase together, then point errors
    int *thread_path_alloc = new int [seqlen];
    int *thread_path = &thread_path_alloc[-start_coord];
    vector<char> phase_alloc(seqlen);
    char *phase = &phase_alloc[0] - start_coord;
    ArgHmmMatrixIter matrix_iter2(model, NULL, trees);
    matrix_iter2.set_internal(internal, minage);
    traceback_phase(trees, &matrix_iter2, forward.get_table(),
                    forward2.get_table(), mix, switch_prob,
                    thread_path, phase);
    vector<int> orient;
    sample_orientations(blocks, point_prob, thread_path, phase, orient11,
                        orient12, orient);

    // swap alleles to match the sampled orientations
    int nswapped = 0;
    for (int k=0; k<blocks.size(); k++) {
        int pos = blocks.sites[k];
        bool is11 = (sequences->seqs[blocks.hap1][pos] == blocks.allele1[k]);
        if (is11 != (orient[k] == 0)) {
            sequences->switch_alleles(pos, blocks.hap1, blocks.hap2);
            nswapped++;
        }
    }

    // sample recombination points and add the thread to the ARG
    vector<int> recomb_pos;
    vector<Spr> recombs;
    sample_recombinations(trees, model, &matrix_iter2, thread_path,
                          recomb_pos, recombs, internal);
    add_arg_thread_path(trees, matrix_iter.states_model, model->ntimes,
                        thread_path, recomb_pos, recombs, model->pop_tree);

    delete [] thread_path_alloc;
    return nswapped;
}

} // namespace argweaver
