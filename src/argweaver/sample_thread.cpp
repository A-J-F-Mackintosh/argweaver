// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "hmm.h"
#include "local_tree.h"
#include "logging.h"
#include "matrices.h"
#include "model.h"
#include "recomb.h"
#include "sample_thread.h"
#include "sequences.h"
#include "sequences.h"
#include "states.h"
#include "thread.h"
#include "trans.h"
#include "total_prob.h"


namespace argweaver {

using namespace std;


//=============================================================================
// Forward algorithm for thread path

// Mix the two read-backed phase layers of a forward column at a
// heterozygous site of the individual being phased (see phase_thread.h).
void mix_phase_layers(double *col_a, double *col_b, int nstates, char mix,
                      double switch_prob)
{
    if (mix == PHASE_MIX_SWITCH) {
        for (int k=0; k<nstates; k++) {
            double a = col_a[k], b = col_b[k];
            col_a[k] = (1.0 - switch_prob) * a + switch_prob * b;
            col_b[k] = switch_prob * a + (1.0 - switch_prob) * b;
        }
    } else if (mix == PHASE_MIX_REDRAW) {
        for (int k=0; k<nstates; k++)
            col_a[k] = col_b[k] = 0.5 * (col_a[k] + col_b[k]);
    }
}


// compute one block of forward algorithm with compressed transition matrices
// NOTE: first column of forward table should be pre-populated
//
// If fw2 is given, the state space is doubled by a read-backed phase state:
// fw/emit are the forward table and emissions for phase state "1-1", and
// fw2/emit2 for "1-2".  mix[i] says how the phase state changes at column i.
void arghmm_forward_block(const ArgModel *model,
                          const LocalTree *tree,
                          const int blocklen, const States &states,
                          const LineageCounts &lineages,
                          const TransMatrix *matrix,
                          const double* const *emit, double **fw,
                          double **fw2, const double* const *emit2,
                          const char *mix, double switch_prob)
{
    const int nlayers = (fw2 == NULL) ? 1 : 2;
    const int nstates = states.size();
    const LocalNode *nodes = tree->nodes;
    const int ntimes = model->ntimes;

    //  handle internal branch resampling special cases
    int minage = matrix->minage;
    int maintree_root = 0;
    if (matrix->internal) {
        maintree_root = nodes[tree->root].child[1];

        if (nstates == 0) {
            // handle fully given case
            for (int i=1; i<blocklen; i++) {
                fw[i][0] = fw[i-1][0];
                if (fw2) {
                    fw2[i][0] = fw2[i-1][0];
                    mix_phase_layers(fw[i], fw2[i], 1, mix[i], switch_prob);
                    double norm = fw[i][0] + fw2[i][0];
                    fw[i][0] /= norm;
                    fw2[i][0] /= norm;
                }
            }
            return;
        }
    }

    // get max time
    int maxtime = 0;
    for (int k=0; k<nstates; k++)
        if (maxtime < states[k].time)
            maxtime = states[k].time;

    int numpath = model->num_pop_paths();
    int numpath_per_time[ntimes];
    int paths_per_time[ntimes][numpath];
    int path_map[states.size()];
    int max_numpath = 1;
    if (numpath > 1) {
        for (int i=0; i < ntimes; i++) {
            numpath_per_time[i]=0;
            for (int j=0; j < numpath; j++)
                paths_per_time[i][j]=0;
        }
        for (unsigned int i=0; i < states.size(); i++) {
            int t = states[i].time;
            int p = states[i].pop_path;
            int j=0;
            for ( ; j < numpath_per_time[t]; j++) {
                if (model->paths_equal(p, paths_per_time[t][j],
                                       minage, t))
                    break;
            }
            path_map[i] = j;
            if (j == numpath_per_time[t])
                paths_per_time[t][numpath_per_time[t]++] = p;
        }
        for (int i=0; i < ntimes; i++)
            if (numpath_per_time[i] > max_numpath)
                max_numpath = numpath_per_time[i];
    } else {
        for (int i=0; i < ntimes; i++) {
            numpath_per_time[i] = 1;
            paths_per_time[i][0] = 0;
        }
        for (unsigned int i=0; i < states.size(); i++)
            path_map[i] = 0;
        max_numpath = 1;
    }

    // get branch ages
    // set ages1[i] to age of each branch
    // set ages2[i] to age of each branch's parent
    // set indexes[i] to index for state (node[i], ages1[i])
    int ages1[tree->nnodes];
    int ages2[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        ages1[i] = max(nodes[i].age, minage);
        if (matrix->internal)
            ages2[i] = (i == maintree_root || i == tree->root) ?
                maxtime : nodes[nodes[i].parent].age;
        else
            ages2[i] = (i == tree->root) ? maxtime : nodes[nodes[i].parent].age;
    }

    // compute ntimes*ntimes and ntime*nstates temp matrices
    double tmatrix[ntimes-1][max_numpath][ntimes-1][max_numpath];
    for (int b=0; b<ntimes-1; b++) {
        for (int pb=0; pb < numpath_per_time[b]; pb++) {
            for (int a=0; a<ntimes-1; a++) {
                for (int pa=0; pa < numpath_per_time[a]; pa++) {
                    tmatrix[b][pb][a][pa] =
                        matrix->get_time(a, b, 0,
                                         paths_per_time[a][pa],
                                         paths_per_time[b][pb],
                                         -1, minage, false);
                    //                    printf("tmatrix %i %i = %e\n", a, b, tmatrix[pa][pb][a][b]);
                    assert(!isnan(tmatrix[b][pb][a][pa]));
                    assert(!isinf(tmatrix[b][pb][a][pa]));
                }
            }
        }
    }

    // take advantage of fact that same branch case is only special
    // if path a and path b are same; otherwise there must be recomb
    // on branch being threaded and same branch case is not special
    double tmatrix2[nstates][ntimes];
    for (int k=0; k<nstates; k++) {
        for (int a=0; a < ntimes; a++) tmatrix2[k][a]=0.0;
        const int b = states[k].time;
        const int node2 = states[k].node;
        const int c = nodes[node2].age;
        const int p = states[k].pop_path;
        const int pc = nodes[node2].pop_path;
        for (int a=ages1[node2]; a <= ages2[node2]; a++) {
            tmatrix2[k][a] =
                matrix->get_time(a, b, c, p, p, pc, minage, true, k) -
                matrix->get_time(a, b, 0, p, p, -1, minage, false);
            /*            printf("tmatrix2\t%i\t%i\t%e\t%e\t%e\n", a, k, tmatrix2[k][a],
                   matrix->get_time(a, b, c, p, p, pc, minage, true, k),
                   matrix->get_time(a, b, 0, p, p, -1, minage, false));*/
            if (isnan(tmatrix2[k][a]) || isinf(tmatrix2[k][a]) || tmatrix2[k][a] < 0) {
                printf("a=%i k=%i b=%i node2=%i c=%i p=%i pc=%i\n",
                       a, k, b, node2, c, p, pc);
                assert(false);
            }
        }
    }

    // there is one more special case for different path, same time, same node
    double tmatrix3[nstates][max_numpath];
    if (max_numpath > 1) {
        for (int k=0; k < nstates; k++) {
            for (int i=0; i <max_numpath; i++) tmatrix3[k][i]=0.0;
            int b = states[k].time;
            const int pb = states[k].pop_path;
            for (int j=0; j < numpath_per_time[b]; j++) {
                int pa = paths_per_time[b][j];
                if (!model->paths_equal(pa, pb, minage, states[k].time)) {
                    tmatrix3[k][j] =
                        ( matrix->get_time(b, b, -1, pa, pb, -1, minage, true, k) -
                          matrix->get_time(b, b, -1, pa, pb, -1, minage, false, k));
                }
            }
        }
    }

    NodeStateLookup state_lookup(states, minage, model->pop_tree);
    int max_idx = ntimes*nstates + max_numpath*nstates;
    int nextState[max_idx];
    int idx=0;
    int age1_state[nstates];
    for (int k=0; k<nstates; k++) {
        const int b = states[k].time;
        const int node2 = states[k].node;
        int age1 = ages1[node2];
        const int age2 = ages2[node2];
        const int path1 = tree->nodes[node2].pop_path;
        const int path2 = states[k].pop_path;
        int j=state_lookup.lookup_idx(node2, age1, path2);
        while (j < 0 && age1 <= age2) {
            age1++;
            j = state_lookup.lookup_idx(node2, age1, path2);
        }
        age1_state[k] = age1;
        for (int a=age1; a <= age2; a++, j++) {
            int j_state = state_lookup.lookup_by_idx(j);
            if (j_state >= 0 &&
                (model->pop_tree == NULL || a >= b ||
                 model->paths_equal(path1,
                                    path2, a, b))) {
                nextState[idx++]=j_state;
            } else nextState[idx++] = -1;
        }
        // this setion accounts for self-recombinations that change paths
        // (same node, same time, different path)
        if (max_numpath > 1) {
            for (int pa=0; pa < numpath_per_time[b]; pa++) {
                int path_a = paths_per_time[b][pa];
                if (!model->paths_equal(path_a, path2, minage, b))
                    nextState[idx++] = state_lookup.lookup(node2, b, path_a);
                else nextState[idx++] = -1;
            }
        }
    }
    assert(idx <= max_idx);


    double tmatrix_fgroups[max_numpath][ntimes];
    double fgroups[max_numpath][ntimes];
    for (int i=1; i<blocklen; i++) {
        // apply the thread transitions to each phase layer
        for (int layer=0; layer<nlayers; layer++) {
            const double *col1 = (layer == 0) ? fw[i-1] : fw2[i-1];
            double *col2 = (layer == 0) ? fw[i] : fw2[i];
            idx = 0;

            // precompute the fgroup sums
            for (int p=0; p < max_numpath; p++)
                fill(fgroups[p], fgroups[p]+ntimes, 0.0);
            for (int j=0; j<nstates; j++) {
                const int a = states[j].time;
                fgroups[path_map[j]][a] += col1[j];
                assert(!isinf(col1[j]));
            }

            // multiply tmatrix and fgroups together
            for (int b=0; b<ntimes-1; b++) {
                for (int pb=0; pb < numpath_per_time[b]; pb++) {
                    double sum = 0.0;
                    for (int a=0; a<ntimes-1; a++) {
                        for (int pa=0; pa < numpath_per_time[a]; pa++) {
                            sum += tmatrix[b][pb][a][pa] * fgroups[pa][a];
                        }
                    }
                    tmatrix_fgroups[pb][b] = sum;
                }
            }

            // fill in one column of forward table (before emissions)
            for (int k=0; k<nstates; k++) {
                const int b = states[k].time;
                const int node2 = states[k].node;
                const int age2 = ages2[node2];
                double sum = tmatrix_fgroups[path_map[k]][b];

                // same branch case
                for (int a=age1_state[k]; a <= age2; a++) {
                    int j_state = nextState[idx++];
                    if (j_state >= 0 && col1[j_state] > 0) {
                        sum += tmatrix2[k][a] * col1[j_state];
                    }
                }
                // this setion accounts for self-recombinations that change paths
                // (same node, same time, different path)
                if (max_numpath > 1) {
                    for (int pa=0; pa < numpath_per_time[b]; pa++) {
                        int j_state = nextState[idx++];
                        if (j_state >= 0 && col1[j_state] > 0) {
                            sum += tmatrix3[k][pa] * col1[j_state];
                        }
                    }
                }
                col2[k] = sum;
            }
        }

        // change of read-backed phase state
        if (nlayers == 2)
            mix_phase_layers(fw[i], fw2[i], nstates, mix[i], switch_prob);

        // apply emissions
        double norm = 0.0;
        for (int layer=0; layer<nlayers; layer++) {
            double *col2 = (layer == 0) ? fw[i] : fw2[i];
            const double *emit_col = (layer == 0) ? emit[i] : emit2[i];
            for (int k=0; k<nstates; k++) {
                col2[k] *= emit_col[k];
                norm += col2[k];
                if (isnan(col2[k]))
                    assert(false);
            }
        }
        assert(norm > 0);
        assert(!isnan(norm));
        assert(!isinf(norm));

        // normalize column for numerical stability
        for (int k=0; k<nstates; k++)
            fw[i][k] /= norm;
        if (nlayers == 2)
            for (int k=0; k<nstates; k++)
                fw2[i][k] /= norm;
    }
}



// compute one block of forward algorithm with compressed transition matrices
// NOTE: first column of forward table should be pre-populated
// This can be used for testing
void arghmm_forward_block_slow(const LocalTree *tree, const int ntimes,
                               const int blocklen, const States &states,
                               const LineageCounts &lineages,
                               const TransMatrix *matrix,
                               const double* const *emit, double **fw)
{
    const int nstates = states.size();
    if (nstates == 0) {
        // handle fully given case
        for (int i=1; i<blocklen; i++)
            fw[i][0] = fw[i-1][0];
        return;
    }
    // get transition matrix
    double **transmat = new_matrix<double>(nstates, nstates);
    for (int k=0; k<nstates; k++)
        for (int j=0; j<nstates; j++)
            transmat[j][k] = matrix->get(tree, states, j, k);

    // fill in forward table
    for (int i=1; i<blocklen; i++) {
        const double *col1 = fw[i-1];
        double *col2 = fw[i];
        double norm = 0.0;

        for (int k=0; k<nstates; k++) {
            double sum = 0.0;
            for (int j=0; j<nstates; j++)
                sum += col1[j] * transmat[j][k];
            col2[k] = sum * emit[i][k];
            norm += col2[k];
        }

        // normalize column for numerical stability
        for (int k=0; k<nstates; k++)
            col2[k] /= norm;
    }

    // cleanup
    delete_matrix<double>(transmat, nstates);
}




// apply the switch matrix transitions to one column (no emissions)
void arghmm_forward_switch_trans(const double *col1, double* col2,
                                 const TransMatrixSwitch *matrix)
{
    // if state space is size zero, we still treat it as size 1
    const int nstates1 = max(matrix->nstates1, 1);
    const int nstates2 = max(matrix->nstates2, 1);
    //    printf("nstates1=%i nstates2=%i\n", nstates1, nstates2);

    // initialize all entries in col2 to 0
    for (int k=0; k<nstates2; k++)
        col2[k] = 0.0;

    // add deterministic transitions
    for (int j=0; j<nstates1; j++) {
        int k = matrix->determ[j];
        if (k != -1 && matrix->recombsrc[j] < 0 && matrix->recoalsrc[j] < 0) {
            col2[k] += col1[j] * matrix->determprob[j];
            assert(!isnan(col2[k]));
        }
    }

    // add recombination and recoalescing transitions
    for (int j=0; j < nstates1; j++) {
        if (matrix->recombsrc[j] >= 0) {
            assert(matrix->recoalsrc[j] < 0);
            for (int k=0; k<nstates2; k++) {
                double val = matrix->get(j, k);
                if (val > 0) {
                    col2[k] += col1[j] * val;
                    assert(!isnan(col2[k]));
                    //                    printf("recombsrc %i %i %e %e\n", j, k, val, col2[k]);
                }
            }
        }
    }
    for (int j=0; j < nstates1; j++) {
        if (matrix->recoalsrc[j] >= 0) {
            assert(matrix->recombsrc[j] < 0);
            for (int k=0; k<nstates2; k++) {
                double val = matrix->get(j, k);
                if (val > 0) {
                    col2[k] += col1[j] * val;
                    assert(!isnan(col2[k]));
                    //                    printf("recoalsrc %i %i %e %e\n", j, k, val, col2[k]);
                }
            }
        }
    }
}


// run forward algorithm for one column of the table
// use switch matrix
void arghmm_forward_switch(const double *col1, double* col2,
                           const TransMatrixSwitch *matrix,
                           const double *emit)
{
    const int nstates2 = max(matrix->nstates2, 1);
    arghmm_forward_switch_trans(col1, col2, matrix);

    double norm = 0.0;
    for (int k=0; k<nstates2; k++) {
        col2[k] *= emit[k];
        norm += col2[k];
    }
    assert(norm != 0.0);
    assert(!isnan(norm));
    assert(!isinf(norm));

#ifdef DEBUG
    // assert that probability is valid
    double top = max_array(col2, nstates2);
    if (top <= 0.0) {
        for (int i=0; i<nstates2; i++) {
            printf("col2[%d] = %f\n", i, col2[i]);
        }
        assert(false);
    }
#endif

    // normalize column for numerical stability
    for (int k=0; k<nstates2; k++)
        col2[k] /= norm;
    if (isnan(norm))
        assert(false);
}



// Run forward algorithm for all blocks
void arghmm_forward_alg(const LocalTrees *trees, const ArgModel *model,
    const Sequences *sequences, ArgHmmMatrixIter *matrix_iter,
    ArgHmmForwardTable *forward, PhaseProbs *phase_pr,
    bool prior_given, bool internal, bool slow)
{
    LineageCounts lineages(model->ntimes, model->num_pops());
    States states;
    ArgModel local_model;
    int mu_idx=0, rho_idx=0;
    LocalTree *tree;
#ifdef DEBUG
    LocalTree *last_tree = NULL;
#endif

    double **fw = forward->get_table();
    // forward algorithm over local trees
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        // get block information
        tree = matrix_iter->get_tree_spr()->tree;
        ArgHmmMatrices &matrices = matrix_iter->ref_matrices(phase_pr);
        int pos = matrix_iter->get_block_start();
        int blocklen = matrices.blocklen;
        model->get_local_model(pos, local_model, &mu_idx, &rho_idx);
        double **emit = matrices.emit;

        // allocate the forward table
        if (pos > trees->start_coord || !prior_given)
            forward->new_block(pos, pos+matrices.blocklen, matrices.nstates2);
        double **fw_block = &fw[pos];

        matrices.states_model.get_coal_states(tree, states);
        lineages.count(tree, model->pop_tree, internal);

        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == trees->start_coord) {
            // calculate prior of first state
            int minage = matrices.states_model.minage;
            if (!prior_given) {
                if (internal) {
                    int subtree_root = tree->nodes[tree->root].child[0];
                    if (subtree_root != -1)
                        minage = max(minage, tree->nodes[subtree_root].age);
                }
                calc_state_priors(states, &lineages, &local_model,
                                  fw[pos], minage);
            }
        } else if (matrices.transmat_switch) {
            // perform one column of forward algorithm with transmat_switch
            arghmm_forward_switch(fw[pos-1], fw[pos],
                matrices.transmat_switch, matrices.emit[0]);
        } else {
            // we are still inside the same ARG block, therefore the
            // state-space does not change and no switch matrix is needed
            fw_block = &fw[pos-1];
            emit--;
            blocklen++;
        }

        int nstates = max(matrices.transmat->nstates, 1);
        double top = max_array(fw_block[0], nstates);
        //        printf("%i %i top =%f nstates=%i\n", count, pos, top, nstates);
        for (int i=0; i < nstates; i++) assert(!isnan(fw_block[0][i]));
        assert(!isnan(top));
        assert(top > 0.0);

        // calculate rest of block
        if (slow)
            arghmm_forward_block_slow(tree, model->ntimes, blocklen,
                                      states, lineages, matrices.transmat,
                                      emit, fw_block);
        else
            arghmm_forward_block(model, tree, blocklen,
                                 states, lineages, matrices.transmat,
                                 emit, fw_block);

        // safety check
        double top2 = max_array(fw[pos + matrices.blocklen - 1], nstates);
        //        printf("%i %i top2=%f\n", count++, pos, top2);
        assert(top2 > 0.0);
#ifdef DEBUG
        last_tree = tree;
#endif
    }
}




//=============================================================================
// Sample thread paths



double sample_hmm_posterior(
    int blocklen, const LocalTree *tree, const States &states,
    const TransMatrix *matrix, const double *const *fw, int *path)
{
    // NOTE: path[blocklen-1] must already be sampled

    const int nstates = max(states.size(), (size_t)1);
    double A[nstates];
    double trans[nstates];
    int last_k = -1;
    double lnl = 0.0;

    // recurse
    for (int i=blocklen-2; i>=0; i--) {
        int k = path[i+1];

        // recompute transition probabilities if state (k) changes
        if (k != last_k) {
            for (int j=0; j<nstates; j++)
                trans[j] = matrix->get(tree, states, j, k);
            last_k = k;
        }

        for (int j=0; j<nstates; j++)
            A[j] = fw[i][j] * trans[j];
        path[i] = sample(A, nstates);
        //lnl += log(A[path[i]]);

        // DEBUG
        assert(trans[path[i]] != 0.0);
    }

    return lnl;
}


int sample_hmm_posterior_step(const TransMatrixSwitch *matrix,
                              const double *col1, int state2)
{
    const int nstates1 = max(matrix->nstates1, 1);
    double A[nstates1];

    for (int j=0; j<nstates1; j++)
        A[j] = col1[j] * matrix->get(j, state2);
    int k = sample(A, nstates1);

    // DEBUG
    assert(matrix->get(k, state2) != 0.0);
    return k;
}


double stochastic_traceback(
    const LocalTrees *trees, const ArgModel *model,
    ArgHmmMatrixIter *matrix_iter,
    double **fw, int *path, bool last_state_given, bool internal)
{
    States states;
    double lnl = 0.0;
    /*    printf("stochastic_traceback last_state_given=%i internal=%i\n",
          (int)last_state_given, (int)internal);*/

    // choose last column first
    matrix_iter->rbegin();
    int pos = trees->end_coord;

    if (!last_state_given) {
        ArgHmmMatrices &mat = matrix_iter->ref_matrices();
        const int nstates = max(mat.nstates2, 1);
        path[pos-1] = sample(fw[pos-1], nstates);
        lnl = fw[pos-1][path[pos-1]];
    }

    // iterate backward through blocks
    for (; matrix_iter->more(); matrix_iter->prev()) {
        ArgHmmMatrices &mat = matrix_iter->ref_matrices();
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        mat.states_model.get_coal_states(tree, states);
        pos -= mat.blocklen;

        lnl += sample_hmm_posterior(mat.blocklen, tree, states,
                                    mat.transmat, &fw[pos], &path[pos]);

        // fill in last col of next block
        if (pos > trees->start_coord) {
            if (mat.transmat_switch) {
                // use switch matrix
                int i = pos - 1;
                path[i] = sample_hmm_posterior_step(
                    mat.transmat_switch, fw[i], path[i+1]);
                lnl += log(fw[i][path[i]] *
                           mat.transmat_switch->get(path[i], path[i+1]));
            } else {
                // use normal matrix
                lnl += sample_hmm_posterior(2, tree, states,
                    mat.transmat, &fw[pos-1], &path[pos-1]);
            }
        }
    }

    return lnl;
}



//=============================================================================
// ARG sampling


// sample the thread of the last chromosome
void sample_arg_thread(const ArgModel *model, Sequences *sequences,
                       LocalTrees *trees, int new_chrom)
{
    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];
    int start_pop = sequences->get_pop(new_chrom);

    assert_trees(trees, model->pop_tree);

    // gives probability of current phasing for each heterozygous site;
    // will only be filled in if model->unphased==true
    PhaseProbs phase_pr(new_chrom, trees->get_num_leaves(),
			sequences, trees, model);
    if (model->unphased)
      printf("treemap = %i %i\n", phase_pr.treemap1, phase_pr.treemap2);

    // build matrices
    ArgHmmMatrixIter matrix_iter(model, sequences, trees, new_chrom);
    matrix_iter.set_start_pop(start_pop);

    // compute forward table
    Timer time;
    arghmm_forward_alg(trees, model, sequences, &matrix_iter, &forward,
		       model->unphased ? &phase_pr : NULL);
    int nstates = get_num_coal_states(trees->front().tree, model->ntimes);
    printTimerLog(time, LOG_LOW,
                  "forward (%3d states, %6d blocks):",
                  nstates, trees->get_num_trees());

    // traceback
    time.start();
    double **fw = forward.get_table();
    ArgHmmMatrixIter matrix_iter2(model, NULL, trees, new_chrom);
    matrix_iter2.set_start_pop(start_pop);
    stochastic_traceback(trees, model, &matrix_iter2, fw, thread_path);
    printTimerLog(time, LOG_LOW,
                  "trace:                              ");

    time.start();

    if (model->unphased)
	phase_pr.sample_phase(thread_path);

    // sample recombination points
    vector<int> recomb_pos;
    vector<Spr> recombs;
    sample_recombinations(trees, model, &matrix_iter2,
                          thread_path, recomb_pos, recombs);
    assert_trees(trees, model->pop_tree);

    // add thread to ARG
    //    matrix_iter.states_model.set_start_pop(start_pop);
    add_arg_thread(trees, matrix_iter.states_model,
                   model->ntimes, thread_path, new_chrom,
                   recomb_pos, recombs, model->pop_tree);
    assert_trees(trees, model->pop_tree);
    printTimerLog(time, LOG_LOW,
                  "add thread:                         ");

    // clean up
    delete [] thread_path_alloc;
}


// sample the thread of the internal branch
void sample_arg_thread_internal(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    int minage, PhaseProbs *phase_pr)
{
    const bool internal = true;

    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];

    // build matrices
    ArgHmmMatrixIter matrix_iter(model, sequences, trees);
    matrix_iter.set_internal(internal, minage);

    if (phase_pr != NULL)
        printLog(LOG_HIGH, "treemap = %i %i\n",
                 phase_pr->treemap1, phase_pr->treemap2);

    // compute forward table
    Timer time;
    arghmm_forward_alg(trees, model, sequences, &matrix_iter, &forward,
                       phase_pr, false, internal);
    int nstates = get_num_coal_states_internal(
           trees->front().tree, model->ntimes, minage);
    printTimerLog(time, LOG_LOW,
                  "forward (%3d states, %6d blocks):",
                  nstates, trees->get_num_trees());

    // traceback
    time.start();
    double **fw = forward.get_table();
    ArgHmmMatrixIter matrix_iter2(model, NULL, trees);
    matrix_iter2.set_internal(internal, minage);
    stochastic_traceback(trees, model, &matrix_iter2, fw, thread_path,
                         false, internal);
    printTimerLog(time, LOG_LOW,
                  "trace:                              ");

    if (phase_pr != NULL)
        phase_pr->sample_phase(thread_path);

    // sample recombination points
    time.start();
    vector<int> recomb_pos;
    vector<Spr> recombs;
    sample_recombinations(trees, model, &matrix_iter2,
                          thread_path, recomb_pos, recombs, internal);

    // add thread to ARG
    add_arg_thread_path(trees, matrix_iter.states_model,
                        model->ntimes, thread_path,
                        recomb_pos, recombs, model->pop_tree);
    printTimerLog(time, LOG_LOW,
                  "add thread:                         ");

    // clean up
    delete [] thread_path_alloc;
}



// sample the thread of the last chromosome, conditioned on a given
// start and end state
/*void cond_sample_arg_thread(const ArgModel *model, const Sequences *sequences,
                            LocalTrees *trees, int new_chrom,
                            State start_state, State end_state)
{
    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    States states;
    double **fw = forward.get_table();
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());

    // fill in first column of forward table
    matrix_list.begin();
    matrix_list.get_coal_states(states);
    forward.new_block(matrix_list.get_block_start(),
                      matrix_list.get_block_end(), states.size());
    int j = find_vector(states, start_state);
    assert(j != -1);
    double *col = fw[trees->start_coord];
    fill(col, col + states.size(), 0.0);
    col[j] = 1.0;

    // compute forward table
    time.start();
    arghmm_forward_alg(trees, model, sequences, &matrix_list, &forward, NULL,
		       true);
    int nstates = get_num_coal_states(trees->front().tree, model->ntimes);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           nstates, trees->get_num_trees());

    // fill in last state of traceback
    matrix_list.rbegin();
    matrix_list.get_coal_states(states);
    thread_path[trees->end_coord-1] = find_vector(states, end_state);
    assert(thread_path[trees->end_coord-1] != -1);

    // traceback
    time.start();
    stochastic_traceback(trees, model, &matrix_list, fw, thread_path, true);
    printf("trace:       %e s\n", time.time());
    assert(fw[trees->start_coord][thread_path[trees->start_coord]] == 1.0);


    // sample recombination points
    time.start();
    vector<int> recomb_pos;
    vector<Spr> recombs;
    sample_recombinations(trees, model, &matrix_list,
                          thread_path, recomb_pos, recombs);

    // add thread to ARG
    add_arg_thread(trees, matrix_list.states_model,
                   model->ntimes, thread_path, new_chrom,
                   recomb_pos, recombs, model->pop_tree);

    printf("add thread:  %e s\n", time.time());

    // clean up
    delete [] thread_path_alloc;
    } */



// sample the thread of the last chromosome, conditioned on a given
// start and end state
void cond_sample_arg_thread_internal(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    const State start_state, const State end_state)
{
    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    States states;
    double **fw = forward.get_table();
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];
    const bool internal = true;
    bool prior_given = true;
    bool last_state_given = true;
    /*    printf("cond_sample_arg_thread_internal start_state=(%i,%i) end_state=(%i,%i)\n",
           start_state.node, start_state.time,
           end_state.node, end_state.time);*/

    assert_trees(trees, model->pop_tree, true);

    // build matrices
    ArgHmmMatrixIter matrix_iter(model, sequences, trees);
    matrix_iter.set_internal(internal);

    // fill in first column of forward table
    matrix_iter.begin();
    matrix_iter.get_coal_states(states);
    forward.new_block(matrix_iter.get_block_start(),
                      matrix_iter.get_block_end(), states.size());

    if (states.size() > 0) {
        if (!start_state.is_null()) {
            // find start state
            LocalTree *first_tree = trees->front().tree;
            int subtree_root = first_tree->nodes[first_tree->root].child[0];
            assert(subtree_root != -1);
            int minage = first_tree->nodes[subtree_root].age;
            int j = find_state(states, start_state, model, minage);
            assert(j != -1);
            double *col = fw[trees->start_coord];
            fill(col, col + states.size(), 0.0);
            col[j] = 1.0;
        } else {
            // open ended, sample start state
            prior_given = false;
        }
    } else {
        // fully specified tree
        fw[trees->start_coord][0] = 1.0;
    }

    // compute forward table
    Timer time;
    arghmm_forward_alg(trees, model, sequences, &matrix_iter, &forward, NULL,
                       prior_given, internal);

    // TODO: Check that we don't need more arguments here!
    int nstates = get_num_coal_states_internal(
         trees->front().tree, model->ntimes);
    printTimerLog(time, LOG_LOW,
                  "forward (%3d states, %6d blocks):",
                  nstates, trees->get_num_trees());

    // fill in last state of traceback
    matrix_iter.rbegin();
    matrix_iter.get_coal_states(states);
    if (states.size() > 0) {
        if (!end_state.is_null()) {
            LocalTree *last_tree = trees->back().tree;
            int subtree_root = last_tree->nodes[last_tree->root].child[0];
            assert(subtree_root != -1);
            int minage = last_tree->nodes[subtree_root].age;
            thread_path[trees->end_coord-1] =
                find_state(states, end_state, model, minage);
            assert(thread_path[trees->end_coord-1] != -1);
        } else {
            // sample end start
            last_state_given = false;
        }
    } else {
        // fully specified tree
        thread_path[trees->end_coord-1] = 0;
    }

    // traceback
    time.start();
    ArgHmmMatrixIter matrix_iter2(model, NULL, trees);
    matrix_iter2.set_internal(internal);
    stochastic_traceback(trees, model, &matrix_iter2, fw, thread_path,
                         last_state_given, internal);
    printTimerLog(time, LOG_LOW,
                  "trace:                              ");
    if (!start_state.is_null())
        assert(fw[trees->start_coord][thread_path[trees->start_coord]] == 1.0);

    // sample recombination points
    time.start();
    vector<int> recomb_pos;
    vector<Spr> recombs;
    sample_recombinations(trees, model, &matrix_iter2,
                          thread_path, recomb_pos, recombs, internal);

    // add thread to ARG
    assert_trees(trees, model->pop_tree, true);
    add_arg_thread_path(trees, matrix_iter.states_model,
                        model->ntimes, thread_path,
                        recomb_pos, recombs, model->pop_tree);
    assert_trees(trees, model->pop_tree);
    printTimerLog(time, LOG_LOW,
                  "add thread:                         ");

    // clean up
    delete [] thread_path_alloc;
}



// resample the threading of one chromosome
void resample_arg_thread(const ArgModel *model, Sequences *sequences,
                         LocalTrees *trees, int chrom)
{
    // remove chromosome from ARG and resample its thread
    remove_arg_thread(trees, chrom, model);
    sample_arg_thread(model, sequences, trees, chrom);
}



//=============================================================================
// C interface

/*extern "C" {


// perform forward algorithm
double **arghmm_forward_alg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, bool prior_given, double *prior,
    bool internal, bool slow)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, trees);
    matrix_list.set_internal(internal);
    matrix_list.setup();
    matrix_list.begin();

    ArgHmmForwardTableOld forward(0, sequences.length());

    // setup prior
    if (prior_given) {
        LocalTree *tree = matrix_list.get_tree_spr()->tree;
        LineageCounts lineages(ntimes);
        States states;
        get_coal_states(tree, model.ntimes, states, internal);

        int start = matrix_list.get_block_start();
        forward.new_block(start, matrix_list.get_block_end(), states.size());
        double **fw = forward.get_table();
        for (unsigned int i=0; i<states.size(); i++)
            fw[0][i] = prior[i];
    }

    arghmm_forward_alg(trees, &model, &sequences, &matrix_list,
                       &forward, NULL, prior_given, internal, slow);

    // steal pointer
    double **fw = forward.detach_table();

    return fw;
}


// perform forward algorithm and sample threading path from posterior
intstate *arghmm_sample_posterior(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, intstate *path=NULL)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    Sequences sequences(seqs, nseqs, seqlen);

    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, &trees);
    matrix_list.setup();

    // compute forward table
    ArgHmmForwardTable forward(0, seqlen);
    arghmm_forward_alg(&trees, &model, &sequences, &matrix_list, &forward);

    // traceback
    int *ipath = new int [seqlen];
    stochastic_traceback(&trees, &model, &matrix_list,
                         forward.get_table(), ipath);

    // convert path
    if (path == NULL)
        path = new intstate [seqlen];

    States states;
    int end = trees.start_coord;
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); ++it) {
        int start = end;
        int end = start + it->blocklen;
        get_coal_states(it->tree, ntimes, states, false);

        for (int i=start; i<end; i++) {
            int istate = ipath[i];
            path[i][0] = states[istate].node;
            path[i][1] = states[istate].time;
        }
    }


    // clean up
    delete [] ipath;

    return path;
}


// sample the thread of an internal branch
void arghmm_sample_arg_thread_internal(LocalTrees *trees,
    double *times, int ntimes, double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int *thread_path)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    const bool internal = true;

    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    thread_path = &thread_path[-trees->start_coord];

    // forward algorithm
    ArgHmmMatrixIter matrix_iter(&model, &sequences, trees);
    matrix_iter.set_internal(internal);
    arghmm_forward_alg(trees, &model, &sequences, &matrix_iter, &forward,
                       NULL, false, internal);

    // traceback
    double **fw = forward.get_table();
    ArgHmmMatrixIter matrix_iter2(&model, NULL, trees);
    matrix_iter2.set_internal(internal);
    stochastic_traceback(trees, &model, &matrix_iter2, fw, thread_path,
                         false, internal);
}


// add one chromosome to an ARG
LocalTrees *arghmm_sample_thread(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    int new_chrom = nseqs -  1;

    sample_arg_thread(&model, &sequences, trees, new_chrom);

    return trees;
}


void delete_path(int *path)
{
    delete [] path;
}


void delete_double_matrix(double **mat, int nrows)
{
    delete_matrix<double>(mat, nrows);
}

void delete_forward_matrix(double **mat, int nrows)
{
    for (int i=0; i<nrows; i++)
        delete [] mat[i];
    delete [] mat;
}



} // extern "C"
*/
} // namespace argweaver
