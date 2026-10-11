//=============================================================================
// Joint re-threading of a haplotype and its read-backed phase
//
// When a leaf (haplotype) of a diploid individual is re-threaded, the
// threading HMM's state (where the new branch joins the tree) is paired with
// the 2-state read-backed phase state of phase_hmm.h ("1-1" or "1-2").  The
// branch moves as in ordinary threading, and independently the phase state
// moves as in the phase HMM: at each heterozygous site of the individual it
// is kept with 1-s and switched with s within a phase block, or drawn afresh
// 50/50 at the start of a block.  Emissions are the likelihood of the data
// with the individual's two alleles in each orientation.  A forward pass and
// stochastic traceback then sample the new thread and the phase path
// together, as an exact Gibbs step.
//
// With point errors (see phase_hmm.h), the emissions of each phase state at
// sites in a block are averaged over the alleles at the site being the
// other way round (probability point_prob), and after the traceback each
// such site is sampled as a point error or not, given its phase state and
// thread state.

#ifndef ARGWEAVER_PHASE_THREAD_H
#define ARGWEAVER_PHASE_THREAD_H

#include "local_tree.h"
#include "model.h"
#include "phase_hmm.h"
#include "sample_thread.h"
#include "sequences.h"

namespace argweaver {

// Forward algorithm over the doubled (thread state, phase state) space, for
// internal re-threading of the leaf in phase_pr, whose individual's
// haplotypes are both in the ARG.  mix[i] and block_site[i] give, for
// column start_coord+i, the change of phase state and the index of the site
// in 'blocks' (or -1).  'forward' receives phase state "1-1" and 'forward2'
// phase state "1-2".  With point errors (point_prob > 0), orient11 and
// orient12, if given, receive the emissions of each orientation at each
// site in a block, for each thread state, for sampling point errors after
// the traceback.  Exposed for testing.
void forward_phase(const LocalTrees *trees, const ArgModel *model,
                   ArgHmmMatrixIter *matrix_iter, PhaseProbs *phase_pr,
                   const Sequences *sequences, const PhaseBlocks &blocks,
                   const vector<int> &block_site,
                   const vector<char> &mix, double switch_prob,
                   double point_prob, ArgHmmForwardTable *forward,
                   ArgHmmForwardTable *forward2,
                   vector<vector<double> > *orient11=NULL,
                   vector<vector<double> > *orient12=NULL);

// Re-thread leaf 'node' (whose haplotype belongs to the individual of
// 'blocks') and resample that individual's phase jointly.  Alleles in
// 'sequences' are swapped to match the sampled phase.  Returns the number of
// sites swapped.
int resample_arg_leaf_phase(const ArgModel *model, Sequences *sequences,
                            LocalTrees *trees, int node,
                            const PhaseBlocks &blocks, double switch_prob,
                            double point_prob);

} // namespace argweaver

#endif // ARGWEAVER_PHASE_THREAD_H
