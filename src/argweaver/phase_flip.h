//=============================================================================
// Flipping the relative phase of two read-backed phase blocks
//
// Swapping an individual's two haplotypes to the right of a boundary between
// phase blocks flips the relative phase of the blocks on either side.  In
// the ARG, it is the same as swapping which leaf is which from that point
// on: to the right of a window around the boundary nothing changes but the
// labels, so the likelihood and prior there stay the same.  Only the join
// inside the window needs repair.
//
// The move (Metropolis-Hastings):
//   1. Pick at random which haplotype is threaded first (u), and which
//      second (v).
//   2. Swap the individual's alleles to the right of the boundary.
//   3. Remove both haplotypes from the ARG.  Re-thread u through the window,
//      pinned on the left to its old state and on the right to v's old
//      state, with v absent.  Then re-thread v, pinned on the left to its
//      old state and on the right to u's old state.  Outside the window
//      the threads are kept (with the labels swapped on the right).
//   4. Accept with the ratio of posteriors times the ratio of the
//      probabilities of proposing the reverse and the forward move.  The
//      reverse move is the same flip from the new state; it re-threads the
//      old threads through the same window.

#ifndef ARGWEAVER_PHASE_FLIP_H
#define ARGWEAVER_PHASE_FLIP_H

#include "local_tree.h"
#include "model.h"
#include "phase_hmm.h"
#include "sequences.h"

namespace argweaver {

// terms of one flip, for testing
struct FlipMoveInfo
{
    int lo, hi;              // pinned columns; threads re-sampled in between
    int first;               // haplotype threaded first
    double dprior, dlik;     // change in log prior and log likelihood
    double lq_fwd, lq_rev;   // log probability of the forward and reverse
                             // proposals
    double log_ratio;
    bool supported;          // false if the old threads could not be
                             // proposed by the move (then rejected)
};

// The window for boundary k (between heterozygous sites k-1 and k of the
// individual): the threads are pinned at columns lo and hi and re-sampled
// in between.  Sites more than 'margin' positions outside the gap between
// the two sites are left out.
void flip_window(const LocalTrees *trees, const PhaseBlocks &blocks, int k,
                 int margin, int *lo, int *hi);

// Propose flipping the phase at boundary k (blocks.same_block(k) must be
// false).  Returns whether the flip was accepted; if not, the ARG and
// sequences are unchanged.  'first', if given, is the haplotype threaded
// first (for testing); otherwise it is chosen at random.
bool flip_phase_boundary(const ArgModel *model, Sequences *sequences,
                         LocalTrees *trees, const PhaseBlocks &blocks, int k,
                         int margin, FlipMoveInfo *info=NULL,
                         int first=-1);

// Propose a flip at every boundary between blocks, in turn, accepting or
// rejecting each one separately.  Returns the number accepted, and gives
// the number tried.
int flip_phase_boundaries(const ArgModel *model, Sequences *sequences,
                          LocalTrees *trees, const PhaseBlocks &blocks,
                          int margin, int *ntried=NULL);

// For testing: remove leaf 'seqid' from 'full', read its thread off 'full'
// and add it back.  'rebuilt' should then equal 'full'.  Returns false if
// the thread cannot be read off.
bool rebuild_leaf_thread(const ArgModel *model, const Sequences *sequences,
                         const LocalTrees *full, int seqid,
                         LocalTrees *rebuilt);

// For testing: accept every flip that the move can propose
extern bool flip_force_accept;

} // namespace argweaver

#endif // ARGWEAVER_PHASE_FLIP_H
