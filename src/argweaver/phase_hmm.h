//=============================================================================
// Phase resampling that respects read-backed phase blocks (e.g. WhatsHap)
//
// For one diploid individual, given the current ARG, the orientation of its
// two haplotypes at each heterozygous site is modelled with a 2-state HMM:
//
//   state 0 ("1-1"): ARG haplotype 1 carries read-backed haplotype 1's allele
//   state 1 ("1-2"): ARG haplotype 1 carries read-backed haplotype 2's allele
//
// Within a phase block the state is kept with probability 1-s and switched
// with probability s (a switch error).  Between blocks, or at sites not in
// any block, the state is chosen afresh with probability 1/2 each.  The
// emission at each site is the likelihood of the local tree with the
// individual's two alleles in each orientation.
//
// Optionally, each site in a block can also be a point error: with
// probability e, independently of the other sites, its alleles are the
// other way round from the state of the block at that site, without the
// state itself switching.  The orientation of the alleles at the site is
// then the state, flipped if there is a point error.

#ifndef ARGWEAVER_PHASE_HMM_H
#define ARGWEAVER_PHASE_HMM_H

#include <vector>

#include "local_tree.h"
#include "model.h"
#include "sequences.h"

namespace argweaver {

using namespace std;


// Read-backed phase blocks for one diploid individual.  All vectors are
// indexed by k, the k-th heterozygous site of this individual.
class PhaseBlocks
{
public:
    PhaseBlocks(int hap1=-1, int hap2=-1) : hap1(hap1), hap2(hap2) {}

    int size() const { return sites.size(); }

    // true if sites k-1 and k are in the same phase block
    bool same_block(int k) const {
        return k > 0 && block[k] >= 0 && block[k] == block[k-1];
    }

    int hap1, hap2;        // sequence indices of the two haplotypes
    vector<int> sites;     // site positions (ARGweaver coordinates), sorted
    vector<int> block;     // phase block id of each site; -1 if unphased
    vector<char> allele1;  // allele on read-backed haplotype 1 at each site
};


// Forward pass of the phase HMM.  On return, f11[k] is the probability of
// state "1-1" at site k given the emissions at sites 0..k (the probability
// of "1-2" is 1 - f11[k]).  Emissions are given as natural logs.
void phase_hmm_forward(const PhaseBlocks &blocks, double switch_prob,
                       const vector<double> &log_emit11,
                       const vector<double> &log_emit12,
                       vector<double> &f11);

// Sample a state path (0 = "1-1", 1 = "1-2") from the posterior of the
// phase HMM, by a forward pass followed by stochastic traceback.
void sample_phase_path(const PhaseBlocks &blocks, double switch_prob,
                       const vector<double> &log_emit11,
                       const vector<double> &log_emit12,
                       vector<int> &path);

// Sample the orientation of the alleles at each site (0 = "1-1",
// 1 = "1-2") from the posterior of the phase HMM with point errors of
// probability 'point_prob' at sites in blocks: sample the state path with
// the emissions of each state averaged over a point error, then sample
// whether each site is a point error given its state.  With point_prob = 0
// this is sample_phase_path.
void sample_phase_orientation(const PhaseBlocks &blocks, double switch_prob,
                              double point_prob,
                              const vector<double> &log_emit11,
                              const vector<double> &log_emit12,
                              vector<int> &orient);

// Compute the log emissions of both orientations at each site in 'blocks',
// using the local trees of the current ARG.  Sequences and trees must use
// the same coordinates.
void calc_phase_emissions(const ArgModel *model, const Sequences *sequences,
                          const LocalTrees *trees, const PhaseBlocks &blocks,
                          vector<double> &log_emit11,
                          vector<double> &log_emit12);

// Resample the phase of one individual given the current ARG, and swap
// alleles in 'sequences' to match.  Returns the number of sites swapped.
// With point errors (point_prob > 0), this also moves single sites within
// a block.
int resample_phase(const ArgModel *model, Sequences *sequences,
                   const LocalTrees *trees, const PhaseBlocks &blocks,
                   double switch_prob, double point_prob=0.0);

// Read phase blocks from a PHASE file (may be gzipped).  The format follows
// the SITES format:
//
//   NAMES   <ind1>  <ind2> ...
//   REGION  <chrom> <start> <end>
//   <pos>   <block1>,<block2>,...
//
// with one block id (a non-negative integer, or '.' if unknown) per
// individual, separated by commas.  Individual <ind> has haplotypes named
// <ind>_1 and <ind>_2 in the sequences.  Positions are 1-based, as in
// SITES, and are mapped to ARGweaver's (compressed) coordinates using
// 'sites_mapping'.  The current phase of the sequences is taken as the
// read-backed phase.
bool read_phase_blocks(const char *filename, const Sequences *sequences,
                       const SitesMapping *sites_mapping,
                       vector<PhaseBlocks> &all_blocks);

} // namespace argweaver

#endif // ARGWEAVER_PHASE_HMM_H
