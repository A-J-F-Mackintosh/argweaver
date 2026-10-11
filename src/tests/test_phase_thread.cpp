// Tests for joint re-threading of a haplotype and its read-backed phase
// (src/argweaver/phase_thread.cpp).
//
// usage: test_phase_thread <sites file> <smc file> <popsize> <maxtime> [nreps]
//
// The first two sequences form the test individual.  Two checks:
//
// 1. Forward pass.  If every heterozygous site is its own phase block (or
//    s = 1/2), the phase state carries no information between sites, and the
//    model is exactly ARGweaver's own unphased threading.  The two phase
//    layers of our forward table must then add up to ARGweaver's unphased
//    forward table.
//
// 2. Sampling.  In the same setting, re-threading with the joint move and
//    with ARGweaver's unphased leaf move must give the same distribution of
//    phase swaps at each site and of the re-threaded branch.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "argweaver/common.h"
#include "argweaver/local_tree.h"
#include "argweaver/matrices.h"
#include "argweaver/model.h"
#include "argweaver/phase_hmm.h"
#include "argweaver/phase_thread.h"
#include "argweaver/sample_arg.h"
#include "argweaver/sample_thread.h"
#include "argweaver/sequences.h"
#include "argweaver/thread.h"

using namespace argweaver;
using namespace std;

static int nfail = 0;

static void check(bool ok, const char *what)
{
    printf("%s  %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) nfail++;
}


static int find_leaf(const LocalTrees *trees, int seqid)
{
    for (int i=0; i<trees->get_num_leaves(); i++)
        if (trees->seqids[i] == seqid) return i;
    return -1;
}


// heterozygous sites of the individual, each its own block or all in one
static PhaseBlocks make_blocks(const Sequences *seqs, bool one_block)
{
    PhaseBlocks blocks(0, 1);
    for (int i=0; i<seqs->length(); i++) {
        if (seqs->seqs[0][i] != seqs->seqs[1][i]) {
            blocks.block.push_back(one_block ? 0 : blocks.size());
            blocks.sites.push_back(i);
            blocks.allele1.push_back(seqs->seqs[0][i]);
        }
    }
    return blocks;
}


// largest difference between ARGweaver's unphased forward table and the sum
// of our two phase layers, after removing the leaf of sequence 0
static double forward_difference(ArgModel *model, Sequences *seqs,
                                 const LocalTrees *trees0,
                                 const PhaseBlocks &blocks, double s,
                                 double e=0.0)
{
    LocalTrees trees;
    trees.copy(*trees0);
    int node = find_leaf(&trees, 0);
    int *removal_path = new int [trees.get_num_trees()];
    sample_arg_removal_leaf_path(&trees, node, removal_path);
    remove_arg_thread_path(&trees, removal_path,
                           model->get_removed_root_time(), model->pop_tree);
    delete [] removal_path;
    const int start = trees.start_coord, len = trees.length();

    model->unphased = true;

    // ARGweaver's unphased forward pass
    PhaseProbs pp1(0, node, seqs, &trees, model);
    ArgHmmMatrixIter iter1(model, seqs, &trees);
    iter1.set_internal(true, 0);
    ArgHmmForwardTable fw1(start, len);
    arghmm_forward_alg(&trees, model, seqs, &iter1, &fw1, &pp1, false, true);

    // ours
    vector<char> mix(len, PHASE_MIX_NONE);
    vector<int> block_site(len, -1);
    for (int k=0; k<blocks.size(); k++) {
        mix[blocks.sites[k] - start] =
            blocks.same_block(k) ? PHASE_MIX_SWITCH : PHASE_MIX_REDRAW;
        block_site[blocks.sites[k] - start] = k;
    }
    PhaseProbs pp2(0, node, seqs, &trees, model);
    ArgHmmMatrixIter iter2(model, seqs, &trees);
    iter2.set_internal(true, 0);
    ArgHmmForwardTable fwa(start, len), fwb(start, len);
    forward_phase(&trees, model, &iter2, &pp2, seqs, blocks, block_site,
                  mix, s, e, &fwa, &fwb);

    model->unphased = false;

    // compare column by column
    double **t1 = fw1.get_table(), **ta = fwa.get_table(),
        **tb = fwb.get_table();
    double worst = 0.0;
    int i = start;
    for (LocalTrees::const_iterator it=trees.begin(); it!=trees.end(); ++it) {
        States states;
        iter1.states_model.get_coal_states(it->tree, states);
        int nstates = max((int) states.size(), 1);
        for (int j=0; j<it->blocklen; j++, i++)
            for (int k=0; k<nstates; k++)
                worst = max(worst, fabs(t1[i][k] - (ta[i][k] + tb[i][k])));
    }
    return worst;
}


// age (time index) of the parent of leaf 'node' at position pos
static int parent_age(const LocalTrees *trees, int node, int pos)
{
    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        end += it->blocklen;
        if (pos < end) {
            const LocalTree *tree = it->tree;
            return tree->nodes[tree->nodes[node].parent].age;
        }
    }
    return -1;
}


// re-thread the leaf of sequence 0 nreps times from the same starting
// point; record how often each site is swapped and the mean parent age
static void sample_stats(ArgModel *model, Sequences *seqs,
                         const LocalTrees *trees0, const PhaseBlocks &blocks,
                         bool joint, int nreps, const vector<int> &positions,
                         vector<double> &swap_freq, vector<double> &mean_age,
                         double s=0.5, double e=0.0)
{
    const int len = seqs->length();
    vector<char> orig0(seqs->seqs[0], seqs->seqs[0] + len);
    vector<char> orig1(seqs->seqs[1], seqs->seqs[1] + len);
    swap_freq.assign(blocks.size(), 0.0);
    mean_age.assign(positions.size(), 0.0);

    for (int r=0; r<nreps; r++) {
        LocalTrees trees;
        trees.copy(*trees0);
        memcpy(seqs->seqs[0], &orig0[0], len);
        memcpy(seqs->seqs[1], &orig1[0], len);
        int node = find_leaf(&trees, 0);
        if (joint) {
            resample_arg_leaf_phase(model, seqs, &trees, node, blocks, s,
                                    e);
        } else {
            model->unphased = true;
            resample_arg_leaf(model, seqs, &trees, node);
            model->unphased = false;
        }
        for (int k=0; k<blocks.size(); k++)
            if (seqs->seqs[0][blocks.sites[k]] != orig0[blocks.sites[k]])
                swap_freq[k] += 1.0 / nreps;
        for (unsigned int p=0; p<positions.size(); p++)
            mean_age[p] += parent_age(&trees, node, positions[p]) /
                double(nreps);
    }
    memcpy(seqs->seqs[0], &orig0[0], len);
    memcpy(seqs->seqs[1], &orig1[0], len);
}


int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr, "usage: test_phase_thread <sites> <smc> <popsize>"
                " <maxtime> [nreps]\n");
        return 1;
    }
    srand(1);
    int nreps = (argc > 5) ? atoi(argv[5]) : 500;
    ArgModel model(20, atof(argv[4]), atof(argv[3]), 1e-8, 1.2e-8);

    Sites sites;
    Sequences seqs;
    LocalTrees trees;
    vector<string> seqnames;
    if (!read_sites(argv[1], &sites)) return 1;
    make_sequences_from_sites(&sites, &seqs);
    seqs.set_age();
    if (!read_local_trees(argv[2], model.times, model.ntimes, &trees,
                          seqnames) ||
        !trees.set_seqids(seqnames, seqs.names)) {
        fprintf(stderr, "cannot read %s\n", argv[2]);
        return 1;
    }
    seqs.pairs.assign(seqs.get_num_seqs(), -1);
    seqs.pairs[0] = 1;
    seqs.pairs[1] = 0;

    PhaseBlocks own = make_blocks(&seqs, false);
    PhaseBlocks one = make_blocks(&seqs, true);
    printf("individual %s/%s: %d heterozygous sites\n",
           seqs.names[0].c_str(), seqs.names[1].c_str(), own.size());

    // 1. forward pass
    double d1 = forward_difference(&model, &seqs, &trees, own, 0.01);
    double d2 = forward_difference(&model, &seqs, &trees, one, 0.5);
    double d3 = forward_difference(&model, &seqs, &trees, one, 0.01, 0.5);
    printf("    every site its own block: max difference %g\n", d1);
    printf("    one block with s = 1/2:   max difference %g\n", d2);
    printf("    one block with s = 0.01, point errors e = 1/2:"
           " max difference %g\n", d3);
    check(d1 < 1e-9 && d2 < 1e-9 && d3 < 1e-9,
          "forward pass matches ARGweaver's unphased forward pass");

    // 2. sampling
    vector<int> positions;
    for (int p=1; p<=4; p++)
        positions.push_back(trees.start_coord + p * trees.length() / 5);
    vector<double> swap_a, age_a;
    sample_stats(&model, &seqs, &trees, own, false, nreps, positions,
                 swap_a, age_a);
    const char *labels[3] = {"joint, every site its own block",
                             "joint, one block with s = 1/2",
                             "joint, one block with s = 0.01, e = 1/2"};
    for (int m=0; m<3; m++) {
        vector<double> swap_b, age_b;
        sample_stats(&model, &seqs, &trees, (m == 0) ? own : one, true,
                     nreps, positions, swap_b, age_b,
                     (m == 2) ? 0.01 : 0.5, (m == 2) ? 0.5 : 0.0);

        // swap frequencies: compare with a two-sample z-score
        double worst_z = 0.0;
        int nflagged = 0;
        for (int k=0; k<own.size(); k++) {
            double p = 0.5 * (swap_a[k] + swap_b[k]);
            double se = sqrt(2.0 * p * (1.0 - p) / nreps);
            if (se > 0) {
                double z = fabs(swap_a[k] - swap_b[k]) / se;
                worst_z = max(worst_z, z);
                if (z > 3.0) nflagged++;
            }
        }
        printf("    %s: phase swap frequency at %d sites, %d repeats:"
               " worst z = %.2f, sites with z > 3: %d\n", labels[m],
               own.size(), nreps, worst_z, nflagged);
        for (int k=0; k<min(5, own.size()); k++)
            printf("      site %d: ARGweaver %.3f, joint %.3f\n",
                   own.sites[k], swap_a[k], swap_b[k]);
        for (unsigned int p=0; p<positions.size(); p++)
            printf("      mean parent age (time index) at %d: ARGweaver"
                   " %.3f, joint %.3f\n", positions[p], age_a[p], age_b[p]);
        // with this many sites, more than one z > 3 would be suspicious
        check(nflagged <= 1, "phase swaps match ARGweaver's unphased"
              " leaf move");
    }

    printf("%s\n", nfail == 0 ? "all tests passed" : "SOME TESTS FAILED");
    return nfail == 0 ? 0 : 1;
}
