// Checks for the block-boundary flip move (src/argweaver/phase_flip.cpp).
//
// usage: test_flip_move <sites> <smc> <popsize> <maxtime> <phase file>
//                       [margin] [nboundaries]
//
// 1. Reading a thread off the ARG: removing a haplotype and adding back the
//    thread read off the ARG gives back the same ARG (same prior and
//    likelihood).
// 2. Each flip: the change in posterior computed over the window matches
//    the change over the whole ARG, and flipping back from the new state
//    (threading the same haplotype first) gives as its reverse proposal
//    probability the forward proposal probability of the first flip.
//    Flips are forced to be accepted for this check.
// 3. Proposal probabilities are normalised (an importance-sampling check
//    between the two orders of threading).
// 4. Acceptance, without forcing: log ratios and acceptance at each
//    boundary of the first individual.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>
#include <algorithm>

#include "argweaver/common.h"
#include "argweaver/compress.h"
#include "argweaver/local_tree.h"
#include "argweaver/model.h"
#include "argweaver/phase_flip.h"
#include "argweaver/phase_hmm.h"
#include "argweaver/sequences.h"
#include "argweaver/total_prob.h"

using namespace argweaver;
using namespace std;

static int nfail = 0;

static void check(bool ok, const char *what)
{
    printf("%s  %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) nfail++;
}


int main(int argc, char **argv)
{
    if (argc < 6) {
        fprintf(stderr, "usage: test_flip_move <sites> <smc> <popsize>"
                " <maxtime> <phase file> [margin] [nboundaries]\n");
        return 1;
    }
    srand(1);
    int margin = (argc > 6) ? atoi(argv[6]) : 1000;
    int nbound = (argc > 7) ? atoi(argv[7]) : 20;
    ArgModel model(20, atof(argv[4]), atof(argv[3]), 1e-8, 1.2e-8);
    model.smc_prime = false;  // as arg-sample by default

    Sites sites;
    Sequences seqs;
    LocalTrees trees;
    vector<string> seqnames;
    if (!read_sites(argv[1], &sites)) return 1;
    SitesMapping sites_mapping;
    find_compress_cols(&sites, 1, &sites_mapping);
    compress_sites(&sites, &sites_mapping);
    make_sequences_from_sites(&sites, &seqs);
    seqs.set_age();
    bool ok = read_local_trees(argv[2], model.times, model.ntimes, &trees,
                               seqnames);
    trees.seqids.resize(seqnames.size());
    if (!ok || !trees.set_seqids(seqnames, seqs.names)) {
        fprintf(stderr, "cannot read %s\n", argv[2]);
        return 1;
    }
    vector<PhaseBlocks> all_blocks;
    if (!read_phase_blocks(argv[5], &seqs, &sites_mapping, all_blocks)) {
        fprintf(stderr, "cannot read %s\n", argv[5]);
        return 1;
    }
    const PhaseBlocks &blocks = all_blocks[0];
    vector<int> bounds;
    for (int k=1; k<blocks.size(); k++)
        if (!blocks.same_block(k)) bounds.push_back(k);
    printf("individual %s/%s: %d heterozygous sites, %d boundaries\n",
           seqs.names[blocks.hap1].c_str(), seqs.names[blocks.hap2].c_str(),
           blocks.size(), (int) bounds.size());

    // 1. reading threads off the ARG
    const double prior0 = calc_arg_prior(&model, &trees);
    const double lik0 = calc_arg_likelihood(&model, &seqs, &trees);
    double worst = 0.0;
    int nread = 0;
    for (int h=0; h<seqs.get_num_seqs(); h++) {
        LocalTrees rebuilt;
        if (!rebuild_leaf_thread(&model, &seqs, &trees, h, &rebuilt))
            continue;
        nread++;
        worst = max(worst, fabs(calc_arg_prior(&model, &rebuilt) - prior0));
        worst = max(worst, fabs(calc_arg_likelihood(&model, &seqs, &rebuilt)
                                - lik0));
    }
    printf("    threads read off: %d of %d; largest change in prior or"
           " likelihood after re-adding: %g\n", nread, seqs.get_num_seqs(),
           worst);
    check(nread == seqs.get_num_seqs() && worst < 1e-6,
          "a thread read off the ARG adds back to the same ARG");

    // 2. forced flips and their reverses
    flip_force_accept = true;
    double worst_post = 0.0, worst_rev = 0.0;
    int nbad = 0;
    int nflip = 0;
    for (unsigned int b=0; b<bounds.size() && (int) b<nbound; b++) {
        const int k = bounds[b];
        const double p0 = calc_arg_prior(&model, &trees);
        const double l0 = calc_arg_likelihood(&model, &seqs, &trees);
        FlipMoveInfo info, info2;
        if (!flip_phase_boundary(&model, &seqs, &trees, blocks, k, margin,
                                 &info, blocks.hap1)) {
            printf("    boundary %d: not proposed (supported %d)\n", k,
                   (int) info.supported);
            continue;
        }
        if (!assert_trees(&trees, NULL)) nbad++;
        const double p1 = calc_arg_prior(&model, &trees);
        const double l1 = calc_arg_likelihood(&model, &seqs, &trees);
        const double dpost = fabs((p1 + l1 - p0 - l0) -
                                  (info.dprior + info.dlik));
        // flip back from the new state
        bool back = flip_phase_boundary(&model, &seqs, &trees, blocks, k,
                                        margin, &info2, blocks.hap1);
        const double drev = back ? fabs(info2.lq_rev - info.lq_fwd) : -1;
        printf("    boundary %3d window %5d-%5d (%4d): dprior %8.2f"
               " dlik %8.2f lq_fwd %8.2f lq_rev %8.2f ratio %8.2f;"
               " whole-ARG diff %.2g, reverse diff %.2g\n",
               k, info.lo, info.hi, info.hi - info.lo, info.dprior,
               info.dlik, info.lq_fwd, info.lq_rev, info.log_ratio, dpost,
               drev);
        nflip++;
        worst_post = max(worst_post, dpost);
        if (back) worst_rev = max(worst_rev, drev);
        else printf("    flip back not proposed\n");
    }
    // again, with the window chosen to end just before a recombination
    // (where the ARG's blocks are cut differently when the window is put
    // back)
    int nstub = 0;
    for (unsigned int b=0; b<bounds.size() && nstub < 5; b++) {
        const int k = bounds[b];
        int bstart = -1, end = trees.start_coord;
        for (LocalTrees::iterator it=trees.begin(); it!=trees.end(); ++it) {
            if (end > blocks.sites[k] + 1) { bstart = end; break; }
            end += it->blocklen;
        }
        if (bstart < 0) continue;
        const int m = bstart - 1 - blocks.sites[k];
        const double p0 = calc_arg_prior(&model, &trees);
        const double l0 = calc_arg_likelihood(&model, &seqs, &trees);
        FlipMoveInfo info, info2;
        if (!flip_phase_boundary(&model, &seqs, &trees, blocks, k, m,
                                 &info, blocks.hap1))
            continue;
        if (!assert_trees(&trees, NULL)) nbad++;
        const double p1 = calc_arg_prior(&model, &trees);
        const double l1 = calc_arg_likelihood(&model, &seqs, &trees);
        const double dpost = fabs((p1 + l1 - p0 - l0) -
                                  (info.dprior + info.dlik));
        bool back = flip_phase_boundary(&model, &seqs, &trees, blocks, k, m,
                                        &info2, blocks.hap1);
        if (!assert_trees(&trees, NULL)) nbad++;
        const double drev = back ? fabs(info2.lq_rev - info.lq_fwd) : 0;
        printf("    boundary %3d, window %d-%d ending before a"
               " recombination at %d: whole-ARG diff %.2g, reverse diff"
               " %.2g\n", k, info.lo, info.hi, bstart, dpost, drev);
        worst_post = max(worst_post, dpost);
        worst_rev = max(worst_rev, drev);
        nstub++;
    }
    flip_force_accept = false;
    check(nstub > 0, "windows ending before a recombination were tried");
    check(nbad == 0, "the ARG is valid after each flip");
    check(nflip > 0 && worst_post < 1e-6,
          "change in posterior over the window = change over the whole ARG");
    check(nflip > 0 && worst_rev < 1e-6,
          "reverse proposal probability matches the forward one");

    // 3. proposal probabilities are normalised: a proposal G' drawn with
    // hap1 threaded first, weighted by the probability of proposing it
    // with hap2 threaded first over the probability with hap1 first, has
    // mean weight 1 (both orders propose from the same set).  The second
    // probability is the reverse proposal probability of flipping back
    // from G' with hap2 first.
    for (int pass=0; pass<2; pass++) {
        const int reps = 400;
        const int margin3 = (pass == 0) ? 0 : 200;
        int kbest = -1, wbest = 1 << 30;
        for (unsigned int b=0; b<bounds.size(); b++) {
            int lo, hi;
            flip_window(&trees, blocks, bounds[b], margin3, &lo, &hi);
            if (hi - lo < wbest) { wbest = hi - lo; kbest = bounds[b]; }
        }
        LocalTrees start;
        start.copy(trees);
        vector<char> s1(seqs.seqs[blocks.hap1],
                        seqs.seqs[blocks.hap1] + seqs.length());
        vector<char> s2(seqs.seqs[blocks.hap2],
                        seqs.seqs[blocks.hap2] + seqs.length());
        flip_force_accept = true;
        double sw = 0.0, sw2 = 0.0;
        int n = 0;
        for (int r=0; r<reps; r++) {
            trees.copy(start);
            memcpy(seqs.seqs[blocks.hap1], &s1[0], seqs.length());
            memcpy(seqs.seqs[blocks.hap2], &s2[0], seqs.length());
            FlipMoveInfo info, info2;
            if (!flip_phase_boundary(&model, &seqs, &trees, blocks, kbest,
                                     margin3, &info, blocks.hap1))
                continue;
            if (!flip_phase_boundary(&model, &seqs, &trees, blocks, kbest,
                                     margin3, &info2, blocks.hap2))
                continue;
            double w = exp(info2.lq_rev - info.lq_fwd);
            sw += w;
            sw2 += w * w;
            n++;
        }
        flip_force_accept = false;
        trees.copy(start);
        memcpy(seqs.seqs[blocks.hap1], &s1[0], seqs.length());
        memcpy(seqs.seqs[blocks.hap2], &s2[0], seqs.length());
        double mean = sw / n;
        double se = sqrt((sw2 / n - mean * mean) / n);
        printf("    boundary %d, window %d: mean weight %.3f (se %.3f)"
               " over %d proposals\n", kbest, wbest, mean, se, n);
        check(n > 0 && fabs(mean - 1.0) < 4 * se + 0.02,
              "proposal probabilities are normalised");
    }

    // 4. unforced
    int naccept = 0, ntry = 0;
    vector<double> ratios;
    for (unsigned int b=0; b<bounds.size(); b++) {
        FlipMoveInfo info;
        bool acc = flip_phase_boundary(&model, &seqs, &trees, blocks,
                                       bounds[b], margin, &info);
        if (!info.supported) continue;
        ntry++;
        naccept += acc;
        ratios.push_back(info.log_ratio);
    }
    sort(ratios.begin(), ratios.end());
    printf("    unforced: accepted %d of %d; log ratio median %.2f,"
           " quartiles %.2f %.2f\n", naccept, ntry,
           ratios.size() ? ratios[ratios.size()/2] : 0.0,
           ratios.size() ? ratios[ratios.size()/4] : 0.0,
           ratios.size() ? ratios[3*ratios.size()/4] : 0.0);

    printf("%s\n", nfail == 0 ? "all tests passed" : "SOME TESTS FAILED");
    return nfail == 0 ? 0 : 1;
}
