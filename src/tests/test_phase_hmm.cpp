// Tests for the read-backed phase HMM (src/argweaver/phase_hmm.cpp).
//
// usage: test_phase_hmm [<sites file> <smc file>]
//
// Without arguments only the HMM itself is tested.  With a sites file and an
// ARG (smc file, uncompressed, same coordinates) the emissions are tested too.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>

#include "argweaver/common.h"
#include "argweaver/local_tree.h"
#include "argweaver/model.h"
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


// The 5-site worked example from phase-hmm-explained.md
static void test_worked_example()
{
    PhaseBlocks blocks;
    int block[] = {0, 0, 0, 1, 1};
    double e11[] = {.9, .8, .3, .2, .4};
    double e12[] = {.1, .2, .7, .8, .6};
    double expect[] = {.9, .9706, .9139, .2, .1475};
    vector<double> log11, log12, f11;
    for (int k=0; k<5; k++) {
        blocks.sites.push_back(k);
        blocks.block.push_back(block[k]);
        blocks.allele1.push_back('A');
        log11.push_back(log(e11[k]));
        log12.push_back(log(e12[k]));
    }

    phase_hmm_forward(blocks, .01, log11, log12, f11);
    bool ok = true;
    for (int k=0; k<5; k++) {
        printf("    site %d: f(1-1) = %.4f (expected %.4f)\n",
               k+1, f11[k], expect[k]);
        if (fabs(f11[k] - expect[k]) > 1e-4) ok = false;
    }
    check(ok, "forward pass matches the worked example");
}


// Exact posterior probability of every path, by listing them all
static void brute_force_posterior(const PhaseBlocks &blocks, double s,
                                  const vector<double> &log11,
                                  const vector<double> &log12,
                                  vector<double> &post)
{
    int n = blocks.size();
    int npaths = 1 << n;
    post.assign(npaths, 0.0);
    double total = 0.0;
    for (int p=0; p<npaths; p++) {
        double prob = 1.0;
        for (int k=0; k<n; k++) {
            int state = (p >> k) & 1;
            if (blocks.same_block(k)) {
                int prev = (p >> (k-1)) & 1;
                prob *= (state == prev) ? 1.0 - s : s;
            } else {
                prob *= 0.5;
            }
            prob *= exp(state == 0 ? log11[k] : log12[k]);
        }
        post[p] = prob;
        total += prob;
    }
    for (int p=0; p<npaths; p++)
        post[p] /= total;
}


// Sampled paths should occur at their exact posterior frequencies
static void test_against_brute_force()
{
    const int ncases = 20;
    const int nsamples = 1000000;
    double worst = 0.0;

    for (int c=0; c<ncases; c++) {
        // random small case: 2 to 7 sites, some unphased, random blocks
        int n = 2 + irand(6);
        double s = frand(0.001, 0.3);
        PhaseBlocks blocks;
        vector<double> log11, log12;
        int b = 0;
        for (int k=0; k<n; k++) {
            if (frand() < 0.3) b++;
            blocks.sites.push_back(k);
            blocks.block.push_back(frand() < 0.15 ? -1 : b);
            blocks.allele1.push_back('A');
            log11.push_back(log(frand(0.01, 1.0)));
            log12.push_back(log(frand(0.01, 1.0)));
        }

        vector<double> post;
        brute_force_posterior(blocks, s, log11, log12, post);

        vector<double> freq(post.size(), 0.0);
        vector<int> path;
        for (int i=0; i<nsamples; i++) {
            sample_phase_path(blocks, s, log11, log12, path);
            int p = 0;
            for (int k=0; k<n; k++)
                p |= path[k] << k;
            freq[p] += 1.0 / nsamples;
        }

        // total variation distance between sampled and exact
        double tv = 0.0;
        for (unsigned int p=0; p<post.size(); p++)
            tv += 0.5 * fabs(freq[p] - post[p]);
        if (tv > worst) worst = tv;
    }
    printf("    worst total variation distance over %d cases: %.4f\n",
           ncases, worst);
    check(worst < 0.01, "sampled paths match the exact posterior");
}


// Emissions should equal the change in the whole ARG likelihood when the
// pair's alleles are swapped at one site
static void test_emissions(const char *sites_file, const char *smc_file)
{
    ArgModel model(20, 200e3, 10000, 1.6e-8, 1.8e-8);

    Sites sites;
    if (!read_sites(sites_file, &sites)) {
        check(false, "read sites file");
        return;
    }
    Sequences sequences;
    make_sequences_from_sites(&sites, &sequences);

    LocalTrees trees;
    vector<string> seqnames;
    if (!read_local_trees(smc_file, model.times, model.ntimes, &trees,
                          seqnames) ||
        !trees.set_seqids(seqnames, sequences.names)) {
        check(false, "read smc file");
        return;
    }

    // first two sequences form the individual; read-backed haplotype 1 is
    // taken to be the current first haplotype
    PhaseBlocks blocks(0, 1);
    for (int i=0; i<sequences.length(); i++) {
        if (sequences.seqs[0][i] != sequences.seqs[1][i]) {
            blocks.sites.push_back(i);
            blocks.block.push_back(0);
            blocks.allele1.push_back(sequences.seqs[0][i]);
        }
    }
    printf("    %s and %s: %d heterozygous sites\n",
           sequences.names[0].c_str(), sequences.names[1].c_str(),
           blocks.size());

    vector<double> log11, log12;
    calc_phase_emissions(&model, &sequences, &trees, blocks, log11, log12);

    double full = calc_arg_likelihood(&model, &sequences, &trees, -1, -1);
    double worst = 0.0;
    int prefer11 = 0;
    for (int k=0; k<blocks.size(); k++) {
        int pos = blocks.sites[k];
        sequences.switch_alleles(pos, 0, 1);
        double full_swapped = calc_arg_likelihood(&model, &sequences,
                                                  &trees, -1, -1);
        sequences.switch_alleles(pos, 0, 1);
        double diff = fabs((full - full_swapped) - (log11[k] - log12[k]));
        if (diff > worst) worst = diff;
        if (log11[k] > log12[k]) prefer11++;
    }
    printf("    largest mismatch in log likelihood ratio: %g\n", worst);
    printf("    sites where the ARG prefers the input phase: %d of %d\n",
           prefer11, blocks.size());
    check(blocks.size() > 0 && worst < 1e-6,
          "emissions match whole-ARG likelihood differences");
}


// With point errors, sampled orientations should occur at their exact
// posterior frequencies: the prior of an orientation vector is the sum over
// all state paths of the path's prior times e for each site in a block
// whose orientation differs from its state, and 1 - e for each that agrees.
static void test_point_errors()
{
    const int ncases = 20;
    const int nsamples = 1000000;
    double worst = 0.0;

    for (int c=0; c<ncases; c++) {
        int n = 2 + irand(6);
        double s = frand(0.001, 0.3);
        double e = (c % 4 == 0) ? 0.5 : frand(0.001, 0.3);
        PhaseBlocks blocks;
        vector<double> log11, log12;
        int b = 0;
        for (int k=0; k<n; k++) {
            if (frand() < 0.3) b++;
            blocks.sites.push_back(k);
            blocks.block.push_back(frand() < 0.15 ? -1 : b);
            blocks.allele1.push_back('A');
            log11.push_back(log(frand(0.01, 1.0)));
            log12.push_back(log(frand(0.01, 1.0)));
        }

        // exact posterior of orientations
        const int nvec = 1 << n;
        vector<double> post(nvec, 0.0);
        double total = 0.0;
        for (int o=0; o<nvec; o++) {
            double prior = 0.0;
            for (int p=0; p<nvec; p++) {
                double prob = 1.0;
                for (int k=0; k<n; k++) {
                    int state = (p >> k) & 1, orient = (o >> k) & 1;
                    if (blocks.same_block(k))
                        prob *= (state == ((p >> (k-1)) & 1)) ? 1.0 - s : s;
                    else
                        prob *= 0.5;
                    if (blocks.block[k] >= 0)
                        prob *= (orient == state) ? 1.0 - e : e;
                    else if (orient != state)
                        prob = 0.0;
                }
                prior += prob;
            }
            double lik = 1.0;
            for (int k=0; k<n; k++)
                lik *= exp(((o >> k) & 1) == 0 ? log11[k] : log12[k]);
            post[o] = prior * lik;
            total += post[o];
        }

        vector<double> freq(nvec, 0.0);
        vector<int> orient;
        for (int i=0; i<nsamples; i++) {
            sample_phase_orientation(blocks, s, e, log11, log12, orient);
            int o = 0;
            for (int k=0; k<n; k++)
                o |= orient[k] << k;
            freq[o] += 1.0 / nsamples;
        }
        double tv = 0.0;
        for (int o=0; o<nvec; o++)
            tv += 0.5 * fabs(freq[o] - post[o] / total);
        if (tv > worst) worst = tv;
    }
    printf("    worst total variation distance over %d cases: %.4f\n",
           ncases, worst);
    check(worst < 0.01,
          "with point errors, sampled orientations match the exact"
          " posterior");
}


int main(int argc, char **argv)
{
    srand(1);
    test_worked_example();
    test_against_brute_force();
    test_point_errors();
    if (argc == 3)
        test_emissions(argv[1], argv[2]);
    printf("%s\n", nfail == 0 ? "all tests passed" : "SOME TESTS FAILED");
    return nfail == 0 ? 0 : 1;
}
