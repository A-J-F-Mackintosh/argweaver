// c++ includes
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// arghmm includes
#include "common.h"
#include "compress.h"
#include "emit.h"
#include "logging.h"
#include "parsing.h"
#include "phase_hmm.h"

namespace argweaver {

using namespace std;


void phase_hmm_forward(const PhaseBlocks &blocks, double switch_prob,
                       const vector<double> &log_emit11,
                       const vector<double> &log_emit12,
                       vector<double> &f11)
{
    const int nsites = blocks.size();
    assert((int) log_emit11.size() == nsites);
    assert((int) log_emit12.size() == nsites);
    f11.resize(nsites);

    for (int k=0; k<nsites; k++) {
        // probability of "1-1" at site k before seeing its emission
        double prior11 = 0.5;
        if (blocks.same_block(k))
            prior11 = f11[k-1] * (1.0 - switch_prob) +
                (1.0 - f11[k-1]) * switch_prob;

        // weigh by the emissions; only their ratio matters, so subtract
        // the larger log emission to avoid underflow
        double top = max(log_emit11[k], log_emit12[k]);
        double p11 = prior11 * exp(log_emit11[k] - top);
        double p12 = (1.0 - prior11) * exp(log_emit12[k] - top);
        f11[k] = p11 / (p11 + p12);
    }
}


void sample_phase_path(const PhaseBlocks &blocks, double switch_prob,
                       const vector<double> &log_emit11,
                       const vector<double> &log_emit12,
                       vector<int> &path)
{
    const int nsites = blocks.size();
    path.resize(nsites);
    if (nsites == 0)
        return;

    vector<double> f11;
    phase_hmm_forward(blocks, switch_prob, log_emit11, log_emit12, f11);

    // last site: draw straight from the forward probabilities
    int k = nsites - 1;
    path[k] = (frand() < f11[k]) ? 0 : 1;

    // walk left, conditioning on the state just drawn to the right
    for (k=nsites-2; k>=0; k--) {
        double p11 = f11[k];
        double p12 = 1.0 - f11[k];
        if (blocks.same_block(k+1)) {
            int next = path[k+1];
            p11 *= (next == 0) ? 1.0 - switch_prob : switch_prob;
            p12 *= (next == 1) ? 1.0 - switch_prob : switch_prob;
        }
        path[k] = (frand() * (p11 + p12) < p11) ? 0 : 1;
    }
}


// log(exp(a) + exp(b))
static inline double log_add(double a, double b)
{
    double top = max(a, b);
    if (top == -INFINITY)
        return top;
    return top + log(exp(a - top) + exp(b - top));
}


void sample_phase_orientation(const PhaseBlocks &blocks, double switch_prob,
                              double point_prob,
                              const vector<double> &log_emit11,
                              const vector<double> &log_emit12,
                              vector<int> &orient)
{
    const int nsites = blocks.size();
    assert(point_prob >= 0.0 && point_prob <= 0.5);
    if (point_prob == 0.0) {
        sample_phase_path(blocks, switch_prob, log_emit11, log_emit12,
                          orient);
        return;
    }

    // emissions of each state, averaged over whether the site is a point
    // error (only sites in a block can be; elsewhere the state is drawn
    // afresh anyway)
    const double lkeep = log(1.0 - point_prob), lerr = log(point_prob);
    vector<double> mix11(log_emit11), mix12(log_emit12);
    for (int k=0; k<nsites; k++) {
        if (blocks.block[k] < 0)
            continue;
        mix11[k] = log_add(lkeep + log_emit11[k], lerr + log_emit12[k]);
        mix12[k] = log_add(lkeep + log_emit12[k], lerr + log_emit11[k]);
    }
    vector<int> path;
    sample_phase_path(blocks, switch_prob, mix11, mix12, path);

    // given the state, is the site a point error?
    orient = path;
    for (int k=0; k<nsites; k++) {
        if (blocks.block[k] < 0)
            continue;
        double own = (path[k] == 0) ? log_emit11[k] : log_emit12[k];
        double other = (path[k] == 0) ? log_emit12[k] : log_emit11[k];
        double p_err = 1.0 / (1.0 + exp(lkeep + own - lerr - other));
        if (frand() < p_err)
            orient[k] = 1 - path[k];
    }
}


void calc_phase_emissions(const ArgModel *model, const Sequences *sequences,
                          const LocalTrees *trees, const PhaseBlocks &blocks,
                          vector<double> &log_emit11,
                          vector<double> &log_emit12)
{
    const int nsites = blocks.size();
    const int nleaves = trees->get_num_leaves();
    log_emit11.resize(nsites);
    log_emit12.resize(nsites);

    // genotype probabilities are not supported yet
    const vector<vector<BaseProbs> > no_base_probs;
    assert(sequences->base_probs.size() == 0);

    // find the leaves holding the two haplotypes
    int leaf1 = -1, leaf2 = -1;
    for (int j=0; j<nleaves; j++) {
        if (trees->seqids[j] == blocks.hap1) leaf1 = j;
        if (trees->seqids[j] == blocks.hap2) leaf2 = j;
    }
    assert(leaf1 != -1 && leaf2 != -1);

    // sequences in leaf order, as they are now and with the pair swapped
    const char *seqs[nleaves];
    const char *swapped[nleaves];
    for (int j=0; j<nleaves; j++)
        seqs[j] = swapped[j] = sequences->seqs[trees->seqids[j]];
    swapped[leaf1] = seqs[leaf2];
    swapped[leaf2] = seqs[leaf1];

    // walk the local trees and the sites together
    int k = 0;
    int end = trees->start_coord;
    int mu_idx = 0, rho_idx = 0;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end() && k < nsites; ++it) {
        int start = end;
        end = start + it->blocklen;
        for (; k < nsites && blocks.sites[k] < end; k++) {
            int pos = blocks.sites[k];
            assert(pos >= start);

            ArgModel local_model;
            model->get_local_model(pos, local_model, &mu_idx, &rho_idx);
            double lk_now = likelihood_tree(it->tree, &local_model, seqs,
                                            no_base_probs, nleaves,
                                            pos, pos+1);
            double lk_swap = likelihood_tree(it->tree, &local_model, swapped,
                                             no_base_probs, nleaves,
                                             pos, pos+1);

            // the current orientation is "1-1" if ARG haplotype 1 carries
            // read-backed haplotype 1's allele
            if (sequences->seqs[blocks.hap1][pos] == blocks.allele1[k]) {
                log_emit11[k] = lk_now;
                log_emit12[k] = lk_swap;
            } else {
                log_emit11[k] = lk_swap;
                log_emit12[k] = lk_now;
            }
        }
    }
    assert(k == nsites);
}


int resample_phase(const ArgModel *model, Sequences *sequences,
                   const LocalTrees *trees, const PhaseBlocks &blocks,
                   double switch_prob, double point_prob)
{
    vector<double> log_emit11, log_emit12;
    vector<int> path;
    calc_phase_emissions(model, sequences, trees, blocks,
                         log_emit11, log_emit12);
    sample_phase_orientation(blocks, switch_prob, point_prob, log_emit11,
                             log_emit12, path);

    // in state "1-1" ARG haplotype 1 should carry allele1, otherwise not
    int nswapped = 0;
    for (int k=0; k<blocks.size(); k++) {
        int pos = blocks.sites[k];
        bool is11 = (sequences->seqs[blocks.hap1][pos] == blocks.allele1[k]);
        if (is11 != (path[k] == 0)) {
            sequences->switch_alleles(pos, blocks.hap1, blocks.hap2);
            nswapped++;
        }
    }
    return nswapped;
}


// split a line on tabs and spaces
static void split_whitespace(const char *line, vector<string> &tokens)
{
    tokens.clear();
    istringstream in(line);
    string token;
    while (in >> token)
        tokens.push_back(token);
}


// parse a block id: a non-negative integer, or '.' for unknown (-1)
static bool parse_block_id(const string &text, int *block)
{
    if (text == ".") {
        *block = -1;
        return true;
    }
    if (text.size() == 0)
        return false;
    char *end;
    long value = strtol(text.c_str(), &end, 10);
    if (*end != '\0' || value < 0 || value > INT_MAX)
        return false;
    *block = (int) value;
    return true;
}


bool read_phase_blocks(const char *filename, const Sequences *sequences,
                       const SitesMapping *sites_mapping,
                       vector<PhaseBlocks> &all_blocks)
{
    all_blocks.clear();
    CompressStream stream(filename);
    if (stream.stream == NULL) {
        printError("cannot read phase file '%s'", filename);
        return false;
    }

    // map from 0-based input positions to ARGweaver's coordinates
    map<int, int> site_coords;
    for (unsigned int i=0; i<sites_mapping->old_sites.size(); i++)
        site_coords[sites_mapping->old_sites[i]] = sites_mapping->new_sites[i];

    vector<string> names;
    map<int, vector<int> > rows;  // ARGweaver coordinate -> block ids
    int nignored = 0;
    int lineno = 0;
    int last_pos = -1;
    bool have_region = false;
    char *line;
    vector<string> tokens, ids;

    while ((line = fgetline(stream.stream))) {
        chomp(line);
        lineno++;
        split_whitespace(line, tokens);
        delete [] line;
        if (tokens.size() == 0)
            continue;

        if (tokens[0] == "NAMES") {
            names.assign(tokens.begin() + 1, tokens.end());
            if (names.size() == 0) {
                printError("phase file line %d: no names", lineno);
                return false;
            }
        } else if (tokens[0] == "REGION") {
            if (tokens.size() != 4) {
                printError("phase file line %d: REGION needs <chrom> <start>"
                           " <end>", lineno);
                return false;
            }
            int start = atoi(tokens[2].c_str()) - 1;
            int end = atoi(tokens[3].c_str());
            if (start != sites_mapping->old_start ||
                end != sites_mapping->old_end) {
                printError("phase file REGION (%s-%s) does not match the"
                           " sites REGION", tokens[2].c_str(),
                           tokens[3].c_str());
                return false;
            }
            have_region = true;
        } else {
            if (names.size() == 0 || !have_region) {
                printError("phase file line %d: NAMES and REGION must come"
                           " before the data", lineno);
                return false;
            }
            if (tokens.size() != 2) {
                printError("phase file line %d: expected a position and a"
                           " comma-separated list of block ids", lineno);
                return false;
            }
            char *end;
            long pos = strtol(tokens[0].c_str(), &end, 10);
            if (*end != '\0' || pos <= last_pos) {
                printError("phase file line %d: positions must be integers"
                           " in increasing order", lineno);
                return false;
            }
            last_pos = pos;

            split(tokens[1].c_str(), ",", ids);
            if (ids.size() != names.size()) {
                printError("phase file line %d: expected %d block ids,"
                           " found %d", lineno, (int) names.size(),
                           (int) ids.size());
                return false;
            }
            vector<int> blocks(ids.size());
            for (unsigned int j=0; j<ids.size(); j++) {
                if (!parse_block_id(ids[j], &blocks[j])) {
                    printError("phase file line %d: block id '%s' is not a"
                               " non-negative integer or '.'", lineno,
                               ids[j].c_str());
                    return false;
                }
            }

            // sites dropped by ARGweaver (or absent from the sites file)
            // are ignored
            map<int, int>::iterator it = site_coords.find(pos - 1);
            if (it == site_coords.end())
                nignored++;
            else
                rows[it->second] = blocks;
        }
    }

    if (names.size() == 0 || !have_region) {
        printError("phase file '%s' has no NAMES or REGION line", filename);
        return false;
    }

    // build phase blocks for each individual
    for (unsigned int j=0; j<names.size(); j++) {
        int hap1 = -1, hap2 = -1;
        for (int i=0; i<sequences->get_num_seqs(); i++) {
            if (sequences->names[i] == names[j] + "_1") hap1 = i;
            if (sequences->names[i] == names[j] + "_2") hap2 = i;
        }
        if (hap1 == -1 || hap2 == -1) {
            printError("phase file: cannot find sequences %s_1 and %s_2",
                       names[j].c_str(), names[j].c_str());
            return false;
        }

        // heterozygous sites missing from the phase file are unphased
        PhaseBlocks blocks(hap1, hap2);
        const char *seq1 = sequences->seqs[hap1];
        const char *seq2 = sequences->seqs[hap2];
        int nphased = 0;
        for (int i=0; i<sequences->length(); i++) {
            if (seq1[i] == seq2[i] || seq1[i] == 'N' || seq2[i] == 'N')
                continue;
            map<int, vector<int> >::iterator it = rows.find(i);
            int block = (it == rows.end()) ? -1 : it->second[j];
            blocks.sites.push_back(i);
            blocks.block.push_back(block);
            blocks.allele1.push_back(seq1[i]);
            if (block >= 0) nphased++;
        }

        int nblocks = 0;
        for (int k=0; k<blocks.size(); k++)
            if (blocks.block[k] >= 0 && !blocks.same_block(k))
                nblocks++;
        printLog(LOG_LOW, "phase blocks %s: %d heterozygous sites, %d in"
                 " %d blocks\n", names[j].c_str(), blocks.size(), nphased,
                 nblocks);
        all_blocks.push_back(blocks);
    }

    if (nignored > 0)
        printLog(LOG_LOW, "phase file: ignored %d positions that are not"
                 " among the sites used by ARGweaver\n", nignored);
    return true;
}

} // namespace argweaver
