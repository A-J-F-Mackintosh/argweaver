// C/C++ includes
#ifdef ARGWEAVER_MPI
#include "mpi.h"
#endif
#include <time.h>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>

// arghmm includes
#include "argweaver/compress.h"
#include "argweaver/ConfigParam.h"
#include "argweaver/emit.h"
#include "argweaver/fs.h"
#include "argweaver/logging.h"
#include "argweaver/mem.h"
#include "argweaver/parsing.h"
#include "argweaver/sample_arg.h"
#include "argweaver/sequences.h"
#include "argweaver/total_prob.h"
#include "argweaver/track.h"
#include "argweaver/est_popsize.h"
#include "argweaver/mcmcmc.h"
#include "argweaver/phase_hmm.h"
#include "argweaver/phase_flip.h"
#include "argweaver/phase_thread.h"
#include "argweaver/coal_records.h"
#include "argweaver/recomb.h"


using namespace argweaver;

// version info
#define VERSION_TEXT "1.0"
#define VERSION_INFO  "\
ARGweaver " VERSION_TEXT " \n\
Matt Rasmussen\n\
Sampler for large ancestral recombination graphs\n\
"

// file extensions
const char *SMC_SUFFIX = ".smc";
const char *SITES_SUFFIX = ".sites";
const char *STATS_SUFFIX = ".stats";
const char *LOG_SUFFIX = ".log";
const char *COAL_RECORDS_SUFFIX = ".cr";

// help categories
const int ADVANCED_OPT = 1;
const int POPMODEL_OPT = 2;
const int EXPERIMENTAL_OPT = 3;

const int EXIT_ERROR = 1;



// parsing command-line options
class Config
{
public:

    Config()
    {
        make_parser();

        resample_region[0] = -1;
        resample_region[1] = -1;
    }

    void make_parser()
    {
        config.clear();
        config.add(new ConfigParamComment("Input/output options"));
        // input/output
        config.add(new ConfigParam<string>
		   ("-o", "--output", "<output prefix>", &out_prefix,
                    "arg-sample",
                    "prefix for all output filenames (default='arg-sample')"));
	config.add(new ConfigParam<string>
		   ("-s", "--sites", "<sites alignment>", &sites_file,
		    "sequence alignment in sites format"));
	config.add(new ConfigParam<string>
		   ("-f", "--fasta", "<fasta alignment>", &fasta_file,
		    "sequence alignment in FASTA format"));
        config.add(new ConfigParam<string>
                   ("", "--vcf", "<.vcf.gz file>", &vcf_file,
                    "sequence alignment in gzipped vcf format. Must also supply"
                    " --region in format chr:start-end, and tabix index file"
                    " (.vcf.gz.tbi) must also be present. tabix program needs to"
                    " be installed and in PATH or path must be provided with"
                    " --tabix-dir. Assumes samples are diploid and unphased;"
                    " each individual will have _1 and _2 appended to its name for its"
                    " two haploid lineages. Indel-type variants are skipped."
                    " Any positions not specified in VCF are assumed to be"
                    " invariant"));
        config.add(new ConfigParam<string>
                   ("", "--vcf-files", "<vcf_file_list.txt>", &vcf_list_file,
                    "Same as --vcf but loads all individuals from multiple .vcf.gz files."
                    " The argument should be a file containing list of files to read"
                    " (one per line). All files should be aligned to same reference"
                    " genome"));
        config.add(new ConfigParam<string>
                   ("", "--rename-seqs", "<name_map_file.txt>", &rename_file,
                    "Used to rename sequences (usually from cryptic names in VCF"
                    " files to more meaningful names). The file should have two"
                    " columns, with the old sequence name, followed by the new"
                    " sequence name. Names can be diploid or haploid, and any old"
                    " names not matching the current sequence set will be ignored."));
        config.add(new ConfigParam<string>
                   ("", "--tabix-dir", "<directory>", &tabix_dir,
                    " path to tabix executable. May be required if using --vcf"
                    " and tabix is not in PATH"));
	config.add(new ConfigParam<string>
		   ("", "--age-file", "<age file>", &age_file,
		    " file giving age for any ancient samples (two-columns, "
		    " first column is sample name, second age in generations)."
		    " Age will be rounded to nearest discrete time point in"
		    " the model, with a maximum of the second-to-oldest point."));

        config.add(new ConfigSwitch
                   ("", "--use-genotype-probs", &use_genotype_probs,
                    "(for VCF input) Use genotype probabilities given by the"
                    " PL, GL, or PP scores in the VCF file. Default: ignore these scores and"
                    " treat assigned genotypes (after any filters) as"
                    " certain.", ADVANCED_OPT));


        config.add(new ConfigParam<string>
                   ("", "--region", "<start>-<end>",
                    &subregion_str, "",
                    "sample ARG for only a region of the sites (optional). Note"
                    "the [chr:] prefix should be added ONLY if alignment is in"
                    "vcf format. If this option is given, the vcf must also"
                    "be indexed with tabix"));
        config.add(new ConfigParam<string>
                   ("", "--subsites", "<subsites file>", &subsites_file,
                    "file listing NAMES from sites file (or sequences from"
                    " fasta) to keep; others will not be used. May be diploid"
                    " or haploid names (i.e., ind will have same effect as"
                    " ind_1 and ind_2)"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<SMC file>", &arg_file, "",
                    "initial ARG file (*.smc) for resampling (optional)"));
        /*        config.add(new ConfigParam<string>
                   ("", "--cr", "<CR file>", &cr_file, "",
                   "initial ARGfile (*.cf) for resampling (optional)"));*/

        config.add(new ConfigSwitch
                   ("", "--unphased", &unphased,
                    "data is unphased (will integrate over phasings)."));
        config.add(new ConfigSwitch
                   ("", "--phased-vcf", &phased_vcf,
                    "Treat data as phased even with VCF input (warning: there is no"
                    " support to read phasing information from VCF. This option is"
                    " for special case where VCF input is for inbred individuals"
                    " and the redundant haplotype for each individiual is removed"
                    " with --subsites"));
        config.add(new ConfigParam<string>
                   ("", "--unphased-file", "<filename>", &unphased_file, "",
                    "use this file to identify haplotype pairs (file should"
                    " have two sequence names per line). By default, will"
                    " first check for the naming convention <ind>_1 <ind>_2 to"
                    " determine pairs. If the convention is not used, will assume"
                    " haplotype pairs are adjacent to each other in sequence"
                    " file. This option should not be used with VCF files."));
        config.add(new ConfigParam<double>
                   ("", "--randomize-phase", "<frac_random>", &randomize_phase,
                    0.0, "randomize phasings at start", ADVANCED_OPT));
        config.add(new ConfigParam<double>
                   ("", "--add-switch-errors", "<switch_error_rate>", &switch_error_rate,
                    0.0, "add switch errors at this rate", ADVANCED_OPT));
        config.add(new ConfigSwitch
                   ("", "--no-sample-phase", &no_sample_phase,
                    "Do not sample phase. Otherwise, if data is unphased, phase"
                    " will be sampled at same frequency as ARGs", ADVANCED_OPT));
        config.add(new ConfigParam<int>
                   ("", "--sample-phase", "<niters>", &sample_phase_step, 0,
                    "output phasings every <niters> samples. The default"
                    " is to sample at same frequency as ARGs (see --sample-step),"
                    " unless --no-sample-phase is given",
                    ADVANCED_OPT));
        config.add(new ConfigParam<string>
                   ("", "--phase-blocks", "<phase file>", &phase_blocks_file,
                    "", "read-backed phase blocks (e.g. from WhatsHap) for"
                    " each individual. Phase is then resampled respecting"
                    " these blocks. The file follows the sites format, with"
                    " comma-separated block ids ('.' if unknown) instead of"
                    " bases, and individuals <ind> named in NAMES, whose"
                    " haplotypes are <ind>_1 and <ind>_2"));
        config.add(new ConfigParam<double>
                   ("", "--phase-switch-prob", "<s>", &phase_switch_prob,
                    0.01, "with read-backed phase blocks: probability of a"
                    " switch error between neighbouring heterozygous sites"
                    " in the same block (default=0.01)"));
        config.add(new ConfigParam<double>
                   ("", "--phase-point-error-prob", "<e>",
                    &phase_point_prob, -1.0, "with read-backed phase"
                    " blocks: probability that a single heterozygous site"
                    " in a block has its alleles the other way round, without"
                    " a switch of the block's phase (default: the same as"
                    " --phase-switch-prob)"));
        config.add(new ConfigParam<int>
                   ("", "--phase-flip-margin", "<bp>", &phase_flip_margin,
                    500, "with read-backed phase blocks: flip moves re-thread"
                    " a window around each boundary between blocks, extending"
                    " this far beyond the two heterozygous sites either side"
                    " of it (default=500)"));
        config.add(new ConfigParam<int>
                   ("", "--phase-flip-burnin", "<niters>", &phase_flip_burnin,
                    1000, "with read-backed phase blocks: number of burn-in"
                    " iterations, during which flip moves are made more often"
                    " (default=1000)"));
        config.add(new ConfigParam<double>
                   ("", "--phase-flip-freq-burnin", "<p>",
                    &phase_flip_freq_burnin, 1.0/3.0,
                    "with read-backed phase blocks: fraction of iterations"
                    " that are flip moves during burn-in (default=1/3)"));
        config.add(new ConfigParam<double>
                   ("", "--phase-flip-freq", "<p>", &phase_flip_freq, 0.1,
                    "with read-backed phase blocks: fraction of iterations"
                    " that are flip moves after burn-in (default=0.1).  A"
                    " flip move tries flipping the relative phase at every"
                    " boundary between blocks of one individual, then"
                    " resamples every individual's phase given the ARG.  The"
                    " other iterations are split evenly between re-threading"
                    " a haplotype jointly with its phase and ARGweaver's"
                    " internal moves"));

        config.add(new ConfigParamComment("Masking options"));
	config.add(new ConfigParam<string>
		   ("", "--maskmap", "<sites mask>",
                    &maskmap, "",
                    "mask map file (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--ind-maskmap", "<ind mask file>",
                    &ind_maskmap, "",
                    "File with two columns giving sample name and file with"
                    " mask for that sample. Sample name can refer to haploid or"
                    " diploid. If they are diploid, will append _1 and _2 to"
                    " names to search for relevant lineages to mask"));
        config.add(new ConfigSwitch
                   ("", "--expand-ind-mask",
                    &expand_ind_mask,
                    "This option is used to deal with Ns in a SITES file in the"
                    " case when the sites file was written by arg-sample run"
                    " with compression. In this case the N should apply to the"
                    " entire compressed site, and not a single base. Hint: to"
                    " run arg-sample and initialize the sites from a previous"
                    " run, use --sites <prev_outroot>[.rep].sites.gz"
                    " --maskmap <prev_outroot>.masked_regions.bed"
                    " --expand-ind-mask. The compression rate needs to be the"
                    " same as in the previous run. No other masking options are"
                    " necessary",
                    ADVANCED_OPT));
        config.add(new ConfigParam<int>
                   ("", "--mask-Ns", "<num>",
                    &maskN, -1,
                    "mask out any site in which at least num lineages have "
                    "missing data (this may increase efficiency)"));
        config.add(new ConfigParam<string>
                   ("", "--mask-cluster", "<numsnp,window>", &mask_cluster, "",
                    "mask any windows of size <window> which have at least"
                    " <numsnp> variants (at least two alleles other than N"
                    " must be observed at each site to count as a variant;"
                    " windows are determined after applying --subsites but"
                    " before any masking or compression is performed)"));
        config.add(new ConfigParam<string>
                   ("", "--vcf-genotype-filter", "<filter string>", &vcf_filter,
                    "String describing filtering for individual genotypes in VCF file."
                    " For example: \"GQ<5;DP<4;DP>30\" will mask genotypes with"
                    " GQ < 5 or DP < 4 or DP > 30. Currently the comparison operator"
                    " can only be < or >. All filtered genotypes will be masked."));
        config.add(new ConfigParam<double>
                   ("", "--vcf-min-qual", "<minQualScore>", &vcf_min_qual, 0.0,
                    "Minimum QUAL score for variants read from VCF. Others will be"
                    " masked. Default=0"));
        config.add(new ConfigParam<double>
                   ("", "--mask-uncertain", "<cutoff>", &mask_uncertain, 0.0,
                    "(for use --use-genotype-probs)"
                    " Mask any genotype if probability of most likely call < cutoff"
                    " according to PL score in VCF file", ADVANCED_OPT));

#ifdef ARGWEAVER_MPI
        config.add(new ConfigSwitch
                   ("", "--mpi", &mpi, "this is an mpi run, add <rank>.sites"
                    " to sites file name and <rank>. to out root, and"
                    " <rank>.smc.gz to --arg option (if given)"));
#endif

        // model parameters
        config.add(new ConfigParamComment("Model parameters"));
        config.add(new ConfigParam<double>
                   ("-m", "--mutrate", "<mutation rate>", &mu, 2.5e-8,
                    "mutations per site per generation (default=2.5e-8)"));
        config.add(new ConfigParam<string>
                   ("-M", "--mutmap", "<mutation rate map file>", &mutmap, "",
                    "mutation map file (optional)"));
        config.add(new ConfigParam<double>
                   ("-r", "--recombrate", "<recombination rate>", &rho, 1.5e-8,
                    "recombination per site per generation (default=1.5e-8)"));
        config.add(new ConfigParam<string>
                   ("-R", "--recombmap", "<recombination rate map file>",
                    &recombmap, "",
                    "recombination map file (optional)"));

        config.add(new ConfigParam<int>
                   ("-t", "--ntimes", "<ntimes>", &ntimes, 20,
                    "number of time points (default=20)"));
        config.add(new ConfigParam<double>
                   ("", "--maxtime", "<maxtime>", &maxtime, 200e3,
                    "maximum time point in generations (default=200e3)"));
        config.add(new ConfigParam<double>
                   ("", "--time-step", "<time>", &time_step, 0,
                    "linear time step in generations (optional)"));
        config.add(new ConfigParam<double>
                   ("", "--delta", "<delta>", &delta, 0.01,
                    "delta value for choosing log times (bigger value-> more"
                    " dense time points at leaves (default=0.01)"));
        config.add(new ConfigParam<string>
                   ("", "--popsize-file", "<popsize filename>", &popsize_file, "",
                    "file containing population sizes for each time span (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--times-file", "<times filename>", &times_file, "",
                    "file containing time points (optional)"));
        config.add(new ConfigParam<string>
                   ("-N", "--popsize", "<population size>", &popsize_str,
                    "10000",
                    "effective population size (default=1e4)"));
        config.add(new ConfigParam<string>
                   ("", "--popsize-file", "<popsize filename>", &popsize_file, "",
                    "Two column file with columns time, popsize."
                    " The time indicates the time (in generations from present) that"
                    " a population changes to the given diploid size."
                    " The rows should be sorted in order of increasing time, with"
                    " the first row having t=0.\n"
                    " If using the multiple population model, a third colum indicates"
                    "   the population number, with the first population numbered zero"));

        config.add(new ConfigSwitch
                   ("", "--smc-prime", &model.smc_prime,
                    "Use SMC' model (this option is generally slower"
                    " but represents a more accurate model). It is"
                    " recommended with the --pop-tree-file option, or for"
                    " ARGs with small numbers of highly diverged samples"));
        // note that --invisible-recombs has value=true,
        // so that providing the switch makes the variable false
        config.add(new ConfigSwitch
                   ("", "--invisible-recombs", &invisible_recombs,
                    "Output invisible recombinations which do not affect tree"
                    " topology (this will provide more accurate count of"
                    " recombination events, but may increase the runtime",
                    ADVANCED_OPT, false));
        config.add(new ConfigSwitch
                   ("", "--smc-orig", &do_nothing,
                    "This option is deprecated and does nothing; it is kept"
                    " for backwards compatibility. The original SMC model is the"
                    " default; see --smc-prime", ADVANCED_OPT));

        // Population model options
        config.add(new ConfigParam<string>
                   ("", "--pop-file", "<population assigment file>", &pop_file,
                    "file assigning each haplotype to a population index",
                    POPMODEL_OPT));
        config.add(new ConfigParam<string>
                   ("-P", "--pop-tree-file", "<population file>", &pop_tree_file,
                    "", "File describing population tree (for multiple populations)"
                    " (Note --smc-prime is recommended with this option)",
                    POPMODEL_OPT));
        config.add(new ConfigParam<int>
                   ("", "--max-migs", "<max_migs>", &max_migrations,
                    1, "For use with --pop-tree-file, do not thread lineages with more"
                    " than this many migrations. Setting this to any value other"
                    " than one is experimental, and may significantly increase"
                    " run-time and cause issues with MCMC mixing",
                    POPMODEL_OPT));
        config.add(new ConfigSwitch
                   ("", "--no-resample-mig", &no_resample_mig,
                    "Do not perform migration-specific resampling",
                    POPMODEL_OPT));
        config.add(new ConfigParam<int>
                   ("", "--start-mig", "<start_mig>", &start_mig_iter,
                    0, "Start allowing migration events at this MCMC iteration"
                    " <start_mig>", POPMODEL_OPT));


        // Population size sampling options
        config.add(new ConfigParamComment("Population size sampling", EXPERIMENTAL_OPT));
        config.add(new ConfigParam<int>
                   ("", "--sample-popsize", "<num>", &sample_popsize_num, 0,
                    "sample population size using Hamiltonian Monte Carlo every"
                    "<num> threading operations (default=0 means do not sample)",
                    EXPERIMENTAL_OPT));
	config.add(new ConfigParam<int>
		   ("", "--popsize-em", "<n>", &popsize_em, 0,
		    "Do EM update of popsizes after every n threading operations",
                    EXPERIMENTAL_OPT));
	config.add(new ConfigParam<double>
		   ("", "--popsize-em-min-event", "<num>", &popsize_em_min_event,
		    2000.0,
		    "Minimum number of events per time interval; time intervals with"
		    " fewer events will be combined with previous time interval for"
		    " EM computations", EXPERIMENTAL_OPT));
        config.add(new ConfigParam<int>
                   ("", "--popsize-config", "<num>", &popsize_config, 0,
                    "Choose configuration for population sizes:\n"
                    "  0: constant population size across time (default)\n"
                    "  1: different population size in each population\n"
                    "  2: most fine-grained; different population size in each time\n"
                    " interval and population",
                    EXPERIMENTAL_OPT));
        config.add(new ConfigSwitch
                   ("", "--sample-popsize-const", &sample_popsize_const,
                    "update popsize but keep constant across times/populations",
                    EXPERIMENTAL_OPT));
        config.add(new ConfigParam<string>
                   ("", "--sample-popsize-config", "<popsize config file>",
                    &popsize_config_file, "",
                    "optional, for use with --sample-popsize.\n"
                    " Overrides popsize-config option, which provides some common configs\n"
                    " This file should contain up to 2-5 entries per line, in the format:"
                    " <param_name> <time_idx> <pop_idx=0> <optimize=1> <init_val=N>\n"
                    " The last three values are optional. All entries with the same\n"
                    " param_name will be constrained to have the same population size\n"
                    " Population sizes are initialized to the default value of N unless\n"
                    " specified here.", EXPERIMENTAL_OPT));
	config.add(new ConfigParam<double>
		   ("", "--epsilon", "<val>", &epsilon,
		    0.01, "(for use with --sample-popsize) epsilon value for"
		    "Hamiltonian population size updates", EXPERIMENTAL_OPT));
        config.add(new ConfigParam<double>
		   ("", "--pseudocount", "<val>", &pseudocount,
		    1.0, "(for use with --sample-popsize) gives weight to prior",
		    EXPERIMENTAL_OPT));
#ifdef ARGWEAVER_MPI
        config.add(new ConfigParam<int>
                   ("", "--mcmcmc", "<int>", &mcmcmc_numgroup,
                    1, "number of mcmcmc threads",
                    EXPERIMENTAL_OPT));
        config.add(new ConfigParam<double>
                   ("", "--mcmcmc-heat", "<val>", &mcmcmc_heat,
                    0.05, "heat interval for each thread in (MC)^3 group",
                    EXPERIMENTAL_OPT));
#endif
        config.add(new ConfigSwitch
                   ("", "--init-popsize-random", &init_popsize_random,
                    "(for use with --sample-popsize). Initialize each"
                    " population size parameter to a random number sampled"
                    " uniformly in [5000, 50000]", EXPERIMENTAL_OPT));


        // sampling
        config.add(new ConfigParamComment("Sampling"));
        config.add(new ConfigParam<int>
                   ("-n", "--iters", "<# of iterations>", &niters, 1000,
                    "(default=1000)"));
        config.add(new ConfigParam<string>
                   ("", "--resample-region", "<start>-<end>",
                    &resample_region_str, "",
                    "region to resample of input ARG (optional)"));
        config.add(new ConfigSwitch
                   ("", "--resume", &resume, "resume a previous run"));
        config.add(new ConfigSwitch
                   ("", "--overwrite", &overwrite,
                    "force an overwrite of a previous run"));
	config.add(new ConfigSwitch
		   ("", "--no-sample-arg", &no_sample_arg,
		    "Do not sample the ARG; hold at initial value "
		    "(usually for use with --sample-popsize)",
		    EXPERIMENTAL_OPT));

        // misc
        config.add(new ConfigParamComment("Miscellaneous"));
        config.add(new ConfigParam<int>
                   ("-c", "--compress-seq", "<compression factor>",
                    &compress_seq, 1,
                    "alignment compression factor (default=1)"));
        config.add(new ConfigParam<int>
                   ("", "--climb", "<# of climb iterations>", &nclimb, 0,
                    "(default=0)", ADVANCED_OPT));
        config.add(new ConfigParam<int>
		   ("", "--num-buildup", "<# of buildup iterations>", &num_buildup,
                    1, "(default=0)", ADVANCED_OPT));
        config.add(new ConfigParam<int>
                   ("", "--sample-step", "<sample step size>", &sample_step,
                    10, "number of iterations between steps (default=10)"));
        config.add(new ConfigSwitch
                   ("", "--no-compress-output", &no_compress_output,
                    "do not gzip output files"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default=current time)"));
        config.add(new ConfigSwitch
                   ("", "--write-sites", &write_sites,
                    "Write sites (after all masking options) to <outroot>.sites.gz"));
        config.add(new ConfigSwitch
                   ("", "--write-sites-only", &write_sites_only,
                    "Same as --write-sites, but exit after reading and writing sites"));
        config.add(new ConfigSwitch
                   ("", "--include-masked-sites", &write_masked_sites,
                    "Relevant for --write-sites, --write-sites-only, or when sampling phase:"
                    " include sites that are masked in all samples in the output. By default"
                    " these are not included, but the set of masked sites is written to a"
                    " bed file named <outroot>.masked_sites.bed"));

        // advanced options
        config.add(new ConfigParamComment("Advanced Options", ADVANCED_OPT));
        config.add(new ConfigSwitch
                   ("", "--gibbs", &gibbs,
                    "use Gibbs sampling", ADVANCED_OPT));
        config.add(new ConfigParam<double>
                   ("", "--prob-path-switch", "<probability>",
                    &prob_path_switch, .1,
                    "removal path switch (default=.1)", ADVANCED_OPT));
        config.add(new ConfigSwitch
                   ("", "--infsites", &infsites,
                    "assume infinite sites model (at most one mutation per site)",
                    ADVANCED_OPT));
        config.add(new ConfigParam<int>
                   ("", "--resample-window", "<window size>",
                    &resample_window, 100000,
                    "sliding window for resampling (default=100000)",
                    ADVANCED_OPT));
        config.add(new ConfigParam<int>
                   ("", "--resample-window-iters", "<iterations>",
                    &resample_window_iters, 10,
                    "number of iterations per sliding window for resampling"
                    " (default=10)", ADVANCED_OPT));


        // help information
        config.add(new ConfigParamComment("Information"));
        config.add(new ConfigParam<int>
                   ("-V", "--verbose", "<verbosity level>",
                    &verbose, LOG_LOW,
                    "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
        config.add(new ConfigSwitch
                   ("-q", "--quiet", &quiet, "suppress logging to stderr"));
        config.add(new ConfigSwitch
                   ("-v", "--version", &version, "display version information"));
        config.add(new ConfigSwitch
                   ("-h", "--help", &help,
                    "display help information"));
        config.add(new ConfigSwitch
                   ("", "--help-advanced", &help_advanced,
                    "display help information about advanced options"));
        config.add(new ConfigSwitch
                   ("", "--help-experimental", &help_experimental,
                    "display help information about experimental options",
                    ADVANCED_OPT));
        config.add(new ConfigSwitch
                   ("-d", "--help-popmodel", &help_popmodel,
                    "display help information for population-based models (ARGweaver-D)"));
    }

    int parse_args(int argc, char **argv)
    {
        // parse arguments
        if (!config.parse(argc, (const char**) argv)) {
            if (argc < 2)
                config.printHelp();
            return EXIT_ERROR;
        }

        // display help
        if (help) {
            config.printHelp();
            return EXIT_ERROR;
        }

        // display debug help
        if (help_advanced) {
            config.printHelp(stderr, ADVANCED_OPT);
            return EXIT_ERROR;
        }

        if (help_popmodel) {
            config.printHelp(stderr, POPMODEL_OPT);
            return EXIT_ERROR;
        }

        if (help_experimental) {
            config.printHelp(stderr, EXPERIMENTAL_OPT);
            return EXIT_ERROR;
        }

        // display version info
        if (version) {
            printf(VERSION_INFO);
            return EXIT_ERROR;
        }
#ifdef ARGWEAVER_MPI
        mcmcmc_group = 0;
        int groupsize = MPI::COMM_WORLD.Get_size() / mcmcmc_numgroup;
        mcmcmc_group = MPI::COMM_WORLD.Get_rank() / groupsize;
        if (mcmcmc_group != 0) {
            char tmp[1000];
            sprintf(tmp, ".%i", mcmcmc_group);
            mcmcmc_prefix = string(tmp);
        }
        printLog(LOG_LOW, "mcmcmc_prefix = %s\n", mcmcmc_prefix.c_str());
        printLog(LOG_LOW, "mcmcmc_group=%i\n", mcmcmc_group);
#endif

        return 0;
    }

    ConfigParser config;

    // input/output
    string fasta_file;
    string sites_file;
    string vcf_file;
    string vcf_list_file;
    string rename_file;
    string vcf_filter;
    double vcf_min_qual;
    string subsites_file;
    string out_prefix;
    string arg_file;
    //    string cr_file;
    string subregion_str;
    string age_file;
    string maskmap;
    string ind_maskmap;


    bool expand_ind_mask;
    int maskN;
    string mask_cluster;
    string tabix_dir;
    double mask_uncertain;
    bool use_genotype_probs;

    // model parameters
    string popsize_str;
    string pop_tree_file;
    int max_migrations;
    string pop_file;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    double time_step;
    double delta;
    string popsize_file;
    string times_file;
    string mutmap;
    string recombmap;
    ArgModel model;
    int popsize_em;
    double popsize_em_min_event;
    bool popsize_prior_neighbor;
    bool init_popsize_random;
    int popsize_config;
    string popsize_config_file;
    int sample_popsize_num;
    bool sample_popsize_const;
    bool invisible_recombs;
    double epsilon;
    double pseudocount;

#ifdef ARGWEAVER_MPI
    double mcmcmc_heat;
    int mcmcmc_group;
    int mcmcmc_numgroup;
    bool mpi;
#endif
    string mcmcmc_prefix;
    bool no_sample_arg;
    bool no_resample_mig;
    int start_mig_iter;

    // search
    int nclimb;
    int num_buildup;
    int niters;
    string resample_region_str;
    int resample_region[2];
    bool resume;
    bool overwrite;
    string resume_stage;
    int resume_iter;
    int resample_window;
    int resample_window_iters;
    bool gibbs;

    // misc
    int compress_seq;
    int sample_step;
    bool no_compress_output;
    int randseed;
    double prob_path_switch;
    bool infsites;
    bool unphased;
    string unphased_file;
    bool phased_vcf;
    double randomize_phase;
    double switch_error_rate;
    bool no_sample_phase;
    int sample_phase_step;
    string phase_blocks_file;
    double phase_switch_prob;
    double phase_point_prob;
    int phase_flip_margin;
    int phase_flip_burnin;
    double phase_flip_freq_burnin;
    double phase_flip_freq;
    vector<PhaseBlocks> phase_blocks;
    bool all_masked;
    bool write_sites;
    bool write_sites_only;
    bool write_masked_sites;

    // help/information
    bool quiet;
    int verbose;
    bool version;
    bool help;
    bool help_advanced;
    bool help_experimental;
    bool help_popmodel;
    bool do_nothing;

    // logging
    FILE *stats_file;
};



bool parse_region(const char *region, int *start, int *end)
{
    int start_idx=0;
    for (int i=0; i < (int)strlen(region); i++)
        if (region[i]==':') {
            start_idx = i+1;
            break;
        }

    return sscanf(&region[start_idx], "%d-%d", start, end) == 2;
}

//=============================================================================
// logging

void set_up_logging(const Config &c, int level, const char *log_mode) {
    Logger *logger;
    /*    printf("set_up_logging %s %s %s %i %i\n",
           c.out_prefix.c_str(), c.mcmcmc_prefix.c_str(), LOG_SUFFIX, c.quiet,
           level);*/
    setLogLevel(level);
    if (c.quiet) {
        // log only to file
        logger = &g_logger;
    } else {
        // log to both stdout and file
        logger = new Logger(NULL, level);
        g_logger.setChain(logger);
    }
    string log_filename = c.out_prefix + c.mcmcmc_prefix + LOG_SUFFIX;
    if (!logger->openLogFile(log_filename.c_str(), log_mode)) {
        printError("Could not open log file '%s'", log_filename.c_str());
        abort();
    }
}

// log the program version and start time
void log_intro(int level)
{
    time_t t = time(NULL);

    printLog(level, "arg-sample " VERSION_TEXT "\n");
    printLog(level, "start time: %s", ctime(&t));  // newline is in ctime
}

// log the command used
void log_prog_commands(int level, int argc, char **argv)
{
    printLog(level, "command:");
    for (int i=0; i<argc; i++) {
        printLog(level, " %s", argv[i]);
    }
    printLog(level, "\n");
}


//=============================================================================
// statistics output

void print_stats_header(Config *config) {
    fprintf(config->stats_file, "stage\titer\tprior\tprior2\tlikelihood\tjoint\t"
            "recombs\tnoncompats\targlen");
    if (config->invisible_recombs)
        fprintf(config->stats_file, "\tinvis_recombs");
    if (config->model.popsize_config.sample) {
        list<PopsizeConfigParam> l = config->model.popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin();
             it != l.end(); ++it) {
            fprintf(config->stats_file, "\t%s", it->name.c_str());
        }
    } else if (config->popsize_em) {
        for (int pop=0; pop < config->model.num_pops(); pop++) {
            for (int i=0; i < config->model.ntimes-1; i++) {
                char str[100];
                if (config->model.num_pops() == 1)
                    sprintf(str, "N%i", i);
                else sprintf(str, "N%i.%i", pop, i);
                fprintf(config->stats_file, "\t%s", str);
            }
        }
    }
    if (config->model.pop_tree != NULL) {
        for (unsigned int i=0; i < config->model.pop_tree->mig_params.size(); i++)
            fprintf(config->stats_file, "\t%s", config->model.pop_tree->mig_params[i].name.c_str());
    }
    fprintf(config->stats_file, "\n");
}


void print_stats(FILE *stats_file, const char *stage, int iter,
                 ArgModel *model,
                 const Sequences *sequences, LocalTrees *trees,
                 const SitesMapping* sites_mapping, const Config *config,
                 const TrackNullValue *maskmap_uncompressed,
                 const vector<int> &invisible_recomb_pos0=vector<int>(),
                 const vector<Spr> &invisible_recombs=vector<Spr>())
{

    // calculate number of recombinations
    int nrecombs = trees->get_num_trees() - 1;

    // calculate number of non-compatiable sites
    int nseqs = trees->get_num_leaves();
    char *seqs[nseqs];
    for (int i=0; i<nseqs; i++)
        seqs[i] = sequences->seqs[trees->seqids[i]];
    int noncompats = count_noncompat(trees, seqs, nseqs, sequences->length());

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;


    // calculate likelihood, prior, and joint probabilities
    // uncompressed local trees
    vector<int> invisible_recomb_pos;
    if (sites_mapping) {
        uncompress_local_trees(trees, sites_mapping);
	uncompress_model(model, sites_mapping, config->compress_seq);
        sites_mapping->uncompress(invisible_recomb_pos0, invisible_recomb_pos);
    }

    double prior = calc_arg_prior(model, trees, NULL, NULL, -1, -1,
                                  invisible_recomb_pos, invisible_recombs);
    double prior2 = calc_arg_prior_recomb_integrate(model, trees, NULL, NULL, NULL);
    double likelihood = config->all_masked ? 0.0 :
        calc_arg_likelihood(model, sequences, trees,
                            sites_mapping,
                            maskmap_uncompressed);
    double joint = prior + likelihood;
    double arglen = get_arglen(trees, model->times);

    // recompress local trees
    if (sites_mapping) {
        compress_local_trees(trees, sites_mapping);
	compress_model(model, sites_mapping, config->compress_seq);
    }

    // output stats
    fprintf(stats_file, "%s\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%f",
            stage, iter,
            prior, prior2, likelihood, joint, nrecombs, noncompats, arglen);
    if (config->invisible_recombs)
        fprintf(stats_file, "\t%i", (int)invisible_recomb_pos.size());
    if (model->popsize_config.sample) {
        list<PopsizeConfigParam> l=model->popsize_config.params;
        for (list<PopsizeConfigParam>::iterator it=l.begin();
             it != l.end(); ++it) {
            set<PopTime>::iterator it2 = it->intervals.begin();
            fprintf(stats_file, "\t%f", model->popsizes[it2->pop][it2->time]);
            it2++;
            // just checking here
            set<PopTime>::iterator it3 = it2;
            while (it3 != it->intervals.end()) {
                assert(model->popsizes[it2->pop][it2->time] ==
                       model->popsizes[it3->pop][it3->time]);
                it3++;
            }
        }
    } else if (config->popsize_em) {
        for (int pop=0; pop < model->num_pops(); pop++)
            for (int i=0; i < model->ntimes-1; i++)
                fprintf(stats_file, "\t%f", model->popsizes[pop][2*i]);
    }
    if (model->pop_tree != NULL) {
        for (unsigned int i=0; i < model->pop_tree->mig_params.size(); i++) {
            MigParam mp = model->pop_tree->mig_params[i];
            fprintf(stats_file, "\t%e",
                    model->pop_tree->mig_matrix[mp.time_idx].get(mp.from_pop, mp.to_pop));
        }
    }
    fprintf(stats_file, "\n");
    fflush(stats_file);

    printLog(LOG_LOW, "\n"
             "prior:      %f\n"
             "likelihood: %f\n"
             "joint:      %f\n"
             "nrecombs:   %d\n"
             "noncompats: %d\n"
             "arglen:     %f\n"
             "max memory: %.1f MB\n\n",
             prior, likelihood, joint, nrecombs, noncompats, arglen, maxrss);

}

//=============================================================================
// sample output


// Returns the iteration-specific ARG filename
string get_out_arg_file(const Config &config, int iter)
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    return config.out_prefix + config.mcmcmc_prefix + iterstr + SMC_SUFFIX;
}

/*string get_out_cr_file(const Config &config, int iter)
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    return config.out_prefix + config.mcmcmc_prefix + iterstr + COAL_RECORDS_SUFFIX;
    }*/

string get_out_sites_file(const Config &config, int iter)
{
    char iterstr[12];
    if (iter < 0) {
        iterstr[0]='\0';
    } else {
        snprintf(iterstr, 12, ".%d", iter);
    }
    string sitesfile = config.out_prefix + config.mcmcmc_prefix + iterstr + SITES_SUFFIX;
    if (!config.no_compress_output)
        sitesfile = sitesfile + ".gz";
    return sitesfile;
}

bool log_sequences(string chrom, const Sequences *sequences,
                   const Config *config,
                   const SitesMapping *sites_mapping, int iter) {
    Sites sites(chrom);
    string out_sites_file = get_out_sites_file(*config, iter);
    make_sites_from_sequences(sequences, &sites);
    if (sites_mapping)
        uncompress_sites(&sites, sites_mapping);
    CompressStream stream(out_sites_file.c_str(), "w");
    if (!stream.stream) {
        printError("cannot write '%s'", out_sites_file.c_str());
        return false;
    }
    write_sites(stream.stream, &sites, config->write_masked_sites);
    return true;
}

bool log_local_trees(const ArgModel *model, const Sequences *sequences,
                     LocalTrees *trees, const SitesMapping* sites_mapping,
                     const Config *config, int iter,
                     const vector<int> &self_recomb_pos0=vector<int>(),
                     const vector<Spr> &self_recombs=vector<Spr>())
{
    string out_arg_file = get_out_arg_file(*config, iter);
    const vector<int> *self_recomb_ptr;
    if (!config->no_compress_output)
        out_arg_file += ".gz";

    // write local trees uncompressed
    vector<int> self_recomb_pos1;
    if (sites_mapping) {
        uncompress_local_trees(trees, sites_mapping);
        sites_mapping->uncompress(self_recomb_pos0, self_recomb_pos1);
        self_recomb_ptr = &self_recomb_pos1;
    } else self_recomb_ptr = &self_recomb_pos0;

    // setup output stream
    CompressStream stream(out_arg_file.c_str(), "w");
    if (!stream.stream) {
        printError("cannot write '%s'", out_arg_file.c_str());
        return false;
    }


    write_local_trees(stream.stream, trees, sequences, model->times,
                      model->pop_tree != NULL,
                      *self_recomb_ptr, self_recombs);

    // testing for now; output coal records version
    /*    string out_cr_file = get_out_cr_file(*config, iter);
    if (!config->no_compress_output)
        out_cr_file += ".gz";
    CompressStream stream2(out_cr_file.c_str(), "w");
    if (!stream2.stream) {
        printError("cannot write '%s'", out_cr_file.c_str());
        return false;
    }

    write_coal_records(stream2.stream, model, trees, sequences); */

    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);

    return true;
}


//=============================================================================


// Read a list of doubles from a file.
bool read_doubles(const char *filename, vector<double> &values) {
    FILE *infile = fopen(filename, "r");
    if (infile == NULL) {
        return false;
    }

    // Read all doubles and sort them.
    double value;
    while (fscanf(infile, "%lf", &value) == 1)
        values.push_back(value);
    fclose(infile);

    return true;
}


/*bool read_init_arg_cr(const char *cr_file, const ArgModel *model,
                   LocalTrees *trees, vector<string> &seqnames)
{
    CompressStream stream(cr_file, "r");
    if (!stream.stream) {
        printError("cannot read '%s'", cr_file);
        return false;
    }
    return read_coal_records(stream.stream, model, trees, seqnames);
    }*/

bool read_init_arg(const char *arg_file, const ArgModel *model,
                   LocalTrees *trees, vector<string> &seqnames)
{
    CompressStream stream(arg_file, "r");
    if (!stream.stream) {
        printError("cannot read '%s'", arg_file);
        return false;
    }
    return read_local_trees(stream.stream, model->times, model->ntimes,
                            trees, seqnames);
}



//=============================================================================
// sampling methods

// build initial arg by sequential sampling
bool seq_sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                    SitesMapping* sites_mapping, Config *config,
                    const TrackNullValue *maskmap_orig)
{

    if (trees->get_num_leaves() < sequences->get_num_seqs()) {
        printLog(LOG_LOW, "Sequentially Sample Initial ARG (%d sequences)\n",
                 sequences->get_num_seqs());
        printLog(LOG_LOW, "------------------------------------------------\n");
        sample_arg_seq(model, sequences, trees, true, config->num_buildup);
        print_stats(config->stats_file, "seq", trees->get_num_leaves(),
                    model, sequences, trees, sites_mapping, config,
                    maskmap_orig);
	return true;
    }
    return false;
}


void climb_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
               SitesMapping* sites_mapping, Config *config,
               const TrackNullValue *maskmap_orig)
{
    if (config->resume)
        return;

    printLog(LOG_LOW, "Climb Search (%d iterations)\n", config->nclimb);
    printLog(LOG_LOW, "-----------------------------\n");
    double recomb_preference = .9;
    for (int i=0; i<config->nclimb; i++) {
        printLog(LOG_LOW, "climb %d\n", i+1);
        resample_arg_climb(model, sequences, trees, recomb_preference);
        print_stats(config->stats_file, "climb", i, model, sequences, trees,
                    sites_mapping, config, maskmap_orig);
    }
    printLog(LOG_LOW, "\n");
}

void mcmcmc_swap(Config *config, ArgModel *model, const Sequences *sequences,
                 const LocalTrees *trees, const SitesMapping *sites_mapping) {
#ifdef ARGWEAVER_MPI
    printLog(LOG_LOW, "mcmcmc_swap model->mc3.max_group=%i\n", model->mc3.max_group);
    if (model->mc3.max_group == 0) return;
    Mc3Config *mc3 = &(model->mc3);
    int swap[2];
    int rank = MPI::COMM_WORLD.Get_rank();
    bool accept;
    double swapstats[4];
    //pick two groups to swap
    if (rank == 0) {
        swap[0] = irand(mc3->max_group + 1);
        swap[1] = irand(mc3->max_group);
        if (swap[1] >= swap[0]) swap[1]++;
        assert(swap[0] >= 0 && swap[0] <= mc3->max_group);
        assert(swap[1] >= 0 && swap[1] <= mc3->max_group);
        assert(swap[0] != swap[1]);
    }
    MPI::COMM_WORLD.Bcast(swap, 2, MPI::INT, 0);
    if (mc3->group == swap[0] || mc3->group == swap[1]) {
        double vals[2];
        vals[0] = calc_arg_prior(model, trees) +
            calc_arg_likelihood(model, sequences, trees, sites_mapping);
        vals[1] = mc3->heat;
        if (mc3->group_comm->Get_rank()==0)
            mc3->group_comm->Reduce(MPI_IN_PLACE, vals, 1, MPI::DOUBLE, MPI_SUM,
                                   0);
        else mc3->group_comm->Reduce(vals, vals, 1, MPI::DOUBLE, MPI_SUM, 0);
        if (mc3->group_comm->Get_rank()==0)
            MPI::COMM_WORLD.Isend(vals, 2, MPI::DOUBLE, 0, 301);
    }
    if (rank == 0) {
        double like[2];
        double heat[2];
        for (int i=0; i < 2; i++) {
            double vals[2];
            MPI::COMM_WORLD.Recv(vals, 2, MPI::DOUBLE, MPI_ANY_SOURCE, 301);
            like[i] = vals[0];
            heat[i] = vals[1];
        }
        double accept_ratio = (heat[0]-heat[1])*like[1] +
            (heat[1] - heat[0])*like[0];
        accept = (accept_ratio >= 0.0 ||
                  frand() < exp(accept_ratio));
        //Note heat[0] may correspond to swap[0] or swap[1]
        swapstats[0] = accept ? 1.0 : 0.0;
        swapstats[1] = heat[0];
        swapstats[2] = heat[1];
        swapstats[3] = accept_ratio;
    }
    MPI::COMM_WORLD.Bcast(swapstats, 4, MPI::DOUBLE, 0);
    accept = ( swapstats[0] > 0.5);
    printLog(LOG_LOW, "swap\t%i\t%i\t%f\t%f\t%f\t%s\n",
             swap[0], swap[1], swapstats[1], swapstats[2], swapstats[3],
             accept ? "accept" : "reject");

    // now all processes know accept. Do the switch
    if (accept && (mc3->group == swap[0] || mc3->group == swap[1])) {
        int other = (mc3->group == swap[1] ? 0 : 1);
        mc3->group = swap[other];
        mc3->heat = 1.0 - mc3->heat_interval * mc3->group;
        //need to switch output files as well, including stats_file, arg output,
        //phase output, log files.  First close all the files.
        fclose(config->stats_file);
        if (config->verbose) {
            Logger *chain = g_logger.getChain();
            chain->closeLogFile();
        }
        if (mc3->group == 0) config->mcmcmc_prefix = "";
        else {
            char tmp[1000];
            sprintf(tmp, ".%i", mc3->group);
            config->mcmcmc_prefix = string(tmp);
        }
    }
    MPI::COMM_WORLD.Barrier();
    //Now reopen in append mode
    if (accept && (mc3->group == swap[0] || mc3->group == swap[1])) {
        string stats_filename = config->out_prefix + config->mcmcmc_prefix
            + STATS_SUFFIX;
        if (!(config->stats_file = fopen(stats_filename.c_str(), "a"))) {
            printError("Error reopening stats file %s in mcmcmc_swap\n",
                    stats_filename.c_str());
            abort();
        }
        if (config->verbose) {
            string log_filename = config->out_prefix + config->mcmcmc_prefix
                + LOG_SUFFIX;
            Logger *chain = g_logger.getChain();
            if (!chain->openLogFile(log_filename.c_str(), "a")) {
                fprintf(stderr, "Error opening %s\n", log_filename.c_str());
                abort();
            }
        }
    }
#endif
}


// One iteration with read-backed phase blocks.  Choose one of three moves:
//   - flip: try flipping the relative phase at each boundary between blocks
//     of one individual, then resample every individual's phase given the
//     ARG (a fraction of iterations set by --phase-flip-freq-burnin during
//     burn-in and --phase-flip-freq after);
//   - re-thread a haplotype jointly with its individual's phase;
//   - ARGweaver's internal moves.
// The last two share the remaining iterations evenly.
void resample_arg_phase_blocks(ArgModel *model, Sequences *sequences,
                               LocalTrees *trees, Config *config, int iter,
                               int window, int niters, double heat)
{
    const vector<PhaseBlocks> &all_blocks = config->phase_blocks;
    const double flip_freq = (iter < config->phase_flip_burnin) ?
        config->phase_flip_freq_burnin : config->phase_flip_freq;

    if (frand() < flip_freq) {
        const PhaseBlocks &blocks = all_blocks[irand(all_blocks.size())];
        int ntried;
        int naccept = flip_phase_boundaries(
            model, sequences, trees, blocks, config->phase_flip_margin,
            &ntried);
        printLog(LOG_LOW, "flip_phase_boundaries %s %s: %d of %d flips"
                 " accepted\n", sequences->names[blocks.hap1].c_str(),
                 sequences->names[blocks.hap2].c_str(), naccept, ntried);

        for (unsigned int j=0; j<all_blocks.size(); j++) {
            const PhaseBlocks &blocks = all_blocks[j];
            int nswapped = resample_phase(
                model, sequences, trees, blocks, config->phase_switch_prob,
                config->phase_point_prob);
            printLog(LOG_LOW, "resample_phase %s %s: %d of %d sites"
                     " swapped\n", sequences->names[blocks.hap1].c_str(),
                     sequences->names[blocks.hap2].c_str(), nswapped,
                     blocks.size());
        }
    } else if (irand(2) == 0) {
        int node = irand(trees->get_num_leaves());
        int seqid = trees->seqids[node];
        const PhaseBlocks *blocks = NULL;
        for (unsigned int j=0; j<all_blocks.size(); j++)
            if (all_blocks[j].hap1 == seqid || all_blocks[j].hap2 == seqid)
                blocks = &all_blocks[j];
        if (blocks) {
            int nswapped = resample_arg_leaf_phase(
                model, sequences, trees, node, *blocks,
                config->phase_switch_prob, config->phase_point_prob);
            printLog(LOG_LOW, "resample_arg_leaf_phase %s: %d of %d sites"
                     " swapped\n", sequences->names[seqid].c_str(),
                     nswapped, blocks->size());
        } else {
            // a haplotype without a partner
            resample_arg_leaf(model, sequences, trees, node);
            printLog(LOG_LOW, "resample_arg_leaf: accept=%f\n", 1.0);
        }
    } else {
        resample_arg_mcmc_all(model, sequences, trees, false, window, niters,
                              heat, config->no_resample_mig);
    }
}


void resample_arg_all(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                      SitesMapping* sites_mapping, Config *config,
                      const TrackNullValue *maskmap_orig)
{
    // setup search options
    bool do_leaf[config->niters+1];
    //int window = 100000;
    //int niters = 10;
    int window = config->resample_window;
    int niters = config->resample_window_iters;
    window /= config->compress_seq;

    vector<int> invisible_recomb_pos;
    vector<Spr> invisible_recombs;
    assert_trees(trees, model->pop_tree);

    // set iteration counter
    int iter = 1;
    if (config->resume)
        iter = config->resume_iter + 1;
    else {
        // save first ARG (iter=0)
        print_stats(config->stats_file, "resample", 0, model, sequences, trees,
                    sites_mapping, config, maskmap_orig);
        log_local_trees(model, sequences, trees, sites_mapping, config, 0,
                        invisible_recomb_pos, invisible_recombs);
        if (config->sample_phase_step > 0)
            log_sequences(trees->chrom, sequences, config, sites_mapping, 0);
    }


    printLog(LOG_LOW, "Resample All Branches (%d iterations)\n",
             config->niters);
    printLog(LOG_LOW, "--------------------------------------\n");

    double frac_leaf = 0.5;

#ifdef ARGWEAVER_MPI
    for (int i=0; i <= config->niters; i++) do_leaf[i] = (frand() < frac_leaf);
    MPI::COMM_WORLD.Bcast(do_leaf, config->niters+1, MPI::BOOL, 0);
#endif

    for (int i=iter; i<=config->niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i);
        Timer timer;
        double heat = model->mc3.heat;
        if (model->pop_tree != NULL && i >= config->start_mig_iter) {
            if (model->pop_tree->max_migrations != config->max_migrations)
                printLog(LOG_LOW, "Changing max_migrations to %i\n", config->max_migrations);
            model->pop_tree->max_migrations = config->max_migrations;
        }

#ifndef ARGWEAVER_MPI
        do_leaf[i] = ( frand() < frac_leaf );
#endif

	if (config->phase_blocks.size() > 0) {
            resample_arg_phase_blocks(model, sequences, trees, config, i,
                                      window, niters, heat);
        } else if ( ! config->no_sample_arg) {
	    if (config->gibbs)
		resample_arg(model, sequences, trees);
	    else
		resample_arg_mcmc_all(model, sequences, trees, do_leaf[i],
				      window, niters, heat,
                                      config->no_resample_mig);
	}



            // TODO: implement popsize updates
            /*	if (config->popsize_em > 0 && i % config->popsize_em == 0)
	    mle_popsize(model, trees, config->popsize_em_min_event);
            else */
        if (model->popsize_config.sample > 0 && i % model->popsize_config.sample == 0) {
            resample_popsizes_mh(model, trees, true, heat);
            //	    update_popsize_hmc(model, trees);
        } /*else {
            printError("Have not implemented popsize update for multipop yet\n");
            //no_update_popsize(model, trees);
            }*/

        printTimerLog(timer, LOG_LOW, "sample time:");

        mcmcmc_swap(config, model, sequences, trees, sites_mapping);

        if (model->smc_prime && config->invisible_recombs) {
            sample_invisible_recombinations(model, trees,
                                            invisible_recomb_pos,
                                            invisible_recombs);
        }

        if (model->pop_tree != NULL) {
            resample_migrates(model, trees,
                              invisible_recombs);
        }


        // logging
        print_stats(config->stats_file, "resample", i, model, sequences, trees,
                    sites_mapping, config, maskmap_orig,
                    invisible_recomb_pos, invisible_recombs);

        // sample saving
        if (i % config->sample_step == 0 && ! config->no_sample_arg)
            log_local_trees(model, sequences, trees, sites_mapping, config, i,
                            invisible_recomb_pos, invisible_recombs);

        if (config->sample_phase_step > 0 && i%config->sample_phase_step == 0)
            log_sequences(trees->chrom, sequences, config, sites_mapping, i);
    }
    printLog(LOG_LOW, "\n");
}


// overall sampling workflow
void sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                SitesMapping* sites_mapping, Config *config,
                const TrackNullValue *maskmap_orig)
{

    if (!config->resume)
        print_stats_header(config);

    // build initial arg by sequential sampling
    bool seq_sample = seq_sample_arg(model, sequences, trees, sites_mapping, config,
                   maskmap_orig);

    // re-thread any ancient samples using internal threading to make use of minage
    if (seq_sample && sequences->ages.size() > 0) {
	for (int i=0; i < (int)sequences->ages.size(); i++) {
	    if (sequences->ages[i] > 0) {
		int mintime = sequences->ages[i];
		for (int j=0; j < trees->get_num_leaves(); j++) {
		    if (trees->seqids[j] == i) {
			printLog(LOG_LOW, "Re-threading ancient sample %s to set sample age %i (%.1f)\n",
				 sequences->names[i].c_str(), mintime,
				 model->times[mintime]);
                        resample_arg_leaf(model, sequences, trees, j);
			break;
		    }
		}
	    }
	}
    }

    if (config->resample_region[0] != -1) {
        // region sampling
        printLog(LOG_LOW, "Resample Region (%d-%d, %d iterations)\n",
                 config->resample_region[0], config->resample_region[1],
                 config->niters);
        printLog(LOG_LOW, "--------------------------------------------\n");

        print_stats(config->stats_file, "resample_region", 0,
                    model, sequences, trees, sites_mapping, config,
                    maskmap_orig);

        for (int i=0; i < config->niters; i++) {
            resample_arg_region(model, sequences, trees,
                                config->resample_region[0],
                                config->resample_region[1],
                                1);

            // logging
            print_stats(config->stats_file, "resample_region", config->niters,
                        model, sequences, trees, sites_mapping, config,
                        maskmap_orig);
            log_local_trees(model, sequences, trees, sites_mapping, config, i);
        }

    } else{
        // climb sampling
        climb_arg(model, sequences, trees, sites_mapping, config,
                  maskmap_orig);
        // resample all branches
        resample_arg_all(model, sequences, trees, sites_mapping, config,
                         maskmap_orig);
    }
}


//=============================================================================

bool parse_status_line(const char* line, Config &config,
                       string &stage, int &iter, string &arg_file,
                       vector<string> header)
{
    bool found=false;
    // parse stage and last iter
    vector<string> tokens;
    split(line, "\t", tokens);
    if (tokens.size() < 2) {
        printError("incomplete line in status file");
        return false;
    }

    string stage2 = tokens[0];
    int iter2;
    if (sscanf(tokens[1].c_str(), "%d", &iter2) != 1) {
        printError("iter column is not an integer");
        return false;
    }

    // NOTE: only resume resample stage for now
    if (stage2 != "resample")
        return true;

    // see if ARG file exists
    string out_arg_file = get_out_arg_file(config, iter2);
    struct stat st;
    if (stat(out_arg_file.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_arg_file;
        found=true;
    }

    // try compress output
    out_arg_file += ".gz";
    if (stat(out_arg_file.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_arg_file;
        found=true;
    }

    if (found)
        config.model.init_params_from_statfile(header, line,
                                               config.popsize_em);
    return true;
}


bool setup_resume(Config &config)
{
    if (!config.resume)
        return true;

    printLog(LOG_LOW, "Resuming previous run\n");

    // open stats file
    string stats_filename = config.out_prefix + config.mcmcmc_prefix
        + STATS_SUFFIX;
    printLog(LOG_LOW, "Checking previous run from stats file: %s\n",
             stats_filename.c_str());

    FILE *stats_file;
    if (!(stats_file = fopen(stats_filename.c_str(), "r"))) {
        printError("could not open stats file '%s'",
                   stats_filename.c_str());
        return false;
    }

    // find last line of stats file that has a written ARG
    char *line = NULL;

    line = fgetline(stats_file);
    if (!line) {
        printError("status file is empty");
        return false;
    }
    chomp(line);
    vector<string> header;
    split(line, "\t", header);
    delete [] line;

    // loop through status lines
    string arg_file = "";
    while ((line = fgetline(stats_file))) {
        if (!parse_status_line(
             line, config, config.resume_stage, config.resume_iter, arg_file,
             header)) {
            delete [] line;
            return false;
        }
        delete [] line;
    }

    if (arg_file == "") {
        printLog(LOG_LOW, "Could not find any previously written ARG files. "
                 "Try disabling resume\n");
        return false;
    }
    config.arg_file = arg_file;

    string sites_file = get_out_sites_file(config, config.resume_iter);
    Sites test_sites;
    if (read_sites(sites_file.c_str(), &test_sites, -1, -1, true)) {
        printLog(LOG_LOW, "Detected phased output sites file. Using %s as input"
                 " and assuming data is unphased\n", sites_file.c_str());
        config.sites_file = sites_file;
        config.unphased=1;
        config.vcf_file = "";
        config.vcf_list_file = "";
        config.subsites_file="";
        config.fasta_file="";

        string maskfile = config.out_prefix + config.mcmcmc_prefix + ".masked_regions.bed";
        CompressStream stream(maskfile.c_str(), "r");
        if (stream.stream != NULL) {
            printLog(LOG_LOW, "Using mask file %s output from prior run\n",
                     maskfile.c_str());
            config.maskmap=maskfile;
            config.maskN=-1;
            config.mask_cluster="";
            config.ind_maskmap="";
            config.mask_uncertain=0.0;
            config.expand_ind_mask=true;
        }
    }

    printLog(LOG_LOW, "resuming at stage=%s, iter=%d, arg=%s\n",
             config.resume_stage.c_str(), config.resume_iter,
             config.arg_file.c_str());

    // clean up
    fclose(stats_file);

    return true;
}


bool check_overwrite(Config &config)
{
    // check for stats file
    string stats_filename = config.out_prefix + config.mcmcmc_prefix
        + STATS_SUFFIX;
    bool exists = !access(stats_filename.c_str(), F_OK);
    if (config.resume || config.overwrite || !exists) {
        return true;
    } else {
        printError("Stats file already exists: %s", stats_filename.c_str());
        printError("To force overwrite use --overwrite option");
        return false;
    }
}


bool ensure_output_dir(const char *outdir)
{
    char *path = strdup(outdir);
    char *dir = dirname(path);
    bool result = true;
    if (!makedirs(dir)) {
        printError("could not make directory for output files '%s'", dir);
        result = false;
    }
    free(path);
    return result;
}


//=============================================================================


int main(int argc, char **argv)
{

#ifdef ARGWEAVER_MPI
    MPI::Init(argc, argv);
#endif

    // parse command line arguments
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
        return ret;

    // ensure output dir
    if (!ensure_output_dir(c.out_prefix.c_str()))
        return EXIT_ERROR;

    // check overwriting
    if (!check_overwrite(c))
        return EXIT_ERROR;

    // setup logging
    set_up_logging(c, c.verbose, (c.resume ? "a" : "w"));

    // try to resume a previous run
    if (!setup_resume(c)) {
        printError("resume failed.");
        if (c.overwrite) {
            c.resume = false;
            printLog(LOG_LOW, "Resume failed.  Sampling will start from scratch"
                     " since overwrite is enabled.\n");
        } else {
            return EXIT_ERROR;
        }
    }

#ifdef ARGWEAVER_MPI
    if (c.mpi) {
        int numcore = MPI::COMM_WORLD.Get_size();
        if (numcore % c.mcmcmc_numgroup != 0) {
            fprintf(stderr, "Error: number of cores should be evenly divisible"
                    " by number of mcmcmc threads");
        }
        int groupsize = numcore / c.mcmcmc_numgroup;
        int sites_num = MPI::COMM_WORLD.Get_rank() % groupsize;
        char tmp[10000];
        sprintf(tmp, "%s%i.sites", c.sites_file.c_str(),
                sites_num);
        c.sites_file = (string)tmp;
        sprintf(tmp, "%s%i", c.out_prefix.c_str(),
                sites_num);
        c.out_prefix = (string)tmp;
        /*        if (c.cr_file != "") {
            sprintf(tmp, "%s%i.cr.gz", c.cr_file.c_str(), sites_num);
            c.cr_file = (string)tmp;
            }*/
	if (c.arg_file != "") {
	    sprintf(tmp, "%s%i.smc.gz", c.arg_file.c_str(), sites_num);
	    c.arg_file = (string)tmp;
	}
    }
#endif



    // log intro
    if (c.resume)
        printLog(LOG_LOW, "RESUME\n");
    log_intro(LOG_LOW);
    log_prog_commands(LOG_LOW, argc, argv);
    Timer timer;


    // init random number generator
    if (c.randseed == 0)
        c.randseed = time(NULL);
#ifdef ARGWEAVER_MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
        for (int i=1; i < MPI::COMM_WORLD.Get_size(); i++) {
            int seed = irand(12581020);
            MPI::COMM_WORLD.Send(&seed, 1, MPI::INT, i, 13);
        }
    } else {
        MPI::COMM_WORLD.Recv(&c.randseed, 1, MPI::INT, 0, 13);
    }
#endif
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n", c.randseed);

    // read sequences
    Sites sites;
    Sequences sequences;
    SitesMapping *sites_mapping = NULL;
    unique_ptr<SitesMapping> sites_mapping_ptr;
    Region seq_region;
    Region seq_region_compress;

    set<string> keep_inds;

    if (c.subsites_file != "") {
        FILE *infile;
        if (!(infile = fopen(c.subsites_file.c_str(), "r"))) {
            printError("Could not open subsites file %s",
                       c.subsites_file.c_str());
            return false;
        }
        char name[10000];
        while (EOF != fscanf(infile, "%s", name))
            keep_inds.insert(string(name));
        fclose(infile);
    }

    if (c.fasta_file != "") {
        // read FASTA file

        if (!read_fasta(c.fasta_file.c_str(), &sequences)) {
            printError("could not read fasta file");
            return EXIT_ERROR;
        }
        seq_region.set("chr", 0, sequences.length());

        printLog(LOG_LOW, "read input sequences (nseqs=%d, length=%d)\n",
                 sequences.get_num_seqs(), sequences.length());

        make_sites_from_sequences(&sequences, &sites);
    } else if (c.sites_file != "") {
        // read sites file

        // parse subregion if given
        int subregion[2] = {-1, -1};
        if (c.subregion_str != "") {
            if (!parse_region(c.subregion_str.c_str(),
                              &subregion[0], &subregion[1])) {
                printError("subregion is not specified as 'start-end'");
                return EXIT_ERROR;
            }
            subregion[0] -= 1; // convert to 0-index
        }

        // read sites
        CompressStream stream(c.sites_file.c_str());
        if (!stream.stream ||
            !read_sites(stream.stream, &sites, subregion[0], subregion[1])) {
            printError("could not read sites file");
            return EXIT_ERROR;
        }
        stream.close();

        printLog(LOG_LOW, "read input sites (chrom=%s, start=%d, end=%d, "
                 "length=%d, nseqs=%d, nsites=%d)\n",
                 sites.chrom.c_str(), sites.start_coord, sites.end_coord,
                 sites.length(), sites.get_num_seqs(),
                 sites.get_num_sites());

        // sanity check for sites
        if (sites.get_num_sites() == 0) {
            printLog(LOG_LOW, "no sites given.  terminating.\n");
            return EXIT_ERROR;
        }
        seq_region.set(sites.chrom, sites.start_coord, sites.end_coord);
    } else if (c.vcf_file != "") {
        if (c.vcf_list_file != "") {
            printLog(LOG_LOW, "Cannot use both --vcf-file and --vcf-list-file. Terminating\n");
            return EXIT_ERROR;
        }
        if (!read_vcf(c.vcf_file, &sites, c.subregion_str,
                      c.vcf_min_qual, c.vcf_filter, c.use_genotype_probs,
                      c.mask_uncertain, false, c.tabix_dir, keep_inds)) {
            printError("Could not read VCF file");
            return EXIT_ERROR;
        }
        printLog(LOG_LOW, "read input sites from VCF (chrom=%s, start=%d, end=%d, length=%d, nseqs=%d, nsites=%d)\n",
                 sites.chrom.c_str(), sites.start_coord, sites.end_coord,
                 sites.length(), sites.get_num_seqs(),
                 sites.get_num_sites());
        if (sites.get_num_sites() == 0) {
            printLog(LOG_LOW, "no sites given.  terminating.\n");
            return EXIT_ERROR;
        }
        seq_region.set(sites.chrom, sites.start_coord, sites.end_coord);
    } else if (c.vcf_list_file != "") {
        vector<string> vcf_files;
        FILE *infile;
        if (!(infile= fopen(c.vcf_list_file.c_str(), "r"))) {
            printError("Could not open vcf_list_file %s\n", c.vcf_list_file.c_str());
            return false;
        }
        char *tmpstr;
        while (NULL != (tmpstr = fgetline(infile))) {
            tmpstr = trim(tmpstr);
            if (strlen(tmpstr) > 0)
                vcf_files.push_back(string(tmpstr));
            delete [] tmpstr;

        }
        if (!read_vcfs(vcf_files, &sites, c.subregion_str,
                       c.vcf_min_qual, c.vcf_filter, c.use_genotype_probs,
                       c.mask_uncertain, c.tabix_dir, keep_inds)) {
            printError("Error reading VCF files\n");
            return EXIT_ERROR;
        }
        seq_region.set(sites.chrom, sites.start_coord, sites.end_coord);
    } else {
        // no input sequence specified
        printError("must specify sequences (use --fasta or --sites)");
        return EXIT_ERROR;
    }

    if (c.rename_file != "")
        sites.rename(c.rename_file);

    if (keep_inds.size() > 0) {
        if (sites.subset(keep_inds)) {
            printError("Error subsetting sites\n");
            return false;
        }
    }

    //read in masks
    vector<TrackNullValue> ind_maskmap;
    for (int i=0; i < sites.get_num_seqs(); i++)
        ind_maskmap.push_back(TrackNullValue());
    if (c.ind_maskmap != "") {
        FILE *infile = fopen(c.ind_maskmap.c_str(), "r");
        if (infile == NULL) {
            fprintf(stderr, "Error opening ind_maskmap %s\n", c.ind_maskmap.c_str());
            return EXIT_ERROR;
        }
        char ind[10000], maskfile[100000];
        while (EOF != fscanf(infile, "%s %s", ind, maskfile)) {
            CompressStream stream(maskfile, "r");
            TrackNullValue indmask;
            if (!stream.stream ||
                !read_track_filter(stream.stream, &indmask, seq_region)) {
                printError("cannot read mask %s for ind %s\n", maskfile, ind);
                return EXIT_ERROR;
            }
            apply_mask_sites(&sites, indmask, ind, &ind_maskmap);
        }
        fclose(infile);
    }

    TrackNullValue maskmap_orig;
    TrackNullValue maskmap = sites.remove_masked();
    if (c.maskmap != "") {
        //read mask
        printLog(LOG_LOW, "Reading %s\n", c.maskmap.c_str());
        TrackNullValue curr_maskmap;
        CompressStream stream(c.maskmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &curr_maskmap, seq_region)) {
            printError("cannot read mask map '%s'",
                       c.maskmap.c_str());
            return EXIT_ERROR;
        }
        maskmap.merge_tracks(curr_maskmap);
    }
    c.all_masked=false;

    if (c.maskN >= 0) {
        TrackNullValue maskNmap = get_n_regions(sites, c.maskN);
        maskmap.merge_tracks(maskNmap);
    }

    if (c.mask_cluster != "") {
        int numsnp, window;
        if (2 != sscanf(c.mask_cluster.c_str(), "%i,%i", &numsnp, &window)) {
            printError("Bad format in mask_cluster agument; expect numsnp,windowSize");
            return EXIT_ERROR;
        }
        TrackNullValue cluster_mask = get_snp_clusters(sites, numsnp, window);
        maskmap.merge_tracks(cluster_mask);
    }

    if (maskmap.size() > 0 && sites.get_num_sites() > 0) {
        int old_num_sites = sites.get_num_sites();
        // this is done so that compression is not affected by missing data
        // sites (because compression does not actually reduce the number of
        // sites and will fail if too many sites exist). After compression,
        // sites are converted to sequences and a compressed mask applied
        sites.remove_overlapping(maskmap);
        int new_num_sites = sites.get_num_sites();
        printLog(LOG_LOW, "Removed %i sites overlapping mask (old=%i , new=%i)\n",
                 old_num_sites - new_num_sites, old_num_sites, new_num_sites);
    }

    // before compression unmask everything; do not want masked sites to
    // prevent compression- reapply masks after compression
    if (sites.get_num_sites() > 0) {
        unmask_inds(&sites, &ind_maskmap);
        int nremove = sites.remove_invariant();
        printLog(LOG_LOW, "%i sites are partially masked but otherwise invariant\n",
                 nremove);
    }

    // compress sequences
    // first remove any sites that fall under mask

    sites_mapping = new SitesMapping();
    sites_mapping_ptr = unique_ptr<SitesMapping>(sites_mapping);

    if (!find_compress_cols(&sites, c.compress_seq, sites_mapping)) {
        printError("unable to compress sequences at given compression level"
                   " (--compress-seq)");
        return EXIT_ERROR;
    }
    compress_sites(&sites, sites_mapping);
    make_sequences_from_sites(&sites, &sequences);
    seq_region_compress.set(seq_region.chrom, 0, sequences.length());

    // compress mask
    maskmap_orig = maskmap;
    if (!c.resume)
        maskmap.write_track_regions(c.out_prefix + c.mcmcmc_prefix + ".masked_regions.bed");
    if (maskmap.size() > 0) {
        // apply mask
        if (sites_mapping)
            compress_mask(maskmap, sites_mapping);
        apply_mask_sequences(&sequences, maskmap);
    }
    if (ind_maskmap.size() > 0) {
        for (int i=0; i < (int)ind_maskmap.size(); i++) {
            compress_mask(ind_maskmap[i], sites_mapping, c.expand_ind_mask);
            apply_mask_sequences(&sequences, ind_maskmap[i],
                                 sites.names[i].c_str());
        }
    }

    // report number of masked sites
    bool *masked = new bool [sequences.length()];
    find_masked_sites(sequences.get_seqs(), sequences.get_num_seqs(),
                      sequences.length(), masked);
    int nmasked = 0;
    for (int i=0; i<sequences.length(); i++)
        nmasked += int(masked[i]);
    delete [] masked;
    printLog(LOG_LOW, "masked %d (%.1f%%) positions\n", nmasked,
             100.0 * nmasked / double(sequences.length()));
    if (nmasked == (int)sequences.length())
        c.all_masked=true;

    // setup phasing options
    if (c.unphased_file != "") {
        if (c.vcf_file != "" || c.vcf_list_file != "")
            printError("Cannot use --unphased-file with VCF input");
        c.model.unphased_file = c.unphased_file;
        c.unphased = true;
    }
    if (c.unphased || c.vcf_file != "" || c.vcf_list_file != "")
        c.model.unphased = true;
    if (c.model.unphased)
        sequences.set_pairs(&c.model);
    if (c.randomize_phase > 0.0) {
        sequences.randomize_phase(c.randomize_phase);
    }
    if (c.switch_error_rate > 0.0)
        sequences.add_switch_errors(c.switch_error_rate);
    if (c.no_sample_phase)
        c.sample_phase_step=0;
    else if (c.sample_phase_step == 0 && c.model.unphased)
        c.sample_phase_step = c.sample_step;

    if (c.write_sites || c.write_sites_only) {
        log_sequences(sites.chrom, &sequences, &c, sites_mapping, -1);
        printLog(LOG_LOW, "Wrote sites\n");
        if (c.write_sites_only) return(0);
    }


    // setup model parameters
    if (c.times_file != "")
        c.model.set_times_from_file(c.times_file);
    else if (c.time_step)
        c.model.set_linear_times(c.time_step, c.ntimes);
    else
        c.model.set_log_times(c.maxtime, c.ntimes, c.delta);
    if (c.pop_tree_file != "") {
        c.model.read_population_tree(c.pop_tree_file);
        c.model.pop_tree->max_migrations = ( c.start_mig_iter <= 0 ? 1 : 0 );
        if (c.pop_file != "") {
            sequences.set_pops_from_file(c.pop_file);
        }
        if (c.model.num_pops() != 1 && !c.model.smc_prime) {
            printLog(LOG_LOW, "Warning: Using original SMC model. --smc-prime option "
                     "is recommended for multiple population models\n");
        }
    }

    c.model.rho = c.rho;
    c.model.mu = c.mu;
    if (c.popsize_file != "") {
        // use population sizes from a file
        c.model.read_population_sizes(c.popsize_file);
    } else {
        // use single constant population size
        c.model.set_popsizes(c.popsize_str);
    }
    const double infsites_penalty = 1e-100; // TODO: make configurable
    if (c.infsites)
        c.model.infsites_penalty = infsites_penalty;

    if (c.age_file != "")
	sequences.set_age(c.age_file, c.model.ntimes, c.model.times);
    else sequences.set_age();

    // setup phasing options
    if (c.unphased_file != "") {
        if (c.vcf_file != "" || c.vcf_list_file != "")
            printError("Cannot use --unphased-file with VCF input");
        c.model.unphased_file = c.unphased_file;
    }
    if (c.unphased || (!c.phased_vcf && (c.vcf_file != "" || c.vcf_list_file != "")))
        c.model.unphased = true;
    if (c.model.unphased)
        sequences.set_pairs(&c.model);
    if (c.randomize_phase > 0.0) {
        sequences.randomize_phase(c.randomize_phase);
    }
    if (c.no_sample_phase)
        c.sample_phase_step=0;
    else if (c.sample_phase_step == 0)
        c.sample_phase_step = c.sample_step;

    // setup read-backed phase blocks
    if (c.phase_blocks_file != "") {
        bool mcmcmc = false;
#ifdef ARGWEAVER_MPI
        mcmcmc = (c.mcmcmc_numgroup > 1);
#endif
        if (c.gibbs || c.no_sample_arg || mcmcmc || c.infsites ||
            sequences.base_probs.size() > 0) {
            printError("read-backed phase blocks cannot yet be used with"
                       " --gibbs, --no-sample-arg, --mcmcmc, --infsites or"
                       " genotype probabilities");
            return EXIT_ERROR;
        }
        if (c.phase_point_prob < 0.0)
            c.phase_point_prob = c.phase_switch_prob;
        if (c.phase_point_prob > 0.5) {
            printError("--phase-point-error-prob must be at most 0.5");
            return EXIT_ERROR;
        }
        // threading treats the current phase as fixed (this overrides
        // --unphased); phase is only changed by the phase block moves
        c.model.unphased = false;
        if (!read_phase_blocks(c.phase_blocks_file.c_str(), &sequences,
                               sites_mapping, c.phase_blocks))
            return EXIT_ERROR;
        // record the diploid pairs given by the phase blocks
        sequences.pairs.assign(sequences.get_num_seqs(), -1);
        for (unsigned int j=0; j<c.phase_blocks.size(); j++) {
            sequences.pairs[c.phase_blocks[j].hap1] = c.phase_blocks[j].hap2;
            sequences.pairs[c.phase_blocks[j].hap2] = c.phase_blocks[j].hap1;
        }
        printLog(LOG_LOW, "read-backed phase blocks for %d individuals,"
                 " switch probability %g, point error probability %g\n",
                 (int) c.phase_blocks.size(), c.phase_switch_prob,
                 c.phase_point_prob);
    }

    if (c.sample_popsize_num > 0) {
	if (c.popsize_em) {
	    printError("Error: cannot use --popsize-em with --sample-popsize\n");
	    return 1;
	}
        if (c.popsize_config_file != "") {
            c.model.popsize_config =
                PopsizeConfig(c.popsize_config_file, c.model.ntimes,
                              c.model.num_pops(), c.model.popsizes);
        } else if (c.popsize_config == 0) {
            c.model.popsize_config =
                PopsizeConfig(c.ntimes, c.model.num_pops(), true, true);
        } else if (c.popsize_config == 1) {
            c.model.set_popsize_config_by_pop_tree();
        } else if (c.popsize_config == 2) {
            c.model.popsize_config =
                PopsizeConfig(c.ntimes, c.model.num_pops(), true, false);
        } else {
            printError("Error: invalid value for --popsize-config\n");
            return 1;
        }
        c.model.popsize_config.numsample = c.sample_popsize_num;
        c.model.popsize_config.neighbor_prior = c.popsize_prior_neighbor;
	c.model.popsize_config.pseudocount = c.pseudocount;
    }
#ifdef ARGWEAVER_MPI
    c.model.mc3 = Mc3Config(c.mcmcmc_group, c.mcmcmc_heat);

    printf("MPI rank=%i size=%i\n",
           MPI::COMM_WORLD.Get_rank(),
           MPI::COMM_WORLD.Get_size());
    MPI::COMM_WORLD.Barrier();
    #endif
    if (c.init_popsize_random)
        c.model.set_popsizes_random();

    // read model parameter maps if given
    if (c.mutmap != "") {
        CompressStream stream(c.mutmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model.mutmap, seq_region)) {
            printError("cannot read mutation rate map '%s'", c.mutmap.c_str());
            return EXIT_ERROR;
        }
    }
    if (c.recombmap != "") {
        CompressStream stream(c.recombmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model.recombmap, seq_region)) {
            printError("cannot read recombination rate map '%s'",
                       c.recombmap.c_str());
            return EXIT_ERROR;
        }
    }


    // make compressed model
    ArgModel model(c.model);
    if (!model.setup_maps(seq_region.chrom, seq_region.start, seq_region.end))
        return EXIT_ERROR;
    compress_model(&model, sites_mapping, c.compress_seq);

    // log original model
    model.log_model();

    // setup init ARG
    LocalTrees *trees = NULL;
    unique_ptr<LocalTrees> trees_ptr;
    if (c.arg_file != "") { // || c.cr_file != "") {
        // init ARG from file

        trees = new LocalTrees();
        trees_ptr = unique_ptr<LocalTrees>(trees);
        vector<string> seqnames;
        /*        if (c.cr_file != "") {
            if (c.arg_file != "") {
                printError("cannot use --arg-file and --cr");
                return EXIT_ERROR;
            }
            if (!read_init_arg_cr(c.cr_file.c_str(), &c.model, trees, seqnames)) {
                printError("Could not read ARG");
                return EXIT_ERROR;
            }
            } else { */
            if (!read_init_arg(c.arg_file.c_str(), &c.model, trees, seqnames)) {
                printError("could not read ARG");
                return EXIT_ERROR;
            }
            //        }

        if (!trees->set_seqids(seqnames, sequences.names)) {
            printError("input ARG's sequence names do not match input"
                       " sequences");
            return EXIT_ERROR;
        }

        printLog(LOG_LOW, "read input ARG (chrom=%s, start=%d, end=%d,"
                 " nseqs=%d)\n",
                 trees->chrom.c_str(), trees->start_coord, trees->end_coord,
                 trees->get_num_leaves());

        // check ARG matches sites/sequences
        if (trees->start_coord != seq_region.start ||
            trees->end_coord != seq_region.end) {
            printError("trees range does not match sites: tree(start=%d,"
                       " end=%d), sites(start=%d, end=%d) [compressed"
                       " coordinates]",
                       trees->start_coord, trees->end_coord,
                       seq_region.start, seq_region.end);
            return EXIT_ERROR;
        }

        // compress input ARG if compression is requested
        if (sites_mapping)
            compress_local_trees(trees, sites_mapping, true);

    } else {
        // create new init ARG
        trees = new LocalTrees(seq_region_compress.start,
                               seq_region_compress.end);
        trees->chrom = seq_region.chrom;
        trees_ptr = unique_ptr<LocalTrees>(trees);
    }


    // check for region sample
    if (c.resample_region_str != "") {
        if (!parse_region(c.resample_region_str.c_str(),
                          &c.resample_region[0], &c.resample_region[1]))
            {
                printError("--resample-region is not specified as 'start-end'");
                return EXIT_ERROR;
            }
        c.resample_region[0] -= 1; // convert to 0-index

        if (sites_mapping) {
            c.resample_region[0] = sites_mapping->compress(c.resample_region[0]);
            c.resample_region[1] = sites_mapping->compress(c.resample_region[1]-1)+1;
        }
    }


    // init stats file
    string stats_filename = c.out_prefix + c.mcmcmc_prefix + STATS_SUFFIX;
    const char *stats_mode = (c.resume ? "a" : "w");
    if (!(c.stats_file = fopen(stats_filename.c_str(), stats_mode))) {
        printError("could not open stats file '%s'", stats_filename.c_str());
        return EXIT_ERROR;
    }

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);

    // sample ARG
    printLog(LOG_LOW, "\n");
    sample_arg(&model, &sequences, trees, sites_mapping, &c, &maskmap_orig);

    // final log message
    maxrss = get_max_memory_usage() / 1000.0;
    printTimerLog(timer, LOG_LOW, "sampling time: ");
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);
    printLog(LOG_LOW, "FINISH\n");

    // clean up
    fclose(c.stats_file);

#ifdef ARGWEAVER_MPI
    MPI_Finalize();
#endif

    return 0;
}
