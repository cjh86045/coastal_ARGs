import shutil
from ast import literal_eval


from utils import setup_logging, run_command, get_sample_names_list, process_samples_in_parallel
from configure import *


class BinAnnotation:
    def __init__(self, input_folder):
        # Input directories from upstream pipelines
        self.acc_output_folder = Path(input_folder, "preprocessing", "acc_output")
        self.binning_output_folder = Path(input_folder, "binning")
        self.dastool_output_folder = Path(self.binning_output_folder, "dastool_out")
        self.checkm_output_folder = Path(self.binning_output_folder, "checkm_out")

        # Output root directory
        self.bin_annotation_folder = Path(input_folder, "bin_annotation")
        self.bin_annotation_folder.mkdir(parents=True, exist_ok=True)

        # Quality bins (one TSV per sample)
        self.quality_bins_output_folder = Path(self.bin_annotation_folder, "quality_bins_results")
        # Count of quality bins per sample
        self.bin_count_output_folder = Path(self.bin_annotation_folder, "bin_count_sample_results")
        # High-quality bins FASTA (one subdir per sample, named <sample>_quality_bins)
        self.quality_bins_fasta_folder = Path(self.bin_annotation_folder, "quality_bins_fasta")
        # GTDB-tk output (single run on merged bins)
        self.gtdb_output_folder = Path(self.bin_annotation_folder, "gtdb_results")
        # Final merged results (per sample + global)
        self.final_results_folder = Path(self.bin_annotation_folder, "final_results_with_species")
        self.sample_merge_output_folder = Path(self.bin_annotation_folder, "sample_merge_results")

        # Directory holding ALL quality bins merged together (with sample prefix)
        self.merged_bins_folder = Path(self.bin_annotation_folder, "merged_quality_bins")

        # Configure logging
        setup_logging(self.bin_annotation_folder, "bin_annotation.log")
        self.logger = logging.getLogger("BinAnnotation")

        # Sample names list
        self.sample_names_list = get_sample_names_list(self.dastool_output_folder, "species")

    # ------------------------------------------------------------------ #
    # CheckM parsing
    # ------------------------------------------------------------------ #
    def get_checkm_results(self, sample):
        """
        Parse CheckM bin_stats_ext.tsv for a sample.
        The stats are stored as a JSON string in the second column.
        """
        checkm_results_file = self.checkm_output_folder / f"{sample}_checkm" / "storage" / "bin_stats_ext.tsv"
        if not checkm_results_file.exists():
            self.logger.info(f"No checkm results file found, skipping {sample}")
            return pd.DataFrame()

        df = pd.read_csv(str(checkm_results_file), sep='\t', header=None, names=['bin_name', 'stats_json'])
        df['stats'] = df['stats_json'].apply(lambda x: literal_eval(x) if isinstance(x, str) else {})

        result = pd.DataFrame({
            'bin_name': df['bin_name'],
            'Completeness': df['stats'].apply(lambda x: x.get('Completeness', '')),
            'Contamination': df['stats'].apply(lambda x: x.get('Contamination', '')),
            'Genome size': df['stats'].apply(lambda x: x.get('Genome size', '')),
            'scaffolds': df['stats'].apply(lambda x: x.get('# scaffolds', '')),
            'contigs': df['stats'].apply(lambda x: x.get('# contigs', '')),
            'Longest scaffold': df['stats'].apply(lambda x: x.get('Longest scaffold', '')),
            'Longest contig': df['stats'].apply(lambda x: x.get('Longest contig', '')),
            'N50 (scaffolds)': df['stats'].apply(lambda x: x.get('N50 (scaffolds)', '')),
            'N50 (contigs)': df['stats'].apply(lambda x: x.get('N50 (contigs)', ''))
        })

        return result

    # ------------------------------------------------------------------ #
    # Step 1: identify high-quality bins
    # ------------------------------------------------------------------ #
    def get_quality_bins(self):
        """Identify high-quality bins based on CheckM results (Completeness > 50%, Contamination < 10%)."""
        self.logger.info("Identifying quality bins based on CheckM results")
        self.quality_bins_output_folder.mkdir(parents=True, exist_ok=True)
        self.bin_count_output_folder.mkdir(parents=True, exist_ok=True)

        bin_count = []
        with tqdm(total=len(self.sample_names_list), desc="Quality bins", unit="sample") as pbar:
            for sample_name in self.sample_names_list:
                das_tool_contig2bin_file = self.dastool_output_folder / f"{sample_name}_dastool" / "output_DASTool_contig2bin.tsv"
                if not das_tool_contig2bin_file.exists():
                    self.logger.info(f"No DASTool contig2bin file found, skipping {sample_name}")
                    bin_count.append(0)
                    pbar.update(1)
                    continue

                df = pd.read_table(str(das_tool_contig2bin_file), header=None, names=['contig_id', 'bin_name'])

                checkm_result_df = self.get_checkm_results(sample_name)
                if checkm_result_df.empty:
                    self.logger.info(f"No checkm result found for {sample_name}, skipping")
                    bin_count.append(0)
                    pbar.update(1)
                    continue

                # Merge bin info with CheckM results
                result = pd.merge(df[['bin_name']].drop_duplicates(),
                                  checkm_result_df,
                                  on="bin_name",
                                  how="left")

                # Filter high-quality bins
                quality_bins = result[(result['Completeness'] > 50) &
                                      (result['Contamination'] < 10)].copy()

                bin_count.append(len(quality_bins))

                quality_bins.to_csv(
                    self.quality_bins_output_folder / f"{sample_name}_quality_bins_results.tsv",
                    sep="\t", index=False
                )
                pbar.update(1)

        # Save quality bin count per sample
        bin_count_df = pd.DataFrame({
            "sample": self.sample_names_list,
            "quality_bin_count": bin_count
        })
        bin_count_df.to_csv(
            self.bin_count_output_folder / "bin_count_sample_results.tsv",
            sep="\t", index=False
        )

    # ------------------------------------------------------------------ #
    # Step 2: extract high-quality bin FASTA files
    # ------------------------------------------------------------------ #
    def extract_quality_bins(self):
        """Extract high-quality bin FASTA files for GTDB-tk analysis."""
        self.logger.info("Extracting quality bins for GTDB-tk analysis")
        self.quality_bins_fasta_folder.mkdir(parents=True, exist_ok=True)

        with tqdm(total=len(self.sample_names_list), desc="Extract bins", unit="sample") as pbar:
            for sample_name in self.sample_names_list:
                bins_file = self.dastool_output_folder / f"{sample_name}_dastool" / "output_DASTool_bins"
                output_dir = self.quality_bins_fasta_folder / f"{sample_name}_quality_bins"
                output_dir.mkdir(parents=True, exist_ok=True)

                quality_bins_result = self.quality_bins_output_folder / f"{sample_name}_quality_bins_results.tsv"
                if not quality_bins_result.exists():
                    self.logger.warning(f"No quality bins result file found, skipping {sample_name}")
                    pbar.update(1)
                    continue

                df = pd.read_csv(str(quality_bins_result), sep='\t')
                bin_names = df['bin_name'].unique()
                for bin_name in bin_names:
                    src_file = bins_file / f"{bin_name}.fa"
                    if src_file.exists():
                        shutil.copy(str(src_file), str(output_dir))
                    else:
                        self.logger.warning(f"Bin file {src_file} not found")

                pbar.update(1)

    # ------------------------------------------------------------------ #
    # Step 3: prepare merged bins (with sample prefix)
    # ------------------------------------------------------------------ #
    def _prepare_merged_bins(self):
        """
        Collect all high-quality bins from all samples into a single directory.
        Each bin is renamed with the sample name as prefix to avoid name collisions
        and to preserve sample-of-origin information.

        Naming rule:  {sample_name}__{original_bin_name}.fa
        Returns the merged bins directory.
        """
        self.merged_bins_folder.mkdir(parents=True, exist_ok=True)

        # If merging was already done and is non-empty, skip
        existing = list(self.merged_bins_folder.glob("*.fa"))
        if existing:
            self.logger.info(f"Merged bins directory already contains {len(existing)} bins. Skipping copy.")
            return self.merged_bins_folder

        self.logger.info("Merging high-quality bins from all samples (adding sample prefix)...")
        with tqdm(total=len(self.sample_names_list), desc="Merge bins", unit="sample") as pbar:
            for sample_name in self.sample_names_list:
                sample_bins_dir = self.quality_bins_fasta_folder / f"{sample_name}_quality_bins"
                if not sample_bins_dir.exists():
                    self.logger.warning(f"No quality bins folder for {sample_name}. Skipping.")
                    pbar.update(1)
                    continue

                for bin_fa in sample_bins_dir.glob("*.fa"):
                    # e.g. sampleA__bin.1.fa
                    new_name = f"{sample_name}__{bin_fa.stem}.fa"
                    shutil.copy(str(bin_fa), str(self.merged_bins_folder / new_name))

                pbar.update(1)

        total = len(list(self.merged_bins_folder.glob("*.fa")))
        self.logger.info(f"Total merged bins: {total}")
        return self.merged_bins_folder

    # ------------------------------------------------------------------ #
    # Step 3 (cont.): run GTDB-tk ONCE on merged bins
    # ------------------------------------------------------------------ #
    def run_gtdbtk(self):
        """
        Run GTDB-tk ONCE on all high-quality bins merged from every sample.
        Bin names are prefixed with sample name (see _prepare_merged_bins),
        so GTDB-tk output can be traced back to the original sample.
        """
        self.gtdb_output_folder.mkdir(parents=True, exist_ok=True)

        # Output directory for the single GTDB-tk run
        output_dir = self.gtdb_output_folder / "gtdbtk_all"

        # Skip if final summary already exists
        summary_file = output_dir / "gtdb_.bac120.summary.tsv"
        if summary_file.exists() and summary_file.stat().st_size > 0:
            self.logger.info("GTDB-tk results already exist. Skipping.")
            return

        # Step 1: prepare merged bins (with sample prefix)
        merged_dir = self._prepare_merged_bins()
        if not list(merged_dir.glob("*.fa")):
            self.logger.error("No merged bins available for GTDB-tk. Aborting.")
            return

        # Step 2: run GTDB-tk a single time on the whole collection
        gtdb_command = [
            "gtdbtk", "classify_wf",
            "--genome_dir", str(merged_dir),
            "--extension", "fa",
            "--prefix", "gtdb_",
            "--out_dir", str(output_dir),
            "--cpus", str(THREAD_CONFIG["gtdbtk"]),
        ]

        self.logger.info(f"Running GTDB-tk on ALL merged bins (single run): {' '.join(gtdb_command)}")
        run_command(gtdb_command, logger=self.logger, log_message=' '.join(gtdb_command))

        if not summary_file.exists():
            self.logger.error("GTDB-tk finished but summary file was not found. Check logs.")

    # ------------------------------------------------------------------ #
    # ARG-containing bin identification
    # ------------------------------------------------------------------ #
    def get_arg_containing_bins(self, sample_name):
        """
        Identify bins that contain at least one ARG-carrying contig,
        by intersecting ACC contig IDs with the contig2bin mapping.
        """
        arg_bins = set()

        acc_file = self.acc_output_folder / f"{sample_name}_acc_contigs.fa"
        if not acc_file.exists():
            self.logger.debug(f"No ACC file found for {sample_name}")
            return arg_bins

        # ARG-carrying contig IDs
        arg_contigs = {rec.id for rec in SeqIO.parse(acc_file, "fasta")}

        contig2bin_file = self.dastool_output_folder / f"{sample_name}_dastool" / "output_DASTool_contig2bin.tsv"
        if not contig2bin_file.exists():
            self.logger.warning(f"No contig2bin file found for {sample_name}")
            return arg_bins

        contig2bin_df = pd.read_csv(str(contig2bin_file), sep='\t',
                                    header=None, names=['contig_id', 'bin_name'])
        arg_bins = set(contig2bin_df[contig2bin_df['contig_id'].isin(arg_contigs)]['bin_name'].unique())

        return arg_bins

    # ------------------------------------------------------------------ #
    # Step 4: merge quality + species + ARG info
    # ------------------------------------------------------------------ #
    def merge_results(self):
        """
        Merge per-sample quality bin info with GTDB-tk species classification.
        GTDB-tk was run once on merged (sample-prefixed) bins, so here we
        split the prefix back into sample_name to reconstruct per-sample results.
        """
        self.logger.info("Merging quality bins with species information")
        self.final_results_folder.mkdir(parents=True, exist_ok=True)
        self.sample_merge_output_folder.mkdir(parents=True, exist_ok=True)

        # Global GTDB-tk summary (single run)
        summary_file = self.gtdb_output_folder / "gtdbtk_all" / "gtdb_.bac120.summary.tsv"
        if not summary_file.exists():
            self.logger.error("GTDB-tk summary not found. Did you run run_gtdbtk()?")
            return

        species_info_df = pd.read_csv(str(summary_file), sep='\t')
        # 'user_genome' looks like "sampleA__bin.1" -> split into sample + bin
        species_info_df[['sample_name', 'bin_name']] = species_info_df['user_genome'].str.split('__', n=1, expand=True)

        sample_merge_rows = []

        with tqdm(total=len(self.sample_names_list), desc="Merge results", unit="sample") as pbar:
            for sample_name in self.sample_names_list:
                quality_bins_file = self.quality_bins_output_folder / f"{sample_name}_quality_bins_results.tsv"
                if not quality_bins_file.exists():
                    self.logger.warning(f"No quality bins file found, skipping {sample_name}")
                    pbar.update(1)
                    continue

                quality_bins_df = pd.read_csv(str(quality_bins_file), sep='\t')

                # Subset GTDB results for this sample (by prefix)
                sample_species = species_info_df[species_info_df['sample_name'] == sample_name].copy()

                final_df = pd.merge(quality_bins_df, sample_species, on='bin_name', how='left')

                # Mark ARG-containing bins
                arg_bins = self.get_arg_containing_bins(sample_name)
                final_df['contains_ARG'] = final_df['bin_name'].isin(arg_bins)

                # Save per-sample final result
                final_output = self.final_results_folder / f"{sample_name}_final_results_with_species.tsv"
                final_df.to_csv(str(final_output), sep='\t', index=False)

                # Collect ARG-containing rows for global summary
                arg_rows = final_df[final_df['contains_ARG']][['sample_name', 'bin_name', 'classification']].copy()
                if not arg_rows.empty:
                    sample_merge_rows.append(arg_rows)

                pbar.update(1)

        if sample_merge_rows:
            sample_merge_df = pd.concat(sample_merge_rows, ignore_index=True)
            sample_merge_df.to_csv(
                str(self.sample_merge_output_folder / "sample_merge_results.tsv"),
                sep="\t", index=False
            )
        else:
            self.logger.warning("No ARG-containing bins found in any sample; summary file not written.")

    # ------------------------------------------------------------------ #
    # Pipeline entry point
    # ------------------------------------------------------------------ #
    def run(self):
        """Run the entire bin annotation pipeline."""
        import time
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("Starting bin annotation pipeline...")
        self.logger.info("=" * 60)

        # Step 1: identify high-quality bins
        self.logger.info("Step 1/4: Identifying high-quality bins...")
        self.get_quality_bins()

        # Step 2: extract high-quality bin FASTA files
        self.logger.info("Step 2/4: Extracting high-quality bins...")
        self.extract_quality_bins()

        # Step 3: run GTDB-tk ONCE on all merged bins
        self.logger.info("Step 3/4: Running GTDB-tk (single run on merged bins)...")
        self.run_gtdbtk()

        # Step 4: merge quality + species + ARG info
        self.logger.info("Step 4/4: Merging results...")
        self.merge_results()

        # Runtime summary
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)

        self.logger.info("=" * 60)
        self.logger.info(f"Bin annotation pipeline completed!")
        self.logger.info(f"Total runtime: {int(hours)}h {int(minutes)}m {int(seconds)}s")
        self.logger.info("=" * 60)


if __name__ == '__main__':
    # Initialize and run the bin annotation pipeline
    annotator = BinAnnotation("")
    annotator.run()