from utils import setup_logging, run_command, get_sample_names_list,process_samples_in_parallel
from configure import *

class ARGsAbundanceProfile:
    def __init__(self, input_folder):
        # Preprocessing output directory
        self.input_folder = pathlib.Path(input_folder, "preprocessing")

        # Required inputs: clean reads and ARG-like ORFs
        self.fastp_output_folder = self.input_folder / "fastp_output"
        self.arg_like_orfs_folder = self.input_folder / "arg_like_orfs_output"

        # Abundance output root directory
        self.abundance_output_folder = pathlib.Path(input_folder, "args_abundance")

        # Phase 1: Estimate cell counts (ncells) by aligning clean reads to KO30 marker database
        self.count_cells_output_folder = pathlib.Path(self.abundance_output_folder, "count_cells_output")

        # Phase 2: Align clean reads to ARG-like ORFs database
        self.diamond_output_folder = pathlib.Path(self.abundance_output_folder, "arg_like_orfs_diamond_output")

        # Phase 3: Calculate and output ARGs abundance for each sample
        self.sample_abundance_output_folder = pathlib.Path(self.abundance_output_folder, "abundance_output")

        # Phase 3: Calculate and output ARGs abundance by antibiotic type
        self.type_abundance_output_folder = pathlib.Path(self.abundance_output_folder, "type_abundance_output")

        # Phase 4: Merge all sample results into a single summary file
        self.merge_results_output_folder = pathlib.Path(self.abundance_output_folder, "merge_results")

        # Create abundance output root directory
        self.abundance_output_folder.mkdir(parents=True, exist_ok=True)

        # Configure logging
        setup_logging(self.abundance_output_folder, "args_abundance_profile.log")
        self.logger = logging.getLogger("ARGsAbundanceProfile")

        # Get sample names list and store as instance variable
        self.sample_names_list = get_sample_names_list(self.fastp_output_folder, "megahit")
        print(self.sample_names_list)

        # Load CARD metadata for ARG annotation
        self.card_metadata = pd.read_table(CARD_STRUCTURE, sep="\t")

    def count_cells(self):
        """
        Phase 1: Align clean reads to KO30 marker database using DIAMOND BLASTx
        to estimate cell counts (ncells) for each sample.
        """
        self.count_cells_output_folder.mkdir(parents=True, exist_ok=True)

        # Check database exists once before processing all samples
        if not KO30_DB.exists():
            self.logger.error(f"KO30 database not found at {KO30_DB}. "
                              f"Please check the database path in configure.py.")
            return

        sample_names_list = list(self.sample_names_list)  # avoid iterating generator issues
        if not sample_names_list:
            self.logger.warning("No samples found for count_cells processing!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting count_cells (KO30 BLASTx)...")

        def process_sample(sample_name):
            """Process a single sample: align clean reads to KO30 marker database"""
            # Input files: clean paired-end reads
            clean_fwd = self.fastp_output_folder / f"clean_{sample_name}_1.fastq.gz"
            clean_rev = self.fastp_output_folder / f"clean_{sample_name}_2.fastq.gz"

            # Fix: use 'or' — skip if either read file is missing or empty
            if not clean_fwd.exists() or not clean_rev.exists():
                self.logger.warning(f"clean_reads not found for {sample_name}: "
                                    f"{clean_fwd.name} or {clean_rev.name}. Skipping.")
                return False

            if clean_fwd.stat().st_size == 0 or clean_rev.stat().st_size == 0:
                self.logger.warning(f"clean_reads are empty for {sample_name}. Skipping.")
                return False

            # Output file: DIAMOND alignment results
            diamond_output_file = self.count_cells_output_folder / f"{sample_name}_count_cells_results.tsv"

            # Skip if output already exists and is not empty
            if diamond_output_file.exists() and diamond_output_file.stat().st_size > 0:
                self.logger.info(f"count_cells output already exists for {sample_name}. Skipping.")
                return True

            # Build DIAMOND BLASTx command
            diamond_command = [
                "diamond", "blastx",
                "--db", str(KO30_DB),
                "--query", str(clean_fwd), str(clean_rev),
                "-o", str(diamond_output_file),
                "--evalue", "3",
                "--id", "45",
                "--query-cover", "0",
                "--max-hsps", "1",
                "--max-target-seqs", "1",
                "--threads", str(THREAD_CONFIG["diamond"]),
                "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
            ]

            # Run DIAMOND command
            try:
                self.logger.info(f"Running DIAMOND BLASTx for sample {sample_name}...")
                run_command(diamond_command, logger=self.logger,
                            log_message=f"Running DIAMOND BLASTx for sample {sample_name}: {' '.join(diamond_command)}")

                # Verify output file was generated successfully
                if diamond_output_file.exists() and diamond_output_file.stat().st_size > 0:
                    self.logger.info(f"count_cells completed successfully for {sample_name}")
                    return True
                else:
                    self.logger.error(f"count_cells failed to produce output for {sample_name}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running count_cells for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names_list), desc="count_cells (KO30)", unit="sample") as pbar:
            results = process_samples_in_parallel(
                self.logger, process_sample, sample_names_list,
                    "diamond", THREAD_CONFIG.get("diamond_pool", 2))
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")


        self.logger.info(f"count_cells (KO30 BLASTx) completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def run_diamond(self):
        """
        Align clean reads to ARG-like ORFs using DIAMOND BLASTx.
        Note: DIAMOND now accepts a .faa file directly as --db,
        so no separate makedb step is required.
        """
        self.diamond_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = list(self.sample_names_list)
        if not sample_names_list:
            self.logger.warning("No samples found for DIAMOND BLASTx alignment!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting DIAMOND BLASTx alignment to ARG-like ORFs...")

        def process_sample(sample_name):
            """Process a single sample: align clean reads to its ARG-like ORFs"""
            # Input files: clean paired-end reads
            clean_fwd = self.fastp_output_folder / f"clean_{sample_name}_1.fastq.gz"
            clean_rev = self.fastp_output_folder / f"clean_{sample_name}_2.fastq.gz"

            # Fix: skip if either read file is missing or empty
            if not clean_fwd.exists() or not clean_rev.exists():
                self.logger.warning(f"clean_reads not found for {sample_name}: "
                                    f"{clean_fwd.name} or {clean_rev.name}. Skipping.")
                return False

            if clean_fwd.stat().st_size == 0 or clean_rev.stat().st_size == 0:
                self.logger.warning(f"clean_reads are empty for {sample_name}. Skipping.")
                return False

            # Input file: ARG-like ORFs (used directly as DIAMOND database)
            arg_like_orfs_file = self.arg_like_orfs_folder / f"{sample_name}_arg_like_orfs.faa"
            if not arg_like_orfs_file.exists() or arg_like_orfs_file.stat().st_size == 0:
                self.logger.warning(f"arg_like_orfs file not found or empty for {sample_name}: "
                                    f"{arg_like_orfs_file}")
                return False

            # Output file: DIAMOND alignment results
            diamond_output_file = self.diamond_output_folder / f"{sample_name}_diamond_results.tsv"

            # Skip if output already exists and is not empty
            if diamond_output_file.exists() and diamond_output_file.stat().st_size > 0:
                self.logger.info(f"diamond blastx output already exists for {sample_name}. Skipping alignment.")
                return True

            # Build DIAMOND BLASTx command
            # Note: --db accepts the .faa file directly (no makedb needed)
            diamond_command = [
                "diamond", "blastx",
                "--db", str(arg_like_orfs_file),  # ARG-like ORFs as database (direct .faa input)
                "--query", str(clean_fwd), str(clean_rev),
                "-o", str(diamond_output_file),
                "--threads", str(THREAD_CONFIG["diamond"]),
                "--evalue", "1e-10",
                "--id", "90",
                "--query-cover", "80",
                "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore",
            ]

            # Run DIAMOND command
            try:
                self.logger.info(f"Running DIAMOND BLASTx for sample {sample_name}...")
                run_command(diamond_command, logger=self.logger,
                            log_message=f"Running DIAMOND BLASTx for sample {sample_name}: {' '.join(diamond_command)}")

                # Verify output file was generated successfully
                if diamond_output_file.exists() and diamond_output_file.stat().st_size > 0:
                    self.logger.info(f"DIAMOND BLASTx completed successfully for {sample_name}")
                    return True
                else:
                    self.logger.error(f"DIAMOND BLASTx failed to produce output for {sample_name}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running DIAMOND BLASTx for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0

        with tqdm(total=len(sample_names_list), desc="DIAMOND BLASTx", unit="sample") as pbar:
            results = process_samples_in_parallel(
                self.logger, process_sample, sample_names_list,
                    "diamond", THREAD_CONFIG.get("diamond_pool", 2))
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"DIAMOND BLASTx alignment completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def get_ncells(self, sample_name):
        """
        Estimate the number of cells (ncells) for a sample based on
        DIAMOND BLASTx alignment results against the KO30 marker database.

        The ncells value is calculated as the total coverage of KO30 marker
        genes divided by 30 (the number of universal single-copy marker genes).

        Args:
            sample_name (str): Name of the sample to process.

        Returns:
            float or None: Estimated cell count, or None if data is unavailable.
        """
        # Check if DIAMOND output exists and is not empty
        diamond_output_file = self.count_cells_output_folder / f"{sample_name}_count_cells_results.tsv"
        if not diamond_output_file.exists() or diamond_output_file.stat().st_size == 0:
            self.logger.warning(f"count_cells_results file not found or empty for sample {sample_name}: "
                                f"{diamond_output_file}")
            return None

        try:
            # Read DIAMOND BLASTx output
            df = pd.read_table(
                str(diamond_output_file),
                header=None,
                names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'evalue', 'bitscore']
            )

            # Merge with KO30 marker gene annotations
            ko30_df = pd.read_table(KO30_STRUCTURE, header=None, names=['sseqid', 'ko30'])
            df = pd.merge(df, ko30_df, on='sseqid', how='left')

            # For each read, keep only the best hit (lowest e-value)
            df = df.sort_values('evalue', ascending=True)
            df = df.drop_duplicates('qseqid', keep='first')

            if len(df) == 0:
                self.logger.warning(f'No marker-like sequences found in file <{sample_name}>.')
                return None

            # Calculate coverage per marker gene:
            # coverage = sum(length / slen) for each KO30 marker
            # Then divide by 30 (number of universal single-copy markers)
            ncells = df.groupby('ko30')[['length', 'slen']].apply(
                lambda x: sum(x['length'] / x['slen'])
            ).sum() / 30

            self.logger.info(f"Estimated ncells for {sample_name}: {ncells:.2f}")
            return ncells

        except Exception as e:
            self.logger.error(f"Error calculating ncells for {sample_name}: {str(e)}")
            return None

    def _calculate_abundance(self, diamond_output_file, abundance_output_file, type_abundance_output_file,
                             sample_name, ncells):
        """
        Calculate ARGs abundance for a single sample based on DIAMOND BLASTx results.

        Abundance is calculated as: (query_length / subject_length) / ncells
        i.e., normalized by estimated cell count.

        Args:
            diamond_output_file (Path): DIAMOND BLASTx alignment results (.tsv)
            abundance_output_file (Path): Output file for gene-level abundance
            type_abundance_output_file (Path): Output file for type-level abundance
            sample_name (str): Sample name
            ncells (float): Estimated cell count for normalization
        """
        self.logger.info(f"Calculating ARGs abundance for sample {sample_name}...")

        # Read DIAMOND output
        blast_df = pd.read_csv(
            str(diamond_output_file), sep="\t", header=None,
            names=["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]
        )

        # For each read, keep only the best hit (lowest e-value, then highest bitscore, longest alignment)
        blast_df = blast_df.sort_values(
            ['qseqid', 'evalue', 'bitscore', 'length'],
            ascending=[True, True, False, False]
        )
        blast_df = blast_df.drop_duplicates(subset='qseqid', keep='first')

        # Filter alignment results (based on pident, evalue, length, qcov)
        blast_df["qcov"] = blast_df["length"] / blast_df["qlen"]
        blast_df = blast_df[
            (blast_df["pident"] >= 80) &
            (blast_df["evalue"] <= 1e-7) &
            (blast_df["length"] >= 8) &
            (blast_df["qcov"] >= 0.25)
            ]

        if len(blast_df) == 0:
            self.logger.warning(f"No valid ARG alignments found for {sample_name} after filtering.")
            # Create empty output files to mark as "processed"
            pd.DataFrame().to_csv(abundance_output_file, sep="\t", index=False)
            pd.DataFrame().to_csv(type_abundance_output_file, sep="\t", index=False)
            return

        # Merge with CARD metadata to get ARG annotations
        merged_df = pd.merge(blast_df, self.card_metadata, how='left', left_on='sseqid', right_on='ARGs')

        # Calculate abundance normalized by cell count
        merged_df["abundance(copies/cell)"] = (merged_df["qlen"] / merged_df["slen"]) / ncells

        # Select relevant columns (risk_index removed)
        result = merged_df[["sseqid", "Type", "Resistance_Mechanism_Type", "abundance(copies/cell)"]]

        # Group by sseqid (ARG gene), summing abundance
        result = result.groupby('sseqid')[
            ['Type', 'Resistance_Mechanism_Type', 'abundance(copies/cell)']
        ].agg({
            'Type': 'first',
            'Resistance_Mechanism_Type': 'first',
            'abundance(copies/cell)': 'sum'
        }).reset_index()

        result["sample_name"] = sample_name

        # Aggregate by antibiotic type
        type_abundance = result.groupby(['Type'], as_index=False).agg({
            'abundance(copies/cell)': 'sum',
            'sample_name': 'first',
        })

        # Write outputs
        result.to_csv(abundance_output_file, sep="\t", index=False)
        type_abundance.to_csv(type_abundance_output_file, sep="\t", index=False)

        self.logger.info(f"ARGs abundance for sample {sample_name} saved to {abundance_output_file}.")

    def calculate_abundance(self):
        """
        Calculate ARGs abundance for all samples.
        For each sample:
          1. Get estimated cell count (ncells) from KO30 marker alignment
          2. Calculate gene-level and type-level abundance
        """
        self.logger.info("Starting ARGs abundance calculation...")

        self.sample_abundance_output_folder.mkdir(parents=True, exist_ok=True)
        self.type_abundance_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = list(self.sample_names_list)
        if not sample_names_list:
            self.logger.warning("No samples found for abundance calculation!")
            return

        self.logger.info(f"{len(sample_names_list)} samples to process for abundance calculation...")

        success_count = 0
        with tqdm(total=len(sample_names_list), desc="Calculate abundance", unit="sample") as pbar:
            for sample_name in sample_names_list:
                # Input: DIAMOND alignment results
                diamond_output_file = self.diamond_output_folder / f"{sample_name}_diamond_results.tsv"
                if not diamond_output_file.exists() or diamond_output_file.stat().st_size == 0:
                    self.logger.warning(f"DIAMOND output not found or empty for {sample_name}. "
                                        f"Skipping abundance calculation.")
                    pbar.update(1)
                    pbar.set_postfix(success=f"{success_count}/{pbar.n}")
                    continue

                # Output files
                abundance_output_file = self.sample_abundance_output_folder / f"{sample_name}_abundance.tsv"
                type_abundance_output_file = self.type_abundance_output_folder / f"{sample_name}_type_abundance.tsv"

                # Skip if both outputs already exist and are not empty
                if (abundance_output_file.exists() and abundance_output_file.stat().st_size > 0 and
                        type_abundance_output_file.exists() and type_abundance_output_file.stat().st_size > 0):
                    self.logger.info(f"Abundance outputs already exist for {sample_name}. Skipping.")
                    success_count += 1
                    pbar.update(1)
                    pbar.set_postfix(success=f"{success_count}/{pbar.n}")
                    continue

                # Get normalization factor (estimated cell count)
                ncells = self.get_ncells(sample_name)
                if ncells is None or ncells == 0:
                    self.logger.error(f"Failed to get ncells for {sample_name}. Skipping abundance calculation.")
                    pbar.update(1)
                    pbar.set_postfix(success=f"{success_count}/{pbar.n}")
                    continue

                try:
                    self._calculate_abundance(
                        diamond_output_file, abundance_output_file,
                        type_abundance_output_file, sample_name, ncells
                    )

                    # Verify outputs were generated
                    if (abundance_output_file.exists() and abundance_output_file.stat().st_size > 0 and
                            type_abundance_output_file.exists() and type_abundance_output_file.stat().st_size > 0):
                        success_count += 1
                    else:
                        self.logger.error(f"Abundance calculation failed to produce outputs for {sample_name}")

                except Exception as e:
                    self.logger.error(f"Error calculating abundance for {sample_name}: {str(e)}")

                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"ARGs abundance calculation completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def merge_results(self):
        """
        Merge all sample-level abundance results into a single summary file.
        """
        self.merge_results_output_folder.mkdir(parents=True, exist_ok=True)

        merged_result = []
        for sample_name in self.sample_names_list:
            abundance_output_file = self.sample_abundance_output_folder / f"{sample_name}_abundance.tsv"
            if not abundance_output_file.exists():
                self.logger.info(f"abundance_output_file not found for sample {sample_name}. Skipping merge.")
                continue

            # Read sample results
            sample_result = pd.read_table(abundance_output_file)

            # Calculate total abundance (copies/cell)
            total_abundance = sample_result['abundance(copies/cell)'].sum()

            # Store result (risk_abundance removed)
            merged_result.append({
                'sample_name': sample_name,
                'total_abundance(copies/cell)': total_abundance,
            })

        # Convert to DataFrame
        merged_df = pd.DataFrame(merged_result)

        # Save as TSV file
        output_tsv = self.merge_results_output_folder / "merged_abundance_results.tsv"
        merged_df.to_csv(output_tsv, sep="\t", index=False)

        self.logger.info(f"Merged results saved to {output_tsv}")

    def run(self):
        """Run the entire ARGs abundance profiling pipeline"""
        import time
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("Starting the entire ARGs abundance profiling pipeline...")
        self.logger.info("=" * 60)

        # Step 1: Estimate cell counts (ncells) via KO30 marker gene alignment
        self.logger.info("Step 1/4: Estimating cell counts (KO30 BLASTx)...")
        self.count_cells()

        # Step 2: Align clean reads to ARG-like ORFs via DIAMOND BLASTx
        self.logger.info("Step 2/4: Running DIAMOND BLASTx against ARG-like ORFs...")
        self.run_diamond()

        # Step 3: Calculate ARGs abundance for each sample
        self.logger.info("Step 3/4: Calculating ARGs abundance...")
        self.calculate_abundance()

        # Step 4: Merge all sample results into a single summary
        self.logger.info("Step 4/4: Merging results...")
        self.merge_results()

        # Calculate total runtime
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)

        self.logger.info("=" * 60)
        self.logger.info(f"ARGs abundance profiling pipeline completed successfully!")
        self.logger.info(f"Total runtime: {int(hours)}h {int(minutes)}m {int(seconds)}s")
        self.logger.info("=" * 60)

if __name__ == '__main__':

    # ARGsAbundanceCalculator
    args_abundance_profile = ARGsAbundanceProfile(input_folder="")

    args_abundance_profile.run()




