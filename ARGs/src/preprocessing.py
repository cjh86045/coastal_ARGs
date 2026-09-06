import shutil
from utils import setup_logging, run_command, get_sample_names_list, process_samples_in_parallel
from configure import *
from Bio import SeqIO



class Preprocessing:
    def __init__(self, input_folder):
        # Input data root directory path: input_folder/data/
        self.data_folder = pathlib.Path(input_folder, "data")

        # Preprocessing pipeline: generate inputs for downstream analysis

        # Create preprocessing output folder
        self.preprocessing_folder = pathlib.Path(input_folder, 'preprocessing')
        self.preprocessing_folder.mkdir(parents=True, exist_ok=True)

        # fastp output directory
        self.fastp_output_folder = self.preprocessing_folder / "fastp_output"

        # MEGAHIT assembly output directory
        self.megahit_output_folder = self.preprocessing_folder / "megahit_output"

        # Prodigal ORF prediction output directory
        self.prodigal_output_folder = self.preprocessing_folder / "prodigal_output"

        # DIAMOND alignment results against ARGs database output directory
        self.ARGs_blast_result_output_folder = self.preprocessing_folder / "ARGs_blast_result_output"

        # ARG-like ORFs output directory
        self.arg_like_orfs_folder = self.preprocessing_folder / "arg_like_orfs_output"

        # ARG-carrying contigs (ACCs) output directory
        self.acc_output_folder = self.preprocessing_folder / "acc_output"

        # ORFs on ARG-carrying contigs output directory
        self.acc_orfs_output_folder = self.preprocessing_folder / "acc_orfs"

        # Configure logging
        setup_logging(self.preprocessing_folder, "preprocessing.log")
        self.logger = logging.getLogger("Preprocessing")

    def run_fastp(self):
        """Run fastp for raw data quality control"""
        self.fastp_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = get_sample_names_list(self.data_folder, flag="fastp")
        if not sample_names_list:
            self.logger.warning("No samples found for fastp processing!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting fastp processing...")

        def process_sample(sample_name):
            """Process a single sample with fastp"""
            # Input files: raw reads
            raw_fwd = pathlib.Path(self.data_folder) / f"{sample_name}_1.fastq.gz"
            raw_rev = pathlib.Path(self.data_folder) / f"{sample_name}_2.fastq.gz"

            # Check input files exist
            if not raw_fwd.exists() or not raw_rev.exists():
                self.logger.warning(f"Input files not found for {sample_name}: "
                                    f"{raw_fwd.name} or {raw_rev.name}. Skipping.")
                return False

            # Check input files are not empty
            if raw_fwd.stat().st_size == 0 or raw_rev.stat().st_size == 0:
                self.logger.warning(f"Input files are empty for {sample_name}. Skipping.")
                return False

            # Output files: clean reads
            clean_fwd = self.fastp_output_folder / f"clean_{sample_name}_1.fastq.gz"
            clean_rev = self.fastp_output_folder / f"clean_{sample_name}_2.fastq.gz"

            # Check if output files already exist and are not empty
            if clean_fwd.exists() and clean_rev.exists() and \
                    clean_fwd.stat().st_size > 0 and clean_rev.stat().st_size > 0:
                self.logger.info(f"Output files already exist for {sample_name}. Skipping.")
                return True

            # Get thread count from config
            threads = THREAD_CONFIG.get("fastp", 4)  # Default to 4 threads
            if threads <= 0:
                self.logger.warning(f"Invalid thread count {threads} for fastp. Using default 4.")
                threads = 4

            # Build fastp command
            fastp_command = [
                "fastp",
                "-i", str(raw_fwd),  # Input file 1 (forward reads)
                "-I", str(raw_rev),  # Input file 2 (reverse reads)
                "-o", str(clean_fwd),  # Output file 1 (clean forward reads)
                "-O", str(clean_rev),  # Output file 2 (clean reverse reads)
                "--json", str(self.fastp_output_folder / f"{sample_name}_fastp.json"),  # JSON report
                "--html", str(self.fastp_output_folder / f"{sample_name}_fastp.html"),  # HTML report
                "--thread", str(threads),  # Number of threads
            ]

            # Run fastp command
            try:
                self.logger.info(f"Running fastp for sample {sample_name}...")
                run_command(fastp_command, logger=self.logger,
                            log_message=f"fastp command: {' '.join(fastp_command)}")

                # Verify output files were generated successfully
                if clean_fwd.exists() and clean_rev.exists() and \
                        clean_fwd.stat().st_size > 0 and clean_rev.stat().st_size > 0:
                    self.logger.info(f"fastp completed successfully for {sample_name}")
                    return True
                else:
                    self.logger.error(f"fastp failed to produce valid output for {sample_name}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running fastp for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names_list), desc="fastp QC", unit="sample") as pbar:
            # 先获取所有结果（返回 (sample_name, success) 元组列表）
            results = process_samples_in_parallel(
                self.logger, process_sample, sample_names_list,
                "fastp", THREAD_CONFIG.get("fastp_pool", 2)
            )

            # 遍历结果
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"fastp processing completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def run_megahit(self):
        """Run MEGAHIT for metagenome assembly"""
        self.megahit_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = get_sample_names_list(self.fastp_output_folder, flag="megahit")
        if not sample_names_list:
            self.logger.warning("No samples found for megahit assembly!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting megahit assembly processing...")

        def process_sample(sample_name):
            """Process a single sample with MEGAHIT assembly"""
            # Input files: clean reads
            clean_fwd = self.fastp_output_folder / f"clean_{sample_name}_1.fastq.gz"
            clean_rev = self.fastp_output_folder / f"clean_{sample_name}_2.fastq.gz"

            # Check input files exist and are not empty
            if not clean_fwd.exists() or not clean_rev.exists():
                self.logger.warning(f"Input files not found for {sample_name}: "
                                    f"{clean_fwd.name} or {clean_rev.name}. Skipping.")
                return False

            if clean_fwd.stat().st_size == 0 or clean_rev.stat().st_size == 0:
                self.logger.warning(f"Input files are empty for {sample_name}. Skipping.")
                return False

            # Output file: final contigs
            final_contigs = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"

            # Case 1: Final output already exists and is not empty -> skip
            if final_contigs.exists() and final_contigs.stat().st_size > 0:
                self.logger.info(f"Final contigs already exist for {sample_name}. Skipping.")
                return True

            # Define temporary output folder
            sample_output_folder = self.megahit_output_folder / f"megahit_{sample_name}"
            tmp_final_contig = sample_output_folder / "final.contigs.fa"

            # Get thread count from config
            threads = THREAD_CONFIG.get("megahit", 8)
            if threads <= 0:
                self.logger.warning(f"Invalid thread count {threads} for megahit. Using default 8.")
                threads = 8

            # Build megahit command
            megahit_command = [
                "megahit",
                "-1", str(clean_fwd),
                "-2", str(clean_rev),
                "-o", str(sample_output_folder),
                "--presets", "meta-large",
                "--num-cpu-threads", str(threads),
            ]

            # Case 2: Resume from interrupted assembly
            if sample_output_folder.exists():
                self.logger.info(f"Found existing MEGAHIT output dir for {sample_name}. "
                                 f"Attempting to resume assembly with --continue...")
                megahit_command.append("--continue")
            else:
                self.logger.info(f"Starting fresh MEGAHIT assembly for {sample_name}...")

            # Run megahit command
            try:
                run_command(megahit_command, logger=self.logger,
                            log_message=f"megahit command: {' '.join(megahit_command)}")

                # Verify final contigs were generated
                if tmp_final_contig.exists() and tmp_final_contig.stat().st_size > 0:
                    # Step 1: Move the result file to the final destination
                    shutil.move(str(tmp_final_contig), str(final_contigs))
                    self.logger.info(f"Moved final.contigs.fa to {final_contigs}")

                    # Step 2: Only after confirming the move succeeded, clean up temp folder
                    if final_contigs.exists() and final_contigs.stat().st_size > 0:
                        shutil.rmtree(str(sample_output_folder))
                        self.logger.info(f"Deleted megahit temporary folder: {sample_output_folder}")
                        return True
                    else:
                        self.logger.error(f"Failed to verify final contigs at {final_contigs}")
                        return False
                else:
                    self.logger.error(f"final.contigs.fa not found or empty in {sample_output_folder}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running megahit for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names_list), desc="MEGAHIT assembly", unit="sample") as pbar:
            # 先获取所有结果（返回 (sample_name, success) 元组列表）
            results = process_samples_in_parallel(
                self.logger, process_sample, sample_names_list,
                "megahit", THREAD_CONFIG.get("megahit_pool", 2)
            )

            # 遍历结果
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"Megahit assembly completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def run_prodigal(self):
        """Run Prodigal for ORF prediction on assembled contigs"""
        self.prodigal_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names = get_sample_names_list(self.megahit_output_folder, flag="prodigal")
        if not sample_names:
            self.logger.warning("No samples found for Prodigal ORF prediction!")
            return

        self.logger.info(f"{len(sample_names)} samples Starting Prodigal ORF prediction...")

        def process_sample(sample_name):
            """Process a single sample with Prodigal"""
            # Input file: assembled contigs
            final_contigs = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"

            # Check input file exists and is not empty
            if not final_contigs.exists():
                self.logger.warning(f"Input contigs not found for {sample_name}: "
                                    f"{final_contigs.name}. Skipping.")
                return False

            if final_contigs.stat().st_size == 0:
                self.logger.warning(f"Input contigs are empty for {sample_name}. Skipping.")
                return False

            # Output files (only gkb and faa are required)
            output_gkb = self.prodigal_output_folder / f"{sample_name}_prodigal.gkb"
            output_faa = self.prodigal_output_folder / f"{sample_name}_prodigal.faa"

            # Check if all required output files already exist and are not empty
            required_outputs = [output_gkb, output_faa]
            if all(f.exists() and f.stat().st_size > 0 for f in required_outputs):
                self.logger.info(f"All Prodigal outputs already exist for {sample_name}. Skipping.")
                return True

            # If some outputs exist but not all, remove incomplete outputs and re-run
            for f in required_outputs:
                if f.exists():
                    self.logger.warning(f"Removing incomplete output file: {f}")
                    f.unlink()

            # Build Prodigal command
            prodigal_command = [
                "prodigal",
                "-i", str(final_contigs),  # Input contigs file
                "-o", str(output_gkb),  # Output gene coordinates (GBK format)
                "-a", str(output_faa),  # Output protein sequences
                "-p", "meta",  # Metagenomic mode
                "-q"  # Quiet mode (suppress stderr output)
            ]

            # Run Prodigal command
            try:
                self.logger.info(f"Running Prodigal for sample {sample_name}...")
                run_command(prodigal_command, logger=self.logger,
                            log_message=f"Prodigal command: {' '.join(prodigal_command)}")

                # Verify all output files were generated successfully
                if all(f.exists() and f.stat().st_size > 0 for f in required_outputs):
                    self.logger.info(f"Prodigal completed successfully for {sample_name}")
                    return True
                else:
                    self.logger.error(f"Prodigal failed to produce all outputs for {sample_name}")
                    for f in required_outputs:
                        if not f.exists():
                            self.logger.error(f"  Missing: {f}")
                        elif f.stat().st_size == 0:
                            self.logger.error(f"  Empty: {f}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running Prodigal for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names), desc="Prodigal ORF prediction", unit="sample") as pbar:
            # 先获取所有结果（返回 (sample_name, success) 元组列表）
            results = process_samples_in_parallel(
                self.logger, process_sample, sample_names,
                "prodigal", THREAD_CONFIG.get("prodigal_pool", 2)
            )

            # 遍历结果
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"Prodigal ORF prediction completed. "
                         f"Successfully processed {success_count}/{len(sample_names)} samples.")

    def run_diamond(self):
        """Run DIAMOND BLASTP to align predicted ORFs against CARD database"""
        self.ARGs_blast_result_output_folder.mkdir(parents=True, exist_ok=True)

        # Check database exists once before processing all samples
        if not CARD_DB.exists():
            self.logger.error(f"CARD database not found at {CARD_DB}. "
                              f"Please check the database path in configure.py.")
            return

        sample_names_list = get_sample_names_list(self.prodigal_output_folder, flag="diamond")
        if not sample_names_list:
            self.logger.warning("No samples found for DIAMOND BLASTP processing!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting DIAMOND BLASTP against CARD database...")

        def process_sample(sample_name):
            """Process a single sample with DIAMOND BLASTP"""
            # Input file: Prodigal predicted ORFs (protein sequences)
            orfs_file = self.prodigal_output_folder / f"{sample_name}_prodigal.faa"

            # Check input file exists and is not empty
            if not orfs_file.exists():
                self.logger.warning(f"Input ORFs file not found for {sample_name}: "
                                    f"{orfs_file.name}. Skipping.")
                return False

            if orfs_file.stat().st_size == 0:
                self.logger.warning(f"Input ORFs file is empty for {sample_name}. Skipping.")
                return False

            # Output file: DIAMOND alignment results
            diamond_output_file = self.ARGs_blast_result_output_folder / f"{sample_name}_card_results.tsv"

            # Check if output already exists and is not empty
            if diamond_output_file.exists() and diamond_output_file.stat().st_size > 0:
                self.logger.info(f"DIAMOND output already exists for {sample_name}. Skipping.")
                return True

            # Get thread count from config
            threads = THREAD_CONFIG.get("diamond", 8)
            if threads <= 0:
                self.logger.warning(f"Invalid thread count {threads} for diamond. Using default 8.")
                threads = 8

            # Build DIAMOND command
            diamond_command = [
                "diamond", "blastp",
                "-d", str(CARD_DB),  # CARD database
                "-q", str(orfs_file),  # Input ORFs file (protein sequences)
                "-o", str(diamond_output_file),  # Output results file
                "--threads", str(threads),  # Number of threads
                "--evalue", "1e-7",  # E-value threshold
                "--id", "70",  # Sequence identity threshold (%)
                "--query-cover", "70",  # Query coverage threshold (%)
                "--more-sensitive",  # More sensitive alignment mode
                "--max-hsps", "1",  # Maximum HSPs per query
                "--max-target-seqs", "1",  # Maximum target sequences per query
                "--outfmt", "6",  # Output format: tabular with query and subject IDs
                "--quiet",  # Suppress non-error output
            ]

            # Run DIAMOND command
            try:
                self.logger.info(f"Running DIAMOND BLASTP for sample {sample_name}...")
                run_command(diamond_command, logger=self.logger,
                            log_message=f"DIAMOND command: {' '.join(diamond_command)}")

                # Verify output file was generated successfully
                if diamond_output_file.exists() and diamond_output_file.stat().st_size > 0:
                    self.logger.info(f"DIAMOND BLASTP completed successfully for {sample_name}")
                    return True
                else:
                    self.logger.error(f"DIAMOND BLASTP failed to produce output for {sample_name}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running DIAMOND BLASTP for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names_list), desc="DIAMOND BLASTP", unit="sample") as pbar:
            results = process_samples_in_parallel(self.logger, process_sample, sample_names_list,
                "diamond", THREAD_CONFIG.get("diamond_pool", 2))
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"DIAMOND BLASTP completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def _extract_arg_like_orfs(self, orfs_file, diamond_output_file, arg_like_orfs_file, sample_name):
        """
        Extract ORF sequences that aligned to the CARD database.

        Parameters:
        orfs_file (Path): Prodigal predicted ORFs (.faa)
        diamond_output_file (Path): DIAMOND BLASTP results (.tsv)
        arg_like_orfs_file (Path): Output ARG-like ORFs file (.faa)
        sample_name (str): Sample name
        """
        self.logger.info(f"Extracting ARG-like ORFs for sample {sample_name}...")

        # Read only the first two columns (qseqid, sseqid)
        df = pd.read_table(str(diamond_output_file), sep="\t", header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
       "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

        if df.empty:
            self.logger.warning(f"DIAMOND output is empty for {sample_name}. No ARG-like ORFs found.")
            # Create an empty output file to mark as "processed"
            Path(arg_like_orfs_file).touch()
            return

        # Build mapping: ORF_ID -> CARD_gene_ID (sseqid)
        diamond_output = dict(zip(df["qseqid"], df["sseqid"]))

        # Extract matching ORF sequences
        arg_like_orfs = []
        for orf in SeqIO.parse(str(orfs_file), "fasta"):
            if orf.id in diamond_output:
                arg_like_orfs.append(SeqIO.SeqRecord(
                    orf.seq,
                    id=diamond_output[orf.id],  # CARD gene ID as sequence ID
                    description=f"{sample_name}|{orf.id}"  # Keep original ORF ID for traceability
                ))

        # Write output
        SeqIO.write(arg_like_orfs, str(arg_like_orfs_file), "fasta")

        self.logger.info(f"Extracted {len(arg_like_orfs)} ARG-like ORFs for {sample_name} "
                         f"-> saved to {arg_like_orfs_file}")

    def extract_arg_like_orfs(self):
        """Extract ARG-like ORF sequences from DIAMOND alignment results"""
        self.arg_like_orfs_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = get_sample_names_list(self.prodigal_output_folder, flag="diamond")
        if not sample_names_list:
            self.logger.warning("No samples found for ARG-like ORF extraction!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting ARG-like ORF extraction...")

        success_count = 0
        with tqdm(total=len(sample_names_list), desc="Extract ARG-like ORFs", unit="sample") as pbar:
            for sample_name in sample_names_list:
                # Input files
                diamond_output_file = self.ARGs_blast_result_output_folder / f"{sample_name}_card_results.tsv"
                orfs_file = self.prodigal_output_folder / f"{sample_name}_prodigal.faa"

                # Check input files exist and are not empty
                if not diamond_output_file.exists() or not orfs_file.exists():
                    self.logger.warning(f"Input files not found for {sample_name}: "
                                        f"{diamond_output_file.name} or {orfs_file.name}. Skipping.")
                    pbar.update(1)
                    pbar.set_postfix(success=f"{success_count}/{pbar.n}")
                    continue

                if diamond_output_file.stat().st_size == 0 or orfs_file.stat().st_size == 0:
                    self.logger.warning(f"Input files are empty for {sample_name}. Skipping.")
                    pbar.update(1)
                    pbar.set_postfix(success=f"{success_count}/{pbar.n}")
                    continue

                # Output file (overwrite if exists)
                arg_like_orfs_file = self.arg_like_orfs_folder / f"{sample_name}_arg_like_orfs.faa"

                # Extract ARG-like ORFs (no skip check, always overwrite)
                try:
                    self.logger.info(f"Extracting ARG-like ORFs for {sample_name}...")
                    self._extract_arg_like_orfs(orfs_file, diamond_output_file, arg_like_orfs_file, sample_name)

                    # Verify output was generated
                    if arg_like_orfs_file.exists() and arg_like_orfs_file.stat().st_size > 0:
                        self.logger.info(f"ARG-like ORF extraction completed for {sample_name}")
                        success_count += 1
                    else:
                        self.logger.error(f"ARG-like ORF extraction failed to produce output for {sample_name}")

                except Exception as e:
                    self.logger.error(f"Error extracting ARG-like ORFs for {sample_name}: {str(e)}")

                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"ARG-like ORF extraction completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def extract_contigs_with_arg_like_orfs(self, contig_file, diamond_output_file, acc_output_file, sample_name):
        """
        Extract contigs containing ARG-like ORFs from the assembled contigs.

        Parameters:
        contig_file (Path): Assembled contigs file (.fa)
        diamond_output_file (Path): DIAMOND BLASTP results (.tsv)
        acc_output_file (Path): Output contigs file (.fa)
        sample_name (str): Sample name
        """
        self.logger.info(f"Extracting contigs containing ARG-like ORFs for sample {sample_name}...")

        # Read DIAMOND output
        df = pd.read_table(str(diamond_output_file), sep="\t", header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
       "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

        if df.empty:
            self.logger.warning(f"DIAMOND output is empty for {sample_name}. No ARG-carrying contigs found.")
            # Create an empty output file to mark as "processed"
            Path(acc_output_file).touch()
            return

        # Extract contig ID from ORF ID (e.g., "k141_1_1" -> "k141_1")
        # Prodigal ORF IDs follow the pattern: {contig_id}_{orf_number}
        df["acc_contig_ids"] = df["qseqid"].str.rsplit("_", n=1).str[0]

        # Create description with CARD gene ID and ORF ID
        df["acc_contig_description"] = df["sseqid"] + "(" + df["qseqid"] + ")"

        # Group by contig ID, join multiple descriptions with "@"
        result = df.groupby('acc_contig_ids')['acc_contig_description'].apply('@'.join).reset_index()
        diamond_output_merge = dict(zip(result["acc_contig_ids"], result["acc_contig_description"]))

        # Extract matching contigs
        arg_carry_contigs = []
        for contig in SeqIO.parse(str(contig_file), "fasta"):
            if contig.id in diamond_output_merge:
                arg_carry_contigs.append(SeqIO.SeqRecord(
                    contig.seq,
                    id=contig.id,
                    description=diamond_output_merge[contig.id]
                ))

        # Write output
        SeqIO.write(arg_carry_contigs, str(acc_output_file), "fasta")

        self.logger.info(f"Extracted {len(arg_carry_contigs)} ARG-carrying contigs for {sample_name} "
                         f"-> saved to {acc_output_file}")

    def extract_arg_carry_contig(self):
        """Extract contigs containing ARG-like ORFs for all samples"""
        self.acc_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = get_sample_names_list(self.megahit_output_folder, flag="prodigal")
        if not sample_names_list:
            self.logger.warning("No samples found for ARG-carrying contig extraction!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting ARG-carrying contig extraction...")

        def process_sample(sample_name):
            """Process a single sample to extract ARG-carrying contigs"""
            # Input files
            final_contigs = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"
            diamond_output_file = self.ARGs_blast_result_output_folder / f"{sample_name}_card_results.tsv"

            # Check input files exist and are not empty
            if not final_contigs.exists() or not diamond_output_file.exists():
                self.logger.warning(f"Input files not found for {sample_name}: "
                                    f"{final_contigs.name} or {diamond_output_file.name}. Skipping.")
                return False

            if final_contigs.stat().st_size == 0 or diamond_output_file.stat().st_size == 0:
                self.logger.warning(f"Input files are empty for {sample_name}. Skipping.")
                return False

            # Output file
            acc_output_file = self.acc_output_folder / f"{sample_name}_acc_contigs.fa"

            # Check if output already exists and is not empty
            if acc_output_file.exists() and acc_output_file.stat().st_size > 0:
                self.logger.info(f"ACC output already exists for {sample_name}. Skipping.")
                return True

            # Extract ARG-carrying contigs
            try:
                self.extract_contigs_with_arg_like_orfs(final_contigs, diamond_output_file, acc_output_file,
                                                        sample_name)

                # Verify output was generated
                if acc_output_file.exists() and acc_output_file.stat().st_size > 0:
                    self.logger.info(f"ARG-carrying contig extraction completed for {sample_name}")
                    return True
                else:
                    self.logger.error(f"ARG-carrying contig extraction failed to produce output for {sample_name}")
                    return False

            except Exception as e:
                self.logger.error(f"Error extracting ARG-carrying contigs for {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names_list), desc="Extract ARG-carrying contigs", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                self.logger, process_sample, sample_names_list,
                "extract_arg_carry_contig", THREAD_CONFIG.get("extract_pool", 2)
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"ARG-carrying contig extraction completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def get_acc_orfs(self):
        """Predict ORFs on ARG-carrying contigs (ACCs) using Prodigal"""
        self.acc_orfs_output_folder.mkdir(parents=True, exist_ok=True)

        sample_names_list = get_sample_names_list(self.acc_output_folder, flag="binning")
        if not sample_names_list:
            self.logger.warning("No samples found for ACC ORF prediction!")
            return

        self.logger.info(f"{len(sample_names_list)} samples Starting ACC ORF prediction...")

        def process_sample(sample_name):
            """Process a single sample to predict ORFs on ACCs"""
            # Input file: ARG-carrying contigs
            acc_file = self.acc_output_folder / f"{sample_name}_acc_contigs.fa"

            # Check input file exists and is not empty
            if not acc_file.exists():
                self.logger.warning(f"ACC file not found for {sample_name}: "
                                    f"{acc_file.name}. Skipping.")
                return False

            if acc_file.stat().st_size == 0:
                self.logger.warning(f"ACC file is empty for {sample_name}. Skipping.")
                return False

            # Output files
            acc_orfs_file = self.acc_orfs_output_folder / f"{sample_name}_acc_orfs.faa"
            acc_orfs_gbk = self.acc_orfs_output_folder / f"{sample_name}_acc_orfs.gbk"

            # Check if all required output files already exist and are not empty
            required_outputs = [acc_orfs_file, acc_orfs_gbk]
            if all(f.exists() and f.stat().st_size > 0 for f in required_outputs):
                self.logger.info(f"ACC ORF outputs already exist for {sample_name}. Skipping.")
                return True

            # If some outputs exist but not all, remove incomplete outputs and re-run
            for f in required_outputs:
                if f.exists():
                    self.logger.warning(f"Removing incomplete output file: {f}")
                    f.unlink()

            # Build Prodigal command
            prodigal_command = [
                "prodigal",
                "-i", str(acc_file),  # Input ACC contigs file
                "-o", str(acc_orfs_gbk),  # Output gene coordinates (GBK format)
                "-a", str(acc_orfs_file),  # Output protein sequences
                "-p", "meta",  # Metagenomic mode
                "-q"  # Quiet mode (suppress stderr output)
            ]

            # Run Prodigal command
            try:
                self.logger.info(f"Running Prodigal for ACC ORF prediction on {sample_name}...")
                run_command(prodigal_command, logger=self.logger,
                            log_message=f"Prodigal command: {' '.join(prodigal_command)}")

                # Verify all output files were generated successfully
                if all(f.exists() and f.stat().st_size > 0 for f in required_outputs):
                    self.logger.info(f"ACC ORF prediction completed successfully for {sample_name}")
                    return True
                else:
                    self.logger.error(f"ACC ORF prediction failed to produce all outputs for {sample_name}")
                    for f in required_outputs:
                        if not f.exists():
                            self.logger.error(f"  Missing: {f}")
                        elif f.stat().st_size == 0:
                            self.logger.error(f"  Empty: {f}")
                    return False

            except Exception as e:
                self.logger.error(f"Error running Prodigal for ACC ORF prediction on {sample_name}: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(sample_names_list), desc="ACC ORF prediction", unit="sample") as pbar:
            results = process_samples_in_parallel(self.logger, process_sample, sample_names_list,
                "get_acc_orfs", THREAD_CONFIG.get("prodigal_pool", 2))
            for sample_name, success in results:
                if success:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"ACC ORF prediction completed. "
                         f"Successfully processed {success_count}/{len(sample_names_list)} samples.")

    def run(self):
        """Run the entire preprocessing pipeline"""
        import time
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("Starting the entire preprocessing pipeline...")
        self.logger.info("=" * 60)

        # Step 1: fastp quality control
        self.logger.info("Step 1/7: Running fastp quality control...")
        self.run_fastp()

        # Step 2: MEGAHIT assembly
        self.logger.info("Step 2/7: Running MEGAHIT assembly...")
        self.run_megahit()

        # Step 3: Prodigal ORF prediction
        self.logger.info("Step 3/7: Running Prodigal ORF prediction...")
        self.run_prodigal()

        # Step 4: DIAMOND alignment against CARD database
        self.logger.info("Step 4/7: Running DIAMOND BLASTP against CARD database...")
        self.run_diamond()

        # Step 5: Extract ARG-like ORFs
        self.logger.info("Step 5/7: Extracting ARG-like ORFs...")
        self.extract_arg_like_orfs()

        # Step 6: Extract ARG-carrying contigs (ACCs)
        self.logger.info("Step 6/7: Extracting ARG-carrying contigs...")
        self.extract_arg_carry_contig()

        # Step 7: Predict ORFs on ACCs
        self.logger.info("Step 7/7: Predicting ORFs on ARG-carrying contigs...")
        self.get_acc_orfs()

        # Calculate total runtime
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)

        self.logger.info("=" * 60)
        self.logger.info(f"Entire preprocessing pipeline completed successfully!")
        self.logger.info(f"Total runtime: {int(hours)}h {int(minutes)}m {int(seconds)}s")
        self.logger.info("=" * 60)


if __name__ == '__main__':
    # Initialize Preprocessing class with input folder path
    preprocessing = Preprocessing(input_folder="")

    preprocessing.run()


