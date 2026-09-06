import shutil
from utils import setup_logging, run_command, get_sample_names_list,process_samples_in_parallel
from configure import *

def has_long_sequences(fasta_file, threshold=1000):
    """
    Check if a FASTA file contains any sequence longer than the given threshold.

    Args:
        fasta_file: Path to the FASTA file.
        threshold: Length threshold in bp (default: 1000).

    Returns:
        bool: True if at least one sequence exceeds the threshold, False otherwise.
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) > threshold:
            return True
    return False

def generate_contigs2bin(bin_files, output_file):
    """
    Convert bin files to a contig2bin table file required by DAS Tool.

    Args:
        bin_files (list): List of bin file paths.
        output_file (Path): Output file path.
    """
    with open(output_file, "w") as out_f:
        for bin_file in bin_files:
            bin_name = bin_file.stem  # Extract bin file name (without extension) as bin-ID
            for record in SeqIO.parse(bin_file, "fasta"):
                out_f.write(f"{record.id}\t{bin_name}\n")


class Binning:
    def __init__(self, input_folder):
        self.input_folder = pathlib.Path(input_folder, "preprocessing")
        self.fastp_output_folder = self.input_folder / "fastp_output"
        self.megahit_output_folder = self.input_folder / "megahit_output"
        self.prodigal_output_folder = self.input_folder / "prodigal_output"
        self.acc_output_folder = self.input_folder / "acc_output"

        self.binning_output_folder = pathlib.Path(input_folder, "binning")
        self.binning_output_folder.mkdir(parents=True, exist_ok=True)

        self.get_depth_output_folder = self.binning_output_folder / "get_depth_out"
        self.get_depth_index_output_folder = self.binning_output_folder /"get_depth_out" / "build_index"
        self.get_depth_bam_output_folder = self.binning_output_folder / "get_depth_out" / "sorted_bam"


        self.metabat_output_folder = self.binning_output_folder / "metabat_out"
        self.maxbin_output_folder = self.binning_output_folder / "maxbin_out"
        self.concoct_output_folder = self.binning_output_folder / "concoct_out"

        self.dastool_output_folder = self.binning_output_folder / "dastool_out"
        self.checkm_output_folder = self.binning_output_folder /  "checkm_out"


        # 配置 logging
        setup_logging(self.binning_output_folder, "binning.log")
        self.logger = logging.getLogger("Binning")
        # 获取样本名称和文件扩展名，并存储在实例变量中
        self.sample_names_list = get_sample_names_list(self.fastp_output_folder, "megahit")

    def filter_needless_sample(self):
        """
        Filter out samples that do not have contigs longer than 1000 bp in the ACC output.
        These samples are skipped in the binning process.
        """
        samples_to_remove = []
        for sample in self.sample_names_list:
            acc_file = self.acc_output_folder / f"{sample}_acc_contigs.fa"
            if not acc_file.exists():
                self.logger.warning(f"ACC contigs file not found for {sample}: {acc_file}. Skipping.")
                samples_to_remove.append(sample)
            elif not has_long_sequences(str(acc_file)):
                self.logger.warning(f"{sample} did not have contigs exceeding 1000bp. Skipping binning process.")
                samples_to_remove.append(sample)

        for sample in samples_to_remove:
            self.sample_names_list.remove(sample)

        self.logger.info(f"Filtered {len(samples_to_remove)} samples. "
                         f"{len(self.sample_names_list)} samples remaining for binning.")

    def generate_bam_and_depth(self):
        """
        Generate BAM files and coverage depth files for each sample.
        """
        self.get_depth_output_folder.mkdir(parents=True, exist_ok=True)
        self.get_depth_index_output_folder.mkdir(parents=True, exist_ok=True)
        self.get_depth_bam_output_folder.mkdir(parents=True, exist_ok=True)
        self.logger.info("Starting BAM and depth files generation...")

        def process_sample(sample_name):
            clean_fwd = self.fastp_output_folder / f"clean_{sample_name}_1.fastq.gz"
            clean_rev = self.fastp_output_folder / f"clean_{sample_name}_2.fastq.gz"
            if not clean_fwd.exists() or not clean_rev.exists():
                self.logger.warning(f"clean_reads not found for {sample_name}: "
                                    f"{clean_fwd.name} or {clean_rev.name}. Skipping.")
                return False

            contigs_file = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"
            if not contigs_file.exists():
                self.logger.warning(f"contigs file not found for sample {sample_name}: {contigs_file}")
                return False

            # --- Step 1: Build Bowtie2 index ---
            index_file = self.get_depth_index_output_folder / f"{sample_name}_index.1.bt2"
            index_prefix = str(index_file).split(".1.bt2")[0]
            if index_file.exists():
                self.logger.info(f"Bowtie2 index for {sample_name} already exists. Skipping.")
            else:
                bowtie2_build_command = [
                    "bowtie2-build",
                    "--quiet",
                    "--threads", str(THREAD_CONFIG["bowtie2"]),
                    str(contigs_file),
                    index_prefix,
                ]
                run_command(bowtie2_build_command, logger=self.logger,
                            log_message=f"Building Bowtie2 index for {sample_name}: {' '.join(bowtie2_build_command)}")

            # --- Step 2: Align reads and generate sorted BAM ---
            sorted_bam_file = self.get_depth_bam_output_folder / f"{sample_name}.sorted.bam"
            if sorted_bam_file.exists() and sorted_bam_file.stat().st_size > 0:
                self.logger.info(f"Bowtie2 sorted bam file for {sample_name} already exists. Skipping.")
            else:
                # Build pipeline command as a single shell string
                get_sorted_bam_command = (
                    f"bowtie2 --quiet -x {index_prefix} "
                    f"-1 {clean_fwd} -2 {clean_rev} "
                    f"--very-sensitive-local -p {THREAD_CONFIG['bowtie2']} | "
                    f"samtools view -bS -@ {THREAD_CONFIG['samtools']} | "
                    f"samtools sort -o {sorted_bam_file} -@ {THREAD_CONFIG['samtools']}"
                )
                run_command(get_sorted_bam_command, logger=self.logger, shell=True,
                            log_message=f"Generating sorted BAM for {sample_name}: {get_sorted_bam_command}")

            # --- Step 3: Index BAM file ---
            bam_index_file = sorted_bam_file.with_suffix(".bam.bai")
            if bam_index_file.exists():
                self.logger.info(f"samtools index for {sample_name} already exists. Skipping.")
            else:
                samtools_index_command = [
                    "samtools", "index",
                    "-@", str(THREAD_CONFIG["samtools"]),
                    str(sorted_bam_file),
                ]
                run_command(samtools_index_command, logger=self.logger,
                            log_message=f"samtools indexing {sample_name}")

            # --- Step 4: Summarize depth ---
            depth_file = self.get_depth_output_folder / f"{sample_name}_depth.txt"
            if depth_file.exists() and depth_file.stat().st_size > 0:
                self.logger.info(f"depth file for {sample_name} already exists. Skipping.")
                return True
            else:
                summarize_depth_command = [
                    "jgi_summarize_bam_contig_depths",
                    "--outputDepth", str(depth_file),
                    str(sorted_bam_file),
                ]
                run_command(summarize_depth_command, logger=self.logger,
                            log_message=f"Generating depth file for {sample_name}: {' '.join(summarize_depth_command)}")

            return True

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(self.sample_names_list), desc="BAM & depth", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                    self.logger, process_sample, self.sample_names_list,
                    "get_depth", THREAD_CONFIG["get_depth_pool"]
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"BAM and depth files generation completed. "
                         f"Successfully processed {success_count}/{len(self.sample_names_list)} samples.")

    def run_metabat(self):
        """
        Run MetaBAT2 binning for all samples in parallel.
        """
        self.metabat_output_folder.mkdir(parents=True, exist_ok=True)
        self.logger.info("Starting MetaBAT2 binning...")

        def process_sample(sample_name):
            try:
                # Check input files
                contigs_file = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"
                if not contigs_file.exists():
                    self.logger.warning(f"Contigs file not found for {sample_name}: {contigs_file}")
                    return False

                depth_file = self.get_depth_output_folder / f"{sample_name}_depth.txt"
                if not depth_file.exists():
                    self.logger.warning(f"Depth file not found for {sample_name}: {depth_file}")
                    return False

                # Check output directory
                sample_output_dir = self.metabat_output_folder / f"{sample_name}_metabat"
                sample_output_dir.mkdir(exist_ok=True)

                # Skip if already processed
                metabat_output = sample_output_dir / "bin.1.fa"
                if metabat_output.exists():
                    self.logger.info(f"MetaBAT2 output for {sample_name} already exists. Skipping.")
                    return True

                # Build command
                metabat_command = [
                    "metabat2",
                    "-i", str(contigs_file),
                    "-a", str(depth_file),
                    "-o", str(sample_output_dir / "bin"),
                    "-m", "1500",
                    "-t", str(THREAD_CONFIG["metabat"])
                ]

                run_command(
                    metabat_command,
                    logger=self.logger,
                    log_message=f"Running MetaBAT2 for {sample_name}: {' '.join(metabat_command)}"
                )
                return True

            except Exception as e:
                self.logger.error(f"Error processing {sample_name} with MetaBAT2: {str(e)}")
                return False

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(self.sample_names_list), desc="MetaBAT2", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                    self.logger, process_sample, self.sample_names_list,
                    "MetaBAT2", THREAD_CONFIG["metabat_pool"]
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"MetaBAT2 binning completed. "
                         f"Successfully processed {success_count}/{len(self.sample_names_list)} samples.")

    def run_maxbin(self):
        """
        Run MaxBin2 binning for all samples in parallel.
        """
        self.maxbin_output_folder.mkdir(parents=True, exist_ok=True)
        self.logger.info("Starting MaxBin2 binning...")

        def process_sample(sample_name):
            contigs_file = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"
            if not contigs_file.exists():
                self.logger.warning(f"Contigs file not found for {sample_name}: {contigs_file}")
                return False

            depth_file = self.get_depth_output_folder / f"{sample_name}_depth.txt"
            if not depth_file.exists():
                self.logger.warning(f"Depth file not found for {sample_name}: {depth_file}")
                return False

            maxbin_output = self.maxbin_output_folder / f"{sample_name}_maxbin" / "bin.001.fasta"
            if maxbin_output.exists():
                self.logger.info(f"MaxBin output for {sample_name} already exists. Skipping.")
                return True

            # Only remove the directory if it exists but the expected output is missing
            # (i.e., a previous run was interrupted)
            if maxbin_output.parent.exists():
                existing_bins = list(maxbin_output.parent.glob("*.fasta"))
                if existing_bins:
                    self.logger.info(f"MaxBin output for {sample_name} has {len(existing_bins)} bins but "
                                     f"bin.001.fasta is missing. Cleaning up and re-running.")
                    shutil.rmtree(maxbin_output.parent)
                maxbin_output.parent.mkdir(parents=True, exist_ok=True)
            else:
                maxbin_output.parent.mkdir(parents=True, exist_ok=True)

            maxbin_output_prefix = str(maxbin_output).split(".001.fasta")[0]
            maxbin_command = [
                "run_MaxBin.pl",
                "-contig", str(contigs_file),
                "-out", maxbin_output_prefix,
                "-abund", str(depth_file),
            ]
            run_command(maxbin_command, logger=self.logger,
                        log_message=f"Running MaxBin for {sample_name}: {' '.join(maxbin_command)}")
            return True

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(self.sample_names_list), desc="MaxBin2", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                    self.logger, process_sample, self.sample_names_list,
                    "MaxBin", THREAD_CONFIG["maxbin_pool"]
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"MaxBin2 binning completed. "
                         f"Successfully processed {success_count}/{len(self.sample_names_list)} samples.")

    def run_concoct(self):
        """
        Run CONCOCT binning for all samples in parallel.
        Pipeline: cut contigs → coverage table → CONCOCT → merge clustering → extract bins
        """
        self.concoct_output_folder.mkdir(parents=True, exist_ok=True)
        self.logger.info("Starting CONCOCT binning...")

        def process_sample(sample_name):
            # Check input files
            contigs_file = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"
            if not contigs_file.exists():
                self.logger.warning(f"Contigs file not found for {sample_name}: {contigs_file}")
                return False

            bam_file = self.get_depth_bam_output_folder / f"{sample_name}.sorted.bam"
            if not bam_file.exists():
                self.logger.warning(f"BAM file not found for {sample_name}: {bam_file}")
                return False

            # Define output directory
            sample_output_dir = self.concoct_output_folder / f"{sample_name}_concoct"
            sample_output_dir.mkdir(exist_ok=True)

            # Skip if final merged clustering already exists
            merged_clustering_file = sample_output_dir / "clustering_merged.csv"
            if merged_clustering_file.exists():
                self.logger.info(f"CONCOCT output for {sample_name} already exists. Skipping.")
                return True

            # --- Step 1: Cut contigs into 10kbp fragments ---
            cutup_fasta = sample_output_dir / "contigs_10K.fa"
            cutup_bed = sample_output_dir / "contigs_10K.bed"
            if not cutup_fasta.exists() or cutup_fasta.stat().st_size == 0:
                cutup_command = (
                    f"cut_up_fasta.py {contigs_file} "
                    f"-c 10000 -o 0 --merge_last "
                    f"-b {cutup_bed} > {cutup_fasta}"
                )
                run_command(cutup_command, logger=self.logger, shell=True,
                            log_message=f"Cutting contigs for {sample_name}: {cutup_command}")

            # --- Step 2: Generate coverage table ---
            coverage_table = sample_output_dir / "coverage_table.tsv"
            if not coverage_table.exists() or coverage_table.stat().st_size == 0:
                coverage_command = (
                    f"concoct_coverage_table.py {cutup_bed} {bam_file} > {coverage_table}"
                )
                run_command(coverage_command, logger=self.logger, shell=True,
                            log_message=f"Generating coverage table for {sample_name}: {coverage_command}")

            # --- Step 3: Run CONCOCT ---
            clustering_file = sample_output_dir / "clustering_gt1000.csv"
            if not clustering_file.exists():
                concoct_command = [
                    "concoct",
                    "--composition_file", str(cutup_fasta),
                    "--coverage_file", str(coverage_table),
                    "-b", str(sample_output_dir),
                    "-t", str(THREAD_CONFIG["concoct"])
                ]
                run_command(concoct_command, logger=self.logger,
                            log_message=f"Running CONCOCT for {sample_name}: {' '.join(concoct_command)}")

            # --- Step 4: Merge sub-contig clustering results ---
            if not merged_clustering_file.exists():
                merge_command = (
                    f"merge_cutup_clustering.py {clustering_file} > {merged_clustering_file}"
                )
                run_command(merge_command, logger=self.logger, shell=True,
                            log_message=f"Merging clustering for {sample_name}: {merge_command}")

            # --- Step 5: Extract final bins ---
            fasta_bins_dir = sample_output_dir / "fasta_bins"
            if not fasta_bins_dir.exists() or not list(fasta_bins_dir.glob("*.fa")):
                if not fasta_bins_dir.exists():
                    fasta_bins_dir.mkdir()
                extract_command = [
                    "extract_fasta_bins.py",
                    str(contigs_file),
                    str(merged_clustering_file),
                    "--output_path", str(fasta_bins_dir)
                ]
                run_command(extract_command, logger=self.logger,
                            log_message=f"Extracting bins for {sample_name}: {' '.join(extract_command)}")

            return True

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(self.sample_names_list), desc="CONCOCT", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                    self.logger, process_sample, self.sample_names_list,
                    "CONCOCT", THREAD_CONFIG["concoct_pool"]
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"CONCOCT binning completed. "
                         f"Successfully processed {success_count}/{len(self.sample_names_list)} samples.")

    def run_dastool(self):
        """
        Integrate MetaBAT2, MaxBin2, and CONCOCT binning results using DAS Tool.
        """
        self.dastool_output_folder.mkdir(parents=True, exist_ok=True)
        self.logger.info("Starting DASTool binning integration...")

        def process_sample(sample_name):
            # Check required input files
            contigs_file = self.megahit_output_folder / f"{sample_name}_final.contigs.fa"
            if not contigs_file.exists():
                self.logger.warning(f"Contigs file not found for {sample_name}: {contigs_file}")
                return False

            proteins_file = self.prodigal_output_folder / f"{sample_name}_prodigal.faa"
            if not proteins_file.exists():
                self.logger.warning(f"Proteins file not found for {sample_name}: {proteins_file}")
                return False

            # Skip if already processed
            dastool_summary = self.dastool_output_folder / f"{sample_name}_dastool" / "DASTool_summary.tsv"
            if dastool_summary.exists():
                self.logger.info(f"DASTool output for {sample_name} already exists. Skipping.")
                return True

            # Create output directory
            sample_output_dir = self.dastool_output_folder / f"{sample_name}_dastool"
            sample_output_dir.mkdir(parents=True, exist_ok=True)

            # --- Step 1: Collect binning results from each tool ---
            tool_dirs = {
                "MetaBAT2": self.metabat_output_folder / f"{sample_name}_metabat",
                "MaxBin2": self.maxbin_output_folder / f"{sample_name}_maxbin",
                "CONCOCT": self.concoct_output_folder / f"{sample_name}_concoct" / "fasta_bins"
            }

            active_tools = {}
            for tool_name, tool_dir in tool_dirs.items():
                if tool_name == "CONCOCT":
                    bin_files = list(tool_dir.glob("*.fa"))
                else:
                    bin_files = list(tool_dir.glob("*.fa")) if tool_name == "MetaBAT2" else list(
                        tool_dir.glob("*.fasta"))

                if not bin_files:
                    self.logger.warning(f"No bins found for {tool_name} in {sample_name}. Skipping this tool.")
                else:
                    active_tools[tool_name] = bin_files

            if len(active_tools) < 2:
                self.logger.warning(
                    f"At least 2 binning tools required for DASTool. "
                    f"Only {len(active_tools)} found for {sample_name}.")
                return False

            # --- Step 2: Generate contig2bin files for each tool ---
            contig2bin_files = []
            tool_labels = []
            for tool_name, bin_files in active_tools.items():
                output_file = sample_output_dir / f"{tool_name.lower()}_contigs2bins.tsv"
                if not output_file.exists():
                    generate_contigs2bin(bin_files, output_file)
                contig2bin_files.append(str(output_file))
                tool_labels.append(tool_name)

            # --- Step 3: Run DAS Tool ---
            dastool_output_prefix = sample_output_dir / "output"
            dastool_command = [
                "DAS_Tool",
                "-i", ",".join(contig2bin_files),
                "-l", ",".join(tool_labels),
                "-c", str(contigs_file),
                "-o", str(dastool_output_prefix),
                "--proteins", str(proteins_file),
                "--write_bins",
                "--search_engine", "diamond",
                "--threads", str(THREAD_CONFIG["dastool"])
            ]

            run_command(" ".join(dastool_command), logger=self.logger, shell=True,
                        log_message=f"Running DASTool for {sample_name}: {' '.join(dastool_command)}")

            return True

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(self.sample_names_list), desc="DASTool", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                    self.logger, process_sample, self.sample_names_list,
                    "DASTool", THREAD_CONFIG["dastool_pool"]
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"DASTool binning integration completed. "
                         f"Successfully processed {success_count}/{len(self.sample_names_list)} samples.")

    def run_checkm(self):
        """
        Run CheckM quality assessment on DAS Tool output bins.
        """
        self.checkm_output_folder.mkdir(parents=True, exist_ok=True)
        self.logger.info("Starting CheckM quality assessment...")

        def process_sample(sample_name):
            # Define file paths
            dastool_bins_file = self.dastool_output_folder / f"{sample_name}_dastool" / "output_DASTool_bins"
            checkm_output_folder = self.checkm_output_folder / f"{sample_name}_checkm"

            # Create output folder
            checkm_output_folder.mkdir(parents=True, exist_ok=True)

            # Check if DAS Tool bins exist
            if not dastool_bins_file.exists():
                self.logger.warning(f"DAS Tool bins for sample {sample_name} not found. Skipping.")
                return False

            # Check if CheckM output already exists
            checkm_stats_file = checkm_output_folder / "storage" / "bin_stats_ext.tsv"
            if checkm_stats_file.exists():
                self.logger.info(f"CheckM output for {sample_name} already exists. Skipping.")
                return True  # Fix: return True (already processed successfully)

            # Build CheckM command
            checkm_command = [
                "checkm", "lineage_wf",
                "-t", str(THREAD_CONFIG["checkm"]),
                "-x", "fa",
                str(dastool_bins_file),
                str(checkm_output_folder)
            ]

            run_command(" ".join(checkm_command), logger=self.logger, shell=True,
                        log_message=f"Running CheckM for {sample_name}: {' '.join(checkm_command)}")
            return True

        # Process all samples in parallel with progress bar
        success_count = 0
        with tqdm(total=len(self.sample_names_list), desc="CheckM", unit="sample") as pbar:
            for result in process_samples_in_parallel(
                    self.logger, process_sample, self.sample_names_list,
                    "CheckM", THREAD_CONFIG["checkm_pool"]
            ):
                if result:
                    success_count += 1
                pbar.update(1)
                pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"CheckM quality assessment completed. "
                         f"Successfully processed {success_count}/{len(self.sample_names_list)} samples.")

    def run(self):
        """Run the entire binning pipeline."""
        import time
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("Starting binning pipeline...")
        self.logger.info("=" * 60)

        # Step 1: Filter samples without long contigs
        self.logger.info("Step 1/6: Filtering samples without long contigs...")
        self.filter_needless_sample()

        # Step 2: Generate BAM and depth files
        self.logger.info("Step 2/6: Generating BAM and depth files...")
        self.generate_bam_and_depth()

        # Step 3: Run MetaBAT2
        self.logger.info("Step 3/6: Running MetaBAT2 binning...")
        self.run_metabat()

        # Step 4: Run MaxBin2
        self.logger.info("Step 4/6: Running MaxBin2 binning...")
        self.run_maxbin()

        # Step 5: Run CONCOCT
        self.logger.info("Step 5/6: Running CONCOCT binning...")
        self.run_concoct()

        # Step 6: Integrate with DASTool and assess with CheckM
        self.logger.info("Step 6/6: Running DASTool integration and CheckM quality assessment...")
        self.run_dastool()
        self.run_checkm()

        # Calculate total runtime
        elapsed_time = time.time() - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)

        self.logger.info("=" * 60)
        self.logger.info(f"Binning pipeline completed successfully!")
        self.logger.info(f"Total runtime: {int(hours)}h {int(minutes)}m {int(seconds)}s")
        self.logger.info("=" * 60)


if __name__ == '__main__':
    binning = Binning("")
    binning.run()