from configure import *

# ====================================================================== #
# Configuration (user should modify as needed)
# ====================================================================== #
EXCEL_PATH = ""
SRA_OUTPUT_DIR = ""
FASTQ_OUTPUT_DIR = ""
SHEET_NAME = ""

# Performance settings
DOWNLOAD_POOL_SIZE = 12   # prefetch parallel workers
DUMP_POOL_SIZE = 6       # fasterq-dump parallel workers
DUMP_THREADS = 24         # threads per fasterq-dump / pigz


class SRAProcessor:
    """Download SRA files via prefetch and convert them to paired-end FASTQ."""

    def __init__(self, excel_path=EXCEL_PATH, sra_dir=SRA_OUTPUT_DIR,
                 fastq_dir=FASTQ_OUTPUT_DIR, sheet_name=SHEET_NAME):
        """
        Initialize the SRA downloader/converter.

        Args:
            excel_path: Path to the Excel file with run IDs.
            sra_dir:    Directory to store downloaded .sra files.
            fastq_dir:  Directory to store converted FASTQ files.
            sheet_name: Sheet name in the Excel file.
        """
        # Store config as instance attributes (IMPORTANT: must be set before _get_runid_list)
        self.excel_path = excel_path
        self.sheet_name = sheet_name
        self.sra_dir = Path(sra_dir)
        self.fastq_dir = Path(fastq_dir)

        self.sra_dir.mkdir(parents=True, exist_ok=True)
        self.fastq_dir.mkdir(parents=True, exist_ok=True)

        # Configure logging
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        self.logger = logging.getLogger("SRAProcessor")

        # Now safe to call _get_runid_list (uses self.excel_path / self.sheet_name)
        self.runid_list = self._get_runid_list()

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def _get_runid_list(self):
        """Read run IDs from the configured Excel sheet."""
        df = pd.read_excel(self.excel_path, sheet_name=self.sheet_name)
        runids = df["runid"].dropna().astype(str).str.strip().unique().tolist()
        self.logger.info(f"Loaded {len(runids)} run IDs from {self.excel_path}")
        return runids

    def _run_cmd(self, cmd, step_name):
        """Run a subprocess command, logging on failure."""
        self.logger.debug(f"{step_name}: {' '.join(map(str, cmd))}")
        subprocess.run(cmd, check=True, capture_output=True)

    # ------------------------------------------------------------------ #
    # Step 1: download SRA
    # ------------------------------------------------------------------ #
    def download_sra(self):
        """Download all SRA files in parallel using prefetch."""

        def _download_single(run_id):
            sra_file = self.sra_dir / run_id / f"{run_id}.sra"

            # Skip if already downloaded and non-empty
            if sra_file.exists() and sra_file.stat().st_size > 0:
                self.logger.info(f"Already downloaded, skipping: {run_id}")
                return run_id, True

            try:
                sra_file.parent.mkdir(parents=True, exist_ok=True)
                cmd = [
                    "prefetch", run_id,
                    "--max-size", "100G",
                    "--output-directory", str(self.sra_dir),
                ]
                self._run_cmd(cmd, "prefetch")

                if sra_file.exists() and sra_file.stat().st_size > 0:
                    self.logger.info(f"Downloaded: {run_id}")
                    return run_id, True
                else:
                    self.logger.warning(f"Download finished but file is missing/empty: {run_id}")
                    return run_id, False

            except subprocess.CalledProcessError as e:
                self.logger.error(f"prefetch failed for {run_id}: {e}")
                return run_id, False
            except Exception as e:
                self.logger.error(f"Unexpected error for {run_id}: {e}")
                return run_id, False

        self.logger.info("=" * 50)
        self.logger.info(f"Downloading {len(self.runid_list)} SRA files...")
        self.logger.info("=" * 50)

        results = []
        with tqdm(total=len(self.runid_list), desc="Download SRA", unit="run") as pbar:
            with ThreadPoolExecutor(max_workers=DOWNLOAD_POOL_SIZE) as executor:
                for result in executor.map(_download_single, self.runid_list):
                    results.append(result)
                    pbar.update(1)

        success_ids = [r[0] for r in results if r[1]]
        failed_ids = [r[0] for r in results if not r[1]]

        self.logger.info("=" * 50)
        self.logger.info(f"Download summary: total={len(self.runid_list)}, "
                         f"success={len(success_ids)}, failed={len(failed_ids)}")
        if failed_ids:
            self.logger.warning(f"Failed run IDs: {failed_ids}")

        return success_ids, failed_ids

    # ------------------------------------------------------------------ #
    # Step 2: convert SRA -> FASTQ
    # ------------------------------------------------------------------ #
    def convert_to_fastq(self, runid_list=None):
        """
        Convert downloaded SRA files to paired-end FASTQ (.fastq.gz).

        Args:
            runid_list: Optional subset of run IDs to process.
                        Defaults to all downloaded run IDs.
        """
        targets = runid_list or self.runid_list

        def _process_single(run_id):
            sra_file = self.sra_dir / run_id / f"{run_id}.sra"
            fwd_gz = self.fastq_dir / f"{run_id}_1.fastq.gz"
            rev_gz = self.fastq_dir / f"{run_id}_2.fastq.gz"

            # Skip if both paired gz files already exist and are non-empty
            if (fwd_gz.exists() and fwd_gz.stat().st_size > 0 and
                    rev_gz.exists() and rev_gz.stat().st_size > 0):
                self.logger.info(f"FASTQ already exists, skipping: {run_id}")
                return run_id, True

            # Check SRA source
            if not sra_file.exists() or sra_file.stat().st_size == 0:
                self.logger.warning(f"SRA file missing or empty: {sra_file}")
                return run_id, False

            try:
                # fasterq-dump: produce paired uncompressed FASTQ
                fasterq_cmd = [
                    "fasterq-dump", "--split-files",
                    "--threads", str(DUMP_THREADS),
                    "--outdir", str(self.fastq_dir),
                    str(sra_file),
                ]
                self._run_cmd(fasterq_cmd, "fasterq-dump")

                fwd = self.fastq_dir / f"{run_id}_1.fastq"
                rev = self.fastq_dir / f"{run_id}_2.fastq"

                # Compress both FASTQ files with pigz
                pigz_cmd = ["pigz", "-p", str(DUMP_THREADS), "--force", str(fwd), str(rev)]
                self._run_cmd(pigz_cmd, "pigz")

                # Verify both .fastq.gz files were created
                if (fwd_gz.exists() and fwd_gz.stat().st_size > 0 and
                        rev_gz.exists() and rev_gz.stat().st_size > 0):
                    self.logger.info(f"Converted: {run_id}")
                    return run_id, True
                else:
                    self.logger.warning(f"Compression finished but .fastq.gz missing: {run_id}")
                    return run_id, False

            except subprocess.CalledProcessError as e:
                self.logger.error(f"Conversion failed for {run_id}: {e}")
                self._cleanup_incomplete(run_id)
                return run_id, False
            except Exception as e:
                self.logger.error(f"Unexpected error for {run_id}: {e}")
                self._cleanup_incomplete(run_id)
                return run_id, False

        self.logger.info("=" * 50)
        self.logger.info(f"Converting {len(targets)} SRA files to FASTQ...")
        self.logger.info("=" * 50)

        success_count = 0
        with tqdm(total=len(targets), desc="Convert FASTQ", unit="run") as pbar:
            with ThreadPoolExecutor(max_workers=DUMP_POOL_SIZE) as executor:
                for result in executor.map(_process_single, targets):
                    if result[1]:
                        success_count += 1
                    pbar.update(1)
                    pbar.set_postfix(success=f"{success_count}/{pbar.n}")

        self.logger.info(f"Conversion done: {success_count}/{len(targets)} succeeded")
        return success_count

    def _cleanup_incomplete(self, run_id):
        """Remove incomplete intermediate files after a failed conversion."""
        for pattern in [f"{run_id}_1.fastq", f"{run_id}_2.fastq",
                        f"{run_id}_1.fastq.gz", f"{run_id}_2.fastq.gz"]:
            f = self.fastq_dir / pattern
            if f.exists():
                self.logger.debug(f"Cleaning up incomplete file: {f}")
                f.unlink()

    # ------------------------------------------------------------------ #
    # Pipeline entry point
    # ------------------------------------------------------------------ #
    def run(self):
        """Execute the full download + conversion pipeline."""
        import time
        start = time.time()

        success_ids, _ = self.download_sra()
        self.convert_to_fastq(runid_list=success_ids)

        elapsed = time.time() - start
        h, rem = divmod(elapsed, 3600)
        m, s = divmod(rem, 60)
        self.logger.info("=" * 50)
        self.logger.info(f"Pipeline finished in {int(h)}h {int(m)}m {int(s)}s")


if __name__ == '__main__':
    processor = SRAProcessor()
    processor.run()