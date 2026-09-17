from utils import setup_logging, run_command, get_sample_names_list, process_samples_in_parallel
from pathlib import Path
from configure import *  # Import configuration from configure.py
import pandas as pd
import logging


class VirulenceFactorAnnotation(object):
    def __init__(self, input_folder):
        self.input_folder = input_folder
        self.checkm_output_folder = Path(self.input_folder, "binning", "checkm_out")
        self.arb_info_output_folder = Path(self.input_folder, "species_annotation", "bin_count_sample_results")

        self.vf_annotation_folder = Path(input_folder, "vf_annotation")
        self.vf_annotation_folder.mkdir(parents=True, exist_ok=True)

        self.vf_annotation_output_folder = Path(self.vf_annotation_folder, "vf_annotation_diamond_results")
        self.vf_number_samples_count_output_folder = Path(self.vf_annotation_folder, "vf_number_samples_count_results")

        setup_logging(self.vf_annotation_folder, "vf_annotation.log")
        self.logger = logging.getLogger("vf_annotation")
        self.sample_names_list = get_sample_names_list(self.checkm_output_folder, "vf")

    def run_diamond_vf(self):
        self.logger.info("Running diamond vf annotation")
        self.vf_annotation_output_folder.mkdir(parents=True, exist_ok=True)

        # Read the combined ARB info result file
        arb_info_result = self.arb_info_output_folder / "bin_count_sample_results.tsv"  # Assuming this is the filename
        if not arb_info_result.exists():
            self.logger.error("No combined ARB info result file found!")
            return

        # Read the entire TSV file
        df_all = pd.read_csv(str(arb_info_result), sep='\t')

        def process_sample(sample_name):
            # Filter data for the current sample from the combined file
            sample_data = df_all[df_all['sample'] == sample_name]
            if sample_data.empty:
                self.logger.warning(f"No data found for sample {sample_name}")
                return None

            # Get all bin names for this sample
            bin_names = sample_data['arg_bin_names'].iloc[0].split(',')

            bins_orf_file = self.checkm_output_folder / f"{sample_name}_checkm" / "bins"
            for bin in bin_names:
                bin_orf_file = bins_orf_file / bin / "genes.faa"
                if not bin_orf_file.exists():
                    self.logger.warning(f"No genes.faa file found for {sample_name}:{bin}")
                    continue

                vf_blastp_output_folder = self.vf_annotation_output_folder / f"{sample_name}_vf_result"
                vf_blastp_output_folder.mkdir(parents=True, exist_ok=True)

                result_file = vf_blastp_output_folder / f"{bin}_vf.tsv"
                if result_file.exists():
                    self.logger.info(f"VF result already exists. Skipping {sample_name}:{bin}")
                    continue

                diamond_command = [
                    "diamond", "blastp",
                    "-d", str(VFDB_PRO),
                    "-q", str(bin_orf_file),
                    "-o", str(result_file),
                    "--threads", str(THREAD_CONFIG["diamond"]),
                    "--evalue", "1e-7",
                    "--id", "70",
                    "--query-cover", "70",
                    "--more-sensitive",
                    "--max-hsps", "1",
                    "--max-target-seqs", "1",
                    "--outfmt", "6",
                    "--quiet",
                ]

                run_command(diamond_command, logger=self.logger,
                            log_message=f"Running DIAMOND for sample {sample_name}: {' '.join(diamond_command)}")
            return True

        process_samples_in_parallel(self.logger, process_sample, self.sample_names_list, "diamond",
                                    THREAD_CONFIG["diamond_pool"])
        self.logger.info("DIAMOND BLASTP completed.")

    def vf_number_samples_count(self):
        self.logger.info("Running vf number samples count")

        self.vf_number_samples_count_output_folder.mkdir(parents=True, exist_ok=True)
        vf_count_list = []
        bins_count_list = []

        for sample_name in self.sample_names_list:
            vf_blastp_output_folder = self.vf_annotation_output_folder / f"{sample_name}_vf_result"
            if not vf_blastp_output_folder.exists():
                self.logger.warning(f"No VF blastp output folder found, skipping {sample_name}")
                vf_count_list.append(0)
                bins_count_list.append(0)
                continue

            vf_count = 0
            bins_count = 0

            for file in vf_blastp_output_folder.glob("*_vf.tsv"):
                try:
                    # Check if the file is empty
                    if file.stat().st_size == 0:
                        self.logger.debug(f"Empty VF result file: {file}")
                        bins_count += 1  # Still count the bin
                        continue

                    # Try to read the file
                    df = pd.read_csv(str(file), sep='\t', header=None)
                    if not df.empty:
                        vf_count += len(df)
                    bins_count += 1

                except Exception as e:
                    self.logger.error(f"Error processing {file}: {str(e)}")
                    bins_count += 1  # Count the bin even if there's an error
                    continue

            vf_count_list.append(vf_count)
            bins_count_list.append(bins_count)
            self.logger.info(f"Sample {sample_name}: {vf_count} VFs found in {bins_count} bins")

        # Create result DataFrame
        result = pd.DataFrame({
            "sample": self.sample_names_list,
            "vf_count": vf_count_list,
            "bins_number": bins_count_list,
            "vf_per_bin": [round(vf / bins, 4) if bins > 0 else 0
                           for vf, bins in zip(vf_count_list, bins_count_list)]
        })

        # Save results
        output_file = self.vf_number_samples_count_output_folder / "vf_number_samples_count.tsv"
        result.to_csv(output_file, sep='\t', index=False)
        self.logger.info(f"VF count results saved to {output_file}")

    def run(self):
        # self.run_diamond_vf()
        self.vf_number_samples_count()


if __name__ == '__main__':
    vf_annotation = VirulenceFactorAnnotation("")
    vf_annotation.run()