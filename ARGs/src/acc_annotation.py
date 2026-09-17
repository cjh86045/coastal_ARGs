from utils import setup_logging, run_command, get_sample_names_list, process_samples_in_parallel
from configure import *


class ACCAnnotation:
    def __init__(self, input_folder):
        """
        Initialize the ACC annotation pipeline.

        Args:
            input_folder: Path to the input folder containing preprocessing results
        """
        self.acc_orfs_output_folder = pathlib.Path(input_folder, "preprocessing", "acc_orfs")
        self.arg_orfs_blast_output_folder = pathlib.Path(input_folder, "preprocessing", "ARGs_blast_result_output")

        self.acc_annotate_output_folder = pathlib.Path(input_folder, "acc_annotation")
        self.acc_annotate_output_folder.mkdir(parents=True, exist_ok=True)

        self.ncbi_blast_output_folder = self.acc_annotate_output_folder / "acc_blastp_NR_output"
        self.acc_arg_mge_output_folder = self.acc_annotate_output_folder / "acc_arg_mge_output"
        self.arg_mge_contig_output_folder = self.acc_annotate_output_folder / "acc_arg_mge_contig_output"
        self.sample_mge_output_folder = self.acc_annotate_output_folder / "acc_sample_mge_output"

        # Configure logging
        setup_logging(self.acc_annotate_output_folder, "acc_annotation.log")
        self.logger = logging.getLogger("acc_annotation")

        self.sample_names_list = get_sample_names_list(self.acc_orfs_output_folder, "acc_orfs")
        self.logger.info(f"{len(self.sample_names_list)} Samples starting ACC Annotation processing:")

    def run_diamond_blastp(self):
        """
        Run DIAMOND BLASTP against NCBI NR database for all samples.
        """
        self.ncbi_blast_output_folder.mkdir(parents=True, exist_ok=True)

        def process_sample(sample_name):
            acc_orfs_file = self.acc_orfs_output_folder / f"{sample_name}_acc_orfs.faa"
            if not acc_orfs_file.exists():
                self.logger.warning(f"acc_orfs file not found for {sample_name} skipping")
                return False

            acc_diamond_output_file = self.ncbi_blast_output_folder / f"{sample_name}_acc_blastp_NR.tsv"
            if acc_diamond_output_file.exists():
                self.logger.info(f"acc_diamond_file exists for {sample_name} skipping")
                return True

            # Build DIAMOND BLASTP command
            diamond_command = [
                "diamond", "blastp",
                "--query", str(acc_orfs_file),  # Input ORFs file
                "--db", str(NCBI_NR_DB),  # NCBI NR database
                "--out", str(acc_diamond_output_file),  # Output result file
                "--evalue", "1e-7",  # E-value threshold
                "--id", "70",  # Sequence identity threshold
                "--query-cover", "70",  # Query coverage threshold
                "--max-target-seqs", "1",  # Limit max target sequences per query
                "--outfmt", "6", "qseqid", "sseqid",  # Output format
                "--threads", str(THREAD_CONFIG["diamond"]),  # Use configured threads
                "--quiet",
            ]
            # Run DIAMOND BLASTP command
            run_command(diamond_command, logger=self.logger,
                        log_message=f"Running DIAMOND BLASTP for sample {sample_name}: {' '.join(diamond_command)}")
            return True

        process_samples_in_parallel(self.logger, process_sample, self.sample_names_list,
                                    "acc_diamond_blastp", THREAD_CONFIG["diamond_pool"])

        self.logger.info("DIAMOND BLASTP against NCBI NR database completed.")

    def _load_prodigal_coords(self, sample_name):
        """
        Parse Prodigal coordinate file to extract ORF positions.

        Args:
            sample_name: Name of the sample

        Returns:
            DataFrame with ORF coordinates or None if file missing
        """
        fasta_file = self.acc_orfs_output_folder / f"{sample_name}_acc_orfs.faa"
        if not fasta_file.exists():
            self.logger.warning(f"FASTA file missing: {fasta_file}")
            return None

        records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            parts = record.description.split(" # ")
            records.append({
                "orfs_name": parts[0],
                "Start": int(parts[1]),
                "End": int(parts[2]),
                "Strand": int(parts[3])
            })
        return pd.DataFrame(records)

    def _merge_blast_results(self, sample_name: str, coord_df: pd.DataFrame):
        """
        Merge NR and ARG BLAST results with coordinates.

        Args:
            sample_name: Name of the sample
            coord_df: DataFrame with ORF coordinates

        Returns:
            Merged DataFrame or None if BLAST results missing
        """
        nr_file = self.ncbi_blast_output_folder / f"{sample_name}_acc_blastp_NR.tsv"
        arg_file = self.arg_orfs_blast_output_folder / f"{sample_name}_card_results.tsv"

        if not nr_file.exists() or not arg_file.exists():
            self.logger.warning(f"BLAST results missing for {sample_name}")
            return None

        nr_df = pd.read_csv(nr_file, sep="\t", names=["orfs_name", "nr_indexes"])
        arg_df = pd.read_csv(arg_file, sep="\t", names=["orfs_name", "args_annotation"])

        return (
            coord_df
            .merge(nr_df, on="orfs_name", how="left")
            .merge(arg_df, on="orfs_name", how="left")
        )

    def _add_annotations(self, df: pd.DataFrame, nr_index) -> pd.DataFrame:
        """
        Add NR descriptions and MGE flags to the merged DataFrame.

        Args:
            df: Merged DataFrame with BLAST results
            nr_index: NCBI NR sequence index for fetching descriptions

        Returns:
            DataFrame with annotations added
        """
        df["final_annotation"] = df["args_annotation"].fillna(df["nr_indexes"])
        df["arg_flag"] = ~df["args_annotation"].isna()

        # Batch query NR descriptions
        mask = ~df["arg_flag"] & df["nr_indexes"].notna()
        unique_ids = df.loc[mask, "nr_indexes"].unique()
        desc_map = {}
        for nr_id in unique_ids:
            try:
                record = nr_index.get(nr_id)
                desc_map[nr_id] = record.description if record else None
            except Exception as e:
                self.logger.warning(f"Failed to fetch {nr_id}: {e}")
                desc_map[nr_id] = None

        df.loc[mask, "nr_description"] = df.loc[mask, "nr_indexes"].map(desc_map)

        # Mark MGEs based on keywords
        df["MGE_flag"] = df["nr_description"].apply(
            lambda x: any(kw.lower() in str(x).lower() for kw in MGEs_keywords) if pd.notna(x) else False
        )

        return df.drop(columns=["nr_indexes", "args_annotation"])

    def _calculate_mge_to_arg_gene_distance(self, merged_df):
        """
        Calculate whether MGEs are within 5 genes of ARGs on the same contig.

        This method replaces the previous distance-based approach with a more
        biologically meaningful gene-interval-based approach. An MGE is considered
        "near" an ARG if there are 5 or fewer genes between them on the same contig.

        Args:
            merged_df: DataFrame containing:
                - orfs_name: ORF name (e.g., k127_1417741_1)
                - contig_name: Contig name (extracted from orfs_name)
                - Start: ORF start position
                - End: ORF end position
                - arg_flag: Whether it's an ARG annotation (boolean)
                - MGE_flag: Whether it's an MGE (boolean)

        Returns:
            DataFrame with two additional columns:
                - mge_near_arg_flag: Whether MGE is within 5 genes of an ARG
                - distance_to_arg_in_genes: Gene interval distance to nearest ARG
        """
        # Sort by contig and start position to establish gene order
        merged_df = merged_df.sort_values(['contig_name', 'Start']).reset_index(drop=True)

        # Add gene order index within each contig
        merged_df['gene_order'] = merged_df.groupby('contig_name').cumcount()

        # Initialize new columns
        merged_df['mge_near_arg_flag'] = False
        merged_df['distance_to_arg_in_genes'] = np.nan

        # Process each contig separately
        for contig_name, contig_df in merged_df.groupby('contig_name'):
            # Get gene order indices of all ARGs in this contig
            arg_orders = contig_df.loc[contig_df['arg_flag'], 'gene_order'].tolist()

            if not arg_orders:  # No ARGs in this contig, skip
                continue

            # Check each MGE in this contig
            for idx, row in contig_df[contig_df['MGE_flag']].iterrows():
                mge_order = row['gene_order']

                # Calculate gene interval distances to all ARGs
                distances = [abs(mge_order - arg_order) for arg_order in arg_orders]
                min_distance = min(distances)

                # Check if within 5 genes (distance ≤ 5, including adjacent)
                if min_distance <= 5:
                    merged_df.at[idx, 'mge_near_arg_flag'] = True
                    merged_df.at[idx, 'distance_to_arg_in_genes'] = min_distance

        # Remove helper column
        merged_df = merged_df.drop(columns=['gene_order'])

        return merged_df

    def get_annotation(self):
        """
        Main annotation process: merge coordinates, NR and ARG results,
        add annotations, and calculate MGE-ARG proximity.
        """
        nr_database_index = None
        try:
            nr_database_index = SeqIO.index_db(str(NCBI_NR_DB_index), str(NCBI_NR_DB_fa), "fasta")
            self.acc_arg_mge_output_folder.mkdir(parents=True, exist_ok=True)

            for sample_name in self.sample_names_list:
                merge_result = self.acc_arg_mge_output_folder / f"{sample_name}_arg_mge.tsv"
                if merge_result.exists():
                    self.logger.info(f"merge_result file already exists for {sample_name} skipping")
                    continue

                coord_df = self._load_prodigal_coords(sample_name)
                if coord_df is None:
                    continue

                merged_df = self._merge_blast_results(sample_name, coord_df)
                if merged_df is None:
                    continue

                merged_df = self._add_annotations(merged_df, nr_database_index)

                # Extract contig name from ORF name
                merged_df["contig_name"] = merged_df["orfs_name"].str.extract(r'^(.*)_\d+$')

                # Calculate MGE-ARG proximity using gene interval approach
                merged_df = self._calculate_mge_to_arg_gene_distance(merged_df)

                # Save results
                merged_df.to_csv(str(merge_result), sep="\t", index=False)

        finally:
            # Close the database index connection
            if nr_database_index:
                nr_database_index.close()

    def get_arg_mge_contigs(self):
        """
        Generate contig-level and sample-level MGE summaries.
        """
        # Create an empty DataFrame for summary information
        summary_data = []

        self.arg_mge_contig_output_folder.mkdir(parents=True, exist_ok=True)
        self.sample_mge_output_folder.mkdir(parents=True, exist_ok=True)

        for sample_name in self.sample_names_list:
            merge_result = self.acc_arg_mge_output_folder / f"{sample_name}_arg_mge.tsv"
            if not merge_result.exists():
                self.logger.info(f"merge_result file not found for {sample_name}, skipping")
                continue

            merge_result_df = pd.read_csv(merge_result, sep="\t")

            # Calculate MGE count for current sample
            mge_count = merge_result_df['MGE_flag'].sum()
            acc_numbers = merge_result_df['orfs_name'].nunique()
            summary_data.append({'sample_name': sample_name, 'MGE_count': mge_count, "acc_numbers": acc_numbers})

            # Filter contigs containing at least one MGE
            mge_contigs = merge_result_df[merge_result_df['MGE_flag'] == True]['contig_name'].unique()
            filtered_df = merge_result_df[merge_result_df['contig_name'].isin(mge_contigs)]

            # Define processing function for each contig
            def process_contig(group_df):
                records = []
                mge_records = []

                for _, row in group_df.iterrows():
                    record = f"{row['final_annotation']}({row['Start']},{row['End']},{row['Strand']})"
                    records.append(record)

                    # Only include MGEs that are within 5 genes of an ARG
                    if row['MGE_flag'] and row['mge_near_arg_flag'] and pd.notna(row['distance_to_arg_in_genes']):
                        mge_records.append(
                            f"MGE:{row['final_annotation']}:{int(row['distance_to_arg_in_genes'])}"
                        )

                contig_str = ",".join(records)
                if mge_records:
                    contig_str += "," + ",".join(mge_records)
                return contig_str

            # Process contig data
            result = (
                filtered_df.groupby('contig_name')
                .agg(annotation_info=('final_annotation', lambda g: process_contig(filtered_df.loc[g.index])))
                .reset_index()
            )
            # Save contig-level detailed results
            arg_mge_contig_file = self.arg_mge_contig_output_folder / f"{sample_name}_arg_mge_contig.tsv"
            result.to_csv(arg_mge_contig_file, sep="\t", index=False)

        # Save summary file
        if summary_data:  # Only save if there's data
            summary_df = pd.DataFrame(summary_data)
            summary_file = self.sample_mge_output_folder / "sample_mge_summary.tsv"
            summary_df.to_csv(summary_file, sep="\t", index=False)
            self.logger.info(f"MGE summary saved to {summary_file}")

    def run(self):
        """
        Run the complete ACC annotation pipeline.
        """
        self.logger.info("Starting the entire ACC ANNOTATION pipeline...")
        self.run_diamond_blastp()
        self.get_annotation()
        self.get_arg_mge_contigs()

        self.logger.info("Entire preprocessing pipeline completed.")


if __name__ == '__main__':
    # Initialize and run the pipeline
    acc_annotation = ACCAnnotation("")
    acc_annotation.run()