import shutil
from configure import *
import pathlib
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


def setup_logging(log_folder, log_name):
    """
    初始化 logging 配置。
    参数:
    log_folder (str): 日志文件存储的文件夹路径。
    log_name (str): 日志文件名。
    """
    if not log_folder.exists():
        log_folder.mkdir()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(str(log_folder / log_name)),
            logging.StreamHandler()
        ]
    )


def get_sample_names_list(path, flag):
    sample_names_list = []
    if flag == "fastp":
        samples = set()
        for sample in path.iterdir():
            samples.add(sample.name.split("_")[0])
        sample_names_list = list(samples)

    elif flag == "megahit":
        for sample in path.iterdir():
            if sample.name.startswith("clean_") and "_1" in sample.name:
                sample_names_list.append(sample.name.split('_')[1])

    elif flag == "prodigal":
        for sample in path.iterdir():
            if sample.name.endswith(".contigs.fa"):
                sample_names_list.append(sample.name.split('_final.contigs.fa')[0])

    elif flag == "diamond":
        for sample in path.iterdir():
            if sample.name.endswith("faa"):
                sample_names_list.append(sample.name.split('_prodigal.faa')[0])

    elif flag == "binning":
        for sample in path.iterdir():
            if sample.name.endswith("_acc_contigs.fa"):
                sample_names_list.append(sample.name.split('_acc_contigs.fa')[0])

    elif flag == "acc_orfs":
        for sample in path.iterdir():
            if sample.name.endswith("_acc_orfs.faa"):
                sample_names_list.append(sample.name.split('_acc_orfs.faa')[0])

    elif flag == "species":
        for sample in path.iterdir():
            if sample.name.endswith("_dastool"):
                sample_names_list.append(sample.name.split('_dastool')[0])

    elif flag == "vf":
        for sample in path.iterdir():
            if sample.name.endswith("_checkm"):
                sample_names_list.append(sample.name.split('_checkm')[0])

    elif flag == "plasme":
        for sample in path.iterdir():
            if sample.name.endswith("_acc_contigs.fa"):
                sample_names_list.append(sample.name.split('_acc_contigs.fa')[0])
    elif flag == "arg_like_orfs":
        for sample in path.iterdir():
            if sample.name.endswith("_arg_like_orfs.faa"):
                sample_names_list.append(sample.name.split('_arg_like_orfs.faa')[0])

    return sample_names_list


def run_command(command, logger=None, shell=False, log_message=None):
    """
    运行命令并捕获异常。

    参数:
    command (list or str): 要执行的命令。
    logger (logging.Logger): 日志记录器实例。
    shell (bool): 是否使用 shell 执行命令。
    log_message (str): 自定义日志信息。
    """
    if logger and log_message:
        logger.info(log_message)

    try:
        result = subprocess.run(
            command,
            check=True,
            shell=shell,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        if result.stderr and logger:
            logger.warning("Command error output:\n{0}".format(result.stderr))
    except subprocess.CalledProcessError as e:
        error_message = "Error: Command failed with exit code {0}\nCommand: {1}\nError output: {2}".format(
            e.returncode, e.cmd, e.stderr
        )
        if logger:
            logger.error(error_message)
        else:
            print(error_message)
        raise e


def process_samples_in_parallel(logger, process_func, sample_names, process_name, max_workers=1):
    """
    Process samples in parallel and monitor progress.

    Args:
        process_func: Function to process a single sample.
        sample_names: List of sample names.
        process_name: Process name for logging.
        max_workers: Maximum number of worker threads.

    Returns:
        list: List of (sample_name, success) tuples.
    """
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_func, sample_name): sample_name
                  for sample_name in sample_names}

        for future in as_completed(futures):
            sample_name = futures[future]
            try:
                success = future.result()
                results.append((sample_name, success))
                if success:
                    logger.info("{0} assembly completed for sample {1}".format(
                        process_name.capitalize(), sample_name))
                else:
                    logger.warning("{0} assembly failed or skipped for sample {1}".format(
                        process_name.capitalize(), sample_name))
            except Exception as e:
                logger.error("Error processing sample {0}: {1}".format(sample_name, str(e)))
                results.append((sample_name, False))

    return results  # ← 关键：必须有 return


def clear_folder_pathlib(folder_path):
    """
    使用 pathlib 清空文件夹
    """
    path = pathlib.Path(folder_path)

    if path.exists() and path.is_dir():
        # 删除文件夹内所有内容
        for item in path.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
        print(f"已清空文件夹: {folder_path}")
    else:
        print(f"文件夹不存在: {folder_path}")