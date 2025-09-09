import logging
import os

# TODO: Change the message for Ea warnings to specify the reaction 
# TODO: Fix the problem of writing in the wrong log file 



class DedupFilter(logging.Filter):
    def __init__(self):
        super().__init__()
        self.seen = set()
    def filter(self, record):
        if record.msg in self.seen:
            return False
        self.seen.add(record.msg)
        return True


def log_warning(log_file_name):
    os.makedirs("warnings", exist_ok=True)

    # use log_file_name as part of the logger name → unique per file
    logger_name = f"{os.path.basename(log_file_name)}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.WARNING)

    if not logger.handlers:
        mode = "w" if logger_name in logging.Logger.manager.loggerDict else "a"
        fh = logging.FileHandler(log_file_name, mode=mode)
        fh.setLevel(logging.WARNING)

        formatter = logging.Formatter(
            "%(asctime)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        fh.setFormatter(formatter)

        fh.addFilter(DedupFilter())
        logger.addHandler(fh)

    return logger


logger_extraction = log_warning("warnings/warning_extraction.log")
logger_retrieval = log_warning("warnings/warning_retrieval.log")
logger_prediction = log_warning("warnings/warning_prediction.log")