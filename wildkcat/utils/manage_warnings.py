import logging
import os


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
    logger = logging.getLogger("warning_logger")
    logger.setLevel(logging.WARNING)

    if not logger.handlers:
        fh = logging.FileHandler(log_file_name, mode="w") 
        fh.setLevel(logging.WARNING)

        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        fh.setFormatter(formatter)

        fh.addFilter(DedupFilter())
        logger.addHandler(fh)

    return logger


logger_extraction = log_warning("warning_extraction.log")
logger_retrieval = log_warning("warning_retrieval.log")
logger_prediction = log_warning("warning_prediction.log")