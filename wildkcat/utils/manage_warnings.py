import logging
# import os

# TODO: Change the message for Ea warnings to specify the reaction 


class DedupFilter(logging.Filter):
    def __init__(self):
        super().__init__()
        self.seen = set()
    def filter(self, record):
        if record.msg in self.seen:
            return False
        self.seen.add(record.msg)
        return True

