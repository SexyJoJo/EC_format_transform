import logging
import logging.config
import json
import os

class Log:
    def __init__(self, logger):
        if not os.path.exists("./log"): # 用于存放日志
            os.mkdir("./log")

        with open(r"log_config.json", "r") as f:
            log_config = json.load(f)
        logging.config.dictConfig(log_config)
        self.logger = logging.getLogger(logger)


trans_log = Log("trans_logger")
