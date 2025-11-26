import logging

logger = logging.getLogger()

fmt = "%(asctime)s - %(levelname)s - %(filename)s[:%(lineno)d] - %(message)s"
datefmt = "%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
f_handler = logging.FileHandler("reactions_logs.log", mode='a', encoding='utf-8')
f_handler.setLevel(logging.INFO)
f_handler.setFormatter(formatter)
logger.addHandler(f_handler)