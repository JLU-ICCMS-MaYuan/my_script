import logging

class qe_log:
    # logging.basicConfig(level = logging.DEBUG,format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # logger = logging.getLogger(__name__)
    # console_handler = logging.StreamHandler()
    # logger.addHandler(console_handler)
    def __init__(self, logging_level):
        logger = logging.getLogger(__name__)

        # 默认只在屏幕上输出INFO及以上级别的日志
        if logging_level == "DEBUG":
            logger.setLevel(logging.DEBUG)
        elif logging_level == "INFO":
            logger.setLevel(logging.INFO)
        elif logging_level == "WARNING":
            logger.setLevel(logging.WARNING)
        elif logging_level == "ERROR":
            logger.setLevel(logging.ERROR)
        elif logging_level == "CRITICAL":
            logger.setLevel(logging.CRITICAL)
        else:
            logger.error(f"无效的日志级别设置：{logging_level}, 采用默认设置 INFO")
            logger.setLevel(logging.INFO)
