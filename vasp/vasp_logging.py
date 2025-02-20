import logging

class vasp_logging:
    def __init__(self, logging_level):
        # 获取根日志记录器
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)  # 根日志记录器保留所有日志信息

        # 创建一个控制台日志处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging_level)  # 设置处理器的日志级别

        # 设置日志格式
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)

        # 添加处理器到根日志记录器
        root_logger.addHandler(console_handler)

        # 如果日志级别无效，记录一条错误日志
        if logging_level not in [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL]:
            root_logger.error(f"无效的日志级别设置：{logging_level}, 采用默认设置 INFO")
            console_handler.setLevel(logging.INFO)