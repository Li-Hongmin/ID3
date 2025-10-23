"""



"""

import logging
import sys
from pathlib import Path
from typing import Optional, Union
from datetime import datetime
import json


class ColoredFormatter(logging.Formatter):

    

    COLORS = {





    }
    RESET = '\033[0m'
    
    def format(self, record):

        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{self.RESET}"
        

        message = super().format(record)
        
        return message


class StructuredFormatter(logging.Formatter):
    """结构化JSON日志格式化器（用于文件输出和分析）"""
    
    def format(self, record):
        log_obj = {
            'timestamp': datetime.utcnow().isoformat(),
            'level': record.levelname,
            'logger': record.name,
            'message': record.getMessage(),
            'module': record.module,
            'function': record.funcName,
            'line': record.lineno,
        }
        

        if hasattr(record, 'context'):
            log_obj['context'] = record.context
            

        if record.exc_info:
            log_obj['exception'] = self.formatException(record.exc_info)
        
        return json.dumps(log_obj, ensure_ascii=False)


def setup_logging(
    level: Union[int, str] = logging.INFO,
    log_file: Optional[str] = None,
    console: bool = True,
    colored: bool = True,
    structured: bool = False,
    name: str = 'id3'
) -> logging.Logger:
    """

    
    Args:






        
    Returns:

    """

    logger = logging.getLogger(name)
    logger.setLevel(level)
    

    logger.handlers.clear()
    

    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        

        if colored and sys.stdout.isatty():
            console_formatter = ColoredFormatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
        else:
            console_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
        
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
    

    if log_file:

        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(level)
        

        if structured:
            file_formatter = StructuredFormatter()
        else:
            file_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(funcName)s:%(lineno)d - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
        
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    

    logger.propagate = False
    
    return logger


def get_logger(name: str, **kwargs) -> logging.Logger:
    """

    
    Args:


        
    Returns:

    """

    if '.' in name and name.startswith('id3'):
        parent_name = 'id3'
        parent_logger = logging.getLogger(parent_name)
        

        if parent_logger.handlers:
            return logging.getLogger(name)
    

    return setup_logging(name=name, **kwargs)


class LogContext:

    
    def __init__(self, logger: logging.Logger, **context):
        self.logger = logger
        self.context = context
        self._old_factory = None
    
    def __enter__(self):
        self._old_factory = logging.getLogRecordFactory()
        
        def record_factory(*args, **kwargs):
            record = self._old_factory(*args, **kwargs)
            record.context = self.context
            return record
        
        logging.setLogRecordFactory(record_factory)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        logging.setLogRecordFactory(self._old_factory)


class PerformanceLogger:
    """性能日志记录器，用于记录执行时间和资源使用"""
    
    def __init__(self, logger: logging.Logger, operation: str):
        self.logger = logger
        self.operation = operation
        self.start_time = None
        
    def __enter__(self):
        self.start_time = datetime.now()
        self.logger.debug(f"Starting {self.operation}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        elapsed = (datetime.now() - self.start_time).total_seconds()
        
        if exc_type is None:
            self.logger.info(f"Completed {self.operation} in {elapsed:.3f}s")
        else:
            self.logger.error(f"Failed {self.operation} after {elapsed:.3f}s: {exc_val}")



def configure_default_logging(debug: bool = False):
    """

    
    Args:

    """
    level = logging.DEBUG if debug else logging.INFO
    

    setup_logging(
        level=level,
        console=True,
        colored=True,
        name='id3'
    )
    

    logging.getLogger('id3.experiments').setLevel(level)
    logging.getLogger('id3.constraints').setLevel(level)
    logging.getLogger('id3.cai').setLevel(level)
    logging.getLogger('id3.optimizers').setLevel(level)
    

    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('torch').setLevel(logging.WARNING)



__all__ = [
    'setup_logging',
    'get_logger',
    'LogContext',
    'PerformanceLogger',
    'configure_default_logging',
    'ColoredFormatter',
    'StructuredFormatter',
]