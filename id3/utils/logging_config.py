"""
Logging Configuration Module

Provides colored console logging, structured JSON logging, and performance monitoring.
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Union
from datetime import datetime
import json


class ColoredFormatter(logging.Formatter):
    """Colored log formatter for console output"""

    COLORS = {
        'DEBUG': '\033[36m',    # Cyan
        'INFO': '\033[32m',     # Green
        'WARNING': '\033[33m',  # Yellow
        'ERROR': '\033[31m',    # Red
        'CRITICAL': '\033[35m'  # Magenta
    }
    RESET = '\033[0m'

    def format(self, record):
        # Add color to log level
        levelname = record.levelname
        if levelname in self.COLORS:
            record.levelname = f"{self.COLORS[levelname]}{levelname}{self.RESET}"

        # Format message
        message = super().format(record)

        return message


class StructuredFormatter(logging.Formatter):
    """Structured JSON log formatter (for file output and analysis)"""

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

        # Add context if available
        if hasattr(record, 'context'):
            log_obj['context'] = record.context

        # Add exception information
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
    Configure logging with console and file output

    Args:
        level: Logging level
        log_file: Optional log file path
        console: Whether to enable console output
        colored: Whether to use colored output
        structured: Whether to use structured JSON format (for file output)
        name: Logger name

    Returns:
        Configured logger instance
    """
    # Get or create logger
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Clear existing handlers to avoid duplication
    logger.handlers.clear()

    # Add console handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)

        # Use colored formatter if supported
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

    # Add file handler
    if log_file:
        # Create log directory if it doesn't exist
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(level)

        # Use structured formatter or standard formatter
        if structured:
            file_formatter = StructuredFormatter()
        else:
            file_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(funcName)s:%(lineno)d - %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )

        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    # Disable propagation to prevent duplicate logging
    logger.propagate = False

    return logger


def get_logger(name: str, **kwargs) -> logging.Logger:
    """
    Get or create a logger instance

    Args:
        name: Logger name
        **kwargs: Additional arguments passed to setup_logging

    Returns:
        Logger instance
    """
    # If this is a submodule logger and parent logger exists, use hierarchical logging
    if '.' in name and name.startswith('id3'):
        parent_name = 'id3'
        parent_logger = logging.getLogger(parent_name)

        # If parent logger is configured, return child logger
        if parent_logger.handlers:
            return logging.getLogger(name)

    # Otherwise configure a new logger
    return setup_logging(name=name, **kwargs)


class LogContext:
    """Context manager for adding contextual information to log records"""

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
    """Performance logger for recording execution time and resource usage"""

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
    Configure default logging for ID3 package

    Args:
        debug: Whether to enable debug level logging
    """
    level = logging.DEBUG if debug else logging.INFO

    # Setup main logger
    setup_logging(
        level=level,
        console=True,
        colored=True,
        name='id3'
    )

    # Configure submodule loggers
    logging.getLogger('id3.experiments').setLevel(level)
    logging.getLogger('id3.constraints').setLevel(level)
    logging.getLogger('id3.cai').setLevel(level)
    logging.getLogger('id3.optimizers').setLevel(level)

    # Suppress verbose third-party library logging
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
