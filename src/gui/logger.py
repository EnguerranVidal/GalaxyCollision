from __future__ import annotations

import os
import threading
from datetime import datetime


class SessionLogger:
    def __init__(self, logFilePath: str, threadName: str):
        self.logFilePath = logFilePath
        self.threadIndex = None
        self.threadName = threadName
        self._lock = threading.Lock()
        os.makedirs(os.path.dirname(self.logFilePath) or ".", exist_ok=True)
        with open(self.logFilePath, "w", encoding="utf-8") as f:
            f.write(f"# session log | threadIndex={self.threadIndex}\n")

    def bindThread(self):
        if self.threadIndex is None:
            self.threadIndex = threading.get_ident()
            with open(self.logFilePath, "a", encoding="utf-8") as f:
                f.write(f"# bound threadName={self.threadName} ident={self.threadIndex}\n")

    def log(self, message: str, level: str = "INFO"):
        if self.threadIndex is None:
            self.bindThread()
        timeString = datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
        line = f"{timeString} | idx={self.threadIndex} | {self.threadName} | {level} | {message}\n"
        with self._lock:
            with open(self.logFilePath, "a", encoding="utf-8") as f:
                f.write(line)

    def info(self, message: str):
        self.log(message, "INFO")

    def debug(self, message: str):
        self.log(message, "DEBUG")

    def warning(self, message: str):
        self.log(message, "WARN")

    def error(self, message: str):
        self.log(message, "ERROR")
